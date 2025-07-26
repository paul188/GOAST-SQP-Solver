#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <goast/SQP/Utils/ObjectFactory.h>

#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <goast/Smoothers.h>
#include <unordered_set>
#include <fstream>
#include <iostream>

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

/** 
 * \brief Simulation of a simple fold 
 * \author Johannssen
 *
 * We optimize this energy by direct optimization via gradient descent or BFGS.
 * 
 */

/**/
int main(int argc, char *argv[])
{
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/paperCrissCross.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial;
    getGeometry( plate, plateGeomInitial);
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );

    std::vector<int> bdryMask;
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        if(coords[0] < 0.04)
        {
            bdryMask.push_back(i);
        }
        if(coords[0] > 0.96)
        {
            bdryMask.push_back(i);
        }
    }

    //Next, set the bending edge weights for the edge at x = 0.5
    VectorType edge_weights = VectorType::Ones(plateTopol.getNumEdges());

    for(int edgeIdx = 0; edgeIdx < plateTopol.getNumEdges(); edgeIdx++){
        int node_i = plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
        int node_j = plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

        VecType coords_i, coords_j;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords_i, node_i);
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords_j, node_j);

        if(coords_i[0] == 0.5 && coords_j[0] == 0.5){
            edge_weights[edgeIdx] = 1.0;
        }
    }

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMask );

    SimpleBendingEnergy<DefaultConfigurator> E_bend( plateTopol, plateGeomRef, true , edge_weights);

    OpenMesh::IO::read_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/bendingFoldSol_withNewton2.ply");
    getGeometry( plate, plateGeomDef );
    RealType energy;
    E_bend.apply(plateGeomDef, energy);

    std::string filename = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/SimpleBendingEnergy.log";
    std::ifstream infile(filename);
    std::vector<double> values;
    std::string line;

    if (!infile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }

    while (std::getline(infile, line)) {
        try {
            double val = std::stod(line);
            values.push_back(val);
        } catch (const std::exception& e) {
            std::cerr << "Skipping invalid line: " << line << std::endl;
        }
    }

    infile.close();

    // Optional: print the result
    std::cout << "Read " << values.size() << " values." << std::endl;
    std::cout<< " Number of edges for comparison: "<<plateTopol.getNumEdges() << std::endl;

    for(int i = 0; i < values.size(); i++)
    {
        std::cout << "Value " << i << ": " << values[i] << std::endl;
    }

    std::ofstream stream;
    stream.open("/home/s24pjoha_hpc/goast_old_old/goast/build/examples/VertexValues.log");

    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        TriMesh::VertexHandle vh(i);
        RealType value_vertex = 0;
        int num_neighbouring_edges = 0;
        for (auto ve_it = plate.ve_iter(vh); ve_it.is_valid(); ++ve_it) {
            OpenMesh::EdgeHandle eh = *ve_it;
            int edgeIdx = eh.idx();
            value_vertex += values[edgeIdx];
            num_neighbouring_edges++;
        }
        stream <<(value_vertex/num_neighbouring_edges) << std::endl;
    }

    stream.close();

    return 0;
}