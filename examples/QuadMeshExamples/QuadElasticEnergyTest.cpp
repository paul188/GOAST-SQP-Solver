#include <goast/Core.h>
#include "goast/Core/Auxiliary.h"
#include "goast/Developability/DevelopabilityCentroidEnergy.h"
#include <iostream>
#include <ctime>
#include <string>
#include <vector>

#include <goast/QuadMesh/QuadTopology.h>
#include <goast/Developability/Developability.h>
//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

int main(int argc, char *argv[])
{

try{

    MyMesh mesh;

    std::cout<<"Hello10"<<std::endl;
    
    OpenMesh::IO::read_mesh(mesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/quadTest.ply");

    VectorType geom;

    QuadMeshTopologySaver quadTopol(mesh);

    size_t num_vertices = quadTopol.getNumVertices();

    QuadMeshTopologySaver::getGeometry(mesh,geom);
    VectorType quadGeomRef = geom;

    // set the boundary
    std::vector<int> bdryMask;

    for(int i = 0; i < num_vertices; i++){
        VecType coords;
        coords[0] = geom[3*i];
        if(coords[0] == 0.0)
        {
            bdryMask.push_back(i);
        }
    }

    std::cout<<"Hello11"<<std::endl;

    // Initialize all the centroid stuff
    TriMesh centroid_base_mesh = quadTopol.makeQuadMeshCentroid();
    MeshTopologySaver centroidTopol(centroid_base_mesh);
    VectorType centroidReferenceGeometry;
    getGeometry(centroid_base_mesh, centroidReferenceGeometry);
    OpenMesh::IO::write_mesh(centroid_base_mesh, "centroid_base_coarse.ply");

    std::cout<<"Hello12"<<std::endl;

    VectorType _factors_elasticity_dev = VectorType::Ones(2);
    _factors_elasticity_dev[0] = 1.0;
    _factors_elasticity_dev[1] = 10000.0;

    VectorType _factors_mem_bend = VectorType::Ones(2);
    _factors_mem_bend[0] = 10000.0;
    _factors_mem_bend[1] = 1.0;
    VectorType test_input = centroidReferenceGeometry;
    /*
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-0.05, 0.05); // range: [-0.05, 0.05]

    
    for (int i = 0; i < num_vertices; i++)
    {
        VectorType coords;
        getXYZCoord(test_input, coords, i);
        coords[0] += distribution(generator); // x
        coords[1] += distribution(generator); // y
        coords[2] += distribution(generator); // z
        setXYZCoord(test_input, coords, i);
    }
        */
    
    std::cout<<"Hello1"<<std::endl;

    RealType dest;
    VectorType Dest;

    QuadElasticEnergy<DefaultConfigurator> quadElasticEnergy(quadTopol, centroidTopol, quadGeomRef, centroidReferenceGeometry,_factors_elasticity_dev, _factors_mem_bend);
    QuadElasticEnergyGradient<DefaultConfigurator> quadElasticEnergyGrad(quadTopol, centroidTopol, quadGeomRef, centroidReferenceGeometry, _factors_elasticity_dev, _factors_mem_bend);

    std::cout<<"Hello2"<<std::endl;

    quadElasticEnergy.apply(test_input, dest);
    quadElasticEnergyGrad.apply(test_input, Dest);

    std::cout<<"Energy: "<<dest<<std::endl;
    std::cout<<"Gradient size: "<<Dest.size()<<std::endl;

    ScalarValuedDerivativeTester<DefaultConfigurator> tester(quadElasticEnergy,quadElasticEnergyGrad,0.01,50);
    //VectorValuedDerivativeTester<DefaultConfigurator> tester(constraint, constraintGrad,0.01, Dest.rows());
    tester.plotAllDirections(test_input,"/lustre/scratch/data/s24pjoha_hpc-results/deriv_test/");

}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}
return 0;
}