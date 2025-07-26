// TO BE DONE

#include <cmath>
#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <iostream>
#include <unordered_set>

#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <goast/Smoothers.h>

#include <goast/SQP/Algorithm/CostFunctional.h>
#include <goast/SQP/Utils/SparseMat.h>
#include <goast/SQP/DOFHandling/FoldDofs.h>

typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;
typedef DefaultConfigurator::SparseMatrixType MatrixType;
typedef DefaultConfigurator::FullMatrixType FullMatrixType;

bool is_near(RealType coord, RealType value, RealType tol = 1e-6)
{
    return std::abs(coord - value) < tol;
}


int main(int argc, char *argv[])
{

try{
    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE FOLD OPTIMIZATION OF A ROTATING CROSS FOLD" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    Eigen::setNbThreads(8);  // Adjust based on CPU cores

    // load flat plate [0,1]^2
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/paperCrissCross.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial;
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );
    getGeometry( plate, plateGeomInitial );

    // bdryMaskRef_1 fixes x,y,z coordinates in the reference geometry
    std::vector<int> bdryMaskRef_1;
    // bdryMaskRef_2 fixes y,z coordinates in the reference geometry
    std::vector<int> bdryMaskRef_2;
    // bdryMaskRef_3 fixes x,z coordinates in the reference geometry
    std::vector<int> bdryMaskRef_3;

    RealType t_init = 0.0;
    RealType tolerance = 1e-4;

    VectorType edge_weights = VectorType::Ones(plateTopol.getNumEdges());

    // Just for checking, we will cover the fold edges red
    plate.request_edge_colors();
    plate.request_vertex_colors();

    for(int edgeIdx = 0; edgeIdx < plateTopol.getNumEdges(); edgeIdx++){
        int node_i = plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
        int node_j = plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

        VecType coords_i, coords_j;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords_i, node_i);
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords_j, node_j);
        //std::cout<<std::setprecision(12)<<"coords0: "<<coords_i[0]<<" "<<coords_j[0]<<" "<<std::endl;
        //std::cout<<std::setprecision(12)<<"coords1: "<<coords_i[1]<<" "<<coords_j[1]<<" "<<std::endl;
        // Set the edge weights for the edge at x = 0.5 to zero -> fold
        if((is_near(coords_i[1],0.5, 1e-4) && is_near(coords_j[1], 0.5, 1e-4)) ||(is_near(coords_i[0],0.5, 1e-4) && is_near(coords_j[0],0.5, 1e-4))){
            edge_weights[edgeIdx] = 0;
            // Set colors for both vertices of the edge
            plate.set_color(plate.vertex_handle(node_i), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
            plate.set_color(plate.vertex_handle(node_j), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
        }
    }

    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomInitial, coords, i);
        // first, fix the outer vertices of the square
        if( (std::abs(coords[0] < tolerance) || std::abs(coords[0] - 1.0) < tolerance) && (std::abs(coords[1] < tolerance) || std::abs(coords[1] - 1.0) < tolerance) )
        {
            bdryMaskRef_1.push_back( i );            continue;
        }

        // fix the middle vertex
        if(is_near(coords[0], 0.5, tolerance) && is_near(coords[1], 0.5, tolerance))
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }

        // Would like the possibility to let the vertices vary along possibly rotated fold lines
        // how to do this ? Difficult
        if(is_near(coords[0], 0.5, tolerance))
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }

        if(is_near(coords[1], 0.5, tolerance))
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }

        if(is_near(coords[0], 1.0) || is_near(coords[0], 0.0))
        {
            bdryMaskRef_3.push_back(i);
            continue;
        }

        if(is_near(coords[1], 1.0) || is_near(coords[1], 0.0))
        {
            bdryMaskRef_2.push_back(i);
            continue;
        }

    }

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef_1 );
    std::vector<int> activeRef_2 = (std::vector<int>){0,1,1};
    std::vector<int> activeRef_3 = (std::vector<int>){1,0,1};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_2 , activeRef_2);
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_3 , activeRef_3);

    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntries(bdryMaskRef_1.begin(), bdryMaskRef_1.end());
    uniqueEntries.insert(bdryMaskRef_2.begin(), bdryMaskRef_2.end());
    uniqueEntries.insert(bdryMaskRef_3.begin(), bdryMaskRef_3.end());

    // Move the unique elements back to bdryMaskRef_1 in order
    bdryMaskRef_1.assign(uniqueEntries.begin(), uniqueEntries.end());
    std::sort(bdryMaskRef_1.begin(), bdryMaskRef_1.end());

    std::vector<int> foldVertices;

    auto foldDofs = FoldDofsCrossInterpolation<DefaultConfigurator>(plateTopol,plateGeomInitial, plateGeomRef, bdryMaskRef_1);
    foldDofs.getFoldVertices(foldVertices);
    auto DfoldDofs = FoldDofsCrossInterpolationGradient<DefaultConfigurator>(plateTopol, bdryMaskRef_1, plateGeomInitial, foldVertices);

    // TEST THE FOLD DOFS ACTION
    VectorType translate_vec = VectorType::Ones(foldDofs.getNumDofs());
    translate_vec[0] = 5.0;
    translate_vec[1] = 5.0;
    VectorType testGeom = plateGeomInitial;
    MatrixType testMat;
    foldDofs.apply(translate_vec,testGeom);
    DfoldDofs.apply(translate_vec,testMat);
    setGeometry(plate, testGeom);
    OpenMesh::IO::write_mesh(plate, "deformed_mesh_test_new.ply");

    // Get the coords of the point with index 0
    VecType coords;
    getXYZCoord<VectorType, VecType>(plateGeomInitial, coords,0);
    std::cout<<"Point zero: "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<std::endl;
    
    VectorValuedDerivativeTester<DefaultConfigurator> tester3(foldDofs, DfoldDofs, 0.005, testMat.rows());
    tester3.plotAllDirections(translate_vec, "deriv_test/");
} 
catch ( BasicException &el ){
      std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
}

return 0;
}
