//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACK

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
#include <goast/SQP.h>

typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;
typedef DefaultConfigurator::SparseMatrixType MatrixType;
typedef DefaultConfigurator::FullMatrixType FullMatrixType;

int main(int argc, char *argv[])
{

try{
    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE FOLD OPTIMIZATION" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    Eigen::setNbThreads(8);  // Adjust based on CPU cores

    // load flat plate [0,1]^2
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/testMesh2.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial;
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );
    getGeometry( plate, plateGeomInitial );

    // determine boundary mask for reference geometry and foldVertices
    // the part of the reference boundary where everything is fixed
    std::vector<int> bdryMaskRef_1;
    std::vector<int> bdryMaskRef_2;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
        if( std::abs(coords[1]) < 1e-6 || std::abs(coords[1] - 4.0) < 1e-10 ){
            bdryMaskRef_1.push_back( i );
        }
        if(std::abs(coords[1] - 2.0) < 1e-6){
            coords[1]+=0.3;
            bdryMaskRef_1.push_back( i );
        }
        if( std::abs(coords[0]) < 1e-6 || std::abs(coords[0] - 1.0) < 6 ){
            bdryMaskRef_2.push_back( i );
        }

        setXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
    }

    setGeometry(plate,plateGeomRef);
    OpenMesh::IO::write_mesh(plate,"reference_mesh.ply");

    // extend all boundary masks to (x,y,z) coordinates

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef_1 );
    std::vector<int> activeRef_2 = (std::vector<int>){1,0,1};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_2 , activeRef_2);

    std::vector<VectorType> def_geometries;
    std::vector<VectorType> ref_geometries;
    std::vector<VectorType> fold_DOFs;

    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntries(bdryMaskRef_1.begin(), bdryMaskRef_1.end());
    uniqueEntries.insert(bdryMaskRef_2.begin(), bdryMaskRef_2.end());

    // Move the unique elements back to bdryMaskRef_1 in order
    bdryMaskRef_1.assign(uniqueEntries.begin(), uniqueEntries.end());
    std::sort(bdryMaskRef_1.begin(), bdryMaskRef_1.end());

     // determine boundary mask for optimization
     // and deform part of boundary
    std::vector<int> bdryMaskOpt;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        if( coords[1] == 0.0 || coords[1] == 4.0 ){
            bdryMaskOpt.push_back( i );
            if(std::abs(coords[1] - 2.0) < 1e-6){
                coords[1]+=0.3;
            }
            if(coords[1] == 0.0){
                coords[1] += 0.5;
                coords[2] += 0.2;
            }
            else{
                coords[1] -= 0.5;
                coords[2] += 0.2;
            }
        }
        
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
    }

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );

    auto foldDofsPtr = std::make_shared<FoldDofsFreeLine<DefaultConfigurator>>(plateTopol,plateGeomInitial, plateGeomRef, bdryMaskRef_1, 2.0, true);

    std::vector<int> foldVertices;
    foldDofsPtr -> getFoldVertices(foldVertices);
    for(int i = 0; i < foldVertices.size(); i++){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, foldVertices[i]);
        std::cout<<"("<<coords[0]<<", "<<coords[1]<<")"<<std::endl;
    }

    auto DfoldDofsPtr = std::make_shared<FoldDofsFreeLineGradient<DefaultConfigurator>>(plateTopol, bdryMaskRef_1, plateGeomInitial, foldVertices, true);

    VectorType edge_weights = VectorType::Zero(plateTopol.getNumEdges());
    foldDofsPtr->getEdgeWeights(edge_weights);

    for(int i = 0; i < foldVertices.size(); i++){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, foldVertices[i]);
        coords[1] += 0.16;
        coords[2]-= 0.5;
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, foldVertices[i]);
    }

    size_t nFoldDOFs = foldDofsPtr->getNumDofs();
    size_t nVertexDOFs = 3*plateTopol.getNumVertices();

    std::cout<<"Number of fold DOFs: "<<nFoldDOFs<<std::endl;
    
    SQPLineSearchParams<DefaultConfigurator> pars;
    CostFunctional<DefaultConfigurator> costFunctional(foldVertices);
    CostFunctionalGradient<DefaultConfigurator> DcostFunctional(plateTopol,foldVertices);
    
    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;

    VectorType factors(2);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;

    auto factory = std::make_shared<ElasticEnergyFactory<DefaultConfigurator>>(factors, plateTopol, edge_weights);
    BoundaryDOFS<DefaultConfigurator> boundaryDOFs(bdryMaskOpt, nVertexDOFs, nFoldDOFs);
    // Create the degrees of freedom object
    std::vector<RealType> deviations;
    ProblemDOFs<DefaultConfigurator> problemDOFs(VectorType::Zero(nFoldDOFs), plateGeomDef, foldDofsPtr, DfoldDofsPtr);
    SQPLineSearchSolver<DefaultConfigurator> solver(pars, costFunctional, DcostFunctional, factory, boundaryDOFs, problemDOFs, 20);
    solver.solve(plateGeomRef, def_geometries, ref_geometries, fold_DOFs);

    std::string filename;

    for(int i = 0; i < def_geometries.size(); i++){
        filename = "deformed/plate_" + std::to_string(i) + ".ply";
        setGeometry(plate, def_geometries[i]);
        OpenMesh::IO::write_mesh(plate,filename);
        setGeometry(plate, ref_geometries[i]);
        filename = "reference/plate_" + std::to_string(i) + ".ply";
        OpenMesh::IO::write_mesh(plate, filename);
    }

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}