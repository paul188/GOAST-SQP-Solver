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

RealType norm_l1(const VectorType &constr) {
            // avoid division by zero
            RealType c_l1 = std::numeric_limits<RealType>::epsilon();

            // l <= c(x) <= u
            c_l1 += (- constr).cwiseMax(0.0).sum();
            c_l1 += constr.cwiseMax(0.0).sum();

            return c_l1;
}

int main(int argc, char *argv[])
{

  using MatrixType = DefaultConfigurator::SparseMatrixType;
  using FullMatrixType = DefaultConfigurator::FullMatrixType;

try{
    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE FOLD OPTIMIZATION" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

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
            //coords[1]+=0.16;
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
            if(coords[1] == 0.0){
                coords[1] += 0.5;
            }
            else{
                coords[1] -= 0.5;
            }
        }

        if(coords[1] == 2.0){
            coords[1] += 0.16;
            coords[2]-= 0.5;
        }
        
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
    }

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );
    
    auto foldDofsPtr = std::make_shared<FoldDofsSimpleLine<DefaultConfigurator>>(plateTopol,plateGeomInitial, plateGeomRef, bdryMaskRef_1);

    std::vector<int> foldVertices;
    foldDofsPtr -> getFoldVertices(foldVertices);

    for(int i = 0; i < foldVertices.size(); i++){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, foldVertices[i]);
        coords[2]-=0.5;
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, foldVertices[i]);
    }

    auto DfoldDofsPtr = std::make_shared<FoldDofsSimpleLineGradient<DefaultConfigurator>>(plateTopol, bdryMaskRef_1, plateGeomInitial, foldVertices);

    VectorType edge_weights = VectorType::Zero(plateTopol.getNumEdges());
    foldDofsPtr->getEdgeWeights(edge_weights);

    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;

    VectorType factors(2);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;

    SimpleBendingEnergy<DefaultConfigurator> E_bend( plateTopol, plateGeomRef, true , edge_weights);
    SimpleBendingGradientDef<DefaultConfigurator> DE_bend( plateTopol, plateGeomRef , edge_weights);
    SimpleBendingHessianDef<DefaultConfigurator> D2E_bend( plateTopol, plateGeomRef , edge_weights);

    typename DefaultConfigurator::RealType energy;

    NonlinearMembraneEnergy<DefaultConfigurator> E_mem( plateTopol, plateGeomRef, true );
    NonlinearMembraneGradientDef<DefaultConfigurator> DE_mem( plateTopol, plateGeomRef );
    NonlinearMembraneHessianDef<DefaultConfigurator> D2E_mem( plateTopol, plateGeomRef );

    AdditionOp<DefaultConfigurator> E_tot( factors, E_mem, E_bend);
    AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_mem, DE_bend);
    AdditionHessian<DefaultConfigurator> D2E_tot(factors, D2E_mem, D2E_bend);

    // set outer optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setGradientIterations( 1000);
    optPars.setBFGSIterations( 1000 );
    optPars.setNewtonIterations( 10000 );
    optPars.setQuietMode( SHOW_TERMINATION_INFO );
    VectorType initialization = plateGeomDef;

    std::cerr<< "Start Newton " <<std::endl;
    initialization = plateGeomDef;
    NewtonMethod<DefaultConfigurator> N( DE_tot, D2E_tot, optPars);
    N.setBoundaryMask( bdryMaskOpt );
    N.solve( initialization, plateGeomDef );

    // saving
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "bendingFoldSol_withNewton2.ply");

    CostFunctional<DefaultConfigurator> costFunctional(foldVertices);
    RealType cost = 0.0;
    costFunctional.apply(plateGeomDef, cost);
    std::cout<<"Cost Functional value: "<<std::setprecision(12)<<cost<<std::endl;
    MyObjectFactory<DefaultConfigurator> factory(factors, plateTopol, edge_weights);
    BoundaryDOFS<DefaultConfigurator> boundaryDOFs(bdryMaskOpt, 3*plateTopol.getNumVertices(), foldVertices.size());
    ProblemDOFs<DefaultConfigurator> problemDOFs(VectorType::Zero(foldVertices.size()), plateGeomDef, foldDofsPtr, DfoldDofsPtr);
    VectorType DE_val;
    factory.produceDE_vertex(problemDOFs,DE_val);
    boundaryDOFs.transformToReducedSpace(DE_val);

    std::cout<<std::setprecision(12)<<norm_l1(DE_val)<<std::endl;

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}