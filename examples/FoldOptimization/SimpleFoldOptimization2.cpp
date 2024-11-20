/* 
 Simulating a fold that results from bending a plate with clamped 
 boundary conditions. Only NonlinearMembraneEnergy and SimpleBendingEnergy are used.
 The gravitational energy is not used in this example.
 We set the weight of the contribution of the middle edge to the bending energy to 0.
 We thus expect a lot of the bending to happen in the middle of the plate, at the fold.
 Imagine a piece of paper that is folded in the middle in both directions.
 So that the fold does not point into any direction, but the edge just doesnt
 resist bending.
*/

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <iostream>
#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include "Arena.h"
#include <goast/Smoothers.h>
#include "CostFunctional.h"
#include "FoldDofs.h"
#include "Merit.h"
#include "SQP.h"
#include "SparseMat.h"
//#include "BoundaryDOFS.h"
#include <chrono>
#include <thread>
#include "ObjectFactory.h"

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

/** 
 * \brief Simulation of the optimization of a simple fold
 * \author Johannssen
 *
 * The cost functional is defined in CostFunctional.h
 * Arena.h is used to manage memory allocation
 * 
 */

/**/

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
    OpenMesh::IO::read_mesh(plate, "../../data/plate/plate4SD.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial;
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );
    getGeometry( plate, plateGeomInitial );

    Arena arena;

    // determine boundary mask for reference geometry and foldVertices
    // the part of the reference boundary where everything is fixed
    std::vector<int> bdryMaskRef_1;
    // the part of the reference boundary where only y,z coordinates are fixed -> x free
    std::vector<int> bdryMaskRef_2;

     // determine boundary mask for optimization
     // and deform part of boundary
    std::vector<int> bdryMaskOpt;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        if( coords[0] == 0.0 || coords[0] == 1.0 ){
            bdryMaskOpt.push_back( i );
        // deform part of boundary
        }
        
        coords[0] *= 0.8;
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
    }

    VectorType edge_weights_2 = VectorType::Ones(plateTopol.getNumEdges());

    for(int edgeIdx = 0; edgeIdx < plateTopol.getNumEdges(); edgeIdx++){
        int node_i = plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
        int node_j = plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

        VecType coords_i, coords_j;
        getXYZCoord<VectorType, VecType>( plateGeomInitial, coords_i, node_i);
        getXYZCoord<VectorType, VecType>( plateGeomInitial, coords_j, node_j);

        if(coords_i[0] == 0.5 && coords_j[0] == 0.5){
            edge_weights_2[edgeIdx] = 0;
        }
    }

    // extend all boundary masks to (x,y,z) coordinates

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef_1 );
    std::vector<int> activeRef_2 = (std::vector<int>){0,1,1};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_2 , activeRef_2);
    // append to the Dirichlet bdry mask of Reference boundary 1
    bdryMaskRef_1.insert(bdryMaskRef_1.end(), bdryMaskRef_2.begin(), bdryMaskRef_2.end());
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );

    // prepare the initial deformed geometry
    //DirichletSmoother<DefaultConfigurator> smoother(plateGeomRef, plateGeomDef, bdryMaskRef_1, plateTopol);
    //smoother.smoothMesh(plateGeomDef);

    OpenMesh::IO::write_mesh(plate, "testPlate1.ply");

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "initialDefPlate.ply");

    auto foldDofsPtr = std::make_unique<FoldDofsSimpleLine<DefaultConfigurator>>(plateTopol,plateGeomRef);

    std::vector<int> foldVertices;
    foldDofsPtr -> getFoldVertices(foldVertices);
    SQPLineSearchParams<DefaultConfigurator> pars;
    CostFunctional<DefaultConfigurator> costFunctional(foldVertices);
    CostFunctionalGradient<DefaultConfigurator> DcostFunctional(plateTopol,foldVertices);
    MyObjectFactory<DefaultConfigurator> factory;
    auto DFoldDofsPtr = std::make_unique<FoldDofsSimpleLineGradient<DefaultConfigurator>>(plateTopol,foldVertices);
    SQPLineSearchSolver<DefaultConfigurator> solver(pars, plateTopol, costFunctional, DcostFunctional, factory, bdryMaskOpt, bdryMaskRef_1, std::move(foldDofsPtr), std::move(DFoldDofsPtr));

    solver.solve(plateGeomDef, plateGeomRef, plateGeomInitial);

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "result_plateGeomDef.ply");

    setGeometry( plate, plateGeomRef );
    OpenMesh::IO::write_mesh(plate, "result_plateGeomRef.ply");
  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}