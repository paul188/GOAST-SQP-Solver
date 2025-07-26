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
#include <goast/SQP/Utils/ObjectFactory.h>

#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <goast/Smoothers.h>
#include <unordered_set>

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

try{
    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE FOLD SIMULATION" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;
    
// load flat plate [0,1]^2
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/paperCrissCross.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial;
    getGeometry( plate, plateGeomInitial);
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "initPlate_0_0.ply");

    // determine boundary mask and deform part of boundary
    std::vector<int> bdryMask, bdryMaskDirichletDef_1, bdryMaskDirichletDef_y;
    RealType angle = M_PI / 3.0;
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        if(coords[0] < 0.04)
        {
            // First translate to x = 0
            RealType x = coords[0];
            coords[0] -= 0.015625;
            RealType x_new = std::cos(angle) * coords[0];
            RealType z_new = -std::sin(angle) * coords[0];
            coords[0] = x_new + 0.015625;
            coords[2] = z_new;
            bdryMask.push_back(i);
            bdryMaskDirichletDef_1.push_back(i);
            coords[0] += 0.1;
        }
        if(coords[0] > 0.96)
        {
            // First, translate x to 0
            RealType x = coords[0];
            coords[0] -= (1.0 - 0.015625);
            RealType x_new = std::cos(-angle) * coords[0];
            RealType z_new = -std::sin(-angle) * coords[0];
            coords[0] = x_new + (1.0 - 0.015625);
            coords[2] = z_new;
            bdryMask.push_back(i);
            bdryMaskDirichletDef_1.push_back(i);
            coords[0] -= 0.1;
        }
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);

        if(std::abs(coords[1]) < 1e-4 || std::abs(coords[1] - 1.0) < 1e-4)
        {
            bdryMaskDirichletDef_y.push_back(i);
        }
    }

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "initPlate_0_1.ply");

    //Next, set the bending edge weights for the edge at x = 0.5
    VectorType edge_weights = VectorType::Ones(plateTopol.getNumEdges());

    for(int edgeIdx = 0; edgeIdx < plateTopol.getNumEdges(); edgeIdx++){
        int node_i = plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
        int node_j = plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

        VecType coords_i, coords_j;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords_i, node_i);
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords_j, node_j);

        if(coords_i[0] == 0.5 && coords_j[0] == 0.5){
            std::cout<<"Edge index: "<<edgeIdx<<std::endl;
            edge_weights[edgeIdx] = 0.0;
        }
    }

    std::vector<int> activeRef_y = (std::vector<int>){0,1,1};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskDirichletDef_y , activeRef_y);
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskDirichletDef_1);

    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntries(bdryMaskDirichletDef_1.begin(), bdryMaskDirichletDef_1.end());
    uniqueEntries.insert(bdryMaskDirichletDef_y.begin(), bdryMaskDirichletDef_y.end());

    std::vector<int> bdryMaskDirichletDef;

    // Move the unique elements back to bdryMaskRef_1 in order
    bdryMaskDirichletDef.assign(uniqueEntries.begin(), uniqueEntries.end());
    std::sort(bdryMaskDirichletDef.begin(), bdryMaskDirichletDef.end());

    DirichletSmoother<DefaultConfigurator> smoother(plateGeomInitial, bdryMaskDirichletDef, plateTopol);
    smoother.apply(plateGeomDef, plateGeomDef);

    // save initialization
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "initPlate_1.ply");

    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords; 
        getXYZCoord<VectorType, VecType>(plateGeomDef, coords, i);
        if(coords[0] < 0.13)
        {
            continue;
        }
        if(coords[0] > 0.86)
        {
            continue;
        }
        coords[2] += (5.0*sqrt(3.0))/4.0*(coords[0]-0.5)*(coords[0]-0.5) - 0.2*sqrt(3.0);
        setXYZCoord<VectorType, VecType>(plateGeomDef, coords, i);
    }
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "initPlate_2.ply");

    int numBdryNodes = bdryMask.size();
    std::cerr << "num of bdry nodes = " << numBdryNodes << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMask );

    std::cerr << "\n\nb) OPTIMIZATION BY DIRECT MINIMIZATION" << std::endl;

    SimpleBendingEnergy<DefaultConfigurator> E_bend( plateTopol, plateGeomRef, true , edge_weights);
    SimpleBendingGradientDef<DefaultConfigurator> DE_bend( plateTopol, plateGeomRef , edge_weights);
    SimpleBendingHessianDef<DefaultConfigurator> D2E_bend( plateTopol, plateGeomRef , edge_weights);

    typename DefaultConfigurator::RealType energy;

    NonlinearMembraneEnergy<DefaultConfigurator> E_mem( plateTopol, plateGeomRef, true );
    NonlinearMembraneGradientDef<DefaultConfigurator> DE_mem( plateTopol, plateGeomRef );
    NonlinearMembraneHessianDef<DefaultConfigurator> D2E_mem( plateTopol, plateGeomRef );

    const VectorType& mass_distribution = VectorType::Ones( plateTopol.getNumVertices() );

    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;

    VectorType factors(2);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;

    AdditionOp<DefaultConfigurator> E_tot( factors, E_mem, E_bend);
    AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_mem, DE_bend);
    AdditionHessian<DefaultConfigurator> D2E_tot(factors, D2E_mem, D2E_bend);

    // set outer optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setGradientIterations( 1000);
    optPars.setBFGSIterations( 1000 );
    optPars.setNewtonIterations( 1000 );
    optPars.setQuietMode( SHOW_TERMINATION_INFO );
    VectorType initialization = plateGeomDef;

    std::cerr<< "Start Newton " <<std::endl;
    initialization = plateGeomDef;
    NewtonMethod<DefaultConfigurator> N( DE_tot, D2E_tot, optPars);
    N.setBoundaryMask( bdryMask );
    N.solve( initialization, plateGeomDef );

    // saving
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "bendingFoldSol_withNewton2.ply");

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}