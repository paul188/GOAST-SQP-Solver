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

#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <cmath>

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

int main(int argc, char *argv[])
{

try{
    
    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "EXAMPLE ON DIRICHLET ENERGY" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    // load flat plate [0,1]^2
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/plate4SD.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef;
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );

    // determine boundary mask and deform part of boundary
    std::vector<int> bdryMask;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        if( coords[0] < 0.05 || coords[0] > 0.95 )
            bdryMask.push_back( i );
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        if( coords[0]  > 0.95 ) {
            coords[0] = 1.0 - coords[0];
            RealType phi = 2*M_1_PI * 0.5;
            RealType coords0 = coords[0];
            RealType coords1 = coords[1];
            coords[0] = cos(phi) * coords0 + sin(phi) * coords1;
            coords[1] = - sin(phi) * coords0 + cos(phi) * coords1;
            // coords[2] = -0.4 + coords[0];
            coords[0] += 0.3;

            coords[0] = 1.0 - coords[0];
            
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        }

    }

    int numBdryNodes = bdryMask.size();
    std::cerr << "num of bdry nodes = " << numBdryNodes << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMask );

    // save initialization
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "initPlate.ply");

    // FIRST CONSIDER SOLUTION OF EULER-LAGRANGE-EQUATION
    std::cerr << "\n\na) OPTIMIZATION BY SOLVING EULER-LAGRANGE EQUATION" << std::endl;
    // assemble and mask stiffness matrix
    std::cerr << "Set up system matrix" << std::endl;
    typename DefaultConfigurator::SparseMatrixType StiffnessMatrix;
    computeStiffnessMatrix<DefaultConfigurator>( plateTopol, plateGeomRef, StiffnessMatrix );
    applyMaskToMajor<typename DefaultConfigurator::SparseMatrixType>( bdryMask, StiffnessMatrix );

    // set up right hand side and mask
    VectorType rhs( plateGeomDef );
    rhs.setZero();
    for( int i = 0; i < bdryMask.size(); i++ )
        rhs[bdryMask[i]] = plateGeomDef[bdryMask[i]];

    // set up linear system and solve
    std::cerr << "Set up linear system and solve" << std::endl;
    LinearSolver<DefaultConfigurator> directSolver;
    VectorType solution;
    directSolver.prepareSolver( StiffnessMatrix );
    directSolver.backSubstitute( rhs, solution );

    // get final energy value
    computeStiffnessMatrix<DefaultConfigurator>( plateTopol, plateGeomRef, StiffnessMatrix );
    VectorType temp = StiffnessMatrix * solution;
    std::cerr << "Final Dirichlet energy = " << 0.5 * temp.dot( solution ) << std::endl;

    // saving
    setGeometry( plate, solution );
    OpenMesh::IO::write_mesh(plate, "solutionDirichlet_EulerLangrange.ply");

    // NOW CONSIDER SOLUTION BY OPTIMIZING DIRICHLET ENERGY
    std::cerr << "\n\nb) OPTIMIZATION BY DIRECT MINIMIZATION" << std::endl;

    // set up energies
    SimpleBendingEnergy<DefaultConfigurator> E_bend( plateTopol, plateGeomRef, true );
    SimpleBendingGradientDef<DefaultConfigurator> DE_bend( plateTopol, plateGeomRef );
    SimpleBendingHessianDef<DefaultConfigurator> D2E_bend( plateTopol, plateGeomRef );
    NonlinearMembraneEnergy<DefaultConfigurator> E_mem( plateTopol, plateGeomRef, true );
    NonlinearMembraneGradientDef<DefaultConfigurator> DE_mem( plateTopol, plateGeomRef );
    NonlinearMembraneHessianDef<DefaultConfigurator> D2E_mem( plateTopol, plateGeomRef );
    typename DefaultConfigurator::RealType energy;

    VectorType factor = VectorType::Ones(2);
    factor[1] *= 1.e-3;

    AdditionOp<DefaultConfigurator> E_tot( factor, E_mem, E_bend );
    AdditionGradient<DefaultConfigurator> DE_tot( factor, DE_mem, DE_bend );
    AdditionHessian<DefaultConfigurator> D2E_tot( factor, D2E_mem, D2E_bend );

    // set outer optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setGradientIterations( 10000 );
    optPars.setBFGSIterations( 10000 );
    optPars.setNewtonIterations( 100000 );
    optPars.setQuietMode( SHOW_TERMINATION_INFO );
    VectorType initialization = solution;
    
    // optmization with gradient descent
    std::cerr << "Start gradient descent... " << std::endl;

    GradientDescent<DefaultConfigurator> GD( E_tot, DE_tot, optPars );
    GD.setBoundaryMask( bdryMask );
    GD.solve( initialization, plateGeomDef );

    // optmization with BFGS
    std::cerr << "Start Quasi-Newton... " << std::endl;
    initialization = plateGeomDef;
    QuasiNewtonBFGS<DefaultConfigurator> QNM( E_tot, DE_tot, optPars );
    QNM.setBoundaryMask( bdryMask );
    QNM.solve( initialization, plateGeomDef );

    // optmization with Newton
    std::cerr << "Start Newton... " << std::endl;
    initialization = plateGeomDef;
    NewtonMethod<DefaultConfigurator> NM( DE_tot, D2E_tot, optPars );
    NM.setBoundaryMask( bdryMask );
    NM.solve( initialization, plateGeomDef );
    

    // saving
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "solutionDirichlet_DirectMinimization.ply");

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}