// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>

#include <goast/Core.h>
#include <goast/DiscreteShells.h>

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

/** 
 * \brief Minimization of Simple Bending Energy + Potential Gravitational Energy
 * \author Johannssen
 *
 * We optimize this energy by direct optimization via gradient descent or BFGS.
 * 
 */

int main(int argc, char *argv[])
{

try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "EXAMPLE ON GRAVITATIONAL POTENTIAL ENERGY TOGETHER WITH BENDING" << std::endl;
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
        if( coords[0] == 0.0 || coords[0] == 1.0 )
            bdryMask.push_back( i );
        // deform part of boundary
        if( coords[0] == 0.0 ) {
            coords[0] += 0.2;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        }
        else{
            coords[0] -= 0.2;
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

    SimpleBendingEnergy<DefaultConfigurator> E_bend( plateTopol, plateGeomRef, true );
    SimpleBendingGradientDef<DefaultConfigurator> DE_bend( plateTopol, plateGeomRef );
    SimpleBendingHessianDef<DefaultConfigurator> D2E_bend( plateTopol, plateGeomRef );

    typename DefaultConfigurator::RealType energy;

    NonlinearMembraneEnergy<DefaultConfigurator> E_mem( plateTopol, plateGeomRef, true );
    NonlinearMembraneGradientDef<DefaultConfigurator> DE_mem( plateTopol, plateGeomRef );
    NonlinearMembraneHessianDef<DefaultConfigurator> D2E_mem( plateTopol, plateGeomRef );

    RealType factor_membrane = 1.0;
    RealType factor_bending = 1.0;
    //RealType factor_gravity = 1.0;

    VectorType factors(2);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;
    //factors[2] = factor_gravity;

    AdditionOp<DefaultConfigurator> E_tot( factors, E_mem, E_bend);
    AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_mem, DE_bend);

    // set outer optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setGradientIterations( 1000);
    optPars.setBFGSIterations( 10000 );
    optPars.setNewtonIterations( 1000 );
    optPars.setQuietMode( SHOW_TERMINATION_INFO );
    VectorType initialization = plateGeomDef;

    std::cerr << "Start gradient descent... " << std::endl;
    GradientDescent<DefaultConfigurator> GD( E_tot, DE_tot,optPars);
    GD.setBoundaryMask( bdryMask );
    GD.solve( initialization, plateGeomDef );

    std::cerr << "Start Quasi-Newton... " << std::endl;
    initialization = plateGeomDef;
    QuasiNewtonBFGS<DefaultConfigurator> QNM( E_tot, DE_tot, optPars);
    QNM.setBoundaryMask( bdryMask );
    QNM.solve( initialization, plateGeomDef );

    // saving
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "gravitationBendingSolution.ply");

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}