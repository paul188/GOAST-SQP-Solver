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
    std::cerr << "EXAMPLE ON GRAVITATIONAL POTENTIAL ENERGY" << std::endl;
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
        // deform part of boundary
        if( coords[0] > 0.95 ) {
            coords[2] += 0.2;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        }
    }
    int numBdryNodes = bdryMask.size();
    std::cerr << "num of bdry nodes = " << numBdryNodes << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMask );

    // save initialization
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "initPlate.ply");

    // NOW CONSIDER SOLUTION BY OPTIMIZING DIRICHLET ENERGY
    std::cerr << "\n\nb) OPTIMIZATION BY DIRECT MINIMIZATION" << std::endl;

    //THIS IS MY CODE
    VectorType uniform_mass_dist = VectorType::Ones( plateTopol.getNumVertices() );

    // set up energies
    GravitationalEnergy<DefaultConfigurator> E( plateTopol, plateGeomDef, true, uniform_mass_dist, 1.0 );


    //OLD CODE 

    // set up energies
    DirichletEnergy<DefaultConfigurator> E( plateTopol, plateGeomRef, true );
    DirichletGradientDef<DefaultConfigurator> DE( plateTopol, plateGeomRef );
    DirichletHessianDef<DefaultConfigurator> D2E( plateTopol, plateGeomRef );
    typename DefaultConfigurator::RealType energy;

    // set outer optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setGradientIterations( 100 );
    optPars.setBFGSIterations( 50 );
    optPars.setNewtonIterations( 10 );
    optPars.setQuietMode( SHOW_TERMINATION_INFO );
    
    // optmization with gradient descent
    std::cerr << "Start gradient descent... " << std::endl;
    VectorType initialization = plateGeomDef;
    GradientDescent<DefaultConfigurator> GD( E, DE, optPars );
    GD.setBoundaryMask( bdryMask );
    GD.solve( initialization, plateGeomDef );

    // optmization with BFGS
    std::cerr << "Start Quasi-Newton... " << std::endl;
    initialization = plateGeomDef;
    QuasiNewtonBFGS<DefaultConfigurator> QNM( E, DE, optPars );
    QNM.setBoundaryMask( bdryMask );
    QNM.solve( initialization, plateGeomDef );

    // saving
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "solutionDirichlet_DirectMinimization.ply");

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}