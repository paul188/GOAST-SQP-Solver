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
    std::cerr << "EXAMPLE ON GRAVITATIONAL POTENTIAL ENERGY TOGETHER WITH BENDING AND MEMBRANE" << std::endl;
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
        if( coords[0] == 0.0 || coords[0] == 1.0 ){
            bdryMask.push_back( i );
        // deform part of boundary
        if( coords[0] == 0.0 ) {
            coords[0] += 0.2;
        }
        else if(coords[0] == 1.0){
            coords[0] -= 0.2;
        }
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        }
    }
    int numBdryNodes = bdryMask.size();
    std::cerr << "num of bdry nodes = " << numBdryNodes << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMask );

    // save initialization
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "initPlate.ply");

    std::cerr << "\n\nb) OPTIMIZATION BY DIRECT MINIMIZATION" << std::endl;

    SimpleBendingEnergy<DefaultConfigurator> E_bend( plateTopol, plateGeomRef, true );
    SimpleBendingGradientDef<DefaultConfigurator> DE_bend( plateTopol, plateGeomRef );
    SimpleBendingHessianDef<DefaultConfigurator> D2E_bend( plateTopol, plateGeomRef );

    NonlinearMembraneEnergy<DefaultConfigurator> E_mem( plateTopol, plateGeomRef, true );
    NonlinearMembraneGradientDef<DefaultConfigurator> DE_mem( plateTopol, plateGeomRef );
    NonlinearMembraneHessianDef<DefaultConfigurator> D2E_mem( plateTopol, plateGeomRef );

    const VectorType& mass_distribution = VectorType::Ones( plateTopol.getNumVertices() );

    VectorType gravity_dir;
    gravity_dir.resize(3);
    //In which direction does the gravitational force act?
    gravity_dir << 0.0, 0.0, -1.0;

    GravitationalEnergy<DefaultConfigurator> E_grav(plateTopol, plateGeomRef, true, mass_distribution, gravity_dir);
    GravitationalEnergyGradientDef<DefaultConfigurator> DE_grav(plateTopol, plateGeomRef, mass_distribution, gravity_dir);

    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;
    RealType factor_gravity = 1.0;

    VectorType factors(3);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;
    factors[2] = factor_gravity;

    AdditionOp<DefaultConfigurator> E_tot( factors, E_mem, E_bend, E_grav);
    AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_mem, DE_bend, DE_grav);
    AdditionHessian<DefaultConfigurator> D2E_tot( factors, D2E_mem, D2E_bend);

    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setNewtonIterations( 1000 );
    optPars.setQuietMode( SHOW_TERMINATION_INFO );
    VectorType initialization = plateGeomDef;

    std::cerr<< "Startting Newton Linesearch..."<<std::endl;
    LineSearchNewton<DefaultConfigurator> NLS( E_tot, DE_tot, D2E_tot, optPars);
    NLS.setBoundaryMask( bdryMask );
    NLS.solve( initialization, plateGeomDef );

    // saving
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "gravitationBendingSolution.ply");

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}