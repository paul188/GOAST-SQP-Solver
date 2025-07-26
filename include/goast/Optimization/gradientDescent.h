// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for gradient descent methods
 * \author Heeren
 *
 * \todo Convert to new interface
 * \todo Documentation!
 */

#ifndef OPTIMIZATION_GRADIENTDESCENT_H
#define OPTIMIZATION_GRADIENTDESCENT_H

#include "optInterface.h"
#include "optParameters.h"
#include "stepsizeControl.h"

/**
 * \brief Simple gradient descent method
 * \tparam ConfiguratorType Container with data types
 * \author Heeren
 */
template<typename ConfiguratorType>
class GradientDescent {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;

  const BaseOp<VectorType, RealType> &_E;
  const BaseOp<VectorType, VectorType> &_DE;

  TIMESTEP_CONTROLLER _timestepController;
  StepsizeControl<ConfiguratorType> _stepsizeControl;

  int _maxIterations;
  RealType _stopEpsilon;
  bool _useNonlinearCG, _useGradientBasedStopping;
  QUIET_MODE _quietMode;
  const std::vector<int> *_bdryMask;

  const BaseOp<VectorType, VectorType> *_Preconditioner;

public:
  GradientDescent( const BaseOp<VectorType, RealType> &E,
                   const BaseOp<VectorType, VectorType> &DE,
                   int MaxIterations,
                   RealType StopEpsilon,
                   TIMESTEP_CONTROLLER TimestepController,
                   bool quiet,
                   RealType sigma = 0.1,
                   RealType tauMin = 1.e-6,
                   RealType tauMax = 4. )
          : _E( E ), _DE( DE ), _timestepController( static_cast<TIMESTEP_CONTROLLER>(TimestepController)),
            _stepsizeControl( _E, _DE, _timestepController, sigma, 0.9, 1., tauMin, tauMax ),
            _maxIterations( MaxIterations ), _stopEpsilon( StopEpsilon ), _useNonlinearCG( false ),
            _useGradientBasedStopping( true ), _quietMode( quiet ? SUPERQUIET : SHOW_ALL ), _bdryMask( nullptr ), _Preconditioner( nullptr ) {}
            
  GradientDescent( const BaseOp<VectorType, RealType> &E,
                   const BaseOp<VectorType, VectorType> &DE,
                   int MaxIterations = 1000,
                   RealType StopEpsilon = 1e-8,
                   TIMESTEP_CONTROLLER TimestepController = ARMIJO,
                   QUIET_MODE quietMode = SUPERQUIET,
                   RealType sigma = 0.1,
                   RealType tauMin = 1.e-6,
                   RealType tauMax = 4. )
          : _E( E ), _DE( DE ), _timestepController( static_cast<TIMESTEP_CONTROLLER>(TimestepController)),
            _stepsizeControl( _E, _DE, _timestepController, sigma, 0.9, 1., tauMin, tauMax ),
            _maxIterations( MaxIterations ), _stopEpsilon( StopEpsilon ), _useNonlinearCG( false ),
            _useGradientBasedStopping( true ), _quietMode( quietMode ), _bdryMask( nullptr ), _Preconditioner( nullptr ) {}

  GradientDescent( const BaseOp<VectorType, RealType> &E,
                   const BaseOp<VectorType, VectorType> &DE,
                   const OptimizationParameters<ConfiguratorType> &optPars )
          : _E( E ), _DE( DE ), _timestepController( static_cast<TIMESTEP_CONTROLLER>(optPars.getGradTimeStepping())),
            _stepsizeControl( _E, _DE, _timestepController, optPars.getSigma(), optPars.getBeta(),
                              optPars.getStartTau(), optPars.getTauMin(), optPars.getTauMax()),
            _maxIterations( optPars.getGradientIterations()), _stopEpsilon( optPars.getStopEpsilon()),
            _useNonlinearCG( false ), _useGradientBasedStopping( true ), _quietMode( optPars.getQuietMode()), _bdryMask( nullptr ),
            _Preconditioner( nullptr ) 
            {}

  void setGradientBasedStopping( bool useGradientBasedStopping = true ) {
    _useGradientBasedStopping = useGradientBasedStopping;
  }

  void setBoundaryMask( const std::vector<int> &Mask ) {
    _bdryMask = &Mask;
    _stepsizeControl.setBoundaryMask( Mask );
  }


  void setNonlinearCG( const bool value ) {
    _useNonlinearCG = value;
  }
  void setPreconditioner( const BaseOp<VectorType, VectorType> &Preconditioner ) {
    _Preconditioner = &Preconditioner;
  }

  // x^{k+1} = x^k - tau * E'[x^k]
  void solve( const VectorType &x_0, VectorType &x_k ) const {

    std::string _path;
    std::string _path_vert_1;
    std::string _path_vert_2;
    std::string _path_edge_vec_1;
    std::string _path_edge_vec_2;
    std::string _path_normal_1;
    std::string _path_normal_2;
    std::string _path_normal_3;
    std::string _path_ruling_0;
    std::string _path_ruling_1;
    std::string _path_ruling_2;
    std::string _path_dev;

    std::ofstream _outFile;
    std::ofstream _outFile_vert_1;
    std::ofstream _outFile_vert_2;
    std::ofstream _outFile_edge_vec_1;
    std::ofstream _outFile_edge_vec_2;
    std::ofstream _outFile_normal_1;
    std::ofstream _outFile_normal_2;
    std::ofstream _outFile_normal_3;
    std::ofstream _outFile_ruling_0;
    std::ofstream _outFile_ruling_1;
    std::ofstream _outFile_ruling_2;
    std::ofstream _outFile_dev;

    _path = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/norm.txt";
    _path_vert_1 = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/vert_1.txt";
    _path_vert_2 = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/vert_2.txt";
    _path_edge_vec_1 = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/edge_vec_1.txt";
    _path_edge_vec_2 = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/edge_vec_2.txt";
    _path_normal_1 = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/normal_1.txt";
    _path_normal_2 = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/normal_2.txt";
    _path_normal_3 = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/normal_3.txt";
    _path_ruling_0 = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/ruling_0.txt";
    _path_ruling_1 = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/ruling_1.txt";
    _path_ruling_2 = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/ruling_2.txt";
    _path_dev = "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/total_norm/dev.txt";
    
    _outFile = std::ofstream(_path, std::ios::app);  // open in append mode
    if (!_outFile) {
        std::cerr << "Failed to open file.\n";
    }

    _outFile_vert_1 = std::ofstream(_path_vert_1, std::ios::app);  // open in append mode
    if (!_outFile_vert_1) {
        std::cerr << "Failed to open file for vert_1.\n";
    }

    _outFile_vert_2 = std::ofstream(_path_vert_2, std::ios::app);  // open in append mode
    if (!_outFile_vert_2) {
        std::cerr << "Failed to open file for vert_2.\n";
    }

    _outFile_edge_vec_1 = std::ofstream(_path_edge_vec_1, std::ios::app);  // open in append mode
    if (!_outFile_edge_vec_1) {
        std::cerr << "Failed to open file for edge_vec_1.\n";
    }

    _outFile_edge_vec_2 = std::ofstream(_path_edge_vec_2, std::ios::app);  // open in append mode
    if (!_outFile_edge_vec_2) {
        std::cerr << "Failed to open file for edge_vec_2.\n";
    }

    _outFile_normal_1 = std::ofstream(_path_normal_1, std::ios::app);  // open in append mode
    if (!_outFile_normal_1) {
        std::cerr << "Failed to open file for normal_1.\n";
    }

    _outFile_normal_2 = std::ofstream(_path_normal_2, std::ios::app);  // open in append mode
    if (!_outFile_normal_2) {
        std::cerr << "Failed to open file for normal_2.\n";
    }

    _outFile_normal_3 = std::ofstream(_path_normal_3, std::ios::app);  // open in append mode
    if (!_outFile_normal_3) {
        std::cerr << "Failed to open file for normal_3.\n";
    }

    _outFile_ruling_0 = std::ofstream(_path_ruling_0, std::ios::app);  // open in append mode
    if (!_outFile_ruling_0) {
        std::cerr << "Failed to open file for ruling_0.\n";
    }

    _outFile_ruling_1 = std::ofstream(_path_ruling_1, std::ios::app);  // open in append mode
    if (!_outFile_ruling_1) {
        std::cerr << "Failed to open file for ruling_1.\n";
    }

    _outFile_ruling_2 = std::ofstream(_path_ruling_2, std::ios::app);  // open in append mode
    if (!_outFile_ruling_2) {
        std::cerr << "Failed to open file for ruling_2.\n";
    }

    _outFile_dev = std::ofstream(_path_dev, std::ios::app);  // open in append mode
    if (!_outFile_dev) {
        std::cerr << "Failed to open file for dev.\n";
    }


    if ( _maxIterations == 0 )
      return;

    x_k = x_0;
    RealType tau = _stepsizeControl.getStartTau();

    // compute initla energy
    RealType energyScalar;
    _E.apply( x_k, energyScalar );                             // just for the stopping criterion
    RealType energy = energyScalar;
    RealType energyNew = energyScalar;

    if ( _quietMode == SHOW_ALL ) {
      std::cout << "=========================================================================================" << std::endl;
      std::cout << "Start gradient descent with " << _maxIterations << " iterations and eps = " << _stopEpsilon << "."
                << std::endl;
      std::cout << "Initial energy: " << energy << std::endl;
      std::cout << "=========================================================================================" << std::endl;
    }

    int iterations = 0;
    bool stoppingCriterion;
    RealType error;

    // store old gradient and direction for nonlinear CG
    VectorType oldDirection( x_k.size()), oldGradient( x_k.size()), currGradient;
    _DE.apply( x_k, currGradient );
    if ( _bdryMask )
      applyMaskToVector( *_bdryMask, currGradient );

    // iteration loop
    do {
      auto t_start = std::chrono::high_resolution_clock::now();
      iterations++;
      // Save the energy at the current position.
      energy = energyNew;

      // compute descent diraction (ie. negative gradient)
      VectorType descentDir( currGradient );
      descentDir *= -1.;

      // apply preconditioner to descent direction if given
      if ( _Preconditioner )
        _Preconditioner->apply( descentDir, descentDir );


      // Fletcher-Reeves nonlinear conjugate gradient method.
      if ( _useNonlinearCG ) {
        const RealType nonlinCGBeta = ((iterations - 1) % 10) ? (descentDir.squaredNorm() / oldGradient.squaredNorm())
                                                              : 0.;
        oldGradient = descentDir;
        if ( nonlinCGBeta > 0. )
          descentDir += nonlinCGBeta * oldDirection;
        oldDirection = descentDir;
      }

      // compute tau
      tau = _stepsizeControl.getStepsize( x_k, currGradient, descentDir, tau, energy );

      // If the nonlinear CG direction is not a descent direction, use the normal descent direction instead.
      if ( _useNonlinearCG && (tau == 0)) {
        descentDir = oldGradient;
        oldDirection.setZero();
        tau = _stepsizeControl.getStepsize( x_k, currGradient, descentDir, tau, energyNew );
      }

      // update position
      x_k += tau * descentDir;

      // Calculate energy and gradient at the new position
      _DE.apply( x_k, currGradient );
      if ( _bdryMask )
        applyMaskToVector( *_bdryMask, currGradient );
      
      _outFile << "here"<<std::endl;
      _outFile_vert_1 <<" here"<<std::endl;
      _outFile_vert_2 <<" here"<<std::endl;
      _outFile_edge_vec_1 <<" here"<<std::endl;
      _outFile_edge_vec_2 <<" here"<<std::endl;
      _outFile_normal_1 <<" here"<<std::endl;
      _outFile_normal_2 <<" here"<<std::endl;
      _outFile_normal_3 <<" here"<<std::endl;
      _outFile_ruling_0 <<" here"<<std::endl;
      _outFile_ruling_1 <<" here"<<std::endl;
      _outFile_ruling_2 <<" here"<<std::endl;
      _outFile_dev <<" here"<<std::endl;
      _E.apply( x_k, energyScalar );
      energyNew = energyScalar;

      // compute error
      error = _useGradientBasedStopping ? currGradient.norm() : energy - energyNew;

      auto t_end = std::chrono::high_resolution_clock::now();

      if ( _quietMode == SHOW_ALL )
        std::cout << "step = " << iterations
                  << std::scientific << " , stepsize = " << tau << ", energy = " << energyNew << ", error = " << error
                  << std::fixed << ", time = " << std::chrono::duration<double, std::milli>( t_end - t_start ).count()
                  << "ms" << std::endl;

      // stoppingCriterion == true means the iteration will be continued
      stoppingCriterion = (error > _stopEpsilon) && (iterations < _maxIterations) && (tau > 0);

    } while ( stoppingCriterion ); // end iteration loop

    if ( _quietMode != SUPERQUIET ) {
      std::cout << "=========================================================================================" << std::endl;
      std::cout << "Finished gradient descent after " << iterations << " steps (max. steps = " << _maxIterations << ", tol = " << _stopEpsilon << ")." << std::endl;
      std::cout << "Final stepsize = " << std::scientific << std::setprecision(15) << tau << ", energy = " << energyNew << ", error = " << error << std::endl;
      std::cout << "=========================================================================================" << std::endl;
    }

  }


};


#endif //GRADIENTDESCENT_H
