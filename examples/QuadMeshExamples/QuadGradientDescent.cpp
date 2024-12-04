#pragma once

#include "optInterface.h"
#include "optParameters.h"
#include "stepsizeControl.h"
#include "QuadMesh.h"

/**
 * \brief Simple gradient descent method adapted for the QuadMesh type
 * \tparam ConfiguratorType Container with data types
 * \author Heeren
 */
template<typename ConfiguratorType>
class QuadGradientDescent {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;

  const BaseOp<VectorType, RealType> &_E1;
  const BaseOp<VectorType, RealType> &_E2;
  const BaseOp<VectorType, VectorType> &_DE1;
  const BaseOp<VectorType, VectorType> &_DE2;

  TIMESTEP_CONTROLLER _timestepController;
  StepsizeControl<ConfiguratorType> _stepsizeControl1;
  StepsizeControl<ConfiguratorType> _stepsizeControl2;

  int _maxIterations;
  RealType _stopEpsilon;
  bool _useNonlinearCG, _useGradientBasedStopping;
  QUIET_MODE _quietMode;
  const std::vector<int> *_bdryMask;

  const BaseOp<VectorType, VectorType> *_Preconditioner;

  public:
  QuadGradientDescent( const BaseOp<VectorType, RealType> &E1,
                   const BaseOp<VectorType, RealType> &E2,
                   const BaseOp<VectorType, VectorType> &DE1,
                   const BaseOp<VectorType, VectorType> &DE2,
                   int MaxIterations,
                   RealType StopEpsilon,
                   TIMESTEP_CONTROLLER TimestepController,
                   bool quiet,
                   RealType sigma = 0.1,
                   RealType tauMin = 1.e-6,
                   RealType tauMax = 4. )
          : _E1( E1 ), _E2( E2 ), _DE1( DE1 ), _DE2( DE2 ), _timestepController( static_cast<TIMESTEP_CONTROLLER>(TimestepController)),
            _stepsizeControl1( _E1, _DE1, _timestepController, sigma, 0.9, 1., tauMin, tauMax ),
            _stepsizeControl2( _E2, _DE2, _timestepController, sigma, 0.9, 1., tauMin, tauMax ),
            _maxIterations( MaxIterations ), _stopEpsilon( StopEpsilon ), _useNonlinearCG( false ),
            _useGradientBasedStopping( true ), _quietMode( quiet ? SUPERQUIET : SHOW_ALL ), _bdryMask( nullptr ), _Preconditioner( nullptr ) {}
            
  QuadGradientDescent( const BaseOp<VectorType, RealType> &E1,
                   const BaseOp<VectorType, RealType> &E2,
                   const BaseOp<VectorType, VectorType> &DE1,
                   const BaseOp<VectorType, VectorType> &DE2,
                   int MaxIterations = 1000,
                   RealType StopEpsilon = 1e-8,
                   TIMESTEP_CONTROLLER TimestepController = ARMIJO,
                   QUIET_MODE quietMode = SUPERQUIET,
                   RealType sigma = 0.1,
                   RealType tauMin = 1.e-6,
                   RealType tauMax = 4. )
          : _E1( E1 ), _E2( E2 ), _DE1( DE1 ), _DE2( DE2 ), _timestepController( static_cast<TIMESTEP_CONTROLLER>(TimestepController)),
            _stepsizeControl1( _E1, _DE1, _timestepController, sigma, 0.9, 1., tauMin, tauMax ),
            _stepsizeControl2( _E2, _DE2, _timestepController, sigma, 0.9, 1., tauMin, tauMax ),
            _maxIterations( MaxIterations ), _stopEpsilon( StopEpsilon ), _useNonlinearCG( false ),
            _useGradientBasedStopping( true ), _quietMode( quietMode ), _bdryMask( nullptr ), _Preconditioner( nullptr ) {}

   QuadGradientDescent( const BaseOp<VectorType, RealType> &E1,
                    const BaseOp<VectorType, RealType> &E2,
                    const BaseOp<VectorType, VectorType> &DE1,
                    const BaseOp<VectorType, VectorType> &DE2,
                    const OptimizationParameters<ConfiguratorType> &optPars )
            : _E1( E1 ), _E2( E2 ), _DE1( DE1 ), _DE2( DE2 ), _timestepController( static_cast<TIMESTEP_CONTROLLER>(optPars.getGradTimeStepping())),
            _stepsizeControl1( _E1, _DE1, _timestepController, optPars.getSigma(), optPars.getBeta(),
                                optPars.getStartTau(), optPars.getTauMin(), optPars.getTauMax()),
            _stepsizeControl2( _E2, _DE2, _timestepController, optPars.getSigma(), optPars.getBeta(),
                                optPars.getStartTau(), optPars.getTauMin(), optPars.getTauMax()),
            _maxIterations( optPars.getGradientIterations()), _stopEpsilon( optPars.getStopEpsilon()),
            _useNonlinearCG( false ), _useGradientBasedStopping( true ), _quietMode( optPars.getQuietMode()), _bdryMask( nullptr ),
            _Preconditioner( nullptr ) {}

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
  void solve( const VectorType &x1_0, const VectorType &x2_0, VectorType &x1_k, VectorType &x2_k ) const {

    if ( _maxIterations == 0 )
      return;

    x1_k = x1_0;
    x2_k = x2_0;
    RealType tau = std::min(_stepsizeControl1.getStartTau(), stepsizeControl2.getStartTau());

    // compute initial energy
    RealType energyScalar1;
    RealType energyScalar2;
    _E1.apply( x1_k, energyScalar1 );
    _E2.apply( x2_k, energyScalar2 );                             // just for the stopping criterion
    RealType energy = 0.5*(energyScalar1 + energyScalar2);
    RealType energyNew = energy;

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
    VectorType oldDirection1( x1_k.size()), oldGradient1( x1_k.size()), currGradient1;
    VectorType oldDirection2( x2_k.size()), oldGradient2( x2_k.size()), currGradient2;

    _DE.apply( x1_k, currGradient1 );
    _DE.apply( x2_k, currGradient2 );

    if ( _bdryMask )
      applyMaskToVector( *_bdryMask, currGradient1 );
      applyMaskToVector( *_bdryMask, currGradient2 );

    // iteration loop
    do {
      auto t_start = std::chrono::high_resolution_clock::now();
      iterations++;
      // Save the energy at the current position.
      energy = energyNew;

      // compute descent diraction (ie. negative gradient)
      VectorType descentDir1( currGradient1 );
      descentDir1 *= -1.;

      VectorType descentDir2( currGradient2 );
      descentDir2 *= -1.;

      // apply preconditioner to descent direction if given
      if ( _Preconditioner )
        _Preconditioner->apply( descentDir1, descentDir1 );
        _Preconditioner->apply( descentDir2, descentDir2 );


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