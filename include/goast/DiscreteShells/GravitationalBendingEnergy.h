//==========================================================================================================
// CLASS FOR THE POTENTIAL ENERGY RESULTING FROM MASS OF VERTICES AND BENDING ENERGY
//==========================================================================================================
/**
 * See SimpleBendingEnergy.h for the bending energy part.
 * See GravitationalEnergy.h for the gravitational energy part.
 * \tparam ConfiguratorType Underlying types for scales, vectors, matrices, etc.
 * \author Johannssen
 */
#ifndef  GRAVITATIONALBENDINGENERGY_HH
#define  GRAVITATIONALBENDINGENERGY_HH

#include "GravitationalEnergy.h"
#include "SimpleBendingEnergy.h"

template<typename ConfiguratorType>
class GravitationalBendingEnergy
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  RealType _factor_bending;
  RealType _factor_gravity;

  // Objects for the simple bending and the gravitational energy
  GravitationalEnergy<ConfiguratorType>& _gravitationalEnergy;
  SimpleBendingEnergy<ConfiguratorType>& _simpleBendingEnergy;

public:
  GravitationalBendingEnergy( GravitationalEnergy<ConfiguratorType>& gravitationalEnergy,
                              SimpleBendingEnergy<ConfiguratorType>& simpleBendingEnergy,
                              RealType factor_gravity = 1.,
                              RealType factor_bending = 1.)
                              : _gravitationalEnergy(gravitationalEnergy),
                              _simpleBendingEnergy(simpleBendingEnergy),
                              _factor_bending(factor_bending),
                              _factor_gravity(factor_gravity){}

    //energy evaluation
  void apply( const VectorType& ActiveGeometry, RealType& Dest ) const override {

    RealType Dest_gravity = 0;
    RealType Dest_bending = 0;

    _gravitationalEnergy.apply(ActiveGeometry, Dest_gravity);
    _simpleBendingEnergy.apply(ActiveGeometry, Dest_bending);

    Dest = _factor_gravity*Dest_gravity + _factor_bending*Dest_bending;
  }
};

//==========================================================================================================
//! \brief First derivative of GravitationalBendingEnergy w.r.t. the deformed configuration
//! \author Johannssen
template<typename ConfiguratorType>
class GravitationalBendingEnergyGradientDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;
	typedef typename ConfiguratorType::VectorType VectorType;
	typedef typename ConfiguratorType::VecType    VecType;
	typedef typename ConfiguratorType::MatType    MatType;

  const RealType _factor_gravity;
  const RealType _factor_bending;
  // Objects for the gradients of the simple bending and gravitational energy
  GravitationalEnergyGradientDef<ConfiguratorType> &_gravitationalEnergyGradient;
  SimpleBendingGradientDef<ConfiguratorType> &_simpleBendingGradient;

public:
    GravitationalBendingEnergyGradientDef(
                           GravitationalEnergyGradientDef<ConfiguratorType>& gravitationalEnergyGradient,
                           SimpleBendingGradientDef<ConfiguratorType>& simpleBendingGradient,
                           RealType factor_gravity = 1.,
                           RealType factor_bending = 1.)
                           : _factor_gravity( factor_gravity ),
                           _factor_bending( factor_bending ),
                           _gravitationalEnergyGradient(gravitationalEnergyGradient),
                           _simpleBendingGradient(simpleBendingGradient){}

    void apply(const VectorType& defShell, VectorType& Dest) const {
        VectorType Dest_gravity;
        _gravitationalEnergyGradient.apply(defShell, Dest_gravity);
        VectorType Dest_bending;
        _simpleBendingGradient.apply(defShell, Dest_bending);
        Dest = _factor_gravity*Dest_gravity + _factor_bending*Dest_bending;
    }
};

//==========================================================================================================
//! \brief First derivative of GravitationalBendingEnergy w.r.t. the undeformed configuration
//! \author Johannssen
template<typename ConfiguratorType>
class GravitationalBendingEnergyGradientUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;
	typedef typename ConfiguratorType::VectorType VectorType;
	typedef typename ConfiguratorType::VecType    VecType;
	typedef typename ConfiguratorType::MatType    MatType;

    const RealType _factor_gravity;
    const RealType _factor_bending;
    // Objects for the gradients of the simple bending and gravitational energy
    GravitationalEnergyGradientUndef<ConfiguratorType> &_gravitationalEnergyGradient;
    SimpleBendingGradientUndef<ConfiguratorType> &_simpleBendingGradient;

public:
    GravitationalBendingEnergyGradientUndef(
                           GravitationalEnergyGradientUndef<ConfiguratorType>& gravitationalEnergyGradient,
                           SimpleBendingGradientUndef<ConfiguratorType>& simpleBendingGradient,
                           RealType factor_gravity = 1.,
                           RealType factor_bending = 1.)
                           :_gravitationalEnergyGradient(gravitationalEnergyGradient),
                            _simpleBendingGradient(simpleBendingGradient),
                           _factor_gravity( factor_gravity ),
                           _factor_bending( factor_bending ){}

    void apply(const VectorType& defShell, VectorType& Dest) const {
        VectorType Dest_gravity;
        _gravitationalEnergyGradient.apply(defShell, Dest_gravity);
        VectorType Dest_bending;
        _simpleBendingGradient.apply(defShell, Dest_bending);
        Dest = _factor_gravity*Dest_gravity + _factor_bending*Dest_bending;
    }
};

#endif