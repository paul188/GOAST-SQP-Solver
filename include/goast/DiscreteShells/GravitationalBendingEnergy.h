//==========================================================================================================
// CLASS FOR THE POTENTIAL ENERGY RESULTING FROM MASS OF VERTICES AND BENDING ENERGY
//==========================================================================================================
/**
 * See SimpleBendingEnergy.h for the bending energy part.
 * See GravitationalEnergy.h for the gravitational energy part.
 * \tparam ConfiguratorType Underlying types for scales, vectors, matrices, etc.
 * \author Johannssen
 */

#include "GravitationalEnergy.h"
#include "SimpleBendingEnergy.h"

template<typename ConfiguratorType>
class GravitationalBendingEnergy
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const MeshTopologySaver& _topology;
  const VectorType&  _inactiveGeometry;
  const bool _activeShellIsDeformed;
  // Stores the forces applied to the vertices
  const VectorType &_mass_distribution; 
  // Factor to scale the potential energy
  const RealType _factor_gravity;
  const RealType _factor_bending;
  // Objects for the simple bending and the gravitational energy
  GravitationalEnergy<ConfiguratorType> &_gravitationalEnergy;
  SimpleBendingEnergy<ConfiguratorType> &_simpleBendingEnergy

public:
  GravitationalBendingEnergy( const MeshTopologySaver& topology,
                                      const VectorType& InactiveGeometry,
		                              const bool ActiveShellIsDeformed,
                                      const VectorType &mass_distribution, 
                                      RealType factor_gravity = 1.,
                                      RealType factor_bending = 1. )
                            : _topology( topology), 
                              _inactiveGeometry(InactiveGeometry), 
                              _activeShellIsDeformed( ActiveShellIsDeformed),
                              _mass_distribution( mass_distribution ),
                              _factor_gravity( factor_gravity ),
                              _factor_bending( factor_bending ){
    
    this->_gravitational_energy = GravitationalEnergy<ConfiguratorType>(_topology,
                                                                _inactiveGeometry,
                                                                _activeShellIsDeformed,
                                                                _mass_distribution,
                                                                _factor_gravity)
    
    this->_simple_bending_energy = SimpleBendingEnergy<ConfiguratorType>(_topology,
                                                                _inactiveGeometry,
		                                                        _activeShellIsDeformed,
                                                                factor_bending)
    
                              }

    //energy evaluation
  void apply( const VectorType& ActiveGeometry, RealType & Dest ) const {

    RealType & Dest_gravity;
    RealType & Dest_bending;


    _gravitational_energy.apply(ActiveGeometry, Dest_gravity);
    _simple_bending_energy.apply(ActiveGeometry, Dest_bending);

    Dest = Dest_gravity + Dest_bending;
  }
};

//==========================================================================================================
//! \brief First derivative of GravitationalBendingEnergy w.r.t. the deformed configuration
//! \author Heeren
template<typename ConfiguratorType>
class GravitationalBendingEnergyGradientDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;
	typedef typename ConfiguratorType::VectorType VectorType;
	typedef typename ConfiguratorType::VecType    VecType;
	typedef typename ConfiguratorType::MatType    MatType;

	const VectorType&  _undefShell;
	const VectorType& _weight;
  const MeshTopologySaver& _topology;
  const VectorType&  _inactiveGeometry;
  const bool _activeShellIsDeformed;
  // Stores the forces applied to the vertices
  const VectorType &_mass_distribution; 
  // Factor to scale the potential energy
  const RealType _factor;

public:
	GravitationalEnergyGradientDef(const MeshTopologySaver& topology,
                           const VectorType& InactiveGeometry,
		                       const bool ActiveShellIsDeformed,
                           const VectorType &mass_distribution, 
                           RealType factor = 1. )
                           : _topology( topology), 
                           _inactiveGeometry(InactiveGeometry), 
                           _activeShellIsDeformed( ActiveShellIsDeformed),
                           _mass_distribution( mass_distribution ),
                           _factor( factor ) {}

	void apply(const VectorType& defShell, VectorType& Dest) const {

  if( defShell.size() != _undefShell.size() ){
      std::cerr << "size of deformed = " << defShell.size() << " vs. size of undeformed = " << _undefShell.size() << std::endl;
      throw BasicException( "GravitationalEnergyGradientDef::apply(): sizes dont match!");
    }

  if (_mass_distribution.size() != _topology.getNumVertices()){
    std::cerr << "size of masses = " << _mass_distribution.size() << " vs. size of vertices = " << _topology.getNumVertices() << std::endl;
    throw BasicException( "GravitationalEnergyGradientDef::apply(): sizes dont match!");
  }

  if( Dest.size() != _undefShell.size() )
    Dest.resize( _undefShell.size() );

  Dest.setZero();

  const TriMesh& mesh = _topology.getGrid();

  for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
    if (!mesh.is_boundary(*v_it)) {
        int vertexIdx = v_it->idx();
        Dest[2*_topology.getNumVertices() + vertexIdx] = VecType(0.,0.,_factor*_mass_distribution[vertexIdx]);
    }
  }
  }
};