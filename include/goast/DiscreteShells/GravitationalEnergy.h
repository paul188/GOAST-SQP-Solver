//==========================================================================================================
// CLASS FOR THE POTENTIAL ENERGY RESULTING FROM MASS OF VERTICES
//==========================================================================================================
/**
 * Integral term of the potential energy resulting from a force field applied to the vertices of a mesh.
 * $\int_{\Omega} F(\phi(x)) dx$ where $\phi$ is the deformation of the mesh.
 * Simple example is the gravitational force, here 
 * F(\phi(x)) = f(x) \cdot \langle (\phi(x)-x), e_3 \rangle
 * up to a constant f ~ (vertex-mass) and 
 * $\langle \phi(x) - x \rangle $ is a difference in height along some axis.
 * \tparam ConfiguratorType Underlying types for scales, vectors, matrices, etc.
 * \author Johannssen
 */
template<typename ConfiguratorType>
class GravitationalEnergy
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
  const RealType _factor;

public:
  GravitationalEnergy( const MeshTopologySaver& topology,
                                      const VectorType& InactiveGeometry,
		                                  const bool ActiveShellIsDeformed,
                                      const VectorType &mass_distribution, 
                                      RealType factor = 1. )
                            : _topology( topology), 
                              _inactiveGeometry(InactiveGeometry), 
                              _activeShellIsDeformed( ActiveShellIsDeformed),
                              _mass_distribution( mass_distribution ),
                              _factor( factor ) {}

  //energy evaluation
  void apply( const VectorType& ActiveGeometry, RealType & Dest ) const {
    
    if( ActiveGeometry.size() != _inactiveGeometry.size() ){
      std::cerr << "size of active = " << ActiveGeometry.size() << " vs. size of inactive = " << _inactiveGeometry.size() << std::endl;
      throw BasicException( "GravitationalEnergy::apply(): sizes dont match!");
    }

    if (_mass_distribution.size() != _topology.getNumVertices()){
      std::cerr << "size of masses = " << _mass_distribution.size() << " vs. size of vertices = " << _topology.getNumVertices() << std::endl;
      throw BasicException( "GravitationalEnergy::apply(): sizes dont match!");
    }

    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry : &_inactiveGeometry;
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    
    Dest = 0.;

    for( int vertexIdx = 0; vertexIdx < _topology.getNumVertices(); ++vertexIdx ){
      //Does this also sum over boundary vertices? That should be excluded, since we consider them
      // under fixed boundary conditions. MIGHT NEED TO BE FIXED
      // But should be fine, the boundary vertices should not move

      VecType coords_def, coords_undef;
      getXYZCoord<VectorType, VecType>( *defShellP, coords_def, vertexIdx);
      getXYZCoord<VectorType, VecType>( *undefShellP, coords_undef, vertexIdx);

      Dest += _mass_distribution[vertexIdx] * ( coords_def[3] - coords_undef[3] );
    }

    Dest *= _factor;
  }
};

//==========================================================================================================
//! \brief First derivative of GravitationalEnergy w.r.t. the deformed configuration
//! \author Heeren
template<typename ConfiguratorType>
class GravitationalEnergyGradientDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

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

//Gradient w.r.t. the undeformed configuration is just the negative of the gradient w.r.t. the deformed configuration
// so VecType(0.,0.,-_factor*_mass_distribution[vertexIdx]);

//All second derivatives are zero, since the potential energy is linear in the deformation