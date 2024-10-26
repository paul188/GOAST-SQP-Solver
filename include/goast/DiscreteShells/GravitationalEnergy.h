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
#ifndef  GRAVITATIONALENERGY_HH
#define  GRAVITATIONALENERGY_HH

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
  // Which direction does the gravity act in (default is z-direction)
  const VectorType &_gravityDirection;

public:
  GravitationalEnergy( const MeshTopologySaver& topology,
                                      const VectorType& InactiveGeometry,
		                                  const bool ActiveShellIsDeformed,
                                      const VectorType &mass_distribution, 
                                      const VectorType &gravityDirection,
                                      RealType factor = 1.)
                            : _topology( topology), 
                              _inactiveGeometry(InactiveGeometry), 
                              _activeShellIsDeformed( ActiveShellIsDeformed),
                              _mass_distribution( mass_distribution ),
                              _factor( factor ),
                              _gravityDirection( gravityDirection ) {}

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

    if( _gravityDirection.size() != 3 ){
      std::cerr << "size of gravity direction = " << _gravityDirection.size() << std::endl;
      throw BasicException( "GravitationalEnergy::apply(): gravity direction must be 3-dimensional!");
    }

    typename DefaultConfigurator::SparseMatrixType MassMatrix;
    computeMassMatrix<DefaultConfigurator>( _topology, _inactiveGeometry, MassMatrix , false);
    // In contrast to Dirichlet Energy, mass acting on boundary will be taken into account
    // even though the boundary is fixed and does not move
    // -> so not calling applyMaskToMajor()

    // Next, initialize the vector of displacements
    Eigen::MatrixXd Displacements(3,_topology.getNumVertices());

    //How much mass is at each vertex
    MatrixType MassDistributionMatrix = MatrixType(_mass_distribution.asDiagonal());

    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry : &_inactiveGeometry;
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;

    for( int vertexIdx = 0; vertexIdx < _topology.getNumVertices(); ++vertexIdx ){
      VecType coords_def, coords_undef;
      getXYZCoord<VectorType, VecType>( *defShellP, coords_def, vertexIdx);
      getXYZCoord<VectorType, VecType>( *undefShellP, coords_undef, vertexIdx);
      for (int i = 0; i < 3; i++){
        Displacements(i,vertexIdx) = coords_def[i] - coords_undef[i];
      }
    }

    Dest = -_factor*((_gravityDirection.transpose()*Displacements*MassMatrix*MassDistributionMatrix).sum());
  }
};

//==========================================================================================================
//! \brief First derivative of GravitationalEnergy w.r.t. the deformed configuration
//! \author Johannssen
template<typename ConfiguratorType>
class GravitationalEnergyGradientDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;
	typedef typename ConfiguratorType::VectorType VectorType;
	typedef typename ConfiguratorType::VecType    VecType;
	typedef typename ConfiguratorType::MatType    MatType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

	const VectorType&  _undefShell;
  const MeshTopologySaver& _topology;
  // Stores the forces applied to the vertices
  const VectorType &_mass_distribution; 
  // Factor to scale the potential energy
  const RealType _factor;
  // Which direction does the gravity act in (default is z-direction)
  const VectorType &_gravityDirection;

public:
	GravitationalEnergyGradientDef(const MeshTopologySaver& topology,
                           const VectorType& undefShell,
                           const VectorType& mass_distribution, 
                           const VectorType &gravityDirection,
                           RealType factor = 1.)
                           : _topology( topology), 
                            _undefShell(undefShell),
                           _mass_distribution( mass_distribution ),
                           _factor( factor ),
                           _gravityDirection( gravityDirection) {}

	void apply(const VectorType& defShell, VectorType& Dest) const {

  if( defShell.size() != _undefShell.size() ){
      std::cerr << "size of deformed = " << defShell.size() << " vs. size of undeformed = " << _undefShell.size() << std::endl;
      throw BasicException( "GravitationalEnergyGradientDef::apply(): sizes dont match!");
    }

  if (_mass_distribution.size() != _topology.getNumVertices()){
    std::cerr << "size of masses = " << _mass_distribution.size() << " vs. size of vertices = " << _topology.getNumVertices() << std::endl;
    throw BasicException( "GravitationalEnergyGradientDef::apply(): sizes dont match!");
  }

  if( _gravityDirection.size() != 3 ){
      std::cerr << "size of gravity direction = " << _gravityDirection.size() << std::endl;
      throw BasicException( "GravitationalEnergy::apply(): gravity direction must be 3-dimensional!");
  }

  if( Dest.size() != _undefShell.size() )
    Dest.resize( _undefShell.size() );

  Dest.setZero();

  MatrixType MassDistributionMatrix = MatrixType(_mass_distribution.asDiagonal());

  typename DefaultConfigurator::SparseMatrixType MassMatrix;
  computeMassMatrix<DefaultConfigurator>( _topology, _undefShell, MassMatrix , false);
  
  // Now calculate the matrix mass_matrix * stiffness_matrix
  MatrixType AreaMassMatrix = MassDistributionMatrix*MassMatrix;

  for (int i = 0; i < 3; i++){
    for(int vertexIdx = 0; vertexIdx < _topology.getNumVertices(); vertexIdx++){
				Dest[i*_topology.getNumVertices() + vertexIdx] = -_factor*_gravityDirection[i]*AreaMassMatrix.row(vertexIdx).sum();
    }
		}
  }
};

//Gradient w.r.t. the undeformed configuration is just the negative of the gradient w.r.t. the deformed configuration

//==========================================================================================================
//! \brief First derivative of GravitationalBendingEnergy w.r.t. the undeformed configuration
//! \author Johannssen
template<typename ConfiguratorType>
class GravitationalEnergyGradientUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;
	typedef typename ConfiguratorType::VectorType VectorType;
	typedef typename ConfiguratorType::VecType    VecType;
	typedef typename ConfiguratorType::MatType    MatType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;


	const VectorType&  _undefShell;
  const MeshTopologySaver& _topology;
  // Stores the forces applied to the vertices
  const VectorType &_mass_distribution; 
  // Factor to scale the potential energy
  const RealType _factor;
  // Which direction does the gravity act in (default is z-direction)
  const VectorType &_gravityDirection;
public:
	GravitationalEnergyGradientUndef(const MeshTopologySaver& topology,
                           const VectorType& undefShell,
                           const VectorType& mass_distribution, 
                           const VectorType &gravityDirection,
                           RealType factor = 1.)
                           : _topology( topology), 
                            _undefShell(undefShell),
                           _mass_distribution( mass_distribution ),
                           _factor( factor ),
                           _gravityDirection( gravityDirection) {}

	void apply(const VectorType& defShell, VectorType& Dest) const {

  if( defShell.size() != _undefShell.size() ){
      std::cerr << "size of deformed = " << defShell.size() << " vs. size of undeformed = " << _undefShell.size() << std::endl;
      throw BasicException( "GravitationalEnergyGradientDef::apply(): sizes dont match!");
    }

  if (_mass_distribution.size() != _topology.getNumVertices()){
    std::cerr << "size of masses = " << _mass_distribution.size() << " vs. size of vertices = " << _topology.getNumVertices() << std::endl;
    throw BasicException( "GravitationalEnergyGradientDef::apply(): sizes dont match!");
  }

  if( _gravityDirection.size() != 3 ){
      std::cerr << "size of gravity direction = " << _gravityDirection.size() << std::endl;
      throw BasicException( "GravitationalEnergy::apply(): gravity direction must be 3-dimensional!");
  }

  if( Dest.size() != _undefShell.size() )
    Dest.resize( _undefShell.size() );

  Dest.setZero();

  typename DefaultConfigurator::SparseMatrixType MassDistributionMatrix = _mass_distribution.asDiagonal();


  typename DefaultConfigurator::SparseMatrixType MassMatrix;
  computeStiffnessMatrix<DefaultConfigurator>( _topology, _undefShell, MassMatrix , false);

  // Now calculate the matrix mass_matrix * stiffness_matrix
  MatrixType AreaMassMatrix = MassDistributionMatrix*MassMatrix;

  for (int i = 0; i < 3; i++){
    for(int vertexIdx = 0; vertexIdx < _topology.getNumVertices(); vertexIdx++){
				Dest[i*_topology.getNumVertices() + vertexIdx] = _factor*_gravityDirection[i]*AreaMassMatrix.row(vertexIdx).sum();
    }
		}
  }
};
// Only the sign changed in the apply method. Everything else is identical when takin gradient w.r.t. undeformed instead of deformed configuration

//All second derivatives are zero, since the potential energy is linear in the deformation
#endif