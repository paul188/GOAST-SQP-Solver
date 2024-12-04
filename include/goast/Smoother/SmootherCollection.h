#include <goast/Core.h>

template<typename ConfiguratorType>
class DirichletSmoother : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    protected:
        typedef DefaultConfigurator::VectorType VectorType;
        typedef DefaultConfigurator::VecType VecType;
        typedef DefaultConfigurator::RealType RealType;

        const VectorType &_plateGeomInitial;
        const std::vector<int> &_bdryMask;
        const MeshTopologySaver& _plateTopol;

    public:
        DirichletSmoother(const VectorType &plateGeomInitial, const std::vector<int> &bdryMask, const MeshTopologySaver &plateTopol)
        : _plateGeomInitial(plateGeomInitial), _bdryMask(bdryMask), _plateTopol(plateTopol) {}

        // smoothing the mesh by solving a Poisson problem with boundary Conditions (bdryMask) with Stiffness matrix 
        // computed w.r.t. the baseGeometry on the active Geometry. The Dirichlet energy is returned
        void apply(const VectorType& ActiveGeometry, VectorType &Dest) const override{
            
            if(ActiveGeometry.size() != _plateGeomInitial.size()){
                std::cerr << "size of active = " << ActiveGeometry.size() << " vs. size of inactive = " << _plateGeomInitial.size() << std::endl;
                throw BasicException( "DirichletSmoother::smoothMesh(): sizes dont match!");
            }

            typename DefaultConfigurator::SparseMatrixType StiffnessMatrix;
            computeStiffnessMatrix<DefaultConfigurator>( _plateTopol, _plateGeomInitial, StiffnessMatrix );
            applyMaskToMajor<typename DefaultConfigurator::SparseMatrixType>( _bdryMask, StiffnessMatrix );

            if(Dest.size() != ActiveGeometry.size()){
                Dest.resize(ActiveGeometry.size());
            }

            // set up right hand side and mask
            VectorType rhs( ActiveGeometry.size() );
            rhs.setZero();
            for( int i = 0; i < _bdryMask.size(); i++ )
                rhs[_bdryMask[i]] = ActiveGeometry[_bdryMask[i]];

            // set up linear system and solve
            LinearSolver<DefaultConfigurator> directSolver;
            directSolver.prepareSolver( StiffnessMatrix );
            directSolver.backSubstitute( rhs, Dest );
        }
};