#include <goast/Core.h>

template<typename ConfiguratorType>
class DirichletSmoother{
    protected:
        typedef DefaultConfigurator::VectorType VectorType;
        typedef DefaultConfigurator::VecType VecType;
        typedef DefaultConfigurator::RealType RealType;

        const VectorType &baseGeometry;
        const VectorType &activeGeometry;
        const std::vector<int> &bdryMask;
        const MeshTopologySaver& plateTopol;

    public:
        DirichletSmoother(const VectorType &baseGeometry, const VectorType &activeGeometry, const std::vector<int> &bdryMask, const MeshTopologySaver &plateTopol)
        : baseGeometry(baseGeometry), activeGeometry(activeGeometry), bdryMask(bdryMask), plateTopol(plateTopol) {}

        void setActiveGeometry(const VectorType &activeGeometry){
            this-> activeGeometry = activeGeometry;
        }

        void setBoundaryMask(const std::vector<int> &bdryMask){
            this-> bdryMask = bdryMask;
        }

        void setBaseGeometry(const VectorType &baseGeometry){
            this-> baseGeometry = baseGeometry;
        }

        void setMeshTopologySaver(const MeshTopologySaver &plateTopol){
            this-> plateTopol = plateTopol;
        }

        // smoothing the mesh by solving a Poisson problem with boundary Conditions (bdryMask) with Stiffness matrix 
        // computed w.r.t. the baseGeometry on the active Geometry. The Dirichlet energy is returned
        RealType smoothMesh(VectorType &resultGeometry){
            typename DefaultConfigurator::SparseMatrixType StiffnessMatrix;
            computeStiffnessMatrix<DefaultConfigurator>( plateTopol, baseGeometry, StiffnessMatrix );
            applyMaskToMajor<typename DefaultConfigurator::SparseMatrixType>( bdryMask, StiffnessMatrix );

            // set up right hand side and mask
            VectorType rhs( activeGeometry );
            rhs.setZero();
            for( int i = 0; i < bdryMask.size(); i++ )
                rhs[bdryMask[i]] = activeGeometry[bdryMask[i]];

            // set up linear system and solve
            LinearSolver<DefaultConfigurator> directSolver;
            directSolver.prepareSolver( StiffnessMatrix );
            directSolver.backSubstitute( rhs, resultGeometry );
            // get final energy value
            computeStiffnessMatrix<DefaultConfigurator>( plateTopol, baseGeometry, StiffnessMatrix );
            VectorType temp = StiffnessMatrix * resultGeometry;

            return (0.5 * temp.dot( resultGeometry ));
        }
};