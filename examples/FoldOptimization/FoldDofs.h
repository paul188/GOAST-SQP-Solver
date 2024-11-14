template <typename ConfiguratorType>
// Domain type template is VectorType but in this case, the vector has one entry
class FoldDofs : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    
    protected:
        const VectorType &_Geometry;
        const std::vector<int> &_foldVertices;
    public:
        FoldDofs(const VectorType& Geometry,const std::vector<int> &foldVertices):
        _Geometry(Geometry), _foldVertices(foldVertices){}

        void apply(const VectorType &t, VectorType &Dest) const override{
            // Translate the vertices of the fold with one parameter t
            if(Dest.size() != _Geometry.size()){
                Dest.resize(_Geometry.size());
            }
            Dest = _Geometry;
            for(int i = 0; i < _foldVertices.size(); i++){
                VecType coords;
                getXYZCoord<VectorType, VecType>( Dest, coords, _foldVertices[i]);
                coords[0] += t[0];
                setXYZCoord<VectorType, VecType>( Dest, coords, _foldVertices[i]);
            }
        }
};

template <typename ConfiguratorType>
class FoldDofsGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    protected:
        const MeshTopologySaver &_plateTopol;
        const std::vector<int> &_foldVertices;

    public:
        FoldDofsGradient(const MeshTopologySaver &plateTopol,const std::vector<int> &foldVertices): 
        _plateTopol(plateTopol),_foldVertices(foldVertices){}
        // Gradient of translating all fold vertices with one parameter t into the x-direction
        // as in translateFoldVerticesAsOne
        // Geometry is the geometry w.r.t. which the smoothing is done, i.e. w.r.t.
        // which the Stiffness matrix is computed
        // Usually not the reference geometry, but the basic mesh
        void apply(const VectorType &Geometry, VectorType& dest) const override{
            // Compute the gradient of the translation of the vertices of the fold with one parameter t
           if(dest.size() != Geometry.size()){
                dest.resize(Geometry.size());
            }
            dest.setZero();

            VectorType indicator_dof = VectorType::Zero(Geometry.size());

            for(int i = 0; i < _foldVertices.size(); i++){
                // Set the x components to 1
                indicator_dof[_foldVertices[i]] = 1;
            }

            typename ConfiguratorType::SparseMatrixType StiffnessMatrix;
            computeStiffnessMatrix<ConfiguratorType>(_plateTopol,Geometry, StiffnessMatrix);
            
            LinearSolver<DefaultConfigurator> directSolver;
            directSolver.prepareSolver( StiffnessMatrix );
            directSolver.backSubstitute( indicator_dof, dest );
        }
};