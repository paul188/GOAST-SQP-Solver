#ifndef FOLD_DOFS_H
#define FOLD_DOFS_H

template <typename ConfiguratorType>
class FoldDofs : public BaseOp<const typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    public:
        virtual size_t getNumDofs() const = 0;
        virtual std::vector<RealType> getFoldDofs() const = 0;
        virtual bool isFoldVertex(const RealType coord_x, const RealType coord_y) const = 0;
        virtual void getFoldVertices(std::vector<int> &foldVertices) const = 0;
        virtual void getEdgeWeights(VectorType &edge_weights) const = 0;
        ~FoldDofs() = default;
};

template<typename ConfiguratorType>
class FoldDofsGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    public:
        virtual void apply(const VectorType &Geometry, VectorType& dest) const = 0;
};


// Domain type template is VectorType but in this case, the vector has one entry
template <typename ConfiguratorType>
class FoldDofsSimpleLine : public FoldDofs<ConfiguratorType> {

    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    
    public:
        FoldDofsSimpleLine(const MeshTopologySaver &plateTopol, const VectorType &Geom):_plateTopol(plateTopol)
        {
            _foldDofs.resize(1);
            _foldDofs[0] = 0;

            for( int i = 0; i < _plateTopol.getNumVertices(); i++ ){
                VecType coords;
                getXYZCoord<VectorType, VecType>( Geom, coords, i);
                if( isFoldVertex(coords[0],coords[1]) ){
                    _foldVertices.push_back( i );
                }
            }

            _edge_weights.resize(_plateTopol.getNumEdges());
            _edge_weights = VectorType::Ones(_plateTopol.getNumEdges());
            for(int edgeIdx = 0; edgeIdx < _plateTopol.getNumEdges(); edgeIdx++){
                int node_i = plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
                int node_j = plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

                VecType coords_i, coords_j;
                getXYZCoord<VectorType, VecType>( Geom, coords_i, node_i);
                getXYZCoord<VectorType, VecType>( Geom, coords_j, node_j);

                if(isFoldVertex(coords_i[0],coords_i[1]) && isFoldVertex(coords_j[0],coords_j[1])){
                    _edge_weights[edgeIdx] = 0;
                }
            }

        }

        void apply(const VectorType &t, VectorType &Dest) const{
            // Translate the vertices of the fold with one parameter t into the x-direction
            if(Dest.size() != 3*_plateTopol.getNumVertices()){
                throw std::invalid_argument("Geometry vector has wrong size");
            }

            for(int i = 0; i < _foldVertices.size(); i++){
                VecType coords;
                getXYZCoord<VectorType, VecType>( Dest, coords, _foldVertices[i]);
                coords[0] += t[0];
                _foldDofs[0] += t[0];
                setXYZCoord<VectorType, VecType>( Dest, coords, _foldVertices[i]);
            }
        }

        bool isFoldVertex(const RealType coord_x, const RealType coord_y) const
        {
            const RealType tolerance = 1e-5;  // Adjust as necessary
            return std::abs(coord_x - 0.5) < tolerance;
        }

        size_t getNumDofs() const {
            return 1;
        }

        std::vector<RealType> getFoldDofs() const {
            return _foldDofs;
        }

        void getFoldVertices(std::vector<int> &foldVertices) const{
            foldVertices = _foldVertices;
        }

        void getEdgeWeights(VectorType &edge_weights) const{
            edge_weights = _edge_weights;
        }

    private:
        const MeshTopologySaver &_plateTopol;
        std::vector<int> _foldVertices;
        VectorType _edge_weights;
        mutable std::vector<RealType> _foldDofs;
};

template <typename ConfiguratorType>
class FoldDofsSimpleLineGradient : public FoldDofsGradient<ConfiguratorType> {
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    protected:
        const MeshTopologySaver &_plateTopol;
        const std::vector<int> &_foldVertices;

    public:
        FoldDofsSimpleLineGradient(const MeshTopologySaver &plateTopol,const std::vector<int> &foldVertices): 
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

#endif