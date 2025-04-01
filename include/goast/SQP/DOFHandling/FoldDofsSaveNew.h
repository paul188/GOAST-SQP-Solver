#pragma once

#include "../Utils/SparseMat.h"
#include <goast/Smoothers.h>
#include <goast/Core.h>

template <typename ConfiguratorType>
class FoldDofs{
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;
    public:
        virtual size_t getNumDofs() const = 0;
        virtual bool isFoldVertex(const RealType coord_x, const RealType coord_y) const = 0;
        virtual bool isFoldEdge(const int edgeIdx) const = 0;
        virtual void getFoldVertices(std::vector<int> &foldVertices) const = 0;
        virtual void getEdgeWeights(VectorType &edge_weights) const = 0;
        ~FoldDofs() = default;
        virtual void apply(const VectorType &t, VectorType& dest) const = 0;
        FoldDofs(const MeshTopologySaver &plateTopol,
                 const VectorType &plateGeomInitial,
                 const VectorType &plateGeomRef_basic,
                 const std::vector<int> &bdryMaskRef)
                :_plateTopol(plateTopol)
                , _plateGeomInitial(plateGeomInitial)
                , _plateGeomRef_basic(plateGeomRef_basic)
                , _bdryMaskRef(bdryMaskRef)
            {}

        void initialize_folds_edges() {
            // fold vertices are specified using the isFoldVertex function
            // relative to the basic mesh -> plateGeomInitial
            for( int i = 0; i < _plateTopol.getNumVertices(); i++ ){
                VecType coords;
                getXYZCoord<VectorType, VecType>( _plateGeomInitial, coords, i);
                if( isFoldVertex(coords[0],coords[1]) ){
                    _foldVertices.push_back( i );
                }
            }

            _edge_weights.resize(_plateTopol.getNumEdges());
            _edge_weights = VectorType::Ones(_plateTopol.getNumEdges());
            for(int edgeIdx = 0; edgeIdx < _plateTopol.getNumEdges(); edgeIdx++){
                if(isFoldEdge(edgeIdx)){
                    _edge_weights[edgeIdx] = 0;
                }
            }
        }

    protected:
        const MeshTopologySaver &_plateTopol;
        const VectorType &_plateGeomInitial;
        const VectorType &_plateGeomRef_basic;
        const std::vector<int> &_bdryMaskRef;
        std::vector<int> _foldVertices;
        VectorType _edge_weights;
};

template<typename ConfiguratorType>
class FoldDofsGradient{
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;

    public:
        virtual void apply(const VectorType &t, MatrixType& dest) const = 0;
        ~FoldDofsGradient() = default;
        FoldDofsGradient(const MeshTopologySaver &plateTopol,
                         const std::vector<int> &bdryMaskRef,
                         const VectorType &plateGeomInitial,
                         const std::vector<int> &foldVertices)
                         : _plateTopol( plateTopol ),
                           _bdryMaskRef( bdryMaskRef ),
                           _plateGeomInitial( plateGeomInitial ),
                           _foldVertices( foldVertices ){}
    protected:
        const MeshTopologySaver &_plateTopol;
        const std::vector<int> &_foldVertices;
        const std::vector<int> &_bdryMaskRef;
        const VectorType &_plateGeomInitial;
};

// Domain type template is VectorType but in this case, the vector has one entry
template <typename ConfiguratorType>
class FoldDofsSimpleLine : public FoldDofs<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;

    public:
        // We start from a certain reference geometry relative to which we translate the foldDofs. This is plateGeomRef_basic
        // This reference geometry should not be updated during the optimization
        // No new object needs to be created during optimization loop
        FoldDofsSimpleLine(const MeshTopologySaver &plateTopol, const VectorType &plateGeomInitial, const VectorType &plateGeomRef_basic, const std::vector<int> &bdryMaskRef)
        : FoldDofs<ConfiguratorType>(plateTopol,plateGeomInitial,plateGeomRef_basic,bdryMaskRef)
        {
            this->initialize_folds_edges();
        }

        void apply(const VectorType &t, VectorType &Dest) const override{
            // Translate the vertices of the fold with one parameter t into the x-direction
            if(Dest.size() != 3*this->_plateTopol.getNumVertices()){
                Dest.resize(3*this->_plateTopol.getNumVertices());
            }

            Dest = this->_plateGeomRef_basic;

            for(int i = 0; i < this->_foldVertices.size(); i++){
                VecType coords;
                getXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
                coords[1] += t[0];
                setXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
            }

            DirichletSmoother<ConfiguratorType> smoother(this->_plateGeomInitial, this->_bdryMaskRef, this->_plateTopol);
            smoother.apply(Dest,Dest);
        }

        bool isFoldVertex(const RealType coord_x, const RealType coord_y) const
        {
            const RealType tolerance = 1e-5;  // Adjust as necessary
            return std::abs(coord_y - 2.0) < tolerance;
        }

        bool isFoldEdge(const int edgeIdx) const
        {
               int node_i = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
               int node_j = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

                VecType coords_i, coords_j;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_i, node_i);
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_j, node_j);

                if(isFoldVertex(coords_i[0],coords_i[1]) && isFoldVertex(coords_j[0],coords_j[1])){
                    return true;
                }
            return false;
        }

        size_t getNumDofs() const {
            return 1;
        }

        void getFoldVertices(std::vector<int> &foldVertices) const{
            foldVertices = this->_foldVertices;
        }

        void getEdgeWeights(VectorType &edge_weights) const{
            edge_weights = this->_edge_weights;
        }
};

template <typename ConfiguratorType>
class FoldDofsSimpleLineGradient : public FoldDofsGradient<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    public:
        FoldDofsSimpleLineGradient(const MeshTopologySaver &plateTopol,
                                   const std::vector<int> &bdryMaskRef,
                                   const VectorType &plateGeomInitial,
                                   const std::vector<int> &foldVertices)
                                   : FoldDofsGradient<ConfiguratorType>( plateTopol, bdryMaskRef, plateGeomInitial, foldVertices ){}

        // Gradient of translating all fold vertices with one parameter t into the x-direction
        // as in translateFoldVerticesAsOne
        // Geometry is the geometry w.r.t. which the smoothing is done, i.e. w.r.t.
        // which the Stiffness matrix is computed
        // Usually not the reference geometry, but the basic mesh
        void apply(const VectorType &t, MatrixType& Dest) const override{
            // Compute the gradient of the translation of the vertices of the fold with one parameter t
           
           VectorType dest;
           
           if(dest.size() != this->_plateGeomInitial.size()){
                dest.resize(this->_plateGeomInitial.size());
            }
            dest.setZero();
        
            VectorType indicator_dof = VectorType::Zero(this->_plateGeomInitial.size());

            for(int i = 0; i < this->_foldVertices.size(); i++){
                // Set the Y components of foldVertices to 1
                indicator_dof[this->_plateTopol.getNumVertices() + this->_foldVertices[i]] = 1;
            }

            typename ConfiguratorType::SparseMatrixType StiffnessMatrix;
            computeStiffnessMatrix<ConfiguratorType>(this->_plateTopol,this->_plateGeomInitial, StiffnessMatrix);
            applyMaskToMajor<typename DefaultConfigurator::SparseMatrixType>( this->_bdryMaskRef, StiffnessMatrix );

            LinearSolver<DefaultConfigurator> directSolver;
            directSolver.prepareSolver( StiffnessMatrix );
            directSolver.backSubstitute( indicator_dof, dest );

            if(!(StiffnessMatrix*dest).isApprox(indicator_dof)){
                std::cerr<<"Gradient in Fold Dofs Gradient is incorrect !\n";
            }
            Dest = convertVecToSparseMat(dest);
        }
};

template <typename ConfiguratorType>
class FoldDofsFreeLine : public FoldDofs<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;

    private:
        bool parallel_x;
        RealType initial_line_pos;

    public:
        // We start from a certain reference geometry relative to which we translate the foldDofs. This is plateGeomRef_basic
        // This reference geometry should not be updated during the optimization
        // No new object needs to be created during optimization loop
        // parallel_x specifies whether the line is oriented parallel to the x-direction (true) or y-direction (false)
        FoldDofsFreeLine(const MeshTopologySaver &plateTopol, const VectorType &plateGeomInitial, const VectorType &plateGeomRef_basic, const std::vector<int> &bdryMaskRef, RealType initial_line_pos, bool parallel_x)
        : FoldDofs<ConfiguratorType>(plateTopol,plateGeomInitial,plateGeomRef_basic,bdryMaskRef)
        {
            this->initialize_folds_edges();
            this->parallel_x = parallel_x;
        }

        void apply(const VectorType &t, VectorType &Dest) const override{
            // Translate each vertex along the fold line with an individual parameter t
            if(Dest.size() != 3*this->_plateTopol.getNumVertices()){
                Dest.resize(3*this->_plateTopol.getNumVertices());
            }

            Dest = this->_plateGeomRef_basic;

            for(int i = 0; i < this->_foldVertices.size(); i++){
                VecType coords;
                getXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
                // This is the difference to FoldDofsSimpleLine
                // -> have one param t[i] for each fold vertex
                if(parallel_x){
                    coords[1] += t[i];
                }else{
                    coords[0] += t[i];
                }
                setXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
            }

            DirichletSmoother<ConfiguratorType> smoother(this->_plateGeomInitial, this->_bdryMaskRef, this->_plateTopol);
            smoother.apply(Dest,Dest);
        }

        // checks whether the vertex is a fold vertex in the geometry plateGeom_basic
        bool isFoldVertex(const RealType coord_x, const RealType coord_y) const
        {
            const RealType tolerance = 1e-5;  // Adjust as necessary
            if(parallel_x){
                return std::abs(coord_y - initial_line_pos) < tolerance;
            }else{
                return std::abs(coord_x - initial_line_pos) < tolerance;
            }
        }

        bool isFoldEdge(const int edgeIdx) const
        {
            int node_i = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
            int node_j = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

            VecType coords_i, coords_j;
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_i, node_i);
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_j, node_j);

            if(isFoldVertex(coords_i[0],coords_i[1]) && isFoldVertex(coords_j[0],coords_j[1])){
                return true;
            }
            return false;
        }

        size_t getNumDofs() const {
            return this->_foldVertices.size();
        }

        void getFoldVertices(std::vector<int> &foldVertices) const{
            foldVertices = this->_foldVertices;
        }

        void getEdgeWeights(VectorType &edge_weights) const{
            edge_weights = this->_edge_weights;
        }
};

template <typename ConfiguratorType>
class FoldDofsFreeLineGradient : public FoldDofsGradient<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;

    private:
        bool parallel_x;

    public:
        FoldDofsFreeLineGradient(const MeshTopologySaver &plateTopol,
                                   const std::vector<int> &bdryMaskRef,
                                   const VectorType &plateGeomInitial,
                                   const std::vector<int> &foldVertices,
                                   bool parallel_x)
                                   : FoldDofsGradient<ConfiguratorType>(plateTopol,bdryMaskRef,plateGeomInitial,foldVertices)
                                   {
                                        this->parallel_x = parallel_x;
                                   }

        void apply(const VectorType &t, MatrixType& Dest) const override{
            // Compute the gradient of the translation of the vertices of the fold with parameters t_i       
           
            VectorType dest;
            dest.resize(this->_plateGeomInitial.size());
           
            if(Dest.rows() != this->_plateGeomInitial.size() || Dest.cols() != this->_foldVertices.size()){
                Dest.resize(this->_plateGeomInitial.size(), this->_foldVertices.size());
            }

            dest.setZero();
            Dest.setZero();

            MatrixType indicator_dof(this->_plateGeomInitial.size(),this->_foldVertices.size());
            indicator_dof.setZero();

            std::vector<Eigen::Triplet<RealType>> triplets;
            if(parallel_x){
                for (int i = 0; i < this->_foldVertices.size(); i++) {
                    VecType coords;
                    getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords, this->_foldVertices[i]);
                    triplets.push_back(Eigen::Triplet<RealType>(this->_plateTopol.getNumVertices() + this->_foldVertices[i], i, 1.0));
                }
            }else{
                for (int i = 0; i < this->_foldVertices.size(); i++) {
                    VecType coords;
                    getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords, this->_foldVertices[i]);
                    triplets.push_back(Eigen::Triplet<RealType>(this->_foldVertices[i], i, 1.0));
                }
            }

            indicator_dof.setFromTriplets(triplets.begin(), triplets.end());

            typename ConfiguratorType::SparseMatrixType StiffnessMatrix;
            computeStiffnessMatrix<ConfiguratorType>(this->_plateTopol,this->_plateGeomInitial, StiffnessMatrix);
            applyMaskToMajor<typename ConfiguratorType::SparseMatrixType>( this->_bdryMaskRef, StiffnessMatrix );

            Eigen::BiCGSTAB<MatrixType> solver;
            solver.compute(StiffnessMatrix);

            std::vector<Eigen::Triplet<RealType>> tripletList;
            // solve the equation L*dest = indicator_dof for each column of indicator_dof
            for(int i = 0; i < this->_foldVertices.size(); i++)
            {
                dest = solver.solve(indicator_dof.col(i));
                assignSparseBlockInplace(Dest, convertVecToSparseMat(dest), 0, i, tripletList);
                dest.setZero();
            }
        }
};

template <typename ConfiguratorType>
class FoldDofsCross : public FoldDofs<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    
    public:
        FoldDofsCross(const MeshTopologySaver &plateTopol, const VectorType &plateGeomInitial, const VectorType &plateGeomRef_basic, const std::vector<int> &bdryMaskRef)
        : FoldDofs<ConfiguratorType>(plateTopol,plateGeomInitial,plateGeomRef_basic,bdryMaskRef)
        {
            for( int i = 0; i < this->_plateTopol.getNumVertices(); i++ ){
                VecType coords;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords, i);
                if(isLine1(coords[0], coords[1])){
                    _line1Vertices.push_back(i);
                }
                if(isLine2(coords[0], coords[1])){
                    _line2Vertices.push_back(i);
                }
            }

            this->initialize_folds_edges();
        }

        void apply(const VectorType &t, VectorType &Dest) const override{
            // Translate each vertex along the fold line with an individual parameter t
            if(Dest.size() != 3*this->_plateTopol.getNumVertices()){
                Dest.resize(3*this->_plateTopol.getNumVertices());
            }

            Dest = this->_plateGeomRef_basic;

            for(int i = 0; i < _line1Vertices.size(); i++){
                VecType coords;
                getXYZCoord<VectorType, VecType>( Dest, coords, _line1Vertices[i]);
                coords[0] += t[0];
                setXYZCoord<VectorType, VecType>( Dest, coords, _line1Vertices[i]);
            }
            
            for(int i = 0; i < _line2Vertices.size(); i++){
                VecType coords;
                getXYZCoord<VectorType, VecType>( Dest, coords, _line2Vertices[i]);
                coords[1] += t[1];
                setXYZCoord<VectorType, VecType>( Dest, coords, _line2Vertices[i]);
            }

            DirichletSmoother<ConfiguratorType> smoother(this->_plateGeomInitial, this->_bdryMaskRef, this-> _plateTopol);
            smoother.apply(Dest,Dest);
        }

        // the cross consists of two lines at x=0.5, resp. y=0.5
        bool isLine1(const RealType coord_x, const RealType coord_y) const
        {
            const RealType tolerance = 1e-5;  // Adjust as necessary
            return std::abs(coord_x - 0.625) < tolerance;
        }

        bool isLine2(const RealType coord_x, const RealType coord_y) const
        {
            const RealType tolerance = 1e-5;  // Adjust as necessary
            return std::abs(coord_y - 0.625) < tolerance;
        }

        bool isFoldVertex(const RealType coord_x, const RealType coord_y) const
        {
            return (isLine1(coord_x,coord_y) || isLine2(coord_x,coord_y));
        }

        bool isFoldEdge(const int edgeIdx) const
        {
            int node_i = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
            int node_j = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

            VecType coords_i, coords_j;
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_i, node_i);
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_j, node_j);

            if(isLine1(coords_i[0],coords_i[1]) && isLine1(coords_j[0],coords_j[1])){
                return true;
            }
            if(isLine2(coords_i[0],coords_i[1]) && isLine2(coords_j[0],coords_j[1])){
                return true;
            }
            return false;
        }

        void getFoldVertices(std::vector<int> &foldVertices) const{
            foldVertices = this->_foldVertices;
        }

        void getLine1Vertices(std::vector<int> &line1Vertices) const{
            line1Vertices = _line1Vertices;
        }

        void getLine2Vertices(std::vector<int> &line2Vertices) const{
            line2Vertices = _line2Vertices;
        }

        void getEdgeWeights(VectorType &edge_weights) const{
            edge_weights = this->_edge_weights;
        }

        size_t getNumDofs() const {
            return 2;
        }

        protected:
            std::vector<int> _line1Vertices;
            std::vector<int> _line2Vertices;
};

template <typename ConfiguratorType>
class FoldDofsCrossGradient : public FoldDofsGradient<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    protected:
        std::vector<int> _line1Vertices;
        std::vector<int> _line2Vertices;

    public:

        FoldDofsCrossGradient(const MeshTopologySaver &plateTopol,
                                   const std::vector<int> &bdryMaskRef,
                                   const VectorType &plateGeomInitial,
                                   const std::vector<int> &foldVertices,
                                   const std::vector<int> &line1Vertices,
                                   const std::vector<int> &line2Vertices)
                                      : FoldDofsGradient<ConfiguratorType>(plateTopol,bdryMaskRef,plateGeomInitial,foldVertices),
                                     _line1Vertices( line1Vertices ),
                                     _line2Vertices( line2Vertices ){}


        void apply(const VectorType &t, MatrixType& Dest) const override{
            // Compute the gradient of the translation of the vertices of the fold with parameters t_i       
           
            VectorType dest;
            dest.resize(this->_plateGeomInitial.size());
           
            if(Dest.rows() != this->_plateGeomInitial.size() || Dest.cols() != 2){
                Dest.resize(this->_plateGeomInitial.size(), 2);
            }

            dest.setZero();
            Dest.setZero();

            MatrixType indicator_dof(this->_plateGeomInitial.size(),2);
            indicator_dof.setZero();

            std::vector<Eigen::Triplet<RealType>> triplets;
            for (int i = 0; i < _line1Vertices.size(); i++) {
                VecType coords;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords, _line1Vertices[i]);
                triplets.push_back(Eigen::Triplet<RealType>(_line1Vertices[i], 0, 1.0));
            }

            for(int i = 0; i < _line2Vertices.size(); i++) {
                VecType coords;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords, _line2Vertices[i]);
                triplets.push_back(Eigen::Triplet<RealType>(this->_plateTopol.getNumVertices() + _line2Vertices[i], 1, 1.0));
            }

            indicator_dof.setFromTriplets(triplets.begin(), triplets.end());

            typename ConfiguratorType::SparseMatrixType StiffnessMatrix;
            computeStiffnessMatrix<ConfiguratorType>(this->_plateTopol, this->_plateGeomInitial, StiffnessMatrix);
            applyMaskToMajor<typename ConfiguratorType::SparseMatrixType>( this->_bdryMaskRef, StiffnessMatrix );

            LinearSolver<DefaultConfigurator> directSolver;
            directSolver.prepareSolver( StiffnessMatrix );

            std::vector<Eigen::Triplet<RealType>> tripletList;
            // solve the equation L*dest = indicator_dof for each column of indicator_dof
            for(int i = 0; i < 2; i++)
            {
                directSolver.backSubstitute( indicator_dof.col(i), dest );
                assignSparseBlockInplace(Dest, convertVecToSparseMat(dest), 0, i, tripletList);
                dest.setZero();
            }
        }

};

template <typename ConfiguratorType>
class FoldDofsArcLine : public FoldDofs<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;
    
    public:
        FoldDofsArcLine(const MeshTopologySaver &plateTopol, const VectorType &plateGeomInitial, const VectorType &plateGeomRef_basic, const std::vector<int> &bdryMaskRef)
          : FoldDofs<ConfiguratorType>(plateTopol,plateGeomInitial,plateGeomRef_basic,bdryMaskRef)
        {
            this->initialize_folds_edges();
        }

        void apply(const VectorType &t, VectorType &Dest) const override{
            if(Dest.size() != 3*this->_plateTopol.getNumVertices()){
                Dest.resize(3*this->_plateTopol.getNumVertices());
            }

            Dest = this->_plateGeomRef_basic;
            for(int i = 0; i < this->_foldVertices.size(); i++)
            {
                VecType coords;
                getXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
                coords[1] = 0.25 + t[0]*(0.25 - std::pow(coords[0] - 0.5, 2.0));
                setXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
            }

            DirichletSmoother<ConfiguratorType> smoother(this->_plateGeomInitial, this->_bdryMaskRef, this->_plateTopol);
            smoother.apply(Dest,Dest);
        }

        bool isFoldVertex(const RealType coord_x, const RealType coord_y) const
        {
            RealType tolerance = 1e-5;
            return (std::abs(coord_y - 0.25) < tolerance);
        }

        bool isFoldEdge(const int edgeIdx) const
        {
            int node_i = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
            int node_j = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

            VecType coords_i, coords_j;
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_i, node_i);
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_j, node_j);

            if(isFoldVertex(coords_i[0],coords_i[1]) && isFoldVertex(coords_j[0],coords_j[1])){
                return true;
            }
            return false;
        }

        void getFoldVertices(std::vector<int> &foldVertices) const{
            foldVertices = this->_foldVertices;
        }

        void getEdgeWeights(VectorType &edge_weights) const{
            edge_weights = this->_edge_weights;
        }

        size_t getNumDofs() const {
            return 1;
        }
};

template<typename ConfiguratorType>
class FoldDofsArcLineGradient : public FoldDofsGradient<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    public:
        FoldDofsArcLineGradient(const MeshTopologySaver &plateTopol,
                                   const std::vector<int> &bdryMaskRef,
                                   const VectorType &plateGeomInitial,
                                   const std::vector<int> &foldVertices)
                                   : FoldDofsGradient<ConfiguratorType>(plateTopol,bdryMaskRef,plateGeomInitial,foldVertices){}

        void apply(const VectorType &t, MatrixType& Dest) const override{
           
            VectorType dest, rhs;
            dest.resize(this->_plateGeomInitial.size());
            rhs.resize(this->_plateGeomInitial.size());
           
            if(Dest.rows() != this->_plateGeomInitial.size() || Dest.cols() != 1){
                Dest.resize(this->_plateGeomInitial.size(), 1);
            }

            rhs.setZero();
            dest.setZero();
            Dest.setZero();

            for(int i = 0; i < this->_foldVertices.size(); i++)
            {
                VecType coords;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords, this->_foldVertices[i]);
                rhs[this->_foldVertices[i] + this->_plateTopol.getNumVertices()] = (0.25 - std::pow(coords[0] - 0.5, 2.0));
            }

            typename ConfiguratorType::SparseMatrixType StiffnessMatrix;
            computeStiffnessMatrix<ConfiguratorType>(this->_plateTopol, this->_plateGeomInitial, StiffnessMatrix);
            applyMaskToMajor(this->_bdryMaskRef, StiffnessMatrix );

            LinearSolver<DefaultConfigurator> directSolver;
            directSolver.prepareSolver( StiffnessMatrix );
            directSolver.backSubstitute( rhs, dest );
            if (dest.array().isNaN().any()) {
                std::cerr << "Warning: Solution contains NaN values!" << std::endl;
            }
            Dest = convertVecToSparseMat(dest);
        }
};