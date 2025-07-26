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
        // is the line parallel to the x axis or y axis
        const bool is_x_parallel;
        // how far away is the line from x-/resp. y-axis
        const RealType start_value;

    public:
        // We start from a certain reference geometry relative to which we translate the foldDofs. This is plateGeomRef_basic
        // This reference geometry should not be updated during the optimization
        // No new object needs to be created during optimization loop
        // parallel_x specifies whether the line is oriented parallel to the x-direction (true) or y-direction (false)
        FoldDofsFreeLine(const MeshTopologySaver &plateTopol, const VectorType &plateGeomInitial, const VectorType &plateGeomRef_basic, const std::vector<int> &bdryMaskRef, const bool is_x_parallel, const RealType start_value) 
        : FoldDofs<ConfiguratorType>(plateTopol,plateGeomInitial,plateGeomRef_basic,bdryMaskRef), is_x_parallel(is_x_parallel), start_value(start_value)
        {
            this->initialize_folds_edges();
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
                if(is_x_parallel){
                    // This is the difference to FoldDofsSimpleLine
                    // -> have one param t[i] for each fold vertex
                    coords[1] += t[i];
                }
                else{
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
            if(is_x_parallel){
                return std::abs(coord_y - start_value) < tolerance;
            }else{
                return std::abs(coord_x - start_value) < tolerance;
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
        // is the line parallel to the x axis or y axis
        const bool is_x_parallel;

    public:
        FoldDofsFreeLineGradient(const MeshTopologySaver &plateTopol,
                                   const std::vector<int> &bdryMaskRef,
                                   const VectorType &plateGeomInitial,
                                   const std::vector<int> &foldVertices,
                                   const bool is_x_parallel)
                                   : FoldDofsGradient<ConfiguratorType>(plateTopol,bdryMaskRef,plateGeomInitial,foldVertices), is_x_parallel(is_x_parallel)
                                   {}

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
            for (int i = 0; i < this->_foldVertices.size(); i++) {
                if(is_x_parallel){
                    triplets.push_back(Eigen::Triplet<RealType>(this->_plateTopol.getNumVertices() + this->_foldVertices[i], i, 1.0));
                }else{
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

template<typename ConfiguratorType>
class FoldDofsSkewedCross : public FoldDofs<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    public:
        FoldDofsSkewedCross(const MeshTopologySaver &plateTopol,
                            const VectorType &plateGeomInitial,
                            const VectorType &plateGeomRef_basic,
                            const std::vector<int> &bdryMaskRef)
        : FoldDofs<ConfiguratorType>(plateTopol,plateGeomInitial,plateGeomRef_basic,bdryMaskRef)
        {
            this->initialize_folds_edges();
        }

        void apply(const VectorType &t, VectorType &Dest) const override{
            if(Dest.size() != 3*this->_plateTopol.getNumVertices()){
                Dest.resize(3*this->_plateTopol.getNumVertices());
            }

            Dest = this->_plateGeomRef_basic;

            // Reparametrise the parameter t such that it is in [-0.3, 0.3]
            RealType w = 0.3*std::tanh(t[0]);

            // Now, apply the translation to the fold vertices
            
            for(int i = 0; i < this->_foldVertices.size(); i++)
            {
                VecType coords, coords_initial;
                getXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_initial, this->_foldVertices[i]);
                // Check which part of the cross this vertex belongs to
                if(std::abs(coords_initial[1] - 0.5) < 1e-4)
                {
                    // first translate by -(0.5, 0.5)
                    coords[1] = 0.5-2*w*(coords_initial[0]-0.5);
                }
                else if(std::abs(coords_initial[0] - 0.5) < 1e-4)
                {
                    coords[0] = 0.5+2*w*(coords_initial[1]-0.5);
                }
                setXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
            }   

            auto exp_smoother = [](RealType x) -> RealType{
                return (1.0 - x);
            };
            

            
            for(int i = 0; i < this->_plateTopol.getNumVertices(); i++)
            {
                VecType coords_initial, coords;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_initial, i);
                getXYZCoord<VectorType, VecType>(this->_plateGeomInitial, coords, i);
                // Apply the exponential smoothing to the coordinates
                // Check in which quadrant we are
                if(coords_initial[0] <= 0.5 && coords_initial[1] <= 0.5)
                {
                    coords[0] = coords_initial[0] - exp_smoother(2*(0.5 - coords_initial[0]))*2*w*(0.5 - coords_initial[1]);
                    coords[1] = coords_initial[1] + exp_smoother(2*(0.5 - coords_initial[1]))*2*w*(0.5 - coords_initial[0]);
                }
                else if(coords_initial[0] <= 0.5 && coords_initial[1] >= 0.5)
                {
                    coords[0] = coords_initial[0] + exp_smoother(2*(0.5 - coords_initial[0]))*2*w*(coords_initial[1] - 0.5);
                    coords[1] = coords_initial[1] + exp_smoother(2*(coords_initial[1] - 0.5))*2*w*(0.5 - coords_initial[0]);
                }
                else if(coords_initial[0] >= 0.5 && coords_initial[1] >= 0.5)
                {
                    coords[0] = coords_initial[0] + exp_smoother(2*(coords_initial[0] - 0.5))*2*w*(coords_initial[1] - 0.5);
                    coords[1] = coords_initial[1] - exp_smoother(2*(coords_initial[1] - 0.5))*2*w*(coords_initial[0] - 0.5);
                }
                else if(coords_initial[0] >= 0.5 && coords_initial[1] <= 0.5)
                {
                    coords[0] = coords_initial[0] - exp_smoother(2*(coords_initial[0] - 0.5))*2*w*(0.5 - coords_initial[1]);
                    coords[1] = coords_initial[1] - exp_smoother(2*(0.5 - coords_initial[1]))*2*w*(coords_initial[0] - 0.5);
                }
                setXYZCoord<VectorType, VecType>( Dest, coords, i);
            }
        }

        bool isFoldVertex(const RealType coord_x, const RealType coord_y) const
        {
            if((std::abs(coord_x - 0.5) < 1e-4) || (std::abs(coord_y - 0.5) < 1e-4)){
                return true;
            }
            return false;
        }

        bool isFoldEdge(const int edgeIdx) const
        {
            int node_i = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
            int node_j = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

            VecType coords_i, coords_j;
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_i, node_i);
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_j, node_j);

            if((std::abs(coords_i[0] - 0.5) < 1e-4 && std::abs(coords_j[0] - 0.5) < 1e-4) || (std::abs(coords_i[1] - 0.5) < 1e-4 && std::abs(coords_j[1] - 0.5) < 1e-4)){
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
class FoldDofsSkewedCrossGradient : public FoldDofsGradient<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    public:
        FoldDofsSkewedCrossGradient(const MeshTopologySaver &plateTopol,
            const std::vector<int> &bdryMaskRef,
            const VectorType &plateGeomInitial,
            const std::vector<int> &foldVertices  
        ): FoldDofsGradient<ConfiguratorType>(plateTopol,bdryMaskRef,plateGeomInitial,foldVertices)
        {}

        void apply(const VectorType &t, MatrixType& Dest) const override{
            VectorType dest, rhs;
            dest.resize(this->_plateGeomInitial.size());
            rhs.resize(this->_plateGeomInitial.size());
           
            if(Dest.rows() != this->_plateGeomInitial.size() || Dest.cols() != 1){
                Dest.resize(this->_plateGeomInitial.size(), 1);
            }

            rhs.setZero();
            Dest.setZero();


            // derivative of the translation of the fold vertices
            // differentiate concatenation of rotation and scaling
            // -> first rotation, then scaling. both are linear operations
            // but at very first, differentiate how the angle changes with t
            // angle = atan(2*t) -> derivative 2/(1 + (2*t^2))

            RealType w = 0.3*std::tanh(t[0]);
            RealType dw_dt = 0.3*(1.0 - std::tanh(t[0])*std::tanh(t[0]));

            std::vector<Eigen::Triplet<RealType>> triplets;
            triplets.reserve(2*this->_plateTopol.getNumVertices());

            for(int i = 0; i < this->_foldVertices.size(); i++)
            {
                VecType coords_initial;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_initial, this->_foldVertices[i]);
                // Check which part of the cross this vertex belongs to
                if(std::abs(coords_initial[1] - 0.5) < 1e-4)
                {
                    triplets.emplace_back(this->_plateTopol.getNumVertices() + this->_foldVertices[i], 0, -2*(coords_initial[0]-0.5)*dw_dt);
                }
                else if(std::abs(coords_initial[0] - 0.5) < 1e-4)
                {
                    triplets.emplace_back(this->_foldVertices[i], 0, 2*(coords_initial[1]-0.5)*dw_dt);
                }
            }
            
            auto exp_smoother = [](RealType x) -> RealType{
                return (1.0 - x);
            };
            
            for(int i = 0; i < this->_plateTopol.getNumVertices(); i++)
            {
                VecType coords_initial, coords;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_initial, i);
                getXYZCoord<VectorType, VecType>(this->_plateGeomInitial, coords, i);

                if(std::abs(coords_initial[1] - 0.5) < 1e-4)
                {
                    continue;
                }
                else if(std::abs(coords_initial[0] - 0.5) < 1e-4)
                {
                    continue;
                }

                // Apply the exponential smoothing to the coordinates
                // Check in which quadrant we are
                if(coords_initial[0] <= 0.5 && coords_initial[1] <= 0.5)
                {
                    triplets.emplace_back(i, 0, - exp_smoother(2*(0.5 - coords_initial[0]))*2*dw_dt*(0.5 - coords_initial[1]));
                    triplets.emplace_back(this->_plateTopol.getNumVertices() + i, exp_smoother(2*(0.5 - coords_initial[1]))*2*dw_dt*(0.5 - coords_initial[0]));
                }
                else if(coords_initial[0] <= 0.5 && coords_initial[1] >= 0.5)
                {
                    triplets.emplace_back(i, 0, exp_smoother(2*(0.5 - coords_initial[0]))*2*dw_dt*(coords_initial[1] - 0.5));
                    triplets.emplace_back(this->_plateTopol.getNumVertices() + i, 0,  exp_smoother(2*(coords_initial[1] - 0.5))*2*dw_dt*(0.5 - coords_initial[0]));
                }
                else if(coords_initial[0] >= 0.5 && coords_initial[1] >= 0.5)
                {
                    triplets.emplace_back(i, 0, exp_smoother(2*(coords_initial[0] - 0.5))*2*dw_dt*(coords_initial[1] - 0.5));
                    triplets.emplace_back(this->_plateTopol.getNumVertices() + i,0, - exp_smoother(2*(coords_initial[1] - 0.5))*2*dw_dt*(coords_initial[0] - 0.5));
                }
                else if(coords_initial[0] >= 0.5 && coords_initial[1] <= 0.5)
                {
                    triplets.emplace_back(i, 0, - exp_smoother(2*(coords_initial[0] - 0.5))*2*dw_dt*(0.5 - coords_initial[1]));
                    triplets.emplace_back(this->_plateTopol.getNumVertices() + i, 0, - exp_smoother(2*(0.5 - coords_initial[1]))*2*dw_dt*(coords_initial[0] - 0.5));
                }
            }

            Dest.setFromTriplets(triplets.begin(), triplets.end());
        }
};   

template<typename ConfiguratorType>
class FoldDofsCrossInterpolation : public FoldDofs<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    public:
        FoldDofsCrossInterpolation(const MeshTopologySaver &plateTopol,
                            const VectorType &plateGeomInitial,
                            const VectorType &plateGeomRef_basic,
                            const std::vector<int> &bdryMaskRef)
        : FoldDofs<ConfiguratorType>(plateTopol,plateGeomInitial,plateGeomRef_basic,bdryMaskRef)
        {
            this->initialize_folds_edges();
        }

        void apply(const VectorType &t, VectorType &Dest) const override{
            if(Dest.size() != 3*this->_plateTopol.getNumVertices()){
                Dest.resize(3*this->_plateTopol.getNumVertices());
            }

            Dest = this->_plateGeomRef_basic;

            // // Different parameters for the translation
            // w_0: controls all the outer vertices
            RealType w_0 = 0.2*std::tanh(t[0]);
            RealType w_1 = 0.2*std::tanh(t[1]);

            // Calculate the baseline polynomial interpolation in the interval x \in [0,0.5]
            RealType a = -16.0*w_1 + 8*w_0;
            RealType b = -2.0*w_0;

            // Now, the interpolation polynomial p(x) depends on the w's via a and b:
            auto p = [a, b](RealType x) -> RealType{
                return (x - 0.5)*(a*x + b);
            };

            // Now, apply the translation to the fold vertices
            
            for(int i = 0; i < this->_foldVertices.size(); i++)
            {
                VecType coords, coords_initial;
                getXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_initial, this->_foldVertices[i]);
                
                
                if((std::abs(coords_initial[1]-0.5) < 1e-4) && (coords_initial[0] <= 0.5))
                {
                    coords[1] = p(coords_initial[0]) + 0.5;
                }
                else if((std::abs(coords_initial[1]- 0.5) < 1e-4) && (coords_initial[0] >= 0.5))
                {
                    coords[1] = -p(1.0 - coords_initial[0]) + 0.5;
                }
                else if(std::abs(coords_initial[0]) - 0.5 < 1e-4 && coords_initial[1] <= 0.5)
                {
                    coords[0] = -p(coords_initial[1]) + 0.5;
                }
                else if(std::abs(coords_initial[0]) - 0.5 < 1e-4 && coords_initial[1] >= 0.5)
                {
                    coords[0] = p(1.0 - coords_initial[1]) + 0.5;
                }
                setXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
            }   

            // Now, use the DirichletSmoother to regularize
            //DirichletSmoother<DefaultConfigurator> smoother(this->_plateGeomInitial,this->_bdryMaskRef, this->_plateTopol);
            //smoother.apply(Dest, Dest);

            // Candidates for mesh smoothing functions:
            auto linear = [](RealType x) -> RealType{
                return 1 - x;
            };
            
            /*
            auto exp_smoother = [](RealType x) -> RealType{
                // Can vary the k
                RealType k = 2.0;
                return (1.0 - std::exp(-k*(1-x)))/(1.0 - std::exp(-k));
            };*/

            auto exp_smoother = [](RealType x) -> RealType{
                return (1.0 - x);
            };

            // Now, smooth all of this using interpolation
            for(int i = 0; i < this->_plateTopol.getNumVertices(); i++)
            {
                VecType coords, coords_initial;
                getXYZCoord<VectorType, VecType>( Dest, coords, i);
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_initial, i);

                // dont want to perturb the folds
                if(std::abs(coords_initial[0] - 0.5) < 1e-4 || std::abs(coords_initial[1] - 0.5) < 1e-4)
                {
                    continue;
                }
                
                // Check in which quadrant we are
                if(coords_initial[0] <= 0.5 && coords_initial[1] <= 0.5)
                {
                    coords[0] = coords_initial[0] - exp_smoother(2*(0.5 - coords_initial[0]))*p(coords_initial[1]);
                    coords[1] = coords_initial[1] + exp_smoother(2*(0.5 - coords_initial[1]))*p(coords_initial[0]);
                }
                else if(coords_initial[0] <= 0.5 && coords_initial[1] >= 0.5)
                {
                    coords[0] = coords_initial[0] + exp_smoother(2*(0.5 - coords_initial[0]))*p(1.0 - coords_initial[1]);
                    coords[1] = coords_initial[1] + exp_smoother(2*(coords_initial[1] - 0.5))*p(coords_initial[0]);
                }
                else if(coords_initial[0] >= 0.5 && coords_initial[1] >= 0.5)
                {
                    coords[0] = coords_initial[0] + exp_smoother(2*(coords_initial[0] - 0.5))*p(1.0 - coords_initial[1]);
                    coords[1] = coords_initial[1] - exp_smoother(2*(coords_initial[1] - 0.5))*p(1.0 - coords_initial[0]);
                }
                else if(coords_initial[0] >= 0.5 && coords_initial[1] <= 0.5)
                {
                    coords[0] = coords_initial[0] - exp_smoother(2*(coords_initial[0] - 0.5))*p(coords_initial[1]);
                    coords[1] = coords_initial[1] - exp_smoother(2*(0.5 - coords_initial[1]))*p(1.0 - coords_initial[0]);
                }
                setXYZCoord<VectorType, VecType>( Dest, coords, i);
            }
        }

        bool isFoldVertex(const RealType coord_x, const RealType coord_y) const
        {
            if((std::abs(coord_x - 0.5) < 1e-4) || (std::abs(coord_y - 0.5) < 1e-4)){
                return true;
            }
            return false;
        }

        bool isFoldEdge(const int edgeIdx) const
        {
            int node_i = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
            int node_j = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

            VecType coords_i, coords_j;
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_i, node_i);
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_j, node_j);

            if((std::abs(coords_i[0] - 0.5) < 1e-4 && std::abs(coords_j[0] - 0.5) < 1e-4) || (std::abs(coords_i[1] - 0.5) < 1e-4 && std::abs(coords_j[1] - 0.5) < 1e-4)){
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
            return 2;
        }
};

template<typename ConfiguratorType>
class FoldDofsCrossInterpolationGradient : public FoldDofsGradient<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    public:
        FoldDofsCrossInterpolationGradient(const MeshTopologySaver &plateTopol,
            const std::vector<int> &bdryMaskRef,
            const VectorType &plateGeomInitial,
            const std::vector<int> &foldVertices 
        ): FoldDofsGradient<ConfiguratorType>(plateTopol,bdryMaskRef,plateGeomInitial,foldVertices)
        {}

        void apply(const VectorType &t, MatrixType& Dest) const override{
           
            if(Dest.rows() != this->_plateGeomInitial.size() || Dest.cols() != 2){
                Dest.resize(this->_plateGeomInitial.size(), 2);
            }

            Dest.setZero();

            // derivative of the translation of the fold vertices
            // differentiate concatenation of rotation and scaling
            // -> first rotation, then scaling. both are linear operations
            // but at very first, differentiate how the angle changes with t
            // angle = atan(2*t) -> derivative 2/(1 + (2*t^2))

            RealType w_0 = 0.2*std::tanh(t[0]);
            RealType w_1 = 0.2*std::tanh(t[1]);
            RealType dw_0_dt = 0.2*(1.0 - std::tanh(t[0])*std::tanh(t[0]));
            RealType dw_1_dt = 0.2*(1.0 - std::tanh(t[1])*std::tanh(t[1]));
            FullMatrixType dw_dt(2,2);
            dw_dt(0,0) = dw_0_dt;
            dw_dt(0,1) = 0.0;
            dw_dt(1,0) = 0.0;
            dw_dt(1,1) = dw_1_dt;
            // Calculate the baseline polynomial interpolation in the interval x \in [0,0.5]
            RealType a = -16.0*w_1 + 8*w_0;
            RealType b = -2.0*w_0;

            FullMatrixType dab_dw(2,2);
            dab_dw(0,0) = 8.0;
            dab_dw(0,1) = -16.0;
            dab_dw(1,0) = -2.0;
            dab_dw(1,1) = 0.0;

            auto exp_smoother = [](RealType x) -> RealType{
                return (1.0 - x);
            };

            // We have different polynomials depending on the quadrant

            // ---------------------- P1 ---------------------------------------

            auto dp_dab = [](RealType x) -> FullMatrixType{
                FullMatrixType dp(1,2);
                dp(0,0) = x*x - 0.5*x;
                dp(0,1) = x - 0.5;
                return dp;
            };  

            auto dp_dt = [dp_dab, dw_dt, dab_dw](RealType x) -> FullMatrixType{
                FullMatrixType dp(1,2);
                dp = dp_dab(x)*dab_dw*dw_dt;
                return dp;
            }; 
            
            MatrixType rhs(3*this->_plateTopol.getNumVertices(), 2);
            std::vector<Eigen::Triplet<RealType>> triplets;

            for(int i = 0; i < this->_plateTopol.getNumVertices(); i++)
            {
                VecType coords, coords_initial;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_initial, i);
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords, i);

                if(coords_initial[0] <= 0.5 && coords_initial[1] <= 0.5)
                {
                    RealType dx_dt_1 = - exp_smoother(2*(0.5 - coords_initial[0]))*dp_dt(coords_initial[1])(0,0);
                    RealType dx_dt_2 = - exp_smoother(2*(0.5 - coords_initial[0]))*dp_dt(coords_initial[1])(0,1);
                    RealType dy_dt_1 = exp_smoother(2*(0.5 - coords_initial[1]))*dp_dt(coords_initial[0])(0,0);
                    RealType dy_dt_2 = exp_smoother(2*(0.5 - coords_initial[1]))*dp_dt(coords_initial[0])(0,1);
                    triplets.push_back(Eigen::Triplet<RealType>(i, 0, dx_dt_1));
                    triplets.push_back(Eigen::Triplet<RealType>(i, 1, dx_dt_2));
                    triplets.push_back(Eigen::Triplet<RealType>(i + this->_plateTopol.getNumVertices(), 0, dy_dt_1));
                    triplets.push_back(Eigen::Triplet<RealType>(i + this->_plateTopol.getNumVertices(), 1, dy_dt_2));
                }
                else if(coords_initial[0] <= 0.5 && coords_initial[1] >= 0.5)
                {
                    RealType dx_dt_1 = exp_smoother(2*(0.5 - coords_initial[0]))*dp_dt(1.0 - coords_initial[1])(0,0);
                    RealType dx_dt_2 = exp_smoother(2*(0.5 - coords_initial[0]))*dp_dt(1.0 - coords_initial[1])(0,1);
                    RealType dy_dt_1 = exp_smoother(2*(coords_initial[1] - 0.5))*dp_dt(coords_initial[0])(0,0);
                    RealType dy_dt_2 = exp_smoother(2*(coords_initial[1] - 0.5))*dp_dt(coords_initial[0])(0,1);
                    triplets.push_back(Eigen::Triplet<RealType>(i, 0, dx_dt_1));
                    triplets.push_back(Eigen::Triplet<RealType>(i, 1, dx_dt_2));
                    triplets.push_back(Eigen::Triplet<RealType>(i + this->_plateTopol.getNumVertices(), 0, dy_dt_1));
                    triplets.push_back(Eigen::Triplet<RealType>(i + this->_plateTopol.getNumVertices(), 1, dy_dt_2));
                }
                else if(coords_initial[0] >= 0.5 && coords_initial[1] >= 0.5)
                {
                    /*
                    if(i == 0)
                    {
                        std::cout<<"Stop here"<<std::endl;
                        std::cout<<dp_dt(1.0 - coords_initial[1])(0,0);
                        std::cout<<"Should be positive"<<std::endl;
                    }*/
                    RealType dx_dt_1 = exp_smoother(2*(coords_initial[0] - 0.5))*dp_dt(1.0 - coords_initial[1])(0,0);
                    RealType dx_dt_2 = exp_smoother(2*(coords_initial[0] - 0.5))*dp_dt(1.0 - coords_initial[1])(0,1);
                    RealType dy_dt_1 = - exp_smoother(2*(coords_initial[1] - 0.5))*dp_dt(1.0 - coords_initial[0])(0,0);
                    RealType dy_dt_2 = - exp_smoother(2*(coords_initial[1] - 0.5))*dp_dt(1.0 - coords_initial[0])(0,1);
                    triplets.push_back(Eigen::Triplet<RealType>(i, 0, dx_dt_1));
                    triplets.push_back(Eigen::Triplet<RealType>(i, 1, dx_dt_2));
                    triplets.push_back(Eigen::Triplet<RealType>(i + this->_plateTopol.getNumVertices(), 0, dy_dt_1));
                    triplets.push_back(Eigen::Triplet<RealType>(i + this->_plateTopol.getNumVertices(), 1, dy_dt_2));
                }
                else if(coords_initial[0] >= 0.5 && coords_initial[1] <= 0.5)
                {
                    RealType dx_dt_1 = - exp_smoother(2*(coords_initial[0] - 0.5))*dp_dt(coords_initial[1])(0,0);
                    RealType dx_dt_2 = - exp_smoother(2*(coords_initial[0] - 0.5))*dp_dt(coords_initial[1])(0,1);
                    RealType dy_dt_1 = - exp_smoother(2*(0.5 - coords_initial[1]))*dp_dt(1.0 - coords_initial[0])(0,0);
                    RealType dy_dt_2 = - exp_smoother(2*(0.5 - coords_initial[1]))*dp_dt(1.0 - coords_initial[0])(0,1);
                    triplets.push_back(Eigen::Triplet<RealType>(i, 0, dx_dt_1));
                    triplets.push_back(Eigen::Triplet<RealType>(i, 1, dx_dt_2));
                    triplets.push_back(Eigen::Triplet<RealType>(i + this->_plateTopol.getNumVertices(), 0, dy_dt_1));
                    triplets.push_back(Eigen::Triplet<RealType>(i + this->_plateTopol.getNumVertices(), 1, dy_dt_2));
                }
            }
            Dest.setFromTriplets(triplets.begin(), triplets.end());
        }
};   

template<typename ConfiguratorType>
class FoldDofsArcInterpolation : public FoldDofs<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    const size_t _num_dofs;
    const std::vector<int> _polDegrees;
    const std::vector<RealType> _x_points;

    Eigen::PartialPivLU<FullMatrixType> _lu;

    public:
        FoldDofsArcInterpolation(const MeshTopologySaver &plateTopol,
            const VectorType &plateGeomInitial,
            const VectorType &plateGeomRef_basic,
            const std::vector<int> &bdryMaskRef,
            std::vector<int> polDegrees,
            std::vector<RealType> x_points)
            : FoldDofs<ConfiguratorType>(plateTopol,plateGeomInitial,plateGeomRef_basic,bdryMaskRef),
            _num_dofs(polDegrees.size()),
            _polDegrees(polDegrees),
            _x_points(x_points)
        {
            this->initialize_folds_edges();

            // Initialize the inerpolationMatrix
            FullMatrixType interpolationMatrix;
            interpolationMatrix.resize(_num_dofs, _num_dofs);
            interpolationMatrix.setZero();

            assert(_x_points.size() == _num_dofs);

            for(size_t i = 0; i < _num_dofs; i++)
            {
                for(size_t j = 0; j < _num_dofs; j++)
                {
                    interpolationMatrix(i,j) = std::pow(_x_points[i] - 0.5, _polDegrees[j]);
                }
            }

            Eigen::PartialPivLU<FullMatrixType> lu(interpolationMatrix); 
            _lu = lu;
        }

        void apply(const VectorType &t, VectorType &Dest) const override{
            if(Dest.size() != 3*this->_plateTopol.getNumVertices()){
                Dest.resize(3*this->_plateTopol.getNumVertices());
            }

            assert(t.size() == _num_dofs);

            // // Different parameters for the translation
            // w_0: controls all the outer vertices
            VectorType w(_num_dofs);
            for(size_t i = 0; i < _num_dofs; i++)
            {
                //w[i] = 0.25 + 0.5*std::tanh(t[i]);
                w[i] = 0.45 + 0.35*std::tanh(t[i]);
            }

            Dest = this->_plateGeomRef_basic;

            VectorType coeffs = _lu.solve(w);

            auto p = [coeffs, this](RealType x) -> RealType{
                RealType y_val = 0.0;
                for(int i = 0; i < coeffs.size(); i++)
                {
                    y_val += coeffs[i] * std::pow(x-0.5, this->_polDegrees[i]);
                }
                return y_val;
            };

            // Now, apply the translation to the fold vertices
            for(int i = 0; i < this->_foldVertices.size(); i++)
            {
                VecType coords, coords_initial;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_initial, this->_foldVertices[i]);
                getXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);

                coords[1] = p(coords_initial[0]);
                setXYZCoord<VectorType, VecType>( Dest, coords, this->_foldVertices[i]);
            }

            DirichletSmoother<ConfiguratorType> smoother(this->_plateGeomInitial, this->_bdryMaskRef, this->_plateTopol);
            smoother.apply(Dest, Dest);
        }

        bool isFoldVertex(const RealType coord_x, const RealType coord_y) const
        {
            if((std::abs(coord_y - 0.25) < 1e-4)){
                return true;
            }
            return false;
        }

        bool isFoldEdge(const int edgeIdx) const
        {
            int node_i = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
            int node_j = this->_plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

            VecType coords_i, coords_j;
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_i, node_i);
            getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_j, node_j);

            if(isFoldVertex(coords_i[0], coords_i[1]) && isFoldVertex(coords_j[0], coords_j[1])){
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
            return _num_dofs;
        }
};

template<typename ConfiguratorType>
class FoldDofsArcInterpolationGradient : public FoldDofsGradient<ConfiguratorType>, public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;

        const std::vector<int> _polDegrees;
        const std::vector<RealType> _x_points;
        const size_t _num_dofs;

        Eigen::PartialPivLU<FullMatrixType> _lu;

        FullMatrixType _interpolationMatrix;
        FullMatrixType _interpolationMatrix_inv;

    public:
        FoldDofsArcInterpolationGradient(const MeshTopologySaver &plateTopol,
            const std::vector<int> &bdryMaskRef,
            const VectorType &plateGeomInitial,
            const std::vector<int> &foldVertices,
            std::vector<int> polDegrees,
            std::vector<RealType> x_points 
        ): FoldDofsGradient<ConfiguratorType>(plateTopol,bdryMaskRef,plateGeomInitial,foldVertices),
           _polDegrees(polDegrees),
           _x_points(x_points),
           _num_dofs(polDegrees.size())
        {
            assert(_x_points.size() == _num_dofs);

            // Initialize the inerpolationMatrix
            FullMatrixType interpolationMatrix;
            interpolationMatrix.resize(_num_dofs, _num_dofs);
            interpolationMatrix.setZero();

            for(size_t i = 0; i < _num_dofs; i++)
            {
                for(size_t j = 0; j < _num_dofs; j++)
                {
                    interpolationMatrix(i,j) = std::pow(_x_points[i] - 0.5, _polDegrees[j]);
                }
            }

            FullMatrixType interpolationMatrix_inv(_num_dofs, _num_dofs);
            interpolationMatrix_inv.setZero();
            interpolationMatrix_inv = interpolationMatrix.fullPivLu().inverse();
            _interpolationMatrix_inv = interpolationMatrix_inv;

            Eigen::PartialPivLU<FullMatrixType> lu(interpolationMatrix);
            _lu = lu;

        }

        void apply(const VectorType &t, MatrixType& Dest) const override{

            if(Dest.rows() != 3*this->_plateTopol.getNumVertices() || Dest.cols() != _num_dofs){
                Dest.resize(3*this->_plateTopol.getNumVertices(), _num_dofs);
            }

            Dest.setZero();
            assert(t.size() == _num_dofs);
            FullMatrixType dw_dt(_num_dofs, _num_dofs);
            dw_dt.setZero();
            for(int i = 0; i < _num_dofs; i++)
            {
                dw_dt(i,i) = 0.35*(1.0 - std::tanh(t[i])*std::tanh(t[i]));
            }

            FullMatrixType da_dt(_num_dofs, _num_dofs);
            da_dt.setZero();
            for(int i = 0; i < _num_dofs; i++)
            {
                VectorType row_i = _lu.solve(dw_dt.col(i));
                da_dt.row(i) = row_i.transpose();
            }

            FullMatrixType dp_da(this->_plateTopol.getNumVertices(), _num_dofs);
            dp_da.setZero();
            std::vector<Eigen::Triplet<RealType>> triplets;
            for(int i = 0; i < this->_foldVertices.size(); i++)
            {
                VecType coords_initial;
                getXYZCoord<VectorType, VecType>( this->_plateGeomInitial, coords_initial, this->_foldVertices[i]);

                for(int j = 0; j < _num_dofs; j++)
                {
                    dp_da(this->_foldVertices[i],j) = std::pow(coords_initial[0] - 0.5, _polDegrees[j]);
                }
            }

            MatrixType dp_dt_block = (dp_da * da_dt.transpose()).sparseView();
            MatrixType dp_dt(3*this->_plateTopol.getNumVertices(), _num_dofs);
            dp_dt.setZero();
            std::vector<Eigen::Triplet<RealType>> triplets_dp_dt;
            assignSparseBlockInplace(dp_dt, dp_dt_block, this->_plateTopol.getNumVertices(), 0, triplets_dp_dt);

            typename ConfiguratorType::SparseMatrixType StiffnessMatrix;
            computeStiffnessMatrix<ConfiguratorType>(this->_plateTopol,this->_plateGeomInitial, StiffnessMatrix);
            applyMaskToMajor<typename ConfiguratorType::SparseMatrixType>( this->_bdryMaskRef, StiffnessMatrix );

            Eigen::BiCGSTAB<MatrixType> solver;
            solver.compute(StiffnessMatrix);

            std::vector<Eigen::Triplet<RealType>> tripletList;
            // solve the equation L*dest = indicator_dof for each column of indicator_dof
            VectorType dest(3*this->_plateTopol.getNumVertices());
            for(int i = 0; i < _num_dofs; i++)
            {
                dest = solver.solve(dp_dt.col(i));
                assignSparseBlockInplace(Dest, convertVecToSparseMat(dest), 0, i, tripletList);
                dest.setZero();
            }

            //assignSparseBlockInplace(Dest, dp_dt, this->_plateTopol.getNumVertices(), 0, triplets);
        }
};   