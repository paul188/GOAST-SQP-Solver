#pragma once

#include <goast/QuadMesh/QuadTopology.h>
#include <cmath>
#include <goast/SQP/Utils/SparseMat.h>
#include <optional>

class SkippingBdryHalfEdgeIterator{
    private:
        const QuadMeshTopologySaver &_quadTopol;
        const std::vector<int> &_bdryHalfEdges;
        const size_t _num_halfedges;
        std::set<int> _noBdryHalfEdges;
        std::map<int,int> _heh_idx_to_nbdry;
        std::set<int>::iterator _it;
    public:
        SkippingBdryHalfEdgeIterator(const QuadMeshTopologySaver &quadTopol)
            : _quadTopol(quadTopol),
              _bdryHalfEdges(quadTopol.getBdryHalfEdges()),
              _num_halfedges(quadTopol.getNumHalfEdges())
        {
            int bdryCounter = 0;
            for(int i = 0; i < _num_halfedges; i++){
                if(i == _bdryHalfEdges[bdryCounter]){
                    _heh_idx_to_nbdry[i] = -1;
                    bdryCounter++;
                }
                else{
                    _heh_idx_to_nbdry[i] = i - bdryCounter;
                    _noBdryHalfEdges.insert(i);
                }
            }

            _it = _noBdryHalfEdges.begin();
        }

        void operator++(){
            _it++;
        }

        // Postfix increment operator
        SkippingBdryHalfEdgeIterator operator++(int) {
            SkippingBdryHalfEdgeIterator temp = *this; // Save the current state
            ++_it; // Use the prefix increment to advance the iterator
            return temp; // Return the saved state
        }

        int idx(){
            return *_it;
        }

        int idx_nobdry(){
            return std::distance(_noBdryHalfEdges.begin(), _it);
        }

        bool valid(){
            return _it != _noBdryHalfEdges.end();
        }

        int to_nbdry(int idx){
            return _heh_idx_to_nbdry[idx];
        }
        
        // reset the iterator
        void reset()
        {
            _it = _noBdryHalfEdges.begin();
        }
};

class SkippingBdryFaceIterator{
    private:
        const QuadMeshTopologySaver &_quadTopol;
        const std::vector<int> &_bdryFaces;
        const size_t _num_faces;
        std::set<int> _noBdryFaces;
        std::map<int,int> _face_idx_to_nbdry;
        std::set<int>::iterator _it;

    public:
        SkippingBdryFaceIterator(const QuadMeshTopologySaver &quadTopol)
            : _quadTopol(quadTopol),
              _bdryFaces(quadTopol.getBdryFaces()),
              _num_faces(quadTopol.getNumFaces())
        {
            int bdryCounter = 0;
            for(int i = 0; i < _num_faces; i++){
                if(i == _bdryFaces[bdryCounter]){
                    _face_idx_to_nbdry[i] = -1;
                    bdryCounter++;
                }
                else{
                    _face_idx_to_nbdry[i] = i - bdryCounter;
                    _noBdryFaces.insert(i);
                }
            }

            _it = _noBdryFaces.begin();
        }

        void operator++(){
            _it++;
        }

        // Postfix increment operator
        SkippingBdryFaceIterator operator++(int) {
            SkippingBdryFaceIterator temp = *this; // Save the current state
            ++_it; // Use the prefix increment to advance the iterator
            return temp; // Return the saved state
        }

        int idx(){
            return *_it;
        }

        int idx_nobdry(){
            return std::distance(_noBdryFaces.begin(), _it);
        }

        bool valid(){
            return (_it != _noBdryFaces.end());
        }

        int to_nbdry(int idx){
            return _face_idx_to_nbdry[idx];
        }

        // reset the iterator
        void reset()
        {
            _it = _noBdryFaces.begin();
        }
};

/* Store the vertex, halfedge and face strips */
template<typename ConfiguratorType>
class StripHandler{
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    
        size_t _num_vertexStrips;
        size_t _num_faceStrips;
        size_t _num_halfedgeStrips;
        std::vector<std::tuple<int,int,int>> _vertex_strips;
        std::vector<std::tuple<int,int,int>> _face_strips;
        std::vector<std::tuple<int,int,int>> _halfedge_strips;

    public:
        StripHandler(QuadMeshTopologySaver &quadTopol)
        {
            _vertex_strips = quadTopol.getVertexStrips();
            _num_vertexStrips = _vertex_strips.size();
            _halfedge_strips = quadTopol.getHalfEdgeStrips();
            _num_halfedgeStrips = _halfedge_strips.size();
            _face_strips = quadTopol.getFaceStrips();
            _num_faceStrips = _face_strips.size();
        }

        std::tuple<int,int,int> vertex_strip(size_t i){
            if(i >= _num_vertexStrips){
                std::cerr << "Index out of bounds for vertex strips"<<std::endl;
            }
            return _vertex_strips[i];
        }

        std::tuple<int,int,int> face_strip(size_t i){
            if(i >= _num_faceStrips){
                std::cerr << "Index out of bounds for face strips"<<std::endl;
            }
            return _face_strips[i];
        }

        std::tuple<int,int,int> halfedge_strip(size_t i){
            if(i >= _num_halfedgeStrips){
                std::cerr << "Index out of bounds for halfedge strips"<<std::endl;
            }
            return _halfedge_strips[i];
        }

        size_t num_vertexStrips() const{
            return _num_vertexStrips;
        }

        size_t num_faceStrips() const{
            return _num_faceStrips;
        }

        size_t num_halfedgeStrips() const{
            return _num_halfedgeStrips;
        }

};

 /* storage class for indices of the variables */
template<typename ConfiguratorType>
class VarsIdx{
     public:
         typedef typename ConfiguratorType::RealType RealType;
         typedef typename ConfiguratorType::VectorType VectorType;
         typedef typename ConfiguratorType::SparseMatrixType MatrixType;
 
     private:
         const size_t _num_vertices;
         const size_t _num_faces;
         const size_t _num_edges;
         const size_t _num_halfedges;
         const size_t _num_bdryHalfEdges;
         const size_t _num_bdryFaces;
         const QuadMeshTopologySaver &_quadTopol;
 
     private:
         std::map<std::string, int> _idx;
 
         /*
         Ordering of the dofs: 
         1. dummy_weights from 0 to num_vertices-1
         2. vertices from num_vertices to 4*num_vertices-1 
         3. weights from 4*num_vertices to 5*num_vertices-1
         4. reweighted_vertices from 5*num_vertices to 8*num_vertices-1
         5. reweighted_edges_1 from 8*num_vertices to (8*num_vertices + 3*num_faces)-1
         6. reweighted_edges_2 from (8*num_vertices + 3*num_faces) to (8*num_vertices + 6*num_faces)-1
         7. normals from (8*num_vertices + 6*num_faces) to (8*num_vertices + 9*num_faces)-1
         8. rulings from (8*num_vertices + 9*num_faces) to (8*num_vertices + 9*num_faces + 3*(num_halfedges - bdryHalfEdges.size()))-1
         9. reweighted_rulings_1 from (8*num_vertices + 9*num_faces + 3*(num_halfedges - bdryHalfEdges.size()) to (8*num_vertices + 9*num_faces + 3*(num_faces - bdryFaces.size()) + 3*(num_halfedges - bdryHalfEdges.size()))-1
         10.reweighted_rulings_2 from (8*num_vertices + 9*num_faces + 3*(num_faces - bdryFaces.size()) + 3*(num_halfedges - bdryHalfEdges.size())) to (8*num_vertices + 9*num_faces + 9*num_faces + 6*(num_faces - bdryFaces.size()) + 3*(num_halfedges- bdryHalfEdges.size()))-1
         */
 
        // We purposely do not include rulings at bdry halfedges 
        // and reweighted rulings at bdry adjacent faces. They cannot be used to express the developability constraint
        // and they are not determined by the constraints and further do not play a role in the fairness energy
 
     public:
         VarsIdx(const QuadMeshTopologySaver &quadTopol)
             : _num_vertices(quadTopol.getNumVertices()),
               _num_faces(quadTopol.getNumFaces()),
               _num_edges(quadTopol.getNumEdges()),
               _num_halfedges(quadTopol.getNumHalfEdges()),
               _num_bdryHalfEdges(quadTopol.getBdryHalfEdges().size()),
               _num_bdryFaces(quadTopol.getBdryFaces().size()),                  
               _quadTopol(quadTopol)
 
         {
             _idx["dummy_weights"] = 0;
             _idx["vertices"] = _num_vertices;
             _idx["weights"] = 4*_num_vertices;
             _idx["reweighted_vertices"] = 5*_num_vertices;
             _idx["reweighted_edges_1"] = 8*_num_vertices;
             _idx["reweighted_edges_2"] = 8*_num_vertices + 3*_num_faces;
             _idx["normals"] = 8*_num_vertices + 6*_num_faces;
             _idx["rulings"] = 8*_num_vertices + 9*_num_faces;
             _idx["reweighted_rulings_1"] = 8*_num_vertices + 9*_num_faces + 3*(_num_halfedges - _num_bdryHalfEdges);
             _idx["reweighted_rulings_2"] = 8*_num_vertices + 12*_num_faces - 3*_num_bdryFaces + 3*(_num_halfedges - _num_bdryHalfEdges);
             _idx["num_dofs"] = 8*_num_vertices + 15*_num_faces - 6*_num_bdryFaces + 3*(_num_halfedges - _num_bdryHalfEdges);
         
         }
 
         RealType dummy_weight(const VectorType &vars, size_t i) const {
             if(i >= _num_vertices){
                 std::cerr << "Index out of bounds for dummy weights"<<std::endl;
             }
             return vars[_idx.at("dummy_weights") + i];
         }
 
         VectorType vertex(const VectorType &vars, size_t i) const {
             if(i >= _num_vertices){
                 std::cerr << "Index out of bounds for vertices"<<std::endl;
             }
             return vars.segment(_idx.at("vertices") + 3*i,3);
         }
 
         RealType vertex_weight(const VectorType &vars, size_t i) const{
             if(i >= _num_vertices){
                 std::cerr << "Index out of bounds for vertex weights"<<std::endl;
             }
             return vars[_idx.at("weights") + i];
         }
 
         VectorType reweighted_vertex(const VectorType &vars, size_t i) const{
             if(i >= _num_vertices){
                 std::cerr << "Index out of bounds for reweighted vertices"<<std::endl;
             }
             return vars.segment(_idx.at("reweighted_vertices") + 3*i,3);
         }
 
         VectorType reweighted_edge_1(const VectorType &vars, size_t i) const{
             if(i >= _num_faces){
                 std::cerr << "Index out of bounds for reweighted edges 1"<<std::endl;
             }
             return vars.segment(_idx.at("reweighted_edges_1") + 3*i, 3);
         }
 
         VectorType reweighted_edge_2(const VectorType &vars, size_t i) const{
             if(i >= _num_faces){
                 std::cerr << "Index out of bounds for reweighted edges 2"<<std::endl;
             }
             return vars.segment(_idx.at("reweighted_edges_2") + 3*i, 3);
         }
 
         VectorType face_normal(const VectorType &vars, size_t i) const{
             if(i >= _num_faces){
                 std::cerr << "Index out of bounds for face normals"<<std::endl;
             }
             return vars.segment(_idx.at("normals") + 3*i, 3);
         }
 
         VectorType ruling(const VectorType &vars, size_t i) const{
             if(i >= _num_halfedges){
                 std::cerr << "Index out of bounds for rulings"<<std::endl;
             }
             return vars.segment(_idx.at("rulings") + 3*i,3);
         }
 
         VectorType reweighted_ruling_1(const VectorType &vars, size_t i) const{
             if(i >= _num_faces){
                 std::cerr << "Index out of bounds for reweighted rulings 1"<<std::endl;
             }
             return vars.segment(_idx.at("reweighted_rulings_1") + 3*i, 3);
         }
 
         VectorType reweighted_ruling_2(const VectorType &vars, size_t i) const{
             if(i >= _num_faces){
                 std::cerr << "Index out of bounds for reweighted rulings 1"<<std::endl;
             }
             return vars.segment(_idx.at("reweighted_rulings_2") + 3*i, 3);
         }
 
         size_t operator[](const std::string& key) const {
             return _idx.at(key);
         }
 };
 
 /* storage class for indices of constraints */
template<typename ConfiguratorType>
class ConstraintIdx{
     public:
         typedef typename ConfiguratorType::RealType RealType;
         typedef typename ConfiguratorType::VectorType VectorType;
         typedef typename ConfiguratorType::SparseMatrixType MatrixType;
         // used to store the boundary indices and corresponding values
         typedef typename std::optional<std::pair<std::vector<int>, VectorType>> OptionalBoundaryData;
 
     private:
         std::map<std::string,size_t> _cons_idx;
     public:
         ConstraintIdx(StripHandler<ConfiguratorType> &stripHandle, const QuadMeshTopologySaver &quadTopol)
         {
 
             // index handling
             size_t num_vertices = quadTopol.getNumVertices();
             size_t num_faces = quadTopol.getNumFaces();
             size_t num_halfedges = 2*quadTopol.getNumEdges();
             size_t num_bdryHalfEdges = quadTopol.getBdryHalfEdges().size();
             size_t num_bdryFaces = quadTopol.getBdryFaces().size();
             size_t num_faceStrips = stripHandle.num_faceStrips();
             size_t num_vertexStrips = stripHandle.num_vertexStrips();
             size_t num_halfedgeStrips = stripHandle.num_halfedgeStrips();
 
             size_t num_constraints_vert_1 = num_vertices;
             size_t num_constraints_vert_2 = 3*num_vertices;
             size_t num_constraints_edge_vec_1 = 3*num_faces;
             size_t num_constraints_edge_vec_2 = 3*num_faces;
             size_t num_constraints_normal_1 = num_faces;
             size_t num_constraints_normal_2 = num_faces;
             size_t num_constraints_normal_3 = num_faces;
             size_t num_constraints_ruling_0 = 3*(num_halfedges - num_bdryHalfEdges);
             size_t num_constraints_ruling_1 = 3*(num_faces - num_bdryFaces);
             size_t num_constraints_ruling_2 = 3*(num_faces - num_bdryFaces);
             size_t num_constraints_dev = 3*(num_faces - num_bdryFaces);
             size_t num_constraints_fair_v =  3*num_vertexStrips;
             size_t num_constraints_fair_n =  3*num_faceStrips;
             size_t num_constraints_fair_r =  3*num_halfedgeStrips;
 
             _cons_idx["vertex_1"] = 0;
             _cons_idx["vertex_2"] = _cons_idx["vertex_1"] + num_constraints_vert_1;
             _cons_idx["edge_vec_1"] = _cons_idx["vertex_2"] + num_constraints_vert_2;
             _cons_idx["edge_vec_2"] = _cons_idx["edge_vec_1"] + num_constraints_edge_vec_1;
             _cons_idx["normal_1"] = _cons_idx["edge_vec_2"] + num_constraints_edge_vec_2;
             _cons_idx["normal_2"] = _cons_idx["normal_1"] + num_constraints_normal_1;
             _cons_idx["normal_3"] = _cons_idx["normal_2"] + num_constraints_normal_2;
             _cons_idx["ruling_0"] = _cons_idx["normal_3"] + num_constraints_normal_3;
             _cons_idx["ruling_1"] = _cons_idx["ruling_0"] + num_constraints_ruling_0;
             _cons_idx["ruling_2"] = _cons_idx["ruling_1"] + num_constraints_ruling_1;
             _cons_idx["dev"] = _cons_idx["ruling_2"] + num_constraints_ruling_2;
             _cons_idx["fair_v"] = _cons_idx["dev"] + num_constraints_dev;
             _cons_idx["fair_n"] = _cons_idx["fair_v"] + num_constraints_fair_v;
             _cons_idx["fair_r"] = _cons_idx["fair_n"] + num_constraints_fair_n;
             _cons_idx["num_cons"] = _cons_idx["fair_r"] + num_constraints_fair_r;
         }
 
         size_t operator[](const std::string& key) const {
             if (key == "vertex_1"||
                 key == "vertex_2"||
                 key == "edge_vec_1"||
                 key == "edge_vec_2"||
                 key == "normal_1"||
                 key == "normal_2"||
                 key == "normal_3"||
                 key == "ruling_0"||
                 key == "ruling_1"||
                 key == "ruling_2"||
                 key == "dev"||
                 key == "fair_v"||
                 key == "fair_n"||
                 key == "fair_r"||
                 key == "num_cons"){
                    return _cons_idx.at(key);
                 }
             throw std::runtime_error("Unknown constraint key: " + key);
         }
 };
 
template<typename ConfiguratorType>
struct constraint_weights{
    typedef typename ConfiguratorType::RealType RealType;

    RealType vertex_1 = 1.0;
    RealType vertex_2 = 1.0;
    RealType edge_vec_1 = 1.0;
    RealType edge_vec_2 = 1.0;
    RealType normal_1 = 1.0;
    RealType normal_2 = 1.0;
    RealType normal_3 = 1.0;
    RealType ruling_0 = 1.0;
    RealType ruling_1 = 1.0;
    RealType ruling_2 = 1.0;
    RealType dev = 1.0;
    RealType fair_v = 1.0;
    RealType fair_n = 1.0;
    RealType fair_r = 1.0;

    RealType operator[](const std::string& key) const {
        if (key == "vertex_1") return vertex_1;
        else if (key == "vertex_2") return vertex_2;
        else if (key == "edge_vec_1") return edge_vec_1;
        else if (key == "edge_vec_2") return edge_vec_2;
        else if (key == "normal_1") return normal_1;
        else if (key == "normal_2") return normal_2;
        else if (key == "normal_3") return normal_3;
        else if (key == "ruling_0") return ruling_0;
        else if (key == "ruling_1") return ruling_1;
        else if (key == "ruling_2") return ruling_2;
        else if (key == "dev") return dev;
        else if (key == "fair_v") return fair_v;
        else if (key == "fair_n") return fair_n;
        else if (key == "fair_r") return fair_r;
        
        throw std::runtime_error("Unknown constraint key: " + key);
    }

    void extend_weights(const ConstraintIdx<ConfiguratorType> &constraintIdx,MatrixType &Dest)
    {
        size_t num_constraints = constraintIdx["num_cons"];
        if(Dest.rows() != num_constraints || Dest.cols() != num_constraints)
        {
            Dest.resize(num_constraints,num_constraints);
        }

        std::vector<Eigen::Triplet<RealType>> triplets;
        triplets.reserve(num_constraints);

        std::vector<std::pair<std::string, std::string>> constraints = {
            {"vertex_1", "vertex_2"}, {"vertex_2", "edge_vec_1"},
            {"edge_vec_1", "edge_vec_2"}, {"edge_vec_2", "normal_1"}, {"normal_1", "normal_2"},
            {"normal_2", "normal_3"}, {"normal_3", "ruling_0"}, {"ruling_0", "ruling_1"},
            {"ruling_1", "ruling_2"}, {"ruling_2", "dev"}, {"dev", "fair_v"},
            {"fair_v", "fair_n"}, {"fair_n", "fair_r"}, {"fair_r", "num_cons"}
        };

        for (const auto& [start, end] : constraints) {
            for (int i = constraintIdx[start]; i < constraintIdx[end]; i++) {
                triplets.emplace_back(i, i, this->operator[](start));  
            }
        }

        Dest.setFromTriplets(triplets.begin(), triplets.end());
    }
};

template <typename ConfiguratorType>
class Constraint : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    private:
        const QuadMeshTopologySaver &_quadTopol;
        size_t _num_vertices;
        size_t _num_faces;
        size_t _num_edges;
        size_t _num_halfedges;
        mutable SkippingBdryFaceIterator _it_face;
        mutable SkippingBdryHalfEdgeIterator _it_he;
        mutable StripHandler<ConfiguratorType> _stripHandle;

        // Serves to access the current vector of input variables at the right positions
        VarsIdx<ConfiguratorType> _vars_idx;
        // Serves to write the constraints into the vector
        // as well as the constraint gradients into the matrix
        ConstraintIdx<ConfiguratorType> _cons_idx;

    public:
        Constraint(const QuadMeshTopologySaver &quadTopol, 
                   StripHandler<ConfiguratorType> &stripHandle,
                   ConstraintIdx<ConfiguratorType> &constraintIdx,
                   VarsIdx<ConfiguratorType> &varsIdx):
              _quadTopol(quadTopol),
              _num_vertices(quadTopol.getNumVertices()),
              _num_faces(quadTopol.getNumFaces()),
              _num_edges(quadTopol.getNumEdges()),
              _num_halfedges(2*_num_edges),
              _stripHandle(stripHandle),
              _it_face(_quadTopol),
              _it_he(_quadTopol),
              _vars_idx(varsIdx),
              _cons_idx(constraintIdx)
        {
        }

    void initialize_vars(const VectorType &geom, VectorType &Dest) const {
        if(Dest.size() != _vars_idx["num_dofs"]){
            Dest.resize(_vars_idx["num_dofs"]);
        }

        Dest.setZero();

        // 1. dummy weights stay initialized at zero
        // 2. vertices are initialized with the input geometry
        Dest.segment(_vars_idx["vertices"],3*_num_vertices) = geom;
        // 3. weights are initialized with 1
        Dest.segment(_vars_idx["weights"],_num_vertices) = VectorType::Ones(_num_vertices);
        // 4. reweighted vertices are initialized with the input geometry
        Dest.segment(_vars_idx["reweighted_vertices"],3*_num_vertices) = geom;
        // 5. initialize reweighted edges:
        for(int i = 0; i < _num_faces; i++){
            // First, obtain first vertex idces of the face
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w_0 = _vars_idx.vertex_weight(Dest, node_0);
            RealType w_1 = _vars_idx.vertex_weight(Dest, node_1);
            RealType w_2 = _vars_idx.vertex_weight(Dest, node_2);
            RealType w_3 = _vars_idx.vertex_weight(Dest, node_3);

            Dest.segment(_vars_idx["reweighted_edges_1"]+3*i,3) = (w_0 + w_1)*(_vars_idx.reweighted_vertex(Dest,node_2)+_vars_idx.reweighted_vertex(Dest, node_3)) - (w_2 + w_3)*(_vars_idx.reweighted_vertex(Dest, node_1) + _vars_idx.reweighted_vertex(Dest, node_0));
            Dest.segment(_vars_idx["reweighted_edges_2"]+3*i,3) = (w_1 + w_2)*(_vars_idx.reweighted_vertex(Dest,node_3)+_vars_idx.reweighted_vertex(Dest, node_0)) - (w_0 + w_3)*(_vars_idx.reweighted_vertex(Dest, node_1) + _vars_idx.reweighted_vertex(Dest, node_2));
        }
        // 6. initialize normals
        for(int i = 0; i < _num_faces; i++){
            VectorType normal = (_vars_idx.reweighted_edge_1(Dest, i).template head<3>()).cross(_vars_idx.reweighted_edge_2(Dest, i).template head<3>());
            Dest.segment(_vars_idx["normals"]+3*i,3) = 1.0/(normal.norm())*normal;
        }

        _it_he.reset();
        // 7. initialize rulings
        for(_it_he; _it_he.valid(); ++_it_he){
            int i = _it_he.idx();
            int i_nobdry = _it_he.idx_nobdry();

            int face_1 = _quadTopol.getFaceOfHalfEdge(i,0);
            int face_2 = _quadTopol.getFaceOfHalfEdge(i,1);

            Dest.segment(_vars_idx["rulings"] + 3*i_nobdry,3) = (_vars_idx.face_normal(Dest, face_1).template head<3>()).cross(_vars_idx.face_normal(Dest, face_2).template head<3>()); 
        }

        //8. Initialize reweighted rulings
        _it_he.reset();
        _it_face.reset();

        for(_it_face; _it_face.valid(); _it_face++){

            int i = _it_face.idx();
            int i_nobdry = _it_face.idx_nobdry();

            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
            int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
            int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
            int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

            RealType w0 = _vars_idx.vertex_weight(Dest, node_0);
            RealType w1 = _vars_idx.vertex_weight(Dest, node_1);
            RealType w2 = _vars_idx.vertex_weight(Dest, node_2);
            RealType w3 = _vars_idx.vertex_weight(Dest, node_3);

            VectorType ruling_01 = _vars_idx.ruling(Dest, _it_he.to_nbdry(heh_01));
            VectorType ruling_12 = _vars_idx.ruling(Dest, _it_he.to_nbdry(heh_12));
            VectorType ruling_23 = _vars_idx.ruling(Dest, _it_he.to_nbdry(heh_23));
            VectorType ruling_30 = _vars_idx.ruling(Dest, _it_he.to_nbdry(heh_30));

            Dest.segment(_vars_idx["reweighted_rulings_1"] + 3*i_nobdry,3) = (w1 + w2)*ruling_12 - (w3 + w0)*(ruling_30);
            Dest.segment(_vars_idx["reweighted_rulings_2"] + 3*i_nobdry,3) = (w0 + w1)*ruling_01 - (w2 + w3)*(ruling_23);
        }
    }
    
    void get_constraint_vertex_1(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _num_vertices; i++){
            Dest[_cons_idx["vertex_1"] + i] = _vars_idx.vertex_weight(vars, i) - (_vars_idx.dummy_weight(vars, i)*_vars_idx.dummy_weight(vars, i)) - 1.0;
        }
    }

    void get_constraint_vertex_2(const VectorType &vars, VectorType &Dest) const{

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _num_vertices; i++){
            Dest.segment(_cons_idx["vertex_2"] + 3*i, 3) = _vars_idx.reweighted_vertex(vars, i) - _vars_idx.vertex_weight(vars, i)*_vars_idx.vertex(vars, i);
        }

    }

    void get_constraint_edge_vec_1(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _num_faces; i++){
            // First, obtain first vertex idces of the face
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w_0 = _vars_idx.vertex_weight(vars, node_0);
            RealType w_1 = _vars_idx.vertex_weight(vars, node_1);
            RealType w_2 = _vars_idx.vertex_weight(vars, node_2);
            RealType w_3 = _vars_idx.vertex_weight(vars, node_3);

            VectorType calc_edge_vec_1 = (w_0 + w_1)*(_vars_idx.reweighted_vertex(vars, node_2)+_vars_idx.reweighted_vertex(vars, node_3)) - (w_2 + w_3)*(_vars_idx.reweighted_vertex(vars, node_1) + _vars_idx.reweighted_vertex(vars, node_0));

            Dest.segment(_cons_idx["edge_vec_1"]+3*i,3) = _vars_idx.reweighted_edge_1(vars, i) - calc_edge_vec_1;
        }
    }

    void get_constraint_edge_vec_2(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _num_faces; i++){
            // First, obtain first vertex idces of the face
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w_0 = _vars_idx.vertex_weight(vars, node_0);
            RealType w_1 = _vars_idx.vertex_weight(vars, node_1);
            RealType w_2 = _vars_idx.vertex_weight(vars, node_2);
            RealType w_3 = _vars_idx.vertex_weight(vars, node_3);

            VectorType calc_edge_vec_2 = (w_1 + w_2)*(_vars_idx.reweighted_vertex(vars, node_3)+_vars_idx.reweighted_vertex(vars, node_0)) - (w_0 + w_3)*(_vars_idx.reweighted_vertex(vars, node_1) + _vars_idx.reweighted_vertex(vars, node_2));

            Dest.segment(_cons_idx["edge_vec_2"] + 3*i,3) = _vars_idx.reweighted_edge_2(vars, i) - calc_edge_vec_2;
        }
    }

    void get_constraint_normal_1(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _num_faces; i++){
            Dest[_cons_idx["normal_1"] + i] = _vars_idx.face_normal(vars, i).dot(_vars_idx.reweighted_edge_1(vars, i));
        }
    }

    void get_constraint_normal_2(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _num_faces; i++){
            Dest[_cons_idx["normal_2"] + i] = _vars_idx.face_normal(vars, i).dot(_vars_idx.reweighted_edge_2(vars, i));
        }
    }

    void get_constraint_normal_3(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _num_faces; i++){
            Dest[_cons_idx["normal_3"] + i] = _vars_idx.face_normal(vars,i).dot(_vars_idx.face_normal(vars,i)) - 1.0;
        }
    }

    void get_constraint_ruling_0(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _it_he.reset();

        for(_it_he; _it_he.valid(); ++_it_he){
            int i = _it_he.idx();
            int i_nobdry = _it_he.idx_nobdry();

            int face_1 = _quadTopol.getFaceOfHalfEdge(i,0);
            int face_2 = _quadTopol.getFaceOfHalfEdge(i,1);

            Dest.segment(_cons_idx["ruling_0"] + 3*i_nobdry,3) = _vars_idx.ruling(vars, i_nobdry) - (_vars_idx.face_normal(vars, face_1).template head<3>()).cross(_vars_idx.face_normal(vars, face_2).template head<3>()); 
        }

    }

    void get_constraint_ruling_1(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _it_he.reset();
        _it_face.reset();

        for(_it_face; _it_face.valid(); _it_face++){

            int i = _it_face.idx();
            int i_nobdry = _it_face.idx_nobdry();

            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
            int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
            int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
            int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

            RealType w0 = _vars_idx.vertex_weight(vars,node_0);
            RealType w1 = _vars_idx.vertex_weight(vars, node_1);
            RealType w2 = _vars_idx.vertex_weight(vars, node_2);
            RealType w3 = _vars_idx.vertex_weight(vars, node_3);

            VectorType ruling_01 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_01));
            VectorType ruling_12 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_12));
            VectorType ruling_23 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_23));
            VectorType ruling_30 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_30));

            Dest.segment(_cons_idx["ruling_1"] + 3*i_nobdry,3) = _vars_idx.reweighted_ruling_1(vars, i_nobdry) - (w1 + w2)*ruling_12 + (w3 + w0)*(ruling_30);
        }
    }

    void get_constraint_ruling_2(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _it_he.reset();
        _it_face.reset();

        for(_it_face; _it_face.valid(); _it_face++){

            int i = _it_face.idx();
            int i_nobdry = _it_face.idx_nobdry();

            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
            int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
            int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
            int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

            RealType w0 = _vars_idx.vertex_weight(vars, node_0);
            RealType w1 = _vars_idx.vertex_weight(vars, node_1);
            RealType w2 = _vars_idx.vertex_weight(vars, node_2);
            RealType w3 = _vars_idx.vertex_weight(vars, node_3);

            VectorType ruling_01 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_01));
            VectorType ruling_12 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_12));
            VectorType ruling_23 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_23));
            VectorType ruling_30 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_30));

            Dest.segment(_cons_idx["ruling_2"] + 3*i_nobdry,3) = _vars_idx.reweighted_ruling_2(vars, i_nobdry) - (w0 + w1)*ruling_01 + (w2 + w3)*(ruling_23);
        }
    }

    void get_constraint_dev(const VectorType &vars, VectorType &Dest) const{

        size_t num_constraints = _cons_idx["num_cons"];

        std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _it_face.reset();

        for(_it_face; _it_face.valid(); _it_face++){
            int i = _it_face.idx();
            int i_nobdry = _it_face.idx_nobdry();
            Dest.segment(_cons_idx["dev"] + 3*i_nobdry,3) = (_vars_idx.reweighted_ruling_1(vars, i_nobdry). template head<3>()).cross(_vars_idx.reweighted_ruling_2(vars, i_nobdry). template head<3>());
        }

        VectorType constraint_dev = Dest.segment(_cons_idx["dev"], 3*(_num_faces - bdryFaces.size()));
        //std::cout<<"Developability constraint: "<<constraint_dev.norm()<<std::endl;
    }

    void get_constraint_fair_v(const VectorType &vars, VectorType &Dest) const
    {

        size_t num_constraints = _cons_idx["num_cons"];
        
        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _stripHandle.num_vertexStrips(); i++){
            auto triple = _stripHandle.vertex_strip(i);
            int node_0 = std::get<0>(triple);
            int node_1 = std::get<1>(triple);
            int node_2 = std::get<2>(triple);
            Dest.segment(_cons_idx["fair_v"] + 3*i,3) = _vars_idx.vertex(vars, node_0) - 2*_vars_idx.vertex(vars, node_1) + _vars_idx.vertex(vars, node_2);
        }
    }

    void get_constraint_fair_n(const VectorType &vars, VectorType &Dest) const
    {
        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints)
        {
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _stripHandle.num_faceStrips(); i++)
        {
            auto triple = _stripHandle.face_strip(i);
            int face_0 = std::get<0>(triple);
            int face_1 = std::get<1>(triple);
            int face_2 = std::get<2>(triple);

            VectorType normal_0 = _vars_idx.face_normal(vars, face_0);
            VectorType normal_1 = _vars_idx.face_normal(vars, face_1);
            VectorType normal_2 = _vars_idx.face_normal(vars, face_2);

            Dest.segment(_cons_idx["fair_n"] + 3*i,3) = normal_0 - 2*normal_1 + normal_2;
        }

    }

    void get_constraint_fair_r(const VectorType &vars, VectorType &Dest) const
    {
        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints)
        {
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _stripHandle.num_halfedgeStrips(); i++)
        {
            auto triple = _stripHandle.halfedge_strip(i);
            int he_0 = std::get<0>(triple);
            int he_1 = std::get<1>(triple);
            int he_2 = std::get<2>(triple);

            VectorType ruling_0 = _vars_idx.ruling(vars, _it_he.to_nbdry(he_0));
            VectorType ruling_1 = _vars_idx.ruling(vars, _it_he.to_nbdry(he_1));
            VectorType ruling_2 = _vars_idx.ruling(vars, _it_he.to_nbdry(he_2)); 
            Dest.segment(_cons_idx["fair_r"] + 3*i,3) = ruling_0 - 2*ruling_1 + ruling_2;
        }
    }

    void apply(const VectorType &vars, VectorType &Dest) const override
{
    if(vars.size() != _vars_idx["num_dofs"]){
        std::cout<<"num dofs: "<<_vars_idx["num_dofs"]<<std::endl;
        std::cout<<"vars size: "<<vars.size()<<std::endl;
        throw std::runtime_error("The size of the input vector does not match the number of degrees of freedom.");
    }
    Dest.setZero();
    
    get_constraint_vertex_1(vars, Dest);
    get_constraint_vertex_2(vars, Dest);
    get_constraint_edge_vec_1(vars, Dest);
    get_constraint_edge_vec_2(vars, Dest);
    get_constraint_normal_1(vars, Dest);
    get_constraint_normal_2(vars, Dest);
    get_constraint_normal_3(vars, Dest);
    get_constraint_ruling_0(vars, Dest);
    get_constraint_ruling_1(vars, Dest);
    get_constraint_ruling_2(vars, Dest);
    get_constraint_dev(vars, Dest);
    get_constraint_fair_v(vars, Dest);
    get_constraint_fair_n(vars, Dest);
    get_constraint_fair_r(vars, Dest);
}
};

template <typename ConfiguratorType>
class ConstraintGrad : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    private:
        typedef typename std::vector<Eigen::Triplet<RealType>> TripletVectorType;
        QuadMeshTopologySaver &_quadTopol;
        size_t _num_vertices;
        size_t _num_faces;
        size_t _num_edges;
        size_t _num_halfedges;
        size_t _num_bdry_faces;
        size_t _num_bdry_halfedges;
        mutable SkippingBdryFaceIterator _it_face;
        mutable SkippingBdryHalfEdgeIterator _it_he;
        mutable StripHandler<ConfiguratorType> _stripHandle;

        // Serves to access the current vector of input variables at the right positions
        VarsIdx<ConfiguratorType> _vars_idx;
        ConstraintIdx<ConfiguratorType> _cons_idx;

    public:
        ConstraintGrad(QuadMeshTopologySaver &quadTopol, 
                    StripHandler<ConfiguratorType> &stripHandle,
                    ConstraintIdx<ConfiguratorType> &constraintIdx,
                    VarsIdx<ConfiguratorType> &varsIdx):
          _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges()),
          _num_bdry_faces(quadTopol.getNumBdryFaces()),
          _num_bdry_halfedges(quadTopol.getNumBdryHalfEdges()),
          _stripHandle(stripHandle),
          _it_face(quadTopol),
          _it_he(quadTopol),
          _vars_idx(varsIdx),
          _cons_idx(constraintIdx)
        {
        }

        void get_constraintgrad_vertex_1(const VectorType &vars, MatrixType &Dest) const {

            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows()!= num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }
            
            for(int i = 0; i < _num_vertices; i++){
                RealType dummy_i = _vars_idx.dummy_weight(vars, i);
                // Assign values directly
                Dest.coeffRef(_cons_idx["vertex_1"] + i, i + _vars_idx["dummy_weights"]) = -2.0 * dummy_i;
                Dest.coeffRef(_cons_idx["vertex_1"] + i, i + _vars_idx["weights"]) = 1.0;
            }
        }

        void get_constraintgrad_vertex_1_triplet(const VectorType &vars, TripletVectorType &triplets) const {
            
            for(int i = 0; i < _num_vertices; i++){
                RealType dummy_i = _vars_idx.dummy_weight(vars, i);
                // Assign values directly
                triplets.emplace_back(_cons_idx["vertex_1"] + i, i + _vars_idx["dummy_weights"], -2.0 * dummy_i);
                triplets.emplace_back(_cons_idx["vertex_1"] + i, i + _vars_idx["weights"], 1.0);
            }
        }

        void get_constraintgrad_vertex_2(const VectorType &vars, MatrixType &Dest) const {

            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            auto Id = MatrixType(3,3);
            Id.setIdentity();

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _num_vertices; i++){
                RealType w_i = _vars_idx.vertex_weight(vars, i);
                VectorType vec = -_vars_idx.vertex(vars, i);
                MatrixType vec_mat = vectorToSparseMat(vec);
                assignSparseBlockInplace(Dest,-w_i*Id,_cons_idx["vertex_2"] + 3*i,3*i + _vars_idx["vertices"],triplet);
                assignSparseBlockInplace(Dest,vec_mat,_cons_idx["vertex_2"] + 3*i,i + _vars_idx["weights"],triplet);
                assignSparseBlockInplace(Dest,Id,_cons_idx["vertex_2"] +3*i,3*i + _vars_idx["reweighted_vertices"],triplet);
            }

        }

        void get_constraintgrad_vertex_2_triplet(const VectorType &vars, TripletVectorType &triplet) const {

            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_vertices; i++){
                RealType w_i = _vars_idx.vertex_weight(vars, i);
                VectorType vec = -_vars_idx.vertex(vars, i);
                assignSparseMatBlockTriplet(-w_i*Id,_cons_idx["vertex_2"] + 3*i,3*i + _vars_idx["vertices"],triplet);
                assignSparseVecBlockTriplet(vec,_cons_idx["vertex_2"] + 3*i,i + _vars_idx["weights"],triplet);
                assignSparseMatBlockTriplet(Id,_cons_idx["vertex_2"] +3*i,3*i + _vars_idx["reweighted_vertices"],triplet);
            }

        }

        void get_constraintgrad_edge_vec_1(const VectorType &vars, MatrixType &Dest)const {

            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            // From average edge vector 1
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                RealType w0 = _vars_idx.vertex_weight(vars, node_0);
                RealType w1 = _vars_idx.vertex_weight(vars, node_1);
                RealType w2 = _vars_idx.vertex_weight(vars, node_2);
                RealType w3 = _vars_idx.vertex_weight(vars, node_3);

                VectorType v0_tilde = _vars_idx.reweighted_vertex(node_0);
                VectorType v1_tilde = _vars_idx.reweighted_vertex(node_1);
                VectorType v2_tilde = _vars_idx.reweighted_vertex(node_2);
                VectorType v3_tilde = _vars_idx.reweighted_vertex(node_3);

                assignSparseBlockInplace(Dest,Id,_cons_idx["edge_vec_1"] + 3*i,3*i+_vars_idx["reweighted_edges_1"],triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                MatrixType sum_1 = vectorToSparseMat((v2_tilde + v3_tilde));
                MatrixType sum_2 = vectorToSparseMat((v1_tilde + v0_tilde));
                assignSparseBlockInplace(Dest, -sum_1,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -sum_1,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, sum_2,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, sum_2,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseBlockInplace(Dest, (w2+w3)*Id,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseBlockInplace(Dest, (w2+w3)*Id,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseBlockInplace(Dest, -(w0+w1)*Id,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseBlockInplace(Dest, -(w0+w1)*Id,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_3,triplet);
            }

        }

        void get_constraintgrad_edge_vec_1_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            // From average edge vector 1
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                RealType w0 = _vars_idx.vertex_weight(vars, node_0);
                RealType w1 = _vars_idx.vertex_weight(vars, node_1);
                RealType w2 = _vars_idx.vertex_weight(vars, node_2);
                RealType w3 = _vars_idx.vertex_weight(vars, node_3);

                VectorType v0_tilde = _vars_idx.reweighted_vertex(vars, node_0);
                VectorType v1_tilde = _vars_idx.reweighted_vertex(vars, node_1);
                VectorType v2_tilde = _vars_idx.reweighted_vertex(vars, node_2);
                VectorType v3_tilde = _vars_idx.reweighted_vertex(vars, node_3);

                assignSparseMatBlockTriplet(Id,_cons_idx["edge_vec_1"] + 3*i,3*i+_vars_idx["reweighted_edges_1"],triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                VectorType sum_1 = (v2_tilde + v3_tilde);
                VectorType sum_2 = (v1_tilde + v0_tilde);
                assignSparseVecBlockTriplet(-sum_1,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["weights"] + node_0,triplet);
                assignSparseVecBlockTriplet(-sum_1,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["weights"] + node_1,triplet);
                assignSparseVecBlockTriplet(sum_2,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["weights"] + node_2,triplet);
                assignSparseVecBlockTriplet(sum_2,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseMatBlockTriplet((w2+w3)*Id,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseMatBlockTriplet((w2+w3)*Id,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseMatBlockTriplet(-(w0+w1)*Id,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseMatBlockTriplet(-(w0+w1)*Id,_cons_idx["edge_vec_1"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_3,triplet);
            }

        }

        void get_constraintgrad_edge_vec_2(const VectorType &vars, MatrixType &Dest)const {

            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;
            // From average edge vector 1
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                RealType w0 = _vars_idx.vertex_weight(vars, node_0);
                RealType w1 = _vars_idx.vertex_weight(vars, node_1);
                RealType w2 = _vars_idx.vertex_weight(vars, node_2);
                RealType w3 = _vars_idx.vertex_weight(vars, node_3);

                VectorType v0_tilde = _vars_idx.reweighted_vertex(vars, node_0);
                VectorType v1_tilde = _vars_idx.reweighted_vertex(vars, node_1);
                VectorType v2_tilde = _vars_idx.reweighted_vertex(vars, node_2);
                VectorType v3_tilde = _vars_idx.reweighted_vertex(vars, node_3);

                assignSparseBlockInplace(Dest,Id,_cons_idx["edge_vec_2"] + 3*i,3*i+_vars_idx["reweighted_edges_2"],triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                MatrixType sum_1 = vectorToSparseMat((v2_tilde + v1_tilde));
                MatrixType sum_2 = vectorToSparseMat((v0_tilde + v3_tilde));
    
                assignSparseBlockInplace(Dest, sum_1,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -sum_2,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, -sum_2,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, sum_1,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseBlockInplace(Dest,-(w1+w2)*Id,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseBlockInplace(Dest,(w0+w3)*Id,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseBlockInplace(Dest,(w0+w3)*Id,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseBlockInplace(Dest,-(w1+w2)*Id,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_3,triplet);
            }

        }

        void get_constraintgrad_edge_vec_2_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            // From average edge vector 1
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                RealType w0 = _vars_idx.vertex_weight(vars, node_0);
                RealType w1 = _vars_idx.vertex_weight(vars, node_1);
                RealType w2 = _vars_idx.vertex_weight(vars, node_2);
                RealType w3 = _vars_idx.vertex_weight(vars, node_3);

                VectorType v0_tilde = _vars_idx.reweighted_vertex(vars, node_0);
                VectorType v1_tilde = _vars_idx.reweighted_vertex(vars, node_1);
                VectorType v2_tilde = _vars_idx.reweighted_vertex(vars, node_2);
                VectorType v3_tilde = _vars_idx.reweighted_vertex(vars, node_3);

                assignSparseMatBlockTriplet(Id,_cons_idx["edge_vec_2"] + 3*i,3*i+_vars_idx["reweighted_edges_2"],triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                VectorType sum_1 = v2_tilde + v1_tilde;
                VectorType sum_2 = v0_tilde + v3_tilde;
    
                assignSparseVecBlockTriplet(sum_1,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["weights"] + node_0,triplet);
                assignSparseVecBlockTriplet(-sum_2,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["weights"] + node_1,triplet);
                assignSparseVecBlockTriplet(-sum_2,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["weights"] + node_2,triplet);
                assignSparseVecBlockTriplet(sum_1,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseMatBlockTriplet(-(w1+w2)*Id,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseMatBlockTriplet((w0+w3)*Id,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseMatBlockTriplet((w0+w3)*Id,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseMatBlockTriplet(-(w1+w2)*Id,_cons_idx["edge_vec_2"] + 3*i,_vars_idx["reweighted_vertices"] + 3*node_3,triplet);
            }

        }

        void get_constraintgrad_normal_1(const VectorType &vars, MatrixType &Dest)const {

            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs ){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                MatrixType vec_1 = vectorToSparseMat(_vars_idx.reweighted_edge_1(vars, i),true);
                MatrixType vec_2 = vectorToSparseMat(_vars_idx.face_normal(vars, i),true);
                assignSparseBlockInplace(Dest,vec_1,_cons_idx["normal_1"] + i,_vars_idx["normals"] + 3*i,triplet);
                assignSparseBlockInplace(Dest,vec_2,_cons_idx["normal_1"] + i,_vars_idx["reweighted_edges_1"] + 3*i,triplet);
            }
        }

        void get_constraintgrad_normal_1_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                assignSparseVecBlockTriplet(_vars_idx.reweighted_edge_1(vars, i),_cons_idx["normal_1"] + i,_vars_idx["normals"] + 3*i,triplet,true);
                assignSparseVecBlockTriplet(_vars_idx.face_normal(vars, i),_cons_idx["normal_1"] + i,_vars_idx["reweighted_edges_1"] + 3*i,triplet,true);
            }
        }

        void get_constraintgrad_normal_2(const VectorType &vars, MatrixType &Dest)const {

            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs ){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                MatrixType vec_1 = vectorToSparseMat(_vars_idx.reweighted_edge_2(vars, i),true);
                MatrixType vec_2 = vectorToSparseMat(_vars_idx.face_normal(vars, i),true);
                assignSparseBlockInplace(Dest,vec_1,_cons_idx["normal_2"] + i,_vars_idx["normals"] + 3*i,triplet);
                assignSparseBlockInplace(Dest,vec_2,_cons_idx["normal_2"] + i,_vars_idx["reweighted_edges_2"] + 3*i,triplet);
            }
        }

        void get_constraintgrad_normal_2_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                assignSparseVecBlockTriplet(_vars_idx.reweighted_edge_2(vars, i),_cons_idx["normal_2"] + i,_vars_idx["normals"] + 3*i,triplet,true);
                assignSparseVecBlockTriplet(_vars_idx.face_normal(vars, i),_cons_idx["normal_2"] + i,_vars_idx["reweighted_edges_2"] + 3*i,triplet,true);
            }
        }

        void get_constraintgrad_normal_3(const VectorType &vars, MatrixType &Dest)const {

            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs ){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                MatrixType vec_1 = vectorToSparseMat(_vars_idx.face_normal(vars, i),true);
                assignSparseBlockInplace(Dest,2.0*vec_1,_cons_idx["normal_3"] + i,_vars_idx["normals"] + 3*i,triplet);
            }
        }

        void get_constraintgrad_normal_3_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                assignSparseVecBlockTriplet(2.0*_vars_idx.face_normal(vars, i),_cons_idx["normal_3"] + i,_vars_idx["normals"] + 3*i,triplet,true);
            }
        }
        
        void get_constraintgrad_ruling_0(const VectorType &vars, MatrixType &Dest)const {

            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;

            _it_he.reset();

            for(_it_he; _it_he.valid();_it_he++)
            {
                int i = _it_he.idx();
                int i_nobdry = _it_he.idx_nobdry();

                size_t faceIdx_1 = _quadTopol.getFaceOfHalfEdge(i,0);
                size_t faceIdx_2 = _quadTopol.getFaceOfHalfEdge(i,1);
                Eigen::Vector3d normal_1 = _vars_idx.face_normal(vars, faceIdx_1).template head<3>();
                Eigen::Vector3d normal_2 = _vars_idx.face_normal(vars, faceIdx_2).template head<3>();

                // First, set the derivative w.r.t. the ruling vector
                auto Id = MatrixType(3,3);
                Id.setIdentity();

                assignSparseBlockInplace(Dest,Id,_cons_idx["ruling_0"] + 3*(i_nobdry),_vars_idx["rulings"] + 3*i_nobdry,triplet);

                // Now, set the derivatives w.r.t. the cross product
                // First, differentiate w.r.t. normal_1
                Eigen::Vector3d e;
                MatrixType vec_prod;
                e.setZero();
                e[0] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_1,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_1 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_1 + 2,triplet);

                // Now, differentiate w.r.t. normal_2
                e[2] = 0.0;
                e[0] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_2,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_2 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_2 + 2,triplet);
            }
        }

        void get_constraintgrad_ruling_0_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            _it_he.reset();

            for(_it_he; _it_he.valid();_it_he++)
            {
                int i = _it_he.idx();
                int i_nobdry = _it_he.idx_nobdry();

                size_t faceIdx_1 = _quadTopol.getFaceOfHalfEdge(i,0);
                size_t faceIdx_2 = _quadTopol.getFaceOfHalfEdge(i,1);
                Eigen::Vector3d normal_1 = _vars_idx.face_normal(vars, faceIdx_1).template head<3>();
                Eigen::Vector3d normal_2 = _vars_idx.face_normal(vars, faceIdx_2).template head<3>();

                // First, set the derivative w.r.t. the ruling vector
                auto Id = MatrixType(3,3);
                Id.setIdentity();

                assignSparseMatBlockTriplet(Id,_cons_idx["ruling_0"] + 3*(i_nobdry),_vars_idx["rulings"] + 3*i_nobdry,triplet);

                // Now, set the derivatives w.r.t. the cross product
                // First, differentiate w.r.t. normal_1
                Eigen::Vector3d e;
                VectorType vec_prod;
                e.setZero();
                e[0] = 1.0;
                vec_prod = -e.cross(normal_2);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_1,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = -e.cross(normal_2);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_1 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = -e.cross(normal_2);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_1 + 2, triplet);

                // Now, differentiate w.r.t. normal_2
                e[2] = 0.0;
                e[0] = 1.0;
                vec_prod = e.cross(normal_1);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_2,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = e.cross(normal_1);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_2 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = e.cross(normal_1);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_vars_idx["normals"] + 3*faceIdx_2 + 2,triplet);
            }
        }

        void get_constraintgrad_ruling_1(const VectorType& vars, MatrixType &Dest)const {

            size_t num_dofs = _vars_idx["num_dofs"];

            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            _it_he.reset();
            _it_face.reset();

            for(_it_face; _it_face.valid(); _it_face++){

                int i = _it_face.idx();
                int i_nobdry = _it_face.idx_nobdry();

                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                auto node_0_vert = _vars_idx.vertex(vars, node_0);
                auto node_1_vert = _vars_idx.vertex(vars, node_1);
                auto node_2_vert = _vars_idx.vertex(vars, node_2);
                auto node_3_vert = _vars_idx.vertex(vars, node_3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                RealType w0 = _vars_idx.vertex_weight(vars, node_0);
                RealType w1 = _vars_idx.vertex_weight(vars, node_1);
                RealType w2 = _vars_idx.vertex_weight(vars, node_2);
                RealType w3 = _vars_idx.vertex_weight(vars, node_3);

                VectorType ruling_01 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_01));
                VectorType ruling_12 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_12));
                VectorType ruling_23 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_23));
                VectorType ruling_30 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_30));

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}

                assignSparseBlockInplace(Dest,Id,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["reweighted_rulings_1"] + 3*i_nobdry,triplet);

                MatrixType ruling_30_mat = vectorToSparseMat(ruling_30);
                MatrixType ruling_12_mat = vectorToSparseMat(ruling_12);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseBlockInplace(Dest, ruling_30_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -ruling_12_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, -ruling_12_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, ruling_30_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseBlockInplace(Dest,-(w1 + w2)*Id,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["rulings"] + 3*_it_he.to_nbdry(heh_12),triplet);
                assignSparseBlockInplace(Dest,(w3 + w0)*Id,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["rulings"] + 3*_it_he.to_nbdry(heh_30),triplet);
            }
        }

        void get_constraintgrad_ruling_1_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            auto Id = MatrixType(3,3);
            Id.setIdentity();

            _it_he.reset();
            _it_face.reset();

            for(_it_face; _it_face.valid(); _it_face++){

                int i = _it_face.idx();
                int i_nobdry = _it_face.idx_nobdry();

                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                auto node_0_vert = _vars_idx.vertex(vars, node_0);
                auto node_1_vert = _vars_idx.vertex(vars, node_1);
                auto node_2_vert = _vars_idx.vertex(vars, node_2);
                auto node_3_vert = _vars_idx.vertex(vars, node_3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                RealType w0 = _vars_idx.vertex_weight(vars, node_0);
                RealType w1 = _vars_idx.vertex_weight(vars, node_1);
                RealType w2 = _vars_idx.vertex_weight(vars, node_2);
                RealType w3 = _vars_idx.vertex_weight(vars, node_3);

                VectorType ruling_01 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_01));
                VectorType ruling_12 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_12));
                VectorType ruling_23 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_23));
                VectorType ruling_30 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_30));

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}

                assignSparseMatBlockTriplet(Id,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["reweighted_rulings_1"] + 3*i_nobdry,triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseVecBlockTriplet(ruling_30,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["weights"] + node_0,triplet);
                assignSparseVecBlockTriplet(-ruling_12,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["weights"] + node_1,triplet);
                assignSparseVecBlockTriplet(-ruling_12,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["weights"] + node_2,triplet);
                assignSparseVecBlockTriplet(ruling_30,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseMatBlockTriplet(-(w1 + w2)*Id,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["rulings"] + 3*_it_he.to_nbdry(heh_12),triplet);
                assignSparseMatBlockTriplet((w3 + w0)*Id,_cons_idx["ruling_1"] + 3*i_nobdry,_vars_idx["rulings"] + 3*_it_he.to_nbdry(heh_30),triplet);
            }
        }

        void get_constraintgrad_ruling_2(const VectorType& vars, MatrixType &Dest)const {

            std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

            size_t num_dofs = _vars_idx["num_dofs"];

            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            _it_face.reset();

            for(_it_face; _it_face.valid(); _it_face++){

                int i = _it_face.idx();
                int i_nbdry = _it_face.idx_nobdry();

                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                RealType w0 = _vars_idx.vertex_weight(vars, node_0);
                RealType w1 = _vars_idx.vertex_weight(vars, node_1);
                RealType w2 = _vars_idx.vertex_weight(vars, node_2);
                RealType w3 = _vars_idx.vertex_weight(vars, node_3);

                VectorType ruling_01 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_01));
                VectorType ruling_23 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_23));

                MatrixType ruling_01_mat = vectorToSparseMat(ruling_01);
                MatrixType ruling_23_mat = vectorToSparseMat(ruling_23);

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}
                assignSparseBlockInplace(Dest,Id,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["reweighted_rulings_2"] + 3*i_nbdry,triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseBlockInplace(Dest, -ruling_01_mat,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -ruling_01_mat,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, ruling_23_mat,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, ruling_23_mat,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseBlockInplace(Dest,-(w0 + w1)*Id,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["rulings"] + 3*_it_he.to_nbdry(heh_01),triplet);
                assignSparseBlockInplace(Dest,(w2 + w3)*Id,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["rulings"] + 3*_it_he.to_nbdry(heh_23),triplet);
            }
        }

        void get_constraintgrad_ruling_2_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            auto Id = MatrixType(3,3);
            Id.setIdentity();

            _it_face.reset();

            for(_it_face; _it_face.valid(); _it_face++){

                int i = _it_face.idx();
                int i_nbdry = _it_face.idx_nobdry();

                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                RealType w0 = _vars_idx.vertex_weight(vars, node_0);
                RealType w1 = _vars_idx.vertex_weight(vars, node_1);
                RealType w2 = _vars_idx.vertex_weight(vars, node_2);
                RealType w3 = _vars_idx.vertex_weight(vars, node_3);

                VectorType ruling_01 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_01));
                VectorType ruling_23 = _vars_idx.ruling(vars, _it_he.to_nbdry(heh_23));

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}
                assignSparseMatBlockTriplet(Id,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["reweighted_rulings_2"] + 3*i_nbdry,triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseVecBlockTriplet(-ruling_01,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["weights"] + node_0,triplet);
                assignSparseVecBlockTriplet(-ruling_01,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["weights"] + node_1,triplet);
                assignSparseVecBlockTriplet(ruling_23,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["weights"] + node_2,triplet);
                assignSparseVecBlockTriplet(ruling_23,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseMatBlockTriplet(-(w0 + w1)*Id,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["rulings"] + 3*_it_he.to_nbdry(heh_01),triplet);
                assignSparseMatBlockTriplet((w2 + w3)*Id,_cons_idx["ruling_2"] + 3*i_nbdry,_vars_idx["rulings"] + 3*_it_he.to_nbdry(heh_23),triplet);
            }
        }
        
        void get_constraintgrad_dev(const VectorType &vars, MatrixType &Dest)const {

            std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

            size_t num_constraints = _cons_idx["num_cons"];
            size_t num_dofs = _vars_idx["num_dofs"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;

            _it_face.reset();

            for(_it_face; _it_face.valid(); _it_face++)
            {

                int i = _it_face.idx();
                int i_nobdry = _it_face.idx_nobdry();

                // obtain the reweighted ruling vectors
                Eigen::Vector3d reweighted_ruling_1 = _vars_idx.reweighted_ruling_1(vars, i_nobdry).template head<3>();
                Eigen::Vector3d reweighted_ruling_2 = _vars_idx.reweighted_ruling_2(vars, i_nobdry).template head<3>();

                // Derivatives w.r.t. the first ruling vector \Tilde{r}_{13,f}
                Eigen::Vector3d e;
                e.setZero();
                e[0] = 1.0; 
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_1"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_1"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_1"] + 2,triplet);
            
                // Derivatives w.r.t. the second ruling vector \Tilde{r}_{02,f}
                e[2] = 0.0;
                e[0] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_2"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_2"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_2"] + 2,triplet);
            }
        }

        void get_constraintgrad_dev_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            _it_face.reset();

            for(_it_face; _it_face.valid(); _it_face++)
            {

                int i = _it_face.idx();
                int i_nobdry = _it_face.idx_nobdry();

                // obtain the reweighted ruling vectors
                Eigen::Vector3d reweighted_ruling_1 = _vars_idx.reweighted_ruling_1(vars, i_nobdry).template head<3>();
                Eigen::Vector3d reweighted_ruling_2 = _vars_idx.reweighted_ruling_2(vars, i_nobdry).template head<3>();

                // Derivatives w.r.t. the first ruling vector \Tilde{r}_{13,f}
                Eigen::Vector3d e;
                e.setZero();
                e[0] = 1.0; 
                assignSparseVecBlockTriplet(e.cross(reweighted_ruling_2),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_1"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                assignSparseVecBlockTriplet(e.cross(reweighted_ruling_2),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_1"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseVecBlockTriplet(e.cross(reweighted_ruling_2),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_1"] + 2,triplet);
            
                // Derivatives w.r.t. the second ruling vector \Tilde{r}_{02,f}
                e[2] = 0.0;
                e[0] = 1.0;
                assignSparseVecBlockTriplet(reweighted_ruling_1.cross(e),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_2"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                assignSparseVecBlockTriplet(reweighted_ruling_1.cross(e),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_2"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseVecBlockTriplet(reweighted_ruling_1.cross(e),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _vars_idx["reweighted_rulings_2"] + 2,triplet);
            }
        }

        void get_constraintgrad_fair_v(const VectorType &vars, MatrixType &Dest)const {

            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs)
            {
                Dest.resize(num_constraints, num_dofs);
            }

            MatrixType Id(3,3);
            Id.setIdentity();
            std::vector<Eigen::Triplet<double>> triplet;

            for(int i = 0; i < _stripHandle.num_vertexStrips(); i++){
                auto triple = _stripHandle.vertex_strip(i);
                int node_0 = std::get<0>(triple);
                int node_1 = std::get<1>(triple);
                int node_2 = std::get<2>(triple);
                
                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_v"] + 3*i,3*node_0 + _vars_idx["vertices"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,_cons_idx["fair_v"] + 3*i,3*node_1 + _vars_idx["vertices"], triplet);
                assignSparseBlockInplace(Dest, Id, _cons_idx["fair_v"] + 3*i,3*node_2 + _vars_idx["vertices"], triplet);
            }
        }

        void get_constraintgrad_fair_v_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            MatrixType Id(3,3);
            Id.setIdentity();

            for(int i = 0; i < _stripHandle.num_vertexStrips(); i++){
                auto triple = _stripHandle.vertex_strip(i);
                int node_0 = std::get<0>(triple);
                int node_1 = std::get<1>(triple);
                int node_2 = std::get<2>(triple);
                
                assignSparseMatBlockTriplet(Id,_cons_idx["fair_v"] + 3*i,3*node_0 + _vars_idx["vertices"],triplet);
                assignSparseMatBlockTriplet(-2*Id,_cons_idx["fair_v"] + 3*i,3*node_1 + _vars_idx["vertices"], triplet);
                assignSparseMatBlockTriplet(Id, _cons_idx["fair_v"] + 3*i,3*node_2 + _vars_idx["vertices"], triplet);
            }
        }

        void get_constraintgrad_fair_n(const VectorType &vars, MatrixType &Dest)const {

            size_t num_constraints = _cons_idx["num_cons"];
            size_t num_dofs = _vars_idx["num_dofs"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            MatrixType Id(3,3);
            Id.setIdentity();

            std::vector<Eigen::Triplet<double>> triplet;

            for(int i = 0; i < _stripHandle.num_faceStrips(); i++){
                auto triple = _stripHandle.face_strip(i);
                int normal_0 = std::get<0>(triple);
                int normal_1 = std::get<1>(triple);
                int normal_2 = std::get<2>(triple);

                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_n"] + 3*i,3*normal_0 + _vars_idx["normals"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,_cons_idx["fair_n"] + 3*i,3*normal_1 + _vars_idx["normals"],triplet);
                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_n"] + 3*i,3*normal_2 + _vars_idx["normals"],triplet);
            }
        }

        void get_constraintgrad_fair_n_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            MatrixType Id(3,3);
            Id.setIdentity();

            for(int i = 0; i < _stripHandle.num_faceStrips(); i++){
                auto triple = _stripHandle.face_strip(i);
                int normal_0 = std::get<0>(triple);
                int normal_1 = std::get<1>(triple);
                int normal_2 = std::get<2>(triple);

                assignSparseMatBlockTriplet(Id,_cons_idx["fair_n"] + 3*i,3*normal_0 + _vars_idx["normals"],triplet);
                assignSparseMatBlockTriplet(-2*Id,_cons_idx["fair_n"] + 3*i,3*normal_1 + _vars_idx["normals"],triplet);
                assignSparseMatBlockTriplet(Id,_cons_idx["fair_n"] + 3*i,3*normal_2 + _vars_idx["normals"],triplet);
            }
        }

        void get_constraintgrad_fair_r(const VectorType &vars, MatrixType &Dest)const {

            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs)
            {
                Dest.resize(num_constraints, num_dofs);
            }

            MatrixType Id(3,3);
            Id.setIdentity();

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _stripHandle.num_halfedgeStrips(); i++)
            {
                auto triple = _stripHandle.halfedge_strip(i);
                int heh_0 = std::get<0>(triple);
                int heh_1 = std::get<1>(triple);
                int heh_2 = std::get<2>(triple);

                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_0) + _vars_idx["rulings"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_1) + _vars_idx["rulings"],triplet);
                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_2) + _vars_idx["rulings"],triplet);
            }
        }

        void get_constraintgrad_fair_r_triplet(const VectorType &vars, TripletVectorType &triplet)const {

            MatrixType Id(3,3);
            Id.setIdentity();

            for(int i = 0; i < _stripHandle.num_halfedgeStrips(); i++)
            {
                auto triple = _stripHandle.halfedge_strip(i);
                int heh_0 = std::get<0>(triple);
                int heh_1 = std::get<1>(triple);
                int heh_2 = std::get<2>(triple);

                assignSparseMatBlockTriplet(Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_0) + _vars_idx["rulings"],triplet);
                assignSparseMatBlockTriplet(-2*Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_1) + _vars_idx["rulings"],triplet);
                assignSparseMatBlockTriplet(Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_2) + _vars_idx["rulings"],triplet);
            }
        }


        void apply(const VectorType &vars, MatrixType &Dest) const override
        {
            size_t num_dofs = _vars_idx["num_dofs"];
            // Change num_constraints back later !
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            Dest.setZero();
            
            TripletVectorType triplet;
            size_t num_vertex_1_triplets = 2*_num_vertices;
            size_t num_vertex_2_triplets = 21*_num_vertices;
            size_t num_edge_vec_1_triplets = 57*_num_faces;
            size_t num_edge_vec_2_triplets = 57*_num_faces;
            size_t num_normal_1_triplets = 6*_num_faces;
            size_t num_normal_2_triplets = 6*_num_faces;
            size_t num_normal_3_triplets = 3*_num_faces;
            size_t num_ruling_0_triplets = 27*(_num_faces - _num_bdry_faces);
            size_t num_ruling_1_triplets = 39*(_num_faces - _num_bdry_faces);
            size_t num_ruling_2_triplets = 39*(_num_faces - _num_bdry_faces);
            size_t num_dev_triplets = 18*(_num_faces - _num_bdry_faces);
            size_t num_fair_v_triplets = 27*_stripHandle.num_vertexStrips();
            size_t num_fair_n_triplets = 27*_stripHandle.num_faceStrips();
            size_t num_fair_r_triplets = 27*_stripHandle.num_halfedgeStrips();

            triplet.reserve(num_vertex_1_triplets +
                            num_vertex_2_triplets +
                            num_edge_vec_1_triplets +
                            num_edge_vec_2_triplets +
                            num_normal_1_triplets +
                            num_normal_2_triplets +
                            num_normal_3_triplets +
                            num_ruling_0_triplets +
                            num_ruling_1_triplets +
                            num_ruling_2_triplets +
                            num_dev_triplets +
                            num_fair_v_triplets +
                            num_fair_n_triplets +
                            num_fair_r_triplets);

            // 2*num_vertices for vertex_2 constraint
            get_constraintgrad_vertex_1_triplet(vars, triplet);
            get_constraintgrad_vertex_2_triplet(vars, triplet);
            get_constraintgrad_edge_vec_1_triplet(vars, triplet);
            get_constraintgrad_edge_vec_2_triplet(vars, triplet);
            get_constraintgrad_normal_1_triplet(vars, triplet);
            get_constraintgrad_normal_2_triplet(vars, triplet);
            get_constraintgrad_normal_3_triplet(vars, triplet);
            get_constraintgrad_ruling_0_triplet(vars, triplet);
            get_constraintgrad_ruling_1_triplet(vars, triplet);
            get_constraintgrad_ruling_2_triplet(vars,triplet);
            get_constraintgrad_dev_triplet(vars, triplet);
            get_constraintgrad_fair_v_triplet(vars, triplet);
            get_constraintgrad_fair_n_triplet(vars, triplet);
            get_constraintgrad_fair_r_triplet(vars, triplet);

            Dest.setFromTriplets(triplet.begin(), triplet.end());
        }

    };

template <typename ConfiguratorType>
class ConstraintHessian : public BaseOp<typename ConfiguratorType::VectorType, GenericTensor<typename ConfiguratorType::SparseMatrixType>>{
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;
        typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
        typedef GenericTensor<typename ConfiguratorType::SparseMatrixType> GenericTensorType;

    private:
        typedef typename std::vector<Eigen::Triplet<RealType>> TripletVectorType;
        QuadMeshTopologySaver &_quadTopol;
        size_t _num_vertices;
        size_t _num_faces;
        size_t _num_edges;
        size_t _num_halfedges;
        size_t _num_bdry_faces;
        size_t _num_bdry_halfedges;
        mutable SkippingBdryFaceIterator _it_face;
        mutable SkippingBdryHalfEdgeIterator _it_he;
        mutable StripHandler<ConfiguratorType> _stripHandle;

        // Serves to access the current vector of input variables at the right positions
        VarsIdx<ConfiguratorType> _vars_idx;
        ConstraintIdx<ConfiguratorType> _cons_idx;

    public:
        ConstraintHessian(QuadMeshTopologySaver &quadTopol, 
            StripHandler<ConfiguratorType> &stripHandle,
            ConstraintIdx<ConfiguratorType> &constraintIdx,
            VarsIdx<ConfiguratorType> &varsIdx)
            : _quadTopol(quadTopol),
              _num_vertices(quadTopol.getNumVertices()),
              _num_faces(quadTopol.getNumFaces()),
              _num_edges(quadTopol.getNumEdges()),
              _num_halfedges(2*_num_edges),
              _num_bdry_faces(quadTopol.getBdryFaces().size()),
              _num_bdry_halfedges(quadTopol.getBdryHalfEdges().size()),
              _stripHandle(stripHandle),
              _it_face(_quadTopol),
              _it_he(_quadTopol),
              _vars_idx(varsIdx),
              _cons_idx(constraintIdx)
              {}

        void apply(const VectorType &vars, GenericTensorType &Dest) const override
        {
            std::vector<TripletVectorType> triplet_vector;
            Dest.resize(3*_quadTopol.getNumVertices(), _vars_idx["num_dofs"], _vars_idx["num_dofs"]);
            triplet_vector.resize(3*_quadTopol.getNumVertices());
            for(int i = 0; i < Dest.size(); i++)
            {
                Dest[i].resize(_vars_idx["num_dofs"], _vars_idx["num_dofs"]);
                Dest[i].setZero();
                triplet_vector[i].reserve(2);
            }
            for(int i = 0; i < 3*_num_vertices; i++)
            {
                // \partial_{}
                RealType weight = _vars_idx.vertex_weight(vars, i / 3);
                triplet_vector[i].push_back(Eigen::Triplet<RealType>(_vars_idx["vertices"] + i, _vars_idx["weights"] + i/3, -1.0));
                triplet_vector[i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + i/3, _vars_idx["vertices"] + i, -1.0));
            }

            for(int i = 0; i < Dest.size(); i++)
            {
                Dest[i].setFromTriplets(triplet_vector[i].begin(), triplet_vector[i].end());
            }
        }

        void get_constrainthess_vertex_1_triplet(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            for(int i = 0; i < _num_vertices; i++){
                // \partial_{\omega_i} \partial \omega_{i}
                triplet_vector[_cons_idx["vertex_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["dummy_weights"] + i, _vars_idx["dummy_weights"] + i, -2.0));
            }
        }

        void get_constrainthess_vertex_2_triplet(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            for(int i = 0; i < 3*_num_vertices; i++)
            {
                // \partial_{}
                RealType weight = _vars_idx.vertex_weight(vars, i / 3);
                triplet_vector[_cons_idx["vertex_2"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["vertices"] + i, _vars_idx["weights"] + i/3, -1.0));
                triplet_vector[_cons_idx["vertex_2"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + i/3, _vars_idx["vertices"] + i, -1.0));
            }
        }

        void get_constrainthess_edge_vec_1_triplet(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            for(int i = 0; i < 3*_num_faces; i++)
            {
                int faceIdx = i/3;
                int component = i % 3;

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(faceIdx,0);
                int node_1 = _quadTopol.getNodeOfQuad(faceIdx,1);
                int node_2 = _quadTopol.getNodeOfQuad(faceIdx,2);
                int node_3 = _quadTopol.getNodeOfQuad(faceIdx,3);

                // -(w_0 + w_1)*(\Tilde{v}_2 + \Tilde{v}_3) + (w_2 + w_3)*(\Tilde{v}_0 + \Tilde{v}_1)

                // \partial_{w_0} \partial_{\Tilde{v_2}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_0, _vars_idx["reweighted_vertices"] + 3*node_2 + component, -1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_2 + component, _vars_idx["weights"] + node_0, -1.0));
                // \partial_{w_0} \partial_{\Tilde{v_3}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_0, _vars_idx["reweighted_vertices"] + 3*node_3 + component, -1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_3 + component, _vars_idx["weights"] + node_0, -1.0));
                // \partial_{w_1} \partial_{\Tilde{v_2}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_1, _vars_idx["reweighted_vertices"] + 3*node_2 + component, -1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_2 + component, _vars_idx["weights"] + node_1, -1.0));
                // \partial_{w_1} \partial_{\Tilde{v_3}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_1, _vars_idx["reweighted_vertices"] + 3*node_3 + component, -1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_3 + component, _vars_idx["weights"] + node_1, -1.0));

                // Now, the other term
                // \partial_{w_2} \partial_{\Tilde{v_1}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_2, _vars_idx["reweighted_vertices"] + 3*node_1 + component, 1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_1 + component, _vars_idx["weights"] + node_2, 1.0));
                // \partial_{w_2} \partial_{\Tilde{v_0}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_2, _vars_idx["reweighted_vertices"] + 3*node_0 + component, 1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_0 + component, _vars_idx["weights"] + node_2, 1.0));
                // \partial_{w_3} \partial_{\Tilde{v_1}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_3, _vars_idx["reweighted_vertices"] + 3*node_1 + component, 1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_1 + component, _vars_idx["weights"] + node_3, 1.0));
                // \partial_{w_3} \partial_{\Tilde{v_0}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_3, _vars_idx["reweighted_vertices"] + 3*node_0 + component, 1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_0 + component, _vars_idx["weights"] + node_3, 1.0));
            }
        }

        void get_constrainthess_edge_vec_2_triplet(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            for(int i = 0; i < 3*_num_faces; i++)
            {
                int faceIdx = i/3;
                int component = i % 3;

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(faceIdx,0);
                int node_1 = _quadTopol.getNodeOfQuad(faceIdx,1);
                int node_2 = _quadTopol.getNodeOfQuad(faceIdx,2);
                int node_3 = _quadTopol.getNodeOfQuad(faceIdx,3);

                // \partial_{w_1} \partial_{\Tilde{v_0}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_1, _vars_idx["reweighted_vertices"] + 3*node_0 + component, -1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_0 + component, _vars_idx["weights"] + node_1, -1.0));
                //\partial_{w_1} \partial_{\Tilde{v_3}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_1, _vars_idx["reweighted_vertices"] + 3*node_3 + component, -1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_3 + component, _vars_idx["weights"] + node_1, -1.0));
                // \partial_{w_2} \partial_{\Tilde{v_0}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_2, _vars_idx["reweighted_vertices"] + 3*node_0 + component, -1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_0 + component, _vars_idx["weights"] + node_2, -1.0));
                // \partial_{w_2} \partial_{\Tilde{v_3}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_2, _vars_idx["reweighted_vertices"] + 3*node_3 + component, -1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_3 + component, _vars_idx["weights"] + node_2, -1.0));

                // Now, the other term
                // \partial_{w_0} \partial_{\Tilde{v_2}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_0, _vars_idx["reweighted_vertices"] + 3*node_2 + component, 1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_2 + component, _vars_idx["weights"] + node_0, 1.0));
                // \partial_{w_0} \partial_{\Tilde{v_1}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_0, _vars_idx["reweighted_vertices"] + 3*node_1 + component, 1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_1 + component, _vars_idx["weights"] + node_0, 1.0));
                // \partial_{w_3} \partial_{\Tilde{v_2}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_3, _vars_idx["reweighted_vertices"] + 3*node_2 + component, 1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_2 + component, _vars_idx["weights"] + node_3, 1.0));
                // \partial_{w_3} \partial_{\Tilde{v_1}}
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_3, _vars_idx["reweighted_vertices"] + 3*node_1 + component, 1.0));
                triplet_vector[_cons_idx["edge_vec_1"] + i].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_vertices"] + 3*node_1 + component, _vars_idx["weights"] + node_3, 1.0));
            }
        }

        void get_constrainthess_normal_1_triplet(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            for(int i = 0; i < 3*_num_faces; i++)
            {
                int faceIdx = i /3;
                int component = i % 3;

                triplet_vector[_cons_idx["normal_1"]].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + i, _vars_idx["reweighted_edges_1"] + i, 1.0));
                triplet_vector[_cons_idx["normal_1"]].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_edges_1"] + i, _vars_idx["normals"] + i, 1.0));
            }
        }

        void get_constrainthess_normal_2_triplet(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            for(int i = 0; i < 3*_num_faces; i++)
            {
                int faceIdx = i /3;
                int component = i % 3;

                triplet_vector[_cons_idx["normal_2"]].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + i, _vars_idx["reweighted_edges_2"] + i, 1.0));
                triplet_vector[_cons_idx["normal_2"]].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_edges_2"] + i, _vars_idx["normals"] + i, 1.0));
            }
        }

        void get_constrainthess_normal_3_triplet(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            for(int i = 0; i < 3*_num_faces; i++)
            {
                int faceIdx = i /3;
                int component = i % 3;

                triplet_vector[_cons_idx["normal_3"]].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + i, _vars_idx["normals"] + i, 1.0));
                triplet_vector[_cons_idx["normal_3"]].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + i, _vars_idx["normals"] + i, 1.0));
            }
        }

        void get_constrainthess_ruling_0_triplet(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const{
            size_t num_dofs = _vars_idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            _it_he.reset();

            for(_it_he; _it_he.valid();_it_he++)
            {
                int i = _it_he.idx();
                int i_nobdry = _it_he.idx_nobdry();

                int face1 = _quadTopol.getFaceOfHalfEdge(i, 0);
                int face2 = _quadTopol.getFaceOfHalfEdge(i, 1);

                // first component r_x = (n_1)_y * (n_2)_z - (n_1)_z * (n_2)_y -> WARNING ALL HAS NEGATIVE SIGN IN CONSTRAINT
                // differentiate w.r.t. first summand
                triplet_vector[_cons_idx["ruling_0"] + 3*i].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face1 + 1, _vars_idx["normals"] + 3*face2 + 2, -1.0));
                triplet_vector[_cons_idx["ruling_0"] + 3*i].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face2 + 2, _vars_idx["normals"] + 3*face1 + 1, -1.0));
            
                // differentiate w.r.t. second summand
                triplet_vector[_cons_idx["ruling_0"] + 3*i].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face1 + 2, _vars_idx["normals"] + 3*face2 + 1, 1.0));
                triplet_vector[_cons_idx["ruling_0"] + 3*i].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face2 + 1, _vars_idx["normals"] + 3*face1 + 2, 1.0));
            
                // second component r_y = (n_1)_z * (n_2)_x - (n_1)_x * (n_2)_z
                // differentiate w.r.t. first summand
                triplet_vector[_cons_idx["ruling_0"] + 3*i + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face1 + 2, _vars_idx["normals"] + 3*face2 + 0, -1.0));
                triplet_vector[_cons_idx["ruling_0"] + 3*i + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face2 + 0, _vars_idx["normals"] + 3*face1 + 2, -1.0));

                // differentiate w.r.t. second summand
                triplet_vector[_cons_idx["ruling_0"] + 3*i + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face1 + 0, _vars_idx["normals"] + 3*face2 + 2, 1.0));
                triplet_vector[_cons_idx["ruling_0"] + 3*i + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face2 + 2, _vars_idx["normals"] + 3*face1 + 0, 1.0));
            
                // third component r_z = (n_1)_x * (n_2)_y - (n_1)_y * (n_2)_x
                // differentiate w.r.t. first summand
                triplet_vector[_cons_idx["ruling_0"] + 3*i + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face1 + 0, _vars_idx["normals"] + 3*face2 + 1, -1.0));
                triplet_vector[_cons_idx["ruling_0"] + 3*i + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face2 + 1, _vars_idx["normals"] + 3*face1 + 0, -1.0));
                // differentiate w.r.t. second summand
                triplet_vector[_cons_idx["ruling_0"] + 3*i + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face1 + 1, _vars_idx["normals"] + 3*face2 + 0, 1.0));
                triplet_vector[_cons_idx["ruling_0"] + 3*i + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["normals"] + 3*face2 + 0, _vars_idx["normals"] + 3*face1 + 1, 1.0));
            }
        }

        void get_constrainthess_ruling_1_triplet(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            _it_he.reset();
            _it_face.reset();

            for(_it_face; _it_face.valid(); _it_face++){

                int i = _it_face.idx();
                int i_nobdry = _it_face.idx_nobdry();

                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                // -[(w_1 + w_2)*r_{12} + (w_3 + w_0)*(-r_{30})]
                // \partial_{w_1} \partial r_{12}
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_1, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12), -1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_1, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12) + 1, -1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_1, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12) + 2, -1.0));
                
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12), _vars_idx["weights"] + node_1, -1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12) + 1, _vars_idx["weights"] + node_1, -1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12) + 2, _vars_idx["weights"] + node_1, -1.0));
                // \partial_{w_2} \partial r_{12}
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_2, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12), -1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_2, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12) + 1, -1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_2, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12) + 2, -1.0));

                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12), _vars_idx["weights"] + node_2, -1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12) + 1, _vars_idx["weights"] + node_2, -1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_12) + 2, _vars_idx["weights"] + node_2, -1.0));
                // \partial_{w_3} \partial r_{30}
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_3, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30), 1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_3, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30) + 1, 1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_3, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30) + 2, 1.0));
                
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30), _vars_idx["weights"] + node_3, 1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30) + 1, _vars_idx["weights"] + node_3, 1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30) + 2, _vars_idx["weights"] + node_3, 1.0));
                // \partial_{w_0} \partial r_{30}
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_0, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30), 1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_0, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30) + 1, 1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_0, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30) + 2, 1.0));
                
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30), _vars_idx["weights"] + node_0, 1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30) + 1, _vars_idx["weights"] + node_0, 1.0));
                triplet_vector[_cons_idx["ruling_1"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_30) + 2, _vars_idx["weights"] + node_0, 1.0));
            }
        }

        void get_constrainthess_ruling_2_triplet(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            _it_he.reset();
            _it_face.reset();

            for(_it_face; _it_face.valid(); _it_face++){

                int i = _it_face.idx();
                int i_nobdry = _it_face.idx_nobdry();

                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                // -[(w_0 + w_1)*r_{01} + (w_2 + w_3)*(-r_{23})]
                // \partial_{w_0} \partial r_{01}
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_0, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01), -1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_0, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01) + 1, -1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_0, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01) + 2, -1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01), _vars_idx["weights"] + node_0, -1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01) + 1, _vars_idx["weights"] + node_0, -1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01) + 2, _vars_idx["weights"] + node_0, -1.0));
                // \partial_{w_1} \partial r_{01}
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_1, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01), -1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_1, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01) + 1, -1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_1, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01) + 2, -1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01), _vars_idx["weights"] + node_1, -1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01) + 1, _vars_idx["weights"] + node_1, -1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_01) + 2, _vars_idx["weights"] + node_1, -1.0));
                // \partial_{w_2} \partial r_{23}
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_2, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23), 1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_2, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23) + 1, 1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_2, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23) + 2, 1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23), _vars_idx["weights"] + node_2, 1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23) + 1, _vars_idx["weights"] + node_2, 1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23) + 2, _vars_idx["weights"] + node_2, 1.0));
                // \partial_{w_3} \partial r_{23}
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_3, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23), 1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_3, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23) + 1, 1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["weights"] + node_3, _vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23) + 2, 1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23), _vars_idx["weights"] + node_3, 1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23) + 1, _vars_idx["weights"] + node_3, 1.0));
                triplet_vector[_cons_idx["ruling_2"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings"] + 3*_it_he.to_nbdry(heh_23) + 2, _vars_idx["weights"] + node_3, 1.0));
            }
        }

        void get_constrainthess_dev(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            _it_face.reset();

            for(_it_face; _it_face.valid(); _it_face++){
                int i = _it_face.idx();
                int i_nobdry = _it_face.idx_nobdry();

                // differentiate w.r.t. components of vector product
                // (rwr_1)_{y}*(rwr_2)_{z} - (rwr_1)_{z}*(rwr_2)_{y}
                // first term
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 1, _vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 2, 1.0));
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 2, _vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 1, 1.0));
                // second term
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 2, _vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 1, -1.0));
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 1, _vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 2, -1.0));
                // (rwr_1)_{z}*(rwr_2)_{x} - (rwr_1)_{x}*(rwr_2)_{z}
                // first term
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 2, _vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 0, 1.0));
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 0, _vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 2, 1.0));
                // second term
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 0, _vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 2, -1.0));
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry + 1].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 2, _vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 0, -1.0));
                // (rwr_1)_{x}*(rwr_2)_{y} - (rwr_1)_{y}*(rwr_2)_{x}
                // first term
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 0, _vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 1, 1.0));
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 1, _vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 0, 1.0));
                // second term
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 1, _vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 0, -1.0));
                triplet_vector[_cons_idx["dev"] + 3*i_nobdry + 2].push_back(Eigen::Triplet<RealType>(_vars_idx["reweighted_rulings_2"] + 3*i_nobdry + 0, _vars_idx["reweighted_rulings_1"] + 3*i_nobdry + 1, -1.0));
            }
        }

        void get_constrainthess_fair_v(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
           // the fairness constraint is linear, so dont need to do anything here
        }

        void get_constrainthess_fair_n(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            // the fairness constraint is linear, so dont need to do anything here
        }

        void get_constrainthess_fair_r(const VectorType &vars, std::vector<TripletVectorType> triplet_vector) const
        {
            // the fairness constraint is linear, so dont need to do anything here
        }
        // The hessian of the constraint as 3-dimensional tensor, but without
        void apply_all(const VectorType &vars, GenericTensorType &Dest) const
        {
            std::vector<TripletVectorType> triplet_vector(_cons_idx["num_cons"]);
            get_constrainthess_vertex_1_triplet(vars, triplet_vector);
            get_constrainthess_vertex_2_triplet(vars, triplet_vector);
            get_constrainthess_edge_vec_1_triplet(vars, triplet_vector);
            get_constrainthess_edge_vec_2_triplet(vars, triplet_vector);
            get_constrainthess_normal_1_triplet(vars, triplet_vector);
            get_constrainthess_normal_2_triplet(vars, triplet_vector);
            get_constrainthess_normal_3_triplet(vars, triplet_vector);
            get_constrainthess_ruling_0_triplet(vars, triplet_vector);
            get_constrainthess_ruling_1_triplet(vars, triplet_vector);
            get_constrainthess_ruling_2_triplet(vars, triplet_vector);
            get_constrainthess_dev(vars, triplet_vector);
            get_constrainthess_fair_v(vars, triplet_vector);
            get_constrainthess_fair_n(vars, triplet_vector);
            get_constrainthess_fair_r(vars, triplet_vector);

            for(int i = 0; i < _cons_idx["num_cons"]; i++){
                Dest.setFromTriplets(triplet_vector);
            }
        }
    };

/*
template <typename ConfiguratorType>
class VarsIdxReduced{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    const int _num_vertices;
    const int _num_faces;
    const int _num_edges;
    const int _num_halfedges;
    const int _num_bdryHalfEdges;
    const int _num_bdryFaces;
    const QuadMeshTopologySaver &_quadTopol;

    std::map<std::string, int> _idx;

public:
    VarsIdxReduced(const QuadMeshTopologySaver &quadTopol)
    : _num_vertices(quadTopol.getNumVertices()),
      _num_faces(quadTopol.getNumFaces()),
      _num_edges(quadTopol.getNumEdges()),
      _num_halfedges(quadTopol.getNumHalfEdges()),
      _num_bdryHalfEdges(quadTopol.getBdryHalfEdges().size()),
      _num_bdryFaces(quadTopol.getBdryFaces().size()),                  
      _quadTopol(quadTopol)

{
    _idx["vertices"] = 0;
    _idx["dummy_weights"] = 3*_num_vertices;
    _idx["weights"] = 4*_num_vertices;
    _idx["num_dofs"] = 5*_num_vertices;
}

RealType dummy_weight(const VectorType &vars, size_t i) const {
    if(i >= _num_vertices){
        std::cerr << "Index out of bounds for dummy weights"<<std::endl;
    }
    return vars[_idx.at("dummy_weights") + i];
}

VectorType vertex(const VectorType &vars, size_t i) const {
    if(i >= _num_vertices){
        std::cerr << "Index out of bounds for vertices"<<std::endl;
    }
    return vars.segment(_idx.at("vertices") + 3*i,3);
}

RealType vertex_weight(const VectorType &vars, size_t i) const{
    if(i >= _num_vertices){
        std::cerr << "Index out of bounds for vertex weights"<<std::endl;
    }
    return vars[_idx.at("weights") + i];
}
};*/

template<typename ConfiguratorType>
void export_normals_2(const QuadMeshTopologySaver &quadTopol, const VectorType normal_pos, VectorType normals, const std::string &filepath)
{
    std::ofstream stream;
    stream.open(filepath);
    // Write all the normal positions
    for(int i = 0; i < quadTopol.getNumFaces(); i++)
    {
        auto middle_pos = normal_pos.segment(3*i, 3);
        stream << middle_pos[0] << " " << middle_pos[1] << " " << middle_pos[2] << std::endl;
    }
    // Write all the normals
    for(int i = 0; i < quadTopol.getNumFaces(); i++)
    {
        auto normal = normals.segment(3*i, 3);
        stream << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
    }
    stream.close();
}


template <typename ConfiguratorType>
class ConstraintSqrdReduced : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    private:
        const QuadMeshTopologySaver &_quadTopol;
        size_t _num_vertices;
        size_t _num_faces;
        size_t _num_edges;
        size_t _num_halfedges;
        mutable SkippingBdryFaceIterator _it_face;

        mutable VectorType normals;
        mutable VectorType rw_rulings;

    public:
        ConstraintSqrdReduced(const QuadMeshTopologySaver &quadTopol):
                _quadTopol(quadTopol),
                _num_vertices(quadTopol.getNumVertices()),
                _num_faces(quadTopol.getNumFaces()),
                _num_edges(quadTopol.getNumEdges()),
                _num_halfedges(2*_num_edges),
                _it_face(_quadTopol)
        {
            normals.resize(3*_num_faces);
            normals.setZero();
        }
    
        void apply(const VectorType &vertices, RealType &Dest) const override
        {

            Dest = 0.0;

            normals.setZero();

            // Iterate over all the faces
            for(int i = 0; i < _num_faces; i++)
            {

                // Get the indices of the vertices of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                // NEW 
                // Get the weights of the vertices
                
                RealType w_0 = 0.5;//_vars_idx.vertex_weight(vars, node_0);
                RealType w_1 = 0.5;//_vars_idx.vertex_weight(vars, node_1);
                RealType w_2 = 0.5;//_vars_idx.vertex_weight(vars, node_2);
                RealType w_3 = 0.5;//_vars_idx.vertex_weight(vars, node_3);

                VectorType v_0 = vertices.segment(3*node_0, 3);
                VectorType v_1 = vertices.segment(3*node_1, 3);
                VectorType v_2 = vertices.segment(3*node_2, 3);
                VectorType v_3 = vertices.segment(3*node_3, 3);

                // Calculate the reweighted vertices
                VectorType rw_v_0 = v_0 * w_0;
                VectorType rw_v_1 = v_1 * w_1;
                VectorType rw_v_2 = v_2 * w_2;
                VectorType rw_v_3 = v_3 * w_3;

                // Calculate the reweighted edges
                Eigen::Vector3d rw_e_02 = ((rw_v_2 + rw_v_3) - (rw_v_0 + rw_v_1)).template head<3>();
                Eigen::Vector3d rw_e_13 = ((rw_v_0 + rw_v_3) - (rw_v_1 + rw_v_2)).template head<3>();

                // Calculate the face normals with consistent orientation
                Eigen::Vector3d face_normal = rw_e_02.cross(rw_e_13);
                face_normal = face_normal / face_normal.norm();
                
                normals.segment(3*i, 3) = face_normal;
            }   

            /*
            for(int i = 0; i < _num_faces; i++)
            {
                std::cout<< "Face " << i << ": Normal = " << normals.segment(3*i, 3).transpose() << std::endl;
            }*/

            // Next, compute the reweighted rulings
            _it_face.reset();
            for(_it_face; _it_face.valid(); _it_face++)
            {
                int i = _it_face.idx();
                int i_nbdry = _it_face.idx_nobdry();

                // First, iterate over edges of the face
                int he_0 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int he_1 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int he_2 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int he_3 = _quadTopol.getHalfEdgeOfQuad(i,3);

                int face_0 = _quadTopol.getFaceOfHalfEdge(he_0, 1);
                int face_1 = _quadTopol.getFaceOfHalfEdge(he_1, 1);
                int face_2 = _quadTopol.getFaceOfHalfEdge(he_2, 1);
                int face_3 = _quadTopol.getFaceOfHalfEdge(he_3, 1);

                VectorType ruling_01, ruling_12, ruling_23, ruling_30;
                
                Eigen::Vector3d n_i      = normals.segment(3 * i, 3).template head<3>();
                Eigen::Vector3d n_face0  = normals.segment(3 * face_0, 3).template head<3>();
                Eigen::Vector3d n_face1  = normals.segment(3 * face_1, 3).template head<3>();
                Eigen::Vector3d n_face2  = normals.segment(3 * face_2, 3).template head<3>();
                Eigen::Vector3d n_face3  = normals.segment(3 * face_3, 3).template head<3>();

                ruling_01 = n_i.cross(n_face0);
                ruling_12 = n_i.cross(n_face1);
                ruling_23 = n_i.cross(n_face2);
                ruling_30 = n_i.cross(n_face3);

                // Calculate the reweighted rulings
                Eigen::Vector3d rw_ruling_13 = ruling_12 - ruling_30;
                Eigen::Vector3d rw_ruling_02 = ruling_01 - ruling_23;

                // Now, finally calculate the developability constraint
                VectorType rw_rulings_cross = rw_ruling_13.cross(rw_ruling_02);
                Dest += rw_rulings_cross.dot(rw_rulings_cross);
            }
        }
};

FullMatrixType operator%(const Eigen::Vector3d &v, const FullMatrixType &v_2)
{
    assert(v_2.rows() == 3 && v_2.cols() == 3);
    FullMatrixType result(3,3);
    Eigen::Vector3d v_2_col0 = v_2.col(0);
    Eigen::Vector3d v_2_col1 = v_2.col(1);
    Eigen::Vector3d v_2_col2 = v_2.col(2);
    result.setZero();
    result.col(0) = v.cross(v_2_col0);
    result.col(1) = v.cross(v_2_col1);
    result.col(2) = v.cross(v_2_col2);
    return result;
}

FullMatrixType operator%(const FullMatrixType &v_2, const Eigen::Vector3d &v)
{

    assert((v_2.rows() == 3) && (v_2.cols() % 3 == 0) );
    FullMatrixType result(3,3);
    Eigen::Vector3d v_2_col0 = v_2.col(0);
    Eigen::Vector3d v_2_col1 = v_2.col(1);
    Eigen::Vector3d v_2_col2 = v_2.col(2);
    result.setZero();
    result.col(0) = v_2_col0.cross(v);
    result.col(1) = v_2_col1.cross(v);
    result.col(2) = v_2_col2.cross(v);
    return result;
}

FullMatrixType operator&(const Eigen::Vector3d &v, const FullMatrixType &v_2)
{
    assert((v_2.rows() == 3) && (v_2.cols() % 3 == 0) );

    FullMatrixType result(3, v_2.cols());
    result.setZero();

    size_t num_vertex_cols = v_2.cols() / 3;
    for(int i = 0; i < num_vertex_cols; i++)
    {
        Eigen::Vector3d v_2_col0 = v_2.col(3*i);
        Eigen::Vector3d v_2_col1 = v_2.col(3*i + 1);
        Eigen::Vector3d v_2_col2 = v_2.col(3*i + 2);

        result.col(3*i)     = v.cross(v_2_col0);
        result.col(3*i + 1) = v.cross(v_2_col1);
        result.col(3*i + 2) = v.cross(v_2_col2);
    }
    return result;
}

FullMatrixType operator&(const FullMatrixType &v_2, const Eigen::Vector3d &v)
{
    assert((v_2.rows() == 3) && (v_2.cols() % 3 == 0) );
    assert(v_2.cols() % 3 == 0);

    FullMatrixType result(3, v_2.cols());
    result.setZero();

    size_t num_vertex_cols = v_2.cols() / 3;
    for(int i = 0; i < num_vertex_cols; i++)
    {
        Eigen::Vector3d v_2_col0 = v_2.col(3*i);
        Eigen::Vector3d v_2_col1 = v_2.col(3*i + 1);
        Eigen::Vector3d v_2_col2 = v_2.col(3*i + 2);

        result.col(3*i)     = v_2_col0.cross(v);
        result.col(3*i + 1) = v_2_col1.cross(v);
        result.col(3*i + 2) = v_2_col2.cross(v);
    }
    return result;
}

template <typename ConfiguratorType>
class ConstraintSqrdReducedGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;
        typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

        const QuadMeshTopologySaver &_quadTopol;
        size_t _num_vertices;
        size_t _num_faces;
        size_t _num_edges;
        size_t _num_halfedges;
        size_t _num_bdryHalfEdges;
        size_t _num_bdryFaces;

        //VarsIdxReduced<ConfiguratorType> &_vars_idx;
        mutable SkippingBdryFaceIterator _it_face;

        mutable VectorType normals;
        mutable FullMatrixType dnormal_dv;

    public:
        ConstraintSqrdReducedGradient(const QuadMeshTopologySaver &quadTopol):
            _quadTopol(quadTopol),
            _num_vertices(quadTopol.getNumVertices()),
            _num_faces(quadTopol.getNumFaces()),
            _num_edges(quadTopol.getNumEdges()),
            _num_halfedges(2*_num_edges),
            _num_bdryHalfEdges(quadTopol.getBdryHalfEdges().size()),
            _num_bdryFaces(quadTopol.getBdryFaces().size()),
            _it_face(_quadTopol)
            {
                normals.resize(3*_num_faces);
                normals.setZero();

                dnormal_dv.resize(3*_num_faces, 3*_num_vertices);
                dnormal_dv.setZero();
            }

        void apply(const VectorType &vertices,  VectorType &Dest) const override
        {
            size_t num_dofs = 3*_quadTopol.getNumVertices();

            //Dest.resize(3*(_num_faces - _quadTopol.getNumBdryFaces()), num_dofs);

            Dest.resize(num_dofs);
            Dest.setZero();

            // We only have a derivative of developability constraint w.r.t. vertices
            
            // First, calculate all normal vectors
            normals.setZero();
            dnormal_dv.setZero();

            // Iterate over all the faces
            for(int i = 0; i < _num_faces; i++)
            {

               // std::cout<< "Face " << i << ": ";

                // Get the indices of the vertices of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                // Get the weights of the vertices
                RealType w_0 = 0.5;//_vars_idx.vertex_weight(vars, node_0);
                RealType w_1 = 0.5;//_vars_idx.vertex_weight(vars, node_1);
                RealType w_2 = 0.5;//_vars_idx.vertex_weight(vars, node_2);
                RealType w_3 = 0.5;//_vars_idx.vertex_weight(vars, node_3);

                Eigen::Vector3d v_0 = vertices.segment(3*node_0, 3).template head<3>();
                Eigen::Vector3d v_1 = vertices.segment(3*node_1, 3).template head<3>();
                Eigen::Vector3d v_2 = vertices.segment(3*node_2, 3).template head<3>();
                Eigen::Vector3d v_3 = vertices.segment(3*node_3, 3).template head<3>();

                /*
                std::cout<<"v_0: "<<v_0.transpose()<<std::endl;
                std::cout<<"v_1: "<<v_1.transpose()<<std::endl;
                std::cout<<"v_2: "<<v_2.transpose()<<std::endl;
                std::cout<<"v_3: "<<v_3.transpose()<<std::endl;
                */

                // Calculate the reweighted vertices
                Eigen::Vector3d rw_v_0 = v_0 * w_0;
                Eigen::Vector3d rw_v_1 = v_1 * w_1;
                Eigen::Vector3d rw_v_2 = v_2 * w_2;
                Eigen::Vector3d rw_v_3 = v_3 * w_3;

                // Calculate the reweighted edges
                Eigen::Vector3d rw_e_02 = (rw_v_2 + rw_v_3) - (rw_v_0 + rw_v_1);
                Eigen::Vector3d rw_e_13 = (rw_v_0 + rw_v_3) - (rw_v_1 + rw_v_2);

                /*
                std::cout<<"rw_e_02: "<<rw_e_02.transpose()<<std::endl;
                std::cout<<"rw_e_13: "<<rw_e_13.transpose()<<std::endl;
                */

                // Calculate the face normals with consistent orientation
                Eigen::Vector3d face_normal = rw_e_02.cross(rw_e_13);
                face_normal = face_normal / face_normal.norm();

                //std::cout<<"face_normal: "<<face_normal.transpose()<<std::endl;
                
                normals[3*i] = face_normal[0];
                normals[3*i + 1] = face_normal[1];
                normals[3*i + 2] = face_normal[2];

                // Calculate the normal derivatives for the face
                FullMatrixType A(3,3);
                A(0,0) = rw_e_02[0];
                A(0,1) = rw_e_02[1];
                A(0,2) = rw_e_02[2];
                A(1,0) = rw_e_13[0];
                A(1,1) = rw_e_13[1];
                A(1,2) = rw_e_13[2];
                A(2,0) = 2*face_normal[0];
                A(2,1) = 2*face_normal[1];
                A(2,2) = 2*face_normal[2];

                if(!v_0.allFinite())
                {
                    std::cout << "v_0 contains NaN or Inf" << std::endl;
                }
                if(!v_1.allFinite())
                {
                    std::cout << "v_1 contains NaN or Inf" << std::endl;
                }
                if(!v_2.allFinite())
                {
                    std::cout << "v_2 contains NaN or Inf" << std::endl;
                }
                if(!v_3.allFinite())
                {
                    std::cout << "v_3 contains NaN or Inf" << std::endl;
                }

                /*
                for(int k= 0; k < 3; k++)
                {
                    for(int m = 0; m < 3; m++)
                    {
                        std::cout<< "A(" << k << "," << m << ") = " << A(k,m) << " ";
                    }
                    std::cout<<std::endl;
                }*/

                //std::cout<<"Face i"<<i<<": "<<A.determinant()<<std::endl;
                if (!rw_e_02.allFinite()) {
                    std::cout << "rw_e_02 contains NaN or Inf" << std::endl;
                }
                if (!rw_e_13.allFinite()) {
                    std::cout << "rw_e_02 contains NaN or Inf" << std::endl;
                }
                if (!face_normal.allFinite()) {
                    std::cout << "rw_e_02 contains NaN or Inf" << std::endl;
                }

                // Calculate dn_dv_0, ..., dn_dv_3
                // First, calculate dn_dv_0
                FullMatrixType rhs_0(3,3);
                rhs_0.setZero();
                rhs_0(0,0) = -2*face_normal[0];
                rhs_0(0,1) = -2*face_normal[1];
                rhs_0(0,2) = -2*face_normal[2];
                rhs_0(1,0) = 2*face_normal[0];
                rhs_0(1,1) = 2*face_normal[1];
                rhs_0(1,2) = 2*face_normal[2];

                FullMatrixType dn_dv_0 = A.colPivHouseholderQr().solve(-rhs_0);

                // Now, calculate dn_dv_1
                FullMatrixType rhs_1(3,3);
                rhs_1.setZero();
                rhs_1(0,0) = -2*face_normal[0];
                rhs_1(0,1) = -2*face_normal[1];
                rhs_1(0,2) = -2*face_normal[2];
                rhs_1(1,0) = -2*face_normal[0];
                rhs_1(1,1) = -2*face_normal[1];
                rhs_1(1,2) = -2*face_normal[2];

                // Calculate dn_dv_2
                FullMatrixType dn_dv_1 = A.colPivHouseholderQr().solve(-rhs_1);

                // Now, calculate dn_dv_2
                FullMatrixType rhs_2(3,3);
                rhs_2.setZero();
                rhs_2(0,0) = 2*face_normal[0];
                rhs_2(0,1) = 2*face_normal[1];
                rhs_2(0,2) = 2*face_normal[2];
                rhs_2(1,0) = -2*face_normal[0];
                rhs_2(1,1) = -2*face_normal[1];
                rhs_2(1,2) = -2*face_normal[2];
                FullMatrixType dn_dv_2 = A.colPivHouseholderQr().solve(-rhs_2);

                // Now, calculate dn_dv_3
                FullMatrixType rhs_3(3,3);
                rhs_3.setZero();
                rhs_3(0,0) = 2*face_normal[0];
                rhs_3(0,1) = 2*face_normal[1];
                rhs_3(0,2) = 2*face_normal[2];
                rhs_3(1,0) = 2*face_normal[0];
                rhs_3(1,1) = 2*face_normal[1];
                rhs_3(1,2) = 2*face_normal[2];
                FullMatrixType dn_dv_3 = A.colPivHouseholderQr().solve(-rhs_3);

                dnormal_dv.block(3*i, 3*node_0, 3, 3) = dn_dv_0;
                dnormal_dv.block(3*i, 3*node_1, 3, 3) = dn_dv_1;
                dnormal_dv.block(3*i, 3*node_2, 3, 3) = dn_dv_2;
                dnormal_dv.block(3*i, 3*node_3, 3, 3) = dn_dv_3;
            }   

            if(dnormal_dv.hasNaN())
            {
                std::cerr << "dnormal_dv contains NaN values!" << std::endl;
            }

            _it_face.reset();
            for(_it_face; _it_face.valid(); _it_face++)
            {
                int i = _it_face.idx();
                int i_nbdry = _it_face.idx_nobdry();

                // Calculate the derivatives of all adjacent ruling vectors
                int he_0 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int he_1 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int he_2 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int he_3 = _quadTopol.getHalfEdgeOfQuad(i,3);

                int face_0 = _quadTopol.getFaceOfHalfEdge(he_0, 1);
                int face_1 = _quadTopol.getFaceOfHalfEdge(he_1, 1);
                int face_2 = _quadTopol.getFaceOfHalfEdge(he_2, 1);
                int face_3 = _quadTopol.getFaceOfHalfEdge(he_3, 1);

                // All nodes of the original face
                int node_i_0 = _quadTopol.getNodeOfQuad(i, 0);
                int node_i_1 = _quadTopol.getNodeOfQuad(i, 1);
                int node_i_2 = _quadTopol.getNodeOfQuad(i, 2);
                int node_i_3 = _quadTopol.getNodeOfQuad(i, 3);

                // Nodes of face zero without nodes of face i
                int start_node_0, start_node_1, start_node_2, start_node_3;
                for(int j = 0; j < 4; j++)
                {
                    if(_quadTopol.getNodeOfQuad(face_0, j) == node_i_0)
                    {
                        start_node_0 = j;
                    }
                    if(_quadTopol.getNodeOfQuad(face_1, j) == node_i_1)
                    {
                        start_node_1 = j;
                    }
                    if(_quadTopol.getNodeOfQuad(face_2, j) == node_i_2)
                    {
                        start_node_2 = j;
                    }
                    if(_quadTopol.getNodeOfQuad(face_3, j) == node_i_3)
                    {
                        start_node_3 = j;
                    }
                }

                int sign_face_0 = 1;
                if(_quadTopol.getNodeOfQuad(face_0, ((start_node_0 + 1)%4 + 4)%4) == node_i_1)
                {
                    sign_face_0 = -1;
                }

                // ensure that the node indices always stay in [0,3]
                int node_0_0 = _quadTopol.getNodeOfQuad(face_0, ((start_node_0 + sign_face_0*1)%4 + 4)%4);
                int node_0_1 = _quadTopol.getNodeOfQuad(face_0, ((start_node_0 + sign_face_0*2)%4 + 4)%4);

                int sign_face_1 = 1;
                if(_quadTopol.getNodeOfQuad(face_1, ((start_node_1 + 1)%4 + 4)%4) == node_i_2)
                {
                    sign_face_1 = -1;
                }

                // Nodes of face one without nodes of face i
                int node_1_0 = _quadTopol.getNodeOfQuad(face_1, ((start_node_1 + sign_face_1*1)%4 + 4)%4);
                int node_1_1 = _quadTopol.getNodeOfQuad(face_1, ((start_node_1 + sign_face_1*2)%4 + 4)%4);

                int sign_face_2 = 1;
                if(_quadTopol.getNodeOfQuad(face_2, ((start_node_2 + 1)%4 + 4)%4) == node_i_3)
                {
                    sign_face_2 = -1;
                }

                // All nodes of face two
                int node_2_0 = _quadTopol.getNodeOfQuad(face_2, ((start_node_2 + sign_face_2*1)%4 + 4)%4);
                int node_2_1 = _quadTopol.getNodeOfQuad(face_2, ((start_node_2 + sign_face_2*2)%4 + 4)%4);

                int sign_face_3 = 1;
                if(_quadTopol.getNodeOfQuad(face_3, ((start_node_3 + 1)%4 + 4)%4) == node_i_0)
                {
                    sign_face_3 = -1;
                }

                // All nodes of face three
                int node_3_0 = _quadTopol.getNodeOfQuad(face_3, ((start_node_3 + sign_face_3*1)%4 + 4)%4);
                int node_3_1 = _quadTopol.getNodeOfQuad(face_3, ((start_node_3 + sign_face_3*2)%4 + 4)%4);

                // Halfedge 0. First index is the face index, then the local index.
                // dr_he_0_dv_i_0 -> derivative of halfedge 0 w.r.t. vertex 0 of face i
                Eigen::Vector3d normal_i = normals.segment(3*i, 3).template head<3>();
                Eigen::Vector3d normal_0 = normals.segment(3*face_0, 3).template head<3>();
                Eigen::Vector3d normal_1 = normals.segment(3*face_1, 3).template head<3>();
                Eigen::Vector3d normal_2 = normals.segment(3*face_2, 3). template head<3>();
                Eigen::Vector3d normal_3 = normals.segment(3*face_3, 3).template head<3>();

                // Halfedge 0
                FullMatrixType dr_he_01_dv_i_0 = normal_i%(dnormal_dv.block(3*face_0, 3*node_i_0, 3, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_0, 3,3))%normal_0;

                FullMatrixType dr_he_01_dv_i_1 = normal_i%(dnormal_dv.block(3*face_0, 3*node_i_1, 3,3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_1, 3, 3))%normal_0;

                FullMatrixType dr_he_01_dv_i_2 = //0 + //normals.segment(3*i, 3).cross(dnormal_dv.block(3*face_0, 3*node_i_2, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_2, 3, 3))%normal_0;

                FullMatrixType dr_he_01_dv_i_3 = //0 + //normals.segment(3*i, 3).cross(dnormal_dv.block(3*face_0, 3*node_i_3, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_3, 3, 3))%normal_0;

                FullMatrixType dr_he_01_dv_0_0 = normal_i%(dnormal_dv.block(3*face_0, 3*node_0_0, 3,3));
                                                //0; //normals.segment(3*face_0, 3).cross(dnormal_dv.block(3*i, 3*node_0_0, 3));

                FullMatrixType dr_he_01_dv_0_1 = normal_i%(dnormal_dv.block(3*face_0, 3*node_0_1, 3,3));
                                                //0; //normals.segment(3*face_0, 3).cross(dnormal_dv.block(3*i, 3*node_0_1, 3));

                // Halfedge 1
                FullMatrixType dr_he_12_dv_i_0 = //0 + //normals.segment(3*i, 3).cross(dnormal_dv.block(3*face_1, 3*node_i_0, 3)) + 
                                                (dnormal_dv.block(3*i, 3*node_i_0, 3, 3))%normal_1;

                FullMatrixType dr_he_12_dv_i_1 = normal_i%(dnormal_dv.block(3*face_1, 3*node_i_1, 3, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_1, 3, 3))%normal_1;

                FullMatrixType dr_he_12_dv_i_2 = normal_i%(dnormal_dv.block(3*face_1, 3*node_i_2, 3, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_2, 3, 3))%normal_1;

                FullMatrixType dr_he_12_dv_i_3 = //0 + //normals.segment(3*i, 3).cross(dnormal_dv.block(3*face_1, 3*node_i_3, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_3, 3, 3))%normal_1;

                FullMatrixType dr_he_12_dv_1_0 = normal_i%(dnormal_dv.block(3*face_1, 3*node_1_0, 3, 3));
                                                //0; //normals.segment(3*face_1, 3).cross(dnormal_dv.block(3*i, 3*node_1_0, 3));

                FullMatrixType dr_he_12_dv_1_1 = normal_i%(dnormal_dv.block(3*face_1, 3*node_1_1, 3, 3));
                                                //0; //normals.segment(3*face_1, 3).cross(dnormal_dv.block(3*i, 3*node_1_1, 3));
                
                // Halfedge 2
                FullMatrixType dr_he_23_dv_i_0 = //0 + //normals.segment(3*i, 3).cross(dnormal_dv.block(3*face_2, 3*node_i_0, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_0, 3, 3))%normal_2;

                FullMatrixType dr_he_23_dv_i_1 = //0 + //normals.segment(3*i, 3).cross(dnormal_dv.block(3*face_2, 3*node_i_1, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_1, 3, 3))%normal_2;

                FullMatrixType dr_he_23_dv_i_2 = normal_i%(dnormal_dv.block(3*face_2, 3*node_i_2, 3, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_2, 3, 3))%normal_2;

                FullMatrixType dr_he_23_dv_i_3 = normal_i%(dnormal_dv.block(3*face_2, 3*node_i_3, 3, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_3, 3, 3))%normal_2;

                FullMatrixType dr_he_23_dv_2_0 = normal_i%(dnormal_dv.block(3*face_2, 3*node_2_0, 3, 3));
                                                //0; //normals.segment(3*face_2, 3).cross(dnormal_dv.block(3*i, 3*node_2_0, 3));

                FullMatrixType dr_he_23_dv_2_1 = normal_i%(dnormal_dv.block(3*face_2, 3*node_2_1, 3, 3));
                                                //0; //normals.segment(3*face_2, 3).cross(dnormal_dv.block(3*i, 3*node_2_1, 3));
                
                // Halfedge 3
                FullMatrixType dr_he_30_dv_i_0 = normal_i%(dnormal_dv.block(3*face_3, 3*node_i_0, 3, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_0, 3, 3))%normal_3;

                FullMatrixType dr_he_30_dv_i_1 = //0 + //normals.segment(3*i, 3).cross(dnormal_dv.block(3*face_3, 3*node_i_1, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_1, 3, 3))%normal_3;

                FullMatrixType dr_he_30_dv_i_2 = //0 + //normals.segment(3*i, 3).cross(dnormal_dv.block(3*face_3, 3*node_i_2, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_2, 3, 3))%normal_3;

                FullMatrixType dr_he_30_dv_i_3 = normal_i%(dnormal_dv.block(3*face_3, 3*node_i_3, 3, 3)) +
                                                (dnormal_dv.block(3*i, 3*node_i_3, 3, 3))%normal_3;

                FullMatrixType dr_he_30_dv_3_0 = normal_i%(dnormal_dv.block(3*face_3, 3*node_3_0, 3, 3));
                                                //0; //normals.segment(3*face_3, 3).cross(dnormal_dv.block(3*i, 3*node_3_0, 3));

                FullMatrixType dr_he_30_dv_3_1 = normal_i%(dnormal_dv.block(3*face_3, 3*node_3_1, 3, 3));
                                                //0; //normals.segment(3*face_3, 3).cross(dnormal_dv.block(3*i, 3*node_3_1, 3));

                // Now, we consider the derivatives of the reweighted rulings
                /*
                FullMatrixType drw_13_dv_i_0 = 2*dr_he_12_dv_i_0 - 2*dr_he_30_dv_i_0;
                FullMatrixType drw_13_dv_i_1 = 2*dr_he_12_dv_i_1 - 2*dr_he_30_dv_i_1;
                FullMatrixType drw_13_dv_i_2 = 2*dr_he_12_dv_i_2 - 2*dr_he_30_dv_i_2;
                FullMatrixType drw_13_dv_i_3 = 2*dr_he_12_dv_i_3 - 2*dr_he_30_dv_i_3;
                FullMatrixType drw_13_dv_1_0 = 2*dr_he_12_dv_1_0;
                FullMatrixType drw_13_dv_1_1 = 2*dr_he_12_dv_1_1;
                FullMatrixType drw_13_dv_3_0 = -2*dr_he_30_dv_3_0;
                FullMatrixType drw_13_dv_3_1 = -2*dr_he_30_dv_3_1;

                FullMatrixType drw_02_dv_i_0 = 2*dr_he_01_dv_i_0 - 2*dr_he_23_dv_i_0;
                FullMatrixType drw_02_dv_i_1 = 2*dr_he_01_dv_i_1 - 2*dr_he_23_dv_i_1;
                FullMatrixType drw_02_dv_i_2 = 2*dr_he_01_dv_i_2 - 2*dr_he_23_dv_i_2;
                FullMatrixType drw_02_dv_i_3 = 2*dr_he_01_dv_i_3 - 2*dr_he_23_dv_i_3;
                FullMatrixType drw_02_dv_0_0 = 2*dr_he_01_dv_0_0;
                FullMatrixType drw_02_dv_0_1 = 2*dr_he_01_dv_0_1;
                FullMatrixType drw_02_dv_2_0 = -2*dr_he_23_dv_2_0;
                FullMatrixType drw_02_dv_2_1 = -2*dr_he_23_dv_2_1;
                */

                FullMatrixType dr_he_01(3, num_dofs);
                dr_he_01.setZero();
                dr_he_01.block(0, 3*node_i_0, 3, 3) = dr_he_01_dv_i_0;
                dr_he_01.block(0, 3*node_i_1, 3, 3) = dr_he_01_dv_i_1;
                dr_he_01.block(0, 3*node_i_2, 3, 3) = dr_he_01_dv_i_2;
                dr_he_01.block(0, 3*node_i_3, 3, 3) = dr_he_01_dv_i_3;
                dr_he_01.block(0, 3*node_0_0, 3, 3) = dr_he_01_dv_0_0;
                dr_he_01.block(0, 3*node_0_1, 3, 3) = dr_he_01_dv_0_1;
                
                FullMatrixType dr_he_12(3, num_dofs);
                dr_he_12.setZero();
                dr_he_12.block(0, 3*node_i_0, 3, 3) = dr_he_12_dv_i_0;
                dr_he_12.block(0, 3*node_i_1, 3, 3) = dr_he_12_dv_i_1;
                dr_he_12.block(0, 3*node_i_2, 3, 3) = dr_he_12_dv_i_2;
                dr_he_12.block(0, 3*node_i_3, 3, 3) = dr_he_12_dv_i_3;
                dr_he_12.block(0, 3*node_1_0, 3, 3) = dr_he_12_dv_1_0;
                dr_he_12.block(0, 3*node_1_1, 3, 3) = dr_he_12_dv_1_1;

                FullMatrixType dr_he_23(3, num_dofs);
                dr_he_23.setZero();
                dr_he_23.block(0, 3*node_i_0, 3, 3) = dr_he_23_dv_i_0;
                dr_he_23.block(0, 3*node_i_1, 3, 3) = dr_he_23_dv_i_1;
                dr_he_23.block(0, 3*node_i_2, 3, 3) = dr_he_23_dv_i_2;
                dr_he_23.block(0, 3*node_i_3, 3, 3) = dr_he_23_dv_i_3;
                dr_he_23.block(0, 3*node_2_0, 3, 3) = dr_he_23_dv_2_0;
                dr_he_23.block(0, 3*node_2_1, 3, 3) = dr_he_23_dv_2_1;

                FullMatrixType dr_he_30(3, num_dofs);
                dr_he_30.setZero();
                dr_he_30.block(0, 3*node_i_0, 3, 3) = dr_he_30_dv_i_0;
                dr_he_30.block(0, 3*node_i_1, 3, 3) = dr_he_30_dv_i_1;
                dr_he_30.block(0, 3*node_i_2, 3, 3) = dr_he_30_dv_i_2;
                dr_he_30.block(0, 3*node_i_3, 3, 3) = dr_he_30_dv_i_3;
                dr_he_30.block(0, 3*node_3_0, 3, 3) = dr_he_30_dv_3_0;
                dr_he_30.block(0, 3*node_3_1, 3, 3) = dr_he_30_dv_3_1;

                FullMatrixType drw_13(3, num_dofs);
                drw_13 = dr_he_12 - dr_he_30;

                FullMatrixType drw_02(3, num_dofs);
                drw_02 = dr_he_01 - dr_he_23;

                // Now, finally calculate the developability constraint gradient
                // Compute the reweighted rulings themselves
                Eigen::Vector3d ruling_01 = normal_i.cross(normal_0);
                Eigen::Vector3d ruling_12 = normal_i.cross(normal_1);
                Eigen::Vector3d ruling_23 = normal_i.cross(normal_2);
                Eigen::Vector3d ruling_30 = normal_i.cross(normal_3);

                // Calculate the reweighted rulings
                Eigen::Vector3d rw_ruling_13 = ruling_12 - ruling_30;
                Eigen::Vector3d rw_ruling_02 = ruling_01 - ruling_23;

                Eigen::Vector3d rw_rulings_cross = rw_ruling_13.cross(rw_ruling_02);

                //VectorType store = (factor_1.transpose() * rw_rulings_cross);
                //std::cout<<"store norm: "<<store.norm()<<std::endl;

                Dest += factor_1.transpose() * rw_rulings_cross;

                /*
                Dest.segment(3*node_i_0, 3) += 2*((drw_13_dv_i_0%(rw_ruling_02) + rw_ruling_13%(drw_02_dv_i_0)).transpose()*rw_rulings_cross);
                // Derivative w.r.t. i_1
                Dest.segment(3*node_i_1, 3) += 2*((drw_13_dv_i_1%(rw_ruling_02) + rw_ruling_13%(drw_02_dv_i_1)).transpose()*rw_rulings_cross);
                // Derivative w.r.t. i_2
                Dest.segment(3*node_i_2, 3) += 2*((drw_13_dv_i_2%(rw_ruling_02) + rw_ruling_13%(drw_02_dv_i_2)).transpose()*rw_rulings_cross);
                // Derivative w.r.t. i_3
                Dest.segment(3*node_i_3, 3) += 2*((drw_13_dv_i_3%(rw_ruling_02) + rw_ruling_13%(drw_02_dv_i_3)).transpose()*rw_rulings_cross);
                // Derivative w.r.t. 0_0
                Dest.segment(3*node_0_0, 3) += 2*((rw_ruling_13%(drw_02_dv_0_0)).transpose()*rw_rulings_cross);
                // Derivative w.r.t. 0_1
                Dest.segment(3*node_0_1, 3) += 2*((rw_ruling_13%(drw_02_dv_0_1)).transpose()*rw_rulings_cross);
                // Derivative w.r.t. 1_0
                Dest.segment(3*node_1_0, 3) += 2*((drw_13_dv_1_0%(rw_ruling_02)).transpose()*rw_rulings_cross);
                // Derivative w.r.t. 1_1
                Dest.segment(3*node_1_1, 3) += 2*((drw_13_dv_1_1%(rw_ruling_02)).transpose()*rw_rulings_cross);
                // Derivative w.r.t. 2_0
                Dest.segment(3*node_2_0, 3) += 2*((rw_ruling_13%(drw_02_dv_2_0)).transpose()*rw_rulings_cross);
                // Derivative w.r.t. 2_1
                Dest.segment(3*node_2_1, 3) += 2*((rw_ruling_13%(drw_02_dv_2_1)).transpose()*rw_rulings_cross);
                // Derivative w.r.t. 3_0
                Dest.segment(3*node_3_0, 3) += 2*(((drw_13_dv_3_0)%rw_ruling_02).transpose()*rw_rulings_cross);
                // Derivative w.r.t. 3_1
                Dest.segment(3*node_3_1, 3) += 2*(((drw_13_dv_3_1)%rw_ruling_02).transpose()*rw_rulings_cross);
                */
            }
        }

};

template<typename ConfiguratorType>
class EdgeLengthQuadEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>{
    protected:
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

        const QuadMeshTopologySaver &_quadTopol;
        const VectorType &_inactiveGeometry;

    public:
        EdgeLengthQuadEnergy(const QuadMeshTopologySaver &quadTopol,
                        const VectorType &inactiveGeometry
                    ) : _quadTopol(quadTopol),
                        _inactiveGeometry(inactiveGeometry)
        {}

        void apply(const VectorType &quad_geom, RealType &Dest) const override
        {
            Dest = 0.0;
            // penalize differences in edge lengths
            for(int i = 0; i < _quadTopol.getNumEdges(); i++)
            {
                int node_0 = _quadTopol.getAdjacentNodeOfEdge(i, 0);
                int node_1 = _quadTopol.getAdjacentNodeOfEdge(i, 1);

                // Get the positions of the nodes
                Eigen::Vector3d pos_0 = quad_geom.segment(3*node_0, 3);
                Eigen::Vector3d pos_1 = quad_geom.segment(3*node_1, 3);

                Eigen::Vector3d pos_0_inactive = _inactiveGeometry.segment(3*node_0, 3);
                Eigen::Vector3d pos_1_inactive = _inactiveGeometry.segment(3*node_1, 3);

                // Calculate the squared edge lengths
                RealType edge_length = (pos_0 - pos_1).norm();
                RealType edge_length_inactive = (pos_0_inactive - pos_1_inactive).norm();

                RealType log_diff = (std::log(edge_length) - std::log(edge_length_inactive));
                Dest += log_diff*log_diff;
            }

            // penalize differences in face diagonals
            for(int i = 0; i < _quadTopol.getNumFaces(); i++)
            {
                int node_0 = _quadTopol.getNodeOfQuad(i, 0);
                int node_1 = _quadTopol.getNodeOfQuad(i, 1);
                int node_2 = _quadTopol.getNodeOfQuad(i, 2);
                int node_3 = _quadTopol.getNodeOfQuad(i, 3);

                // Get the positions of the nodes
                Eigen::Vector3d pos_0 = quad_geom.segment(3*node_0, 3);
                Eigen::Vector3d pos_1 = quad_geom.segment(3*node_1, 3);
                Eigen::Vector3d pos_2 = quad_geom.segment(3*node_2, 3);
                Eigen::Vector3d pos_3 = quad_geom.segment(3*node_3, 3);

                Eigen::Vector3d pos_0_inactive = _inactiveGeometry.segment(3*node_0, 3);
                Eigen::Vector3d pos_1_inactive = _inactiveGeometry.segment(3*node_1, 3);
                Eigen::Vector3d pos_2_inactive = _inactiveGeometry.segment(3*node_2, 3);
                Eigen::Vector3d pos_3_inactive = _inactiveGeometry.segment(3*node_3, 3);

                // Calculate the squared face diagonals
                RealType diag_0 = (pos_0 - pos_2).norm();
                RealType diag_1 = (pos_1 - pos_3).norm();

                RealType diag_0_inactive = (pos_0_inactive - pos_2_inactive).norm();
                RealType diag_1_inactive = (pos_1_inactive - pos_3_inactive).norm();

                RealType log_diff_diag_0 = (std::log(diag_0) - std::log(diag_0_inactive));
                RealType log_diff_diag_1 = (std::log(diag_1) - std::log(diag_1_inactive));

                Dest += log_diff_diag_0*log_diff_diag_0 + log_diff_diag_1*log_diff_diag_1;
            }
        }
};


template<typename ConfiguratorType>
class EdgeLengthQuadEnergyGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    protected:
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

        const QuadMeshTopologySaver &_quadTopol;
        const VectorType &_inactiveGeometry;

    public:
        EdgeLengthQuadEnergyGradient(const QuadMeshTopologySaver &quadTopol,
                                const VectorType &inactiveGeometry
                            ) : _quadTopol(quadTopol),
                                _inactiveGeometry(inactiveGeometry)
        {}

        void apply(const VectorType &quad_geom, VectorType &Dest) const override
        {
            Dest.setZero();
            if(Dest.size() != quad_geom.size())
            {
                Dest.resize(quad_geom.size());
            }
            // differentiate penalizations in edge length
            for(int i = 0; i < _quadTopol.getNumEdges(); i++)
            {
                int node_0 = _quadTopol.getAdjacentNodeOfEdge(i, 0);
                int node_1 = _quadTopol.getAdjacentNodeOfEdge(i, 1);

                Eigen::Vector3d pos_0 = quad_geom.segment(3*node_0, 3);
                Eigen::Vector3d pos_1 = quad_geom.segment(3*node_1, 3);

                Eigen::Vector3d pos_0_inactive = _inactiveGeometry.segment(3*node_0, 3);
                Eigen::Vector3d pos_1_inactive = _inactiveGeometry.segment(3*node_1, 3);

                RealType edge_length = (pos_0 - pos_1).norm();
                RealType edge_length_sqrd = (pos_0 - pos_1).squaredNorm();
                RealType edge_length_inactive = (pos_0_inactive - pos_1_inactive).norm();

                RealType log_diff = (std::log(edge_length) - std::log(edge_length_inactive));

                Dest.segment(3*node_0, 3) += 2*log_diff*(1.0/edge_length_sqrd)*(pos_0 - pos_1);
                Dest.segment(3*node_1, 3) += 2*log_diff*(1.0/edge_length_sqrd)*(pos_1 - pos_0);
            }

            // differentiate penalizations in face diagonals
            for(int i = 0; i < _quadTopol.getNumFaces(); i++)
            {
                int node_0 = _quadTopol.getNodeOfQuad(i, 0);
                int node_1 = _quadTopol.getNodeOfQuad(i, 1);
                int node_2 = _quadTopol.getNodeOfQuad(i, 2);
                int node_3 = _quadTopol.getNodeOfQuad(i, 3);

                Eigen::Vector3d pos_0 = quad_geom.segment(3*node_0, 3);
                Eigen::Vector3d pos_1 = quad_geom.segment(3*node_1, 3);
                Eigen::Vector3d pos_2 = quad_geom.segment(3*node_2, 3);
                Eigen::Vector3d pos_3 = quad_geom.segment(3*node_3, 3);

                Eigen::Vector3d pos_0_inactive = _inactiveGeometry.segment(3*node_0, 3);
                Eigen::Vector3d pos_1_inactive = _inactiveGeometry.segment(3*node_1, 3);
                Eigen::Vector3d pos_2_inactive = _inactiveGeometry.segment(3*node_2, 3);
                Eigen::Vector3d pos_3_inactive = _inactiveGeometry.segment(3*node_3, 3);

                // Calculate the squared face diagonals
                RealType diag_0 = (pos_0 - pos_2).norm();
                RealType diag_1 = (pos_1 - pos_3).norm();

                RealType diag_0_inactive = (pos_0_inactive - pos_2_inactive).norm();
                RealType diag_1_inactive = (pos_1_inactive - pos_3_inactive).norm();

                RealType log_diff_diag_0 = (std::log(diag_0) - std::log(diag_0_inactive));
                RealType log_diff_diag_1 = (std::log(diag_1) - std::log(diag_1_inactive));

                Dest.segment(3*node_0, 3) += 2*log_diff_diag_0*(1.0/(diag_0*diag_0))*(pos_0 - pos_2);
                Dest.segment(3*node_2, 3) += 2*log_diff_diag_0*(1.0/(diag_0*diag_0))*(pos_2 - pos_0);
                Dest.segment(3*node_1, 3) += 2*log_diff_diag_1*(1.0/(diag_1*diag_1))*(pos_1 - pos_3);
                Dest.segment(3*node_3, 3) += 2*log_diff_diag_1*(1.0/(diag_1*diag_1))*(pos_3 - pos_1);
            }
        }

};

// Now, write a function to easily export the Gauss map image as plot
// Needs to be able to access 
/*
template<typename ConfiguratorType>
void plot_gauss_image(QuadMeshTopologySaver &quadTopol, VectorView<ConfiguratorType> &view, std::string &filepath)
{
    std::ofstream stream;
    stream.open(filepath);
    // Write all the normals
    for(int i = 0; i < quadTopol.getNumFaces(); i++)
    {
        auto normal = view.face_normal(i);
        stream << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
    }
}

template<typename ConfiguratorType>
void plot_gauss_image(QuadMeshTopologySaver &quadTopol, Constraint<ConfiguratorType> &constraint, std::string &filepath)
{
    std::ofstream stream;
    stream.open(filepath);
    // Write all the normals
    for(int i = 0; i < quadTopol.getNumFaces(); i++)
    {
        auto normal = constraint._view.face_normal(i);
        stream << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
    }
}

*/
template<typename ConfiguratorType>
void export_normals(QuadMeshTopologySaver &quadTopol, VectorType &vars, VarsIdx<ConfiguratorType> &vars_idx, std::string &filepath)
{
    std::ofstream stream;
    stream.open(filepath);
    // We want to calculate the average of face vertices to plot the normals at
    for(int i = 0; i < quadTopol.getNumFaces(); i++)
    {
        auto vertex_0_idx = quadTopol.getNodeOfQuad(i,0);
        auto vertex_1_idx = quadTopol.getNodeOfQuad(i,1);
        auto vertex_2_idx = quadTopol.getNodeOfQuad(i,2);
        auto vertex_3_idx = quadTopol.getNodeOfQuad(i,3);

        auto vertex_0 = vars_idx.vertex(vars, vertex_0_idx);
        auto vertex_1 = vars_idx.vertex(vars,vertex_1_idx);
        auto vertex_2 = vars_idx.vertex(vars, vertex_2_idx);
        auto vertex_3 = vars_idx.vertex(vars, vertex_3_idx);

        auto middle_pos = 0.25*(vertex_0 + vertex_1 + vertex_2 + vertex_3);
        std::cout<<"Middle pos: " << middle_pos[0] << " " << middle_pos[1] << " " << middle_pos[2] << std::endl;
        stream << middle_pos[0] << " " << middle_pos[1] << " " << middle_pos[2] << std::endl;
    }
    // Write all the normals
    for(int i = 0; i < quadTopol.getNumFaces(); i++)
    {
        auto normal = vars_idx.face_normal(vars, i);
        stream << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
    }
}

template<typename ConfiguratorType>
void export_rw_rulings(QuadMeshTopologySaver &quadTopol, VectorType &vars, VarsIdx<ConfiguratorType> &vars_idx, std::string &filepath)
{
    std::ofstream stream;
    stream.open(filepath);

    // Now, for each face, first want again the middle point of the face
    SkippingBdryFaceIterator it_face(quadTopol);
    for(it_face; it_face.valid(); it_face ++)
    {
        int i = it_face.idx();

        auto vertex_0_idx = quadTopol.getNodeOfQuad(i,0);
        auto vertex_1_idx = quadTopol.getNodeOfQuad(i,1);
        auto vertex_2_idx = quadTopol.getNodeOfQuad(i,2);
        auto vertex_3_idx = quadTopol.getNodeOfQuad(i,3);

        auto vertex_0 = vars_idx.vertex(vars, vertex_0_idx);
        auto vertex_1 = vars_idx.vertex(vars, vertex_1_idx);
        auto vertex_2 = vars_idx.vertex(vars, vertex_2_idx);
        auto vertex_3 = vars_idx.vertex(vars, vertex_3_idx);

        auto middle_pos = 0.25*(vertex_0 + vertex_1 + vertex_2 + vertex_3);
        stream << middle_pos[0] << " " << middle_pos[1] << " " << middle_pos[2] << std::endl;
    }

    it_face.reset();

    for(it_face; it_face.valid(); it_face++)
    {
        int i = it_face.idx();
        int i_nobdry = it_face.idx_nobdry();

        auto rw_ruling_1 = vars_idx.reweighted_ruling_1(vars, i_nobdry);
        auto rw_ruling_2 = vars_idx.reweighted_ruling_2(vars, i_nobdry);

        auto rw_ruling_1_vec = rw_ruling_1.template head<3>();
        auto rw_ruling_2_vec = rw_ruling_2.template head<3>();

        stream << rw_ruling_1_vec[0] << " " << rw_ruling_1_vec[1] << " " << rw_ruling_1_vec[2] << std::endl;
        stream << rw_ruling_2_vec[0] << " " << rw_ruling_2_vec[1] << " " << rw_ruling_2_vec[2] << std::endl;
    }
}

template<typename ConfiguratorType>
void export_rulings(QuadMeshTopologySaver &quadTopol, Constraint<ConfiguratorType> &constraint, std::string &filepath)
{
    std::ofstream stream;
    stream.open(filepath);

    // Now, for each halfedge we have a ruling vector. But the ruling vectors for
    // two different halfedges should be identical, just with different signs.
    // so is there some way to export them while respecting the orientation of the mesh?
    SkippingBdryHalfEdgeIterator it_he(quadTopol);
    for(it_he; it_he.valid(); it_he++)
    {
        int i = it_he.idx();
        int i_nobdry = it_he.idx_nobdry();

        auto ruling = constraint._view.ruling(i_nobdry);
        auto ruling_vec = ruling.template head<3>();

        auto fh_0 = quadTopol.getFaceOfHalfEdge(i,0);
        auto fh_1 = quadTopol.getFaceOfHalfEdge(i,1);

        // Need to do this since we cannot access nodes directly from halfedges
        // -> need to go over faces and then edges

        std::vector<int> fh_0_edges;
        std::vector<int> fh_1_edges;
        for(int i = 0; i < 4; i++)
        {
            fh_0_edges.push_back(quadTopol.getEdgeOfQuad(fh_0,i));
            fh_1_edges.push_back(quadTopol.getEdgeOfQuad(fh_1,i));
        }

        std::sort(fh_0_edges.begin(), fh_0_edges.end());
        std::sort(fh_1_edges.begin(), fh_1_edges.end());

        // Check which edge appears in both faces
        int common_edge = -1;
        for(int i = 0; i < 4; i++)
        {
            if(std::find(fh_1_edges.begin(), fh_1_edges.end(), fh_0_edges[i]) != fh_1_edges.end())
            {
                common_edge = fh_0_edges[i];
                break;
            }
        }

        if(common_edge == -1)
        {
            std::cout << "Error: no common edge found for halfedge " << i << std::endl;
            continue;
        }

        // Now, obtain the endpoints of the common edge
        int node_0 = quadTopol.getAdjacentNodeOfEdge(common_edge,0);
        int node_1 = quadTopol.getAdjacentNodeOfEdge(common_edge,1);

        auto node_0_vec = constraint._view.vertex(node_0);
        auto node_1_vec = constraint._view.vertex(node_1);
        auto middle_pos = 0.5*(node_0_vec + node_1_vec);

        stream << middle_pos[0] << " " << middle_pos[1] << " " << middle_pos[2] << std::endl;
    }

    it_he.reset();

    for(it_he; it_he.valid(); it_he++)
    {
        int i = it_he.idx();
        int i_nobdry = it_he.idx_nobdry();

        auto ruling = constraint._view.ruling(i_nobdry);
        auto ruling_vec = ruling.template head<3>();

        stream << ruling_vec[0] << " " << ruling_vec[1] << " " << ruling_vec[2] << std::endl;
    }
}

template<typename ConfiguratorType> 
void export_rw_edges(QuadMeshTopologySaver &quadTopol, Constraint<ConfiguratorType> &constraint, std::string &filepath)
{
    std::ofstream stream;
    stream.open(filepath);

    // First, determine the positions of all the reweighted edges
    for(int i = 0; i < quadTopol.getNumFaces(); i++)
    {
        // First, output the average edge vector e_{0,2} = m_{e2} - m_{e0}
        auto w_0 = constraint._view.vertex_weight(quadTopol.getNodeOfQuad(i,0));
        auto w_1 = constraint._view.vertex_weight(quadTopol.getNodeOfQuad(i,1));
        auto w_2 = constraint._view.vertex_weight(quadTopol.getNodeOfQuad(i,2));
        auto w_3 = constraint._view.vertex_weight(quadTopol.getNodeOfQuad(i,3));

        auto v_0 = constraint._view.vertex(quadTopol.getNodeOfQuad(i,0));
        auto v_1 = constraint._view.vertex(quadTopol.getNodeOfQuad(i,1));
        auto v_2 = constraint._view.vertex(quadTopol.getNodeOfQuad(i,2));
        auto v_3 = constraint._view.vertex(quadTopol.getNodeOfQuad(i,3));

        auto middle_pos = 0.25*(v_0 + v_1 + v_2 + v_3);
        //auto pos_rw_edge_02 = 1.0/(w_0 + w_2)*(1.0/w_0*v_0 + 1.0/w_2*v_2);
        stream << middle_pos[0] << " " << middle_pos[1] << " " << middle_pos[2] << std::endl;
        // Now, output the average edge vector e_{1,3} = m_{e3} - m_{e1}
        //auto pos_rw_edge_13 = 1.0/(w_1 + w_3)*(1.0/w_1*v_1 + 1.0/w_3*v_3);
        stream << middle_pos[0] << " " << middle_pos[1] << " " << middle_pos[2] << std::endl;
    }

    // Second, determine the vectors of the rw edges
    for(int i = 0; i < quadTopol.getNumFaces(); i++)
    {
         // First, output the average edge vector e_{0,2} = m_{e2} - m_{e0}
         auto w_0 = constraint._view.vertex_weight(quadTopol.getNodeOfQuad(i,0));
         auto w_1 = constraint._view.vertex_weight(quadTopol.getNodeOfQuad(i,1));
         auto w_2 = constraint._view.vertex_weight(quadTopol.getNodeOfQuad(i,2));
         auto w_3 = constraint._view.vertex_weight(quadTopol.getNodeOfQuad(i,3));
 
         auto v_0 = constraint._view.vertex(quadTopol.getNodeOfQuad(i,0));
         auto v_1 = constraint._view.vertex(quadTopol.getNodeOfQuad(i,1));
         auto v_2 = constraint._view.vertex(quadTopol.getNodeOfQuad(i,2));
         auto v_3 = constraint._view.vertex(quadTopol.getNodeOfQuad(i,3));

         auto rw_edge_vec_02 = constraint._view.reweighted_edge_1(i);
         //auto rw_edge_vec_02_scaled = 1.0/((w_0 + w_1)*(w_2 + w_3))*rw_edge_vec_02;

         stream << rw_edge_vec_02[0] << " " << rw_edge_vec_02[1] << " " << rw_edge_vec_02[2] << std::endl;

         auto rw_edge_vec_13 = constraint._view.reweighted_edge_2(i);
         //auto rw_edge_vec_13_scaled = 1.0/((w_1 + w_2)*(w_0 + w_3))*rw_edge_vec_13;

         stream << rw_edge_vec_13[0] << " " << rw_edge_vec_13[1] << " " << rw_edge_vec_13[2] << std::endl;
    }
}

template <typename ConfiguratorType>
void print_info(QuadMeshTopologySaver &quadTopol, Constraint<ConfiguratorType> &constraint)
{
    // just print out coordinates and weights corresponding to quad nr. 0
    int node_0 = quadTopol.getNodeOfQuad(0,0);
    int node_1 = quadTopol.getNodeOfQuad(0,1);
    int node_2 = quadTopol.getNodeOfQuad(0,2);
    int node_3 = quadTopol.getNodeOfQuad(0,3);

    auto node_0_vec = constraint._view.vertex(node_0);
    auto node_1_vec = constraint._view.vertex(node_1);
    auto node_2_vec = constraint._view.vertex(node_2);
    auto node_3_vec = constraint._view.vertex(node_3);

    std::cout<<"Node 0: " << node_0_vec[0]<<"; "<<node_0_vec[1]<<"; "<<node_0_vec[2] << std::endl;
    std::cout<<"Node 1: " << node_1_vec[0]<<"; "<<node_1_vec[1]<<"; "<<node_1_vec[2] << std::endl;
    std::cout<<"Node 2: " << node_2_vec[0]<<"; "<<node_2_vec[1]<<"; "<<node_2_vec[2] << std::endl;
    std::cout<<"Node 3: " << node_3_vec[0]<<"; "<<node_3_vec[1]<<"; "<<node_3_vec[2] << std::endl;

    // Now, print out the weights of the vertices
    auto w_0 = constraint._view.vertex_weight(node_0);
    auto w_1 = constraint._view.vertex_weight(node_1);
    auto w_2 = constraint._view.vertex_weight(node_2);
    auto w_3 = constraint._view.vertex_weight(node_3);

    std::cout<<"Weight 0: " << w_0 << std::endl;
    std::cout<<"Weight 1: " << w_1 << std::endl;
    std::cout<<"Weight 2: " << w_2 << std::endl;
    std::cout<<"Weight 3: " << w_3 << std::endl;

    // Now, print the corresponding reweighted edges of the quad
    int faceIdx = 0;
    auto rw_edge_1 = constraint._view.reweighted_edge_1(faceIdx);
    auto rw_edge_2 = constraint._view.reweighted_edge_2(faceIdx);

    std::cout<<"Reweighted edge 1: " << rw_edge_1[0]<<"; "<<rw_edge_1[1]<<"; "<<rw_edge_1[2] << std::endl;
    std::cout<<"Reweighted edge 2: " << rw_edge_2[0]<<"; "<<rw_edge_2[1]<<"; "<<rw_edge_2[2] << std::endl;

    // Now, lets do the reweighted vertices:
    auto rw_vertex_0 = constraint._view.reweighted_vertex(node_0);
    auto rw_vertex_1 = constraint._view.reweighted_vertex(node_1);
    auto rw_vertex_2 = constraint._view.reweighted_vertex(node_2);
    auto rw_vertex_3 = constraint._view.reweighted_vertex(node_3);

    std::cout<<"Reweighted vertex 0: " << rw_vertex_0[0]<<"; "<<rw_vertex_0[1]<<"; "<<rw_vertex_0[2] << std::endl;
    std::cout<<"Reweighted vertex 1: " << rw_vertex_1[0]<<"; "<<rw_vertex_1[1]<<"; "<<rw_vertex_1[2] << std::endl;
    std::cout<<"Reweighted vertex 2: " << rw_vertex_2[0]<<"; "<<rw_vertex_2[1]<<"; "<<rw_vertex_2[2] << std::endl;
    std::cout<<"Reweighted vertex 3: " << rw_vertex_3[0]<<"; "<<rw_vertex_3[1]<<"; "<<rw_vertex_3[2] << std::endl;

    // Now, calculate what reweighted edges should be
    auto rw_edge_1_calculated = -(w_0 + w_1)*(w_2*node_2_vec + w_3*node_3_vec) + (w_2 + w_3)*(w_1*node_1_vec + w_0*node_0_vec);
    auto rw_edge_2_calculated = -(w_1 + w_2)*(w_0*node_0_vec + w_3*node_3_vec) + (w_0 + w_3)*(w_1*node_1_vec + w_2*node_2_vec);

    std::cout<<"Reweighted edge calculated 1: "<< rw_edge_1_calculated[0]<<"; "<<rw_edge_1_calculated[1]<<"; "<<rw_edge_1_calculated[2] << std::endl;
    std::cout<<"Reweighted edge calculated 2: "<< rw_edge_2_calculated[0]<<"; "<<rw_edge_2_calculated[1]<<"; "<<rw_edge_2_calculated[2] << std::endl;

    // check are the constraints actually fulfilled well at this point ?
}