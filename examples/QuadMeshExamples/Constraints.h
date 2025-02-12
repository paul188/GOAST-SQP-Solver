#include <goast/Quads.h>
#include <goast/QuadMesh/QuadTopology.h>
#include <cmath>
#include <goast/SQP/Utils/SparseMat.h>

template<typename ConfiguratorType>
struct constraint_weights{
    typedef typename ConfiguratorType::RealType RealType;

    RealType vertex_1 = 1.0;
    RealType vertex_2 = 1.0;
    RealType bdry_opt = 1.0;
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
        else if (key == "bdry_opt") return bdry_opt;
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
};

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

 /* view of the input vector during optimization, allows easy get-access */
template<typename ConfiguratorType>
class VectorView{
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
        VectorType _vec;
        const QuadMeshTopologySaver &_quadTopol;
        // map to convert halfedge index to non boundary halfedge index
        std::map<int,int> _heh_idx_to_nbdry;

    public:
        std::map<std::string, size_t> _idx;

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

        // VectorView with uninitialized vector
        VectorView(const QuadMeshTopologySaver &quadTopol, std::optional<const VectorType> vec = std::nullopt)
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


        void set_vector(const VectorType &vec){
            if(vec.size() != _idx["num_dofs"]){
                std::cerr << "Vector size does not match the number of dofs"<<std::endl;
            }
            _vec = vec;
        }

        RealType dummy_weight(size_t i){
            if(i >= _num_vertices){
                std::cerr << "Index out of bounds for dummy weights"<<std::endl;
            }
            return _vec[_idx["dummy_weights"] + i];
        }

        VectorType vertex(size_t i){
            if(i >= _num_vertices){
                std::cerr << "Index out of bounds for vertices"<<std::endl;
            }
            return _vec.segment(_idx["vertices"] + 3*i,3);
        }

        RealType vertex_weight(size_t i){
            if(i >= _num_vertices){
                std::cerr << "Index out of bounds for vertex weights"<<std::endl;
            }
            return _vec[_idx["weights"] + i];
        }

        VectorType reweighted_vertex(size_t i){
            if(i >= _num_vertices){
                std::cerr << "Index out of bounds for reweighted vertices"<<std::endl;
            }
            return _vec.segment(_idx["reweighted_vertices"] + 3*i,3);
        }

        VectorType reweighted_edge_1(size_t i){
            if(i >= _num_faces){
                std::cerr << "Index out of bounds for reweighted edges 1"<<std::endl;
            }
            return _vec.segment(_idx["reweighted_edges_1"] + 3*i, 3);
        }

        VectorType reweighted_edge_2(size_t i){
            if(i >= _num_faces){
                std::cerr << "Index out of bounds for reweighted edges 2"<<std::endl;
            }
            return _vec.segment(_idx["reweighted_edges_2"] + 3*i, 3);
        }

        VectorType face_normal(size_t i){
            if(i >= _num_faces){
                std::cerr << "Index out of bounds for face normals"<<std::endl;
            }
            return _vec.segment(_idx["normals"] + 3*i, 3);
        }

        VectorType ruling(size_t i){
            if(i >= _num_halfedges){
                std::cerr << "Index out of bounds for rulings"<<std::endl;
            }
            return _vec.segment(_idx["rulings"] + 3*i,3);
        }

        VectorType reweighted_ruling_1(size_t i){
            if(i >= _num_faces){
                std::cerr << "Index out of bounds for reweighted rulings 1"<<std::endl;
            }
            return _vec.segment(_idx["reweighted_rulings_1"] + 3*i, 3);
        }

        VectorType reweighted_ruling_2(size_t i){
            if(i >= _num_faces){
                std::cerr << "Index out of bounds for reweighted rulings 1"<<std::endl;
            }
            return _vec.segment(_idx["reweighted_rulings_2"] + 3*i, 3);
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

/* storage for indices of constraints */
template<typename ConfiguratorType>
class ConstraintView{
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;
        // used to store the boundary indices and corresponding values
        typedef typename std::optional<std::pair<std::vector<int>, VectorType>> OptionalBoundaryData;

        std::map<std::string,size_t> _cons_idx;
    public:
        ConstraintView(StripHandler<ConfiguratorType> &stripHandle, const QuadMeshTopologySaver &quadTopol, OptionalBoundaryData bdryData = std::nullopt)
        {
            
            size_t num_bdryOpt = bdryData.has_value() ? bdryData.value().first.size() : 0;

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
            _cons_idx["bdry_opt"] = _cons_idx["vertex_2"] + num_constraints_vert_2;
            _cons_idx["edge_vec_1"] = _cons_idx["bdry_opt"] + num_bdryOpt;
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

        void extend_weights(const constraint_weights<ConfiguratorType> &weights, MatrixType &Dest)
        {
            size_t num_constraints = _cons_idx["num_cons"];
            if(Dest.rows() != num_constraints || Dest.cols() != num_constraints)
            {
                Dest.resize(num_constraints,num_constraints);
            }

            std::vector<Eigen::Triplet<RealType>> triplets;
            triplets.reserve(num_constraints);

            std::vector<std::pair<std::string, std::string>> constraints = {
                {"vertex_1", "vertex_2"}, {"vertex_2", "bdry_opt"}, {"bdry_opt", "edge_vec_1"},
                {"edge_vec_1", "edge_vec_2"}, {"edge_vec_2", "normal_1"}, {"normal_1", "normal_2"},
                {"normal_2", "normal_3"}, {"normal_3", "ruling_0"}, {"ruling_0", "ruling_1"},
                {"ruling_1", "ruling_2"}, {"ruling_2", "dev"}, {"dev", "fair_v"},
                {"fair_v", "fair_n"}, {"fair_n", "fair_r"}, {"fair_r", "num_cons"}
            };

            for (const auto& [start, end] : constraints) {
                for (int i = _cons_idx[start]; i < _cons_idx[end]; i++) {
                    triplets.emplace_back(i, i, weights[start]);  
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
    // used to store the boundary indices and corresponding values
    typedef typename std::optional<std::pair<std::vector<int>, VectorType>> OptionalBoundaryData;

    private:
        const QuadMeshTopologySaver &_quadTopol;
        size_t _num_vertices;
        size_t _num_faces;
        size_t _num_edges;
        size_t _num_halfedges;
        mutable SkippingBdryFaceIterator _it_face;
        mutable SkippingBdryHalfEdgeIterator _it_he;
        mutable StripHandler<ConfiguratorType> _stripHandle;

        mutable std::map<std::string,size_t> _cons_idx;

        // Serves to access the current vector of input variables at the right positions
        mutable VectorView<ConfiguratorType> _view;
        // Serves to write the constraints into the vector
        // as well as the constraint gradients into the matrix
        mutable ConstraintView<ConfiguratorType> _cons_view;

        // boundary data, i.e. boundary indices in _bdryMask
        // and the prescribed point vectors in _bdryVals
        // the boundary mask is ordered as: first x, then y, then z. All x indices ordered according to point indices.
        OptionalBoundaryData _bdryData;

    public:
        Constraint(const QuadMeshTopologySaver &quadTopol, 
                   StripHandler<ConfiguratorType> &stripHandle,
                   OptionalBoundaryData bdryData = std::nullopt):
              _quadTopol(quadTopol),
              _num_vertices(quadTopol.getNumVertices()),
              _num_faces(quadTopol.getNumFaces()),
              _num_edges(quadTopol.getNumEdges()),
              _num_halfedges(2*_num_edges),
              _stripHandle(stripHandle),
              _cons_view(_stripHandle, _quadTopol, bdryData),
              _view(_quadTopol),
              _it_face(_quadTopol),
              _it_he(_quadTopol),
              _bdryData(bdryData)
        {
            _cons_idx = _cons_view._cons_idx;

            if(_bdryData.has_value())
            {
                if(3*_bdryData.value().first.size() != _bdryData.value().second.size())
                {
                    std::cerr << "Boundary data does not match in size"<<std::endl;
                }
            }
        }
    
    void get_constraint_vertex_1(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _view.set_vector(vars);

        for(int i = 0; i < _num_vertices; i++){
            Dest[_cons_idx["vertex_1"] + i] = _view.vertex_weight(i) - (_view.dummy_weight(i)*_view.dummy_weight(i)) - 1.0;
        }
    }

    void get_constraint_vertex_2(const VectorType &vars, VectorType &Dest)const {

        _view.set_vector(vars);

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _num_vertices; i++){
            Dest.segment(_cons_idx["vertex_2"] + 3*i, 3) = _view.reweighted_vertex(i) - _view.vertex_weight(i)*_view.vertex(i); 
        }

    }

    void get_constraint_bdry_opt(const VectorType &vars, VectorType &Dest) const{
        
        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        if(!_bdryData.has_value())
        {
            // boundary mask is not set, so just return
            return; 
        }

        _view.set_vector(vars);

        for(int i = 0; i < _bdryData.value().first.size(); i++){
            auto vertex_now = _view.vertex(_bdryData.value().first[i]);
            auto vertex_pos = _bdryData.value().second.segment(3*i,3);
            auto vec = vertex_now - vertex_pos;
            Dest[_cons_idx["bdry_opt"] + i] = vec.dot(vec);
        }
    }

    void get_constraint_edge_vec_1(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            // First, obtain first vertex idces of the face
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w_0 = _view.vertex_weight(node_0);
            RealType w_1 = _view.vertex_weight(node_1);
            RealType w_2 = _view.vertex_weight(node_2);
            RealType w_3 = _view.vertex_weight(node_3);

            Dest.segment(_cons_idx["edge_vec_1"]+3*i,3) = _view.reweighted_edge_1(i) - (w_0 + w_1)*(_view.reweighted_vertex(node_2)+_view.reweighted_vertex(node_3)) + (w_2 + w_3)*(_view.reweighted_vertex(node_1) + _view.reweighted_vertex(node_0));
        }
    }

    void get_constraint_edge_vec_2(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            // First, obtain first vertex idces of the face
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w_0 = _view.vertex_weight(node_0);
            RealType w_1 = _view.vertex_weight(node_1);
            RealType w_2 = _view.vertex_weight(node_2);
            RealType w_3 = _view.vertex_weight(node_3);

            Dest.segment(_cons_idx["edge_vec_2"] + 3*i,3) = _view.reweighted_edge_2(i) - (w_1 + w_2)*(_view.reweighted_vertex(node_3)+_view.reweighted_vertex(node_0)) + (w_0 + w_3)*(_view.reweighted_vertex(node_1) + _view.reweighted_vertex(node_2));
        }
    }

    void get_constraint_normal_1(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            Dest[_cons_idx["normal_1"] + i] = _view.face_normal(i).dot(_view.reweighted_edge_1(i));
        }
    }

    void get_constraint_normal_2(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            Dest[_cons_idx["normal_2"] + i] = _view.face_normal(i).dot(_view.reweighted_edge_2(i));
        }
    }

    void get_constraint_normal_3(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            Dest[_cons_idx["normal_3"] + i] = _view.face_normal(i).dot(_view.face_normal(i)) - 1.0;
        }
    }

    void get_constraint_ruling_0(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        _view.set_vector(vars);

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _it_he.reset();

        for(_it_he; _it_he.valid(); ++_it_he){
            int i = _it_he.idx();
            int i_nobdry = _it_he.idx_nobdry();

            int face_1 = _quadTopol.getFaceOfHalfEdge(i,0);
            int face_2 = _quadTopol.getFaceOfHalfEdge(i,1);

            Dest.segment(_cons_idx["ruling_0"] + 3*i_nobdry,3) = _view.ruling(i_nobdry) - (_view.face_normal(face_1).template head<3>()).cross(_view.face_normal(face_2).template head<3>()); 
        }

    }

    void get_constraint_ruling_1(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        _view.set_vector(vars);

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

            RealType w0 = _view.vertex_weight(node_0);
            RealType w1 = _view.vertex_weight(node_1);
            RealType w2 = _view.vertex_weight(node_2);
            RealType w3 = _view.vertex_weight(node_3);

            VectorType ruling_01 = _view.ruling(_it_he.to_nbdry(heh_01));
            VectorType ruling_12 = _view.ruling(_it_he.to_nbdry(heh_12));
            VectorType ruling_23 = _view.ruling(_it_he.to_nbdry(heh_23));
            VectorType ruling_30 = _view.ruling(_it_he.to_nbdry(heh_30));

            Dest.segment(_cons_idx["ruling_1"] + 3*i_nobdry,3) = _view.reweighted_ruling_1(i_nobdry) - (w1 + w2)*ruling_12 + (w3 + w0)*(ruling_30);
        }
    }

    void get_constraint_ruling_2(const VectorType &vars, VectorType &Dest)const {

        size_t num_constraints = _cons_idx["num_cons"];

        _view.set_vector(vars);

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

            RealType w0 = _view.vertex_weight(node_0);
            RealType w1 = _view.vertex_weight(node_1);
            RealType w2 = _view.vertex_weight(node_2);
            RealType w3 = _view.vertex_weight(node_3);

            VectorType ruling_01 = _view.ruling(_it_he.to_nbdry(heh_01));
            VectorType ruling_12 = _view.ruling(_it_he.to_nbdry(heh_12));
            VectorType ruling_23 = _view.ruling(_it_he.to_nbdry(heh_23));
            VectorType ruling_30 = _view.ruling(_it_he.to_nbdry(heh_30));

            Dest.segment(_cons_idx["ruling_2"] + 3*i_nobdry,3) = _view.reweighted_ruling_2(i_nobdry) - (w0 + w1)*ruling_01 + (w2 + w3)*(ruling_23);
        }
    }

    void get_constraint_dev(const VectorType &vars, VectorType &Dest) const{

        size_t num_constraints = _cons_idx["num_cons"];

        _view.set_vector(vars);

        std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        _it_face.reset();

        for(_it_face; _it_face.valid(); _it_face++){
            int i = _it_face.idx();
            int i_nobdry = _it_face.idx_nobdry();
            Dest.segment(_cons_idx["dev"] + 3*i_nobdry,3) = (_view.reweighted_ruling_1(i_nobdry). template head<3>()).cross(_view.reweighted_ruling_2(i_nobdry). template head<3>());
        }
    }

    void get_constraint_fair_v(const VectorType &vars, VectorType &Dest) const
    {

        size_t num_constraints = _cons_idx["num_cons"];

        _view.set_vector(vars);
        
        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(int i = 0; i < _stripHandle.num_vertexStrips(); i++){
            auto triple = _stripHandle.vertex_strip(i);
            int node_0 = std::get<0>(triple);
            int node_1 = std::get<1>(triple);
            int node_2 = std::get<2>(triple);
            Dest.segment(_cons_idx["fair_v"] + 3*i,3) = _view.vertex(node_0) - 2*_view.vertex(node_1) + _view.vertex(node_2);
        }
    }

    void get_constraint_fair_n(const VectorType &vars, VectorType &Dest) const
    {
        size_t num_constraints = _cons_idx["num_cons"];

        _view.set_vector(vars);

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

            VectorType normal_0 = _view.face_normal(face_0);
            VectorType normal_1 = _view.face_normal(face_1);
            VectorType normal_2 = _view.face_normal(face_2);

            Dest.segment(_cons_idx["fair_n"] + 3*i,3) = normal_0 - 2*normal_1 + normal_2;
        }

    }

    void get_constraint_fair_r(const VectorType &vars, VectorType &Dest) const
    {
        size_t num_constraints = _cons_idx["num_cons"];

        _view.set_vector(vars);

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

            VectorType ruling_0 = _view.ruling(_it_he.to_nbdry(he_0));
            VectorType ruling_1 = _view.ruling(_it_he.to_nbdry(he_1));
            VectorType ruling_2 = _view.ruling(_it_he.to_nbdry(he_2));
            if(_it_he.to_nbdry(he_0) == -1)
            {
                std::cout<<"Halfedge index he_0: "<<he_0<<std::endl;
            }
            if(_it_he.to_nbdry(he_1) == -1)
            {
                std::cout<<"Halfedge index he_1: "<<he_1<<std::endl;
            }
            if(_it_he.to_nbdry(he_2) == -1)
            {
                std::cout<<"Halfedge index he_2: "<<he_2<<std::endl;
            }

            Dest.segment(_cons_idx["fair_r"] + 3*i,3) = ruling_0 - 2*ruling_1 + ruling_2;
        }
    }

    void apply(const VectorType &vars, VectorType &Dest) const override
    {
        Dest.setZero();
        get_constraint_vertex_1(vars, Dest);
        get_constraint_vertex_2(vars, Dest);
        get_constraint_bdry_opt(vars, Dest);
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
     // used to store the boundary indices and corresponding values
    typedef typename std::optional<std::pair<std::vector<int>, VectorType>> OptionalBoundaryData;

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

        OptionalBoundaryData _bdryData;

        mutable std::map<std::string,size_t> _cons_idx;

        // Serves to access the current vector of input variables at the right positions
        mutable VectorView<ConfiguratorType> _view;
        // Serves to write the constraints into the vector
        // as well as the constraint gradients into the matrix
        mutable ConstraintView<ConfiguratorType> _cons_view;
        // no copy is made here
        //using _cons_idx = decltype(_cons_view._cons_idx)&;

    public:
        ConstraintGrad(QuadMeshTopologySaver &quadTopol, 
                   StripHandler<ConfiguratorType> &stripHandle,
                   OptionalBoundaryData bdryData = std::nullopt):
          _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges()),
          _num_bdry_faces(quadTopol.getNumBdryFaces()),
          _num_bdry_halfedges(quadTopol.getNumBdryHalfEdges()),
          _view(quadTopol),
          _stripHandle(stripHandle),
          _cons_view(stripHandle, quadTopol, bdryData),
          _it_face(quadTopol),
          _it_he(quadTopol),
          _bdryData(bdryData)
        {
            _cons_idx = _cons_view._cons_idx;

            if(_bdryData.has_value())
            {
                if(3*_bdryData.value().first.size() != _bdryData.value().second.size())
                {
                    std::cerr << "Boundary data does not match in size"<<std::endl;
                }
            }
        }

        void get_constraintgrad_vertex_1(const VectorType &vars, MatrixType &Dest) const {

            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows()!= num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }
            
            for(int i = 0; i < _num_vertices; i++){
                RealType dummy_i = _view.dummy_weight(i);
                // Assign values directly
                Dest.coeffRef(_cons_idx["vertex_1"] + i, i + _view._idx["dummy_weights"]) = -2.0 * dummy_i;
                Dest.coeffRef(_cons_idx["vertex_1"] + i, i + _view._idx["weights"]) = 1.0;
            }
        }

        void get_constraintgrad_vertex_1_triplet(TripletVectorType &triplets) const {
            
            for(int i = 0; i < _num_vertices; i++){
                RealType dummy_i = _view.dummy_weight(i);
                // Assign values directly
                triplets.emplace_back(_cons_idx["vertex_1"] + i, i + _view._idx["dummy_weights"], -2.0 * dummy_i);
                triplets.emplace_back(_cons_idx["vertex_1"] + i, i + _view._idx["weights"], 1.0);
            }
        }

        void get_constraintgrad_vertex_2(const VectorType &vars, MatrixType &Dest) const {
                    
            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            auto Id = MatrixType(3,3);
            Id.setIdentity();

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _num_vertices; i++){
                RealType w_i = _view.vertex_weight(i);
                VectorType vec = -_view.vertex(i);
                MatrixType vec_mat = vectorToSparseMat(vec);
                assignSparseBlockInplace(Dest,-w_i*Id,_cons_idx["vertex_2"] + 3*i,3*i + _view._idx["vertices"],triplet);
                assignSparseBlockInplace(Dest,vec_mat,_cons_idx["vertex_2"] + 3*i,i + _view._idx["weights"],triplet);
                assignSparseBlockInplace(Dest,Id,_cons_idx["vertex_2"] +3*i,3*i + _view._idx["reweighted_vertices"],triplet);
            }

        }

        void get_constraintgrad_vertex_2_triplet(TripletVectorType &triplet) const {

            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_vertices; i++){
                RealType w_i = _view.vertex_weight(i);
                VectorType vec = -_view.vertex(i);
                assignSparseMatBlockTriplet(-w_i*Id,_cons_idx["vertex_2"] + 3*i,3*i + _view._idx["vertices"],triplet);
                assignSparseVecBlockTriplet(vec,_cons_idx["vertex_2"] + 3*i,i + _view._idx["weights"],triplet);
                assignSparseMatBlockTriplet(Id,_cons_idx["vertex_2"] +3*i,3*i + _view._idx["reweighted_vertices"],triplet);
            }

        }

        void get_constraintgrad_bdry_opt(const VectorType &vars, MatrixType &Dest) const{
            
            size_t num_dofs = _view._idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            if(!_bdryData.has_value())
            {
                // boundary mask is not set, so just return
                return; 
            }

            _view.set_vector(vars);

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _bdryData.value().first.size(); i++){
                // Set the derivatives w.r.t. the normals
                auto vertex_now = _view.vertex(_bdryData.value().first[i]);
                auto vertex_pos = _bdryData.value().second.segment(3*i,3);
                VectorType vec = vertex_now - vertex_pos;
                MatrixType vec_1 = vectorToSparseMat(vec,true);
                assignSparseBlockInplace(Dest,2.0*vec_1, _cons_idx["bdry_opt"] + i,_view._idx["vertex"] + 3*_bdryData.value().first[i],triplet);
            }
        }

        void get_constraintgrad_bdry_opt_triplet(TripletVectorType &triplet) const{
            
            if(!_bdryData.has_value())
            {
                // boundary mask is not set, so just return
                return; 
            }

            for(int i = 0; i < _bdryData.value().first.size(); i++){
                // Set the derivatives w.r.t. the normals
                auto vertex_now = _view.vertex(_bdryData.value().first[i]);
                auto vertex_pos = _bdryData.value().second.segment(3*i,3);
                VectorType vec = vertex_now - vertex_pos;
                assignSparseVecBlockTriplet(2.0*vec, _cons_idx["bdry_opt"] + i,_view._idx["vertex"] + 3*_bdryData.value().first[i],triplet);
            }
        }

        void get_constraintgrad_edge_vec_1(const VectorType &vars, MatrixType &Dest)const {

            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
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

                RealType w0 = _view.vertex_weight(node_0);
                RealType w1 = _view.vertex_weight(node_1);
                RealType w2 = _view.vertex_weight(node_2);
                RealType w3 = _view.vertex_weight(node_3);

                VectorType v0_tilde = _view.reweighted_vertex(node_0);
                VectorType v1_tilde = _view.reweighted_vertex(node_1);
                VectorType v2_tilde = _view.reweighted_vertex(node_2);
                VectorType v3_tilde = _view.reweighted_vertex(node_3);

                assignSparseBlockInplace(Dest,Id,_cons_idx["edge_vec_1"] + 3*i,3*i+_view._idx["reweighted_edges_1"],triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                MatrixType sum_1 = vectorToSparseMat((v2_tilde + v3_tilde));
                MatrixType sum_2 = vectorToSparseMat((v1_tilde + v0_tilde));
                assignSparseBlockInplace(Dest, -sum_1,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -sum_1,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, sum_2,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, sum_2,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseBlockInplace(Dest, (w2+w3)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseBlockInplace(Dest, (w2+w3)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseBlockInplace(Dest, -(w0+w1)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseBlockInplace(Dest, -(w0+w1)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_3,triplet);
            }

        }

        void get_constraintgrad_edge_vec_1_triplet(TripletVectorType &triplet)const {

            // From average edge vector 1
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                RealType w0 = _view.vertex_weight(node_0);
                RealType w1 = _view.vertex_weight(node_1);
                RealType w2 = _view.vertex_weight(node_2);
                RealType w3 = _view.vertex_weight(node_3);

                VectorType v0_tilde = _view.reweighted_vertex(node_0);
                VectorType v1_tilde = _view.reweighted_vertex(node_1);
                VectorType v2_tilde = _view.reweighted_vertex(node_2);
                VectorType v3_tilde = _view.reweighted_vertex(node_3);

                assignSparseMatBlockTriplet(Id,_cons_idx["edge_vec_1"] + 3*i,3*i+_view._idx["reweighted_edges_1"],triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                VectorType sum_1 = (v2_tilde + v3_tilde);
                VectorType sum_2 = (v1_tilde + v0_tilde);
                assignSparseVecBlockTriplet(-sum_1,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_0,triplet);
                assignSparseVecBlockTriplet(-sum_1,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_1,triplet);
                assignSparseVecBlockTriplet(sum_2,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_2,triplet);
                assignSparseVecBlockTriplet(sum_2,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseMatBlockTriplet((w2+w3)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseMatBlockTriplet((w2+w3)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseMatBlockTriplet(-(w0+w1)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseMatBlockTriplet(-(w0+w1)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_3,triplet);
            }

        }

        void get_constraintgrad_edge_vec_2(const VectorType &vars, MatrixType &Dest)const {

            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
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

                RealType w0 = _view.vertex_weight(node_0);
                RealType w1 = _view.vertex_weight(node_1);
                RealType w2 = _view.vertex_weight(node_2);
                RealType w3 = _view.vertex_weight(node_3);

                VectorType v0_tilde = _view.reweighted_vertex(node_0);
                VectorType v1_tilde = _view.reweighted_vertex(node_1);
                VectorType v2_tilde = _view.reweighted_vertex(node_2);
                VectorType v3_tilde = _view.reweighted_vertex(node_3);

                assignSparseBlockInplace(Dest,Id,_cons_idx["edge_vec_2"] + 3*i,3*i+_view._idx["reweighted_edges_2"],triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                MatrixType sum_1 = vectorToSparseMat((v2_tilde + v1_tilde));
                MatrixType sum_2 = vectorToSparseMat((v0_tilde + v3_tilde));
    
                assignSparseBlockInplace(Dest, sum_1,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -sum_2,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, -sum_2,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, sum_1,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseBlockInplace(Dest,-(w1+w2)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseBlockInplace(Dest,(w0+w3)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseBlockInplace(Dest,(w0+w3)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseBlockInplace(Dest,-(w1+w2)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_3,triplet);
            }

        }

        void get_constraintgrad_edge_vec_2_triplet(TripletVectorType &triplet)const {

            // From average edge vector 1
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                RealType w0 = _view.vertex_weight(node_0);
                RealType w1 = _view.vertex_weight(node_1);
                RealType w2 = _view.vertex_weight(node_2);
                RealType w3 = _view.vertex_weight(node_3);

                VectorType v0_tilde = _view.reweighted_vertex(node_0);
                VectorType v1_tilde = _view.reweighted_vertex(node_1);
                VectorType v2_tilde = _view.reweighted_vertex(node_2);
                VectorType v3_tilde = _view.reweighted_vertex(node_3);

                assignSparseMatBlockTriplet(Id,_cons_idx["edge_vec_2"] + 3*i,3*i+_view._idx["reweighted_edges_2"],triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                VectorType sum_1 = v2_tilde + v1_tilde;
                VectorType sum_2 = v0_tilde + v3_tilde;
    
                assignSparseVecBlockTriplet(sum_1,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_0,triplet);
                assignSparseVecBlockTriplet(-sum_2,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_1,triplet);
                assignSparseVecBlockTriplet(-sum_2,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_2,triplet);
                assignSparseVecBlockTriplet(sum_1,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseMatBlockTriplet(-(w1+w2)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseMatBlockTriplet((w0+w3)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseMatBlockTriplet((w0+w3)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseMatBlockTriplet(-(w1+w2)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_3,triplet);
            }

        }

        void get_constraintgrad_normal_1(const VectorType &vars, MatrixType &Dest)const {

            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs ){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                MatrixType vec_1 = vectorToSparseMat(_view.reweighted_edge_1(i),true);
                MatrixType vec_2 = vectorToSparseMat(_view.face_normal(i),true);
                assignSparseBlockInplace(Dest,vec_1,_cons_idx["normal_1"] + i,_view._idx["normals"] + 3*i,triplet);
                assignSparseBlockInplace(Dest,vec_2,_cons_idx["normal_1"] + i,_view._idx["reweighted_edges_1"] + 3*i,triplet);
            }
        }

        void get_constraintgrad_normal_1_triplet(TripletVectorType &triplet)const {

            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                assignSparseVecBlockTriplet(_view.reweighted_edge_1(i),_cons_idx["normal_1"] + i,_view._idx["normals"] + 3*i,triplet,true);
                assignSparseVecBlockTriplet(_view.face_normal(i),_cons_idx["normal_1"] + i,_view._idx["reweighted_edges_1"] + 3*i,triplet,true);
            }
        }

        void get_constraintgrad_normal_2(const VectorType &vars, MatrixType &Dest)const {

            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs ){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                MatrixType vec_1 = vectorToSparseMat(_view.reweighted_edge_2(i),true);
                MatrixType vec_2 = vectorToSparseMat(_view.face_normal(i),true);
                assignSparseBlockInplace(Dest,vec_1,_cons_idx["normal_2"] + i,_view._idx["normals"] + 3*i,triplet);
                assignSparseBlockInplace(Dest,vec_2,_cons_idx["normal_2"] + i,_view._idx["reweighted_edges_2"] + 3*i,triplet);
            }
        }

        void get_constraintgrad_normal_2_triplet(TripletVectorType &triplet)const {

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                assignSparseVecBlockTriplet(_view.reweighted_edge_2(i),_cons_idx["normal_2"] + i,_view._idx["normals"] + 3*i,triplet,true);
                assignSparseVecBlockTriplet(_view.face_normal(i),_cons_idx["normal_2"] + i,_view._idx["reweighted_edges_2"] + 3*i,triplet,true);
            }
        }

        void get_constraintgrad_normal_3(const VectorType &vars, MatrixType &Dest)const {

            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs ){
                Dest.resize(num_constraints, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                MatrixType vec_1 = vectorToSparseMat(_view.face_normal(i),true);
                assignSparseBlockInplace(Dest,2.0*vec_1,_cons_idx["normal_3"] + i,_view._idx["normals"] + 3*i,triplet);
            }
        }

        void get_constraintgrad_normal_3_triplet(TripletVectorType &triplet)const {

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                assignSparseVecBlockTriplet(2.0*_view.face_normal(i),_cons_idx["normal_3"] + i,_view._idx["normals"] + 3*i,triplet,true);
            }
        }
        
        void get_constraintgrad_ruling_0(const VectorType &vars, MatrixType &Dest)const {

            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
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
                Eigen::Vector3d normal_1 = _view.face_normal(faceIdx_1).template head<3>();
                Eigen::Vector3d normal_2 = _view.face_normal(faceIdx_2).template head<3>();

                // First, set the derivative w.r.t. the ruling vector
                auto Id = MatrixType(3,3);
                Id.setIdentity();

                assignSparseBlockInplace(Dest,Id,_cons_idx["ruling_0"] + 3*(i_nobdry),_view._idx["rulings"] + 3*i_nobdry,triplet);

                // Now, set the derivatives w.r.t. the cross product
                // First, differentiate w.r.t. normal_1
                Eigen::Vector3d e;
                MatrixType vec_prod;
                e.setZero();
                e[0] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_1,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_1 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_1 + 2,triplet);

                // Now, differentiate w.r.t. normal_2
                e[2] = 0.0;
                e[0] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_2,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_2 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_2 + 2,triplet);
            }
        }

        void get_constraintgrad_ruling_0_triplet(TripletVectorType &triplet)const {

            _it_he.reset();

            for(_it_he; _it_he.valid();_it_he++)
            {
                int i = _it_he.idx();
                int i_nobdry = _it_he.idx_nobdry();

                size_t faceIdx_1 = _quadTopol.getFaceOfHalfEdge(i,0);
                size_t faceIdx_2 = _quadTopol.getFaceOfHalfEdge(i,1);
                Eigen::Vector3d normal_1 = _view.face_normal(faceIdx_1).template head<3>();
                Eigen::Vector3d normal_2 = _view.face_normal(faceIdx_2).template head<3>();

                // First, set the derivative w.r.t. the ruling vector
                auto Id = MatrixType(3,3);
                Id.setIdentity();

                assignSparseMatBlockTriplet(Id,_cons_idx["ruling_0"] + 3*(i_nobdry),_view._idx["rulings"] + 3*i_nobdry,triplet);

                // Now, set the derivatives w.r.t. the cross product
                // First, differentiate w.r.t. normal_1
                Eigen::Vector3d e;
                VectorType vec_prod;
                e.setZero();
                e[0] = 1.0;
                vec_prod = -e.cross(normal_2);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_1,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = -e.cross(normal_2);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_1 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = -e.cross(normal_2);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_1 + 2, triplet);

                // Now, differentiate w.r.t. normal_2
                e[2] = 0.0;
                e[0] = 1.0;
                vec_prod = e.cross(normal_1);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_2,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = e.cross(normal_1);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_2 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = e.cross(normal_1);
                assignSparseVecBlockTriplet(vec_prod,_cons_idx["ruling_0"]  + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_2 + 2,triplet);
            }
        }

        void get_constraintgrad_ruling_1(const VectorType& vars, MatrixType &Dest)const {

            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];

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

                auto node_0_vert = _view.vertex(node_0);
                auto node_1_vert = _view.vertex(node_1);
                auto node_2_vert = _view.vertex(node_2);
                auto node_3_vert = _view.vertex(node_3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                RealType w0 = _view.vertex_weight(node_0);
                RealType w1 = _view.vertex_weight(node_1);
                RealType w2 = _view.vertex_weight(node_2);
                RealType w3 = _view.vertex_weight(node_3);

                VectorType ruling_01 = _view.ruling(_it_he.to_nbdry(heh_01));
                VectorType ruling_12 = _view.ruling(_it_he.to_nbdry(heh_12));
                VectorType ruling_23 = _view.ruling(_it_he.to_nbdry(heh_23));
                VectorType ruling_30 = _view.ruling(_it_he.to_nbdry(heh_30));

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}

                assignSparseBlockInplace(Dest,Id,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["reweighted_rulings_1"] + 3*i_nobdry,triplet);

                MatrixType ruling_30_mat = vectorToSparseMat(ruling_30);
                MatrixType ruling_12_mat = vectorToSparseMat(ruling_12);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseBlockInplace(Dest, ruling_30_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -ruling_12_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, -ruling_12_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, ruling_30_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseBlockInplace(Dest,-(w1 + w2)*Id,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_12),triplet);
                assignSparseBlockInplace(Dest,(w3 + w0)*Id,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_30),triplet);
            }
        }

        void get_constraintgrad_ruling_1_triplet(TripletVectorType &triplet)const {

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

                auto node_0_vert = _view.vertex(node_0);
                auto node_1_vert = _view.vertex(node_1);
                auto node_2_vert = _view.vertex(node_2);
                auto node_3_vert = _view.vertex(node_3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                RealType w0 = _view.vertex_weight(node_0);
                RealType w1 = _view.vertex_weight(node_1);
                RealType w2 = _view.vertex_weight(node_2);
                RealType w3 = _view.vertex_weight(node_3);

                VectorType ruling_01 = _view.ruling(_it_he.to_nbdry(heh_01));
                VectorType ruling_12 = _view.ruling(_it_he.to_nbdry(heh_12));
                VectorType ruling_23 = _view.ruling(_it_he.to_nbdry(heh_23));
                VectorType ruling_30 = _view.ruling(_it_he.to_nbdry(heh_30));

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}

                assignSparseMatBlockTriplet(Id,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["reweighted_rulings_1"] + 3*i_nobdry,triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseVecBlockTriplet(ruling_30,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_0,triplet);
                assignSparseVecBlockTriplet(-ruling_12,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_1,triplet);
                assignSparseVecBlockTriplet(-ruling_12,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_2,triplet);
                assignSparseVecBlockTriplet(ruling_30,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseMatBlockTriplet(-(w1 + w2)*Id,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_12),triplet);
                assignSparseMatBlockTriplet((w3 + w0)*Id,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_30),triplet);
            }
        }
        
        void get_constraintgrad_ruling_2(const VectorType& vars, MatrixType &Dest)const {

            _view.set_vector(vars);

            std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

            size_t num_dofs = _view._idx["num_dofs"];

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

                RealType w0 = _view.vertex_weight(node_0);
                RealType w1 = _view.vertex_weight(node_1);
                RealType w2 = _view.vertex_weight(node_2);
                RealType w3 = _view.vertex_weight(node_3);

                VectorType ruling_01 = _view.ruling(_it_he.to_nbdry(heh_01));
                VectorType ruling_23 = _view.ruling(_it_he.to_nbdry(heh_23));

                MatrixType ruling_01_mat = vectorToSparseMat(ruling_01);
                MatrixType ruling_23_mat = vectorToSparseMat(ruling_23);

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}
                assignSparseBlockInplace(Dest,Id,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["reweighted_rulings_2"] + 3*i_nbdry,triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseBlockInplace(Dest, -ruling_01_mat,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -ruling_01_mat,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, ruling_23_mat,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, ruling_23_mat,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseBlockInplace(Dest,-(w0 + w1)*Id,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_01),triplet);
                assignSparseBlockInplace(Dest,(w2 + w3)*Id,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_23),triplet);
            }
        }

        void get_constraintgrad_ruling_2_triplet(TripletVectorType &triplet)const {

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

                RealType w0 = _view.vertex_weight(node_0);
                RealType w1 = _view.vertex_weight(node_1);
                RealType w2 = _view.vertex_weight(node_2);
                RealType w3 = _view.vertex_weight(node_3);

                VectorType ruling_01 = _view.ruling(_it_he.to_nbdry(heh_01));
                VectorType ruling_23 = _view.ruling(_it_he.to_nbdry(heh_23));

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}
                assignSparseMatBlockTriplet(Id,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["reweighted_rulings_2"] + 3*i_nbdry,triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseVecBlockTriplet(-ruling_01,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["weights"] + node_0,triplet);
                assignSparseVecBlockTriplet(-ruling_01,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["weights"] + node_1,triplet);
                assignSparseVecBlockTriplet(ruling_23,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["weights"] + node_2,triplet);
                assignSparseVecBlockTriplet(ruling_23,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseMatBlockTriplet(-(w0 + w1)*Id,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_01),triplet);
                assignSparseMatBlockTriplet((w2 + w3)*Id,_cons_idx["ruling_2"] + 3*i_nbdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_23),triplet);
            }
        }
        
        void get_constraintgrad_dev(const VectorType &vars, MatrixType &Dest)const {

            _view.set_vector(vars);

            std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

            size_t num_constraints = _cons_idx["num_cons"];
            size_t num_dofs = _view._idx["num_dofs"];

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
                Eigen::Vector3d reweighted_ruling_1 = _view.reweighted_ruling_1(i_nobdry).template head<3>();
                Eigen::Vector3d reweighted_ruling_2 = _view.reweighted_ruling_2(i_nobdry).template head<3>();

                // Derivatives w.r.t. the first ruling vector \Tilde{r}_{13,f}
                Eigen::Vector3d e;
                e.setZero();
                e[0] = 1.0; 
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_1"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_1"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_1"] + 2,triplet);
            
                // Derivatives w.r.t. the second ruling vector \Tilde{r}_{02,f}
                e[2] = 0.0;
                e[0] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_2"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_2"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_2"] + 2,triplet);
            }
        }

        void get_constraintgrad_dev_triplet(TripletVectorType &triplet)const {

            _it_face.reset();

            for(_it_face; _it_face.valid(); _it_face++)
            {

                int i = _it_face.idx();
                int i_nobdry = _it_face.idx_nobdry();

                // obtain the reweighted ruling vectors
                Eigen::Vector3d reweighted_ruling_1 = _view.reweighted_ruling_1(i_nobdry).template head<3>();
                Eigen::Vector3d reweighted_ruling_2 = _view.reweighted_ruling_2(i_nobdry).template head<3>();

                // Derivatives w.r.t. the first ruling vector \Tilde{r}_{13,f}
                Eigen::Vector3d e;
                e.setZero();
                e[0] = 1.0; 
                assignSparseVecBlockTriplet(e.cross(reweighted_ruling_2),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_1"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                assignSparseVecBlockTriplet(e.cross(reweighted_ruling_2),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_1"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseVecBlockTriplet(e.cross(reweighted_ruling_2),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_1"] + 2,triplet);
            
                // Derivatives w.r.t. the second ruling vector \Tilde{r}_{02,f}
                e[2] = 0.0;
                e[0] = 1.0;
                assignSparseVecBlockTriplet(reweighted_ruling_1.cross(e),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_2"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                assignSparseVecBlockTriplet(reweighted_ruling_1.cross(e),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_2"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseVecBlockTriplet(reweighted_ruling_1.cross(e),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_2"] + 2,triplet);
            }
        }

        void get_constraintgrad_fair_v(const VectorType &vars, MatrixType &Dest)const {
            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
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
                
                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_v"] + 3*i,3*node_0 + _view._idx["vertices"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,_cons_idx["fair_v"] + 3*i,3*node_1 + _view._idx["vertices"], triplet);
                assignSparseBlockInplace(Dest, Id, _cons_idx["fair_v"] + 3*i,3*node_2 + _view._idx["vertices"], triplet);
            }
        }

        void get_constraintgrad_fair_v_triplet(TripletVectorType &triplet)const {

            MatrixType Id(3,3);
            Id.setIdentity();

            for(int i = 0; i < _stripHandle.num_vertexStrips(); i++){
                auto triple = _stripHandle.vertex_strip(i);
                int node_0 = std::get<0>(triple);
                int node_1 = std::get<1>(triple);
                int node_2 = std::get<2>(triple);
                
                assignSparseMatBlockTriplet(Id,_cons_idx["fair_v"] + 3*i,3*node_0 + _view._idx["vertices"],triplet);
                assignSparseMatBlockTriplet(-2*Id,_cons_idx["fair_v"] + 3*i,3*node_1 + _view._idx["vertices"], triplet);
                assignSparseMatBlockTriplet(Id, _cons_idx["fair_v"] + 3*i,3*node_2 + _view._idx["vertices"], triplet);
            }
        }

        void get_constraintgrad_fair_n(const VectorType &vars, MatrixType &Dest)const {
            _view.set_vector(vars);

            size_t num_constraints = _cons_idx["num_cons"];
            size_t num_dofs = _view._idx["num_dofs"];

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

                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_n"] + 3*i,3*normal_0 + _view._idx["normals"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,_cons_idx["fair_n"] + 3*i,3*normal_1 + _view._idx["normals"],triplet);
                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_n"] + 3*i,3*normal_2 + _view._idx["normals"],triplet);
            }
        }

        void get_constraintgrad_fair_n_triplet(TripletVectorType &triplet)const {

            MatrixType Id(3,3);
            Id.setIdentity();

            for(int i = 0; i < _stripHandle.num_faceStrips(); i++){
                auto triple = _stripHandle.face_strip(i);
                int normal_0 = std::get<0>(triple);
                int normal_1 = std::get<1>(triple);
                int normal_2 = std::get<2>(triple);

                assignSparseMatBlockTriplet(Id,_cons_idx["fair_n"] + 3*i,3*normal_0 + _view._idx["normals"],triplet);
                assignSparseMatBlockTriplet(-2*Id,_cons_idx["fair_n"] + 3*i,3*normal_1 + _view._idx["normals"],triplet);
                assignSparseMatBlockTriplet(Id,_cons_idx["fair_n"] + 3*i,3*normal_2 + _view._idx["normals"],triplet);
            }
        }

        void get_constraintgrad_fair_r(const VectorType &vars, MatrixType &Dest)const {
            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
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

                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_0) + _view._idx["rulings"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_1) + _view._idx["rulings"],triplet);
                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_2) + _view._idx["rulings"],triplet);
            }
        }

        void get_constraintgrad_fair_r_triplet(TripletVectorType &triplet)const {

            MatrixType Id(3,3);
            Id.setIdentity();

            for(int i = 0; i < _stripHandle.num_halfedgeStrips(); i++)
            {
                auto triple = _stripHandle.halfedge_strip(i);
                int heh_0 = std::get<0>(triple);
                int heh_1 = std::get<1>(triple);
                int heh_2 = std::get<2>(triple);

                assignSparseMatBlockTriplet(Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_0) + _view._idx["rulings"],triplet);
                assignSparseMatBlockTriplet(-2*Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_1) + _view._idx["rulings"],triplet);
                assignSparseMatBlockTriplet(Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_2) + _view._idx["rulings"],triplet);
            }
        }


        void apply(const VectorType &vars, MatrixType &Dest) const override
        {
            size_t num_dofs = _view._idx["num_dofs"];
            // Change num_constraints back later !
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            Dest.setZero();
            _view.set_vector(vars);
            
            TripletVectorType triplet;
            size_t num_vertex_1_triplets = 2*_num_vertices;
            size_t num_vertex_2_triplets = 21*_num_vertices;
            size_t num_bdry_opt_triplets = _bdryData.has_value() ? 3*_bdryData.value().first.size() : 0;
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
                            num_bdry_opt_triplets +
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
            get_constraintgrad_vertex_1_triplet(triplet);
            get_constraintgrad_vertex_2_triplet(triplet);
            get_constraintgrad_bdry_opt_triplet(triplet);
            get_constraintgrad_edge_vec_1_triplet(triplet);
            get_constraintgrad_edge_vec_2_triplet(triplet);
            get_constraintgrad_normal_1_triplet(triplet);
            get_constraintgrad_normal_2_triplet(triplet);
            get_constraintgrad_normal_3_triplet(triplet);
            get_constraintgrad_ruling_0_triplet(triplet);
            get_constraintgrad_ruling_1_triplet(triplet);
            get_constraintgrad_ruling_2_triplet(triplet);
            get_constraintgrad_dev_triplet(triplet);
            get_constraintgrad_fair_v_triplet(triplet);
            get_constraintgrad_fair_n_triplet(triplet);
            get_constraintgrad_fair_r_triplet(triplet);

            Dest.setFromTriplets(triplet.begin(), triplet.end());
        }

        /*
        What has been achieved here ?
        1. vertex_1 done
        2. vertex_2 done
        3. edge_vec_1 done
        4. edge_vec_2 done
        5. normal_1 done
        6. normal_2 done
        7. normal_3 done
        8. ruling_0 done
        9. ruling_1 done
        10. ruling_2 done
        11. dev done
        12. fair_v done
        13. fair_n
        14. fair_r
        void get_constraintgrads(const VectorType &vars, MatrixType &Dest){
            _view.set_vector(vars);

            size_t num_dofs = _view._idx["num_dofs"];
            size_t num_constraints = _cons_idx["num_cons"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            Dest.setZero();

            std::vector<Eigen::Triplet<RealType>> triplet;
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            // All vertices constraints (vertex_1, vertex_2)
            for(int i = 0; i < _num_vertices; i++){
                RealType dummy_i = _view.dummy_weight(i);
                RealType w_i = _view.vertex_weight(i);
                VectorType vec = -_view.vertex(i);
                MatrixType vec_mat = vectorToSparseMat(vec);

                // Assign values directly
                Dest.coeffRef(i + _cons_idx["vertex_1"], i + _view._idx["dummy_weights"]) = -2.0 * dummy_i;
                Dest.coeffRef(i + _cons_idx["vertex_1"], i + _view._idx["weights"]) = 1.0;
                assignSparseBlockInplace(Dest,-w_i*Id,_cons_idx["vertex_2"] + 3*i,3*i + _view._idx["vertices"],triplet);
                assignSparseBlockInplace(Dest,vec_mat,_cons_idx["vertex_2"] + 3*i,i + _view._idx["weights"],triplet);
                assignSparseBlockInplace(Dest,Id,_cons_idx["vertex_2"] +3*i,3*i + _view._idx["reweighted_vertices"],triplet);
            }

            // All face constraints (edge_vec_1, edge_vec_2, normal_1, normal_2, normal_3)
            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                RealType w0 = _view.vertex_weight(node_0);
                RealType w1 = _view.vertex_weight(node_1);
                RealType w2 = _view.vertex_weight(node_2);
                RealType w3 = _view.vertex_weight(node_3);

                VectorType v0_tilde = _view.reweighted_vertex(node_0);
                VectorType v1_tilde = _view.reweighted_vertex(node_1);
                VectorType v2_tilde = _view.reweighted_vertex(node_2);
                VectorType v3_tilde = _view.reweighted_vertex(node_3);

                assignSparseBlockInplace(Dest,Id,_cons_idx["edge_vec_1"] + 3*i,3*i+_view._idx["reweighted_edges_1"],triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                MatrixType sum_1 = vectorToSparseMat((v2_tilde + v3_tilde));
                MatrixType sum_2 = vectorToSparseMat((v1_tilde + v0_tilde));
                assignSparseBlockInplace(Dest, -sum_1,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -sum_1, _cons_idx["edge_vec_1"]+ 3*i,_view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, sum_2,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, sum_2,_cons_idx["edge_vec_1"] + 3*i,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseBlockInplace(Dest, (w2+w3)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseBlockInplace(Dest, (w2+w3)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseBlockInplace(Dest, -(w0+w1)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseBlockInplace(Dest, -(w0+w1)*Id,_cons_idx["edge_vec_1"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_3,triplet);

                assignSparseBlockInplace(Dest,Id,_cons_idx["edge_vec_2"] + 3*i,3*i+_view._idx["reweighted_edges_2"],triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                MatrixType sum_3 = vectorToSparseMat((v2_tilde + v1_tilde));
                MatrixType sum_4 = vectorToSparseMat((v0_tilde + v3_tilde));
    
                assignSparseBlockInplace(Dest, sum_3,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -sum_4,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, -sum_4,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, sum_3,_cons_idx["edge_vec_2"] + 3*i,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseBlockInplace(Dest,-(w1+w2)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseBlockInplace(Dest,(w0+w3)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseBlockInplace(Dest,(w0+w3)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseBlockInplace(Dest,-(w1+w2)*Id,_cons_idx["edge_vec_2"] + 3*i,_view._idx["reweighted_vertices"] + 3*node_3,triplet);
            
                 // Set the derivatives w.r.t. the normals
                MatrixType rw_edge_1_mat = vectorToSparseMat(_view.reweighted_edge_1(i),true);
                MatrixType normal_mat = vectorToSparseMat(_view.face_normal(i),true);
                assignSparseBlockInplace(Dest,rw_edge_1_mat,_cons_idx["normal_1"] + i,_view._idx["normals"] + 3*i,triplet);
                assignSparseBlockInplace(Dest,normal_mat,_cons_idx["normal_1"] + i,_view._idx["reweighted_edges_1"] + 3*i,triplet);
            
                MatrixType rw_edge_2_mat = vectorToSparseMat(_view.reweighted_edge_2(i),true);
                assignSparseBlockInplace(Dest,rw_edge_2_mat,_cons_idx["normal_2"] + i,_view._idx["normals"] + 3*i,triplet);
                assignSparseBlockInplace(Dest,normal_mat,_cons_idx["normal_2"] + i,_view._idx["reweighted_edges_2"] + 3*i,triplet);
            
                assignSparseBlockInplace(Dest,2.0*normal_mat,_cons_idx["normal_3"] + i,_view._idx["normals"] + 3*i,triplet);
            }

            _it_he.reset();

            // ruling_0 constraint
            for(_it_he; _it_he.valid();_it_he++)
            {
                int i = _it_he.idx();
                int i_nobdry = _it_he.idx_nobdry();

                size_t faceIdx_1 = _quadTopol.getFaceOfHalfEdge(i,0);
                size_t faceIdx_2 = _quadTopol.getFaceOfHalfEdge(i,1);
                Eigen::Vector3d normal_1 = _view.face_normal(faceIdx_1).template head<3>();
                Eigen::Vector3d normal_2 = _view.face_normal(faceIdx_2).template head<3>();

                // First, set the derivative w.r.t. the ruling vector
                auto Id = MatrixType(3,3);
                Id.setIdentity();

                assignSparseBlockInplace(Dest,Id,_cons_idx["ruling_0"] + 3*(i_nobdry),_view._idx["rulings"] + 3*i_nobdry,triplet);

                // Now, set the derivatives w.r.t. the cross product
                // First, differentiate w.r.t. normal_1
                Eigen::Vector3d e;
                MatrixType vec_prod;
                e.setZero();
                e[0] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"] + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_1,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"] + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_1 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"] + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_1 + 2,triplet);

                // Now, differentiate w.r.t. normal_2
                e[2] = 0.0;
                e[0] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"] + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_2,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"] + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_2 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,_cons_idx["ruling_0"] + 3*(i_nobdry),_view._idx["normals"] + 3*faceIdx_2 + 2,triplet);
            }

            _it_face.reset();

            // constraints ruling_1, ruling_2
            for(_it_face; _it_face.valid(); _it_face++){

                int i = _it_face.idx();
                int i_nobdry = _it_face.idx_nobdry();

                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                auto node_0_vert = _view.vertex(node_0);
                auto node_1_vert = _view.vertex(node_1);
                auto node_2_vert = _view.vertex(node_2);
                auto node_3_vert = _view.vertex(node_3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                RealType w0 = _view.vertex_weight(node_0);
                RealType w1 = _view.vertex_weight(node_1);
                RealType w2 = _view.vertex_weight(node_2);
                RealType w3 = _view.vertex_weight(node_3);

                VectorType ruling_01 = _view.ruling(_it_he.to_nbdry(heh_01));
                VectorType ruling_12 = _view.ruling(_it_he.to_nbdry(heh_12));
                VectorType ruling_23 = _view.ruling(_it_he.to_nbdry(heh_23));
                VectorType ruling_30 = _view.ruling(_it_he.to_nbdry(heh_30));

                // obtain the reweighted ruling vectors
                Eigen::Vector3d reweighted_ruling_1 = _view.reweighted_ruling_1(i).template head<3>();
                Eigen::Vector3d reweighted_ruling_2 = _view.reweighted_ruling_2(i).template head<3>();

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}

                assignSparseBlockInplace(Dest,Id,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["reweighted_rulings_1"] + 3*i_nobdry,triplet);

                MatrixType ruling_30_mat = vectorToSparseMat(ruling_30);
                MatrixType ruling_12_mat = vectorToSparseMat(ruling_12);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseBlockInplace(Dest, ruling_30_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -ruling_12_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, -ruling_12_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, ruling_30_mat,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseBlockInplace(Dest,-(w1 + w2)*Id,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_12),triplet);
                assignSparseBlockInplace(Dest,(w3 + w0)*Id,_cons_idx["ruling_1"] + 3*i_nobdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_30),triplet);

                MatrixType ruling_01_mat = vectorToSparseMat(ruling_01);
                MatrixType ruling_23_mat = vectorToSparseMat(ruling_23);

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}
                assignSparseBlockInplace(Dest,Id,_cons_idx["ruling_2"] + 3*i_nobdry,_view._idx["reweighted_rulings_2"] + 3*i_nobdry,triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseBlockInplace(Dest, -ruling_01_mat,_cons_idx["ruling_2"] + 3*i_nobdry,_view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -ruling_01_mat,_cons_idx["ruling_2"] + 3*i_nobdry,_view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, ruling_23_mat,_cons_idx["ruling_2"] + 3*i_nobdry,_view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, ruling_23_mat,_cons_idx["ruling_2"] + 3*i_nobdry,_view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseBlockInplace(Dest,-(w0 + w1)*Id,_cons_idx["ruling_2"] + 3*i_nobdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_01),triplet);
                assignSparseBlockInplace(Dest,(w2 + w3)*Id,_cons_idx["ruling_2"]+ 3*i_nobdry,_view._idx["rulings"] + 3*_it_he.to_nbdry(heh_23),triplet);
            
                // Derivatives w.r.t. the first ruling vector \Tilde{r}_{13,f}
                Eigen::Vector3d e;
                e.setZero();
                e[0] = 1.0; 
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_1"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_1"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_1"] + 2,triplet); 
            
                // Derivatives w.r.t. the second ruling vector \Tilde{r}_{02,f}
                e[2] = 0.0;
                e[0] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_2"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_2"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),_cons_idx["dev"] + 3*i_nobdry,3*i_nobdry + _view._idx["reweighted_rulings_2"] + 2,triplet);   
            }

            for(int i = 0; i < _stripHandle.num_vertexStrips(); i++){
                auto triple = _stripHandle.vertex_strip(i);
                int node_0 = std::get<0>(triple);
                int node_1 = std::get<1>(triple);
                int node_2 = std::get<2>(triple);
                
                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_v"] + 3*i,3*node_0 + _view._idx["vertices"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,_cons_idx["fair_v"] + 3*i,3*node_1 + _view._idx["vertices"], triplet);
                assignSparseBlockInplace(Dest, Id, _cons_idx["fair_v"] + 3*i,3*node_2 + _view._idx["vertices"], triplet);
            }

            for(int i = 0; i < _stripHandle.num_faceStrips(); i++){
                auto triple = _stripHandle.face_strip(i);
                int normal_0 = std::get<0>(triple);
                int normal_1 = std::get<1>(triple);
                int normal_2 = std::get<2>(triple);

                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_n"] + 3*i,3*normal_0 + _view._idx["normals"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,_cons_idx["fair_n"] + 3*i,3*normal_1 + _view._idx["normals"],triplet);
                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_n"] + 3*i,3*normal_2 + _view._idx["normals"],triplet);
            }

            for(int i = 0; i < _stripHandle.num_halfedgeStrips(); i++)
            {
                auto triple = _stripHandle.halfedge_strip(i);
                int heh_0 = std::get<0>(triple);
                int heh_1 = std::get<1>(triple);
                int heh_2 = std::get<2>(triple);

                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_0) + _view._idx["rulings"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_1) + _view._idx["rulings"],triplet);
                assignSparseBlockInplace(Dest,Id,_cons_idx["fair_r"] + 3*i,3*_it_he.to_nbdry(heh_2) + _view._idx["rulings"],triplet);
            }
        }
*/
};