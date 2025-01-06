#include <goast/Quads.h>
#include <goast/QuadMesh/QuadTopology.h>
#include <cmath>
#include <goast/SQP/Utils/SparseMat.h>

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

        int id_nobdry(){
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
        VectorView(QuadMeshTopologySaver &quadTopol, std::optional<VectorType> vec = std::nullopt)
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
        
            if(vec.has_value()){
                if(vec.value().size() != _idx["num_dofs"]){
                    std::cerr << "Vector size does not match the number of dofs"<<std::endl;
                }
                _vec = vec.value();
            }
            else{
                _vec = VectorType::Zero(_idx["num_dofs"]);
            }

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

template <typename ConfiguratorType>
class Constraint : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    private:
        QuadMeshTopologySaver &_quadTopol;
        size_t _num_vertices;
        size_t _num_faces;
        size_t _num_edges;
        size_t _num_halfedges;

    public:
        Constraint(QuadMeshTopologySaver &quadTopol)
            : _quadTopol(quadTopol),
              _num_vertices(quadTopol.getNumVertices()),
              _num_faces(quadTopol.getNumFaces()),
              _num_edges(quadTopol.getNumEdges()),
              _num_halfedges(2*_num_edges)
        {
            VectorView<ConfiguratorType> view;
            view.set_vector(vars);

            std::vector<int> bdryHalfEdges = _quadTopol.getBdryHalfEdges();
            std::vector<int> bdryFaces = _quadTopol.getBdryFaces();
            size_t num_bdryHalfEdges = bdryHalfEdges.size();
            size_t num_bdryFaces = bdryFaces.size();

            std::vector<std::tuple<int,int,int>> vertex_strips = _quadTopol.getVertexStrips();
            size_t num_vertex_strips = vertex_strips.size();

            std::vector<std::tuple<int,int,int>> face_strips = _quadTopol.getFaceStrips();
            size_t num_face_strips = face_strips.size();

            // only strips of halfedges not located at the boundary
            std::vector<std::tuple<int,int,int>> halfedge_strips = _quadTopol.getHalfEdgeStrips();
            size_t num_halfedge_strips = halfedge_strips.size();

            size_t num_constraints_vert_1 = 3*_num_vertices;
            size_t num_constraints_vert_2 = _num_vertices;
            size_t num_constraints_edge_vec_1 = 3*_num_faces;
            size_t num_constraints_edge_vec_2 = 3*_num_faces;
            size_t num_constraints_normal_1 = _num_faces;
            size_t num_constraints_normal_2 = _num_faces;
            size_t num_constraints_normal_3 = _num_faces;
            size_t num_constraints_ruling_0 = 3*(_num_halfedges - num_bdryHalfEdges);
            size_t num_constraints_ruling_1 = 3*(_num_faces - num_bdryFaces);
            size_t num_constraints_ruling_2 = 3*(_num_faces - num_bdryFaces);
            size_t num_constraints_dev = 3*(_num_faces - num_bdryFaces);
            size_t num_constraints_fair_v = 3*num_vertex_strips;
            size_t num_constraints_fair_n = 3*num_face_strips;
            size_t num_constraints_fair_r = 3*num_halfedge_strips;

            // Index handling of the constraints
            std::map<std::string,size_t> cons_idx;

            cons_idx["vertex_1"] = 0;
            cons_idx["vertex_2"] = cons_idx["vertex_1"] + num_constraints_vert_1;
            cons_idx["edge_vec_1"] = cons_idx["vertex_2"] + num_constraints_vert_2;
            cons_idx["edge_vec_2"] = cons_idx["edge_vec_1"] + num_constraints_edge_vec_1;
            cons_idx["normal_1"] = cons_idx["edge_vec_2"] + num_constraints_edge_vec_2;
            cons_idx["normal_2"] = cond_idx["normal_1"] + num_constraints_normal_1;
            cons_idx["normal_3"] = cons_idx["normal_2"] + num_constraints_normal_2;
            cons_idx["ruling_0"] = cons_idx["normal_3"] + num_constraints_normal_3;
            cons_idx["ruling_1"] = cons_idx["ruling_0"] + num_constraints_ruling_0;
            cons_idx["ruling_2"] = cons_idx["ruling_1"] + num_constraints_ruling_1;
            cons_idx["dev"] = cons_idx["ruling_2"] + num_constraints_ruling_2;
            cons_idx["fair_v"] = cons_idx["dev"] + num_constraints_dev;
            cons_idx["fair_n"] = cons_idx["fair_v"] + num_constraints_fair_v;
            cons_idx["fair_r"] = cons_idx["fair_n"] + num_constraints_fair_n;
            cons_idx["num_cons"] = cons_idx["fair_r"] + num_constraints_fair_r;
        }

    void get_constraint_vertex_1(const VectorType &vars, VectorType &Dest){
        VectorView<ConfiguratorType> view(_quadTopol);

        view.set_vector(vars);

        if(Dest.size() != 3*_num_vertices){
            Dest.resize(3*_num_vertices);
        }

        for(int i = 0; i < _num_vertices; i++){
            Dest.segment(3*i,3) = view.reweighted_vertex(i) - view.vertex_weight(i)*view.vertex(i); 
        }
    }

    void get_constraint_vertex_2(const VectorType &vars, VectorType &Dest){

        if(Dest.size() != _num_vertices){
            Dest.resize(_num_vertices);
        }

        VectorView<ConfiguratorType> view(_quadTopol);


        for(int i = 0; i < _num_vertices; i++){
            Dest[i] = view.vertex_weight(i) - (view.dummy_weight(i)*view.dummy_weight(i)) - 1.0;
        }
    }

    void get_constraint_edge_vec_1(const VectorType &vars, VectorType &Dest){

        if(Dest.size() != 3*_num_faces){
            Dest.resize(3*_num_faces);
        }

        VectorView<ConfiguratorType> view(_quadTopol);

        view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            // First, obtain first vertex idces of the face
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w_0 = view.vertex_weight(node_0);
            RealType w_1 = view.vertex_weight(node_1);
            RealType w_2 = view.vertex_weight(node_2);
            RealType w_3 = view.vertex_weight(node_3);

            Dest.segment(3*i,3) = view.reweighted_edge_1(i) - (w_0 + w_1)*(view.reweighted_vertex(node_2)+view.reweighted_vertex(node_3)) + (w_2 + w_3)*(view.reweighted_vertex(node_1) + view.reweighted_vertex(node_0));
        }
    }

    void get_constraint_edge_vec_2(const VectorType &vars, VectorType &Dest){

        if(Dest.size() != 3*_num_faces){
            Dest.resize(3*_num_faces);
        }

        VectorView<ConfiguratorType> view(_quadTopol);
        view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            // First, obtain first vertex idces of the face
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w_0 = view.vertex_weight(node_0);
            RealType w_1 = view.vertex_weight(node_1);
            RealType w_2 = view.vertex_weight(node_2);
            RealType w_3 = view.vertex_weight(node_3);

            Dest.segment(3*i,3) = view.reweighted_edge_2(i) - (w_1 + w_2)*(view.reweighted_vertex(node_3)+view.reweighted_vertex(node_0)) + (w_0 + w_3)*(view.reweighted_vertex(node_1) + view.reweighted_vertex(node_2));
        }
    }

    void get_constraint_normal_1(const VectorType &vars, VectorType &Dest){

        if(Dest.size() != _num_faces){
            Dest.resize(_num_faces);
        }

        VectorView<ConfiguratorType> view(_quadTopol);
        view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            Dest[i] = view.face_normal(i).dot(view.reweighted_edge_1(i));
        }
    }

    void get_constraint_normal_2(const VectorType &vars, VectorType &Dest){

        if(Dest.size() != _num_faces){
            Dest.resize(_num_faces);
        }

        VectorView<ConfiguratorType> view(_quadTopol);
        view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            Dest[i] = view.face_normal(i).dot(view.reweighted_edge_2(i));
        }
    }

    void get_constraint_normal_3(const VectorType &vars, VectorType &Dest){

        if(Dest.size() != _num_faces){
            Dest.resize(_num_faces);
        }

        VectorView<ConfiguratorType> view(_quadTopol);
        view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            Dest[i] = view.face_normal(i).dot(view.face_normal(i)) - 1.0;
        }
    }

    void get_constraint_ruling_0(const VectorType &vars, VectorType &Dest){

        VectorView<ConfiguratorType> view(_quadTopol);
        view.set_vector(vars);

        size_t num_constraints = 3*(_num_halfedges - _quadTopol.getBdryHalfEdges().size());

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        for(SkippingBdryHalfEdgeIterator it(_quadTopol); it.valid(); ++it){
            int i = it.idx();
            int i_nobdry = it.idx_nobdry();

            int face_1 = _quadTopol.getFaceOfHalfEdge(i,0);
            int face_2 = _quadTopol.getFaceOfHalfEdge(i,1);

            Dest.segment(3*i_nobdry,3) = view.ruling(i_nobdry) - (view.face_normal(face_1).template head<3>()).cross(view.face_normal(face_2).template head<3>()); 
        }

    }

    void get_constraint_ruling_1(const VectorType &vars, VectorType &Dest){

        VectorView<ConfiguratorType> view(_quadTopol);
        view.set_vector(vars);

        std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

        size_t num_constraints = 3*(_num_faces - bdryFaces.size());

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        SkippingBdryHalfEdgeIterator it_he(_quadTopol);

        for(SkippingBdryFaceIterator it(_quadTopol); it.valid(); it++){

            int i = it.idx();
            int i_nobdry = it.id_nobdry();

            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
            int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
            int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
            int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

            RealType w0 = view.vertex_weight(node_0);
            RealType w1 = view.vertex_weight(node_1);
            RealType w2 = view.vertex_weight(node_2);
            RealType w3 = view.vertex_weight(node_3);

            VectorType ruling_01 = view.ruling(it_he.to_nbdry(heh_01));
            VectorType ruling_12 = view.ruling(it_he.to_nbdry(heh_12));
            VectorType ruling_23 = view.ruling(it_he.to_nbdry(heh_23));
            VectorType ruling_30 = view.ruling(it_he.to_nbdry(heh_30));

            Dest.segment(3*i_nobdry,3) = view.reweighted_ruling_1(i_nobdry) - (w1 + w2)*ruling_12 + (w3 + w0)*(ruling_30);
        }
    }

    void get_constraint_ruling_2(const VectorType &vars, VectorType &Dest){

        VectorView<ConfiguratorType> view(_quadTopol);

        view.set_vector(vars);

        std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

        size_t num_constraints = 3*(_num_faces - bdryFaces.size());

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        SkippingBdryHalfEdgeIterator it_he(_quadTopol);

        for(SkippingBdryFaceIterator it(_quadTopol); it.valid(); it++){

            int i = it.idx();
            int i_nobdry = it.id_nobdry();

            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
            int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
            int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
            int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

            RealType w0 = view.vertex_weight(node_0);
            RealType w1 = view.vertex_weight(node_1);
            RealType w2 = view.vertex_weight(node_2);
            RealType w3 = view.vertex_weight(node_3);

            VectorType ruling_01 = view.ruling(it_he.to_nbdry(heh_01));
            VectorType ruling_12 = view.ruling(it_he.to_nbdry(heh_12));
            VectorType ruling_23 = view.ruling(it_he.to_nbdry(heh_23));
            VectorType ruling_30 = view.ruling(it_he.to_nbdry(heh_30));

            Dest.segment(3*i_nobdry,3) = view.reweighted_ruling_2(i_nobdry) - (w0 + w1)*ruling_01 + (w2 + w3)*(ruling_23);
        }
    }

    void get_constraint_dev(const VectorType &vars, VectorType &Dest){

        VectorView<ConfiguratorType> view(_quadTopol);

        view.set_vector(vars);

        std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

        size_t num_constraints = 3*(_num_faces - bdryFaces.size());

        if(Dest.size() != num_constraints){
            Dest.resize(num_constraints);
        }

        Dest.setZero();

        for(SkippingBdryFaceIterator it(_quadTopol); it.valid(); it++){
            int i = it.idx();
            int i_nobdry = it.id_nobdry();
            Dest.segment(3*i_nobdry,3) = (view.reweighted_ruling_1(i_nobdry). template head<3>()).cross(view.reweighted_ruling_2(i_nobdry). template head<3>());
        }
    }

    void get_constraint_fair_v(const VectorType &vars, VectorType &Dest)
    {
        VectorView<ConfiguratorType> view(_quadTopol);
        view.set_vector(vars);
        auto vertex_strips = _quadTopol.getVertexStrips();

        int num_strips = vertex_strips.size();
        
        if(Dest.size() != 3*num_strips){
            Dest.resize(3*num_strips);
        }

        Dest.setZero();

        for(int i = 0; i < num_strips; i++){
            int node_0 = std::get<0>(vertex_strips[i]);
            int node_1 = std::get<1>(vertex_strips[i]);
            int node_2 = std::get<2>(vertex_strips[i]);
            Dest.segment(3*i,3) = view.vertex(node_0) - 2*view.vertex(node_1) + view.vertex(node_2);
        }
    }

    void get_constraint_fair_n(const VectorType &vars, VectorType &Dest)
    {
        VectorView<ConfiguratorType> view(_quadTopol);
        view.set_vector(vars);

        auto face_strips = _quadTopol.getFaceStrips();

        size_t num_strips = face_strips.size();

        if(Dest.size() != 3*num_strips)
        {
            Dest.resize(3*num_strips);
        }

        for(int i = 0; i < num_strips; i++)
        {
            int face_0 = std::get<0>(face_strips[i]);
            int face_1 = std::get<1>(face_strips[i]);
            int face_2 = std::get<2>(face_strips[i]);

            VectorType normal_0 = view.face_normal(face_0);
            VectorType normal_1 = view.face_normal(face_1);
            VectorType normal_2 = view.face_normal(face_2);

            Dest.segment(3*i,3) = normal_0 - 2*normal_1 + normal_2;
        }

    }

    void get_constraint_fair_r(const VectorType &vars, VectorType &Dest)
    {
        VectorView<ConfiguratorType> view(_quadTopol);
        view.set_vector(vars);

        auto halfedge_strips = _quadTopol.getHalfEdgeStrips();

        SkippingBdryHalfEdgeIterator he_it(_quadTopol);

        size_t num_strips = halfedge_strips.size();

        if(Dest.size() != 3*num_strips)
        {
            Dest.resize(3*num_strips);
        }

        Dest.setZero();

        for(int i = 0; i < num_strips; i++)
        {
            int he_0 = std::get<0>(halfedge_strips[i]);
            int he_1 = std::get<1>(halfedge_strips[i]);
            int he_2 = std::get<2>(halfedge_strips[i]);

            VectorType ruling_0 = view.ruling(he_it.to_nbdry(he_0));
            VectorType ruling_1 = view.ruling(he_it.to_nbdry(he_1));
            VectorType ruling_2 = view.ruling(he_it.to_nbdry(he_2));

            Dest.segment(3*i,3) = ruling_0 - 2*ruling_1 + ruling_2;
        }
    }
    
    void get_constraints(const VectorType &vars, VectorType &Dest)
    {
        _view.set_vector(vars);

        if(Dest.size() != cons_idx["num_cons"])
        {
            Dest.resize(cons_idx["num_cons"]);
        }

        Dest.setZero();

        // First, iterate over all vertices
        for(int i = 0; i < _num_vertices; i++)
        {
            // constraint vertex_1
            Dest.segment(cons_idx["vertex_1"] + 3*i,3) = view.reweighted_vertex(i) - view.vertex_weight(i)*view.vertex(i); 
            // constraint vertex_2
            Dest[cons_idx["vertex_2"] + i] = view.vertex_weight(i) - (view.dummy_weight(i)*view.dummy_weight(i)) - 1.0;
        }

        // Next, iterate over all faces
        for(int i = 0; i < _num_faces; i++){
              // First, obtain first vertex idces of the face
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w_0 = view.vertex_weight(node_0);
            RealType w_1 = view.vertex_weight(node_1);
            RealType w_2 = view.vertex_weight(node_2);
            RealType w_3 = view.vertex_weight(node_3);

            // edge vec 1 constraint
            Dest.segment(cons_idx["edge_vec_1"] + 3*i,3) = view.reweighted_edge_1(i) - (w_0 + w_1)*(view.reweighted_vertex(node_2)+view.reweighted_vertex(node_3)) + (w_2 + w_3)*(view.reweighted_vertex(node_1) + view.reweighted_vertex(node_0));
            // edge vec 2 constraint
            Dest.segment(cons_idx["edge_vec_2"] + 3*i,3) = view.reweighted_edge_2(i) - (w_1 + w_2)*(view.reweighted_vertex(node_3)+view.reweighted_vertex(node_0)) + (w_0 + w_3)*(view.reweighted_vertex(node_1) + view.reweighted_vertex(node_2));
            // normal 1 constraint
            Dest[cons_idx["normal_1"] + i] = view.face_normal(i).dot(view.reweighted_edge_1(i));
            // normal 2 constraint
            Dest[cons_idx["normal_2"] + i] = view.face_normal(i).dot(view.reweighted_edge_2(i));
            // normal 3 constraint
            Dest[cons_idx["normal_3"] + i] = view.face_normal(i).dot(view.face_normal(i)) - 1.0;
        }

        // Next, iterate over all nonboundary halfedges
        SkippingBdryHalfEdgeIterator he_it(_quadTopol);
        for(he_it; he_it.valid(); he_it++)
        {
            int i = it.idx();
            int i_nobdry = it.idx_nobdry();

            int face_1 = _quadTopol.getFaceOfHalfEdge(i,0);
            int face_2 = _quadTopol.getFaceOfHalfEdge(i,1);

            // Constraint ruling0
            Dest.segment(cons_idx["ruling_0"]+3*i_nobdry,3) = view.ruling(i_nobdry) - (view.face_normal(face_1).template head<3>()).cross(view.face_normal(face_2).template head<3>()); 
        }

        he_it.reset();

        // Next, iterate over all nonboundary faces
        SkippingBdryFaceIterator face_it;
        for(face_it; face_it.valid(); face_it++)
        {
            int i = it.idx();
            int i_nobdry = it.id_nobdry();

            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
            int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
            int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
            int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

            RealType w0 = view.vertex_weight(node_0);
            RealType w1 = view.vertex_weight(node_1);
            RealType w2 = view.vertex_weight(node_2);
            RealType w3 = view.vertex_weight(node_3);

            VectorType ruling_01 = view.ruling(it_he.to_nbdry(heh_01));
            VectorType ruling_12 = view.ruling(it_he.to_nbdry(heh_12));
            VectorType ruling_23 = view.ruling(it_he.to_nbdry(heh_23));
            VectorType ruling_30 = view.ruling(it_he.to_nbdry(heh_30));

            // ruling 1 constraint
            Dest.segment(cons_idx["ruling_1"] + 3*i_nobdry,3) = view.reweighted_ruling_1(i_nobdry) - (w1 + w2)*ruling_12 + (w3 + w0)*(ruling_30);
            // ruling 2 constraint
            Dest.segment(cons_idx["ruling_2"] + 3*i_nobdry,3) = view.reweighted_ruling_2(i_nobdry) - (w0 + w1)*ruling_01 + (w2 + w3)*(ruling_23);
            // developability constraint
            Dest.segment(cons_idx["dev"] + 3*i_nobdry,3) = (view.reweighted_ruling_1(i_nobdry). template head<3>()).cross(view.reweighted_ruling_2(i_nobdry). template head<3>());
        }

        for(int i = 0; i < num_vertex_strips; i++){
            int node_0 = std::get<0>(vertex_strips[i]);
            int node_1 = std::get<1>(vertex_strips[i]);
            int node_2 = std::get<2>(vertex_strips[i]);
            Dest.segment(cons_idx["fair_v"] + 3*i,3) = view.vertex(node_0) - 2*view.vertex(node_1) + view.vertex(node_2);
        }

        for(int i = 0; i < num_face_strips; i++)
        {
            int face_0 = std::get<0>(face_strips[i]);
            int face_1 = std::get<1>(face_strips[i]);
            int face_2 = std::get<2>(face_strips[i]);

            VectorType normal_0 = view.face_normal(face_0);
            VectorType normal_1 = view.face_normal(face_1);
            VectorType normal_2 = view.face_normal(face_2);

            Dest.segment(cons_idx["fair_n"] + 3*i,3) = normal_0 - 2*normal_1 + normal_2;
        }

        for(int i = 0; i < num_halfedge_strips; i++)
        {
            int he_0 = std::get<0>(halfedge_strips[i]);
            int he_1 = std::get<1>(halfedge_strips[i]);
            int he_2 = std::get<2>(halfedge_strips[i]);

            VectorType ruling_0 = view.ruling(he_it.to_nbdry(he_0));
            VectorType ruling_1 = view.ruling(he_it.to_nbdry(he_1));
            VectorType ruling_2 = view.ruling(he_it.to_nbdry(he_2));

            Dest.segment(cons_idx["fair_r"] + 3*i,3) = ruling_0 - 2*ruling_1 + ruling_2;
        }

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
        size_t _num_vertices;
        size_t _num_faces;
        size_t _num_edges;
        size_t _num_halfedges;
        QuadMeshTopologySaver &_quadTopol;

    public:
        ConstraintGrad(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
        {


        }

        //void get_constraintgrad_vertex_1(const VectorType &vars, MatrixType &Dest){
        void get_constraintgrad_vertex_1(const VectorType &vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            size_t num_dofs = view._idx["num_dofs"];

            if(Dest.rows() != 3*_num_vertices || Dest.cols() != num_dofs){
                Dest.resize(3*_num_vertices, num_dofs);
            }

            Dest.setZero();
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_vertices; i++){
                RealType w_i = view.vertex_weight(i);
                std::vector<Eigen::Triplet<RealType>> triplet;
                VectorType vec = -view.vertex(i);
                MatrixType vec_mat = vectorToSparseMat(vec);
                assignSparseBlockInplace(Dest,-w_i*Id,3*i,3*i + view._idx["vertices"],triplet);
                assignSparseBlockInplace(Dest,vec_mat,3*i,3*i + view._idx["weights"],triplet);
                assignSparseBlockInplace(Dest,Id,3*i,3*i + view._idx["reweighted_vertices"],triplet);
            }

        }

        void get_constraintgrad_vertex_2(const VectorType &vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            size_t num_dofs = view._idx["num_dofs"];

            if(Dest.rows()!= _num_vertices || Dest.cols() != num_dofs){
                Dest.resize(_num_vertices, num_dofs);
            }
            
            for(int i = 0; i < _num_vertices; i++){
                RealType dummy_i = view.dummy_weight(i);
                // Assign values directly
                Dest.coeffRef(i, i + view._idx["dummy_weights"]) = -2.0 * dummy_i;
                Dest.coeffRef(i, i + view._idx["weights"]) = 1.0;
            }
        }

        void get_constraintgrad_edge_vec_1(const VectorType &vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            size_t num_dofs = view._idx["num_dofs"];

            if(Dest.rows() != 3*_num_faces || Dest.cols() != num_dofs){
                Dest.resize(3*_num_faces, num_dofs);
            }

            Dest.setZero();

            // From average edge vector 1
            auto Id = MatrixType(3,3);
            Id.setIdentity();
            std::vector<Eigen::Triplet<RealType>> triplet;
            assignSparseBlockInplace(Dest,Id,0,view._idx["reweighted_edges_1"],triplet);

            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                RealType w0 = view.vertex_weight(node_0);
                RealType w1 = view.vertex_weight(node_1);
                RealType w2 = view.vertex_weight(node_2);
                RealType w3 = view.vertex_weight(node_3);

                VectorType v0_tilde = view.reweighted_vertex(node_0);
                VectorType v1_tilde = view.reweighted_vertex(node_1);
                VectorType v2_tilde = view.reweighted_vertex(node_2);
                VectorType v3_tilde = view.reweighted_vertex(node_3);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                MatrixType sum_1 = vectorToSparseMat((v2_tilde + v3_tilde));
                MatrixType sum_2 = vectorToSparseMat((v1_tilde + v0_tilde));
                assignSparseBlockInplace(Dest, -sum_1,3*i,view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -sum_1,3*i,view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, sum_2,3*i,view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, sum_2,3*i,view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseBlockInplace(Dest, (w2+w3)*Id,3*i,view._idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseBlockInplace(Dest, (w2+w3)*Id,3*i,view._idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseBlockInplace(Dest, -(w0+w1)*Id,3*i,view._idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseBlockInplace(Dest, -(w0+w1)*Id,3*i,view._idx["reweighted_vertices"] + 3*node_3,triplet);
            }

        }

        void get_constraintgrad_edge_vec_2(const VectorType &vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            size_t num_dofs = view._idx["num_dofs"];

            if(Dest.rows() != 3*_num_vertices || Dest.cols() != num_dofs){
                Dest.resize(3*_num_faces, num_dofs);
            }

            Dest.setZero();

            std::vector<Eigen::Triplet<RealType>> triplet;
            // From average edge vector 1
            auto Id = MatrixType(3,3);
            Id.setIdentity();
            assignSparseBlockInplace(Dest,Id,0,view._idx["reweighted_edges_2"],triplet);

            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                RealType w0 = view.vertex_weight(node_0);
                RealType w1 = view.vertex_weight(node_1);
                RealType w2 = view.vertex_weight(node_2);
                RealType w3 = view.vertex_weight(node_3);

                VectorType v0_tilde = view.reweighted_vertex(node_0);
                VectorType v1_tilde = view.reweighted_vertex(node_1);
                VectorType v2_tilde = view.reweighted_vertex(node_2);
                VectorType v3_tilde = view.reweighted_vertex(node_3);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                MatrixType sum_1 = vectorToSparseMat((v2_tilde + v1_tilde));
                MatrixType sum_2 = vectorToSparseMat((v0_tilde + v3_tilde));
    
                assignSparseBlockInplace(Dest, sum_1,3*i,view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -sum_2,3*i,view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, -sum_2,3*i,view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, sum_1,3*i,view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                assignSparseBlockInplace(Dest,-(w1+w2)*Id,3*i,view._idx["reweighted_vertices"] + 3*node_0,triplet);
                assignSparseBlockInplace(Dest,(w0+w3)*Id,3*i,view._idx["reweighted_vertices"] + 3*node_1,triplet);
                assignSparseBlockInplace(Dest,(w0+w3)*Id,3*i,view._idx["reweighted_vertices"] + 3*node_2,triplet);
                assignSparseBlockInplace(Dest,-(w1+w2)*Id,3*i,view._idx["reweighted_vertices"] + 3*node_3,triplet);
            }

        }

        void get_constraintgrad_normal_1(const VectorType &vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            size_t num_dofs = view._idx["num_dofs"];

            if(Dest.rows() != _num_faces || Dest.cols() != num_dofs ){
                Dest.resize(_num_faces, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                MatrixType vec_1 = vectorToSparseMat(view.reweighted_edge_1(i),true);
                MatrixType vec_2 = vectorToSparseMat(view.face_normal(i),true);
                assignSparseBlockInplace(Dest,vec_1,i,view._idx["normals"] + 3*i,triplet);
                assignSparseBlockInplace(Dest,vec_2,i,view._idx["reweighted_edges_1"] + 3*i,triplet);
            }
        }

        void get_constraintgrad_normal_2(const VectorType &vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            size_t num_dofs = view._idx["num_dofs"];

            if(Dest.rows() != _num_faces || Dest.cols() != num_dofs ){
                Dest.resize(_num_faces, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                MatrixType vec_1 = vectorToSparseMat(view.reweighted_edge_2(i),true);
                MatrixType vec_2 = vectorToSparseMat(view.face_normal(i),true);
                assignSparseBlockInplace(Dest,vec_1,i,view._idx["normals"] + 3*i,triplet);
                assignSparseBlockInplace(Dest,vec_2,i,view._idx["reweighted_edges_2"] + 3*i,triplet);
            }
        }

        void get_constraintgrad_normal_3(const VectorType &vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            size_t num_dofs = view._idx["num_dofs"];

            if(Dest.rows() != _num_faces || Dest.cols() != num_dofs ){
                Dest.resize(_num_faces, num_dofs);
            }

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(int i = 0; i < _num_faces; i++){
                // Set the derivatives w.r.t. the normals
                MatrixType vec_1 = vectorToSparseMat(view.face_normal(i),true);
                assignSparseBlockInplace(Dest,2.0*vec_1,i,view._idx["normals"] + i,triplet);
            }
        }
        
        void get_constraintgrad_ruling_0(const VectorType &vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            std::vector<int> bdryHalfEdges = _quadTopol.getBdryHalfEdges();

            size_t num_dofs = view._idx["num_dofs"];

            size_t num_constraints = 3*(_num_halfedges - bdryHalfEdges.size());

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            Dest.setZero();

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(SkippingBdryHalfEdgeIterator it(_quadTopol); it.valid();it++)
            {
                int i = it.idx();
                int i_nobdry = it.idx_nobdry();

                size_t faceIdx_1 = _quadTopol.getFaceOfHalfEdge(i,0);
                size_t faceIdx_2 = _quadTopol.getFaceOfHalfEdge(i,1);
                Eigen::Vector3d normal_1 = view.face_normal(faceIdx_1).template head<3>();
                Eigen::Vector3d normal_2 = view.face_normal(faceIdx_2).template head<3>();

                // First, set the derivative w.r.t. the ruling vector
                auto Id = MatrixType(3,3);
                Id.setIdentity();

                assignSparseBlockInplace(Dest,Id,3*(i_nobdry),view._idx["rulings"] + 3*i_nobdry,triplet);

                // Now, set the derivatives w.r.t. the cross product
                // First, differentiate w.r.t. normal_1
                Eigen::Vector3d e;
                MatrixType vec_prod;
                e.setZero();
                e[0] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,3*(i_nobdry),view._idx["normals"] + 3*faceIdx_1,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,3*(i_nobdry),view._idx["normals"] + 3*faceIdx_1 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = convertVecToSparseMat(-e.cross(normal_2));
                assignSparseBlockInplace(Dest,vec_prod,3*(i_nobdry),view._idx["normals"] + 3*faceIdx_1 + 2,triplet);

                // Now, differentiate w.r.t. normal_2
                e[2] = 0.0;
                e[0] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,3*(i_nobdry),view._idx["normals"] + 3*faceIdx_2,triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,3*(i_nobdry),view._idx["normals"] + 3*faceIdx_2 + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                vec_prod = convertVecToSparseMat(e.cross(normal_1));
                assignSparseBlockInplace(Dest,vec_prod,3*(i_nobdry),view._idx["normals"] + 3*faceIdx_2 + 2,triplet);
            }
        }

        void get_constraintgrad_ruling_1(const VectorType& vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);

            view.set_vector(vars);

            std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

            size_t num_dofs = view._idx["num_dofs"];

            size_t num_constraints = 3*(_num_faces - bdryFaces.size());

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            Dest.setZero();

            std::vector<Eigen::Triplet<RealType>> triplet;
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            SkippingBdryHalfEdgeIterator it_he(_quadTopol);

            for(SkippingBdryFaceIterator it(_quadTopol); it.valid(); it++){

                int i = it.idx();
                int i_nobdry = it.id_nobdry();

                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                auto node_0_vert = view.vertex(node_0);
                auto node_1_vert = view.vertex(node_1);
                auto node_2_vert = view.vertex(node_2);
                auto node_3_vert = view.vertex(node_3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                RealType w0 = view.vertex_weight(node_0);
                RealType w1 = view.vertex_weight(node_1);
                RealType w2 = view.vertex_weight(node_2);
                RealType w3 = view.vertex_weight(node_3);

                VectorType ruling_01 = view.ruling(it_he.to_nbdry(heh_01));
                VectorType ruling_12 = view.ruling(it_he.to_nbdry(heh_12));
                VectorType ruling_23 = view.ruling(it_he.to_nbdry(heh_23));
                VectorType ruling_30 = view.ruling(it_he.to_nbdry(heh_30));

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}

                assignSparseBlockInplace(Dest,Id,3*i_nobdry,view._idx["reweighted_rulings_1"] + 3*i_nobdry,triplet);

                MatrixType ruling_30_mat = vectorToSparseMat(ruling_30);
                MatrixType ruling_12_mat = vectorToSparseMat(ruling_12);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseBlockInplace(Dest, ruling_30_mat,3*i_nobdry,view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -ruling_12_mat,3*i_nobdry,view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, -ruling_12_mat,3*i_nobdry,view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, ruling_30_mat,3*i_nobdry,view._idx["weights"] + node_3,triplet);

                std::cout<<"Compare: "<<std::endl;
                std::cout<<view._idx["rulings"] + 3*heh_12<<"; "<<view._idx["rulings"] + 3*it_he.to_nbdry(heh_12)<<std::endl;
                std::cout<<view._idx["rulings"] + 3*heh_30<<"; "<<view._idx["rulings"] + 3*it_he.to_nbdry(heh_30)<<std::endl;

                //Derivatives w.r.t. the rulings
                assignSparseBlockInplace(Dest,-(w1 + w2)*Id,3*i_nobdry,view._idx["rulings"] + 3*it_he.to_nbdry(heh_12),triplet);
                assignSparseBlockInplace(Dest,(w3 + w0)*Id,3*i_nobdry,view._idx["rulings"] + 3*it_he.to_nbdry(heh_30),triplet);
            }
        }

        void get_constraintgrad_ruling_2(const VectorType& vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);

            view.set_vector(vars);

            std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

            size_t num_dofs = view._idx["num_dofs"];

            size_t num_constraints = 3*(_num_faces - bdryFaces.size());

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            Dest.setZero();

            std::vector<Eigen::Triplet<RealType>> triplet;
            auto Id = MatrixType(3,3);
            Id.setIdentity();

            SkippingBdryHalfEdgeIterator it_he(_quadTopol);

            for(SkippingBdryFaceIterator it(_quadTopol); it.valid(); it++){

                int i = it.idx();
                int i_nbdry = it.id_nobdry();

                int node_0 = _quadTopol.getNodeOfQuad(i,0);
                int node_1 = _quadTopol.getNodeOfQuad(i,1);
                int node_2 = _quadTopol.getNodeOfQuad(i,2);
                int node_3 = _quadTopol.getNodeOfQuad(i,3);

                int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
                int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
                int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
                int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

                RealType w0 = view.vertex_weight(node_0);
                RealType w1 = view.vertex_weight(node_1);
                RealType w2 = view.vertex_weight(node_2);
                RealType w3 = view.vertex_weight(node_3);

                VectorType ruling_01 = view.ruling(it_he.to_nbdry(heh_01));
                VectorType ruling_23 = view.ruling(it_he.to_nbdry(heh_23));

                MatrixType ruling_01_mat = vectorToSparseMat(ruling_01);
                MatrixType ruling_23_mat = vectorToSparseMat(ruling_23);

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}
                assignSparseBlockInplace(Dest,Id,3*i_nbdry,view._idx["reweighted_rulings_2"] + 3*i_nbdry,triplet);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                assignSparseBlockInplace(Dest, -ruling_01_mat,3*i_nbdry,view._idx["weights"] + node_0,triplet);
                assignSparseBlockInplace(Dest, -ruling_01_mat,3*i_nbdry,view._idx["weights"] + node_1,triplet);
                assignSparseBlockInplace(Dest, ruling_23_mat,3*i_nbdry,view._idx["weights"] + node_2,triplet);
                assignSparseBlockInplace(Dest, ruling_23_mat,3*i_nbdry,view._idx["weights"] + node_3,triplet);

                //Derivatives w.r.t. the rulings
                assignSparseBlockInplace(Dest,-(w0 + w1)*Id,3*i_nbdry,view._idx["rulings"] + 3*it_he.to_nbdry(heh_01),triplet);
                assignSparseBlockInplace(Dest,(w2 + w3)*Id,3*i_nbdry,view._idx["rulings"] + 3*it_he.to_nbdry(heh_23),triplet);
            }
        }

        void get_constraintgrad_dev(const VectorType &vars, MatrixType &Dest){

            VectorView<ConfiguratorType> view(_quadTopol);

            view.set_vector(vars);

            std::vector<int> bdryFaces = _quadTopol.getBdryFaces();

            size_t num_constraints = 3*(_num_faces - bdryFaces.size());
            size_t num_dofs = view._idx["num_dofs"];

            if(Dest.rows() != num_constraints || Dest.cols() != num_dofs){
                Dest.resize(num_constraints, num_dofs);
            }

            Dest.setZero();

            std::vector<Eigen::Triplet<RealType>> triplet;

            for(SkippingBdryFaceIterator it(_quadTopol); it.valid(); it++)
            {

                int i = it.idx();
                int i_nobdry = it.id_nobdry();

                // obtain the reweighted ruling vectors
                Eigen::Vector3d reweighted_ruling_1 = view.reweighted_ruling_1(i).template head<3>();
                Eigen::Vector3d reweighted_ruling_2 = view.reweighted_ruling_2(i).template head<3>();

                // Derivatives w.r.t. the first ruling vector \Tilde{r}_{13,f}
                Eigen::Vector3d e;
                e.setZero();
                e[0] = 1.0; 
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),3*i_nobdry,3*i_nobdry + view._idx["reweighted_rulings_1"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),3*i_nobdry,3*i_nobdry + view._idx["reweighted_rulings_1"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(e.cross(reweighted_ruling_2)),3*i_nobdry,3*i_nobdry + view._idx["reweighted_rulings_1"] + 2,triplet);
            
                // Derivatives w.r.t. the second ruling vector \Tilde{r}_{02,f}
                e[2] = 0.0;
                e[0] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),3*i_nobdry,3*i_nobdry + view._idx["reweighted_rulings_2"],triplet);
                e[0] = 0.0;
                e[1] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),3*i_nobdry,3*i_nobdry + view._idx["reweighted_rulings_2"] + 1,triplet);
                e[1] = 0.0;
                e[2] = 1.0;
                assignSparseBlockInplace(Dest,vectorToSparseMat(reweighted_ruling_1.cross(e)),3*i_nobdry,3*i_nobdry + view._idx["reweighted_rulings_2"] + 2,triplet);
            }
        }

        void get_constraintgrad_fair_v(const VectorType &vars, MatrixType &Dest){
            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            size_t num_dofs = view._idx["num_dofs"];

            auto vertex_strips = _quadTopol.getVertexStrips();

            if(Dest.rows() != 3*vertex_strips.size() || Dest.cols() != num_dofs)
            {
                Dest.resize(3*vertex_strips.size(), num_dofs);
            }

            MatrixType Id(3,3);
            Id.setIdentity();
            std::vector<Eigen::Triplet<double>> triplet;

            for(int i = 0; i < vertex_strips.size(); i++){
                int node_0 = std::get<0>(vertex_strips[i]);
                int node_1 = std::get<1>(vertex_strips[i]);
                int node_2 = std::get<2>(vertex_strips[i]);
                
                assignSparseBlockInplace(Dest,Id,3*i,3*node_0 + view._idx["vertices"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,3*i,3*node_1 + view._idx["vertices"], triplet);
                assignSparseBlockInplace(Dest, Id, 3*i,3*node_2 + view._idx["vertices"], triplet);
            }
        }

        void get_constraintgrad_fair_n(const VectorType &vars, MatrixType &Dest){
            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            auto normal_strips = _quadTopol.getFaceStrips();

            int num_strips = normal_strips.size();
            int num_dofs = view._idx["num_dofs"];

            if(Dest.rows() != 3*num_strips || Dest.cols() != num_dofs){
                Dest.resize(3*num_strips, num_dofs);
            }

            MatrixType Id(3,3);
            Id.setIdentity();

            std::vector<Eigen::Triplet<double>> triplet;

            for(int i = 0; i < num_strips; i++){
                int normal_0 = std::get<0>(normal_strips[i]);
                int normal_1 = std::get<1>(normal_strips[i]);
                int normal_2 = std::get<2>(normal_strips[i]);

                assignSparseBlockInplace(Dest,Id,3*i,3*normal_0 + view._idx["normals"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,3*i,3*normal_1 + view._idx["normals"],triplet);
                assignSparseBlockInplace(Dest,Id,3*i,3*normal_2 + view._idx["normals"],triplet);
            }
        }

        void get_constraintgrad_fair_r(const VectorType &vars, MatrixType &Dest){
            VectorView<ConfiguratorType> view(_quadTopol);
            view.set_vector(vars);

            auto he_strips = _quadTopol.getHalfEdgeStrips();

            int num_strips = he_strips.size();
            int num_dofs = view._idx["num_dofs"];

            if(Dest.rows() != 3*num_strips || Dest.cols() != num_dofs)
            {
                Dest.resize(3*num_strips, num_dofs);
            }

            MatrixType Id(3,3);
            Id.setIdentity();

            std::vector<Eigen::Triplet<RealType>> triplet;
            SkippingBdryHalfEdgeIterator he_it(_quadTopol);

            for(int i = 0; i < num_strips; i++)
            {
                int heh_0 = std::get<0>(he_strips[i]);
                int heh_1 = std::get<1>(he_strips[i]);
                int heh_2 = std::get<2>(he_strips[i]);

                assignSparseBlockInplace(Dest,Id,3*i,3*he_it.to_nbdry(heh_0) + view._idx["rulings"],triplet);
                assignSparseBlockInplace(Dest,-2*Id,3*i,3*he_it.to_nbdry(heh_1) + view._idx["rulings"],triplet);
                assignSparseBlockInplace(Dest,Id,3*i,3*he_it.to_nbdry(heh_2) + view._idx["rulings"],triplet);
            }
        }
};