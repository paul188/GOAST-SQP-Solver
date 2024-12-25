#include <goast/Quads.h>
#include <goast/QuadMesh/QuadTopology.h>
#include <SQP/Utils/SparseMat.h>
#include <cmath>

template<typename ConfiguratorType>
struct Variables{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    FullMatrixType _vertices;
    VectorType _vertex_weights;
    FullMatrixType _reweighted_vertices;
    VectorType _dummy_weights;
    FullMatrixType _rulings;
    FullMatrixType _face_normals;
    FullMatrixType _reweighted_edges_1;  // Tilde{e}_{0,2}
    FullMatrixType _reweighted_edges_2;  // Tilde{e}_{1,3}
    FullMatrixType _reweighted_rulings_1;   // Tilde{r}_{1,3}
    FullMatrixType _reweighted_rulings_2;   // Tilde{r}_{0,2}

    /* Initialize like follows:
       -> vertices to quad vertices
       -> vertex_weights to 2.0
       -> reweighted_vertices to vertices*vertex_weights
       -> dummy_weights to 1.0
       -> rulings to 0.0
       -> face normals to 0.0
       -> reweighted_edges to 0.0
         -> reweighted_rulings to 0.0
    */
    Variables(QuadMesh &qm, QuadMeshTopologySaver &quadTopol)
    {
        int nVertices = quadTopol.getNumVertices();
        int nFaces = quadTopol.getNumFaces();
        int nEdges = quadTopol.getNumEdges();

        _vertices.resize(3,nVertices);
        _vertex_weights.resize(nVertices);
        _dummy_weights.resize(nVertices);
        _rulings.resize(3,2*nEdges);
        _face_normals.resize(3,nFaces);
        _reweighted_edges_1.resize(3,nFaces);
        _reweighted_edges_2.resize(3,nFaces);
        _reweighted_rulings_1.resize(3,nFaces);
        _reweighted_rulings_2.resize(3,nFaces);

        getGeometryMat<FullMatrixType>(qm._mesh, _vertices);

        _vertex_weights.setConstant(2.0);
        _dummy_weights.setConstant(1.0);

        _reweighted_vertices = 2.0*_vertices;

        _rulings.setZero();
        _face_normals.setZero();
        _reweighted_edges_1.setZero();
        _reweighted_edges_2.setZero();
        _reweighted_rulings_1.setZero();
        _reweighted_rulings_2.setZero();
    }

    setQuadMeshGeometry(QuadMesh &qm, QuadMeshTopologySaver &quadTopol, const Variables& vars){
        for(int i = 0; i < quadTopol.getNumVertices(); i++){
            qm._mesh.set_point(qm._mesh.vertex_handle(i), Point(vars._vertices[0,i], vars._vertices[1,i], vars._vertices[2,i]));
        }
    }
};

template <typename ConfiguratorType>
class Constraint{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    public:
    // How to structure the constraints Dest ?
    // First the x, then the y and then the z components of the constraints
    // This way, should be compatible with structure of the geometry
    void get_constraint_vertex_1(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        size_t num_vertices = qm.num_vertices();

        // Check that all sizes are correct:
        if( vars._vertices.size() != num_vertices ){
            std::cerr << "Wrong number of vertices in the variables to calculate constraints"<<std::endl;
        }

        if( vars._vertex_weights.size() != num_vertices ){
            std::cerr << "Wrong number of vertex weights in the variables to calculate constraints"<<std::endl;
        }

        if( vars._reweighted_vertices.size() != num_vertices ){
            std::cerr << "Wrong number of reweighted vertices in the variables to calculate constraints"<<std::endl;
        }

        if(Dest.size() != 3*num_vertices){
            Dest.resize(3*num_vertices);
        }

        for(int i = 0; i < num_vertices; i++){
            for(int j = 0; j < 3; j++){
                Dest[i*3 + j] = vars._reweighted_vertices[j,i] - vars._vertex_weights[i]*vars._vertices[j,i];
            }
        }
    }

    void get_constraint_vertex_2(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        size_t num_vertices = qm.num_vertices();

        if(vars._dummy_weights.size() != num_vertices ){
            std::cerr << "Wrong number of dummy weights in the variables to calculate constraints"<<std::endl;
        }

        if(Dest.size() != num_vertices){
            Dest.resize(num_vertices);
        }

        for(int i = 0; i < num_vertices; i++){
            Dest[i] = vars._vertex_weights[i] - (vars._dummy_weights[i]*vars._dummy_weights[i]) - 1.0;
        }
    }

    void get_constraint_edge_vec_1(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        size_t num_faces = qm.num_faces();

        if(vars._reweighted_edges.size() != 2*num_faces){
            std::cerr << "Wrong number of reweighted edges in the variables to calculate constraints"<<std::endl;
        }

        if(Dest.size() != 3*num_faces){
            Dest.resize(3*num_faces);
        }

        for(int i = 0; i < num_faces; i++){
                // First, obtain first vertex idces of the face
            int node_0 = quadTopol.getNodeOfQuad(i,0);
            int node_1 = quadTopol.getNodeOfQuad(i,1);
            int node_2 = quadTopol.getNodeOfQuad(i,2);
            int node_3 = quadTopol.getNodeOfQuad(i,3);

            VectorType face_constraint = vars._reweighted_edges_1[:,i] - (vars._vertex_weights[node_0] + vars._vertex_weights[node_1])*(vars._reweighted_vertices[:,node_2]+vars._reweighted_vertices[:,node_3]) + (vars._vertex_weights[node_2] + vars._vertex_weights[node_3])*(vars._reweighted_vertices[:,node_1] + vars._reweighted_vertices[:,node_0]);

            Dest.segment(3*i,3) = face_constraint;
        }
    }

    void get_constraint_edge_vec_2(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        size_t num_faces = qm.num_faces();

        if(vars._reweighted_edges.size() != 2*num_faces){
            std::cerr << "Wrong number of rulings in the variables to calculate constraints"<<std::endl;
        }

        if(Dest.size() != 3*num_faces){
            Dest.resize(3*num_faces);
        }

        for(int i = 0; i < num_faces; i++){
            // First, obtain first vertex idces of the face
            int node_0 = quadTopol.getNodeOfQuad(i,0);
            int node_1 = quadTopol.getNodeOfQuad(i,1);
            int node_2 = quadTopol.getNodeOfQuad(i,2);
            int node_3 = quadTopol.getNodeOfQuad(i,3);

            VectorType face_constraint = vars._reweighted_edges_2[:,i] - (vars._vertex_weights[1] + vars._vertex_weights[node_2])*(vars._reweighted_vertices[:,0]+vars._reweighted_vertices[:,3]) - (vars._vertex_weights[node_0] + vars._vertex_weights[node_3])*(vars._reweighted_vertices[:,2] + vars._reweighted_vertices[:,1]);

            Dest.segment(3*i,3) = face_constraint;
        }
    }

    void get_constraint_normal_1(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        size_t num_faces = qm.num_faces();

        if(vars._face_normals.size() != num_faces){
            std::cerr << "Wrong number of face normals in the variables to calculate constraints"<<std::endl;
        }

        if(Dest.size() != num_faces){
            Dest.resize(num_faces);
        }

        for(int i = 0; i < num_faces; i++){
            Dest[i] = vars._face_normals.col(i).dot(vars._reweighted_edges_1.col(i));
        }
    }

    void get_constraint_normal_2(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        size_t num_faces = qm.num_faces();

        if(vars._face_normals.size() != num_faces){
            std::cerr << "Wrong number of face normals in the variables to calculate constraints"<<std::endl;
        }

        if(Dest.size() != num_faces){
            Dest.resize(num_faces);
        }

        for(int i = 0; i < num_faces; i++){
            Dest[i] = vars._face_normals.col(i).dot(vars._reweighted_edges_2.col(i));
        }
    }

    void get_constraint_normal_3(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        size_t num_faces = qm.num_faces();

        if(vars._face_normals.size() != num_faces){
            std::cerr << "Wrong number of face normals in the variables to calculate constraints"<<std::endl;
        }

        if(Dest.size() != num_faces){
            Dest.resize(num_faces);
        }

        for(int i = 0; i < num_faces; i++){
            Dest[i] = vars._face_normals.col(i).dot(vars._face_normals.col(i)) - 1.0;
        }
    }

    void get_constraint_ruling_0(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        
        size_t num_faces = quadTopol.getNumFaces();
        
        if(Dest.size() != 3*num_faces){
            Dest.resize(3*num_faces);
        }
        
        for(const auto& heh : qm._mesh.halfedge){
            // first, obtain the indices of the halfedges 
            size_t i = heh.idx();
            auto fh1 = qm._mesh.face_handle(heh);
            auto fh2 = qm._mesh.opposite_face_handle(heh);
            Dest.segment(3*i,3) = vars._rulings.col(i) - vars._normals.col(fh1).cross(vars._normals.col(fh2)); 
        }
    }

    void get_constraint_ruling_1(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        
        size_t num_faces = quadTopol.getNumFaces();

        if(Dest.size() != 3*num_faces){
            Dest.resize(3*num_faces);
        }

        for(const auto& fh : qm._mesh.faces){
            int faceIndex = fh.idx();

            int node_0 = quadTopol.getNodeOfQuad(faceIndex,0);
            int node_1 = quadTopol.getNodeOfQuad(faceIndex,1);
            int node_2 = quadTopol.getNodeOfQuad(faceIndex,2);
            int node_3 = quadTopol.getNodeOfQuad(faceIndex,3);

            int heh_01 = quadTopol.getHalfEdgeOfQuad(faceIndex,0);
            int heh_12 = quadTopol.getHalfEdgeOfQuad(faceIndex,1);
            int heh_23 = quadTopol.getHalfEdgeOfQuad(faceIndex,2);
            int heh_30 = quadTopol.getHalfEdgeOfQuad(faceIndex,3);

            RealType w0 = vars._vertex_weights[node_0];
            RealType w1 = vars._vertex_weights[node_1];
            RealType w2 = vars._vertex_weights[node_2];
            RealType w3 = vars._vertex_weights[node_3];

            Vec3d ruling_01 = vars._rulings.col(heh_01);
            Vec3d ruling_12 = vars._rulings.col(heh_12);
            Vec3d ruling_23 = vars._rulings.col(heh_23);
            Vec3d ruling_30 = vars._rulings.col(heh_30);

            Dest.segment(3*faceIndex,3) = vars._reweighted_rulings_1.col(faceIndex) - (w1 + w2)*ruling_12 + (w3 + w0)*(ruling_30);
        }
    }

    void get_constraint_ruling_2(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        
        size_t num_faces = quadTopol.getNumFaces();

        if(Dest.size() != 3*num_faces){
            Dest.resize(3*num_faces);
        }

        for(const auto& fh : qm._mesh.faces){
            int faceIndex = fh.idx();

            int node_0 = quadTopol.getNodeOfQuad(faceIndex,0);
            int node_1 = quadTopol.getNodeOfQuad(faceIndex,1);
            int node_2 = quadTopol.getNodeOfQuad(faceIndex,2);
            int node_3 = quadTopol.getNodeOfQuad(faceIndex,3);

            int heh_01 = quadTopol.getHalfEdgeOfQuad(faceIndex,0);
            int heh_12 = quadTopol.getHalfEdgeOfQuad(faceIndex,1);
            int heh_23 = quadTopol.getHalfEdgeOfQuad(faceIndex,2);
            int heh_30 = quadTopol.getHalfEdgeOfQuad(faceIndex,3);

            RealType w0 = vars._vertex_weights[node_0];
            RealType w1 = vars._vertex_weights[node_1];
            RealType w2 = vars._vertex_weights[node_2];
            RealType w3 = vars._vertex_weights[node_3];

            Vec3d ruling_01 = vars._rulings.col(heh_01);
            Vec3d ruling_12 = vars._rulings.col(heh_12);
            Vec3d ruling_23 = vars._rulings.col(heh_23);
            Vec3d ruling_30 = vars._rulings.col(heh_30);

            Dest.segment(3*faceIndex,3) = vars._reweighted_rulings_2.col(faceIndex) - (w0 + w1)*ruling_01 + (w2 + w3)*(ruling_23);
        }
    }

    void get_constraint_developable(const QuadMesh &qm, const QuadMeshTopologySaver &quadTopol, const Variables &vars, VectorType &Dest){
        
        size_t num_faces = quadTopol.getNumFaces();

        if(Dest.size() != num_faces){
            Dest.resize(num_faces);
        }

        for(int i = 0; i < quadTopol.getNumFaces(); i++){
            Dest[i] = vars.reweighted_rulings_1.col(i).cross(vars.reweighted_rulings_2.col(i));
        }
    }
};

template <typename ConfiguratorType>
class ConstraintGrad{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;

    /*
    Ordering of the dofs: 
    1. dummy_weights from 0 to num_vertices
    2. vertices from num_vertices to 4*num_vertices
    3. weights from 4*num_vertices to 5*num_vertices
    4. reweighted_vertices from 5*num_vertices to 8*num_vertices
    5. reweighted_edges_1 from 8*num_vertices to (8*num_vertices + 3*num_faces)
    6. reweighted_edges_2 from (8*num_vertices + 3*num_faces) to (8*num_vertices + 6*num_faces)
    7. normals from (8*num_vertices + 6*num_faces) to (8*num_vertices + 9*num_faces)
    8. rulings from (8*num_vertices + 9*num_faces) to (8*num_vertices + 9*num_faces + 2*num_edges)
    9. reweighted_rulings from (8*num_vertices + 9*num_faces + 2*num_edges) to 
    */

    private:
        typedef std::pair<size_t,size_t> Bounds;
        const Bounds _dummy_weights_bounds;
        const Bounds _vertices_bounds;
        const Bounds _weights_bounds;
        const Bounds _reweighted_vertices_bounds;
        const Bounds _reweighted_edges_1_bounds;
        const Bounds _reweighted_edges_2_bounds;
        const Bounds _normals_bounds;
        const Bounds _rulings_bounds;
        const Bounds _reweighted_rulings_1_bounds;
        const Bounds _reweighted_rulings_2_bounds;
        const size_t _numDofs;
        const size_t _num_vertices;
        const size_t _num_edges;
        const size_t _num_halfedges;
        const size_t _num_faces;

    public:
        ConstraintGrad(size_t num_vertices, size_t num_edges, size_t num_faces)
                        :_dummy_weights_bounds(0,num_vertices),
                         _vertices_bounds(num_vertices,4*num_vertices),
                         _weights_bounds(4*num_vertices,5*num_vertices),
                         _reweighted_vertices_bounds(5*num_vertices, 8*num_vertices),
                         _reweighted_edges_1_bounds(8*num_vertices,8*num_vertices + 3*num_faces),
                         _reweighted_edges_2_bounds(8*num_vertices + 3*num_faces, 8*num_vertices + 6*num_faces),
                         _normals_bounds(8*num_vertices + 6*num_faces, 8*num_vertices + 9*num_faces),
                         _rulings_bounds(8*num_vertices + 9*num_faces, 8*num_vertices + 9*num_faces + 2*num_edges),
                         _reweighted_rulings_1_bounds(8*num_vertices + 9*num_faces + 2*num_edges, 8*num_vertices + 12*num_faces + 2*num_edges),
                         _reweighted_rulings_2_bounds(8*num_vertices + 12*num_faces + 2*num_edges, 8*num_vertices + 15*num_faces + 2*num_edges),
                         _numDofs(8*num_vertices + 15*num_faces + 2*num_edges),
                         _num_vertices(num_vertices),
                         _num_edges(num_edges),
                         _num_halfedges(_num_halfedges),
                         _num_faces(num_faces)
                         {}

        void get_constraintgrad_vertex_1(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            // Check that all sizes are correct:
            if( vars._vertices.size() != _num_vertices ){
                std::cerr << "Wrong number of vertices in the variables to calculate constraints"<<std::endl;
            }

            if( vars._vertex_weights.size() != _num_vertices ){
                std::cerr << "Wrong number of vertex weights in the variables to calculate constraints"<<std::endl;
            }

            if( vars._reweighted_vertices.size() != _num_vertices ){
                std::cerr << "Wrong number of reweighted vertices in the variables to calculate constraints"<<std::endl;
            }

            if(Dest.rows() != 3*_num_vertices || Dest.cols() != _numDofs){
                Dest.resize(3*_num_vertices, _numDofs);
            }

            Dest.setZero();

            for(int i = 0; i < 3*_num_vertices;i++){
                int vertexIndex = std::floor(i/3);
                int coordIndex = i%3;
                RealType w_i = vars._vertex_weights[vertexIndex];

                Dest(i, i + _vertices_bounds.first) = - w_i;
                Dest(i, vertexIndex + _weights_bounds.first) = vars._vertices[coordIndex,vertexIndex];
                Dest(i, i + _reweighted_vertices_bounds.first) = 1.0;
            }
        }

        void get_constraintgrad_vertex_2(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            if(Dest.rows()!= _num_vertices || Dest.cols() != _numDofs){
                Dest.resize(_num_vertices, _numDofs);
            }
            
            for(int i = 0; i < _num_vertices; i++){
                RealType dummy_i = vars._dummy_weights[i];
                Dest[i,i] = -2.0*dummy_i;
                Dest[i, 4*_num_vertices + i] = 1.0;
            }
        }

        void get_constraintgrad_edge_vec_1(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            if(Dest.rows() != 3*_num_vertices || Dest.cols() != _numDofs){
                Dest.resize(3, _numDofs);
            }

            Dest.setZero();

            // From average edge vector 1
            Dest.block(0, _reweighted_edges_1_bounds.first, 3*_num_faces, 3*_num_faces).setIdentity();

            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = quadTopol.getNodeOfQuad(i,0);
                int node_1 = quadTopol.getNodeOfQuad(i,1);
                int node_2 = quadTopol.getNodeOfQuad(i,2);
                int node_3 = quadTopol.getNodeOfQuad(i,3);

                RealType w0 = vars._vertex_weights[node_0];
                RealType w1 = vars._vertex_weights[node_1];
                RealType w2 = vars._vertex_weights[node_2];
                RealType w3 = vars._vertex_weights[node_3];

                VectorType v0_tilde = vars._reweighted_vertices[:,node_0];
                VectorType v1_tilde = vars._reweighted_vertices[:,node_1];
                VectorType v2_tilde = vars._reweighted_vertices[:,node_2];
                VectorType v3_tilde = vars._reweighted_vertices[:,node_3];

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                Dest.block(_weights_bounds.first + node_0, 3*i, 3, 1) = -(v2_tilde + v3_tilde);
                Dest.block(_weights_bounds.first + node_1, 3*i, 3, 1) = -(v2_tilde + v3_tilde);
                Dest.block(_weights_bounds.first + node_2, 3*i, 3, 1) = (v1_tilde + v0_tilde);
                Dest.block(_weights_bounds.first + node_3, 3*i, 3, 1) = (v1_tilde + v0_tilde);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                Dest.block(_reweighted_vertices_bounds.first + 3*node_0,3*i,3,1) = (w2+w3)*VectorType::Ones(3);
                Dest.block(_reweighted_vertices_bounds.first + 3*node_1,3*i,3,1) = (w2+w3)*VectorType::Ones(3);
                Dest.block(_reweighted_vertices_bounds.first + 3*node_2,3*i,3,1) = -(w0+w1)*VectorType::Ones(3);
                Dest.block(_reweighted_vertices_bounds.first + 3*node_3,3*i,3,1) = -(w0+w1)*VectorType::Ones(3);
            }

        }

        void get_constraintgrad_edge_vec_2(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            if(Dest.rows() != 3*_num_vertices || Dest.cols() != _numDofs){
                Dest.resize(3, _numDofs);
            }

            Dest.setZero();

            // From average edge vector 1
            Dest.block(0, _reweighted_edges_2_bounds.first, 3*num_faces, 3*num_faces).setIdentity();

            for(int i = 0; i < _num_faces; i++){

                // First, obtain first vertex idces of the face
                int node_0 = quadTopol.getNodeOfQuad(i,0);
                int node_1 = quadTopol.getNodeOfQuad(i,1);
                int node_2 = quadTopol.getNodeOfQuad(i,2);
                int node_3 = quadTopol.getNodeOfQuad(i,3);

                RealType w0 = vars._vertex_weights[node_0];
                RealType w1 = vars._vertex_weights[node_1];
                RealType w2 = vars._vertex_weights[node_2];
                RealType w3 = vars._vertex_weights[node_3];

                VectorType v0_tilde = vars._reweighted_vertices[:,node_0];
                VectorType v1_tilde = vars._reweighted_vertices[:,node_1];
                VectorType v2_tilde = vars._reweighted_vertices[:,node_2];
                VectorType v3_tilde = vars._reweighted_vertices[:,node_3];

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                Dest.block(_weights_bounds.first + node_0, 3*i, 3, 1) = (v2_tilde + v1_tilde);
                Dest.block(_weights_bounds.first + node_1, 3*i, 3, 1) = -(v0_tilde + v3_tilde);
                Dest.block(_weights_bounds.first + node_2, 3*i, 3, 1) = -(v0_tilde + v3_tilde);
                Dest.block(_weights_bounds.first + node_3, 3*i, 3, 1) = (v2_tilde + v1_tilde);

                //Derivatives w.r.t. the vertex dofs of reweighted vertices Tilde{v}_0, ..., Tilde{v}_3
                Dest.block(_reweighted_vertices_bounds.first + 3*node_0,3*i,3,1) = -(w1+w2)*VectorType::Ones(3);
                Dest.block(_reweighted_vertices_bounds.first + 3*node_1,3*i,3,1) = (w0+w3)*VectorType::Ones(3);
                Dest.block(_reweighted_vertices_bounds.first + 3*node_2,3*i,3,1) = (w0+w3)*VectorType::Ones(3);
                Dest.block(_reweighted_vertices_bounds.first + 3*node_3,3*i,3,1) = -(w1+w2)*VectorType::Ones(3);
            }

        }

        void get_constraintgrad_normal_1(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            if(Dest.rows() != _num_faces || Dest.cols() != _numDofs ){
                Dest.resize(_num_faces, _numDofs);
            }

            for(const auto& fh : qm._mesh.faces()){
                size_t i = fh.idx();

                // Set the derivatives w.r.t. the normals
                Dest.block(3*i,_normals_bounds.first + 3*i,3,3) = vars._reweighted_edges_1.col(i).asDiagonal();
                Dest.block(3*i,_reweighted_edges_1_bounds.first + 3*i,3,3) = vars._normals.col(i).asDiagonal();
            }
        }

        void get_constraintgrad_normal_2(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            if(Dest.rows() != _num_faces || Dest.cols() != _numDofs ){
                Dest.resize(_num_faces, _numDofs);
            }

            for(const auto& fh : qm._mesh.faces()){
                size_t i = fh.idx();

                // Set the derivatives w.r.t. the normals
                Dest.block(3*i,_normals_bounds.first + 3*i,3,3) = vars._reweighted_edges_2.col(i).asDiagonal();
                Dest.block(3*i,_reweighted_edges_2_bounds.first + 3*i,3,3) = vars._normals.col(i).asDiagonal();
            }
        }

        void get_constraintgrad_normal_3(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            if(Dest.rows() != _num_faces || Dest.cols() != _numDofs ){
                Dest.resize(_num_faces, _numDofs);
            }

            for(const auto& fh : qm._mesh.faces()){
                size_t i = fh.idx();

                // Set the derivatives w.r.t. the normals
                Dest.block(i,_normals_bounds.first + i,1,3) = 2.0*vars._normals.col(i).transpose();
            }
        }
        
        void get_constraintgrad_ruling_0(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            if(Dest.rows() != 3*_num_halfedges || Dest.cols() != _numDofs){
                Dest.resize(3*_num_halfedges, _numDofs);
            }

            Dest.setZero();

            for(const auto& heh : qm._mesh.halfedges){
                size_t i = heh.idx();
                size_t faceIdx_1 = qm._mesh.face_handle(heh);
                size_t faceIdx_2 = qm._mesh.opposite_face_handle(heh);
                VectorType normal_1 = vars._normals.col(faceIdx_1);
                VectorType normal_2 = vars._normals.col(faceIdx_2);

                // First, set the derivative w.r.t. the ruling vector
                Dest.block(3*i, _rulings_bounds.first + 3*i, 3, 3) = FullMatrixType::Identity(3,3);

                // Now, set the derivatives w.r.t. the cross product
                // First, differentiate w.r.t. normal_1
                Eigen::Vector3d e;
                e.setZero();
                e[0] = 1.0; 
                Dest.block(3*i, _normals_bounds.first + 3*faceIdx_1, 3, 1) = -e.cross(normal_2);
                e[0] = 0.0;
                e[1] = 1.0;
                Dest.block(3*i, _normals_bounds.first + 3*faceIdx_1 + 1, 3, 1) = -e.cross(normal_2);
                e[1] = 0.0;
                e[2] = 1.0;
                Dest.block(3*i, _normals_bounds.first + 3*faceIdx_1 + 2, 3, 1) = -e.cross(normal_2);

                // Now, differentiate w.r.t. normal_2
                e[2] = 0.0;
                e[0] = 1.0;
                Dest.block(3*i, _normals_bounds.first + 3*faceIdx_2, 3, 1) = -normal_1.cross(e);
                e[0] = 0.0;
                e[1] = 1.0;
                Dest.block(3*i, _normals_bounds.first + 3*faceIdx_2 + 1, 3, 1) = -normal_1.cross(e);
                e[1] = 0.0;
                e[2] = 1.0;
                Dest.block(3*i, _normals_bounds.first + 3*faceIdx_2 + 2, 3, 1) = -normal_1.cross(e);
            }

        }

         void get_constraintgrad_ruling_1(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            if(Dest.rows() != 3*_num_faces || Dest.cols() != _numDofs){
                Dest.resize(3*_num_faces, _numDofs);
            }

            Dest.setZero();

            for(const auto& fh: qm._mesh.faces())
            {
                int faceIndex = fh.idx();

                int node_0 = quadTopol.getNodeOfQuad(faceIndex,0);
                int node_1 = quadTopol.getNodeOfQuad(faceIndex,1);
                int node_2 = quadTopol.getNodeOfQuad(faceIndex,2);
                int node_3 = quadTopol.getNodeOfQuad(faceIndex,3);

                int heh_01 = quadTopol.getHalfEdgeOfQuad(faceIndex,0);
                int heh_12 = quadTopol.getHalfEdgeOfQuad(faceIndex,1);
                int heh_23 = quadTopol.getHalfEdgeOfQuad(faceIndex,2);
                int heh_30 = quadTopol.getHalfEdgeOfQuad(faceIndex,3);

                RealType w0 = vars._vertex_weights[node_0];
                RealType w1 = vars._vertex_weights[node_1];
                RealType w2 = vars._vertex_weights[node_2];
                RealType w3 = vars._vertex_weights[node_3];

                Vec3d ruling_01 = vars._rulings.col(heh_01);
                Vec3d ruling_12 = vars._rulings.col(heh_12);
                Vec3d ruling_23 = vars._rulings.col(heh_23);
                Vec3d ruling_30 = vars._rulings.col(heh_30);

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}
                Dest.block(3*faceIndex,_reweighted_rulings_1_bounds.first + 3*faceIndex,3,3) = FullMatrixType::Identity(3,3);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                Dest.block(3*faceIndex,_weights_bounds.first + node_0, 3, 1) = (ruling_30);
                Dest.block(3*faceIndex,_weights_bounds.first + node_1, 3, 1) = -(ruling_12);
                Dest.block(3*faceIndex,_weights_bounds.first + node_2, 3, 1) = -(ruling_12);
                Dest.block(3*faceIndex,_weights_bounds.first + node_3, 3, 1) = (ruling_30);

                //Derivatives w.r.t. the rulings
                Dest.block(_rulings_bounds.first + 3*heh_12, 3*faceIndex, 3, 3) = -(w1 + w2)*FullMatrixType::Identity(3);
                Dest.block(_rulings_bounds.first + 3*heh_30, 3*faceIndex, 3, 3) = (w3 + w0)*FullMatrixType::Identity(3);
            }
        }

        void get_constraintgrad_ruling_2(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            if(Dest.rows() != 3*_num_faces || Dest.cols() != _numDofs){
                Dest.resize(3*_num_faces, _numDofs);
            }

            Dest.setZero();

            for(const auto& fh: qm._mesh.faces())
            {
                int faceIndex = fh.idx();

                int node_0 = quadTopol.getNodeOfQuad(faceIndex,0);
                int node_1 = quadTopol.getNodeOfQuad(faceIndex,1);
                int node_2 = quadTopol.getNodeOfQuad(faceIndex,2);
                int node_3 = quadTopol.getNodeOfQuad(faceIndex,3);

                int heh_01 = quadTopol.getHalfEdgeOfQuad(faceIndex,0);
                int heh_12 = quadTopol.getHalfEdgeOfQuad(faceIndex,1);
                int heh_23 = quadTopol.getHalfEdgeOfQuad(faceIndex,2);
                int heh_30 = quadTopol.getHalfEdgeOfQuad(faceIndex,3);

                RealType w0 = vars._vertex_weights[node_0];
                RealType w1 = vars._vertex_weights[node_1];
                RealType w2 = vars._vertex_weights[node_2];
                RealType w3 = vars._vertex_weights[node_3];

                Vec3d ruling_01 = vars._rulings.col(heh_01);
                Vec3d ruling_12 = vars._rulings.col(heh_12);
                Vec3d ruling_23 = vars._rulings.col(heh_23);
                Vec3d ruling_30 = vars._rulings.col(heh_30);

                // Differentiate w.r.t. the reweighted ruling r_{1,3,f}
                Dest.block(3*faceIndex,_reweighted_rulings_2_bounds.first + 3*faceIndex,3,3) = FullMatrixType::Identity(3,3);

                // Derivatives w.r.t. the weights w0, w1, w2, w3
                Dest.block(3*faceIndex,_weights_bounds.first + node_0, 3, 1) = -(ruling_01);
                Dest.block(3*faceIndex,_weights_bounds.first + node_1, 3, 1) = -(ruling_01);
                Dest.block(3*faceIndex,_weights_bounds.first + node_2, 3, 1) = (ruling_23);
                Dest.block(3*faceIndex,_weights_bounds.first + node_3, 3, 1) = (ruling_23);

                //Derivatives w.r.t. the rulings
                Dest.block(_rulings_bounds.first + 3*heh_01, 3*faceIndex, 3, 3) = -(w0 + w1)*FullMatrixType::Identity(3);
                Dest.block(_rulings_bounds.first + 3*heh_23, 3*faceIndex, 3, 3) = (w2 + w0)*FullMatrixType::Identity(3);
            }
        }

        void get_constraintgrad_dev(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

            if(Dest.rows() != _num_faces || Dest.cols() != _numDofs){
                Dest.resize(_num_faces, _numDofs);
            }

            Dest.setZero();

            for(const auto& fh: qm._mesh.faces()){

                int i = fh.idx();

                // Derivatives w.r.t. the first ruling vector \Tilde{r}_{13,f}
                Eigen::Vector3d e;
                e.setZero();
                e[0] = 1.0; 
                Dest.block(i,_reweighted_rulings_1_bounds.first,3,1) = e.cross(vars._reweighted_rulings_2);
                e[0] = 0.0;
                e[1] = 1.0;
                Dest.block(i,_reweighted_rulings_1_bounds.first + 1,3,1) = e.cross(vars._reweighted_rulings_2);
                e[1] = 0.0;
                e[2] = 1.0;
                Dest.block(i,_reweighted_rulings_1_bounds.first + 2,3,1) = e.cross(vars._reweighted_rulings_2);
            
                // Derivatives w.r.t. the second ruling vector \Tilde{r}_{02,f}
                e[2] = 0.0;
                e[0] = 1.0;
                Dest.block(i, _reweighted_rulings_2_bounds.first,3,1) = vars._reweighted_rulings_1.cross(e);
                e[0] = 0.0;
                e[1] = 1.0;
                Dest.block(i, _reweighted_rulings_2_bounds.first + 1,3,1) = vars._reweighted_rulings_1.cross(e);
                e[1] = 0.0;
                e[2] = 1.0;
                Dest.block(i, _reweighted_rulings_2_bounds.first + 2,3,1) = vars._reweighted_rulings_1.cross(e);
            }
        }
};