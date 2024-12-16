#include <goast/Quads.h>
#include <goast/QuadMesh/QuadTopology.h>

template<typename ConfiguratorType>
struct Variables{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    VectorType _vertices;
    VectorType _vertex_weights;
    VectorType _reweighted_vertices;
    VectorType _dummy_weights;
    VectorType _rulings;
    VectorType _face_normals;
    VectorType _reweighted_edges;
    VectorType _reweighted_rulings;

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
        _rulings.resize(3,nEdges);
        _face_normals.resize(3,nFaces);
        _reweighted_edges.resize(3,2*nFaces);
        _reweighted_rulings.resize(3,2*nFaces);

        getGeometryMat<FullMatrixType>(qm._mesh, _vertices);

        _vertex_weights.setConstant(2.0);
        _dummy_weights.setConstant(1.0);

        _reweighted_vertices = _vertices.cwiseProduct(_vertex_weights);

        _rulings.setZero();
        _face_normals.setZero();
        _reweighted_edges.setZero();
        _reweighted_rulings.setZero();
    }
};

template <typename ConfiguratorType>
class Constraint{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    public:
        get_constraint(const QuadMesh &qm, const Variables &vars, VectorType &Dest){
            size_t num_vertices = qm.num_vertices();
            size_t num_edges = qm.num_edges();
            size_t num_faces = qm.num_faces();

            // Check that all sizes are correct:
            if( vars.vertices.size() != num_vertices ){
                std::cerr << "Wrong number of vertices in the variables to calculate constraints"<<std::endl;
            }

            if( vars.vertex_weights.size() != num_vertices ){
                std::cerr << "Wrong number of vertex weights in the variables to calculate constraints"<<std::endl;
            }

            if( vars.reweighted_vertices.size() != num_vertices ){
                std::cerr << "Wrong number of reweighted vertices in the variables to calculate constraints"<<std::endl;
            }

            if(vars.dummy_weights.size() != num_vertices ){
                std::cerr << "Wrong number of dummy weights in the variables to calculate constraints"<<std::endl;
            }

            if(vars.face_normals.size() != num_faces ){
                std::cerr << "Wrong number of face normals in the variables to calculate constraints"<<std::endl;
            }



            //Dest.resize(2*qm.num_vertices() + )
        }

};

template <typename ConfiguratorType>
class ConstraintGrad{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;

    public:
        get_constraint_grad(const QuadMesh &qm, const Variables &vars, MatrixType &Dest){

        }

};