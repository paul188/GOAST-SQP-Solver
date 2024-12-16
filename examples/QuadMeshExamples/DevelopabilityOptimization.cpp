#include <goast/Core.h>
#include <goast/Quads.h>

template<typename ConfiguratorType>
struct Variables{
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;

    VectorType vertices;
    VectorType vertex_weights;
    VectorType reweighted_vertices;
    VectorType dummy_weights;
    std::vector<Vec3d> vertex_normals;
    std::vector<Vec3d> reweighted_edges;
    std::vector<Vec3d> reweighted_rulings;
};

class Optimizer{
    public:
        // need initial geometry and boundary mask
        void solve(QuadMesh &mesh, const std::vector<int> &bdryMask){
            
        }
};