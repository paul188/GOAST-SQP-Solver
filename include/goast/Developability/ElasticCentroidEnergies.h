#include <goast/DiscreteShells.h>

class ElasticDevelopabilityCentroidEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>
{

private:
    typedef typename DefaultConfigurator::RealType   RealType;
    typedef typename DefaultConfigurator::VectorType VectorType;
    typedef typename DefaultConfigurator::VecType    VecType;
    typedef typename DefaultConfigurator::MatType    MatType;

    NonlinearMembraneEnergy<ConfiguratorType> _E_mem;
    SimpleBendingEnergy<ConfiguratorType> _E_bend;

    const QuadMeshTopologySaver& _quadTopol;
    const MeshTopologySaver &_centroidTopol;

    const VectorType &_centroidBaseGeometry;

    

    const VectorType _factors_elasticity_dev;
    const VectorType _factors_mem_bend;

public:
    ElasticDevelopabilityCentroidEnergy(const QuadMeshTopologySaver &quadTopol,
                                        const MeshTopologySaver &centroidTopol,
                                        const VectorType &centroidBaseGeometry,
                                        VectorType &factors_elasticity_dev,
                                        VectorType &factors_mem_bend,,
                                        const std::vector<int> &centroidVerticesIndices
                                    )
        : _quadTopol(quadTopol), 
          _centroidTopol(centroidTopol),
          _centroidBaseGeometry(centroidBaseGeometry),
          _factors_elasticity_dev(factors_elasticity_dev),
          _factors_mem_bend(factors_mem_bend),
          _E_bend(centroidTopol, centroidBaseGeometry),
          _E_mem(centroidTopol, centroidBaseGeometry)
        {
            
        }
}