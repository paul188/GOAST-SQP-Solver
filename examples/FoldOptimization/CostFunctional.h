#ifndef COSTFUNCTIONAL_H
#define COSTFUNCTIONAL_H

#include <goast/Core.h>

template<typename ConfiguratorType>
class CostFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>{
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;

        const std::vector<int> &_foldVertices;
    
    public:
        CostFunctional(const std::vector<int> &foldVertices) : _foldVertices(foldVertices) {}

        void apply(const VectorType &activeGeometry, RealType &dest) const{
            for(int i = 0; i < _foldVertices.size(); i++){
                VecType coords;
                getXYZCoord<VectorType, VecType>(activeGeometry, coords, _foldVertices[i]);
                dest -= coords[2];
            }
        }
};

template<typename ConfiguratorType>
class CostFunctionalGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;

        const std::vector<int> &_foldVertices;
        const MeshTopologySaver &_topology;
    
    public:
        CostFunctionalGradient(const MeshTopologySaver& topology, const std::vector<int> &foldVertices) 
        : _foldVertices(foldVertices), _topology(topology) {}

        void apply(const VectorType &activeGeometry, VectorType &dest) const{

        // First component, i.e. derivative w.r.t. t is zero.
           if(dest.size() != activeGeometry.size()+1){
                dest.resize(activeGeometry.size()+1);
           }
           dest.setZero();
           for(int i = 0; i < _foldVertices.size(); i++){
                dest[2*_topology.getNumVertices() + _foldVertices[i]+1] = -1;
           }
        }
};

#endif //COSTFUNCTIONAL_H
