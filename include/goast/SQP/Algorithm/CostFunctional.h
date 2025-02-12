#pragma once

#include <goast/Core.h>

/*
 \brief CostFunctional class
 
 This class is used to define the cost functional for the optimization problem.
 The cost functional is the sum of the z-components of the top of the _vertices
 
 \author Johannssen
 
*/

template<typename ConfiguratorType>
class CostFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>{
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;

        const std::vector<int> &_vertices;
    
    public:
        CostFunctional(const std::vector<int> &vertices) : _vertices(vertices) {}

        void apply(const VectorType &activeGeometry, RealType &dest) const{
            dest = 0;
            for(int i = 0; i < _vertices.size(); i++){
                VecType coords;
                getXYZCoord<VectorType, VecType>(activeGeometry, coords, _vertices[i]);
                dest += coords[2];
            }
        }
};

template<typename ConfiguratorType>
class CostFunctionalGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;

        const std::vector<int> &_vertices;
        const MeshTopologySaver &_topology;
    
    public:
        CostFunctionalGradient(const MeshTopologySaver& topology, const std::vector<int> &vertices) : _vertices(vertices), _topology(topology){}

        void apply(const VectorType &activeGeometry, VectorType &dest) const{

        // First component, i.e. derivative w.r.t. t is zero.
           if(dest.size() != activeGeometry.size()){
                dest.resize(activeGeometry.size());
           }
           dest.setZero();
           for(int i = 0; i < _vertices.size(); i++){
                dest[2*_topology.getNumVertices() + _vertices[i]] = 1;
           }
        }
};
