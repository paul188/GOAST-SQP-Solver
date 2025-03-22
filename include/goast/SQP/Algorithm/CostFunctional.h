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

//NOT WORKING YET

// Modified cost functional, specially for the ArcFoldOptimization problem
// added penalty that ensures flexible Dirichlet conditions

template<typename ConfiguratorType>
class PenaltyCostFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>{
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;

        const std::vector<int> &_vertices;

        std::pair<std::vector<int>, std::vector<int>> &_penalty_vertices_1;
        std::pair<std::vector<int>, std::vector<int>> &_penalty_vertices_2;
        int _num_penalties;
    
    public:
        PenaltyCostFunctional(const std::vector<int> &vertices,
                              std::pair<std::vector<int>, std::vector<int>> &penalty_vertices_1,
                              std::pair<std::vector<int>, std::vector<int>> &penalty_vertices_2) : 
                              _vertices(vertices),
                              _penalty_vertices_1(penalty_vertices_1),
                              _penalty_vertices_2(penalty_vertices_2) {
            assert(_penalty_vertices_1.first.size() == _penalty_vertices_1.second.size());
            assert(_penalty_vertices_2.first.size() == _penalty_vertices_2.second.size());
            assert(_penalty_vertices_1.first.size() == _penalty_vertices_2.first.size());
            _num_penalties = _penalty_vertices_1.first.size();
        }

        void apply(const VectorType &activeGeometry, RealType &dest) const{
            dest = 0;
            VectorType _penalty_vec_1(_num_penalties);
            VectorType _penalty_vec_2(_num_penalties);

            for(int i = 0; i < _vertices.size(); i++){
                VecType coords;
                getXYZCoord<VectorType, VecType>(activeGeometry, coords, _vertices[i]);
                dest += coords[2];
            }

            for(int i = 0; i < _num_penalties; i++){
                VecType coords_1_1, coords_1_2, coords_2_1, coords_2_2;
                getXYZCoord<VectorType, VecType>(activeGeometry, coords_1_1, _penalty_vertices_1.first[i]);
                getXYZCoord<VectorType, VecType>(activeGeometry, coords_1_2, _penalty_vertices_1.second[i]);
                getXYZCoord<VectorType, VecType>(activeGeometry, coords_2_1, _penalty_vertices_2.first[i]);
                getXYZCoord<VectorType, VecType>(activeGeometry, coords_2_2, _penalty_vertices_1.first[i]);
                _penalty_vec_1[i] = coords_1_1[2] - coords_1_2[2];
                _penalty_vec_2[i] = coords_2_1[2] - coords_2_2[2];
            }

            RealType penalty = _penalty_vec_1.dot(_penalty_vec_2);
            if(penalty < 0)
            {
                dest -= _penalty_vec_1.dot(_penalty_vec_2);
            }
        }
};

template<typename ConfiguratorType>
class PenaltyCostFunctionalGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::VecType VecType;

        const std::vector<int> &_vertices;
        const MeshTopologySaver &_topology;

        std::pair<std::vector<int>, std::vector<int>> &_penalty_vertices_1;
        std::pair<std::vector<int>, std::vector<int>> &_penalty_vertices_2;

        int _num_penalties;
    
    public:
        PenaltyCostFunctionalGradient(const MeshTopologySaver& topology, 
                                      const std::vector<int> &vertices,
                                      std::pair<std::vector<int>, std::vector<int>> &penalty_vertices_1,
                                      std::pair<std::vector<int>, std::vector<int>> &penalty_vertices_2) : 
                                      _vertices(vertices),
                                      _topology(topology),
                                      _penalty_vertices_1(penalty_vertices_1),
                                      _penalty_vertices_2(penalty_vertices_2) {
                    assert(_penalty_vertices_1.first.size() == _penalty_vertices_1.second.size());
                    assert(_penalty_vertices_2.first.size() == _penalty_vertices_2.second.size());
                    assert(_penalty_vertices_1.first.size() == _penalty_vertices_2.first.size());
                    _num_penalties = _penalty_vertices_1.first.size();
        }

        void apply(const VectorType &activeGeometry, VectorType &dest) const{
            VectorType _penalty_vec_1(_num_penalties);
            VectorType _penalty_vec_2(_num_penalties);

            // First component, i.e. derivative w.r.t. t is zero.
           if(dest.size() != activeGeometry.size()){
                dest.resize(activeGeometry.size());
           }
           dest.setZero();
           for(int i = 0; i < _vertices.size(); i++){
                dest[2*_topology.getNumVertices() + _vertices[i]] = 1;
           }

           for(int i = 0; i < _num_penalties; i++){
                VecType coords_1_1, coords_1_2, coords_2_1, coords_2_2;
                getXYZCoord<VectorType, VecType>(activeGeometry, coords_1_1, _penalty_vertices_1.first[i]);
                getXYZCoord<VectorType, VecType>(activeGeometry, coords_1_2, _penalty_vertices_1.second[i]);
                getXYZCoord<VectorType, VecType>(activeGeometry, coords_2_1, _penalty_vertices_2.first[i]);
                getXYZCoord<VectorType, VecType>(activeGeometry, coords_2_2, _penalty_vertices_2.second[i]);
                _penalty_vec_1[i] = coords_1_1[2] - coords_1_2[2];
                _penalty_vec_2[i] = coords_2_1[2] - coords_2_2[2];
            }

            RealType penalty = _penalty_vec_1.dot(_penalty_vec_2);
            if(penalty < 0)
            {
                for(int i = 0; i < _num_penalties; i++)
                {
                    dest[2*_topology.getNumVertices() + _penalty_vertices_1.first[i]] -= _penalty_vec_2[i];
                    dest[2*_topology.getNumVertices() + _penalty_vertices_1.second[i]] += _penalty_vec_2[i];
                    dest[2*_topology.getNumVertices() + _penalty_vertices_2.first[i]] -= _penalty_vec_1[i];
                    dest[2*_topology.getNumVertices() + _penalty_vertices_2.second[i]] += _penalty_vec_1[i];
                }
            }
        }
};