#pragma once

#include <goast/Core.h>

#include <goast/Developability/Constraints.h>
#include <goast/Developability/BoundaryDOFS_quad.h>
#include <goast/DiscreteShells.h>
#include <fstream>
#include <iostream>
#include <memory>


template<typename ConfiguratorType>
class CoordConverter
{
    protected:
        typedef typename DefaultConfigurator::RealType   RealType;
        typedef typename DefaultConfigurator::VectorType VectorType;
        typedef typename DefaultConfigurator::VecType    VecType;
        typedef typename DefaultConfigurator::MatType    MatType;

        const QuadMeshTopologySaver &_quadTopol;
        const MeshTopologySaver &_centroidTopol;

        std::map<int,int> _centroidToQuadIdxMap;
        std::vector<int> _noQuadCentroidVertices;

        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> _perm;
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> _perm2;

        const size_t _num_centroid_vertices;
        const size_t _num_quad_vertices;
        const size_t _num_noquad_vertices;


    public:
        CoordConverter(const QuadMeshTopologySaver &quadTopol,
                       const MeshTopologySaver &centroidTopol,
                       const VectorType &quadBaseGeometry,
                       const VectorType &centroidBaseGeometry)
            : _quadTopol(quadTopol), 
              _centroidTopol(centroidTopol),
              _num_centroid_vertices(centroidTopol.getNumVertices()),
              _num_quad_vertices(quadTopol.getNumVertices()),
              _num_noquad_vertices(centroidTopol.getNumVertices() - quadTopol.getNumVertices())
        {
            // First, assemble centroidToQuadIdxMap and noQuadCentroidVertices
            bool foundMatch = false;
            for(int i = 0; i < _num_centroid_vertices; i++){
                VecType coords_tri;
                getXYZCoord<VectorType, VecType>(centroidBaseGeometry,coords_tri, i);
                for(int j = 0; j < _num_quad_vertices; j++)
                {
                    VecType coords_quad;
                    coords_quad[0] = quadBaseGeometry[3*j];
                    coords_quad[1] = quadBaseGeometry[3*j+1];
                    coords_quad[2] = quadBaseGeometry[3*j+2];
    
                    if((coords_quad - coords_tri).norm() < 1e-4)
                    {
                        _centroidToQuadIdxMap[i] = j;
                        foundMatch = true;
                        break;
                    }
                    
                }
                if(!foundMatch)
                {
                    _noQuadCentroidVertices.push_back(i);
                    _centroidToQuadIdxMap[i] = -1;
                }
                
                // reset foundMatch
                foundMatch = false;
            }

            Eigen::VectorXi permVec(3*_num_centroid_vertices);

            for(int i = 0; i  < _num_centroid_vertices; ++i)
            {
                permVec[i] = i + 2*i;
                permVec[i + _num_centroid_vertices] = i + 2*i + 1;
                permVec[i + 2*_num_centroid_vertices] = i + 2*i + 2;
            }

            _perm.indices() = permVec;

            // Now, there is still the permutation matrix to convert from centroid topology to quadTopology
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> p2;
            
            Eigen::VectorXi permVec2(3*_num_centroid_vertices);

            // Put the noQuadCentroid vertices at the beginning. So that we can later segment only the last part
            for(int i = 0; i < _num_noquad_vertices; ++i)
            {
                permVec2[3*_noQuadCentroidVertices[i]] = 3*i;
                permVec2[3*_noQuadCentroidVertices[i] + 1] = 3*i + 1;
                permVec2[3*_noQuadCentroidVertices[i] + 2] = 3*i +2;
            }

            int counterQuad = 0;
            for(int i = 0; i < _num_centroid_vertices; i++){
                if(_centroidToQuadIdxMap[i] == -1)
                {
                    // Skip this
                    continue;
                }
                else
                {
                    permVec2[3*i] = 3*_num_noquad_vertices + 3*counterQuad;
                    permVec2[3*i + 1] = 3*_num_noquad_vertices + 3*counterQuad + 1;
                    permVec2[3*i + 2] = 3*_num_noquad_vertices + 3*counterQuad + 2;
                    counterQuad++;
                }
            }

            _perm2.indices() = permVec2;

        }

        void convertCentroidToQuad(const VectorType &centroidGeom, VectorType &quadGeom) const
        {
            if(quadGeom.size() != 3*_num_quad_vertices)
            {
                quadGeom.resize(3*_num_quad_vertices);
            }

            VectorType reorderedCentroidGeom = _perm2*_perm*centroidGeom;

            VectorType quadVertices = reorderedCentroidGeom.segment(3*_num_noquad_vertices, 3*_num_quad_vertices);

            // Now, map the indices to the actual indices corresponding to the quad topology
            int noQuadCounter = 0;
            for(int i = 0; i < _num_centroid_vertices; i++)
            {
                if(noQuadCounter < _num_noquad_vertices){
                    if(_centroidToQuadIdxMap.at(i) == -1)
                    {
                        noQuadCounter++;
                        continue;
                    }
                }
                quadGeom[3*_centroidToQuadIdxMap.at(i)] = quadVertices[3*(i - noQuadCounter)];
                quadGeom[3*_centroidToQuadIdxMap.at(i) + 1] = quadVertices[3*(i - noQuadCounter) + 1];
                quadGeom[3*_centroidToQuadIdxMap.at(i) + 2] = quadVertices[3*(i - noQuadCounter) + 2];
            }

        }

        // This method fills the missing points with zeros -> useful for gradient calculations
        void convertQuadToCentroid(const VectorType &quadGeom, VectorType &centroidGeom) const
        {
            if(centroidGeom.size() != 3*_num_centroid_vertices)
            {
                centroidGeom.resize(3*_num_centroid_vertices);
            }

            VectorType reorderedQuadGeom(3*_num_quad_vertices);

            int noQuadCounter = 0;
            for(int i = 0; i < _num_centroid_vertices; i++)
            {
                if(noQuadCounter < _num_noquad_vertices){
                    if(_centroidToQuadIdxMap.at(i) == -1)
                    {
                        noQuadCounter++;
                        continue;
                    }
                }
                reorderedQuadGeom[3*(i - noQuadCounter)] = quadGeom[3*_centroidToQuadIdxMap.at(i)];
                reorderedQuadGeom[3*(i - noQuadCounter) + 1] = quadGeom[3*_centroidToQuadIdxMap.at(i) + 1];
                reorderedQuadGeom[3*(i - noQuadCounter) + 2] = quadGeom[3*_centroidToQuadIdxMap.at(i) + 2];
            }

            centroidGeom.segment(0, 3*_num_noquad_vertices) = VectorType::Zero(3*_num_noquad_vertices);
            centroidGeom.segment(3*_num_noquad_vertices, centroidGeom.size() - 3*_num_noquad_vertices) = reorderedQuadGeom;

            centroidGeom = (_perm2*_perm).inverse()*centroidGeom;
        }

        void convertQuadToCentroid2(const VectorType &quadGeom, VectorType &centroidGeom, TriMesh& centroid_mesh) const
        {
            if(centroidGeom.size() != 3*_num_centroid_vertices)
            {
                centroidGeom.resize(3*_num_centroid_vertices);
            }


            VectorType reorderedQuadGeom(3*_num_quad_vertices);
            int noQuadCounter = 0;
            for(int i = 0; i < _num_centroid_vertices; i++)
            {
                if(noQuadCounter < _num_noquad_vertices){
                    if(_centroidToQuadIdxMap.at(i) == -1)
                    {
                        noQuadCounter++;
                        continue;
                    }
                }
                reorderedQuadGeom[3*(i - noQuadCounter)] = quadGeom[3*_centroidToQuadIdxMap.at(i)];
                reorderedQuadGeom[3*(i - noQuadCounter) + 1] = quadGeom[3*_centroidToQuadIdxMap.at(i) + 1];
                reorderedQuadGeom[3*(i - noQuadCounter) + 2] = quadGeom[3*_centroidToQuadIdxMap.at(i) + 2];
            }

            centroidGeom.segment(0, 3*_num_noquad_vertices) = VectorType::Zero(3*_num_noquad_vertices);
            centroidGeom.segment(3*_num_noquad_vertices, centroidGeom.size() - 3*_num_noquad_vertices) = reorderedQuadGeom;

            centroidGeom = (_perm2*_perm).inverse()*centroidGeom;

            // Now, fill the zero centroid entries
            for(int i = 0; i < _noQuadCentroidVertices.size(); i++)
            {
                int idx = _noQuadCentroidVertices[i];
                VecType coords;
                coords[0] = 0.0;
                coords[1] = 0.0;
                coords[2] = 0.0;
                VertexHandle vh = centroid_mesh.vertex_handle(idx);
                for (TriMesh::VertexVertexIter vv_it = centroid_mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
                    VecType coords_neighbour;
                    getXYZCoord<VectorType, VecType>(centroidGeom,coords_neighbour, vv_it->idx());
                    coords += 0.25*coords_neighbour;
                }
                setXYZCoord<VectorType, VecType>(centroidGeom, coords, idx);
            }
        }

        std::vector<int> getNoQuadCentroidVertices() const
        {
            return _noQuadCentroidVertices;
        }
};

template<typename ConfiguratorType>
class BijectiveCoordConverter
{
    protected:
        typedef typename DefaultConfigurator::RealType   RealType;
        typedef typename DefaultConfigurator::VectorType VectorType;
        typedef typename DefaultConfigurator::VecType    VecType;
        typedef typename DefaultConfigurator::MatType    MatType;

        const QuadMeshTopologySaver &_quadTopol;
        const MeshTopologySaver &_centroidTopol;

        std::map<int,int> _centroidToQuadIdxMap;
        std::map<int,int> _quadToCentroidIdxMap;

        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> _perm;

        const VectorType &_quadBaseGeometry;
        const VectorType &_centroidBaseGeometry;

    public:
        BijectiveCoordConverter(QuadMeshTopologySaver &quadTopol,
                        MeshTopologySaver &centroidTopol,
                       const VectorType &quadBaseGeometry,
                       const VectorType &centroidBaseGeometry)
            : _quadTopol(quadTopol), 
              _centroidTopol(centroidTopol),
              _quadBaseGeometry(quadBaseGeometry),
              _centroidBaseGeometry(centroidBaseGeometry)
        {
            // First, assemble centroidToQuadIdxMap
            bool foundMatch = false;
            for(int i = 0; i < _centroidTopol.getNumVertices(); i++){
                VecType coords_tri;
                getXYZCoord<VectorType, VecType>(_centroidBaseGeometry,coords_tri, i);
                for(int j = 0; j < _quadTopol.getNumVertices(); j++)
                {
                    VecType coords_quad;
                    coords_quad[0] = _quadBaseGeometry[3*j];
                    coords_quad[1] = _quadBaseGeometry[3*j+1];
                    coords_quad[2] = _quadBaseGeometry[3*j+2];
    
                    if((coords_quad - coords_tri).norm() < 1e-4)
                    {
                        _centroidToQuadIdxMap[i] = j;
                        _quadToCentroidIdxMap[j] = i;
                        foundMatch = true;
                        break;
                    }
                    
                }
                if(!foundMatch)
                {
                    throw std::runtime_error("Error: No bijective mapping between centroid and quad mesh possible, as there are centroid vertices that do not correspond to any quad vertex.");
                }
                
                // reset foundMatch
                foundMatch = false;
            }

            Eigen::VectorXi permVec(3*_centroidTopol.getNumVertices());

            for(int i = 0; i  < _quadTopol.getNumVertices(); ++i)
            {
                permVec[i] = _quadToCentroidIdxMap[i];
                permVec[i + 1] = _quadToCentroidIdxMap[i] + _centroidTopol.getNumVertices();
                permVec[i + 2] = _quadToCentroidIdxMap[i] + 2*_centroidTopol.getNumVertices();
            }

            _perm.indices() = permVec;
        }

        void convertQuadToCentroid(const VectorType &quadGeom, VectorType &centroidGeom) const
        {
            
        }
};

template<typename ConfiguratorType>
class QuadTriConverter
{
    protected:
        int _num_vertices;
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> _perm;

    public:
        QuadTriConverter(int num_vertices) 
        :_num_vertices(num_vertices)
        {
            Eigen::VectorXi permVec(3*_num_vertices);
            for(int i = 0; i < _num_vertices; i++)
            {
                permVec[3*i] = i;
                permVec[3*i + 1] = i + _num_vertices;
                permVec[3*i + 2] = i + 2*_num_vertices;
            }

            _perm.indices() = permVec;
        }

        void convertTriToQuad(const VectorType &triGeom, VectorType &quadGeom) const
        {
            if(triGeom.size() != quadGeom.size())
            {
                quadGeom.resize(triGeom.size());
            }
            quadGeom = _perm.inverse()*triGeom;
        }

        void convertQuadToTri(const VectorType &quadGeom, VectorType &triGeom) const
        {
            if(quadGeom.size() != triGeom.size())
            {
                triGeom.resize(quadGeom.size());
            }
            triGeom = _perm*quadGeom;
        }

        void convertTriHessianToQuad(const typename ConfiguratorType::SparseMatrixType &triHessian, typename ConfiguratorType::SparseMatrixType &quadHessian) const
        {
            if(triHessian.rows() != quadHessian.rows() || triHessian.cols() != quadHessian.cols())
            {
                quadHessian.resize(triHessian.rows(), triHessian.cols());
            }
            quadHessian = _perm.transpose()*triHessian*_perm;
        }

};

// Combines the developability constraint with the centroid bending and membrane energies
template<typename ConfiguratorType>
class ElasticDevelopabilityCentroidEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> 
{

protected:
    typedef typename DefaultConfigurator::RealType   RealType;
    typedef typename DefaultConfigurator::VectorType VectorType;
    typedef typename DefaultConfigurator::VecType    VecType;
    typedef typename DefaultConfigurator::MatType    MatType;

    NonlinearMembraneEnergy<ConfiguratorType> _E_mem;
    SimpleBendingEnergy<ConfiguratorType> _E_bend;

    const QuadMeshTopologySaver& _quadTopol;

    const MatrixType _W;

    const Constraint<ConfiguratorType> &_constraint;

    const VectorType _factors_elasticity_dev;
    const VectorType _factors_mem_bend;

    const VarsIdx<ConfiguratorType> _variablesIdx;
    const QuadTriConverter<ConfiguratorType> _quadTriConverter;

    mutable int counter;

    int _num_vertices;

    std::string _path;
    std::string _path_vert_1;
    std::string _path_vert_2;
    std::string _path_edge_vec_1;
    std::string _path_edge_vec_2;
    std::string _path_normal_1;
    std::string _path_normal_2;
    std::string _path_normal_3;
    std::string _path_ruling_0;
    std::string _path_ruling_1;
    std::string _path_ruling_2;
    std::string _path_dev;

    mutable std::ofstream _outFile;
    mutable std::ofstream _outFile_vert_1;
    mutable std::ofstream _outFile_vert_2;
    mutable std::ofstream _outFile_edge_vec_1;
    mutable std::ofstream _outFile_edge_vec_2;
    mutable std::ofstream _outFile_normal_1;
    mutable std::ofstream _outFile_normal_2;
    mutable std::ofstream _outFile_normal_3;
    mutable std::ofstream _outFile_ruling_0;
    mutable std::ofstream _outFile_ruling_1;
    mutable std::ofstream _outFile_ruling_2;
    mutable std::ofstream _outFile_dev;

public:

    ElasticDevelopabilityCentroidEnergy(const QuadMeshTopologySaver &quadTopol,
                                        const MeshTopologySaver &triangleTopol,
                                        const VectorType &baseGeometryTri,
                                        VectorType factors_elasticity_dev,
                                        VectorType factors_mem_bend,
                                        Constraint<ConfiguratorType> &constraint,
                                        VarsIdx<ConfiguratorType> variablesIdx,
                                        MatrixType W)
        : _quadTopol(quadTopol), 
          _constraint(constraint),
          _factors_elasticity_dev(factors_elasticity_dev),
          _factors_mem_bend(factors_mem_bend),
          _variablesIdx(variablesIdx),
          _W(W),
          _E_bend(triangleTopol, baseGeometryTri, true),
          _E_mem(triangleTopol, baseGeometryTri, true),
          _num_vertices(quadTopol.getNumVertices()),
          _quadTriConverter(quadTopol.getNumVertices())
        {
            counter = 0;
        }

        // centroidVars contains in the first components the centroid vertices 
        // instead of the quad vertices in the normal Vars in the quad examples
        // Further, the centroid vertices are ordered component major as in all the TriMesh examples
        // Will convert this to point major format for quad meshes via permutation
        // The components after the centroid vertices are the other quad optimization variables
        void apply(const VectorType &vars, RealType &energy) const override
        {
            // First, evaluate the elastic energies on the centroid vertices
            VectorType quad_vertices = vars.segment(_variablesIdx["vertices"],3*_num_vertices);
            VectorType tri_vertices(3*_num_vertices);
            _quadTriConverter.convertQuadToTri(quad_vertices, tri_vertices);
            
            RealType elasticCentroidEnergy_mem, elasticCentroidEnergy_bend, elasticCentroidEnergy;
            if(!(_factors_elasticity_dev[0] == 0.0))
            {
                _E_bend.apply(tri_vertices, elasticCentroidEnergy_bend);
                _E_mem.apply(tri_vertices, elasticCentroidEnergy_mem);
                elasticCentroidEnergy = _factors_mem_bend[0]*elasticCentroidEnergy_mem + _factors_mem_bend[1]*elasticCentroidEnergy_bend;
            }
            else{
                elasticCentroidEnergy = 0.0;
            }
            // Calculate the constraint
            VectorType constraint(_W.rows());
            if(!_factors_elasticity_dev[1] == 0.0)
            {
                _constraint.apply(vars, constraint);
            }
            else{
                constraint.setZero();
            }

            energy =  _factors_elasticity_dev[0]*elasticCentroidEnergy + _factors_elasticity_dev[1]*constraint.dot(_W*constraint);
            std::cout<<"energy elastic: "<<_factors_elasticity_dev[0]*elasticCentroidEnergy<<"; energy constraint: "<<_factors_elasticity_dev[1]*constraint.dot(_W*constraint)<<"; total energy: "<<energy<<std::endl;
        }
};

// Now, the Gradient of the DevelopabilityCentroidEnergy
template<typename ConfiguratorType>
class ElasticDevelopabilityCentroidGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> 
{
protected:
    typedef typename DefaultConfigurator::RealType   RealType;
    typedef typename DefaultConfigurator::VectorType VectorType;
    typedef typename DefaultConfigurator::VecType    VecType;
    typedef typename DefaultConfigurator::MatType    MatType;

    NonlinearMembraneGradient<ConfiguratorType> _DE_mem;
    SimpleBendingGradient<ConfiguratorType> _DE_bend;

    const QuadMeshTopologySaver& _quadTopol;
    const MeshTopologySaver &_triangleTopol;

    const MatrixType _W;

    const Constraint<ConfiguratorType> &_constraint;
    const ConstraintsGradient<ConfiguratorType> &_DConstraint;

    const QuadTriConverter<ConfiguratorType> _quadTriConverter;

    VectorType _factors;

    const VectorType _factors_elasticity_dev;
    const VectorType _factors_mem_bend;

    const VarsIdx<ConfiguratorType> _variablesIdx;

    int _num_vertices;

public:
    ElasticDevelopabilityCentroidGradient(const QuadMeshTopologySaver &quadTopol,
                                           const MeshTopologySaver &triangleTopol,
                                           const VectorType &baseGeometryTri,
                                           VectorType factors_elasticity_dev,
                                           VectorType factors_mem_bend,
                                           Constraint<ConfiguratorType> &constraint,
                                           ConstraintGrad<ConfiguratorType> &DConstraint,
                                           VarsIdx<ConfiguratorType> &variablesIdx,
                                           MatrixType W)
        : _quadTopol(quadTopol), 
          _triangleTopol(triangleTopol),
          _DConstraint(DConstraint),
          _constraint(constraint),
          _factors_elasticity_dev(factors_elasticity_dev),
          _factors_mem_bend(factors_mem_bend),
          _W(W),
          _variablesIdx(variablesIdx),
          _DE_bend(triangleTopol, baseGeometryTri),
          _DE_mem(triangleTopol, baseGeometryTri),
          _num_vertices(quadTopol.getNumVertices()),
          _quadTriConverter(quadTopol.getNumVertices())
        {}

        void apply(const VectorType &vars, VectorType &gradient) const override
        {
            VectorType quad_vertices = vars.segment(_variablesIdx["vertices"],3*_num_vertices);
            VectorType tri_vertices(3*_num_vertices);
            _quadTriConverter.convertQuadToTri(quad_vertices, tri_vertices);

            VectorType DelasticTriEnergy_mem, DelasticTriEnergy_bend, DelasticTriEnergy;
            _DE_bend.apply(tri_vertices, DelasticTriEnergy_bend);
            _DE_mem.apply(tri_vertices, DelasticTriEnergy_mem);

            DelasticTriEnergy = _factors_mem_bend[0]*DelasticTriEnergy_mem + _factors_mem_bend[1]*DelasticTriEnergy_bend;

            // convert this to the quad format
            VectorType DelasticQuadEnergy(3*_num_vertices);
            _quadTriConverter.convertTriToQuad(DelasticTriEnergy, DelasticQuadEnergy);

            VectorType constraint;
            _constraint.apply(vars, constraint);

            // Calculate the constraint
            MatrixType constraintGrad;
            _DConstraint.apply(vars, constraintGrad);
            
            VectorType constraintGrad_sqrd_1 = (constraintGrad.transpose()*_W) * constraint;
            VectorType constraintGrad_sqrd_2 = ((_W*constraintGrad).transpose()*constraint).transpose();
            VectorType constraintGrad_sqrd = constraintGrad_sqrd_1 + constraintGrad_sqrd_2;

            constraintGrad_sqrd *= _factors_elasticity_dev[1];

            constraintGrad_sqrd.segment(_variablesIdx["vertices"], 3*_num_vertices) += _factors_elasticity_dev[0]*DelasticQuadEnergy;

            gradient = constraintGrad_sqrd;
            //std::cout<<"gradient energy elastic: "<< _factors_elasticity_dev[0]*final_elastic_energies.norm() <<"; "<<" gradient energy constraint: "<< _factors_elasticity_dev[1]*final_constraints.norm() <<"; total gradient: "<<gradient.norm()<<std::endl;
        }
};

template<typename ConfiguratorType>
class ElasticDevelopabilityCentroidHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> 
{
protected:
    typedef typename DefaultConfigurator::RealType   RealType;
    typedef typename DefaultConfigurator::VectorType VectorType;
    typedef typename DefaultConfigurator::VecType    VecType;
    typedef typename DefaultConfigurator::MatType    MatType;
    typedef typename DefaultConfigurator::SparseMatrixType MatrixType;

    NonlinearMembraneHessian<ConfiguratorType> _D2E_mem;
    SimpleBendingHessian<ConfiguratorType> _D2E_bend;

    const QuadMeshTopologySaver& _quadTopol;
    const MeshTopologySaver &_centroidTopol;

    const VectorType &_centroidBaseGeometry;

    const MatrixType _W;

    const Constraint<ConfiguratorType> &_constraint;
    const ConstraintsGradient<ConfiguratorType> &_DConstraint;

    const CoordConverter<ConfiguratorType> _coordConverter;
    const VectorType _factors_elasticity_dev;
    const VectorType _factors_mem_bend;

public:
    ElasticDevelopabilityCentroidHessian(const QuadMeshTopologySaver &quadTopol,
                                           const MeshTopologySaver &centroidTopol,
                                           const VectorType &quadBaseGeometry,
                                           const VectorType &centroidBaseGeometry,
                                           VectorType &factors_elasticity_dev,
                                           VectorType &factors_mem_bend,
                                           Constraint<ConfiguratorType> &constraint,
                                           ConstraintGrad<ConfiguratorType> &DConstraint,
                                           MatrixType W)
        : _quadTopol(quadTopol), 
          _centroidTopol(centroidTopol),
          _centroidBaseGeometry(centroidBaseGeometry),
          _DConstraint(DConstraint),
          _constraint(constraint),
          _coordConverter(quadTopol, centroidTopol, quadBaseGeometry, centroidBaseGeometry),
          _factors_elasticity_dev(factors_elasticity_dev),
          _factors_mem_bend(factors_mem_bend),
          _W(W),
          _D2E_bend(centroidTopol, centroidBaseGeometry),
          _D2E_mem(centroidTopol, centroidBaseGeometry)
        {}

    void apply(const VectorType &centroidVars, MatrixType &hessian) const override
    {
        // First, evaluate the elastic energies on the centroid vertices
        VectorType centroidVertices = centroidVars.segment(0, 3*_centroidTopol.getNumVertices());
        MatrixType D2elasticCentroidEnergy_mem, D2elasticCentroidEnergy_bend;
        _D2E_bend.apply(centroidVertices, D2elasticCentroidEnergy_bend);
        _D2E_mem.apply(centroidVertices, D2elasticCentroidEnergy_mem);

        MatrixType D2elasticCentroidEnergy = _factors_mem_bend[0]*D2elasticCentroidEnergy_mem + _factors_mem_bend[1]*D2elasticCentroidEnergy_bend;

        VectorType quadVertices(3*_quadTopol.getNumVertices());
        _coordConverter.convertCentroidToQuad(centroidVertices, quadVertices);

        // Now, calculate the Quad constraints using this vars vector
        VectorType vars(centroidVars.size() - 3*(_centroidTopol.getNumVertices() - _quadTopol.getNumVertices()));
        vars.segment(0, 3*_quadTopol.getNumVertices()) = quadVertices;
        vars.segment(3*_quadTopol.getNumVertices(), vars.size() - 3*_quadTopol.getNumVertices()) = centroidVars.segment(3*_centroidTopol.getNumVertices(), vars.size() - 3*_quadTopol.getNumVertices());   

        MatrixType Hessian(centroidVars.size(), centroidVars.size());
        std::vector<Eigen::Triplet<double>> tripletList;
        Hessian.setZero();
        assignSparseBlockInplace(Hessian, D2elasticCentroidEnergy, 0, 0, tripletList);
        hessian = Hessian;
    }
};

template<typename ConfiguratorType>
class QuadElasticEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>
{
    protected:
        typedef typename DefaultConfigurator::RealType   RealType;
        typedef typename DefaultConfigurator::VectorType VectorType;
        typedef typename DefaultConfigurator::VecType    VecType;
        typedef typename DefaultConfigurator::MatType    MatType;

        NonlinearMembraneEnergy<ConfiguratorType> _E_mem;
        SimpleBendingEnergy<ConfiguratorType> _E_bend;

        const QuadMeshTopologySaver _quadTopol;

        const ConstraintSqrdReduced<ConfiguratorType> _constraint;

        const VectorType _factors_elasticity_dev;
        const VectorType _factors_mem_bend;

        const QuadTriConverter<ConfiguratorType> _quadTriConverter;

    public:
        QuadElasticEnergy(const QuadMeshTopologySaver &quadTopol,
                            const MeshTopologySaver &triangleTopol,
                            const VectorType &baseGeomTri,
                            VectorType factors_elasticity_dev,
                            VectorType factors_mem_bend)
            : _quadTopol(quadTopol),
            _constraint(quadTopol),
            _factors_elasticity_dev(factors_elasticity_dev),
            _factors_mem_bend(factors_mem_bend),
            _E_bend(triangleTopol, baseGeomTri, true),
            _E_mem(triangleTopol, baseGeomTri, true),
            _quadTriConverter(quadTopol.getNumVertices())
            {}

        void apply(const VectorType &quad_vertices, RealType &energy) const override
        {
            // First, evaluate the elastic energies on the centroid vertices
            RealType membraneEnergy = 0.0;
            RealType bendingEnergy = 0.0;

            VectorType tri_vertices(3*_quadTopol.getNumVertices());
            _quadTriConverter.convertQuadToTri(quad_vertices, tri_vertices);
            
            if(!(_factors_elasticity_dev[0] == 0.0))
            {
                if(!(_factors_mem_bend[0] == 0.0)){
                    _E_mem.apply(tri_vertices, membraneEnergy);
                }
                if(!(_factors_mem_bend[1] == 0.0)){
                    _E_bend.apply(tri_vertices, bendingEnergy);
                }
            }
            RealType elasticEnergies = _factors_mem_bend[0]*membraneEnergy + _factors_mem_bend[1]*bendingEnergy;

            RealType developabilityEnergy = 0.0;
            if(!(_factors_elasticity_dev[1] == 0.0))
            {
                _constraint.apply(quad_vertices, developabilityEnergy);
            }

            energy = _factors_elasticity_dev[0]*elasticEnergies + _factors_elasticity_dev[1]*developabilityEnergy;
            std::cout<<"energy elastic: "<<_factors_elasticity_dev[0]*elasticEnergies<<"; energy constraint: "<<_factors_elasticity_dev[1]*developabilityEnergy<<"; total energy: "<<energy<<std::endl;
        }
};

template<typename ConfiguratorType>
class QuadElasticEnergyGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>
{
    protected:
        typedef typename DefaultConfigurator::RealType   RealType;
        typedef typename DefaultConfigurator::VectorType VectorType;
        typedef typename DefaultConfigurator::VecType    VecType;

        QuadMeshTopologySaver _quadTopol;
        MeshTopologySaver _triangleTopol;
        
        VectorType _factors_elasticity_dev;
        VectorType _factors_mem_bend;

        NonlinearMembraneGradientDef<ConfiguratorType> _DE_mem;
        SimpleBendingGradientDef<ConfiguratorType> _DE_bend;
        
        ConstraintSqrdReducedGradient<ConfiguratorType> _DConstraint;

        QuadTriConverter<ConfiguratorType> _quadTriConverter;

    public:
        QuadElasticEnergyGradient(const QuadMeshTopologySaver &quadTopol,
                                    const MeshTopologySaver &triangleTopol,
                                    const VectorType &triangleBaseGeometry,
                                    VectorType factors_elasticity_dev,
                                    VectorType factors_mem_bend)
            : _quadTopol(quadTopol),
            _triangleTopol(triangleTopol),
            _factors_elasticity_dev(factors_elasticity_dev),
            _factors_mem_bend(factors_mem_bend),
            _DE_bend(triangleTopol, triangleBaseGeometry, true),
            _DE_mem(triangleTopol, triangleBaseGeometry, true),
            _DConstraint(quadTopol),
                _quadTriConverter(quadTopol.getNumVertices())
            {
                std::cout<<"basegeomtri size: "<<triangleBaseGeometry.size()<<std::endl;
            }

        void apply(const VectorType &quad_vertices, VectorType &gradient) const override
        {
            // First, evaluate the elastic energies on the centroid vertices
            VectorType DmembraneEnergy = VectorType::Zero(quad_vertices.size());
            VectorType DbendingEnergy = VectorType::Zero(quad_vertices.size());

            VectorType tri_vertices(3*_quadTopol.getNumVertices());
            _quadTriConverter.convertQuadToTri(quad_vertices, tri_vertices);

            
            if(!(_factors_elasticity_dev[0] == 0.0))
            {
                if(!(_factors_mem_bend[0] == 0.0)){
                    _DE_mem.apply(tri_vertices, DmembraneEnergy);
                }
                if(!(_factors_mem_bend[1] == 0.0)){
                    _DE_bend.apply(tri_vertices, DbendingEnergy);
                }
            }

            VectorType DelasticEnergies = _factors_mem_bend[0]*DmembraneEnergy + _factors_mem_bend[1]*DbendingEnergy;
            VectorType DelasticEnergies_quad(3*_quadTopol.getNumVertices());
            _quadTriConverter.convertTriToQuad(DelasticEnergies, DelasticEnergies_quad);

            VectorType DConstraint = VectorType::Zero(3*_quadTopol.getNumVertices());

            if(!(_factors_elasticity_dev[1] == 0.0))
            {
                _DConstraint.apply(quad_vertices, DConstraint);
            }

            gradient = _factors_elasticity_dev[0]*DelasticEnergies_quad + _factors_elasticity_dev[1]*DConstraint;
        }
};

template<typename ConfiguratorType>
class QuadElasticEnergyHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>
{
    protected:
        typedef typename DefaultConfigurator::RealType   RealType;
        typedef typename DefaultConfigurator::VectorType VectorType;
        typedef typename DefaultConfigurator::SparseMatrixType MatrixType;

        QuadMeshTopologySaver _quadTopol;
        MeshTopologySaver _triangleTopol;
        VectorType _factors_elasticity_dev;
        VectorType _factors_mem_bend;
        NonlinearMembraneHessianDef<ConfiguratorType> _D2E_mem;
        SimpleBendingHessianDef<ConfiguratorType> _D2E_bend;
        ConstraintSqrdReducedHessian<ConfiguratorType> _D2Constraint;

        QuadTriConverter<ConfiguratorType> _quadTriConverter;

    public:
        QuadElasticEnergyHessian(const QuadMeshTopologySaver &quadTopol,
                                    const MeshTopologySaver &triangleTopol,
                                    const VectorType &triangleBaseGeometry,
                                    VectorType factors_elasticity_dev,
                                    VectorType factors_mem_bend)
            : _quadTopol(quadTopol),
            _triangleTopol(triangleTopol),
            _factors_elasticity_dev(factors_elasticity_dev),
            _factors_mem_bend(factors_mem_bend),
            _D2E_bend(triangleTopol, triangleBaseGeometry, true),
            _D2E_mem(triangleTopol, triangleBaseGeometry, true),
            _D2Constraint(quadTopol),
            _quadTriConverter(quadTopol.getNumVertices())
            {}

        void apply(const VectorType &quad_vertices, MatrixType &hessian) const override
        {
            VectorType tri_vertices(3*_quadTopol.getNumVertices());
            _quadTriConverter.convertQuadToTri(quad_vertices, tri_vertices);

            MatrixType D2elasticCentroidEnergy_mem(3*_quadTopol.getNumVertices(), 3*_quadTopol.getNumVertices());
            MatrixType D2elasticCentroidEnergy_bend(3*_quadTopol.getNumVertices(), 3*_quadTopol.getNumVertices());
            MatrixType D2Constraint(3*_quadTopol.getNumVertices(), 3*_quadTopol.getNumVertices());
            if(_factors_elasticity_dev[0] != 0.0)
            {
                if(_factors_mem_bend[0] != 0.0){
                    _D2E_mem.apply(tri_vertices, D2elasticCentroidEnergy_mem);
                }
                if(_factors_mem_bend[1] != 0.0){
                    _D2E_bend.apply(tri_vertices, D2elasticCentroidEnergy_bend);
                }
            }
            else{
                D2elasticCentroidEnergy_mem.setZero();
                D2elasticCentroidEnergy_bend.setZero();
            }

            _quadTriConverter.convertTriHessianToQuad(D2elasticCentroidEnergy_mem, D2elasticCentroidEnergy_mem);
            _quadTriConverter.convertTriHessianToQuad(D2elasticCentroidEnergy_bend, D2elasticCentroidEnergy_bend);

            if(_factors_elasticity_dev[1] != 0.0)
            {
                _D2Constraint.apply(quad_vertices, D2Constraint);
            }
            else{
                D2Constraint.setZero();
            }
            hessian = _factors_elasticity_dev[0]*(_factors_mem_bend[0]*D2elasticCentroidEnergy_mem + _factors_mem_bend[1]*D2elasticCentroidEnergy_bend) + _factors_elasticity_dev[1]*D2Constraint;

        }

};

template<typename ConfiguratorType>
class CombinedQuadEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>
{
    protected:
        typedef typename DefaultConfigurator::RealType   RealType;
        typedef typename DefaultConfigurator::VectorType VectorType;
        typedef typename DefaultConfigurator::VecType    VecType;
        typedef typename DefaultConfigurator::MatType    MatType;
        typedef typename DefaultConfigurator::SparseMatrixType MatrixType;


        QuadMeshTopologySaver _quadTopol;
        MeshTopologySaver _centroidTopol;

        std::shared_ptr<ConstraintSqrdReduced<DefaultConfigurator>> _constraint_ptr;

        std::shared_ptr<SimpleBendingEnergy<DefaultConfigurator>> _E_bend_ptr;
        std::shared_ptr<NonlinearMembraneEnergy<DefaultConfigurator>> _E_mem_ptr;

        const VectorType _factors_mem_bend;
        const VectorType _factors_elasticity_dev;

        VectorType _quadBaseGeometry;
        VectorType _centroidBaseGeometry;

        mutable TriMesh  _centroidMesh;
        mutable MyMesh _quadMesh;

        CoordConverter<ConfiguratorType> _coordConverter;

    public:
        CombinedQuadEnergy(const QuadMeshTopologySaver quadTopol,
                            const MeshTopologySaver centroidTopol,
                            VectorType &quadBaseGeometry,
                            VectorType &centroidBaseGeometry,
                            MyMesh quadMesh,
                            TriMesh centroidMesh,
                            VectorType factors_elasticity_dev,
                            VectorType factors_mem_bend,
                            CoordConverter<ConfiguratorType> coordConverter,
                            std::shared_ptr<SimpleBendingEnergy<DefaultConfigurator>> E_bend_ptr,
                            std::shared_ptr<NonlinearMembraneEnergy<DefaultConfigurator>> E_mem_ptr,
                            std::shared_ptr<ConstraintSqrdReduced<DefaultConfigurator>> constraint_ptr):
          _quadBaseGeometry(quadBaseGeometry),
          _factors_elasticity_dev(factors_elasticity_dev),
          _factors_mem_bend(factors_mem_bend),
          _quadTopol(quadTopol),
          _centroidTopol(centroidTopol),
            _quadMesh(quadMesh),
            _centroidMesh(centroidMesh),
            _coordConverter(coordConverter),
            _E_bend_ptr(std::move(E_bend_ptr)),
            _E_mem_ptr(std::move(E_mem_ptr)),
            _constraint_ptr(std::move(constraint_ptr))
        {
        }

        void apply(const VectorType &quadGeom, RealType &energy) const override
        {
            VectorType centroidGeom(3*_centroidTopol.getNumVertices());
            _coordConverter.convertQuadToCentroid2(quadGeom, centroidGeom, _centroidMesh);

            // Evaluate the elastic energies on the centroid vertices
            RealType membraneEnergy = 0.0;
            RealType bendingEnergy = 0.0;
            if(!(_factors_elasticity_dev[0] == 0.0))
            {
                if(!(_factors_mem_bend[0] == 0.0)){
                    _E_mem_ptr->apply(centroidGeom, membraneEnergy);
                }
                if(!(_factors_mem_bend[1] == 0.0)){
                    _E_bend_ptr->apply(centroidGeom, bendingEnergy);
                }
            }
            RealType elasticEnergies = _factors_mem_bend[0]*membraneEnergy + _factors_mem_bend[1]*bendingEnergy;

            RealType developabilityEnergy = 0.0;
            if(!(_factors_elasticity_dev[1] == 0.0))
            {
                _constraint_ptr->apply(quadGeom, developabilityEnergy);
            }

            energy = _factors_elasticity_dev[0]*elasticEnergies + _factors_elasticity_dev[1]*developabilityEnergy;
            std::cout<<"developability energy: "<<_factors_elasticity_dev[1]*developabilityEnergy<<"; elastic energies: "<<elasticEnergies<<"; total energy: "<<energy<<std::endl;
        }
};


template<typename ConfiguratorType>
class CombinedQuadEnergyGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>
{
    protected:
        typedef typename DefaultConfigurator::RealType   RealType;
        typedef typename DefaultConfigurator::VectorType VectorType;
        typedef typename DefaultConfigurator::VecType    VecType;
        typedef typename DefaultConfigurator::MatType    MatType;
        typedef typename DefaultConfigurator::SparseMatrixType MatrixType;

        QuadMeshTopologySaver _quadTopol;
        MeshTopologySaver _centroidTopol;

        std::shared_ptr<ConstraintSqrdReducedGradient<DefaultConfigurator>> _DConstraint_ptr;
        std::shared_ptr<SimpleBendingGradientDef<DefaultConfigurator>> _DE_bend_ptr;
        std::shared_ptr<NonlinearMembraneGradientDef<DefaultConfigurator>> _DE_mem_ptr;

        const VectorType _factors_mem_bend;
        const VectorType _factors_elasticity_dev;

        VectorType &_quadBaseGeometry;
        VectorType &_centroidBaseGeometry;

        mutable TriMesh  _centroidMesh;
        mutable MyMesh _quadMesh;

        CoordConverter<ConfiguratorType> _coordConverter;

    public:
        CombinedQuadEnergyGradient(const QuadMeshTopologySaver quadTopol,
                                    const MeshTopologySaver centroidTopol,
                                    VectorType &quadBaseGeometry,
                                    VectorType &centroidBaseGeometry,
                                    MyMesh &quadMesh,
                                    TriMesh &centroidMesh,
                                    VectorType factors_elasticity_dev,
                                    VectorType factors_mem_bend,
                                    CoordConverter<ConfiguratorType> coordConverter,
                                    std::shared_ptr<SimpleBendingGradientDef<DefaultConfigurator>> DE_bend_ptr,
                                    std::shared_ptr<NonlinearMembraneGradientDef<DefaultConfigurator>> DE_mem_ptr,
                                    std::shared_ptr<ConstraintSqrdReducedGradient<DefaultConfigurator>> DConstraint_ptr):
            _quadBaseGeometry(quadBaseGeometry),
            _centroidBaseGeometry(centroidBaseGeometry),
            _factors_elasticity_dev(factors_elasticity_dev),
            _factors_mem_bend(factors_mem_bend),
            _quadTopol(quadTopol),
            _centroidTopol(centroidTopol),
            _quadMesh(quadMesh),
            _centroidMesh(centroidMesh),
            _coordConverter(coordConverter),
            _DConstraint_ptr(std::move(DConstraint_ptr)),
            _DE_bend_ptr(std::move(DE_bend_ptr)),
            _DE_mem_ptr(std::move(DE_mem_ptr))
        {}

        void apply(const VectorType &quadGeom, VectorType &gradient) const override
        {
            // Convert the quad geometry to centroid geometry
             // Convert the quad geometry to centroid geometry
            VectorType centroidGeom(3*_centroidTopol.getNumVertices());
            _coordConverter.convertQuadToCentroid2(quadGeom, centroidGeom, _centroidMesh);

            // Evaluate the elastic energies on the centroid vertices
            VectorType DmembraneEnergy = VectorType::Zero(centroidGeom.size());
            VectorType DbendingEnergy = VectorType::Zero(centroidGeom.size());
            
            if(!(_factors_elasticity_dev[0] == 0.0))
            {
                if(!(_factors_mem_bend[0] == 0.0)){
                    _DE_mem_ptr->apply(centroidGeom, DmembraneEnergy);
                }
                if(!(_factors_mem_bend[1] == 0.0)){
                    _DE_bend_ptr->apply(centroidGeom, DbendingEnergy);
                }
            }
            VectorType DelasticEnergies = _factors_mem_bend[0]*DmembraneEnergy + _factors_mem_bend[1]*DbendingEnergy;
            VectorType DelasticEnergies_quad(3*_quadTopol.getNumVertices());
            _coordConverter.convertCentroidToQuad(DelasticEnergies, DelasticEnergies_quad);

            VectorType Dconstraint = VectorType::Zero(3*_quadTopol.getNumVertices());
            if(!(_factors_elasticity_dev[1] == 0.0))
            {
                _DConstraint_ptr->apply(quadGeom, Dconstraint);
            }

            gradient = _factors_elasticity_dev[0]*DelasticEnergies_quad + _factors_elasticity_dev[1]*Dconstraint;
        }
};

/*
template<typename ConfiguratorType>
class CombinedQuadEnergyHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>
{
    protected:
        typedef typename DefaultConfigurator::RealType   RealType;
        typedef typename DefaultConfigurator::VectorType VectorType;
        typedef typename DefaultConfigurator::VecType    VecType;
        typedef typename DefaultConfigurator::MatType    MatType;
        typedef typename DefaultConfigurator::SparseMatrixType MatrixType;

        QuadMeshTopologySaver _quadTopol;
        MeshTopologySaver _centroidTopol;

        std::shared_ptr<ConstraintSqrdReducedHessian<DefaultConfigurator>> _D2Constraint_ptr;
        std::shared_ptr<SimpleBendingHessian<DefaultConfigurator>> _D2E_bend_ptr;
        std::shared_ptr<NonlinearMembraneHessian<DefaultConfigurator>> _D2E_mem_ptr;

        const VectorType _factors_mem_bend;
        const VectorType _factors_elasticity_dev;

        VectorType &_quadBaseGeometry;
        VectorType &_centroidBaseGeometry;

        mutable TriMesh  _centroidMesh;
        mutable MyMesh _quadMesh;

        CoordConverter<ConfiguratorType> _coordConverter;

    public:
        CombinedQuadEnergyHessian(const QuadMeshTopologySaver quadTopol,
                                    const MeshTopologySaver centroidTopol,
                                    VectorType &quadBaseGeometry,
                                    VectorType &centroidBaseGeometry,
                                    MyMesh &quadMesh,
                                    TriMesh &centroidMesh,
                                    VectorType factors_elasticity_dev,
                                    VectorType factors_mem_bend,
                                    CoordConverter<ConfiguratorType> coordConverter,
                                    std::shared_ptr<SimpleBendingHessian<DefaultConfigurator>> D2E_bend_ptr,
                                    std::shared_ptr<NonlinearMembraneHessian<DefaultConfigurator>> D2E_mem_ptr,
                                    std::shared_ptr<ConstraintSqrdReducedHessian<DefaultConfigurator>> D2Constraint_ptr):
            _quadBaseGeometry(quadBaseGeometry),
            _centroidBaseGeometry(centroidBaseGeometry),
            _factors_elasticity_dev(factors_elasticity_dev),
            _factors_mem_bend(factors_mem_bend),
            _quadTopol(quadTopol),
            _centroidTopol(centroidTopol),
            _quadMesh(quadMesh),
            _centroidMesh(centroidMesh),
            _coordConverter(coordConverter),
            _D2Constraint_ptr(std::move(D2Constraint_ptr)),
            _D2E_bend_ptr(std::move(D2E_bend_ptr)),
            _D2E_mem_ptr(std::move(D2E_mem_ptr))
        {}

        void apply(const VectorType &quadGeom, MatrixType &hessian) const override
        {
            // Convert the quad geometry to centroid geometry
            VectorType centroidGeom(3*_centroidTopol.getNumVertices());
            _coordConverter.convertQuadToCentroid2(quadGeom, centroidGeom, _centroidMesh);

            // Evaluate the elastic energies on the centroid vertices
            MatrixType D2membraneEnergy = MatrixType::Zero(centroidGeom.size(), centroidGeom.size());
            MatrixType D2bendingEnergy = MatrixType::Zero(centroidGeom.size(), centroidGeom.size());
            
            if(!(_factors_elasticity_dev[0] == 0.0))
            {
                if(!(_factors_mem_bend[0] == 0.0)){
                    _D2E_mem_ptr->apply(centroidGeom, D2membraneEnergy);
                }
                if(!(_factors_mem_bend[1] == 0.0)){
                    _D2E_bend_ptr->apply(centroidGeom, D2bendingEnergy);
                }
            }
            MatrixType D2elasticEnergies = _factors_mem_bend[0]*D2membraneEnergy + _factors_mem_bend[1]*D2bendingEnergy;

            // Now, need to convert this hessian from centroid to quad coordinates
            // This is done by converting each column of the hessian
            MatrixType D2elasticEnergies_quad = MatrixType::Zero(3*_quadTopol.getNumVertices(), 3*_quadTopol.getNumVertices());
            for(int i=0; i<D2elasticEnergies.cols(); i++)
            {
                VectorType col = D2elasticEnergies.col(i);
                VectorType col_quad(3*_quadTopol.getNumVertices());
                _coordConverter.convertCentroidToQuad(col, col_quad);
                D2elasticEnergies_quad.col(i) = col_quad;
            }

            MatrixType D2constraint = MatrixType::Zero(3*_quadTopol.getNumVertices(), 3*_quadTopol.getNumVertices());
            if(!(_factors_elasticity_dev[1] == 0.0))
            {
                _D2Constraint_ptr->apply(quadGeom, D2constraint);
            }

            hessian = _factors_elasticity_dev[0]*D2elasticEnergies
*/