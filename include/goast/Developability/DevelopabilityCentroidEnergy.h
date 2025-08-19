#pragma once

#include <goast/Core.h>

#include <goast/Developability/Constraints.h>
#include <goast/Developability/BoundaryDOFS_quad.h>
#include <goast/DiscreteShells.h>
#include <fstream>
#include <iostream>


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

        std::vector<int> getNoQuadCentroidVertices() const
        {
            return _noQuadCentroidVertices;
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
    const MeshTopologySaver &_centroidTopol;

    const VectorType &_centroidBaseGeometry;

    const MatrixType _W;

    const Constraint<ConfiguratorType> &_constraint;

    const CoordConverter<ConfiguratorType> _coordConverter;

    const VectorType _factors_elasticity_dev;
    const VectorType _factors_mem_bend;

    const VarsIdx<ConfiguratorType> _variablesIdx;

    mutable int counter;

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
                                        const MeshTopologySaver &centroidTopol,
                                        const VectorType &quadBaseGeometry,
                                        const VectorType &centroidBaseGeometry,
                                        VectorType factors_elasticity_dev,
                                        VectorType factors_mem_bend,
                                        Constraint<ConfiguratorType> &constraint,
                                        VarsIdx<ConfiguratorType> variablesIdx,
                                        MatrixType W)
        : _quadTopol(quadTopol), 
          _centroidTopol(centroidTopol),
          _centroidBaseGeometry(centroidBaseGeometry),
          _constraint(constraint),
          _coordConverter(quadTopol, centroidTopol, quadBaseGeometry, centroidBaseGeometry),
          _factors_elasticity_dev(factors_elasticity_dev),
          _factors_mem_bend(factors_mem_bend),
          _variablesIdx(variablesIdx),
          _W(W),
          _E_bend(centroidTopol, centroidBaseGeometry, true),
          _E_mem(centroidTopol, centroidBaseGeometry, true)
        {
            counter = 0;
        }

        // centroidVars contains in the first components the centroid vertices 
        // instead of the quad vertices in the normal Vars in the quad examples
        // Further, the centroid vertices are ordered component major as in all the TriMesh examples
        // Will convert this to point major format for quad meshes via permutation
        // The components after the centroid vertices are the other quad optimization variables
        void apply(const VectorType &centroidVars, RealType &energy) const override
        {
            // First, evaluate the elastic energies on the centroid vertices
            VectorType centroidVertices = centroidVars.segment(_variablesIdx["vertices"],3*_centroidTopol.getNumVertices());
            
            RealType elasticCentroidEnergy_mem, elasticCentroidEnergy_bend, elasticCentroidEnergy;
            if(!(_factors_elasticity_dev[0] == 0.0))
            {
                _E_bend.apply(centroidVertices, elasticCentroidEnergy_bend);
                _E_mem.apply(centroidVertices, elasticCentroidEnergy_mem);
                elasticCentroidEnergy = _factors_mem_bend[0]*elasticCentroidEnergy_mem + _factors_mem_bend[1]*elasticCentroidEnergy_bend;
            }
            else{
                elasticCentroidEnergy = 0.0;
            }
            
            VectorType quadVertices(3*_quadTopol.getNumVertices());
            _coordConverter.convertCentroidToQuad(centroidVertices, quadVertices);

            //printVectorToFile(quadVertices, "/home/paul_johannssen/Desktop/masterarbeit/goast/goast/build/examples/quad_vertices_in_dev_energy.txt");

            // Now, calculate the Quad constraints using this vars vector
            VectorType vars(centroidVars.size() - 3*(_centroidTopol.getNumVertices() - _quadTopol.getNumVertices()));

            vars.segment(0, _variablesIdx["vertices"]) = centroidVars.segment(0, _variablesIdx["vertices"]);
            vars.segment(_variablesIdx["vertices"], 3*_quadTopol.getNumVertices()) = quadVertices;
            vars.segment(_variablesIdx["vertices"] + 3*_quadTopol.getNumVertices(), vars.size() - (3*_quadTopol.getNumVertices() + _variablesIdx["vertices"])) = centroidVars.segment(_variablesIdx["vertices"] + 3*_centroidTopol.getNumVertices(), centroidVars.size() - (_variablesIdx["vertices"] + 3*_centroidTopol.getNumVertices()));   

            // Calculate the constraint
            VectorType constraint(_W.rows());
            if(!_factors_elasticity_dev[1] == 0.0)
            {
                _constraint.apply(vars, constraint);
            }
            else{
                constraint.setZero();
            }

            // now, output all the constraints
            RealType product = constraint.dot(_W*constraint);
            _outFile.open(_path, std::ios::app);  // open in append mode
            _outFile << product << "\n";  // write scalar followed by newline
            _outFile.close();

            // Output the individual components of the constraint for debugging
             size_t num_vertices = _quadTopol.getNumVertices();
             size_t num_faces = _quadTopol.getNumFaces();
             size_t num_halfedges = 2*_quadTopol.getNumEdges();
             size_t num_bdryHalfEdges = _quadTopol.getBdryHalfEdges().size();
             size_t num_bdryFaces = _quadTopol.getBdryFaces().size();
 
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
 
             int start_vertex_0 = 0;
             int start_vertex_1 = start_vertex_0 + num_constraints_vert_1;
             int start_edge_vec_1 = start_vertex_1 + num_constraints_vert_2;
             int start_edge_vec_2 = start_edge_vec_1 + num_constraints_edge_vec_1;
             int start_normal_1 = start_edge_vec_2 + num_constraints_edge_vec_2;
             int start_normal_2 = start_normal_1 + num_constraints_normal_1;
             int start_normal_3 = start_normal_2 + num_constraints_normal_2;
             int start_ruling_0 = start_normal_3 + num_constraints_normal_3;
             int start_ruling_1 = start_ruling_0 + num_constraints_ruling_0;
             int start_ruling_2 = start_ruling_1 + num_constraints_ruling_1;
             int start_dev = start_ruling_2 + num_constraints_ruling_2;

            RealType constraint_vertex_1 = constraint.segment(start_vertex_0, num_constraints_vert_1).norm();
            RealType constraint_vertex_2 = constraint.segment(start_vertex_1, num_constraints_vert_2).norm();             

            RealType constraint_edge_vec_1 = constraint.segment(start_edge_vec_1, num_constraints_edge_vec_1).norm();
            RealType constraint_edge_vec_2 = constraint.segment(start_edge_vec_2, num_constraints_edge_vec_2).norm();

            RealType constraint_normal_1 = constraint.segment(start_normal_1, num_constraints_normal_1).norm();
            RealType constraint_normal_2 = constraint.segment(start_normal_2, num_constraints_normal_2).norm();

            RealType constraint_normal_3 = constraint.segment(start_normal_3, num_constraints_normal_3).norm();
            RealType constraint_ruling_0 = constraint.segment(start_ruling_0, num_constraints_ruling_0).norm();
            RealType constraint_ruling_1 = constraint.segment(start_ruling_1, num_constraints_ruling_1).norm();
            RealType constraint_ruling_2 = constraint.segment(start_ruling_2, num_constraints_ruling_2).norm();
            RealType constraint_dev = constraint.segment(start_dev, num_constraints_dev).norm();

            //printVectorToFile(vars, "/home/paul_johannssen/Desktop/masterarbeit/goast/goast/build/examples/init_2.txt");
            //printVectorToFile(constraint, "/home/paul_johannssen/Desktop/masterarbeit/goast/goast/build/examples/constraint.txt");

            energy =  _factors_elasticity_dev[0]*elasticCentroidEnergy + _factors_elasticity_dev[1]*constraint.dot(_W*constraint);
            std::cout<<"energy elastic: "<<_factors_elasticity_dev[0]*elasticCentroidEnergy<<"; energy constraint: "<<_factors_elasticity_dev[1]*constraint.dot(_W*constraint)<<"; total energy: "<<energy<<std::endl;
        }

        void apply2(const VectorType &centroidVars, RealType &energy)
        {
            // First, evaluate the elastic energies on the centroid vertices
            VectorType centroidVertices = centroidVars.segment(_variablesIdx["vertices"],3*_centroidTopol.getNumVertices());
            
            RealType elasticCentroidEnergy_mem, elasticCentroidEnergy_bend, elasticCentroidEnergy;
            if(!(_factors_elasticity_dev[0] == 0.0))
            {
                _E_bend.apply(centroidVertices, elasticCentroidEnergy_bend);
                _E_mem.apply(centroidVertices, elasticCentroidEnergy_mem);
                elasticCentroidEnergy = _factors_mem_bend[0]*elasticCentroidEnergy_mem + _factors_mem_bend[1]*elasticCentroidEnergy_bend;
            }
            else{
                elasticCentroidEnergy = 0.0;
            }
            
            VectorType quadVertices(3*_quadTopol.getNumVertices());
            _coordConverter.convertCentroidToQuad(centroidVertices, quadVertices);

            //printVectorToFile(quadVertices, "/home/paul_johannssen/Desktop/masterarbeit/goast/goast/build/examples/quad_vertices_in_dev_energy.txt");

            // Now, calculate the Quad constraints using this vars vector
            VectorType vars(centroidVars.size() - 3*(_centroidTopol.getNumVertices() - _quadTopol.getNumVertices()));

            vars.segment(0, _variablesIdx["vertices"]) = centroidVars.segment(0, _variablesIdx["vertices"]);
            vars.segment(_variablesIdx["vertices"], 3*_quadTopol.getNumVertices()) = quadVertices;
            vars.segment(_variablesIdx["vertices"] + 3*_quadTopol.getNumVertices(), vars.size() - (3*_quadTopol.getNumVertices() + _variablesIdx["vertices"])) = centroidVars.segment(_variablesIdx["vertices"] + 3*_centroidTopol.getNumVertices(), centroidVars.size() - (_variablesIdx["vertices"] + 3*_centroidTopol.getNumVertices()));   

            // Calculate the constraint
            VectorType constraint(_W.rows());
            if(!_factors_elasticity_dev[1] == 0.0)
            {
                _constraint.apply(vars, constraint);
            }
            else{
                constraint.setZero();
            }

            // now, output all the constraints
            RealType product = constraint.dot(_W*constraint);
            _outFile.open(_path, std::ios::app);  // open in append mode
            _outFile << product << "\n";  // write scalar followed by newline
            _outFile.close();

            //printVectorToFile(vars, "/home/paul_johannssen/Desktop/masterarbeit/goast/goast/build/examples/init_2.txt");
            //printVectorToFile(constraint, "/home/paul_johannssen/Desktop/masterarbeit/goast/goast/build/examples/constraint.txt");

            energy =  _factors_elasticity_dev[0]*elasticCentroidEnergy + _factors_elasticity_dev[1]*constraint.dot(_W*constraint);
            //std::cout<<"energy elastic: "<<_factors_elasticity_dev[0]*elasticCentroidEnergy<<"; energy constraint: "<<_factors_elasticity_dev[1]*constraint.dot(_W*constraint)<<"; total energy: "<<energy<<std::endl;
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
    const MeshTopologySaver &_centroidTopol;

    const VectorType &_centroidBaseGeometry;

    const MatrixType _W;

    const Constraint<ConfiguratorType> &_constraint;
    const ConstraintsGradient<ConfiguratorType> &_DConstraint;

    VectorType _factors;

    const CoordConverter<ConfiguratorType> _coordConverter;

    const VectorType _factors_elasticity_dev;
    const VectorType _factors_mem_bend;

    const VarsIdx<ConfiguratorType> _variablesIdx;

    // We are going to set the gradient of noQuadCentroidVertices to the average of the surrounding quad vertices
    const std::vector<int> _noQuadCentroidVertices;

public:
    ElasticDevelopabilityCentroidGradient(const QuadMeshTopologySaver &quadTopol,
                                           const MeshTopologySaver &centroidTopol,
                                           const VectorType &quadBaseGeometry,
                                           const VectorType &centroidBaseGeometry,
                                           VectorType factors_elasticity_dev,
                                           VectorType factors_mem_bend,
                                           Constraint<ConfiguratorType> &constraint,
                                           ConstraintGrad<ConfiguratorType> &DConstraint,
                                           VarsIdx<ConfiguratorType> &variablesIdx,
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
          _variablesIdx(variablesIdx),
          _DE_bend(centroidTopol, centroidBaseGeometry),
          _DE_mem(centroidTopol, centroidBaseGeometry)
        {}

        void apply(const VectorType &centroidVars, VectorType &gradient) const override
        {
            // First, evaluate the elastic energies on the centroid vertices
            VectorType centroidVertices = centroidVars.segment(_variablesIdx["vertices"], 3*_centroidTopol.getNumVertices());
            VectorType DelasticCentroidEnergy_mem, DelasticCentroidEnergy_bend, DelasticCentroidEnergy;
            _DE_bend.apply(centroidVertices, DelasticCentroidEnergy_bend);
            _DE_mem.apply(centroidVertices, DelasticCentroidEnergy_mem);

            DelasticCentroidEnergy = _factors_mem_bend[0]*DelasticCentroidEnergy_mem + _factors_mem_bend[1]*DelasticCentroidEnergy_bend;

            VectorType quadVertices(3*_quadTopol.getNumVertices());
            _coordConverter.convertCentroidToQuad(centroidVertices, quadVertices);

            // Now, calculate the Quad constraints using this vars vector
            VectorType vars(centroidVars.size() - 3*(_centroidTopol.getNumVertices() - _quadTopol.getNumVertices()));

            vars.segment(0, _variablesIdx["vertices"]) = centroidVars.segment(0, _variablesIdx["vertices"]);
            vars.segment(_variablesIdx["vertices"], 3*_quadTopol.getNumVertices()) = quadVertices;
            vars.segment(_variablesIdx["vertices"] + 3*_quadTopol.getNumVertices(), vars.size() - (3*_quadTopol.getNumVertices() + _variablesIdx["vertices"])) = centroidVars.segment(_variablesIdx["vertices"] + 3*_centroidTopol.getNumVertices(), centroidVars.size() - (_variablesIdx["vertices"] + 3*_centroidTopol.getNumVertices()));   
            VectorType constraint;
            _constraint.apply(vars, constraint);

            // Calculate the constraint
            MatrixType constraintGrad;
            _DConstraint.apply(vars, constraintGrad);
            
            VectorType constraintGrad_sqrd_1 = (constraintGrad.transpose()*_W) * constraint;
            VectorType constraintGrad_sqrd_2 = ((_W*constraintGrad).transpose()*constraint).transpose();
            VectorType constraintGrad_sqrd = constraintGrad_sqrd_1 + constraintGrad_sqrd_2;

            // Now, need to convert the first quad vertex components of constraintGrad_sqrd to centroid components

            VectorType constraintGrad_sqrd_quad_vertices = constraintGrad_sqrd.segment(_variablesIdx["vertices"], 3*_quadTopol.getNumVertices());
            VectorType constraintGrad_sqrd_centroid = VectorType::Zero(3*_centroidTopol.getNumVertices());
            _coordConverter.convertQuadToCentroid(constraintGrad_sqrd_quad_vertices, constraintGrad_sqrd_centroid);
            VectorType final_constraints(centroidVars.size());
            final_constraints.segment(0, _variablesIdx["vertices"]) = constraintGrad_sqrd.segment(0, _variablesIdx["vertices"]);
            final_constraints.segment(_variablesIdx["vertices"], 3*_centroidTopol.getNumVertices()) = constraintGrad_sqrd_centroid;
            final_constraints.segment(_variablesIdx["vertices"] + 3*_centroidTopol.getNumVertices(), final_constraints.size() - (_variablesIdx["vertices"] + 3*_centroidTopol.getNumVertices())) = constraintGrad_sqrd.segment(_variablesIdx["vertices"] + 3*_quadTopol.getNumVertices(), constraintGrad_sqrd.size() - (_variablesIdx["vertices"] + 3*_quadTopol.getNumVertices()));

            VectorType final_elastic_energies = VectorType::Zero(centroidVars.size());
            final_elastic_energies.segment(_variablesIdx["vertices"], 3*_centroidTopol.getNumVertices()) = DelasticCentroidEnergy;

            // Now, transfer this to the component major coordinates
            //std::cout<<"factors elasticity dev: "<<_factors_elasticity_dev[0]<<"; "<<_factors_elasticity_dev[1]<<std::endl;
            //std::cout<<"norm1: "<<_factors_elasticity_dev[0]*final_elastic_energies.norm()<<"; norm2: "<<_factors_elasticity_dev[1]*final_constraints.norm()<<std::endl;
            gradient = _factors_elasticity_dev[0]*final_elastic_energies + _factors_elasticity_dev[1]*final_constraints;
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

    QuadMeshTopologySaver _quadTopol;
    MeshTopologySaver _centroidTopol;

    VectorType _quadBaseGeometry;
    VectorType _centroidBaseGeometry;

    VectorType _factors_elasticity_dev;
    VectorType _factors_mem_bend;

    NonlinearMembraneEnergy<ConfiguratorType> _E_mem;
    SimpleBendingEnergy<ConfiguratorType> _E_bend;

    ConstraintSqrdReduced<ConfiguratorType> _constraint;
    EdgeLengthQuadEnergy<ConfiguratorType> _constraintEdgeLength;

    CoordConverter<ConfiguratorType> _coordConverter;

public:
    QuadElasticEnergy(const QuadMeshTopologySaver &quadTopol,
                        const MeshTopologySaver &centroidTopol,
                        const VectorType &quadBaseGeometry,
                        const VectorType &centroidBaseGeometry,
                        VectorType factors_elasticity_dev,
                        VectorType factors_mem_bend)
        : _quadTopol(quadTopol),
        _centroidTopol(centroidTopol),
        _quadBaseGeometry(quadBaseGeometry),
        _centroidBaseGeometry(centroidBaseGeometry),
        _factors_elasticity_dev(factors_elasticity_dev),
        _factors_mem_bend(factors_mem_bend),
        _E_bend(centroidTopol, centroidBaseGeometry, true),
        _E_mem(centroidTopol, centroidBaseGeometry, true),
        _constraint(quadTopol),
        _constraintEdgeLength(quadTopol,quadBaseGeometry),
        _coordConverter(quadTopol, centroidTopol, quadBaseGeometry, centroidBaseGeometry)
        {}

    void apply(const VectorType &centroidGeom, RealType &energy) const override
    {
        // First, evaluate the elastic energies on the centroid vertices
        RealType membraneEnergy = 0.0;
        RealType bendingEnergy = 0.0;
        
        if(!(_factors_elasticity_dev[0] == 0.0))
        {
            if(!(_factors_mem_bend[0] == 0.0)){
                _E_mem.apply(centroidGeom, membraneEnergy);
            }
            if(!(_factors_mem_bend[1] == 0.0)){
                _E_bend.apply(centroidGeom, bendingEnergy);
            }
        }
        RealType elasticEnergies = _factors_mem_bend[0]*membraneEnergy + _factors_mem_bend[1]*bendingEnergy;
        //std::cout<<"elastic energies: "<<elasticEnergies<<std::endl;
        //std::cout<<"Factors mem bend: "<<_factors_mem_bend[0]<<"; "<<_factors_mem_bend[1]<<std::endl;

        VectorType quadGeom(3*_quadTopol.getNumVertices());
        _coordConverter.convertCentroidToQuad(centroidGeom, quadGeom);

        RealType developabilityEnergy = 0.0;
        if(!(_factors_elasticity_dev[1] == 0.0))
        {
            _constraint.apply(quadGeom, developabilityEnergy);
        }
        RealType constraintEdge = 0.0;
        _constraintEdgeLength.apply(quadGeom, constraintEdge);
        //std::cout<<"developability energy: "<<developabilityEnergy<<std::endl;
        //std::cout<<"Factors elasticity dev: "<<_factors_elasticity_dev[0]<<"; "<<_factors_elasticity_dev[1]<<std::endl;
        energy = _factors_elasticity_dev[0]*elasticEnergies + _factors_elasticity_dev[1] + developabilityEnergy + 0.0*constraintEdge;
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
    MeshTopologySaver _centroidTopol;

    VectorType _quadBaseGeometry;
    VectorType _centroidBaseGeometry;
    
    VectorType _factors_elasticity_dev;
    VectorType _factors_mem_bend;

    NonlinearMembraneGradientDef<ConfiguratorType> _DE_mem;
    SimpleBendingGradientDef<ConfiguratorType> _DE_bend;
    
    ConstraintSqrdReducedGradient<ConfiguratorType> _DConstraint;
    EdgeLengthQuadEnergyGradient<ConfiguratorType> _DconstraintEdgeLength;

    CoordConverter<ConfiguratorType> _coordConverter;

public:
    QuadElasticEnergyGradient(const QuadMeshTopologySaver &quadTopol,
                                const MeshTopologySaver &centroidTopol,
                                const VectorType &quadBaseGeometry,
                                const VectorType &centroidBaseGeometry,
                                VectorType factors_elasticity_dev,
                                VectorType factors_mem_bend)
        : _quadTopol(quadTopol),
          _centroidTopol(centroidTopol),
          _quadBaseGeometry(quadBaseGeometry),
          _centroidBaseGeometry(centroidBaseGeometry),
          _factors_elasticity_dev(factors_elasticity_dev),
          _factors_mem_bend(factors_mem_bend),
          _DE_bend(centroidTopol, centroidBaseGeometry, true),
          _DE_mem(centroidTopol, centroidBaseGeometry, true),
          _DConstraint(quadTopol),
          _DconstraintEdgeLength(quadTopol, quadBaseGeometry),
            _coordConverter(quadTopol, centroidTopol, quadBaseGeometry, centroidBaseGeometry)
          {}

    void apply(const VectorType &centroidGeom, VectorType &gradient) const override
    {
        // First, evaluate the elastic energies on the centroid vertices
        VectorType DmembraneEnergy = VectorType::Zero(centroidGeom.size());
        VectorType DbendingEnergy = VectorType::Zero(centroidGeom.size());
        
        if(!(_factors_elasticity_dev[0] == 0.0))
        {
            if(!(_factors_mem_bend[0] == 0.0)){
                _DE_mem.apply(centroidGeom, DmembraneEnergy);
            }
            if(!(_factors_mem_bend[1] == 0.0)){
                _DE_bend.apply(centroidGeom, DbendingEnergy);
            }
        }
        VectorType DelasticEnergies = _factors_mem_bend[0]*DmembraneEnergy + _factors_mem_bend[1]*DbendingEnergy;

        VectorType Dconstraint = VectorType::Zero(3*_quadTopol.getNumVertices());
        VectorType quadGeom(3*_quadTopol.getNumVertices());
        _coordConverter.convertCentroidToQuad(centroidGeom, quadGeom);

        if(!(_factors_elasticity_dev[1] == 0.0))
        {
            _DConstraint.apply(quadGeom, Dconstraint);
        }

        VectorType DconstraintEdgeLength = VectorType::Zero(3*_quadTopol.getNumVertices());
        _DconstraintEdgeLength.apply(quadGeom, DconstraintEdgeLength);

        //std::cout<<"Dconstraint norm: "<<Dconstraint.norm()<<std::endl;

        // Convert the derivative vector from quad to centroid coordinates
        VectorType Dconstraint_centroid(3*_centroidTopol.getNumVertices());
        _coordConverter.convertQuadToCentroid(Dconstraint, Dconstraint_centroid);

        VectorType DconstraintEdgeLength_centroid(3*_centroidTopol.getNumVertices());
        _coordConverter.convertQuadToCentroid(DconstraintEdgeLength, DconstraintEdgeLength_centroid);

        gradient = _factors_elasticity_dev[0]*DelasticEnergies + _factors_elasticity_dev[1]*Dconstraint_centroid + 0.0*DconstraintEdgeLength_centroid;
    }
};