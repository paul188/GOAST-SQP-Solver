#pragma once

#include "Constraints.h"

template<typename ConfiguratorType>
class BoundaryDOFS_quad
{
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;
        typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    private:
        std::vector<int> _extendedBdryMaskOpt; // The boundary mask
        const size_t _extendedBdryMaskSize; // The size of the boundary mask
        const size_t _num_vertices; // The number of vertex DOFs
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> _perm;
        VarsIdx<DefaultConfigurator> _varsIdx;

    public:
        BoundaryDOFS_quad(const std::vector<int> &bdryMaskOpt, const size_t num_vertices, 
                          const VarsIdx<DefaultConfigurator> &varsIdx)
                          : _num_vertices(num_vertices), _varsIdx(varsIdx), _extendedBdryMaskSize(3*bdryMaskOpt.size())
        {

            for(int i = 0; i < bdryMaskOpt.size(); ++i)
            {
                _extendedBdryMaskOpt.push_back(3*bdryMaskOpt[i]);
                _extendedBdryMaskOpt.push_back(3*bdryMaskOpt[i] + 1);
                _extendedBdryMaskOpt.push_back(3*bdryMaskOpt[i] + 2);
            }

            const size_t size_permVec = 3*_num_vertices;
            Eigen::VectorXi permVec(size_permVec);
            permVec = Eigen::VectorXi::LinSpaced(size_permVec,0, size_permVec - 1);
               
            // put the boundary indices to the beginning of the vertices section of vars
            for(int i = 0; i < _extendedBdryMaskOpt.size(); ++i)
            {
                permVec[i] = _extendedBdryMaskOpt[i];
            }

            size_t bdryCounter = 0;
            size_t dofCounter = 0;
            size_t posCounter = _extendedBdryMaskSize;
            while(true)
            {
                if(dofCounter == 3*_num_vertices)
                {
                    break;
                }
                if(dofCounter == _extendedBdryMaskOpt[bdryCounter])
                {
                    bdryCounter++;
                    dofCounter++;
                    continue;
                }
                else{
                    permVec[posCounter] = dofCounter;
                    posCounter++;
                    dofCounter++;
                }
            }

            // Now, initialize the permutation matrix
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> p;
            p.indices() = permVec;

            _perm = p.transpose();
        }

        void transformToReducedSpace(const VectorType &vars, VectorType &Dest)
        {
            assert(vars.size() == _varsIdx["num_dofs"]);
            Dest = VectorType::Zero(_varsIdx["num_dofs"] - _extendedBdryMaskSize);
            // First permute, such that the boundary values are at the beginning  
            VectorType temp = _perm*vars.segment(_varsIdx["vertices"], 3*_num_vertices);
            VectorType temp2 = temp.segment(_extendedBdryMaskSize, temp.size() - _extendedBdryMaskSize);

            size_t begin_segment_0 = 0;
            size_t end_segment_0 = _varsIdx["vertices"];
            size_t begin_segment_1 = _varsIdx["vertices"];
            size_t end_segment_1 = begin_segment_1 + 3*_num_vertices - _extendedBdryMaskSize;
            size_t begin_segment_2 = end_segment_1;
            size_t end_segment_2 = _varsIdx["num_dofs"] - _extendedBdryMaskSize;

            Dest.segment(begin_segment_0, (end_segment_0 - begin_segment_0)) = vars.segment(0, _varsIdx["vertices"]);
            Dest.segment(begin_segment_1, (end_segment_1 - begin_segment_1)) = temp2.segment(0, (end_segment_1 - begin_segment_1));
            Dest.segment(begin_segment_2, (end_segment_2 - begin_segment_2)) = vars.segment(_varsIdx["vertices"] + 3*_num_vertices, _varsIdx["num_dofs"] - (_varsIdx["vertices"] + 3*_num_vertices));
        }

        void InverseTransform(const VectorType &vars, VectorType &Dest)
        {
            assert(vars.size() == _varsIdx["num_dofs"] - _extendedBdryMaskSize);
            VectorType temp = VectorType::Zero(3*_num_vertices);
            temp.segment(_extendedBdryMaskSize, temp.size() - _extendedBdryMaskSize) = vars.segment(_varsIdx["vertices"], 3*_num_vertices - _extendedBdryMaskSize);
            VectorType temp2 = _perm.inverse()*temp;
            // First permute, such that the nonboundary values are at the beginning  
            Dest = VectorType::Zero(_varsIdx["num_dofs"]);
            // First, fill the first entries until the vertices
            Dest.segment(0, _varsIdx["vertices"]) = vars.segment(0, _varsIdx["vertices"]);
            Dest.segment(_varsIdx["vertices"], 3*_num_vertices) = temp2.segment(0, 3*_num_vertices);
            // Now, fill the vertex values without the bdry values
            Dest.segment(_varsIdx["vertices"] + 3*_num_vertices, _varsIdx["num_dofs"] - (3*_num_vertices + _varsIdx["vertices"])) = vars.segment(_varsIdx["vertices"] + 3*_num_vertices - _extendedBdryMaskSize, vars.size() - (_varsIdx["vertices"] + 3*_num_vertices - _extendedBdryMaskSize));
        }

        void transformColsToReducedSpace(const MatrixType &mat, MatrixType &Dest)
        {
            assert(mat.cols() == _varsIdx["num_dofs"]);
            Dest = MatrixType(mat.rows(), _varsIdx["num_dofs"] - _extendedBdryMaskSize);
            Dest.setZero();
            MatrixType temp = extractSparseBlock(mat, 0, _varsIdx["vertices"], mat.rows(), 3*_num_vertices);
            MatrixType temp2 = temp * _perm.transpose();
            temp = extractSparseBlock(temp2,0, _extendedBdryMaskSize, mat.rows(), 3*_num_vertices - _extendedBdryMaskSize);
            // Now the columns with the boundary vertex entries are at the beginning of the vertex entries
            MatrixType block_1 = extractSparseBlock(mat, 0, 0, mat.rows(), _varsIdx["vertices"]);
            MatrixType block_2 = extractSparseBlock(mat, 0, _varsIdx["vertices"] + 3*_num_vertices, mat.rows(), mat.cols() - (_varsIdx["vertices"] + 3*_num_vertices));
            std::vector<Eigen::Triplet<double>> triplets;
            assignSparseBlockInplace(Dest, block_1, 0, 0, triplets);
            assignSparseBlockInplace(Dest, temp, 0, _varsIdx["vertices"], triplets);
            assignSparseBlockInplace(Dest, block_2, 0, (_varsIdx["vertices"] + 3*_num_vertices - _extendedBdryMaskSize), triplets);
        }

        int getNumReducedDOFs()
        {
            return _varsIdx["num_dofs"] - _extendedBdryMaskSize;
        }

        /*
        void InverseTransformColsToReducedSpace(const MatrixType &mat, MatrixType &Dest)
        {

        }
        */

};