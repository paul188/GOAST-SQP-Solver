#ifndef BOUNDARY_DOFS_H
#define BOUNDARY_DOFS_H

#include <goast/Core.h>
#include "SparseMat.h"

// We dont want any translations of the deformed geometry at fixed boundary DOFs
// To account for that, we need to reduce matrices and vectors to the free DOFs
// And we also need functions to transform back in order to be able to add to the deformed geometry
/// @tparam ConfiguratorType 
template <typename ConfiguratorType>
class BoundaryDOFS
{
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;
        typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

        // size_t vertexDOFs: Number of all vertex DOFs with the boundary DOFs (so not only effective DOFs)
        BoundaryDOFS(const std::vector<int> &bdryMaskOpt, const size_t &vertexDOFs, const size_t &foldDofs){
            _bdryMaskOpt = bdryMaskOpt;
            _bdryMaskSize = bdryMaskOpt.size();
            _vertexDOFs = vertexDOFs;
            _foldDOFs = foldDofs;

            size_t numBdryDOFs = _bdryMaskOpt.size();
            // We will create two permutation matrices.
            // The first one will be used to reduce the vertexDOFs on vectors/matrices containing only vertexDOFs
            // The second one will be used to reduce the vertexDOFs on vectors/matrices containing vertexDOFs and foldDOFs

            // Now, calculate the first permutation matrix -> we want to move the DOFs of bdryMaskOpt to the beginning
            Eigen::VectorXi permVec(vertexDOFs);

            // First the boundary DOFs go to the beginning in both cases
            for(int i = 0; i < numBdryDOFs; i++){
                permVec[i] =  bdryMaskOpt[i];
            }

            // Now, insert the vertex dofs without the boundary DOFs
            size_t counter_permVecPos = 0;
            size_t counter_DOFNum = 0;
            size_t counter_bdryMask = 0;
            while(true)
            {
                if(counter_permVecPos + numBdryDOFs  == vertexDOFs){
                    break;
                }

                // Check if the degree of freedom is a boundary degree of freedom
                if(counter_DOFNum == bdryMaskOpt[counter_bdryMask]){
                    // If yes, do not add to the permutation vector
                    // But increment the counter for the boundary mask
                    // Do not increment counter of permVec
                    counter_DOFNum++;
                    counter_bdryMask++;
                }
                else{
                    permVec[counter_permVecPos + numBdryDOFs] = counter_DOFNum;
                    counter_DOFNum++;
                    counter_permVecPos++;
                }
            }

            // Now, initialize the permutation matrix
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> p;
            p.indices() = permVec;

            _perm = p.transpose();


            // ###################### Now second permutation matrix #########################

            size_t totalDOFs = _vertexDOFs + _foldDOFs;
            // Now, calculate the second permutation matrix
            Eigen::VectorXi permVecWithFoldDOFs(totalDOFs);
            // Have to add the number of fold DOFs to the boundary mask
            for(int i = 0; i < numBdryDOFs; i++){
                permVecWithFoldDOFs[i] = _foldDOFs + bdryMaskOpt[i];
            }
            
            // Just take all the foldDOFs identically afterwards
            for(int i = 0; i < _foldDOFs; i++){
                permVecWithFoldDOFs[i + numBdryDOFs] = i;
            }

            // Now, insert the vertex dofs without the boundary DOFs
            counter_permVecPos = 0;
            counter_DOFNum = 0;
            counter_bdryMask = 0;
            while(true)
            {
                if(counter_permVecPos + numBdryDOFs + _foldDOFs == totalDOFs){
                    break;
                }

                // Check if the degree of freedom is a boundary degree of freedom
                if(counter_DOFNum == bdryMaskOpt[counter_bdryMask]){
                    // If yes, do not add to the permutation vector
                    // But increment the counter for the boundary mask
                    // Do not increment counter of permVec
                    counter_DOFNum++;
                    counter_bdryMask++;
                }
                else{
                    permVecWithFoldDOFs[counter_permVecPos + numBdryDOFs + _foldDOFs] = counter_DOFNum + _foldDOFs;
                    counter_DOFNum++;
                    counter_permVecPos++;
                }
            }

            // Now, initialize the permutation matrix
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> pWithFoldDOFs;
            pWithFoldDOFs.indices() = permVecWithFoldDOFs;

            _permWithFoldDOFs = pWithFoldDOFs.transpose();
        }

        size_t getNumBdryDOFs() const{
            return _bdryMaskOpt.size();
        }

        // TRANSFORMATIONS FOR VECTORS ##############################################################

        void transformWithFoldDofsToReducedSpace(VectorType &vec) const{
            // Permute the entries of the vector such that the boundary DOFs are at the beginning
            size_t k = vec.size();
            assert(k == _vertexDOFs + _foldDOFs);
            vec = _permWithFoldDOFs * vec;
            vec = vec.segment(_bdryMaskSize, k - _bdryMaskSize);
        }

        void InverseTransformWithFoldDofs(VectorType &vec) const{
            // Add the boundary DOFs to the beginning
            VectorType vecFull = VectorType::Zero(_vertexDOFs + _foldDOFs);
            vecFull.segment(_bdryMaskSize, vec.size()) = vec;
            vec = _permWithFoldDOFs.inverse() * vecFull;
        }

        // Transform a vector to the reduced space
        void transformToReducedSpace(VectorType &vec) const{
            // Permute the entries of the vector such that the boundary DOFs are at the beginning
            size_t k = vec.size();
            assert(k == _vertexDOFs);
            vec = _perm * vec;
            vec = vec.segment(_bdryMaskSize, k - _bdryMaskSize);
        }

        void inverseTransform(VectorType &vec) const{
            // Add the boundary DOFs to the beginning
            VectorType vecFull = VectorType::Zero(_vertexDOFs);
            vecFull.segment(_bdryMaskSize, vec.size()) = vec;
            vec = _perm.inverse() * vecFull;
        }

        // TRANSFORMATIONS FOR DENSE MATRICES ###########################################################

        void transformRowToReducedSpace(FullMatrixType &mat) const{
            assert(mat.rows() == _vertexDOFs);
            mat = _perm * mat;
            mat = mat.block(_bdryMaskSize,0,_vertexDOFs - _bdryMaskSize, mat.cols());
        }

        void inverseTransformRow(FullMatrixType &mat) const{
            FullMatrixType matFull = FullMatrixType::Zero(_vertexDOFs, mat.cols());
            matFull.block(_bdryMaskSize,0,_vertexDOFs - _bdryMaskSize, mat.cols()) = mat;
            mat = _perm.inverse() * matFull;
        }

        void transformColToReducedSpace(FullMatrixType &mat) const{
            assert(mat.cols() == _vertexDOFs);
            mat = mat * _perm.transpose();
            mat = mat.block(0,_bdryMaskSize, mat.rows(), _vertexDOFs - _bdryMaskSize);
        }

        void inverseTransformCol(FullMatrixType &mat) const{
            assert(mat.cols() == _vertexDOFs - _bdryMaskSize);
            FullMatrixType matFull = FullMatrixType::Zero(mat.rows(), _vertexDOFs);
            matFull.block(0,_bdryMaskSize,mat.rows(), _vertexDOFs - _bdryMaskSize) = mat;
            // inverse of tranposed of a permutation matrix is the permutation matrix itself
            mat = matFull * _perm;
        }

        // TRANSFORMATIONS FOR SPARSE MATRICES #########################################################

        void transformRowToReducedSpace(MatrixType &mat) const{
            assert(mat.rows() == _vertexDOFs);
            mat = _perm * mat;
            mat = extractSparseBlock(mat, _bdryMaskSize, 0, _vertexDOFs - _bdryMaskSize, mat.cols());
        }

        void inverseTransformRow(MatrixType &mat) const{
            assert(mat.rows() == _vertexDOFs - _bdryMaskSize);
            MatrixType matFull(_vertexDOFs, mat.cols());
            matFull.setZero();
            matFull = assignSparseBlock(matFull, mat, _bdryMaskSize, 0);
            mat = _perm.inverse() * matFull;
        }

        void transformColToReducedSpace(MatrixType &mat) const{
            assert(mat.cols() == _vertexDOFs);
            mat = mat * _perm.transpose();
            mat = extractSparseBlock(mat, 0, _bdryMaskSize, mat.rows(), _vertexDOFs - _bdryMaskSize);
        }

        void inverseTransformCol(MatrixType &mat) const{
            MatrixType matFull(mat.rows(), _vertexDOFs);
            matFull.setZero();
            matFull = assignSparseBlock(matFull, mat, 0, _bdryMaskSize);
            mat = matFull * _perm;
        }

        void transformRowColToReducedSpace(MatrixType &mat) const{
            std::cout<<_vertexDOFs<<std::endl;
            std::cout<<mat.rows()<<" "<<mat.cols()<<std::endl;
            assert(mat.rows() == _vertexDOFs && mat.cols() == _vertexDOFs);
            transformRowToReducedSpace(mat);
            transformColToReducedSpace(mat);
        }

        void inverseTransformRowCol(MatrixType &mat) const{
            inverseTransformRow(mat);
            inverseTransformCol(mat);
        }

    protected:
        std::vector<int> _bdryMaskOpt;
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> _perm;
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> _permWithFoldDOFs;
        size_t _vertexDOFs;
        size_t _bdryMaskSize;
        size_t _foldDOFs;
};

#endif