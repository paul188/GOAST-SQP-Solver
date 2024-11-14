#include <goast/Core.h>
#include "SparseMat.h"

// We dont want any translations of the deformed geometry at fixed boundary DOFs
// To account for that, we need to reduce matrices and vectors to the free DOFs
// And we also need functions to transform back in order to be able to add to the deformed geometry
template <typename ConfiguratorType>
class BoundaryDOFS
{
    public:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;
        typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

        // size_t foldDOFS: number of fold DOFs
        BoundaryDOFS(const std::vector<int> &bdryMaskOpt, const size_t &vertexDOFs, size_t foldDOFs){
            _bdryMaskOpt = bdryMaskOpt;

            size_t totalDOFs = vertexDOFs + foldDOFs;
            size_t numBdryDOFs = _bdryMaskOpt.size();
            // Now, calculate the permutation matrix -> we want to move the DOFs of bdryMaskOpt to the beginning
            Eigen::VectorXi permVec(totalDOFs);
            // Have to add the number of fold DOFs to the boundary mask
            for(int i = 0; i < numBdryDOFs; i++){
                permVec[i] = foldDOFs + bdryMaskOpt[i];
            }
            
            // Just take all the foldDOFs identically afterwards
            for(int i = 0; i < foldDOFs; i++){
                permVec[i + numBdryDOFs] = i;
            }

            // Now, insert the vertex dofs without the boundary DOFs
            size_t counter_permVecPos = 0;
            size_t counter_DOFNum = 0;
            size_t counter_bdryMask = 0;
            while(true)
            {
                if(counter_permVecPos + numBdryDOFs + foldDOFs == totalDOFs){
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
                    permVec[counter_permVecPos + numBdryDOFs + foldDOFs] = counter_DOFNum + foldDOFs;
                    counter_DOFNum++;
                    counter_permVecPos++;
                }
            }

            // Now, initialize the permutation matrix
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> p;
            p.indices() = permVec;

            _perm = p.transpose();
        }

        size_t getNumBdryDOFs() const{
            return _bdryMaskOpt.size();
        }

        // TRANSFORMATIONS FOR VECTORS ##############################################################

        // Transform a vector to the reduced space
        void transformToReducedSpace(VectorType &vec) const{
            // Permute the entries of the vector such that the boundary DOFs are at the beginning
            vec = _perm * vec;
            vec = vec.segment(_bdryMaskOpt.size(), vec.size() - _bdryMaskOpt.size());
        }

        void inverseTransform(VectorType &vec) const{
            // Add the boundary DOFs to the beginning
            VectorType vecFull = VectorType::Zero(_perm.indices().size());
            vecFull.segment(_bdryMaskOpt.size(), vec.size()) = vec;
            vec = _perm.inverse() * vecFull;
        }

        // TRANSFORMATIONS FOR DENSE MATRICES ###########################################################

        void transformRowToReducedSpace(FullMatrixType &mat) const{
            mat = _perm * mat;
            mat = mat.block(_bdryMaskOpt.size(),0,mat.rows() - _bdryMaskOpt.size(), mat.cols());
        }

        void inverseTransformRow(FullMatrixType &mat) const{
            FullMatrixType matFull = FullMatrixType::Zero(_perm.indices().size(), mat.cols());
            matFull.block(_bdryMaskOpt.size(),0,mat.rows(), mat.cols()) = mat;
            mat = _perm.inverse() * matFull;
        }

        void transformColToReducedSpace(FullMatrixType &mat) const{
            mat = mat * _perm.transpose();
            mat = mat.block(0,_bdryMaskOpt.size(),mat.rows(), mat.cols() - _bdryMaskOpt.size());
        }

        void inverseTransformCol(FullMatrixType &mat) const{
            FullMatrixType matFull = FullMatrixType::Zero(mat.rows(), _perm.indices().size());
            matFull.block(0,_bdryMaskOpt.size(),mat.rows(), mat.cols()) = mat;
            // inverse of tranposed of a permutation matrix is the permutation matrix itself
            mat = matFull * _perm;
        }

        // TRANSFORMATIONS FOR SPARSE MATRICES #########################################################

        void transformRowToReducedSpace(MatrixType &mat) const{
            mat = _perm * mat;
            mat = extractSparseBlock(mat, _bdryMaskOpt.size(), 0, mat.rows() - _bdryMaskOpt.size(), mat.cols());
        }

        void inverseTransformRow(MatrixType &mat) const{
            MatrixType matFull = MatrixType::Zero(_perm.indices().size(), mat.cols());
            matFull = assignSparseBlock(matFull, mat, _bdryMaskOpt.size(), 0);
            mat = _perm.inverse() * matFull;
        }

        void transformColToReducedSpace(MatrixType &mat) const{
            mat = mat * _perm.transpose();
            mat = extractSparseBlock(mat, 0, _bdryMaskOpt.size(), mat.rows(), mat.cols() - _bdryMaskOpt.size());
        }

        void inverseTransformCol(MatrixType &mat) const{
            MatrixType matFull = MatrixType::Zero(mat.rows(), _perm.indices().size());
            matFull = assignSparseBlock(matFull, mat, 0, _bdryMaskOpt.size());
            mat = matFull * _perm;
        }

    protected:
        std::vector<int> _bdryMaskOpt;
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> _perm;
};