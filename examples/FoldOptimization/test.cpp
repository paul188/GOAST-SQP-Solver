#include "BoundaryDOFS.h"
#include <goast/Core.h>
#include "SparseMat.h"

typedef typename DefaultConfigurator::SparseMatrixType MatrixType;
typedef typename DefaultConfigurator::FullMatrixType FullMatrixType;

int main(){
    std::vector<int> vertexDOFS = {0,1,2,3,4};
    std::vector<int> bdryMask ={2,4};
    size_t foldDOFS = 1;
    BoundaryDOFS<DefaultConfigurator> bdryDOFS(bdryMask, vertexDOFS.size(), foldDOFS);

    MatrixType A(3,3);

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.push_back(Eigen::Triplet<double>(0,0,1));
    tripletList.push_back(Eigen::Triplet<double>(0,2,2));
    tripletList.push_back(Eigen::Triplet<double>(1,1,3));
    tripletList.push_back(Eigen::Triplet<double>(2,0,4));
    tripletList.push_back(Eigen::Triplet<double>(2,2,5));
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    MatrixType B(2,2);
    std::vector<Eigen::Triplet<double>> tripletListB;
    tripletListB.push_back(Eigen::Triplet<double>(0,0,1));
    tripletListB.push_back(Eigen::Triplet<double>(0,1,2));
    tripletListB.push_back(Eigen::Triplet<double>(1,0,3));
    tripletListB.push_back(Eigen::Triplet<double>(1,1,12));
    B.setFromTriplets(tripletListB.begin(), tripletListB.end());

    MatrixType Anew = assignSparseBlock(A, B, 0, 0);

    FullMatrixType A_temp = Anew.toDense();

    for(int i = 0; i < A_temp.rows(); i++){
        for(int j = 0; j < A_temp.cols(); j++){
            std::cout << A_temp(i,j) << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}