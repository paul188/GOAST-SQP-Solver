#include "BoundaryDOFS.h"
#include <goast/Core.h>
#include "SparseMat.h"

typedef typename DefaultConfigurator::SparseMatrixType MatrixType;
typedef typename DefaultConfigurator::FullMatrixType FullMatrixType;

int main(){
    

    /*
    MatrixType A(6,6);
    A.setZero();

    MatrixType B(3,3);
    B.setIdentity();

    assignSparseBlockInplace(A,B,0,0);

    MatrixType C(1,3);
    typedef Eigen::Triplet<double> Tri;
    std::vector<Tri> tripletVec;
    tripletVec.push_back(Tri(0,0,2));
    tripletVec.push_back(Tri(0,1,3));
    tripletVec.push_back(Tri(0,2,4));
    C.setFromTriplets(tripletVec.begin(), tripletVec.end());

    assignSparseBlockInplace(A,C,0,3);

    FullMatrixType A_temp = A.toDense();

    for(int i = 0; i < A_temp.rows(); i++){
        for(int j = 0; j < A_temp.cols(); j++){

            std::cout<<A_temp(i,j)<<" ";
        }
        std::cout<<std::endl;
    } 
    */

    return 0;
}