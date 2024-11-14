#ifndef SPARSEMAT_H
#define SPARSEMAT_H

#include <goast/Core.h>
#include <iostream>

typedef typename DefaultConfigurator::SparseMatrixType MatrixType;

void Sparsity( const MatrixType &matrix ){
    int numEntries = 0;
    for( int i = 0; i < matrix.rows(); i++ ){
        for( MatrixType::InnerIterator it(matrix,i); it; ++it ){
            numEntries++;
        }
    }
    std::cerr << "Relative number of nonzero entries: " << numEntries / ((int)(matrix.cols()*matrix.rows())) << std::endl;
}

// method to extract a block from a sparse matrix
typedef Eigen::Triplet<double> Tri;
MatrixType extractSparseBlock(MatrixType M, int ibegin, int jbegin, int icount, int jcount) {
    assert(ibegin + icount <= M.rows());
    assert(jbegin + jcount <= M.cols());
    
    int Mj, Mi, i, j, currOuterIndex, nextOuterIndex;
    std::vector<Tri> tripletList;
    tripletList.reserve(M.nonZeros());

    for (i = 0; i < icount; i++) {
        Mi = i + ibegin;
        currOuterIndex = M.outerIndexPtr()[Mi];
        nextOuterIndex = M.outerIndexPtr()[Mi + 1];

        for (int a = currOuterIndex; a < nextOuterIndex; a++) {
            Mj = M.innerIndexPtr()[a];

            if (Mj < jbegin) continue;
            if (Mj >= jbegin + jcount) break;

            j = Mj - jbegin;
            tripletList.push_back(Tri(i, j, M.valuePtr()[a]));
        }
    }

    MatrixType matS(icount, jcount);
    matS.setFromTriplets(tripletList.begin(), tripletList.end());
    return matS;
}

MatrixType assignSparseBlock(MatrixType M, MatrixType Mblock, int ibegin, int jbegin) {
    assert(ibegin + Mblock.rows() <= M.rows());
    assert(jbegin + Mblock.cols() <= M.cols());

    M.sortInnerIndices(0, M.outerSize());
    Mblock.sortInnerIndices(0, Mblock.outerSize());

    int Mj,Mi, i, j, currOuterIndex, nextOuterIndex, currOuterIndexBlock, nextOuterIndexBlock;
    std::vector<Tri> tripletList;
    tripletList.reserve(M.nonZeros() + Mblock.nonZeros());

    // Assign the rows above the block matrix
    for(int i = 0; i < M.rows() && i < ibegin; i++)
    {
        currOuterIndex = M.outerIndexPtr()[i];
        nextOuterIndex = M.outerIndexPtr()[i + 1];

        for(int a = currOuterIndex; a < nextOuterIndex ; a++){
            Mj = M.innerIndexPtr()[a];
            tripletList.push_back(Tri(i, Mj, M.valuePtr()[a]));
        }
    }

    // Assign the rows below the block matrix
    for(int i = ibegin + Mblock.rows(); i < M.rows(); i++)
    {
        currOuterIndex = M.outerIndexPtr()[i];
        nextOuterIndex = M.outerIndexPtr()[i + 1];

        for(int a = currOuterIndex; a < nextOuterIndex ; a++){
            Mj = M.innerIndexPtr()[a];
            tripletList.push_back(Tri(i, Mj, M.valuePtr()[a]));
        }
    }

    // Assign the rows mixed with the block matrix
    for(int i = ibegin; i < ibegin + Mblock.rows(); i++)
    {
        currOuterIndex = M.outerIndexPtr()[i];
        nextOuterIndex = M.outerIndexPtr()[i + 1];

        currOuterIndexBlock = Mblock.outerIndexPtr()[i - ibegin];
        nextOuterIndexBlock = Mblock.outerIndexPtr()[i - ibegin + 1];

        // Assign the columns before the block matrix but in the same rows
        for(int a = currOuterIndex; M.innerIndexPtr()[a] < Mblock.innerIndexPtr()[currOuterIndexBlock] && a < nextOuterIndex; a++)
        {
            Mj = M.innerIndexPtr()[j];
            tripletList.push_back(Tri(i, Mj, M.valuePtr()[a]));
        }

        // Assign the block matrix
        for(int a = currOuterIndexBlock; a < nextOuterIndexBlock; a++)
        {
            Mj = Mblock.innerIndexPtr()[a];
            tripletList.push_back(Tri(i, Mj+jbegin, Mblock.valuePtr()[a]));
        }

        // Assign the columns after the block matrix but in the same rows
        for(int a = nextOuterIndexBlock - 1; M.innerIndexPtr()[a] > Mblock.innerIndexPtr()[nextOuterIndexBlock - 1] && a >= currOuterIndex; a++)
        {
            Mj = M.innerIndexPtr()[a];
            tripletList.push_back(Tri(i, Mj, M.valuePtr()[a]));
        }
        for(int i = 0; i < tripletList.size(); i++){
            std::cout<<"("<<tripletList[i].row()<<","<<tripletList[i].col()<<","<<tripletList[i].value()<<")"<<std::endl;
        }
    }
    MatrixType Mnew(M.rows(), M.cols());
    Mnew.setFromTriplets(tripletList.begin(), tripletList.end());
    return Mnew;
}

#endif // SPARSEMAT_H