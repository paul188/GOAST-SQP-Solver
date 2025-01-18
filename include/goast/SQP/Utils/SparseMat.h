#ifndef SPARSEMAT_H
#define SPARSEMAT_H

#include <goast/Core.h>
#include <iostream>

typedef typename DefaultConfigurator::SparseMatrixType MatrixType;
typedef typename DefaultConfigurator::VectorType VectorType;
typedef typename DefaultConfigurator::FullMatrixType FullMatrixType;

void Sparsity( const MatrixType &matrix ){
    int numEntries = 0;
    for( int i = 0; i < matrix.rows(); i++ ){
        for( MatrixType::InnerIterator it(matrix,i); it; ++it ){
            numEntries++;
        }
    }
    std::cerr << "Relative number of nonzero entries: " << numEntries / ((int)(matrix.cols()*matrix.rows())) << std::endl;
}

void computeDiagonalPreconditioner(const MatrixType &A, Eigen::DiagonalMatrix<double, Eigen::Dynamic> &D, Eigen::DiagonalMatrix<double, Eigen::Dynamic> &Dinv){
    D.resize(A.rows());
    Dinv.resize(A.rows());
    for(int i = 0; i < A.rows(); i++){
        D.diagonal()[i] = 1.0/(A.row(i).norm() + 1e-6);
        Dinv.diagonal()[i] = (A.row(i).norm() + 1e-6);
    }
}

void calculateIdfit(size_t rows, size_t cols, FullMatrixType &Dest)
{
    if(Dest.rows() != rows || Dest.cols() != cols)
    {
        Dest.resize(rows, cols);
    }

    Dest.setZero();

    if(rows > cols)
    {
        size_t num_submatrices = rows / cols;
        size_t rest = rows % cols;
        for(int i = 0; i < num_submatrices; i++){
            Dest.block(i*cols,0, cols, cols) = FullMatrixType::Identity(cols, cols);
        }
        Dest.block(num_submatrices*cols, 0, rest, rest) = FullMatrixType::Identity(rest, rest);
    }
    else if(rows < cols)
    {
        size_t num_submatrices = cols / rows;
        size_t rest = cols % rows;
        for(int i = 0; i < num_submatrices; i++){
            Dest.block(0,i*rows, rows, rows) = FullMatrixType::Identity(rows, rows);
        }
        Dest.block(0, num_submatrices*rows, rest, rest) = FullMatrixType::Identity(rest, rest);
    }
    else
    {
        // rows = cols
        Dest.setIdentity(rows, cols);
    }
}

MatrixType convertVecToSparseMat(const VectorType& vec){
    MatrixType sparseMatrix(vec.size(), 1); // Single column sparse matrix
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i] != 0.0) {
            sparseMatrix.insert(i, 0) = vec[i];
        }
    }

    //sparseMatrix.makeCompressed();
    return sparseMatrix;
}

bool is_posdef(MatrixType H) {
    auto H_dense = H.toDense();
    Eigen::LLT<FullMatrixType> llt(H_dense);
    if (llt.info() == Eigen::NumericalIssue) {
        return false;
    }
    return true;
}

inline bool replaceNanWithZero(MatrixType& x) {
    Eigen::Map<Eigen::Array<typename MatrixType::Scalar, Eigen::Dynamic, 1>> values(x.valuePtr(), x.nonZeros());
    if (values.hasNaN()) {
        std::cout << "NaN values detected in the matrix" << std::endl;

        for (int k = 0; k < x.outerSize(); ++k) {
            for (MatrixType::InnerIterator it(x, k); it; ++it) {
                if (std::isnan(it.value())) {
                    it.valueRef() = 0.0;  // Replace NaN with 0
                }
            }
        }
        
        return true;
    }
    return false;
}

// method to extract a block from a sparse matrix
typedef Eigen::Triplet<double> Tri;
MatrixType extractSparseBlock(MatrixType &M, int ibegin, int jbegin, int icount, int jcount) {
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
    //matS.makeCompressed();
    return matS;
}

/*
MatrixType assignSparseBlock(MatrixType &M, MatrixType &Mblock, int ibegin, int jbegin) {
    assert(ibegin + Mblock.rows() <= M.rows());
    assert(jbegin + Mblock.cols() <= M.cols());

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
            std::cout<<"("<<i<<","<<Mj<<", "<< M.valuePtr()[a]<<")"<<std::endl;
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
            std::cout<<"("<<i<<","<<Mj<<", "<< M.valuePtr()[a]<<")"<<std::endl;
        }
    }

    // Assign the rows mixed with the block matrix
    for(int i = ibegin; i < ibegin + Mblock.rows(); i++)
    {
        currOuterIndex = M.outerIndexPtr()[i];
        nextOuterIndex = M.outerIndexPtr()[i + 1];

        currOuterIndexBlock = Mblock.outerIndexPtr()[i - ibegin];
        nextOuterIndexBlock = Mblock.outerIndexPtr()[i - ibegin + 1];

        if(!(currOuterIndex == nextOuterIndex)) // Check if there are even any nonzero elements in the big matrix -> if not innerIndexPtr will be invalid
        {
            // Assign the columns before the block matrix but in the same rows
            for(int a = currOuterIndex; M.innerIndexPtr()[a] < jbegin && a < nextOuterIndex; a++)
            {
                Mj = M.innerIndexPtr()[j];
                tripletList.push_back(Tri(i, Mj, M.valuePtr()[a]));
                std::cout<<"("<<i<<","<<Mj<<", "<< M.valuePtr()[a]<<")"<<std::endl;
            }
        }

        // Assign the block matrix
        for(int a = currOuterIndexBlock; a < nextOuterIndexBlock; a++)
        {
            Mj = Mblock.innerIndexPtr()[a];
            tripletList.push_back(Tri(i, Mj+jbegin, Mblock.valuePtr()[a]));
            std::cout<<"("<<i<<","<<Mj+jbegin<<", "<< Mblock.valuePtr()[a]<<")"<<std::endl;
        }

        if(!(currOuterIndex == nextOuterIndex)){ // Check if row of larger matrix is zero
            
            // Assign the columns after the block matrix but in the same rows
            for(int a = nextOuterIndex - 1; M.innerIndexPtr()[a] > jbegin + Mblock.cols() && a >= currOuterIndex; a--)
            {
                Mj = M.innerIndexPtr()[a];
                tripletList.push_back(Tri(i, Mj, M.valuePtr()[a]));
                std::cout<<"("<<i<<","<<Mj<<", "<< M.valuePtr()[a]<<")"<<std::endl;
            }
        }
    }
    MatrixType Mnew(M.rows(), M.cols());
    Mnew.setFromTriplets(tripletList.begin(), tripletList.end());
    return Mnew;
}

// more efficient than assignSparseBlock
void assignSparseBlockInplace(MatrixType &M, MatrixType &Mblock, int ibegin, int jbegin){
    assert(ibegin + Mblock.rows() <= M.rows());
    assert(jbegin + Mblock.cols() <= M.cols());

    int Mj,Mi, i, currOuterIndex, nextOuterIndex, currOuterIndexBlock, nextOuterIndexBlock;
    std::vector<Tri> tripletList;
    tripletList.reserve(M.nonZeros() + Mblock.nonZeros());

    if(ibegin > 0)
    {
        // Assign the rows above the block matrix
        for(int i = 0; i < M.rows() && i < ibegin; i++)
        {
            currOuterIndex = M.outerIndexPtr()[i];
            nextOuterIndex = M.outerIndexPtr()[i + 1];

            for(int a = currOuterIndex; a < nextOuterIndex ; a++){
                Mj = M.innerIndexPtr()[a];
                tripletList.push_back(Tri(i, Mj, M.valuePtr()[a]));
                if (i < 0 || i >= M.rows() || Mj < 0 || Mj >= M.cols() ) {
                        std::cerr << "Invalid triplet1: (" << i << ", " << Mj <<  ")" << std::endl;
                        std::cerr <<"Current outer Index: "<<currOuterIndex<<std::endl;
                        std::cerr<<"Next outer index: "<<nextOuterIndex<<std::endl;
                    abort();
                }
            }
        }
    }

    if(ibegin + Mblock.rows() < M.rows())
    {
        // Assign the rows below the block matrix
        for(int i = ibegin + Mblock.rows(); i < M.rows(); i++)
        {
            currOuterIndex = M.outerIndexPtr()[i];
            nextOuterIndex = M.outerIndexPtr()[i + 1];

            for(int a = currOuterIndex; a < nextOuterIndex ; a++){
                Mj = M.innerIndexPtr()[a];
                tripletList.push_back(Tri(i, Mj, M.valuePtr()[a]));
                if (i < 0 || i >= M.rows() || Mj < 0 || Mj >= M.cols()) {
                        std::cerr << "Invalid triplet2: (" << i << ", " << Mj  << ")" << std::endl;
                    abort();
                }
            }
        }
    }

    // Assign the rows mixed with the block matrix
    for(int i = ibegin; i < ibegin + Mblock.rows(); i++)
    {
        std::cout<<"BOUND FOR I: "<<ibegin + Mblock.rows()<<std::endl;
        currOuterIndex = M.outerIndexPtr()[i];
        nextOuterIndex = M.outerIndexPtr()[i + 1];

        currOuterIndexBlock = Mblock.outerIndexPtr()[i - ibegin];
        nextOuterIndexBlock = Mblock.outerIndexPtr()[i - ibegin + 1];

        if(!(currOuterIndex == nextOuterIndex)) // Check if there are even any nonzero elements in the big matrix -> if not innerIndexPtr will be invalid
        {
            // Assign the columns before the block matrix but in the same rows
            for(int a = currOuterIndex; a < nextOuterIndex && M.innerIndexPtr()[a] < jbegin ; a++)
            {
                Mj = M.innerIndexPtr()[a];
                tripletList.push_back(Tri(i, Mj, M.valuePtr()[a]));
                if (i < 0 || i >= M.rows() || Mj < 0 || Mj >= M.cols() ) {
                    std::cerr << "Invalid triplet3: (" << i << ", " << Mj << ", " << Mj + jbegin << ")" << std::endl;
                abort();
                }
            }
        }

        // Assign the block matrix
        for(int a = currOuterIndexBlock; a < nextOuterIndexBlock; a++)
        {
            Mj = Mblock.innerIndexPtr()[a];
            tripletList.push_back(Tri(i, Mj+jbegin, Mblock.valuePtr()[a]));
            if (i < 0 || i >= M.rows() || Mj + jbegin >= M.cols()) {
                    std::cerr << "Invalid triplet4: (" << i << ", " << Mj + jbegin << ")" << std::endl;
                abort();
            }
        }

        if(!(currOuterIndex == nextOuterIndex)){ // Check if row of larger matrix is zero
            
            // Assign the columns after the block matrix but in the same rows
            for(int a = nextOuterIndex - 1; a >= currOuterIndex && M.innerIndexPtr()[a] > jbegin + Mblock.cols(); a--)
            {
                Mj = M.innerIndexPtr()[a];
                tripletList.push_back(Tri(i, Mj, M.valuePtr()[a]));
                if (i < 0 || i >= M.rows() || Mj < 0 || Mj >= M.cols()) {
                    std::cerr << "Invalid triplet5: (" << i << ", " << Mj << ", " << Mj + jbegin << ")" << std::endl;
                abort();
                }
            }
        }
    }

    M.setFromTriplets(tripletList.begin(), tripletList.end());
}
*/

void printVectorToFile(const VectorType& vec, const std::string& filename, int precision = 5) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Write vector size
    file << "Vector size: " << vec.size() << "\n\n";

    // Set formatting
    file << std::fixed << std::setprecision(precision);

    // Print vector
    for (size_t i = 0; i < vec.size(); ++i) {
        file << std::setw(8) << vec[i] << "\n";
    }

    file.close();
    std::cout << "Vector successfully printed to " << filename << "\n";
}



void printSparseMatrix(const Eigen::SparseMatrix<double>& matrix, const std::string& filename, int precision = 3) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Write matrix dimensions
    file << "Matrix dimensions: " << matrix.rows() << " x " << matrix.cols() << "\n\n";

    // Convert sparse matrix to dense for readability
    Eigen::MatrixXd denseMatrix = Eigen::MatrixXd(matrix);

    // Set formatting
    file << std::fixed << std::setprecision(precision);

    // Print matrix
    for (int i = 0; i < denseMatrix.rows(); ++i) {
        for (int j = 0; j < denseMatrix.cols(); ++j) {
            file << std::setw(8) << denseMatrix(i, j) << " ";
        }
        file << "\n";
    }

    file.close();
    std::cout << "Sparse matrix successfully printed to " << filename << "\n";
}

void printMatrixToFile(const Eigen::MatrixXd& matrix, const std::string& filename) {
    // Open the file
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write the matrix to the file
    file << matrix.format(Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n"));
    
    // Close the file
    file.close();

    std::cout << "Matrix written to " << filename << std::endl;
}

void assignSparseBlockWithTriplet(std::vector<Tri> &tripletList, const MatrixType &Block, int ibegin, int jbegin)
{
    // First, assemble Triplet List of the Block
    std::vector<Tri> tripletListBlock;
    tripletListBlock.reserve(Block.nonZeros());
    for (int k = 0; k < Block.outerSize(); ++k) {
        for (typename MatrixType::InnerIterator it(Block, k); it; ++it) {
            // Each element can be accessed by (it.row(), it.col(), it.value())
            tripletListBlock.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        }
    }

}

void assignSparseBlockInplace(MatrixType &Base, const MatrixType &Block, int ibegin, int jbegin,
                            std::vector<Tri> &BaseTriplet)
{

    assert(ibegin + Block.rows() <= Base.rows() && jbegin + Block.cols() <= Base.cols());

    // First, collect Eigen Triplets of Base:
    if(BaseTriplet.size() == 0){
        BaseTriplet.reserve(Base.nonZeros() + Block.nonZeros());

        for (int k = 0; k < Base.outerSize(); ++k) {
            for (typename MatrixType::InnerIterator it(Base, k); it; ++it) {
                // If we are outside the range of the block matrix, we keep the entries
                if(it.row() < ibegin || it.row() >= ibegin + Block.rows() || it.col() < jbegin || it.col() >= jbegin + Block.cols())
                {
                    BaseTriplet.push_back(Tri(it.row(), it.col(), it.value()));
                }
            }
        }

    }
    else if(BaseTriplet.size() == Base.nonZeros()){
        // BaseTriplet has already been initialized
        // First, remove all triplets that are in the range of the Block matrix
        for(int i = 0; i < BaseTriplet.size(); i++)
        {
            Tri it = BaseTriplet[i];
            if(it.row() >= ibegin && it.row() < ibegin + Block.rows() && it.col() >= jbegin && it.col() < jbegin + Block.cols())
            {
                    BaseTriplet.erase(BaseTriplet.begin()+i);
            }
        }
    }

    for (int k = 0; k < Block.outerSize(); ++k) {
        for (typename MatrixType::InnerIterator it(Block, k); it; ++it) {
            // Push back all entries of the block matrix
            BaseTriplet.push_back(Tri(it.row() + ibegin, it.col() + jbegin, it.value()));
        }
    }

    Base.setFromTriplets(BaseTriplet.begin(), BaseTriplet.end());
}

bool is_posdef_eigen(FullMatrixType H) {
    Eigen::EigenSolver<FullMatrixType> eigensolver(H);
    for (int i = 0; i < eigensolver.eigenvalues().rows(); i++) {
        double v = eigensolver.eigenvalues()(i).real();
        if (v <= 0) {
            return false;
        }
    }
    return true;
}

MatrixType vectorToSparseMat(const VectorType &vec, bool asRowVector = false) {
    // Create a sparse matrix with appropriate dimensions
    MatrixType sparseMat(asRowVector ? 1 : vec.size(), asRowVector ? vec.size() : 1);

    // Reserve space for non-zero entries
    sparseMat.reserve(vec.size()); // At most, there could be `vec.size()` non-zero elements.

    // Add non-zero elements to the sparse matrix
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i] != 0) {
            if (asRowVector) {
                sparseMat.insert(0, i) = vec[i]; // Insert as a row vector
            } else {
                sparseMat.insert(i, 0) = vec[i]; // Insert as a column vector
            }
        }
    }

    sparseMat.makeCompressed(); // Compress the sparse matrix for efficient storage
    return sparseMat;
}
#endif // SPARSEMAT_H