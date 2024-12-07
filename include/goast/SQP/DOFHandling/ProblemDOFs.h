#include <goast/SQP/DOFHandling/FoldDofs.h>

// This contains the vertex as well as the fold dofs and methods to translate them
template <typename ConfiguratorType>
class ProblemDOFs{
    protected:
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;

    public:
        ProblemDOFs(const VectorType &foldDOFs, 
                    const VectorType &vertexDOFs, 
                    std::shared_ptr<FoldDofs<ConfiguratorType>> foldDofsPtr,
                    std::shared_ptr<FoldDofsGradient<ConfiguratorType>> DfoldDofsPtr)
                     : _foldDOFs(foldDOFs), 
                       _vertexDOFs(vertexDOFs), 
                       _foldDofsPtr(std::move(foldDofsPtr)),
                       _DfoldDofsPtr(std::move(DfoldDofsPtr)) {}

        ProblemDOFs(const ProblemDOFs& other): 
                    _foldDOFs(other._foldDOFs),
                    _vertexDOFs(other._vertexDOFs),
                    _foldDofsPtr(other._foldDofsPtr),
                    _DfoldDofsPtr(other._DfoldDofsPtr){}
        
        ProblemDOFs& operator+=(const VectorType &other) {
            _foldDOFs += other.segment(0, _foldDOFs.size());
            _vertexDOFs += other.segment(_foldDOFs.size(), _vertexDOFs.size());
            return *this;
        }

        ProblemDOFs& operator-=(const VectorType &other) {
            _foldDOFs -= other.segment(0, _foldDOFs.size());
            _vertexDOFs -= other.segment(_foldDOFs.size(), _vertexDOFs.size());
            return *this;
        }

        ProblemDOFs& operator*=(const RealType &factor) {
            _foldDOFs *= factor;
            _vertexDOFs *= factor;
            return *this;
        }

        // apply the degrees of freedom to the respective geometries
        void apply(VectorType &plateGeomRef, VectorType &plateGeomDef) const{
            plateGeomDef += _vertexDOFs;
            _foldDofsPtr->apply(_foldDOFs, plateGeomRef);
        }

        VectorType getReferenceGeometry() const{
            VectorType temp;
            _foldDofsPtr->apply(_foldDOFs, temp);
            return temp;
        }

        MatrixType getReferenceGeometryGradient() const{
            MatrixType temp;
            _DfoldDofsPtr->apply(_foldDOFs, temp);
            return temp;
        }

        VectorType getFoldDOFs() const{
            return _foldDOFs;
        }

        VectorType getVertexDOFs() const{
            return _vertexDOFs;
        }

        void setFoldDOFs(const VectorType &foldDOFs){
            _foldDOFs = foldDOFs;
        }

        void setVertexDOFs(const VectorType &vertexDOFs){
            _vertexDOFs = vertexDOFs;
        }

    private:
        VectorType _foldDOFs;
        VectorType _vertexDOFs;
        std::shared_ptr<FoldDofs<ConfiguratorType>> _foldDofsPtr;
        std::shared_ptr<FoldDofsGradient<ConfiguratorType>> _DfoldDofsPtr; 
};
