#include <goast/SQP/DOFHandling/FoldDofs.h>

// This contains the vertex as well as the fold dofs and methods to translate them
template <typename ConfiguratorType>
class ProblemDOFs{
    protected:
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;

        VectorType _foldDOFs;
        VectorType _vertexDOFs;
        VectorType _baseDeformedGeometry; // base deformed geometry, used to apply the vertex dofs often the result of Newton minimization
        std::shared_ptr<FoldDofs<ConfiguratorType>> _foldDofsPtr; // pointer to the fold dofs object
        std::shared_ptr<FoldDofsGradient<ConfiguratorType>> _DfoldDofsPtr; // pointer to the fold dofs gradient object


    public:
        ProblemDOFs(const VectorType &foldDOFs_initial, 
                    const VectorType &vertexDOFs_initial, 
                    const VectorType &baseDeformedGeometry,
                    std::shared_ptr<FoldDofs<ConfiguratorType>> foldDofsPtr,
                    std::shared_ptr<FoldDofsGradient<ConfiguratorType>> DfoldDofsPtr)
                     : _foldDOFs(foldDOFs_initial), 
                       _vertexDOFs(vertexDOFs_initial), 
                       _baseDeformedGeometry(baseDeformedGeometry),
                       _foldDofsPtr(std::move(foldDofsPtr)),
                       _DfoldDofsPtr(std::move(DfoldDofsPtr)) {}

        ProblemDOFs(const ProblemDOFs& other): 
                    _foldDOFs(other._foldDOFs),
                    _vertexDOFs(other._vertexDOFs),
                    _baseDeformedGeometry(other._baseDeformedGeometry),
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
        void apply(VectorType &geomRefDest, VectorType &geomDefDest) const{
            geomDefDest = _baseDeformedGeometry + _vertexDOFs;
            _foldDofsPtr->apply(_foldDOFs, geomRefDest);
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

        VectorType getDeformedGeometry() const{
            return _baseDeformedGeometry + _vertexDOFs;
        }

        void setFoldDOFs(const VectorType &foldDOFs){
            _foldDOFs = foldDOFs;
        }

        void setVertexDOFs(const VectorType &vertexDOFs){
            _vertexDOFs = vertexDOFs;
        }
};
