#ifndef QP_H
#define QP_H

template <typename ConfiguratorType>
class BaseQP {
    public:
        BaseQP();
        ~BaseQP();
};

template <typename ConfiguratorType>
class LineSearchQP : public BaseQP <ConfiguratorType> {
    public:
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;

        LineSearchQP(const MatrixType &D2E, const VectorType &DF, const VectorType &Constraints, const MatrixType &DConstraints)
        {
            _D2E = D2E;
            _DF = DF;
            _Constraints = Constraints;
            _DConstraints = DConstraints;
        }

        void solveQP(VectorType &sol)
        {
            
        }
        ~LineSearchQP();
    protected:
        MatrixType _D2E;
        VectorType _DF;
        VectorType _Constraints;
        MatrixType _DConstraints;
};

#endif