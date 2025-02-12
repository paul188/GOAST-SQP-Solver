#include "Levenberg_Marquardt.h"
#include "Constraints.h"
#include <goast/Core.h>
#include <random>
#include <iostream>
#include <utility>

template<typename ConfiguratorType>
class CostFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>
{
    public:
        using RealType = typename ConfiguratorType::RealType;
        using VectorType = typename ConfiguratorType::VectorType;
        using MatrixType = typename ConfiguratorType::SparseMatrixType;
    private:
        const VectorType &_points;
    public:
        CostFunctional(VectorType &points):
        _points(points){}
        
        void apply(const VectorType &p, VectorType &Dest) const override
        {
            Dest.resize(_points.size()/2);
            for(int i = 0; i < Dest.size(); i++)
            {
                Dest[i] = (p[0] + p[1]*_points[i]) - _points[i + (Dest.size())];
            }
        }
};

template<typename ConfiguratorType>
class CostFunctionalGrad : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>
{
    public:
        using RealType = typename ConfiguratorType::RealType;
        using VectorType = typename ConfiguratorType::VectorType;
        using MatrixType = typename ConfiguratorType::SparseMatrixType;
    private:
        const VectorType &_points;
    public:
        CostFunctionalGrad(VectorType &points):
        _points(points){}
        
        void apply(const VectorType &p, MatrixType &Dest) const override
        {
            if(Dest.rows() != _points.size()/2 || Dest.cols() != 2)
            {
                Dest.resize(_points.size()/2, 2);
            }

            for(int i = 0; i < Dest.rows(); i++)
            {
                Dest.coeffRef(i,0) = 1.0;
                Dest.coeffRef(i,1) = _points[i];
            }

        }
};

int main()
{

    using VectorType = typename DefaultConfigurator::VectorType;
    using MatrixType = typename DefaultConfigurator::SparseMatrixType;

    LevenbergMarquardtParams<DefaultConfigurator> pars;
    std::random_device rd; // Seed for the random number engine
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_real_distribution<> dis(0.0, 1.0); // Uniform distribution in [0, 1]

    // Vector to store the random points
    VectorType points;
    points.resize(20);

    // Generate the points
    for (int i = 0; i < (points.size()/2); ++i) {
        double x = dis(gen);
        double y = dis(gen);
        std::cout<<"Point: "<<x<<" "<<y<<std::endl;
        points[i] = x;
        points[i + (points.size()/2)] = y;
    }

    CostFunctional<DefaultConfigurator> costFunctional(points);
    CostFunctionalGrad<DefaultConfigurator> gradCostFunctional(points);
    LMAlgorithm<DefaultConfigurator> lm(pars, costFunctional, gradCostFunctional, 2, (points.size()/2));
    VectorType Dest;
    Dest.resize(2);
    VectorType init = VectorType::Ones(2);
    MatrixType W(points.size()/2, points.size()/2);
    W.setIdentity();
    lm.solve(init, Dest, W);

    std::cout<<"Line parameters: "<<Dest[0]<<" "<<Dest[1]<<std::endl;
    return 0;
}