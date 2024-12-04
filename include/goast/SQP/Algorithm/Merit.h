#ifndef MERIT_H
#define MERIT_H

#include <goast/Core.h>

template <typename ConfiguratorType>
class L1Merit{
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;

        const RealType _mu;

    public:
        L1Merit(const RealType mu):_mu( mu ){}

        void apply(const RealType costFunctional_val, const VectorType &DE_val, RealType &dest) const{
                dest = costFunctional_val + _mu*DE_val.template lpNorm<1>();
        }
};

// This is the Gradient of the l1 merit under the condition that x that is supplied 
// in the constructor has been generated from the SQP algorithm
template <typename ConfiguratorType>
class L1MeritGrad{
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;

        const RealType _mu;
        const MatrixType &_A_k;

    public:
        L1MeritGrad(RealType mu, const MatrixType &A_k) : _mu(mu), _A_k(A_k){}

        void apply(const VectorType &d_k, const VectorType& DE_val, const VectorType &DCostFunctional_val, RealType &dest) const{
            
            std::cout<<"Norm: "<<std::endl;
            VectorType temp = _A_k*d_k + DE_val;
            std::cout<< _mu*(temp.template lpNorm<1>()) << std::endl;
            std::cout<<"Norm end"<<std::endl;
            dest = DCostFunctional_val.dot(d_k) - _mu*DE_val.template lpNorm<1>();
        }
};

#endif