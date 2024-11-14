#ifndef MERIT_H
#define MERIT_H

#include <goast/Core.h>

template <typename ConfiguratorType>
class L1Merit {
    protected:
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::VectorType VectorType;

        const RealType _mu;

    public:
        L1Merit(const RealType mu) : _mu(mu) {}

        void apply(const VectorType &constraint, const RealType &costFunctional, RealType &dest) const{
            dest = costFunctional + constraint.template lpNorm<1>();
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

    public:
        L1MeritGrad(RealType mu) : _mu(mu) {}

        void apply(const VectorType &DCostFunctional, const VectorType& p_k, const VectorType &Constraint, RealType &dest) const{
            dest = DCostFunctional.dot(p_k) - _mu*Constraint.template lpNorm<1>();
        }
};

#endif