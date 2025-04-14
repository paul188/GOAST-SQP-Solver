#pragma once
#include <goast/Core.h>

/*
Levenberg-Marquardt Algorithm implemented according to Methods for Non-Linear Least Squares Problems (2nd ed.) by Madsen, K., Nielsen, H. B., & Tingleff, O.
*/

template <typename ConfiguratorType>
class LevenbergMarquardtParams{
    public:
        using RealType = typename ConfiguratorType::RealType;
        RealType eps_1 = 1e-12; // Convergence criterion for gradient
        RealType eps_2 = 1e-8; // Convergence criterion for step
        RealType tau = 1e-6; // Parameter tau, guess how accurate initial value is
        RealType mu = 1.0; // Regularization parameter
        RealType rho = 0.0; // Parameter determining quality of quadratic approximation
        RealType nu = 2.0; // Parameter for increasing mu
        size_t maxIter = 200; // maximum number of iterations
        size_t iter = 0; // current iteration
        bool found = false; // flag to indicate if solution was found
};

template <typename ConfiguratorType>
class LMAlgorithm{
    public:
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    private:
        LevenbergMarquardtParams<ConfiguratorType> _pars;
        BaseOp<VectorType, VectorType> &_f;
        BaseOp<VectorType, MatrixType> &_Df;
        size_t _num_constraints;
        size_t _num_dofs;
    
    public:
        LMAlgorithm(const LevenbergMarquardtParams<ConfiguratorType> &pars,
                    BaseOp<VectorType, VectorType> &f,
                    BaseOp<VectorType, MatrixType> &Df,
                    size_t num_dofs,
                    size_t num_constraints) :
                    _pars(pars),
                    _f(f),
                    _Df(Df),
                    _num_constraints(num_constraints),
                    _num_dofs(num_dofs){}

        void solve(const VectorType &plateGeomInit, VectorType &plateGeomDest, MatrixType &W){
            VectorType x_k = plateGeomInit;
            VectorType f_k;
            _f.apply(x_k, f_k);
            RealType F_k = 0.5*f_k.dot(W*f_k);
            MatrixType Df_k;
            _Df.apply(x_k, Df_k);
            VectorType g_k = Df_k.transpose()*W*f_k;
            MatrixType Hf_k = Df_k.transpose()*W*Df_k;

            VectorType x_kplus1 = VectorType::Zero(_num_dofs);
            VectorType f_kplus1 = VectorType::Zero(_num_constraints);
            RealType F_kplus1 = 0.0;
            VectorType step = VectorType::Zero(_num_dofs);

            // compute maximum absolute value on diagonal
            RealType max_abs_diag = 0.0;
            for (int k = 0; k < Hf_k.outerSize(); ++k) {
                for (typename MatrixType::InnerIterator it(Hf_k, k); it; ++it) {
                    if (it.row() == it.col()) { // Check if diagonal element
                        max_abs_diag = std::max(max_abs_diag, std::abs(it.value()));
                    }
                }
            }

            _pars.found = (g_k.template lpNorm<Eigen::Infinity>() <= _pars.eps_1);

            _pars.mu = _pars.tau*max_abs_diag;
            MatrixType Id(_num_dofs, _num_dofs);
            Id.setIdentity();

            while(_pars.iter < _pars.maxIter && !_pars.found){
                std::cout << "\rIteration: "<< std::setw(4) << _pars.iter
                << " | F_k: " << std::setprecision(12)<< std::setw(10) << 0.5*f_k.dot(W*f_k)
                << " | g_k: " << std::setprecision(12)<<std::setw(10)<<g_k.template lpNorm<Eigen::Infinity>()<<std::endl;

                F_k = 0.5*f_k.dot(W*f_k);
                g_k = Df_k.transpose()*W*f_k;
                Hf_k = Df_k.transpose()*W*Df_k;

                LinearSolver<ConfiguratorType> directSolver;
                directSolver.prepareSolver(Hf_k + _pars.mu*Id);
                directSolver.backSubstitute(-g_k, step);

                if(step.norm() <= _pars.eps_2*(x_k.norm() + _pars.eps_2)){
                    _pars.found = true;
                    plateGeomDest = x_k;
                    return;
                }

                x_kplus1 = x_k + step;
                _f.apply(x_kplus1, f_kplus1);
                F_kplus1 = 0.5*f_kplus1.dot(W*f_kplus1);

                RealType DeltaL = 0.5*step.transpose().dot(_pars.mu*step - g_k);
                _pars.rho = (F_k - F_kplus1)/DeltaL;

                if(_pars.rho > 0)
                {
                    x_k = x_kplus1;
                    _pars.mu = _pars.mu*std::max(1.0/3.0, 1 - std::pow(2*_pars.rho - 1, 3));
                    _pars.nu = 2.0;
                    _f.apply(x_k, f_k);
                    _Df.apply(x_k, Df_k);
                    g_k = Df_k.transpose()*W*f_k;

                    if(g_k.template lpNorm<Eigen::Infinity>() <= _pars.eps_1)
                    {
                        _pars.found = true;
                        plateGeomDest = x_k;
                        return;
                    }

                    F_k = 0.5*f_k.dot(W*f_k);
                    Hf_k = Df_k.transpose()*W*Df_k;
                    _pars.iter++;
                }
                else{
                    _pars.mu = _pars.mu*_pars.nu;
                    _pars.nu = 2*_pars.nu;
                }
                _pars.iter++;
            }
            plateGeomDest = x_k;
        }
    
};