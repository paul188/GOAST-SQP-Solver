template <typename ConfiguratorType>
struct SQP_Parameters{
    typedef typename ConfiguratorType::RealType RealType;
    constexpr static RealType tau = 0.5; // parameter tau for the stepsize control, how much smaller to make alpha after each failed iteration
    constexpr static RealType eta = 0.25; // criterion eta for when we accept a stepsize
    constexpr static RealType rho = 0.5; // parameter rho for the stepsize control update
    constexpr static RealType eps = 1e-6; // convergence criterion
    constexpr static RealType theta = 0.2; // parameter theta for the BFGS secant update
    constexpr static size_t maxiter = 10^4; // maximum number of iterations
    RealType mu = 1.0; // penalty parameter for merit functions to be determined
    RealType alpha = 1.0; //stepsize to be updated each iteration
    size_t iter = 0; // iteration counter
};

template <typename ConfiguratorType>
void BFGS_update(const typename ConfiguratorType::VectorType &y_k, const typename ConfiguratorType::VectorType &s_k, const typename ConfiguratorType::RealType theta,typename ConfiguratorType::FullMatrixType &B_k){
    typename ConfiguratorType::RealType sigma = 0;
    typename ConfiguratorType::VectorType y_k2 = y_k;
    auto Bs = B_k*s_k;
    auto sBs = s_k.dot(Bs);
    auto sy = s_k.dot(y_k);
    if(sy < theta*sBs){
        sigma = (1-theta)*sBs/(sBs - sy);
        y_k2 = sigma*y_k + (1-sigma)*Bs;
    }

    B_k += y_k2*y_k2.transpose()/(sy) - Bs*(Bs.transpose())/(sBs);
}