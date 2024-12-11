#pragma once

#include <goast/Core.h>
#include <goast/SQP/Algorithm/CostFunctional.h>
#include <goast/SQP/Algorithm/Merit.h>
#include <goast/SQP/DOFHandling/BoundaryDOFS.h>
#include <goast/SQP/DOFHandling/FoldDofs.h>
#include <goast/SQP/Utils/ObjectFactory.h>

// Basic Parameters that all SQP Solvers should contain
template<typename ConfiguratorType>
class SQPBaseParams{
    public:
        using RealType = typename ConfiguratorType::RealType;
        RealType eps = 1e-6; // Convergence criterion
        size_t maxIter = 2000; // maximum number of iterations
        size_t iter = 0; // current iteration
};

// Specific parameters needed for SQP Line Search
template <typename ConfiguratorType>
class SQPLineSearchParams : public SQPBaseParams<ConfiguratorType>{
    public:
        using RealType = typename ConfiguratorType::RealType;
        static constexpr RealType tau = 0.5;
        static constexpr RealType eta = 0.5;
        static constexpr RealType rho = 0.5; // Merit control parameter
        static constexpr RealType theta = 0.2;
        
        RealType mu = 1.0; // Merit control parameter
        RealType alpha = 1.0; // Current step size
};

template <typename ConfiguratorType>
class SQPBaseSolver{
    public:
        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;

        SQPBaseSolver(const SQPBaseParams<ConfiguratorType> &pars) : _pars(pars) {}
    protected:
        SQPBaseParams<ConfiguratorType> _pars;
};

template <typename ConfiguratorType>
class SQPLineSearchSolver : SQPBaseSolver<ConfiguratorType>{

    public: 

        typedef typename ConfiguratorType::VectorType VectorType;
        typedef typename ConfiguratorType::RealType RealType;
        typedef typename ConfiguratorType::SparseMatrixType MatrixType;
        typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
        typedef typename ConfiguratorType::VecType VecType;

        // pars: Parameters for the SQP Solver
        // costFunctional: Cost functional
        // DcostFunctional: Derivative of the cost functional
        // factory: Object factory to create energies, their derivatives and hessians on demand
        // FoldDOFs: object to handle degrees of freedom of the fold
        SQPLineSearchSolver(const SQPLineSearchParams<ConfiguratorType> &pars,
                            const CostFunctional<ConfiguratorType> &costFunctional,
                            const CostFunctionalGradient<ConfiguratorType> &DcostFunctional,
                            const MyObjectFactory<ConfiguratorType> &factory,
                            const BoundaryDOFS<ConfiguratorType> &boundaryDOFs,
                            ProblemDOFs<ConfiguratorType> &problemDOFs,
                            size_t BFGS_reset = 15) 
                            : SQPBaseSolver<ConfiguratorType>(pars),
                            _costFunctional(costFunctional),
                            _DcostFunctional(DcostFunctional),
                            _factory(factory),
                            _boundaryDOFs(boundaryDOFs),
                            _problemDOFs(problemDOFs),
                            _BFGS_reset(BFGS_reset){}

        // plateGeomRef_basic is the starting reference geometry of the plate
        // before we translate it with the fold dofs
        // plateGeomDef contains the initial deformed geometry. But since is also contains 
        // our degrees of freedom as the vertexDOFs, we will modify it during the run
        void solve(const VectorType& plateGeomRef_basic, std::vector<VectorType> &def_geometries, std::vector<VectorType> &ref_geometries, std::vector<VectorType> &fold_DOFs)
        {

            size_t nAllVertexDOFs = _boundaryDOFs.getNumVertexDOFs();
            size_t nFoldDOFs = _boundaryDOFs.getNumFoldDOFs();
            size_t nAllDOFs = nAllVertexDOFs + nFoldDOFs;
            size_t bdryMaskOptSize = _boundaryDOFs.getNumBdryDOFs();
            size_t nEffectiveDOFs = nAllDOFs - bdryMaskOptSize;
            size_t nEffectiveVertexDOFs = nAllVertexDOFs - bdryMaskOptSize;

            // Initialize the approximate Hessian
            MatrixType B_k(nEffectiveDOFs, nEffectiveDOFs);
            B_k.setIdentity();

            // Create the constraint, i.e. that derivative of the energy w.r.t. the vertex DOFs is zero
            VectorType Constraint;
            _factory.produceDE_vertex(_problemDOFs, Constraint);

            // Create Hessian only w.r.t. the vertex DOFs
            MatrixType D2E_vertex;
            _factory.produceD2E_vertex(_problemDOFs, D2E_vertex);

            // Create Hessian w.r.t. the vertex and fold DOFs
            MatrixType D2E_mix;
            _factory.produceD2E_mix(_problemDOFs, D2E_mix);

            VectorType DCostFunctional_val;
            _DcostFunctional.apply(_problemDOFs.getVertexDOFs(), DCostFunctional_val);
            _boundaryDOFs.addZeroFoldDOFs(DCostFunctional_val);

            RealType costFunctional_val;
            _costFunctional.apply(_problemDOFs.getVertexDOFs(), costFunctional_val);

            // Now, reduce everything to the effective space
            _boundaryDOFs.transformToReducedSpace(Constraint);
            _boundaryDOFs.transformWithFoldDofsToReducedSpace(DCostFunctional_val);
            _boundaryDOFs.transformRowColToReducedSpace(D2E_vertex);
            _boundaryDOFs.transformRowToReducedSpace(D2E_mix);

            // Reserve memory for the large matrix and r.h.s for the linear system
            // We will do a transformation before initializing A_k s.th. we only neednReducedDOFs x nReducedDOFs matrix
            MatrixType A_k(nEffectiveDOFs + nEffectiveVertexDOFs, nEffectiveDOFs + nEffectiveVertexDOFs);
            A_k.setZero();
            VectorType rhs_k = VectorType::Zero(nEffectiveDOFs + nEffectiveVertexDOFs);

            // Create solution vector for linear system that contains as subvector the Lagrange multipliers
            VectorType sol_k = VectorType::Zero(nEffectiveDOFs + nEffectiveVertexDOFs);

            // Create the subvectors of the solution vector
            VectorType d_k = VectorType::Zero(nEffectiveDOFs);
            VectorType lambda_k = VectorType::Zero(nEffectiveVertexDOFs);;
            VectorType lambda_k_prev = VectorType::Zero(nEffectiveVertexDOFs);
            VectorType lambda_diff = VectorType::Zero(nEffectiveVertexDOFs);

            while(_pars.iter < _pars.maxIter){ //&&(norm_l1(Constraint) > _pars.eps || norm_l1(DCostFunctional_val) > _pars.eps)){

                // Solve the QP Problem and write the sol into sol_k
                // solveQPProblem(sol_k, A_k, rhs_k, B_k, D2E_val, DE_dt_sparse, DCostFunctional_val, DE_val, nEffectiveDOFs, nFoldDOFs, nEffectiveVertexDOFs);

                std::cout << "\rIteration: " << std::setw(4) << _pars.iter
                << " | DE_val: " << std::setprecision(12)<<std::setw(10) << std::setprecision(4) << norm_l1(Constraint)
                << " | CostFunctional_val: " << std::setprecision(12)<<std::setw(10) << std::setprecision(4) << costFunctional_val
                << " | Fold Dofs: "<< std::setprecision(12) << _problemDOFs.getFoldDOFs()[0]
                << std::endl;  // flush to ensure immediate output

                std::vector<Eigen::Triplet<double>> ATriplet;
                assignSparseBlockInplace(A_k, B_k, 0, 0, ATriplet);
                assignSparseBlockInplace(A_k, -D2E_mix.transpose(),0,nEffectiveDOFs, ATriplet);
                assignSparseBlockInplace(A_k, -D2E_vertex, nFoldDOFs,nEffectiveDOFs, ATriplet);
                assignSparseBlockInplace(A_k, D2E_mix, nEffectiveDOFs,0, ATriplet);
                assignSparseBlockInplace(A_k, D2E_vertex, nEffectiveDOFs,nFoldDOFs, ATriplet);

                // Construct the r.h.s for the linear system
                // the derivative of the cost functional
                rhs_k.segment(0,nEffectiveDOFs) = -DCostFunctional_val;
                // then the derivative of the energy, i.e. the constraint
                rhs_k.segment(nEffectiveDOFs,nEffectiveVertexDOFs) = -Constraint;

                if(A_k.rows() <= 200){
                    // Use a dense solver instead for such a small system -> sparse solvers could fail here
                    FullMatrixType A_k_dense = A_k.toDense();
                    Eigen::JacobiSVD<FullMatrixType> svd(A_k_dense, Eigen::ComputeThinU | Eigen::ComputeThinV);
                    sol_k = svd.solve(rhs_k);

                    VectorType singularValues = svd.singularValues();

                    if(!(A_k*sol_k).isApprox(rhs_k)){
                        //std::cerr << "Error: Linear dense system not solved correctly" << std::endl;
                        //std::cout<<"Residual: " << (A_k*sol_k - rhs_k).norm()<<std::endl;
                    }
                }
                else{
                    LinearSolver<ConfiguratorType> solver;
                    solver.prepareSolver(A_k);
                    solver.backSubstitute(rhs_k, sol_k);

                    if(!(A_k*sol_k).isApprox(rhs_k))
                        std::cerr << "Error: Linear sparse system not solved correctly" << std::endl;
                }

                A_k.setZero();

                d_k = sol_k.segment(0, nEffectiveDOFs);
                lambda_k = sol_k.segment(nEffectiveDOFs, nEffectiveVertexDOFs);
                lambda_diff = lambda_k - lambda_k_prev;

                sol_k.segment(nEffectiveDOFs, nEffectiveVertexDOFs) = lambda_diff;

                // Now, determine the stepsize

                RealType comparisonVal = ((DCostFunctional_val.dot(d_k) + 0.5*d_k.transpose().dot(B_k*d_k))/((1-_pars.rho)*Constraint.template lpNorm<1>()));

                if(_pars.mu < comparisonVal){
                    _pars.mu = (1+1e-6)*comparisonVal; // set mu to a value that is larger than the comparison value
                }

                VectorType d_k_full = d_k;
                // Set the boundary DOFs to zero to be able to add this to the problemDOFs
                _boundaryDOFs.InverseTransformWithFoldDofs(d_k_full);

                ProblemDOFs<ConfiguratorType> lineSearchDOFs(_problemDOFs);

                _pars.alpha = line_search(costFunctional_val,
                                          DCostFunctional_val,
                                          d_k,
                                          d_k_full,
                                          Constraint,
                                          B_k,
                                          lineSearchDOFs);

                // Apply the steps
                _problemDOFs += _pars.alpha*d_k_full;
                lambda_k = lambda_k_prev + _pars.alpha*lambda_diff;

                //Evaluate J[x_k+1]
                RealType costFunctional_kplus1_val;
                // can use the same costFunctional object, since foldVertices dont change
                _costFunctional.apply(_problemDOFs.getVertexDOFs(), costFunctional_kplus1_val);

                // Evaluate grad J[x_k+1]
                VectorType DCostFunctional_kplus1_val;
                _DcostFunctional.apply(_problemDOFs.getVertexDOFs(), DCostFunctional_kplus1_val);
                _boundaryDOFs.addZeroFoldDOFs(DCostFunctional_kplus1_val);
                _boundaryDOFs.transformWithFoldDofsToReducedSpace(DCostFunctional_kplus1_val);

                VectorType Constraint_kplus1;
                _factory.produceDE_vertex(_problemDOFs, Constraint_kplus1);
                _boundaryDOFs.transformToReducedSpace(Constraint_kplus1);

                // Update B_k using BFGS update
                VectorType s_k = _pars.alpha*d_k;
                // Want to calculate \grad_{x} L(x_{k+1}) - \grad_{x} L(x_{k})
                VectorType DiffGradxCostFunctional;
                DiffGradxCostFunctional = DCostFunctional_kplus1_val - DCostFunctional_val;

                VectorType DiffGradxEnergy_kplus1(nEffectiveVertexDOFs);
                MatrixType temp(nEffectiveVertexDOFs,nEffectiveDOFs);
                // Recalculate D2E_val

                MatrixType D2E_vertex_kplus1;
                _factory.produceD2E_vertex(_problemDOFs, D2E_vertex_kplus1);
                _boundaryDOFs.transformRowColToReducedSpace(D2E_vertex_kplus1);

                // Recalculate DE/dt:
                MatrixType D2E_mix_kplus1;
                _factory.produceD2E_mix(_problemDOFs, D2E_mix_kplus1);
               
                _boundaryDOFs.transformRowToReducedSpace(D2E_mix_kplus1);
                
                std::vector<Eigen::Triplet<double>> temp_tripletList;
                assignSparseBlockInplace(temp, D2E_mix_kplus1, 0,0, temp_tripletList);
                assignSparseBlockInplace(temp, D2E_vertex_kplus1, 0,nFoldDOFs, temp_tripletList);

                DiffGradxEnergy_kplus1 = lambda_k.transpose()*(temp);
                temp_tripletList.clear();
                temp.setZero();

                VectorType y_k = DiffGradxCostFunctional + DiffGradxEnergy_kplus1;

                BFGS_update(y_k, s_k, B_k);

                // After all the calculations, we can update the values
                Constraint = Constraint_kplus1;
                DCostFunctional_val = DCostFunctional_kplus1_val;
                costFunctional_val = costFunctional_kplus1_val;
                D2E_vertex = D2E_vertex_kplus1;
                D2E_mix = D2E_mix_kplus1;
                lambda_k_prev = lambda_k;
                
                if(_pars.iter %100 == 0){
                    def_geometries.push_back(_problemDOFs.getVertexDOFs());
                    ref_geometries.push_back(_problemDOFs.getReferenceGeometry());
                    fold_DOFs.push_back(_problemDOFs.getFoldDOFs());
                }
                _pars.iter++;
            }
        }

        void BFGS_update(const VectorType &y_k, const VectorType &s_k, MatrixType &B_k){

            if(_pars.iter > 0 && (_pars.iter % _BFGS_reset == 0)){
                B_k.setIdentity();
            }

            RealType sigma = 0;
            VectorType y_k2 = y_k;
            auto Bs = B_k*s_k;
            auto sBs = s_k.dot(Bs);
            auto sy = s_k.dot(y_k);
            if(sy < _pars.theta*sBs){
                sigma = (1-_pars.theta)*sBs/(sBs - sy);
                y_k2 = sigma*y_k + (1-sigma)*Bs;
            }

            RealType yk2_sk = s_k.dot(y_k2);

            MatrixType yk2_sparse = convertVecToSparseMat(y_k2);
            MatrixType Bs_sparse = convertVecToSparseMat(Bs);

            MatrixType outer_yk2 = yk2_sparse * yk2_sparse.transpose() / yk2_sk;
            MatrixType outer_Bs = Bs_sparse * Bs_sparse.transpose() / sBs;

            B_k += outer_yk2 - outer_Bs;

            if (is_nan_sparse(B_k)) {
                std::cout << "Hessian is NaN\n";
            }
        }     

        RealType norm_l1(const VectorType &constr) const {
            // avoid division by zero
            RealType c_l1 = std::numeric_limits<RealType>::epsilon();

            // l <= c(x) <= u
            c_l1 += (- constr).cwiseMax(0.0).sum();
            c_l1 += constr.cwiseMax(0.0).sum();

            return c_l1;
        }

        RealType line_search(const RealType CostFunctional_val,
                            const VectorType &CostFunctionalGradient_val,
                            const VectorType &d_k, 
                            const VectorType &d_k_full,
                            const VectorType &Constraint,
                            const MatrixType &B_k,
                            ProblemDOFs<ConfiguratorType> &lineSearchDOFs) {

            RealType mu, phi_l1, Dp_phi_l1;
            const RealType tau = _pars.tau;  // line search step decrease, 0 < tau < settings.tau

            RealType constr_l1 = norm_l1(Constraint);

            RealType dBd = d_k.dot(B_k * d_k);
            RealType Cd = CostFunctionalGradient_val.dot(d_k);

            if(dBd > 0){
                mu = (Cd + 0.5 * dBd) / ((1 - _pars.rho) * constr_l1);
            }else{
                mu = Cd / ((1 - _pars.rho) * constr_l1);
            }

            if(mu < 0){
                mu*=-1;
            }

            phi_l1 = CostFunctional_val + (mu * constr_l1);
            Dp_phi_l1 = Cd - (mu * constr_l1);

            RealType alpha = 1.0;
            int i;
            
            RealType CostFunctional_val_linesearch;
            VectorType Constraint_linesearch;

            while(alpha > 1e-12) {

                // Update deformed and reference geometry
                lineSearchDOFs += alpha*d_k_full;

                _costFunctional.apply(lineSearchDOFs.getVertexDOFs(),CostFunctional_val_linesearch);
                
                _factory.produceDE_vertex(lineSearchDOFs,Constraint_linesearch);
                _boundaryDOFs.transformToReducedSpace(Constraint_linesearch);

                RealType phi_l1_step = CostFunctional_val_linesearch + mu * norm_l1(Constraint_linesearch);

                //undo the translation
                lineSearchDOFs -= alpha*d_k_full;

                if (phi_l1_step <= (phi_l1 + alpha * _pars.eta * Dp_phi_l1)) {
                    break;
                } else {
                    alpha *= _pars.tau;
                }
            }
            return alpha;
        }

        private:
            SQPLineSearchParams<ConfiguratorType> _pars;
            CostFunctional<ConfiguratorType> _costFunctional;
            CostFunctionalGradient<ConfiguratorType> _DcostFunctional;
            MyObjectFactory<ConfiguratorType> _factory;
            BoundaryDOFS<ConfiguratorType> _boundaryDOFs;
            ProblemDOFs<ConfiguratorType> _problemDOFs;
            size_t _BFGS_reset;

};