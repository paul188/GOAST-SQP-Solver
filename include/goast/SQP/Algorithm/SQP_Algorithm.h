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
        RealType eps = 1e-12; // Convergence criterion
        size_t maxIter = 40000; // maximum number of iterations
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

        static constexpr RealType beta = 0.9; // Parameter for momentum method
        static constexpr RealType sigma = 1.0; // Parameter for momentum method
        
        RealType mu = 1.0; // Merit control parameter
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
                            size_t BFGS_reset = 20,
                            // dont really need second order correction here as our constraint is easy to fulfill
                            // if other constraints added, this might become useful
                            bool secondOrderCorrection = false,
                            // The algorithm will make big trial steps in direction of the d_k_fold_dofs to cover a larger search space more quickly
                            // The normal SQP algorithm is probably a little more reliable, but makes extremely small steps ~1e-8
                            bool bigStep = true) 
                            : SQPBaseSolver<ConfiguratorType>(pars),
                            _costFunctional(costFunctional),
                            _DcostFunctional(DcostFunctional),
                            _factory(factory),
                            _boundaryDOFs(boundaryDOFs),
                            _problemDOFs(problemDOFs),
                            _BFGS_reset(BFGS_reset),
                            _secondOrderCorrection(secondOrderCorrection),
                            _bigStep(bigStep){}

        // plateGeomRef_basic is the starting reference geometry of the plate
        // before we translate it with the fold dofs
        // plateGeomDef contains the initial deformed geometry. But since is also contains 
        // our degrees of freedom as the vertexDOFs, we will modify it during the run
        void solve(const VectorType& plateGeomRef_basic, std::vector<VectorType> &def_geometries, std::vector<VectorType> &ref_geometries, std::vector<VectorType> &fold_DOFs, std::optional<std::pair<VectorType,VectorType>> foldDOF_bounds = std::nullopt)
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
            VectorType rhs_k = VectorType::Zero(nEffectiveDOFs + nEffectiveVertexDOFs);

            // Reserve the memory for the Jacobian of the constraint
            MatrixType C_k(nEffectiveVertexDOFs, nEffectiveDOFs);

            // Create solution vector for linear system that contains as subvector the Lagrange multipliers
            VectorType sol_k = VectorType::Zero(nEffectiveDOFs + nEffectiveVertexDOFs);

            // Create the subvectors of the solution vector
            VectorType d_k = VectorType::Zero(nEffectiveDOFs);
            VectorType lambda_kplus1 = VectorType::Zero(nEffectiveVertexDOFs);
            VectorType lambda_k = VectorType::Zero(nEffectiveVertexDOFs);
            VectorType lambda_diff = VectorType::Zero(nEffectiveVertexDOFs);

            // difference of the gradient of the lagrangian w.r.t. x between the last iteration and this one
            // Set initially to 1.0 so that convergence criterion is not immediately met
            VectorType DiffGradxLagrange = VectorType::Ones(nEffectiveDOFs);

            // DiffConstraint as derivative w.r.t. lambda of the Lagrangian
            VectorType DiffConstraint =  VectorType::Ones(nEffectiveVertexDOFs);

            size_t nBigSteps = 0;
            size_t terminationCounter = 0;
            RealType variationFoldDOFs = 0;

            // Unfortunately, the convergence criterion doesnt quite work yet -> can see improvements even long after "convergence"
            while(_pars.iter < _pars.maxIter && !convergenceCriterion(d_k, DiffGradxLagrange, DiffConstraint, Constraint, variationFoldDOFs, terminationCounter, nFoldDOFs)){

                // Solve the QP Problem and write the sol into sol_k

                //printSparseMatrix(B_k,"B_k_" + std::to_string(_pars.iter) + ".txt");
            
                std::cout << "\rIteration: "<< std::setw(4) << _pars.iter
                << " | CostFunctional_val: " << std::setprecision(12)<< std::setw(10) << costFunctional_val
                << " | DE_val: " << std::setprecision(12)<<std::setw(10) << std::setprecision(4) << norm_l1(Constraint)
                << " | CostFunctional_val: " << std::setprecision(12)<<std::setw(10) << std::setprecision(4) << costFunctional_val
                << " | Fold Dofs: "<< std::setprecision(12) << _problemDOFs.getFoldDOFs()[0]
                << " | D_x L: "<< std::setprecision(12) << norm_l1(DiffGradxLagrange)
                << " | D_lambda L: "<< std::setprecision(12) << norm_l1(DiffConstraint)
                << " | d_k: "<< std::setprecision(12) << norm_l1(d_k)
                << std::endl;  // flush to ensure immediate output

                A_k.setZero();
                C_k.setZero();

                std::vector<Eigen::Triplet<RealType>> CTriplet;
                assignSparseBlockInplace(C_k, D2E_mix,0,0, CTriplet);
                assignSparseBlockInplace(C_k, D2E_vertex, 0, nFoldDOFs, CTriplet);

                std::vector<Eigen::Triplet<RealType>> ATriplet;
                assignSparseBlockInplace(A_k, B_k, 0, 0, ATriplet);
                assignSparseBlockInplace(A_k, -C_k.transpose(),0,nEffectiveDOFs, ATriplet);
                assignSparseBlockInplace(A_k, C_k, nEffectiveDOFs,0, ATriplet);

                // Construct the r.h.s for the linear system
                // the derivative of the cost functional
                rhs_k.segment(0,nEffectiveDOFs) = -DCostFunctional_val;
                // then the derivative of the energy, i.e. the constraint
                rhs_k.segment(nEffectiveDOFs,nEffectiveVertexDOFs) = -Constraint;

                if(A_k.rows() <= 20000){
                    // Use a dense solver instead for such a small system -> sparse solvers could fail here
                    FullMatrixType A_k_dense = A_k.toDense();

                    Eigen::FullPivLU<FullMatrixType> lu_decomp(A_k_dense);
                    sol_k = lu_decomp.solve(rhs_k);

                    /*
                    Eigen::JacobiSVD<FullMatrixType> svd(A_k_dense, Eigen::ComputeThinU | Eigen::ComputeThinV);
                    sol_k = svd.solve(rhs_k);
                    */

                    if(!(A_k*sol_k).isApprox(rhs_k)){
                        //std::cerr << "Error: Linear dense system not solved correctly" << std::endl;
                        //std::cout<<"Residual: " << (A_k*sol_k - rhs_k).norm()<<std::endl;
                    }
                }
                else{
                    Eigen::BiCGSTAB<MatrixType> solver;
                    solver.compute(A_k);
                    sol_k = solver.solve(rhs_k);
                    //LinearSolver<ConfiguratorType> solver;
                    //solver.prepareSolver(A_k);
                    //solver.backSubstitute(rhs_k, sol_k);

                    if(!(A_k*sol_k).isApprox(rhs_k)){
                        std::cerr << "Error: Linear sparse system not solved correctly" << std::endl;
                        std::cout<< "Error norm: "<<norm_l1(A_k*sol_k - rhs_k)<<std::endl;
                    }
                }

                d_k = sol_k.segment(0, nEffectiveDOFs);
                lambda_kplus1 = sol_k.segment(nEffectiveDOFs, nEffectiveVertexDOFs);
                lambda_diff = lambda_kplus1 - lambda_k;

                // Now, determine the stepsize

                RealType comparisonVal = ((DCostFunctional_val.dot(d_k) + 0.5*d_k.transpose().dot(B_k*d_k))/((1-_pars.rho)*Constraint.template lpNorm<1>()));

                if(_pars.mu < comparisonVal){
                    _pars.mu = (1+1e-6)*comparisonVal; // set mu to a value that is larger than the comparison value
                }

                VectorType d_k_full = d_k;
                // Set the boundary DOFs to zero to be able to add this to the problemDOFs
                _boundaryDOFs.InverseTransformWithFoldDofs(d_k_full);

                ProblemDOFs<ConfiguratorType> lineSearchDOFs(_problemDOFs);

                // Algorithm modifies d_k_full and lambda_diff to be the correct steps
                line_search(costFunctional_val,
                            DCostFunctional_val,
                            d_k,
                            d_k_full,
                            lambda_diff,
                            lambda_k,
                            Constraint,
                            C_k,
                            B_k,
                            A_k,
                            rhs_k,
                            nEffectiveDOFs,
                            nEffectiveVertexDOFs,
                            nFoldDOFs,
                            lineSearchDOFs);

                if(norm_l1(Constraint) < 1e-4 && _bigStep && nBigSteps < 15){
                    std::cout<<"Big Step: "<<std::endl;
                    makeBigStep(d_k_full, foldDOF_bounds, nFoldDOFs, nBigSteps);
                }

                // Apply the steps
                _problemDOFs += d_k_full;
                lambda_kplus1 = lambda_k + lambda_diff;

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
                d_k = d_k_full;
                _boundaryDOFs.transformWithFoldDofsToReducedSpace(d_k);
                VectorType s_k = d_k;

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
                
                std::vector<Eigen::Triplet<RealType>> temp_tripletList;
                assignSparseBlockInplace(temp, D2E_mix_kplus1 - D2E_mix, 0,0, temp_tripletList);
                assignSparseBlockInplace(temp, D2E_vertex_kplus1 - D2E_vertex, 0,nFoldDOFs, temp_tripletList);

                DiffGradxEnergy_kplus1 = lambda_kplus1.transpose()*(temp);
                std::cout<<"Dimension test: "<<std::endl;
                std::cout<<temp.rows()<<", "<<temp.cols()<<std::endl;
                std::cout<<lambda_kplus1.size()<<std::endl;
                temp_tripletList.clear();
                temp.setZero();

                DiffGradxLagrange = DiffGradxCostFunctional + DiffGradxEnergy_kplus1;

                BFGS_update(DiffGradxLagrange, s_k, B_k);

                // Difference in constraint to check for convergence
                DiffConstraint = Constraint_kplus1 - Constraint;

                // After all the calculations, we can update the values
                Constraint = Constraint_kplus1;
                DCostFunctional_val = DCostFunctional_kplus1_val;
                costFunctional_val = costFunctional_kplus1_val;
                D2E_vertex = D2E_vertex_kplus1;
                D2E_mix = D2E_mix_kplus1;
                lambda_k = lambda_kplus1;
                
                if(_pars.iter %100 == 0){
                    def_geometries.push_back(_problemDOFs.getVertexDOFs());
                    ref_geometries.push_back(_problemDOFs.getReferenceGeometry());
                    fold_DOFs.push_back(_problemDOFs.getFoldDOFs());
                }
                _pars.iter++;
            }
        }

        void  BFGS_update(const VectorType &y_k, const VectorType &s_k, MatrixType &B_k){

            if(_pars.iter > 0 && (_pars.iter % _BFGS_reset == 0)){
                B_k.setIdentity();
                return;
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

            // Extremely small values lead to nans in B_k
            if(replaceNanWithZero(B_k)){
                MatrixType Id(B_k.rows(),B_k.cols());
                Id.setIdentity();
                B_k += Id*1e-6;
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

        void line_search(const RealType CostFunctional_val,
                            const VectorType &CostFunctionalGradient_val,
                            const VectorType &d_k, 
                            VectorType &d_k_full,
                            VectorType &lambda_diff,
                            const VectorType &lambda_k,
                            const VectorType &Constraint,
                            const MatrixType &C_k,
                            const MatrixType &B_k,
                            const MatrixType &A_k,
                            VectorType &rhs_k,
                            size_t nEffectiveDOFs,
                            size_t nEffectiveVertexDOFs,
                            size_t nFoldDOFs,
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
            
            RealType CostFunctional_val_linesearch;
            VectorType Constraint_linesearch;
            RealType constr_l1_linesearch;

            while(alpha > 1e-13) {

                // Update deformed and reference geometry
                lineSearchDOFs += alpha*d_k_full;

                _costFunctional.apply(lineSearchDOFs.getVertexDOFs(),CostFunctional_val_linesearch);
                
                _factory.produceDE_vertex(lineSearchDOFs,Constraint_linesearch);
                _boundaryDOFs.transformToReducedSpace(Constraint_linesearch);

                constr_l1_linesearch = norm_l1(Constraint_linesearch);
                RealType phi_l1_step = CostFunctional_val_linesearch + mu * constr_l1_linesearch;

                if (phi_l1_step <= (phi_l1 + alpha * _pars.eta * Dp_phi_l1)) {
                    d_k_full *= alpha;
                    lambda_diff *= alpha;
                    return;
                } 
                 // If increase in merit together with increase in constraint use second order correction
                else if ((constr_l1_linesearch > constr_l1) && (alpha == 1.0) && _secondOrderCorrection)
                {
                    // Compute second order correction
                    VectorType d = Constraint_linesearch - C_k*d_k;
                    VectorType sol_k;

                    // Solve the new quadratic problem -> we can reuse the global quadratic problem
                    // Just need to modify the lower segment of the rhs
                    rhs_k.segment(nEffectiveDOFs,nEffectiveVertexDOFs) = -d;

                    // Now, solve the quadratic problem
                    if(A_k.rows() <= 20000){
                        // Use a dense solver instead for such a small system -> sparse solvers could fail here
                        FullMatrixType A_k_dense = A_k.toDense();

                        Eigen::FullPivLU<FullMatrixType> lu_decomp(A_k_dense);
                        sol_k = lu_decomp.solve(rhs_k);

                        if(!(A_k*sol_k).isApprox(rhs_k)){
                            //std::cerr << "Error: Linear dense system not solved correctly" << std::endl;
                            //std::cout<<"Residual: " << (A_k*sol_k - rhs_k).norm()<<std::endl;
                        }
                    }
                    else{
                        Eigen::BiCGSTAB<MatrixType> solver;
                        solver.compute(A_k);
                        sol_k = solver.solve(rhs_k);
                        //LinearSolver<ConfiguratorType> solver;
                        //solver.prepareSolver(A_k);
                        //solver.backSubstitute(rhs_k, sol_k);

                        if(!(A_k*sol_k).isApprox(rhs_k))
                            std::cerr << "Error: Linear sparse system not solved correctly" << std::endl;
                    }

                    VectorType d_k_correction = sol_k.segment(0, nEffectiveDOFs);
                    _boundaryDOFs.InverseTransformWithFoldDofs(d_k_correction);

                    ProblemDOFs<ConfiguratorType> lineSearchDOFs_corrected = lineSearchDOFs;
                    lineSearchDOFs_corrected += d_k_correction;

                    VectorType DE_val_corrected;
                    _factory.produceDE_vertex(lineSearchDOFs_corrected, DE_val_corrected);

                    RealType CostFunctional_val_corrected;
                    _costFunctional.apply(lineSearchDOFs_corrected.getVertexDOFs(), CostFunctional_val_corrected);
                    phi_l1_step = CostFunctional_val_corrected + mu * norm_l1(DE_val_corrected);

                    if(phi_l1_step <= (phi_l1 + _pars.eta * Dp_phi_l1)) {
                        d_k_full += d_k_correction;
                        lambda_diff = sol_k.segment(nEffectiveDOFs, nEffectiveVertexDOFs) - lambda_k;
                        return;
                    }
                    else {
                        //undo the translation
                        lineSearchDOFs -= alpha*d_k_full;
                        alpha *= _pars.tau;
                    }
                }
                else {
                    //undo the translation
                    lineSearchDOFs -= alpha*d_k_full;
                    alpha *= _pars.tau;
                }
            }
            d_k_full *= alpha;
            lambda_diff *= alpha;
        }

        void makeBigStep(VectorType &d_k_full, std::optional<std::pair<VectorType,VectorType>> &foldDOF_bounds, size_t nFoldDOFs, size_t &nBigSteps){

            if(!foldDOF_bounds.has_value()){
                foldDOF_bounds = std::make_pair(VectorType::Zero(nFoldDOFs), VectorType::Ones(nFoldDOFs));
            }

            assert(foldDOF_bounds.value().first.size() == foldDOF_bounds.value().second.size() == nFoldDOFs);

            VectorType distToLower = _problemDOFs.getFoldDOFs() - foldDOF_bounds.value().first;
            VectorType distToUpper = foldDOF_bounds.value().second - _problemDOFs.getFoldDOFs();

            // Check if any value is outside the bounds
            if ((distToLower.array() < 0).any() || (distToUpper.array() < 0).any()) {
                throw std::runtime_error("Error: Vector x is outside the specified bounds.");
            }

            RealType maxAlpha = std::numeric_limits<RealType>::infinity();

            for (int i = 0; i < nFoldDOFs; ++i) {
                if (d_k_full[i] > 1e-8) {
                    maxAlpha = std::min(maxAlpha, (distToUpper[i] / d_k_full[i]));
                } else if (d_k_full[i] < -(1e-8)) {
                    maxAlpha = std::min(maxAlpha, (distToLower[i] / -d_k_full[i]));
                }
                // If step[i] == 0, skip since it doesn't affect bounds
            }

            if(!std::isinf(maxAlpha)){
                // Want to travel halfway to the bounds -> but scaled down exponentially with number of big steps
                d_k_full *= 0.05*maxAlpha*std::pow(0.1, nBigSteps);
                std::cout<<"Big Step: "<<d_k_full[0]<<std::endl;
                nBigSteps++;
            }else{
                // step not really suitable for scaling, too small
                return;
            }
        }

        bool convergenceCriterion(const VectorType &d_k, const VectorType &DiffGradxLagrange, const VectorType &DiffConstraint, const VectorType &Constraint, RealType &variationFoldDOFs, size_t &terminationCounter, size_t nFoldDOFs) const {
            if(!((norm_l1(d_k) < _pars.eps) && (norm_l1(DiffGradxLagrange) < _pars.eps) && (norm_l1(DiffConstraint) < _pars.eps) && norm_l1(Constraint) < _pars.eps)){
                // termination criterion is met for another iteration)
                terminationCounter = _pars.iter;
                variationFoldDOFs = 0;
            }

            if((_pars.iter - terminationCounter) >= 10  && variationFoldDOFs < 1e-7){
                std::cout<<"Termination criterion met"<<std::endl;
                return true;
            }

            variationFoldDOFs += norm_l1(d_k.segment(0, nFoldDOFs))/nFoldDOFs;
            return false;
        }

        private:
            SQPLineSearchParams<ConfiguratorType> _pars;
            CostFunctional<ConfiguratorType> _costFunctional;
            CostFunctionalGradient<ConfiguratorType> _DcostFunctional;
            MyObjectFactory<ConfiguratorType> _factory;
            BoundaryDOFS<ConfiguratorType> _boundaryDOFs;
            ProblemDOFs<ConfiguratorType> _problemDOFs;
            size_t _BFGS_reset;
            bool _secondOrderCorrection;
            bool _bigStep;

};

