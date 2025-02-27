#pragma once
#include <goast/Core.h>
#include <goast/SQP/Algorithm/CostFunctional.h>
#include <goast/SQP/Algorithm/Merit.h>
#include <goast/SQP/DOFHandling/BoundaryDOFS.h>
#include <goast/SQP/DOFHandling/FoldDofs.h>
#include <goast/SQP/Utils/ObjectFactory.h>
#include <goast/SQP/Utils/SparseMat.h>
#include <goast/SQP/Utils/ScipySolver.h>

// Basic Parameters that all SQP Solvers should contain
template<typename ConfiguratorType>
class SQPBaseParams{
    public:
        using RealType = typename ConfiguratorType::RealType;
        RealType eps = 1e-12; // Convergence criterion
        size_t maxIter = 500; // maximum number of iterations
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
    private:
        SQPLineSearchParams<ConfiguratorType> _pars;
        CostFunctional<ConfiguratorType> _costFunctional;
        CostFunctionalGradient<ConfiguratorType> _DcostFunctional;
        std::shared_ptr<EnergyFactory<ConfiguratorType>> _factory;
        BoundaryDOFS<ConfiguratorType> _boundaryDOFs;
        ProblemDOFs<ConfiguratorType> _problemDOFs;
        size_t _BFGS_reset;
        bool _secondOrderCorrection;
        ScipySolver _scipy_solver;

    public:
        // pars: Parameters for the SQP Solver
        // costFunctional: Cost functional
        // DcostFunctional: Derivative of the cost functional
        // factory: Object factory to create energies, their derivatives and hessians on demand
        // FoldDOFs: object to handle degrees of freedom of the fold
        SQPLineSearchSolver(const SQPLineSearchParams<ConfiguratorType> &pars,
                            const CostFunctional<ConfiguratorType> &costFunctional,
                            const CostFunctionalGradient<ConfiguratorType> &DcostFunctional,
                            const std::shared_ptr<EnergyFactory<ConfiguratorType>> factory,
                            const BoundaryDOFS<ConfiguratorType> &boundaryDOFs,
                            ProblemDOFs<ConfiguratorType> &problemDOFs,
                            size_t BFGS_reset = 10,
                            bool secondOrderCorrection = true) 
                            : SQPBaseSolver<ConfiguratorType>(pars),
                            _costFunctional(costFunctional),
                            _DcostFunctional(DcostFunctional),
                            _factory(std::move(factory)),
                            _boundaryDOFs(boundaryDOFs),
                            _problemDOFs(problemDOFs),
                            _BFGS_reset(BFGS_reset),
                            _scipy_solver(boundaryDOFs.getNumVertexDOFs() - boundaryDOFs.getNumFoldDOFs()),
                            _secondOrderCorrection(secondOrderCorrection){}

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
            FullMatrixType B_k(nEffectiveDOFs, nEffectiveDOFs);
            B_k.setIdentity();

            // Approximation for the inverse of the Hessian
            FullMatrixType B_k_inv(nEffectiveDOFs, nEffectiveDOFs);
            B_k_inv.setIdentity();

            //FullMatrixType A_k(nEffectiveDOFs + nEffectiveVertexDOFs, nEffectiveDOFs + nEffectiveVertexDOFs);

            // Reserve the memory for the Jacobian of the constraint
            FullMatrixType C_k( nEffectiveVertexDOFs, nEffectiveDOFs );

            FullMatrixType C_k_inv( nEffectiveDOFs, nEffectiveVertexDOFs );

            // Create the constraint, i.e. that derivative of the energy w.r.t. the vertex DOFs is zero
            VectorType Constraint;
            _factory->produceDE_vertex(_problemDOFs, Constraint);

            // Create Hessian only w.r.t. the vertex DOFs
            MatrixType D2E_vertex;
            _factory->produceD2E_vertex(_problemDOFs, D2E_vertex);

            // Create Hessian w.r.t. the vertex and fold DOFs
            MatrixType D2E_mix;
            _factory->produceD2E_mix(_problemDOFs, D2E_mix);

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

            // Reserve memory for the r.h.s for the linear system
            VectorType rhs_k = VectorType::Zero(nEffectiveDOFs + nEffectiveVertexDOFs);

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

            FullMatrixType mult(nEffectiveVertexDOFs, nEffectiveVertexDOFs);
            FullMatrixType B_invC(nEffectiveDOFs, nEffectiveVertexDOFs);

            size_t terminationCounter = 0;

            // Unfortunately, the convergence criterion doesnt quite work yet -> can see improvements even long after "convergence"
            while(_pars.iter < _pars.maxIter){

                std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
            
                std::cout << "\rIteration: "<< std::setw(4) << _pars.iter
                << " | CostFunctional_val: " << std::setprecision(12)<< std::setw(10) << costFunctional_val
                << " | DE_val: " << std::setprecision(12)<<std::setw(10) << std::setprecision(4) << norm_l1(Constraint);
                for(int i = 0; i < _problemDOFs.getFoldDOFs().size(); i++){
                    std::cout<<" | Fold DOFs " <<i<<" : "<< std::setprecision(12) << _problemDOFs.getFoldDOFs()[i];
                }
                std::cout << " | D_x L: "<< std::setprecision(12) << norm_l1(DiffGradxLagrange)
                << " | D_lambda L: "<< std::setprecision(12) << norm_l1(DiffConstraint)
                << " | d_k: "<< std::setprecision(12) << norm_l1(d_k)
                << std::endl;  // flush to ensure immediate output

                C_k.setZero();
                C_k.block(0,0, D2E_mix.rows(), D2E_mix.cols()) = D2E_mix;
                C_k.block(0,nFoldDOFs, D2E_vertex.rows(), D2E_vertex.cols()) = D2E_vertex;          

                // Construct the r.h.s for the linear system
                // the derivative of the cost functional
                rhs_k.segment(0,nEffectiveDOFs) = -DCostFunctional_val;
                // then the derivative of the energy, i.e. the constraint
                rhs_k.segment(nEffectiveDOFs,nEffectiveVertexDOFs) = -Constraint;

                // Solve the quadratic subproblem, store the solution in d_k and lambda_kplus1
                // In the process, invert C_k and store the 
                solveQP(nFoldDOFs, nEffectiveVertexDOFs, B_k_inv,C_k,mult,B_invC, rhs_k, d_k, lambda_kplus1);

                lambda_diff = lambda_kplus1 - lambda_k;

                // Now, determine the stepsize

                RealType comparisonVal = ((DCostFunctional_val.dot(d_k) + 0.5*d_k.transpose().dot(B_k*d_k))/((1-_pars.rho)*norm_l1(Constraint)));

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
                            C_k_inv,
                            B_k,
                            D2E_vertex,
                            rhs_k,
                            nEffectiveDOFs,
                            nEffectiveVertexDOFs,
                            nFoldDOFs,
                            lineSearchDOFs);

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
                _factory->produceDE_vertex(_problemDOFs, Constraint_kplus1);
                _boundaryDOFs.transformToReducedSpace(Constraint_kplus1);

                // Want to calculate \grad_{x} L(x_{k+1}) - \grad_{x} L(x_{k})
                VectorType DiffGradxCostFunctional;
                DiffGradxCostFunctional = DCostFunctional_kplus1_val - DCostFunctional_val;

                VectorType DiffGradxEnergy_kplus1(nEffectiveVertexDOFs);
                MatrixType temp(nEffectiveVertexDOFs,nEffectiveDOFs);
                // Recalculate D2E_val

                MatrixType D2E_vertex_kplus1;
                _factory->produceD2E_vertex(_problemDOFs, D2E_vertex_kplus1);
                _boundaryDOFs.transformRowColToReducedSpace(D2E_vertex_kplus1);

                // Recalculate DE/dt:
                MatrixType D2E_mix_kplus1;
                _factory->produceD2E_mix(_problemDOFs, D2E_mix_kplus1);
               
                _boundaryDOFs.transformRowToReducedSpace(D2E_mix_kplus1);
                
                std::vector<Eigen::Triplet<RealType>> temp_tripletList;
                assignSparseBlockInplace(temp, D2E_mix_kplus1 - D2E_mix, 0,0, temp_tripletList);
                assignSparseBlockInplace(temp, D2E_vertex_kplus1 - D2E_vertex, 0,nFoldDOFs, temp_tripletList);

                DiffGradxEnergy_kplus1 = lambda_kplus1.transpose()*(temp);
                temp_tripletList.clear();
                temp.setZero();

                DiffGradxLagrange = DiffGradxCostFunctional + DiffGradxEnergy_kplus1;

                // Update BFGS approx of Hessian and inverse Hessian
                BFGS_update(DiffGradxLagrange, d_k, B_k, B_k_inv);

                // Difference in constraint to check for convergence
                DiffConstraint = Constraint_kplus1 - Constraint;

                // After all the calculations, we can update the values
                Constraint = Constraint_kplus1;
                DCostFunctional_val = DCostFunctional_kplus1_val;
                costFunctional_val = costFunctional_kplus1_val;
                D2E_vertex = D2E_vertex_kplus1;
                D2E_mix = D2E_mix_kplus1;
                lambda_k = lambda_kplus1;

                if(_pars.iter %10 == 0){
                    def_geometries.push_back(_problemDOFs.getVertexDOFs());
                    ref_geometries.push_back(_problemDOFs.getReferenceGeometry());
                    fold_DOFs.push_back(_problemDOFs.getFoldDOFs());
                }
                _pars.iter++;
                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
                std::cout<<"Time for iteration: "<<time_span.count()<<std::endl;
            }
        }

        void  BFGS_update(const VectorType &y_k, const VectorType &s_k, FullMatrixType &B_k, FullMatrixType &B_k_inv){

            if(_pars.iter > 0 && (_pars.iter % _BFGS_reset == 0)){
                std::cout<<"Reset BFGS"<<std::endl;
                B_k.setIdentity();
                B_k_inv.setIdentity();
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
                sy = s_k.dot(y_k2);
            }

            B_k += ((y_k2 * y_k2.transpose()) / (sy)) - ((Bs * Bs.transpose()) / (sBs));
            
            // Now update the inverse of B_k
            // We dont operate on complex numbers, thus we know s_k^{t} y_k = y_k^{t} s_k
            FullMatrixType outer_ys = y_k2 * s_k.transpose();
            B_k_inv += (sy + y_k2.dot(B_k_inv*y_k2))*(s_k * s_k.transpose())/(sy*sy) - (B_k_inv*outer_ys + (outer_ys).transpose()*B_k_inv)/(sy);
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
                            const FullMatrixType &C_k,
                            const FullMatrixType &C_k_inv,
                            const FullMatrixType &B_k,
                            MatrixType &D2E_vertex,
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

                lineSearchDOFs = _problemDOFs;
                // Update deformed and reference geometry
                lineSearchDOFs += alpha*d_k_full;

                _costFunctional.apply(lineSearchDOFs.getVertexDOFs(),CostFunctional_val_linesearch);
                
                _factory->produceDE_vertex(lineSearchDOFs,Constraint_linesearch);
                _boundaryDOFs.transformToReducedSpace(Constraint_linesearch);

                constr_l1_linesearch = norm_l1(Constraint_linesearch);
                RealType phi_l1_step = CostFunctional_val_linesearch + mu * constr_l1_linesearch;

                if (phi_l1_step <= (phi_l1 + alpha * _pars.eta * Dp_phi_l1)) {
                    d_k_full *= alpha;
                    lambda_diff *= alpha;
                    return;
                } 
                 // If increase in merit together with increase in constraint use second order correction
                 // Also since the solving of the QP is expensive, only do it after a certain number of iterations
                else if ((constr_l1_linesearch > constr_l1) && (alpha == 1.0) && _secondOrderCorrection)
                {
                    VectorType e_k = C_k.transpose() * (C_k * C_k.transpose()).ldlt().solve(-Constraint_linesearch);
                    
                    
                    _boundaryDOFs.InverseTransformWithFoldDofs(e_k);
                    lineSearchDOFs += e_k;
                    _costFunctional.apply(lineSearchDOFs.getVertexDOFs(),CostFunctional_val_linesearch);
                
                    _factory->produceDE_vertex(lineSearchDOFs,Constraint_linesearch);
                    _boundaryDOFs.transformToReducedSpace(Constraint_linesearch);

                    constr_l1_linesearch = norm_l1(Constraint_linesearch);
                    phi_l1_step = CostFunctional_val_linesearch + mu * constr_l1_linesearch;

                    // subtract 1e-6 for potential numerical inaccuracies
                    if(phi_l1_step <= (phi_l1 + _pars.eta * Dp_phi_l1) + 1e-6) {
                        std::cout<<"Second order correction applied"<<std::endl;
                        d_k_full = d_k_full + e_k;
                        return;
                    }
                    else{
                        alpha *= _pars.tau;
                        continue;
                    }
                }
                else {
                    alpha *= _pars.tau;
                    continue;
                }
            }
            d_k_full *= alpha;
            lambda_diff *= alpha;
        }

    /*
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
    */

        // In the second order correction, dont want to solve for multipliers
        // -> multiplier parameter optional
        void solveQP(size_t nFoldDOFs,
                     size_t nEffectiveVertexDOFs,
                     const FullMatrixType &B_k_inv,
                     const FullMatrixType &C_k,
                     FullMatrixType &mult,
                     FullMatrixType &B_invC,
                     const VectorType &rhs_k, 
                     VectorType &d_k, 
                     VectorType& lambda_kplus1) {

            size_t nEffectiveDOFs = nFoldDOFs + nEffectiveVertexDOFs;

            B_invC = B_k_inv*C_k.transpose();
            mult = -C_k * B_invC;
            VectorType Br1 = B_k_inv*rhs_k.segment(0, nEffectiveDOFs);
            VectorType rhs = C_k*Br1 - rhs_k.segment(nEffectiveDOFs, nEffectiveVertexDOFs);
            VectorType sol(rhs.size());
            sol = _scipy_solver.solve_with_scipy(mult, rhs);
            d_k = Br1 + B_invC*sol;
            lambda_kplus1 = -sol;
        }

        /*
        void solveWithScipy(FullMatrixType &A, Eigen::VectorXd &b, Eigen::VectorXd &sol)
        {
            Py_Initialize();
            if (!Py_IsInitialized()) {
                std::cerr << "Python initialization failed!" << std::endl;
            return;
            }
            npy_intp dims_A[2] = {A.rows(), A.cols()};
            npy_intp dims_b[1] = {b.size()};
               // Prepare dimensions for the NumPy array

            // Create the NumPy array with the same size and data type
            PyObject* A_py = PyArray_SimpleNew(2, dims_A, NPY_DOUBLE);
            if (A_py == NULL) {
                std::cerr << "Error creating NumPy array!" << std::endl;
                Py_Finalize();
                return;
            }

            // Copy the Eigen matrix data into the NumPy array
            std::memcpy(PyArray_DATA((PyArrayObject*)A_py), A.data(), A.size() * sizeof(double));

            PyObject* b_py = PyArray_SimpleNewFromData(1, dims_b, NPY_DOUBLE, b.data());

            // Call scipy.linalg.solve
            PyObject* result = PyObject_CallMethod(_scipy_linalg, "solve", "OO", A_py, b_py);
            if (result == nullptr) {
                PyErr_Print();
                std::cerr << "Error: Could not solve the system." << std::endl;
                Py_Finalize();
                return;
            }

            // Convert result (numpy array) back to a C++ array or Eigen vector
            double* result_data = (double*) PyArray_DATA((PyArrayObject*)result);
            for (int i = 0; i < b.size(); ++i) {
                sol(i) = result_data[i];
            }

            std::cout<<"Error: "<<(A*sol - b).norm()<<std::endl;
            Py_Finalize();
        }
            */

};

