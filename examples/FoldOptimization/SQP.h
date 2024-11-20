#ifndef SQP_H
#define SQP_H

#include <goast/Core.h>
#include "BoundaryDOFS.h"
#include "ObjectFactory.h"

// Basic Parameters that all SQP Solvers should contain
template<typename ConfiguratorType>
class SQPBaseParams{
    public:
        using RealType = typename ConfiguratorType::RealType;
        RealType eps = 1e-6; // Convergence criterion
        size_t maxIter = 10000; // maximum number of iterations
        size_t iter = 0; // current iteration
};

// Specific parameters needed for SQP Line Search
template <typename ConfiguratorType>
class SQPLineSearchParams : public SQPBaseParams<ConfiguratorType>{
    public:
        using RealType = typename ConfiguratorType::RealType;
        static constexpr RealType tau = 0.5;
        static constexpr RealType eta = 0.25;
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

        // pars: Parameters for the SQP Solver
        // plateTopol: Topology of the plate
        // costFunctional: Cost functional
        // DcostFunctional: Derivative of the cost functional
        // factory: Object factory to create energies, their derivatives and hessians on demand
        // std::vector<int> bdryMaskOpt: Boundary mask for the deformed geometry
        // std::vector<int> bdryMaskRef: Boundary mask for the reference geometry
        // FoldDOFs: object to handle degrees of freedom of the fold
        SQPLineSearchSolver(const SQPLineSearchParams<ConfiguratorType> &pars,
                            const MeshTopologySaver &plateTopol,
                            const CostFunctional<ConfiguratorType> &costFunctional,
                            const CostFunctionalGradient<ConfiguratorType> &DcostFunctional,
                            const MyObjectFactory<ConfiguratorType> &factory,
                            std::vector<int> &bdryMaskOpt,
                            std::vector<int> &bdryMaskRef,
                            std::unique_ptr<FoldDofs<ConfiguratorType>> foldDofsPtr,
                            std::unique_ptr<FoldDofsGradient<ConfiguratorType>> DfoldDofsPtr) 
                            : SQPBaseSolver<ConfiguratorType>(pars),
                            _plateTopol(plateTopol),
                            _costFunctional(costFunctional),
                            _DcostFunctional(DcostFunctional),
                            _factory(factory),
                            _bdryMaskOpt(bdryMaskOpt),
                            _bdryMaskRef(bdryMaskRef),
                            _foldDofsPtr(std::move(foldDofsPtr)),
                            _DfoldDofsPtr(std::move(DfoldDofsPtr)) {}

        void solve(VectorType &plateGeomDef, VectorType &plateGeomRef, VectorType &plateGeomInitial)
        {
            size_t nAllVertexDOFs = 3*_plateTopol.getNumVertices();
            size_t nFoldDOFs = _foldDofsPtr->getNumDofs();
            size_t nAllDOFs = nAllVertexDOFs + nFoldDOFs;
            size_t bdryMaskOptSize = _bdryMaskOpt.size();
            size_t nEffectiveDOFs = nAllDOFs - bdryMaskOptSize;
            size_t nEffectiveVertexDOFs = nAllVertexDOFs - bdryMaskOptSize;

            // Initialize the approximate Hessian
            MatrixType B_k(nAllDOFs, nAllDOFs);
            B_k.setZero();

            // get the fold vertices
            std::vector<int> foldVertices;
            _foldDofsPtr->getFoldVertices(foldVertices);

            // get the edge weights
            VectorType edge_weights;
            _foldDofsPtr->getEdgeWeights(edge_weights);

            // Works until here

            // Create all important energies, their derivatives and hessians
            auto DE = _factory.produceDE(_plateTopol, plateGeomRef, edge_weights);
            SimpleBendingGradientDef<DefaultConfigurator> DE_bend( _plateTopol, plateGeomRef , edge_weights);
            auto D2E = _factory.produceD2E(_plateTopol, plateGeomRef, edge_weights);
            auto D2EMixed = _factory.produceD2E_mix(_plateTopol, plateGeomRef, true, false, edge_weights);

            VectorType DE_val;
            MatrixType D2E_val;
            MatrixType D2E_mix_val;

            DE.apply(plateGeomDef, DE_val);
            D2E.apply(plateGeomDef, D2E_val);
            D2EMixed.apply(plateGeomDef, D2E_mix_val);

            // Should work until here

            VectorType DCostFunctional_val;
            _DcostFunctional.apply(plateGeomDef, DCostFunctional_val);

            RealType costFunctional_val;
            _costFunctional.apply(plateGeomDef, costFunctional_val);

            std::cout<<"Initial Cost Functional: "<<costFunctional_val<<std::endl;

            VectorType DFoldDofs_val;
            _DfoldDofsPtr -> apply(plateGeomInitial,DFoldDofs_val); 

            VectorType DE_dt = D2E_mix_val*DFoldDofs_val;
            MatrixType DE_dt_sparse = convertVecToSparseMat(DE_dt);

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

            BoundaryDOFS<ConfiguratorType> boundaryDOFS(_bdryMaskOpt, nAllVertexDOFs, nFoldDOFs);
            
            // Build up B_k matrix in transformed vertex space
            // approximation of Hessian of Lagrangian w.r.t. fold dofs
            MatrixType HessFoldDofs(nFoldDOFs,nFoldDOFs);
            HessFoldDofs.setIdentity();
            MatrixType HessVertexDofs(nAllVertexDOFs, nAllVertexDOFs);
            HessVertexDofs.setIdentity();
            boundaryDOFS.transformRowColToReducedSpace(HessVertexDofs);
            std::vector<Eigen::Triplet<double>> tripletList_B_k;
            assignSparseBlockInplace(B_k, HessFoldDofs, 0, 0, tripletList_B_k);
            assignSparseBlockInplace(B_k, HessVertexDofs, nFoldDOFs, nFoldDOFs, tripletList_B_k);

            // First, reduce the problem to the actual free degrees of freedom
            boundaryDOFS.transformToReducedSpace(DE_val);

            boundaryDOFS.transformRowColToReducedSpace(D2E_val);
            boundaryDOFS.transformToReducedSpace(DE_dt);
            boundaryDOFS.transformWithFoldDofsToReducedSpace(DCostFunctional_val);

            boundaryDOFS.transformRowToReducedSpace(DE_dt_sparse);

            while(_pars.iter < _pars.maxIter &&(DE_val.norm() > _pars.eps || DCostFunctional_val.norm() > _pars.eps)){

                // Solve the QP Problem and write the sol into sol_k
                // solveQPProblem(sol_k, A_k, rhs_k, B_k, D2E_val, DE_dt_sparse, DCostFunctional_val, DE_val, nEffectiveDOFs, nFoldDOFs, nEffectiveVertexDOFs);

                std::cout << "\rIteration: " << std::setw(4) << _pars.iter
                << " | DE_val: " << std::setw(10) << std::fixed << std::setprecision(4) << DE_val.norm()
                << " | DCostFunctional_val: " << std::setw(10) << std::fixed << std::setprecision(4) << DCostFunctional_val.norm()
                << " | CostFunctional_val: " << std::setw(10) << std::fixed << std::setprecision(4) << costFunctional_val
                << std::endl;  // flush to ensure immediate output

                std::vector<Eigen::Triplet<double>> ATriplet;
                assignSparseBlockInplace(A_k, B_k, 0, 0, ATriplet);
                assignSparseBlockInplace(A_k, (-DE_dt_sparse.transpose()).eval(),0,nEffectiveDOFs, ATriplet);
                assignSparseBlockInplace(A_k, (-D2E_val).eval(), nFoldDOFs,nEffectiveDOFs, ATriplet);
                assignSparseBlockInplace(A_k, DE_dt_sparse, nEffectiveDOFs,0, ATriplet);
                assignSparseBlockInplace(A_k, D2E_val, nEffectiveDOFs,nFoldDOFs, ATriplet);

                 if (A_k.isApprox(A_k.transpose())) {
                    std::cout << "The matrix is symmetric.\n";
                } else {
                    std::cout << "The matrix is not symmetric.\n";
                }

                // Construct the r.h.s for the linear system
                // the derivative of the cost functional
                rhs_k.segment(0,nEffectiveDOFs) = -DCostFunctional_val;
                // then the derivative of the energy, i.e. the constraint
                rhs_k.segment(nEffectiveDOFs,nEffectiveVertexDOFs) = -DE_val;

                // Solve the linear system
                Eigen::BiCGSTAB<MatrixType> solver;
                solver.compute(A_k);
                sol_k = solver.solve(rhs_k);

                for(int i = 0; i < sol_k.size(); i++){
                    std::cout<<sol_k[i]<<std::endl;
                }

                d_k = sol_k.segment(0, nEffectiveDOFs);
                lambda_k = sol_k.segment(nEffectiveDOFs, nEffectiveVertexDOFs);

                // Checked until here _____________________________________________

                lambda_diff = lambda_k - lambda_k_prev;
                
                sol_k.segment(nEffectiveDOFs, nEffectiveVertexDOFs) = lambda_diff;

                // Now, determine the stepsizes
                _pars.alpha = 1.0;

                RealType comparisonVal = ((DCostFunctional_val.dot(d_k) + 0.5*d_k.transpose().dot(B_k*d_k))/((1-_pars.rho)*DE_val.template lpNorm<1>() + 1e-6));

                if(_pars.mu < comparisonVal){
                    _pars.mu = (1+1e-4)*comparisonVal; // set mu to a value that is larger than the comparison value
                }
            
                L1Merit<DefaultConfigurator> merit(_pars.mu);
                L1MeritGrad<DefaultConfigurator> Dmerit(_pars.mu);

                RealType merit_x0;
                merit.apply(DE_val,costFunctional_val,merit_x0);

                RealType Dmerit_x0;
                Dmerit.apply(DCostFunctional_val,d_k,DE_val,Dmerit_x0);

                VectorType d_k_vertex_full = d_k.segment(nFoldDOFs,nEffectiveVertexDOFs);
                
                boundaryDOFS.inverseTransform(d_k_vertex_full);

                while(true)
                {

                    // Now, translate the x_k values by alpha*d_k
                    // First, transform the deformed configuration
                    VectorType plateGeomDef_alpha = plateGeomDef + _pars.alpha*d_k_vertex_full;

                    // Next, transform the reference configuration
                    VectorType plateGeomRef_alpha = plateGeomRef;
                    VectorType t = _pars.alpha*d_k.segment(0,nFoldDOFs);

                    _foldDofsPtr->apply(t,plateGeomRef_alpha);

                    // Calculate Energies and constraints for the new configurations
                    auto DE_alpha = _factory.produceDE(_plateTopol, plateGeomRef_alpha, edge_weights);
                    VectorType DE_alpha_val;
                    DE_alpha.apply(plateGeomDef_alpha, DE_alpha_val);
                    boundaryDOFS.transformToReducedSpace(DE_alpha_val);

                    RealType costFunctional_alpha_val;
                    _costFunctional.apply(plateGeomDef_alpha, costFunctional_alpha_val);

                    RealType merit_x_alpha;
                    merit.apply(DE_alpha_val,costFunctional_alpha_val,merit_x_alpha);

                    if(merit_x_alpha <= merit_x0 + _pars.eta*_pars.alpha*Dmerit_x0){
                        break;
                    }

                    _pars.alpha *= _pars.tau;
                }

                // Apply the steps
                plateGeomDef += _pars.alpha*d_k_vertex_full;
                _foldDofsPtr->apply(_pars.alpha*d_k.segment(0,nFoldDOFs),plateGeomRef);
                lambda_k = lambda_k_prev + _pars.alpha*lambda_diff;

                //Evaluate J[x_k+1]
                RealType costFunctional_kplus1_val;
                // can use the same costFunctional object, since foldVertices dont change
                _costFunctional.apply(plateGeomDef, costFunctional_kplus1_val);

                // Evaluate grad J[x_k+1]
                VectorType DCostFunctional_kplus1_val;
                _DcostFunctional.apply(plateGeomDef, DCostFunctional_kplus1_val);
                boundaryDOFS.transformWithFoldDofsToReducedSpace(DCostFunctional_kplus1_val);

                auto DE_kplus1 = _factory.produceDE(_plateTopol, plateGeomRef, edge_weights);

                VectorType DE_kplus1_val;
                DE_kplus1.apply(plateGeomDef, DE_kplus1_val);
                boundaryDOFS.transformToReducedSpace(DE_kplus1_val);

                // Update B_k using BFGS update
                VectorType s_k = _pars.alpha*d_k;
                // Want to calculate \grad_{x} L(x_{k+1}) - \grad_{x} L(x_{k})
                VectorType DiffGradxCostFunctional_kplus1;
                DiffGradxCostFunctional_kplus1 = DCostFunctional_kplus1_val - DCostFunctional_val;

                VectorType DiffGradxEnergy_kplus1(nEffectiveVertexDOFs);
                MatrixType temp(nEffectiveVertexDOFs,nEffectiveDOFs);
                // Recalculate D2E_val

                auto D2E_kplus1 = _factory.produceD2E(_plateTopol, plateGeomRef, edge_weights);
                MatrixType D2E_kplus1_val;
                D2E_kplus1.apply(plateGeomDef, D2E_kplus1_val);
                boundaryDOFS.transformRowColToReducedSpace(D2E_kplus1_val);

                // Recalculate DE/dt:
                auto D2E_mix_kplus1 = _factory.produceD2E_mix(_plateTopol, plateGeomRef, true, false, edge_weights);
                MatrixType D2E_mix_kplus1_val;
                D2E_mix_kplus1.apply(plateGeomDef, D2E_mix_kplus1_val);

                VectorType DFoldDofs_kplus1_val;
                _DfoldDofsPtr -> apply(plateGeomInitial,DFoldDofs_kplus1_val); 
                VectorType DE_dt_kplus_1 = D2E_mix_kplus1_val*DFoldDofs_kplus1_val;
                MatrixType DE_dt_kplus_1_sparse = convertVecToSparseMat(DE_dt_kplus_1);
                boundaryDOFS.transformRowToReducedSpace(DE_dt_kplus_1_sparse);
                
                std::vector<Eigen::Triplet<double>> temp_tripletList;
                assignSparseBlockInplace(temp, DE_dt_kplus_1_sparse, 0,0, temp_tripletList);
                assignSparseBlockInplace(temp, D2E_val, 0,nFoldDOFs, temp_tripletList);

                DiffGradxEnergy_kplus1 = lambda_k.transpose()*(temp);

                VectorType y_k = DiffGradxCostFunctional_kplus1 + DiffGradxEnergy_kplus1;

                FullMatrixType B_k_dense = B_k.toDense();

                BFGS_update(y_k,s_k,B_k_dense);

                B_k = B_k_dense.sparseView();

                if(DE_val.isApprox(DE_kplus1_val,1e-6) && DCostFunctional_val.isApprox(DCostFunctional_kplus1_val,1e-6)){
                    std::cout<<"damnit"<<std::endl;
                }else{
                    std::cout<<"still works though"<<std::endl;
                }

                // After all the calculations, we can update the values
                DE_val = DE_kplus1_val;
                DCostFunctional_val = DCostFunctional_kplus1_val;
                costFunctional_val = costFunctional_kplus1_val;
                D2E_val = D2E_kplus1_val;
                D2E_mix_val = D2E_mix_kplus1_val;
                DE_dt = DE_dt_kplus_1;
                lambda_k_prev = lambda_k;

                _pars.iter++;
            }

        }

/*
        void solveQPProblem(VectorType &sol_k, MatrixType &A_k, const VectorType &rhs_k, const MatrixType &B_k, MatrixType &D2E_val, MatrixType &DE_dt_sparse, VectorType &DCostFunctional_val, VectorType &DE_val, size_t nEffectiveDOFs, size_t nFoldDOFs, size_t nEffectiveVertexDOFs)
        {
            const VectorType neg_DE_val = -DE_val;
            const VectorType neg_DCostFunctional_val = -DCostFunctional_val;

            std::vector<Eigen::Triplet<double>> ATriplet;
            assignSparseBlockInplace(A_k, B_k, 0, 0, ATriplet);
            assignSparseBlockInplace(A_k, (-DE_dt_sparse.transpose()).eval(),0,nEffectiveDOFs, ATriplet);
            assignSparseBlockInplace(A_k, (-D2E_val).eval(), nFoldDOFs,nEffectiveDOFs, ATriplet);
            assignSparseBlockInplace(A_k, DE_dt_sparse, nEffectiveDOFs,0, ATriplet);
            assignSparseBlockInplace(A_k, D2E_val, nEffectiveDOFs,nFoldDOFs, ATriplet);

            // Construct the r.h.s for the linear system
            // the derivative of the cost functional
            rhs_k.block(0,nEffectiveDOFs) = neg_DCostFunctional_val;
            // then the derivative of the energy, i.e. the constraint
            rhs_k.block(nEffectiveDOFs, 0, neg_DE_val.size(), 1) = neg_DE_val;

            LinearSolver<DefaultConfigurator> directSolver;
            directSolver.prepareSolver( A_k );
            directSolver.backSubstitute( rhs_k, sol_k );
        }*/

        void BFGS_update(const VectorType &y_k, const VectorType &s_k, FullMatrixType &B_k){
            RealType sigma = 0;
            VectorType y_k2 = y_k;
            auto Bs = B_k*s_k;
            auto sBs = s_k.dot(Bs);
            auto sy = s_k.dot(y_k);
            if(sy < _pars.theta*sBs){
                sigma = (1-_pars.theta)*sBs/(sBs - sy);
                y_k2 = sigma*y_k + (1-sigma)*Bs;
            }

            B_k.noalias() += y_k2*y_k2.transpose()/(s_k.dot(y_k2)) - Bs*(Bs.transpose())/(sBs);

        }   

        private:
            SQPLineSearchParams<ConfiguratorType> _pars;
            MeshTopologySaver _plateTopol;
            CostFunctional<ConfiguratorType> _costFunctional;
            CostFunctionalGradient<ConfiguratorType> _DcostFunctional;
            MyObjectFactory<ConfiguratorType> _factory;
            std::vector<int> _bdryMaskOpt;
            std::vector<int> _bdryMaskRef;
            std::vector<int> _foldVertices;
            std::unique_ptr<FoldDofs<ConfiguratorType>> _foldDofsPtr;
            std::unique_ptr<FoldDofsGradient<ConfiguratorType>> _DfoldDofsPtr;

};

#endif