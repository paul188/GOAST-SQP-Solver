/* 
 Simulating a fold that results from bending a plate with clamped 
 boundary conditions. Only NonlinearMembraneEnergy and SimpleBendingEnergy are used.
 The gravitational energy is not used in this example.
 We set the weight of the contribution of the middle edge to the bending energy to 0.
 We thus expect a lot of the bending to happen in the middle of the plate, at the fold.
 Imagine a piece of paper that is folded in the middle in both directions.
 So that the fold does not point into any direction, but the edge just doesnt
 resist bending.
*/

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <iostream>
#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include "Arena.h"
#include <goast/Smoothers.h>
#include "CostFunctional.h"
#include "FoldDofs.h"
#include "Merit.h"
#include "SQP.h"
#include "SparseMat.h"
//#include "BoundaryDOFS.h"
#include <chrono>
#include <thread>

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

/** 
 * \brief Simulation of the optimization of a simple fold
 * \author Johannssen
 *
 * The cost functional is defined in CostFunctional.h
 * Arena.h is used to manage memory allocation
 * 
 */

/**/

int main(int argc, char *argv[])
{

  using MatrixType = DefaultConfigurator::SparseMatrixType;
  using FullMatrixType = DefaultConfigurator::FullMatrixType;

try{
    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE FOLD OPTIMIZATION" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    std::cout << "Eigen version: "
              << EIGEN_WORLD_VERSION << "."
              << EIGEN_MAJOR_VERSION << "."
              << EIGEN_MINOR_VERSION << std::endl;

    // load flat plate [0,1]^2
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/plate4SD.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial;
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );
    getGeometry( plate, plateGeomInitial );

    Arena arena;

    // determine boundary mask for reference geometry and foldVertices
    // the part of the reference boundary where everything is fixed
    std::vector<int> bdryMaskRef_1;
    // the part of the reference boundary where only y,z coordinates are fixed -> x free
    std::vector<int> bdryMaskRef_2;
    std::vector<int> foldVertices;
    VectorType edge_weights = VectorType::Ones(plateTopol.getNumEdges());
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
        if( coords[0] == 0.0 || coords[0] == 1.0 ){
            bdryMaskRef_1.push_back( i );
        }
        if( coords[1] == 0.0 || coords[1] == 1.0 ){
            bdryMaskRef_2.push_back( i );
        }
        // Choose the fold vertices at an imperfect location -> so we can see a change
        if( coords[0] == 0.625 ){
            foldVertices.push_back( i );
            edge_weights[i] = 0;
        }
    }

     // determine boundary mask for optimization
     // and deform part of boundary
    std::vector<int> bdryMaskOpt;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        if( coords[0] == 0.0 || coords[0] == 1.0 ){
            bdryMaskOpt.push_back( i );
        // deform part of boundary
        }
        
        coords[0] *= 0.8;
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
    }

    OpenMesh::IO::write_mesh(plate, "testPlate1.ply");

    // extend all boundary masks to (x,y,z) coordinates

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef_1 );
    std::vector<int> activeRef_2 = (std::vector<int>){0,1,1};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_2 , activeRef_2);
    // append to the Dirichlet bdry mask of Reference boundary 1
    bdryMaskRef_1.insert(bdryMaskRef_1.end(), bdryMaskRef_2.begin(), bdryMaskRef_2.end());
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );

    // prepare the initial deformed geometry
    DirichletSmoother<DefaultConfigurator> smoother(plateGeomRef, plateGeomDef, bdryMaskRef_1, plateTopol);
    smoother.smoothMesh(plateGeomDef);

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "testPlate2.ply");

    //BoundaryDOFS<DefaultConfigurator> boundaryDOFS(bdryMaskOpt, plateTopol);

    // we have 3n vertices and we differentiate after one param t -> (3n +1) x (3n +1)
    FullMatrixType B_k(3*plateTopol.getNumVertices()+1, 3*plateTopol.getNumVertices()+1);
    B_k.setIdentity();
    // Create all important energies, their derivatives and hessians and reserve the memory via arena
    // This is prudent, since the members are constant, thus we need to create new energy objects
    // in each iteration of the optimization. This way, we can use the reserved memory in the arena
    SimpleBendingGradientDef<DefaultConfigurator>* DE_bend = new(arena) SimpleBendingGradientDef<DefaultConfigurator>(plateTopol, plateGeomRef, true);
    SimpleBendingHessianDef<DefaultConfigurator>* D2E_bend = new(arena) SimpleBendingHessianDef<DefaultConfigurator>(plateTopol, plateGeomRef, true);
    SimpleBendingHessianMixed<DefaultConfigurator>* D2E_bend_mix = new(arena) SimpleBendingHessianMixed<DefaultConfigurator>(plateTopol, plateGeomRef, true, false);

    NonlinearMembraneGradientDef<DefaultConfigurator>* DE_membrane = new(arena) NonlinearMembraneGradientDef<DefaultConfigurator>(plateTopol, plateGeomRef, true);
    NonlinearMembraneHessianDef<DefaultConfigurator>* D2E_membrane = new(arena) NonlinearMembraneHessianDef<DefaultConfigurator>(plateTopol, plateGeomRef, true);
    NonlinearMembraneHessianMixed<DefaultConfigurator>* D2E_membrane_mix = new(arena) NonlinearMembraneHessianMixed<DefaultConfigurator>(plateTopol, plateGeomRef, true,false);

    VectorType DE_bend_val;
    VectorType DE_membrane_val;
    VectorType DE_val;
    DE_bend->apply(plateGeomDef,DE_bend_val);
    DE_membrane->apply(plateGeomDef,DE_membrane_val);
    DE_val = DE_bend_val + DE_membrane_val;

    MatrixType D2E_bend_val;
    MatrixType D2E_membrane_val;
    MatrixType D2E_val;
    D2E_bend->apply(plateGeomDef,D2E_bend_val);
    D2E_membrane->apply(plateGeomDef,D2E_membrane_val);
    D2E_val = D2E_bend_val + D2E_membrane_val;

    MatrixType D2E_bend_mix_val;
    MatrixType D2E_membrane_mix_val;
    MatrixType D2E_mix_val;
    D2E_bend_mix->apply(plateGeomDef,D2E_bend_mix_val);
    D2E_membrane_mix->apply(plateGeomDef,D2E_membrane_mix_val);
    D2E_mix_val = D2E_bend_mix_val + D2E_membrane_mix_val;

    // Create the cost functional
    CostFunctional<DefaultConfigurator> costFunctional(foldVertices);
    CostFunctionalGradient<DefaultConfigurator> DcostFunctional(plateTopol,foldVertices);

    VectorType DCostFunctional_val;
    DcostFunctional.apply(plateGeomDef, DCostFunctional_val);

    RealType costFunctional_val;
    costFunctional.apply(plateGeomDef, costFunctional_val);

    // Initialize Fold Translation Object and Gradient
    FoldDofs<DefaultConfigurator> foldDofs(plateGeomRef,foldVertices);
    FoldDofsGradient<DefaultConfigurator> DfoldDofs(plateTopol,foldVertices);
    VectorType DFoldDofs_val;
    DfoldDofs.apply(plateGeomInitial,DFoldDofs_val); 

    VectorType DE_dt = D2E_mix_val*DFoldDofs_val;

    // Reserve memory for the large matrix and r.h.s for the linear system
    FullMatrixType A_k(6*plateTopol.getNumVertices()+1, 6*plateTopol.getNumVertices()+1);
    A_k.setZero();
    VectorType rhs_k = VectorType::Zero(6*plateTopol.getNumVertices()+1);

     // Create solution vector for linear system that contains as subvector the Lagrange multipliers
    VectorType sol_k = VectorType::Zero(6*plateTopol.getNumVertices()+1);

    // Create the subvectors of the solution vector
    VectorType::SegmentReturnType d_k = sol_k.segment(0, 3*plateTopol.getNumVertices()+1);
    VectorType::SegmentReturnType lambda_k = sol_k.segment(3*plateTopol.getNumVertices()+1, 3*plateTopol.getNumVertices());
    VectorType lambda_k_prev = VectorType::Zero(3*plateTopol.getNumVertices());
    VectorType lambda_diff = VectorType::Zero(3*plateTopol.getNumVertices());

    SQP_Parameters<DefaultConfigurator> pars;

    RealType t = 0;

    // In the following, we will always order the derivatives in the following way: (t,\phi)
    // that means in rows and columns, the derivative after the DOF t comes first

    while(pars.iter < pars.maxiter &&(DE_val.norm() > pars.eps || DCostFunctional_val.norm() > pars.eps)){

        t += d_k[0];

        std::cout<<"Iteration "<<pars.iter<<std::endl;
        std::cout<<"Norm1: "<<DE_val.norm()<<std::endl;
        std::cout<<"Norm2: "<<costFunctional_val<<std::endl;
        std::cout<<"t: "<<t<<std::endl;

        setGeometry(plate,  plateGeomDef);
        OpenMesh::IO::write_mesh(plate, "plateGeomDef" + std::to_string(pars.iter) + ".ply");
        setGeometry(plate, plateGeomRef);
        OpenMesh::IO::write_mesh(plate, "plateGeomRef" + std::to_string(pars.iter) + ".ply");

        auto start = std::chrono::high_resolution_clock::now();

        A_k.block(0,0,3*plateTopol.getNumVertices()+1,3*plateTopol.getNumVertices()+1) = B_k;
        
        // Now, the other two submatrices
        A_k.block(0,3*plateTopol.getNumVertices()+1,1,3*plateTopol.getNumVertices()) = -DE_dt.transpose();
        A_k.block(1,3*plateTopol.getNumVertices()+1,3*plateTopol.getNumVertices(),3*plateTopol.getNumVertices()) = -D2E_val;
        A_k.block(1+3*plateTopol.getNumVertices(),0,3*plateTopol.getNumVertices(),1) = DE_dt;
        A_k.block(1+3*plateTopol.getNumVertices(),1,3*plateTopol.getNumVertices(),3*plateTopol.getNumVertices()) = D2E_val;

        // Construct the r.h.s for the linear system
        // the derivative of the cost functional
        rhs_k.segment(0,3*plateTopol.getNumVertices()+1) = -DCostFunctional_val;
        // then the derivative of the energy, i.e. the constraint
        rhs_k.segment(1+3*plateTopol.getNumVertices(),3*plateTopol.getNumVertices()) = -DE_val;

        MatrixType A_k_sparse = A_k.sparseView(1e-8);

        LinearSolver<DefaultConfigurator> directSolver;
        directSolver.prepareSolver( A_k_sparse );
        directSolver.backSubstitute( rhs_k, sol_k );
        auto end = std::chrono::high_resolution_clock::now();

        // Calculate duration in milliseconds
        std::chrono::duration<double, std::milli> duration = end - start;
        std::cout << "Time taken for solving the linear system in iteration "<<pars.iter<<": " << duration.count() << " ms" << std::endl;

        lambda_diff = lambda_k - lambda_k_prev;
        lambda_k_prev = lambda_k;
        sol_k.segment(3*plateTopol.getNumVertices()+1, 3*plateTopol.getNumVertices()) = lambda_diff;

        // Now, determine the stepsizes
        pars.alpha = 1.0;
        pars.mu = 1.0;

        RealType comparisonVal = ((DCostFunctional_val.dot(d_k) + 0.5*d_k.transpose().dot(B_k*d_k))/((1-pars.rho)*DE_val.template lpNorm<1>() + 1e-6));

        if(pars.mu < comparisonVal){
            pars.mu = (1+1e-2)*comparisonVal; // set mu to a value that is larger than the comparison value
        }
       
        L1Merit<DefaultConfigurator> merit(pars.mu);
        L1MeritGrad<DefaultConfigurator> Dmerit(pars.mu);

        RealType merit_x0;
        merit.apply(DE_val,costFunctional_val,merit_x0);

        RealType Dmerit_x0;
        Dmerit.apply(DCostFunctional_val,d_k,DE_val,Dmerit_x0);

        auto start_stepsize = std::chrono::high_resolution_clock::now();

        while(true)
        {

            // Now, translate the x_k values by alpha*d_k
            // First, transform the deformed configuration
            VectorType plateGeomDef_alpha = plateGeomDef + pars.alpha*d_k.segment(1,3*plateTopol.getNumVertices());

            // Next, transform the reference configuration
            VectorType plateGeomRef_alpha;
            VectorType t = pars.alpha*d_k.segment(0,1);
            foldDofs.apply(t,plateGeomRef_alpha);

            // Calculate Energies and constraints for the new configurations
            SimpleBendingGradientDef<DefaultConfigurator> E_bend_alpha(plateTopol, plateGeomRef_alpha, true);
            NonlinearMembraneGradientDef<DefaultConfigurator> E_membrane_alpha(plateTopol, plateGeomRef_alpha, true);

            VectorType DE_bend_alpha_val;
            VectorType DE_membrane_alpha_val;
            VectorType DE_alpha_val;
            E_bend_alpha.apply(plateGeomDef_alpha,DE_bend_alpha_val);
            E_membrane_alpha.apply(plateGeomDef_alpha,DE_membrane_alpha_val);
            DE_alpha_val = DE_bend_alpha_val + DE_membrane_alpha_val;

            CostFunctional<DefaultConfigurator> costFunctional_alpha(foldVertices);

            RealType costFunctional_alpha_val;
            costFunctional_alpha.apply(plateGeomDef_alpha, costFunctional_alpha_val);

            RealType merit_x_alpha;
            merit.apply(DE_alpha_val,costFunctional_alpha_val,merit_x_alpha);

            if(merit_x_alpha <= merit_x0 + pars.eta*pars.alpha*Dmerit_x0){
                break;
            }
            pars.alpha *= pars.tau;
        }

        // Apply the steps
        plateGeomDef += pars.alpha*d_k.segment(1,3*plateTopol.getNumVertices());
        foldDofs.apply(pars.alpha*d_k.segment(0,1),plateGeomRef);
        lambda_k_prev = lambda_k;
        lambda_k += pars.alpha*lambda_diff;

        //Evaluate J[x_k+1]
        RealType costFunctional_kplus1_val;
        // can use the same costFunctional object, since foldVertices dont change
        costFunctional.apply(plateGeomDef, costFunctional_kplus1_val);

        // Evaluate grad J[x_k+1]
        VectorType DCostFunctional_kplus1_val;
        DcostFunctional.apply(plateGeomDef, DCostFunctional_kplus1_val);

        SimpleBendingGradientDef<DefaultConfigurator> DE_bend_kplus1(plateTopol, plateGeomRef, true);
        NonlinearMembraneGradientDef<DefaultConfigurator> DE_membrane_kplus1(plateTopol, plateGeomRef, true);

        VectorType DE_bend_kplus1_val;
        VectorType DE_membrane_kplus1_val;
        VectorType DE_kplus1_val;
        DE_bend_kplus1.apply(plateGeomDef,DE_bend_kplus1_val);
        DE_membrane_kplus1.apply(plateGeomDef,DE_membrane_kplus1_val);
        DE_kplus1_val.resize(DE_membrane_kplus1_val.size());
        DE_kplus1_val = DE_bend_kplus1_val + DE_membrane_kplus1_val;

        // Update B_k using BFGS update
        VectorType s_k = pars.alpha*d_k;
        // Want to calculate \grad_{x} L(x_{k+1}) - \grad_{x} L(x_{k})
        VectorType DiffGradxCostFunctional_kplus1(3*plateTopol.getNumVertices()+1);
        DiffGradxCostFunctional_kplus1 = DCostFunctional_kplus1_val - DCostFunctional_val;

        VectorType DiffGradxEnergy_kplus1(3*plateTopol.getNumVertices());
        FullMatrixType temp(3*plateTopol.getNumVertices(),3*plateTopol.getNumVertices()+1);
        // Recalculate D2E_val

        SimpleBendingHessian<DefaultConfigurator> D2E_bend_kplus1(plateTopol, plateGeomRef, true);
        NonlinearMembraneHessian<DefaultConfigurator> D2E_membrane_kplus1(plateTopol, plateGeomRef, true);
        MatrixType D2E_bend_kplus1_val;
        MatrixType D2E_membrane_kplus1_val;
        MatrixType D2E_kplus1_val;
        D2E_bend_kplus1.apply(plateGeomDef,D2E_bend_kplus1_val);
        D2E_membrane_kplus1.apply(plateGeomDef,D2E_membrane_kplus1_val);
        D2E_kplus1_val = D2E_bend_kplus1_val + D2E_membrane_kplus1_val;

        // Recalculate DE/dt:
        SimpleBendingHessianMixed<DefaultConfigurator> D2E_bend_mix_kplus1(plateTopol, plateGeomRef, true, false);
        NonlinearMembraneHessianMixed<DefaultConfigurator> D2E_membrane_mix_kplus1(plateTopol, plateGeomRef, true, false);
        MatrixType D2E_bend_mix_kplus1_val;
        MatrixType D2E_membrane_mix_kplus1_val;
        MatrixType D2E_mix_kplus1_val;
        D2E_bend_mix_kplus1.apply(plateGeomDef,D2E_bend_mix_kplus1_val);
        D2E_membrane_mix_kplus1.apply(plateGeomDef,D2E_membrane_mix_kplus1_val);
        D2E_mix_kplus1_val = D2E_bend_mix_kplus1_val + D2E_membrane_mix_kplus1_val;

        FoldDofsGradient<DefaultConfigurator> DfoldDofs_kplus1(plateTopol,foldVertices);
        VectorType DFoldDofs_kplus1_val;
        DfoldDofs_kplus1.apply(plateGeomInitial,DFoldDofs_kplus1_val); 
        VectorType DE_dt_kplus_1 = D2E_mix_kplus1_val*DFoldDofs_kplus1_val;

        temp.block(0,0,3*plateTopol.getNumVertices(),1) = DE_dt_kplus_1;
        temp.block(0,1,3*plateTopol.getNumVertices(), 3*plateTopol.getNumVertices()) = D2E_val;

        DiffGradxEnergy_kplus1 = lambda_k.transpose()*(temp);

        VectorType y_k = DiffGradxCostFunctional_kplus1 + DiffGradxEnergy_kplus1;

        BFGS_update<DefaultConfigurator>(y_k,s_k,pars.theta,B_k);

        // After all the calculations, we can update the values
        DE_val = DE_kplus1_val;
        DCostFunctional_val = DCostFunctional_kplus1_val;
        costFunctional_val = costFunctional_kplus1_val;
        D2E_val = D2E_kplus1_val;
        D2E_mix_val = D2E_mix_kplus1_val;
        DE_dt = DE_dt_kplus_1;

        std::cout<<"COST FUNCTIONAL VALUE: "<<costFunctional_val<<std::endl;

        pars.iter++;
    }

    setGeometry(plate,plateGeomDef);
    OpenMesh::IO::write_mesh(plate, "result_foldOptim.ply");

    setGeometry(plate, plateGeomRef);
    OpenMesh::IO::write_mesh(plate,"ReferenceGeomOptim.ply");

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}