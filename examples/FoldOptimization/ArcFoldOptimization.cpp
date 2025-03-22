#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACK

#define OPENBLAS_VERBOSE 1

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <unordered_set>
#include <goast/Smoothers.h>
#include <goast/SQP/DOFHandling/FoldDofs.h>
#include <goast/SQP.h>

#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <typeinfo>
#include <math.h>

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

/** 
 * \brief Optimization of an arc fold as parabola t*(0.25 - (x-0.5)^2) in the reference geometry
 *        The dof t is varied here. Positivity not ensured
 *        In our object factory, we let a force act on the top vertices of the plate
 *        We expect a flipping behaviour like in the bending arc fold examples
 *        Our CostFunctional is the sum of z-components of the top of the plate (y=1.0) 
 * \author Johannssen
 *
 */

bool is_near(RealType coord, RealType value, RealType tol = 1e-6)
{
    return std::abs(coord - value) < tol;
}

int main(int argc, char *argv[])
{

try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "ARC FOLD OPTIMIZATION" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    Eigen::setNbThreads(8);

// load flat plate and prepare the arc crease
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/paperCrissCross.ply");
    MeshTopologySaver plateTopol( plate );
    std::cout<<"num vertices: "<<plateTopol.getNumVertices()<<std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial, plateGeomPenalty; // plateGeomRef is the geometry without the arc crease in the mesh
    getGeometry(plate,plateGeomRef);
    getGeometry(plate,plateGeomDef);
    getGeometry(plate,plateGeomInitial);

    // Fix all coordinates of the boundary
    std::vector<int> bdryMaskRef_all;

    RealType t_0 = 0.25; // Initial value of the parameter t

    // Select the vertices with y = 1.0 (topVertices)
    // We will make a Gravitational Force act on them in order to achieve a flipping effect
    std::vector<int> topVertices;
    std::vector<int> upperVertices;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
        RealType x,y,z;
        x = coords[0];
        y = coords[1];
        z = coords[2];

        if(y > 0.25){
            upperVertices.push_back(i);
        }

        if(is_near(y, 1.0) || is_near(y, 0.0) || is_near(x, 0.0) || is_near(x, 1.0)){
            bdryMaskRef_all.push_back(i);
            if(is_near(y,1.0))
            {
                topVertices.push_back(i);
            }
            continue;
        }

        if(is_near(y,0.25)){
            bdryMaskRef_all.push_back(i);
            coords[1] += t_0*(0.25 - std::pow(coords[0] - 0.5, 2.0));
            setXYZCoord<VectorType, VecType>(plateGeomRef, coords, i);
            continue;
        }
    }

    // output preliminary meshes
    setGeometry(plate, plateGeomRef);
    OpenMesh::IO::write_mesh(plate, "reference_mesh_before_smoothing.ply");

    // ----------------- Smooth reference geometry -----------------------------

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef_all );
    //extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_xz , activeRef_xz);
    
    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntriesRef(bdryMaskRef_all.begin(), bdryMaskRef_all.end());
    //uniqueEntriesRef.insert(bdryMaskRef_xz.begin(), bdryMaskRef_xz.end());

    // Move the unique elements back to bdryMaskRef in order
    std::vector<int> bdryMaskRef = bdryMaskRef_all;

    DirichletSmoother<DefaultConfigurator> smoother(plateGeomInitial, bdryMaskRef, plateTopol);
    smoother.apply(plateGeomRef,plateGeomRef);

    setGeometry(plate,plateGeomRef);
    OpenMesh::IO::write_mesh(plate,"reference_mesh.ply");

    getGeometry(plate,plateGeomDef);

    // Fix all coordinates
    std::vector<int> bdryMaskOpt_all;
    std::vector<int> bdryMaskOptSmooth_all;
    // Fix only (x,y) coordinates
    std::vector<int> bdryMaskOptSmooth_xy;
    // Fix only y coordinates
    std::vector<int> bdryMaskOptSmooth_y;
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>(plateGeomInitial,coords,i);
        RealType x,y,z;
        x = coords[0];
        y = coords[1];
        z = coords[2];

        if(is_near(x,0.0) && (y <= 0.25 + 1e-4))
        {
            bdryMaskOptSmooth_all.push_back(i);
            bdryMaskOpt_all.push_back(i);
            coords[2] += 0.1;
            coords[0] += 0.1;
            setXYZCoord<VectorType, VecType>(plateGeomDef,coords,i);
            continue;
        }

        if(is_near(x,1.0) && (y <= 0.25 + 1e-4))
        {
            bdryMaskOptSmooth_all.push_back(i);
            bdryMaskOpt_all.push_back(i);
            coords[2] += 0.1;
            coords[0] -= 0.1;
            setXYZCoord<VectorType, VecType>(plateGeomDef,coords,i);
            continue;
        }

        if(y >= 0.4)
        {
            bdryMaskOptSmooth_xy.push_back(i);
        }

        if(is_near(y,0.0))
        {
            bdryMaskOptSmooth_y.push_back(i);
        }
    }

    std::vector<int> active_xy = (std::vector<int>){1,1,0};
    std::vector<int> active_y = (std::vector<int>){0,1,0};
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOptSmooth_all );
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskOptSmooth_xy, active_xy);
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskOptSmooth_y, active_y);

    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntriesOptSmooth(bdryMaskOptSmooth_all.begin(), bdryMaskOptSmooth_all.end());
    uniqueEntriesOptSmooth.insert(bdryMaskOptSmooth_xy.begin(), bdryMaskOptSmooth_xy.end());
    uniqueEntriesOptSmooth.insert(bdryMaskOptSmooth_y.begin(), bdryMaskOptSmooth_y.end());

    // Move the unique elements back to bdryMaskRef in order
    std::vector<int> bdryMaskOptSmooth;
    bdryMaskOptSmooth.assign(uniqueEntriesOptSmooth.begin(), uniqueEntriesOptSmooth.end());
    std::sort(bdryMaskOptSmooth.begin(), bdryMaskOptSmooth.end());

    setGeometry(plate, plateGeomDef);
    OpenMesh::IO::write_mesh(plate, "deformed_mesh_before_smoothing.ply");

    DirichletSmoother<DefaultConfigurator> DefSmoother(plateGeomInitial,bdryMaskOptSmooth,plateTopol);
    DefSmoother.apply(plateGeomDef, plateGeomDef);

    setGeometry(plate,plateGeomDef);
    OpenMesh::IO::write_mesh(plate,"deformed_mesh.ply");

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt_all );

    std::unordered_set<int> uniqueEntriesOpt(bdryMaskOpt_all.begin(), bdryMaskOpt_all.end());
    //uniqueEntriesOpt.insert(bdryMaskOpt_2.begin(), bdryMaskOpt_2.end());

    std::vector<int> bdryMaskOpt;
    bdryMaskOpt.assign(uniqueEntriesOpt.begin(), uniqueEntriesOpt.end());
    std::sort(bdryMaskOpt.begin(), bdryMaskOpt.end());
   
    auto foldDofsPtr = std::make_shared<FoldDofsArcLine<DefaultConfigurator>>(plateTopol,plateGeomInitial, plateGeomInitial, bdryMaskRef);

    std::vector<int> foldVertices;
    foldDofsPtr -> getFoldVertices(foldVertices);

    auto DfoldDofsPtr = std::make_shared<FoldDofsArcLineGradient<DefaultConfigurator>>(plateTopol, bdryMaskRef, plateGeomInitial, foldVertices);

    VectorType edge_weights = VectorType::Zero(plateTopol.getNumEdges());
    foldDofsPtr->getEdgeWeights(edge_weights);

    size_t nFoldDOFs = foldDofsPtr->getNumDofs();
    size_t nVertexDOFs = 3*plateTopol.getNumVertices();

    std::cout<<"Number of fold DOFs: "<<nFoldDOFs<<std::endl;
    
    SQPLineSearchParams<DefaultConfigurator> pars;
    CostFunctional<DefaultConfigurator> costFunctional(topVertices);
    CostFunctionalGradient<DefaultConfigurator> DcostFunctional(plateTopol,topVertices);
    
    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;
    RealType factor_gravity = 0.0;

    VectorType factors(3);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;
    factors[2] = factor_gravity;

    VectorType mass_distribution = VectorType::Zero(plateTopol.getNumVertices());

    for (int idx : upperVertices) {
        mass_distribution[idx] = 1;
    }

    VectorType gravity_dir(3);
    gravity_dir << 0.0, 0.0, -1.0;

    /*
    // Initialize with Newton deformation:
    SimpleBendingEnergy<DefaultConfigurator> E_bend( plateTopol, plateGeomRef, true , edge_weights);
    SimpleBendingGradientDef<DefaultConfigurator> DE_bend( plateTopol, plateGeomRef , edge_weights);
    SimpleBendingHessianDef<DefaultConfigurator> D2E_bend( plateTopol, plateGeomRef , edge_weights);

    typename DefaultConfigurator::RealType energy;

    NonlinearMembraneEnergy<DefaultConfigurator> E_mem( plateTopol, plateGeomRef, true );
    NonlinearMembraneGradientDef<DefaultConfigurator> DE_mem( plateTopol, plateGeomRef );
    NonlinearMembraneHessianDef<DefaultConfigurator> D2E_mem( plateTopol, plateGeomRef );

    GravitationalEnergy<DefaultConfigurator> E_grav( plateTopol, plateGeomRef, true, mass_distribution, gravity_dir );
    GravitationalEnergyGradientDef<DefaultConfigurator> DE_grav( plateTopol, plateGeomRef, mass_distribution, gravity_dir );

    AdditionOp<DefaultConfigurator> E_tot( factors, E_mem, E_bend, E_grav);
    AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_mem, DE_bend, DE_grav);
    AdditionHessian<DefaultConfigurator> D2E_tot(factors, D2E_mem, D2E_bend);

    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setNewtonIterations( 1000 );
    optPars.setQuietMode( SHOW_ALL );
    VectorType initialization = plateGeomDef;

    std::cerr<< "Startting Newton Linesearch..."<<std::endl;
    LineSearchNewton<DefaultConfigurator> NLS( E_tot, DE_tot, D2E_tot, optPars);
    NLS.setBoundaryMask( bdryMaskOpt );
    NLS.solve( initialization, plateGeomDef );
    

    setGeometry(plate, plateGeomDef);
    OpenMesh::IO::write_mesh(plate, "deformed_mesh_newton.ply");
    */

    // Instead, try to open the newton version
    OpenMesh::IO::read_mesh(plate, "bendingFoldSol_withNewton_scaled.ply");
    getGeometry(plate, plateGeomDef);

    /*
    // Now, test the vector valued derivatives of the penalty cost functional
    // First, specify the penalty vertices
    std::pair<std::vector<int>, std::vector<int>> penalty_vertices_1;
    std::pair<std::vector<int>, std::vector<int>> penalty_vertices_2;

    // Calculate the resolution in x direction by calculating the x values closest to zero
    RealType x_nearest = 1.0;
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>(plateGeomInitial, coords, i);
        RealType x,y,z;
        x = coords[0];
        y = coords[1];
        z = coords[2];
        if(x < x_nearest && (!is_near(x,0.0)))
        {
            x_nearest = x;
        }
    }

    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>(plateGeomInitial, coords, i);
        RealType x,y,z;
        x = coords[0];
        y = coords[1];
        z = coords[2];

        if(y <= 0.25 + 1e-4)
        {
            if(is_near(x,0.0))
            {
                penalty_vertices_1.first.push_back(i);
                continue;
            }
            if(is_near(x,x_nearest))
            {
                penalty_vertices_1.second.push_back(i);
                continue;
            }
            if(is_near(x,1.0))
            {
                penalty_vertices_2.first.push_back(i);
                continue;
            }
            if(is_near(x,1-x_nearest))
            {
                penalty_vertices_2.second.push_back(i);
                continue;
            }
        }
    }

    OpenMesh::IO::read_mesh(plate, "../../data/plate/testPenalty.ply");
    getGeometry(plate, plateGeomPenalty);

    PenaltyCostFunctional<DefaultConfigurator> penaltyCostFunctional(topVertices, penalty_vertices_1, penalty_vertices_2);
    PenaltyCostFunctionalGradient<DefaultConfigurator> DpenaltyCostFunctional(plateTopol, topVertices, penalty_vertices_1, penalty_vertices_2);
    ScalarValuedDerivativeTester<DefaultConfigurator> tester(penaltyCostFunctional, DpenaltyCostFunctional, 0.005, plateGeomPenalty.size());
    tester.plotAllDirections(plateGeomPenalty, "deriv_test/test");
    // Now start with the SQP iterations

    */

    std::vector<VectorType> def_geometries, ref_geometries, fold_DOFs;
    auto factory = std::make_shared<ElasticGravitationalFactory<DefaultConfigurator>>(factors, plateTopol, edge_weights, gravity_dir, mass_distribution);
    BoundaryDOFS<DefaultConfigurator> boundaryDOFs(bdryMaskOpt, nVertexDOFs, nFoldDOFs);
    // Create the degrees of freedom object
    std::vector<RealType> deviations;
    ProblemDOFs<DefaultConfigurator> problemDOFs(VectorType::Ones(1)*t_0, plateGeomDef, foldDofsPtr, DfoldDofsPtr);
    SQPLineSearchSolver<DefaultConfigurator> solver(pars, costFunctional, DcostFunctional, std::move(factory), boundaryDOFs, problemDOFs, 20);
    solver.solve(plateGeomRef, def_geometries, ref_geometries, fold_DOFs);
    std::string filename;

    for(int i = 0; i < def_geometries.size(); i++){
        filename = "deformed/plate_" + std::to_string(i) + ".ply";
        setGeometry(plate, def_geometries[i]);
        OpenMesh::IO::write_mesh(plate,filename);
        setGeometry(plate, ref_geometries[i]);
        filename = "reference/plate_" + std::to_string(i) + ".ply";
        OpenMesh::IO::write_mesh(plate, filename);
    }

  }catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }

    return 0;
}