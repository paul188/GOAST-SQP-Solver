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
    OpenMesh::IO::read_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/paperCrissCross.ply");
    MeshTopologySaver plateTopol( plate );
    std::cout<<"num vertices: "<<plateTopol.getNumVertices()<<std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial; // plateGeomRef is the geometry without the arc crease in the mesh
    getGeometry(plate,plateGeomRef);
    getGeometry(plate,plateGeomDef);
    getGeometry(plate,plateGeomInitial);

    // Fix all coordinates of the boundary
    std::vector<int> bdryMaskRef;

    RealType t_0 = 0.4;//0.8084239959716796; // Initial value of the parameter t
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
            bdryMaskRef.push_back(i);
            if(is_near(y,1.0))
            {
                topVertices.push_back(i);
            }
            continue;
        }

        if(is_near(y,0.25)){
            bdryMaskRef.push_back(i);
            continue;
        }
    }

    // ----------------- Smooth reference geometry -----------------------------

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef );

    // Fix all coordinates
    std::vector<int> bdryMaskOpt;
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>(plateGeomInitial,coords,i);
        RealType x,y,z;
        x = coords[0];
        y = coords[1];
        z = coords[2];

        if( y <= 0.25 + 1e-4 && (is_near(x,0.0) || is_near(x,1.0)) ){
            bdryMaskOpt.push_back( i );
        }
    }

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );
   
    auto foldDofsPtr = std::make_shared<FoldDofsArcLine<DefaultConfigurator>>(plateTopol,plateGeomInitial, plateGeomInitial, bdryMaskRef);

    std::vector<int> foldVertices;
    foldDofsPtr -> getFoldVertices(foldVertices);

    std::cout<<"Num fold vertices: "<<foldVertices.size()<<std::endl;

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

    // Slightly inaccurate loading of the mesh
    //OpenMesh::IO::read_mesh(plate, "bendingFoldSol_withNewton_scaled_2.ply");
    //getGeometry(plate, plateGeomDef);

    // this is the accurate way
    readVectorFromFile(plateGeomDef, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/bendingFoldSol_withNewton_scaled_2.txt");

    std::vector<VectorType> def_geometries, ref_geometries, fold_DOFs;
    auto factory = std::make_shared<ElasticGravitationalFactory<DefaultConfigurator>>(factors, plateTopol, edge_weights, gravity_dir, mass_distribution);
    BoundaryDOFS<DefaultConfigurator> boundaryDOFs(bdryMaskOpt, nVertexDOFs, nFoldDOFs);
    // Create the degrees of freedom object
    std::vector<RealType> deviations;
    VectorType vertexDOFs_initial = VectorType::Zero(3*plateTopol.getNumVertices());
    ProblemDOFs<DefaultConfigurator> problemDOFs(VectorType::Ones(1)*t_0, vertexDOFs_initial, plateGeomDef, foldDofsPtr, DfoldDofsPtr);
    SQPLineSearchSolver<DefaultConfigurator> solver(pars, costFunctional, DcostFunctional, factory, boundaryDOFs, problemDOFs, 20);
    
    // Test the constraint norm:
    ProblemDOFs<DefaultConfigurator> problemDOFsTest(VectorType::Ones(1)*t_0, vertexDOFs_initial, plateGeomDef, foldDofsPtr, DfoldDofsPtr);
    VectorType DE_test;
    factory->produceDE_vertex(problemDOFsTest, DE_test);
    applyMaskToVector(bdryMaskOpt, DE_test);
    RealType constraint_norm = DE_test.norm();
    std::cout<<"Constraint norm: "<<constraint_norm<<std::endl;

    const std::string filename_def_geometries = "/lustre/scratch/data/s24pjoha_hpc-results/thesis_results/deformed_ArcFoldOptimization/";
    const std::string filename_ref_geometries = "/lustre/scratch/data/s24pjoha_hpc-results/thesis_results/reference_ArcFoldOptimization/";

    solver.solve(plateGeomInitial, filename_def_geometries, filename_ref_geometries, plate);
    std::string filename;

    /*
    for(int i = 0; i < def_geometries.size(); i++){
        filename = "deformed/plate_" + std::to_string(i) + ".ply";
        setGeometry(plate, def_geometries[i]);
        OpenMesh::IO::write_mesh(plate,filename);
        setGeometry(plate, ref_geometries[i]);
        filename = "reference/plate_" + std::to_string(i) + ".ply";
        OpenMesh::IO::write_mesh(plate, filename);
    }*/

  }catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }

    return 0;
}