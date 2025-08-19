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
    std::cerr << "SKEWED CROSS FOLD OPTIMIZATION" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    Eigen::setNbThreads(8);

    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/paperCrissCross_coarse.ply");
    MeshTopologySaver plateTopol( plate );
    std::cout<<"Hello2"<<std::endl;
    std::cout<<"num vertices: "<<plateTopol.getNumVertices()<<std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial; // plateGeomRef is the geometry without the arc crease in the mesh
    getGeometry(plate,plateGeomRef);
    getGeometry(plate,plateGeomDef);
    getGeometry(plate,plateGeomInitial);

    RealType t_0 = 0.3;
    RealType tolerance = 1e-4;

    std::cout<<"Fold vertices: "<<std::endl;
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomInitial, coords, i);
        if(std::abs(coords[0] - 0.5)<tolerance)
        {
            std::cout<<i<<std::endl;
        }
        if(std::abs(coords[1] - 0.5)<tolerance)
        {
            std::cout<<i<<std::endl;
        }
    }

    std::cout<<"Hello1"<<std::endl;

    // bdryMaskRef_1 fixes x,y,z coordinates in the reference geometry
    std::vector<int> bdryMaskRef_1;
    // bdryMaskRef_2 fixes y,z coordinates in the reference geometry
    std::vector<int> bdryMaskRef_2;
    // bdryMaskRef_3 fixes x,z coordinates in the reference geometry
    std::vector<int> bdryMaskRef_3;

    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomInitial, coords, i);
        // first, fix the outer vertices of the square
        if( (std::abs(coords[0] < tolerance) || std::abs(coords[0] - 1.0) < tolerance) && (std::abs(coords[1] < tolerance) || std::abs(coords[1] - 1.0) < tolerance) )
        {
            bdryMaskRef_1.push_back( i );
            continue;
        }
        /*
        // identify four small rectangles at the outer corners of the square at which we want to enforce clamped boundary conditions
        // rectangle 0 -> fix y positions
        if( coords[1] <= 0.16 && coords[0] <= 0.07)
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }
        // rectangle 1 -> fix x positions
        if( coords[1] >= (1-0.07) && coords[0] <= 0.16)
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }
        // rectangle 2 -> fix y positions
        if( coords[1] >= (1-0.16) && coords[0] >= (1-0.07))
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }

        if(coords[0] >= (1-0.16) && coords[1] <= 0.07)
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }*/

        // fix the middle vertex
        if(is_near(coords[0], 0.5, tolerance) && is_near(coords[1], 0.5, tolerance))
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }

        // Would like the possibility to let the vertices vary along possibly rotated fold lines
        // how to do this ? Difficult
        if(is_near(coords[0], 0.5, tolerance))
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }

        if(is_near(coords[1], 0.5, tolerance))
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }

        if(is_near(coords[0], 1.0) || is_near(coords[0], 0.0))
        {
            bdryMaskRef_3.push_back(i);
            continue;
        }

        if(is_near(coords[1], 1.0) || is_near(coords[1], 0.0))
        {
            bdryMaskRef_2.push_back(i);
            continue;
        }

    }

    std::cout<<"Hello2"<<std::endl;

    // store four parts of the boundary that would get squashed together by the rotation
    std::vector<int> scaling_piece_1, scaling_piece_2, scaling_piece_3, scaling_piece_4;
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords; 
        getXYZCoord<VectorType, VecType>( plateGeomInitial, coords, i);

        if((std::abs(coords[0]) < tolerance) && (coords[1] > 0.5 + tolerance) && (coords[1] < 1.0 - tolerance))
        {
            scaling_piece_1.push_back(i);
            bdryMaskRef_1.push_back(i);
            continue;
        }
        if((std::abs(coords[1] - 1.0) < tolerance) && (coords[0] > 0.5 + tolerance) && (coords[0] < 1.0 - tolerance))
        {
            scaling_piece_2.push_back(i);
            bdryMaskRef_1.push_back(i);
            continue;
        }
        if(std::abs(coords[0] - 1.0) < tolerance && (coords[1] < 0.5 - tolerance) && (coords[1] > tolerance))
        {
            scaling_piece_3.push_back(i);
            bdryMaskRef_1.push_back(i);
            continue;
        }
        if(std::abs(coords[1]) < tolerance && (coords[0] < 0.5 - tolerance) && (coords[0] > tolerance))
        {
            scaling_piece_4.push_back(i);
            bdryMaskRef_1.push_back(i);
            continue;
        }
    }
    
    std::cout<<"Hello3"<<std::endl;

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef_1 );
    std::vector<int> activeRef_2 = (std::vector<int>){0,1,1};
    std::vector<int> activeRef_3 = (std::vector<int>){1,0,1};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_2 , activeRef_2);
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_3 , activeRef_3);

    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntries(bdryMaskRef_1.begin(), bdryMaskRef_1.end());
    uniqueEntries.insert(bdryMaskRef_2.begin(), bdryMaskRef_2.end());
    uniqueEntries.insert(bdryMaskRef_3.begin(), bdryMaskRef_3.end());

    // Move the unique elements back to bdryMaskRef_1 in order
    bdryMaskRef_1.assign(uniqueEntries.begin(), uniqueEntries.end());
    std::sort(bdryMaskRef_1.begin(), bdryMaskRef_1.end());

    // determine boundary mask for optimization
    // and deform part of boundary
    std::vector<int> bdryMaskOpt;
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);

        // Neumann boundary conditions
        /*
        if( coords[1] <= 0.16 && coords[0] <= 0.07)
        {
            bdryMaskOpt.push_back(i);
            continue;
        }
        // rectangle 1 -> fix x positions
        if( coords[1] >= (1-0.07) && coords[0] <= 0.16)
        {
            bdryMaskOpt.push_back(i);
            continue;
        }
        // rectangle 2 -> fix y positions
        if( coords[1] >= (1-0.16) && coords[0] >= (1-0.07))
        {
            bdryMaskOpt.push_back(i);
            continue;
        }

        if(coords[0] >= (1-0.16) && coords[1] <= 0.07)
        {
            bdryMaskOpt.push_back(i);
            continue;
        }*/

        // Dirichlet boundary conditions
        if(std::abs(coords[0]) < tolerance && (coords[1] <= 0.16))
        {
            bdryMaskOpt.push_back(i);
            continue;
        }
        if(std::abs(coords[1] - 1.0) < tolerance && (coords[0] <= 0.16))
        {
            bdryMaskOpt.push_back(i);
            continue;
        }
        if(std::abs(coords[0] - 1.0) < tolerance && (coords[1] >= (1-0.16)))
        {
            bdryMaskOpt.push_back(i);
            continue;
        }
        if(std::abs(coords[1]) < tolerance && (coords[0] >= (1-0.16)))
        {
            bdryMaskOpt.push_back(i);
            continue;
        }

    }

    std::cout<<"Hello4"<<std::endl;

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );

    /*
    // Now, read the initial deformed plate
    OpenMesh::IO::read_mesh(plate, "SkewedCrossFoldInitPlate.ply");
    getGeometry(plate, plateGeomDef);
    std::cout<<"Plate number of vertices after: "<<plate.n_vertices();
    */

    auto foldDofsPtr = std::make_shared<FoldDofsSkewedCross<DefaultConfigurator>>( plateTopol, plateGeomInitial, plateGeomInitial, bdryMaskRef_1, scaling_piece_1, scaling_piece_2, scaling_piece_3, scaling_piece_4);
    //FoldDofsSkewedCross<DefaultConfigurator> foldDofs( plateTopol, plateGeomInitial, plateGeomInitial, bdryMaskRef_1 );
    
    std::vector<int> foldVertices;
    foldDofsPtr->getFoldVertices(foldVertices);
    std::cout<<"Number of fold vertices: "<<foldVertices.size()<<std::endl;
    
    auto DfoldDofsPtr = std::make_shared<FoldDofsSkewedCrossGradient<DefaultConfigurator>>( plateTopol, bdryMaskRef_1, plateGeomInitial, foldVertices, scaling_piece_1, scaling_piece_2, scaling_piece_3, scaling_piece_4);

    size_t nFoldDOFs = foldDofsPtr->getNumDofs();
    size_t nVertexDOFs = 3*plateTopol.getNumVertices();

    SQPLineSearchParams<DefaultConfigurator> pars;
    CostFunctional<DefaultConfigurator> costFunctional(foldVertices);
    CostFunctionalGradient<DefaultConfigurator> DcostFunctional(plateTopol,foldVertices);

    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;
    RealType factor_gravity = 0.0;

    VectorType factors(2);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;

    OpenMesh::IO::read_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/deformed_plate_after_optimization_dirichlet.ply");
    getGeometry(plate, plateGeomDef);

    VectorType edge_weights;
    foldDofsPtr->getEdgeWeights(edge_weights);

    std::vector<VectorType> def_geometries, ref_geometries, fold_DOFs;
    auto factory = std::make_shared<ElasticEnergyFactory<DefaultConfigurator>>(factors, plateTopol, edge_weights);
    BoundaryDOFS<DefaultConfigurator> boundaryDOFs(bdryMaskOpt, nVertexDOFs, nFoldDOFs);
    // Create the degrees of freedom object
    std::vector<RealType> deviations;
    VectorType t_initial = VectorType::Ones(1)*t_0;
    ProblemDOFs<DefaultConfigurator> problemDOFs(t_initial, plateGeomDef, VectorType::Zero(3*plateTopol.getNumVertices()),foldDofsPtr, DfoldDofsPtr);
    SQPLineSearchSolver<DefaultConfigurator> solver(pars, costFunctional, DcostFunctional, std::move(factory), boundaryDOFs, problemDOFs, 20);
    const std::string filename_def_geometries = "/lustre/scratch/data/s24pjoha_hpc-results/thesis_results/deformed_SkewedCrossFoldOptimization/";
    const std::string filename_ref_geometries = "/lustre/scratch/data/s24pjoha_hpc-results/thesis_results/reference_SkewedCrossFoldOptimization/";

    solver.solve(plateGeomRef, filename_def_geometries, filename_ref_geometries, plate);
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