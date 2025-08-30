#include <cmath>
#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <iostream>
#include <unordered_set>

#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <goast/Smoothers.h>
#include <goast/SQP.h>

typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;
typedef DefaultConfigurator::SparseMatrixType MatrixType;
typedef DefaultConfigurator::FullMatrixType FullMatrixType;

bool is_near(RealType coord, RealType value, RealType tol = 1e-6)
{
    return std::abs(coord - value) < tol;
}

int main(int argc, char *argv[])
{

try{
    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE FOLD OPTIMIZATION OF A ROTATING CROSS FOLD" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    Eigen::setNbThreads(8);  // Adjust based on CPU cores

    // load flat plate [0,1]^2
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/paperCrissCross_coarse.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial, plateGeomTest;
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );
    getGeometry( plate, plateGeomInitial );
    getGeometry( plate, plateGeomTest );

    // bdryMaskRef_1 fixes x,y,z coordinates in the reference geometry
    std::vector<int> bdryMaskRef_1;
    // bdryMaskRef_2 fixes y,z coordinates in the reference geometry
    std::vector<int> bdryMaskRef_2;
    // bdryMaskRef_3 fixes x,z coordinates in the reference geometry
    std::vector<int> bdryMaskRef_3;

    RealType t_init = 0.0;
    RealType tolerance = 1e-4;

    VectorType edge_weights = VectorType::Ones(plateTopol.getNumEdges());

    // Just for checking, we will cover the fold edges red
    plate.request_edge_colors();
    plate.request_vertex_colors();

    for(int edgeIdx = 0; edgeIdx < plateTopol.getNumEdges(); edgeIdx++){
        int node_i = plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
        int node_j = plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

        VecType coords_i, coords_j;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords_i, node_i);
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords_j, node_j);
        //std::cout<<std::setprecision(12)<<"coords0: "<<coords_i[0]<<" "<<coords_j[0]<<" "<<std::endl;
        //std::cout<<std::setprecision(12)<<"coords1: "<<coords_i[1]<<" "<<coords_j[1]<<" "<<std::endl;
        // Set the edge weights for the edge at x = 0.5 to zero -> fold
        if((is_near(coords_i[1],0.5, 1e-4) && is_near(coords_j[1], 0.5, 1e-4)) ||(is_near(coords_i[0],0.5, 1e-4) && is_near(coords_j[0],0.5, 1e-4))){
            edge_weights[edgeIdx] = 0;
            // Set colors for both vertices of the edge
            plate.set_color(plate.vertex_handle(node_i), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
            plate.set_color(plate.vertex_handle(node_j), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
        }
    }

    // bdryMaskRef_1 fixes x,y,z coordinates in the reference geometry
    std::vector<int> bdryMaskDirichletDef_1;
    // bdryMaskRef_2 fixes y coordinates in the reference geometry
    std::vector<int> bdryMaskDirichletDef_2;
    // bdryMaskRef_3 fixes x coordinates in the reference geometry
    std::vector<int> bdryMaskDirichletDef_3;
    // bdryMaskRef_4 fixes x,y coordinates in the reference geometry
    std::vector<int> bdryMaskDirichletDef_4;

    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomInitial, coords, i);
        // first, fix the outer vertices of the square
        if( (std::abs(coords[0] < tolerance) || std::abs(coords[0] - 1.0) < tolerance) && (std::abs(coords[1] < tolerance) || std::abs(coords[1] - 1.0) < tolerance) )
        {
            bdryMaskRef_1.push_back( i );
            bdryMaskDirichletDef_1.push_back(i);
            continue;
        }
        // identify four small rectangles at the outer corners of the square at which we want to enforce clamped boundary conditions
        // rectangle 0 -> fix y positions
        /*
        if( coords[1] <= 0.16 && coords[0] <= 0.07)
        {
            bdryMaskRef_1.push_back(i);
            bdryMaskDirichletDef_1.push_back(i);
            continue;
        }
        // rectangle 1 -> fix x positions
        if( coords[1] >= (1-0.07) && coords[0] <= 0.16)
        {
            bdryMaskRef_1.push_back(i);
            bdryMaskDirichletDef_1.push_back(i);
            continue;
        }
        // rectangle 2 -> fix y positions
        if( coords[1] >= (1-0.16) && coords[0] >= (1-0.07))
        {
            bdryMaskRef_1.push_back(i);
            bdryMaskDirichletDef_1.push_back(i);
            continue;
        }

        if(coords[0] >= (1-0.16) && coords[1] <= 0.07)
        {
            bdryMaskRef_1.push_back(i);
            bdryMaskDirichletDef_1.push_back(i);
            continue;
        }*/

        // fix the middle vertex
        if(is_near(coords[0], 0.5, tolerance) && is_near(coords[1], 0.5, tolerance))
        {
            bdryMaskRef_1.push_back(i);
            bdryMaskDirichletDef_1.push_back(i);
            continue;
        }

        // Would like the possibility to let the vertices vary along possibly rotated fold lines
        // how to do this ? Difficult
        if(is_near(coords[0], 0.5, tolerance))
        {
            bdryMaskRef_1.push_back(i);
            bdryMaskDirichletDef_4.push_back(i);
            continue;
        }

        if(is_near(coords[1], 0.5, tolerance))
        {
            bdryMaskRef_1.push_back(i);
            bdryMaskDirichletDef_4.push_back(i);
            continue;
        }

        if(is_near(coords[0], 1.0) || is_near(coords[0], 0.0))
        {
            bdryMaskRef_3.push_back(i);
            bdryMaskDirichletDef_3.push_back(i);
            continue;
        }

        if(is_near(coords[1], 1.0) || is_near(coords[1], 0.0))
        {
            bdryMaskRef_2.push_back(i);
            bdryMaskDirichletDef_2.push_back(i);
            continue;
        }

    }

    // store four parts of the boundary that would get squashed together by the rotation
    std::vector<int> scaling_piece_1, scaling_piece_2, scaling_piece_3, scaling_piece_4;
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords; 
        getXYZCoord<VectorType, VecType>( plateGeomInitial, coords, i);

        if((std::abs(coords[0]) < tolerance) && (coords[1] > 0.5 + tolerance) && (coords[1] < 1.0 - tolerance))
        {
            scaling_piece_1.push_back(i);
            //bdryMaskRef_1.push_back(i);
            continue;
        }
        if((std::abs(coords[1] - 1.0) < tolerance) && (coords[0] > 0.5 + tolerance) && (coords[0] < 1.0 - tolerance))
        {
            scaling_piece_2.push_back(i);
            //bdryMaskRef_1.push_back(i);
            continue;
        }
        if(std::abs(coords[0] - 1.0) < tolerance && (coords[1] < 0.5 - tolerance) && (coords[1] > tolerance))
        {
            scaling_piece_3.push_back(i);
            //bdryMaskRef_1.push_back(i);
            continue;
        }
        if(std::abs(coords[1]) < tolerance && (coords[0] < 0.5 - tolerance) && (coords[0] > tolerance))
        {
            scaling_piece_4.push_back(i);
            //bdryMaskRef_1.push_back(i);
            continue;
        }
    }

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

    // Now, rotate the fold
    FoldDofsCrossInterpolation<DefaultConfigurator> foldDofs( plateTopol, plateGeomInitial, plateGeomInitial, bdryMaskRef_1);
    std::vector<int> foldVertices;
    foldDofs.getFoldVertices(foldVertices);
    VectorType t_0(2);
    t_0[0] = 10000.0;
    t_0[1] = 10000.0;
    foldDofs.apply(t_0, plateGeomTest);
    setGeometry( plate, plateGeomTest );
    OpenMesh::IO::write_mesh(plate, "interpolated.ply");
    //FoldDofsSkewedCrossGradient<DefaultConfigurator> DfoldDofs( plateTopol, bdryMaskRef_1, plateGeomInitial, foldVertices, scaling_piece_1, scaling_piece_2, scaling_piece_3, scaling_piece_4);

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
            coords[0] += 0.1;
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
            continue;
        }
        // rectangle 1 -> fix x positions
        if( coords[1] >= (1-0.07) && coords[0] <= 0.16)
        {
            bdryMaskOpt.push_back(i);
            coords[1] -= 0.1;
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
            continue;
        }
        // rectangle 2 -> fix y positions
        if( coords[1] >= (1-0.16) && coords[0] >= (1-0.07))
        {
            bdryMaskOpt.push_back(i);
            coords[0] -= 0.1;
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
            continue;
        }

        if(coords[0] >= (1-0.16) && coords[1] <= 0.07)
        {
            bdryMaskOpt.push_back(i);
            coords[1] += 0.1;
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
            continue;
        }*/

        // Dirichlet boundary conditions
        // Dirichlet boundary conditions
        if(std::abs(coords[0]) < tolerance && (coords[1] <= 0.16))
        {
            bdryMaskOpt.push_back(i);
            coords[0] += 0.1;
            coords[1] += 0.1;
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
            continue;
        }
        if(std::abs(coords[1] - 1.0) < tolerance && (coords[0] <= 0.16))
        {
            bdryMaskOpt.push_back(i);
            coords[0] += 0.1;
            coords[1] -= 0.1;
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
            continue;
        }
        if(std::abs(coords[0] - 1.0) < tolerance && (coords[1] >= (1-0.16)))
        {
            bdryMaskOpt.push_back(i);
            coords[1] -= 0.1;
            coords[0] -= 0.1;
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
            continue;
        }
        if(std::abs(coords[1]) < tolerance && (coords[0] >= (1-0.16)))
        {
            bdryMaskOpt.push_back(i);
            coords[0] -= 0.1;
            coords[1] += 0.1;
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
            continue;
        }
        
    }

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskDirichletDef_1 );
    std::vector<int> activeRef_DirichletDef_2 = (std::vector<int>){0,1,0};
    std::vector<int> activeRef_DirichletDef_3 = (std::vector<int>){1,0,0};
    std::vector<int> activeRef_DirichletDef_4 = (std::vector<int>){1,1,0};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskDirichletDef_2 , activeRef_DirichletDef_2);
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskDirichletDef_3 , activeRef_DirichletDef_3);
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskDirichletDef_4 , activeRef_DirichletDef_4);

    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntries_DirichletDef(bdryMaskDirichletDef_1.begin(), bdryMaskDirichletDef_1.end());
    uniqueEntries_DirichletDef.insert(bdryMaskDirichletDef_2.begin(), bdryMaskDirichletDef_2.end());
    uniqueEntries_DirichletDef.insert(bdryMaskDirichletDef_3.begin(), bdryMaskDirichletDef_3.end());
    uniqueEntries_DirichletDef.insert(bdryMaskDirichletDef_4.begin(), bdryMaskDirichletDef_4.end());

    // Move the unique elements back to bdryMaskRef_1 in order
    bdryMaskDirichletDef_1.assign(uniqueEntries_DirichletDef.begin(), uniqueEntries_DirichletDef.end());
    std::sort(bdryMaskDirichletDef_1.begin(), bdryMaskDirichletDef_1.end());

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh( plate, "deformed_plate_before_smoothing.ply" );
    //DirichletSmoother<DefaultConfigurator> dirichletSmoother( plateGeomDef, bdryMaskDirichletDef_1, plateTopol );
    //dirichletSmoother.apply( plateGeomDef, plateGeomDef );

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh( plate, "deformed_plate_after_smoothing.ply" );

    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;

    VectorType factors(2);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;

    SimpleBendingEnergy<DefaultConfigurator> E_bend( plateTopol, plateGeomRef, true , edge_weights);
    SimpleBendingGradientDef<DefaultConfigurator> DE_bend( plateTopol, plateGeomRef , edge_weights);
    SimpleBendingHessianDef<DefaultConfigurator> D2E_bend( plateTopol, plateGeomRef , edge_weights);

    typename DefaultConfigurator::RealType energy;

    NonlinearMembraneEnergy<DefaultConfigurator> E_mem( plateTopol, plateGeomRef, true );
    NonlinearMembraneGradientDef<DefaultConfigurator> DE_mem( plateTopol, plateGeomRef );
    NonlinearMembraneHessianDef<DefaultConfigurator> D2E_mem( plateTopol, plateGeomRef );

    AdditionOp<DefaultConfigurator> E_tot( factors, E_mem, E_bend);
    AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_mem, DE_bend);
    //Gravitational Energy has Hessian = 0
    AdditionHessian<DefaultConfigurator> D2E_tot(factors, D2E_mem, D2E_bend);

    // set outer optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setNewtonIterations( 1000 );
    optPars.setQuietMode( SHOW_ALL );
    VectorType initialization = plateGeomDef;

    std::cerr<< "Startting Newton Linesearch..."<<std::endl;
    LineSearchNewton<DefaultConfigurator> NLS( E_tot, DE_tot, D2E_tot, optPars);
    NLS.setBoundaryMask( bdryMaskOpt );
    NLS.solve( initialization, plateGeomDef );

    CostFunctional<DefaultConfigurator> costFunctional(foldVertices);
    RealType costFuntional_val = 0.0;
    costFunctional.apply( plateGeomDef, costFuntional_val);
    std::cout<<"Cost functional after optimization: "<<std::setprecision(15)<<costFuntional_val<<std::endl;

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh( plate, "deformed_plate_after_optimization_dirichlet.ply" );
} 
catch ( BasicException &el ){
      std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
}

return 0;
}
