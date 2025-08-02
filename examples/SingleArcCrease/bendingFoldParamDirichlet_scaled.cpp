#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACK

#define OPENBLAS_VERBOSE 1

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
 * \brief Simulation of a simple fold 
 * \author Johannssen
 *
 * We optimize this energy by direct optimization via gradient descent or BFGS.
 * 
 */

RealType angle = M_PI / 3.0;
RealType mesh_delta = 0.015625;

bool is_near(RealType coord, RealType value, RealType tol = 1e-6)
{
    return std::abs(coord - value) < tol;
}

void rotateBoundary(RealType &x, RealType &y, RealType &z, int i, std::vector<int>& bdryMask)
{
    if(y <= 0.25 + 1e-4) {
        if(x < mesh_delta + 1e-4)
        {
            // First translate to x = 0
            x -= mesh_delta/2.0;
            RealType x_new = std::cos(angle) * x;
            RealType z_new = -std::sin(angle) * x;
            x = x_new + mesh_delta/2.0;
            z = z_new;
            bdryMask.push_back(i);
            x += 0.1;
        }
        else if(x > (1.0 - mesh_delta) - 1e-4)
        {
            // First, translate x to 0
            x -= (1.0 - mesh_delta/2.0);
            RealType x_new = std::cos(-angle) * x;
            RealType z_new = -std::sin(-angle) * x;
            x = x_new + (1.0 - mesh_delta/2.0);
            z = z_new;
            bdryMask.push_back(i);
            x -= 0.1;
        }
        else {
            x = 0.8033322*(x-0.5) + 0.5;
            z = 2.2531*(x-0.5)*(x-0.5)-0.3329;
        }
    }
}

int main(int argc, char *argv[])
{

try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE ARC FOLD SIMULATION WITH BENDING AND MEMBRANE ENERGIES" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    Eigen::setNbThreads(8);

// load flat plate and prepare the arc crease
    TriMesh plate;
    //OpenMesh::IO::read_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/paperCrissCross.ply");
    OpenMesh::IO::read_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/refine_1.ply");
    MeshTopologySaver plateTopol( plate );
    std::cout<<"num vertices: "<<plateTopol.getNumVertices()<<std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial; // plateGeomRef is the geometry without the arc crease in the mesh
    getGeometry(plate,plateGeomRef);
    getGeometry(plate,plateGeomDef);
    getGeometry(plate,plateGeomInitial);

    // Fix all coordinates of the boundary
    std::vector<int> bdryMaskRef;

    RealType t_0 = 0.25;//0.8084239959716796; // Initial value of the parameter t
    // Select the vertices with y = 1.0 (topVertices)
    // We will make a Gravitational Force act on them in order to achieve a flipping effect
    std::vector<int> topVertices;
    std::vector<int> upperVertices;
    std::vector<int> foldVertices;
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
            foldVertices.push_back(i);
            coords[1] += t_0*(0.25 - std::pow(coords[0] - 0.5, 2.0));
            setXYZCoord<VectorType, VecType>(plateGeomRef, coords, i);
            continue;
        }
    }

    std::cout<<"FoldVertices: "<<std::endl;
    for(int i = 0; i < foldVertices.size(); i++)
    {
        std::cout<<foldVertices[i]<<",  ";
    }
    std::cout<<"End fold vertices"<<std::endl;

    // output preliminary meshes
    setGeometry(plate, plateGeomRef);
    OpenMesh::IO::write_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/reference_mesh_before_smoothing.ply");

    // ----------------- Smooth reference geometry -----------------------------

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef );

    DirichletSmoother<DefaultConfigurator> smoother(plateGeomInitial, bdryMaskRef, plateTopol);
    smoother.apply(plateGeomRef,plateGeomRef);

    setGeometry(plate,plateGeomRef);
    OpenMesh::IO::write_mesh(plate,"/home/s24pjoha_hpc/goast_old_old/goast/build/examples/reference_mesh.ply");

    getGeometry(plate,plateGeomDef);

    // Fix all coordinates
    std::vector<int> bdryMaskOpt;
    // fix all coordinates of the boundary
    std::vector<int> bdryMaskDirichletDef1;

    // fix only (x,y) coordinates of the boundary
    std::vector<int> bdryMaskDirichletDef2;

    // fix only y coordinates of the boundary
    std::vector<int> bdryMaskDirichletDef3;

    //OpenMesh::IO::write_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/bendingFoldSol_withNewton_scaled_2.ply");
    //getGeometry(plate, plateGeomDef);

    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>(plateGeomInitial,coords,i);
        RealType x,y,z;
        x = coords[0];
        y = coords[1];
        z = coords[2];

        if(coords[1] == 0)
        { 
            //bdryMaskOpt.push_back(i);
        }

        rotateBoundary(x, y, z, i, bdryMaskOpt);

        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
    }

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/plateWithDeformedBoundary.ply");

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskDirichletDef1 );
    
    // want to only fix (x,y) coordinates of the boundary
    std::vector<int> active = (std::vector<int>){1,1,0};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskDirichletDef2 , active);

    std::vector<int> active2 = (std::vector<int>){0,1,0};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskDirichletDef3 , active2);

    // append to the Dirichlet bdry mask
    bdryMaskDirichletDef1.insert(bdryMaskDirichletDef1.end(), bdryMaskDirichletDef2.begin(), bdryMaskDirichletDef2.end());
    bdryMaskDirichletDef1.insert(bdryMaskDirichletDef1.end(), bdryMaskDirichletDef3.begin(), bdryMaskDirichletDef3.end());

    //DirichletSmoother<DefaultConfigurator> dirichletSmoother( plateGeomDef, bdryMaskDirichletDef1, plateTopol );
    //dirichletSmoother.apply( plateGeomDef, plateGeomDef );
    
    //interpolate instead
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords_initial, coords_def;
        getXYZCoord<VectorType, VecType>( plateGeomInitial, coords_initial, i);
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords_def, i);
        RealType x,y,z;

        if( coords_initial[1] > (0.25 + 1e-4)){
            RealType x_new = 0.78438*(coords_initial[0]-0.5) + 0.5;
            RealType parabola_z = 2.2531*(x_new - 0.5)*(x_new - 0.5) - 0.3329;
            //RealType dist_from_parabola = coords_def[1] - t_0*(0.25 - (coords_initial[0] - 0.5) * (coords_initial[0] - 0.5)); 
            //RealType max_dist_from_parabola = 1.0 - t_0*(0.25 - (coords_initial[0] - 0.5) * (coords_initial[0] - 0.5));
            // interpolate from that distance
            coords_def[2] =  parabola_z;//(1.0 - dist_from_parabola / max_dist_from_parabola) * parabola_z + (dist_from_parabola / max_dist_from_parabola) * coords_def[2];
        } 
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords_def, i);
    }
    
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/smoothed_initial_plate.ply");

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );
   
    auto foldDofsPtr = std::make_shared<FoldDofsArcLine<DefaultConfigurator>>(plateTopol,plateGeomInitial, plateGeomInitial, bdryMaskRef);

    VectorType edge_weights = VectorType::Zero(plateTopol.getNumEdges());
    foldDofsPtr->getEdgeWeights(edge_weights);

    size_t nFoldDOFs = foldDofsPtr->getNumDofs();
    size_t nVertexDOFs = 3*plateTopol.getNumVertices();

    // save initialization
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/initPlate_rectangle2.ply");

    SimpleBendingEnergy<DefaultConfigurator> E_bend( plateTopol, plateGeomRef, true , edge_weights);
    SimpleBendingGradientDef<DefaultConfigurator> DE_bend( plateTopol, plateGeomRef , edge_weights);
    SimpleBendingHessianDef<DefaultConfigurator> D2E_bend( plateTopol, plateGeomRef , edge_weights);

    typename DefaultConfigurator::RealType energy;

    NonlinearMembraneEnergy<DefaultConfigurator> E_mem( plateTopol, plateGeomRef, true );
    NonlinearMembraneGradientDef<DefaultConfigurator> DE_mem( plateTopol, plateGeomRef );
    NonlinearMembraneHessianDef<DefaultConfigurator> D2E_mem( plateTopol, plateGeomRef );

    // Add Gravitational Energy to the upper part of the plate in order to make it flip
    // We need the 
    VectorType mass_distribution = VectorType::Zero(plateTopol.getNumVertices() );
    
    for (int idx : upperVertices) {
        mass_distribution[idx] = 1;
    }

    //This time, want a force pushing the upper part of the plate upwards
    VectorType gravity_dir;
    gravity_dir.resize(3);
    gravity_dir << 0.0, 0.0, -1.0;

    GravitationalEnergy<DefaultConfigurator> E_grav( plateTopol, plateGeomRef, true, mass_distribution, gravity_dir );
    GravitationalEnergyGradientDef<DefaultConfigurator> DE_grav( plateTopol, plateGeomRef, mass_distribution, gravity_dir );

    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;
    RealType factor_gravity = 500.0;

    VectorType factors(3);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;
    factors[2] = factor_gravity;

    AdditionOp<DefaultConfigurator> E_tot( factors, E_mem, E_bend, E_grav);
    AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_mem, DE_bend, DE_grav);
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

    // saving
    setGeometry( plate, plateGeomDef );
    // slightly inaccurate saving of the geometry
    OpenMesh::IO::write_mesh(plate, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/bendingFoldSol_withNewton_scaled_2.ply");
    printVectorToFile( plateGeomDef, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/bendingFoldSol_withNewton_scaled_2.txt" , 15);
    printVectorToFile( plateGeomRef, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/final_reference.txt" , 15);

    // Write foldVertices to file
    std::ofstream vertexColorsFile("/home/s24pjoha_hpc/goast_old_old/goast/build/examples/VertexColors.txt");
    if (vertexColorsFile.is_open()) {
        for (const int idx : foldVertices) {
            vertexColorsFile << foldVertices[idx] << std::endl;
        }
        vertexColorsFile.close();
    } else {
        std::cerr << "Failed to open VertexColors.txt for writing." << std::endl;
    }

    }
    catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }

    return 0;
}
