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

#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <typeinfo>

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


int main(int argc, char *argv[])
{
try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE ARC FOLD SIMULATION WITH BENDING AND MEMBRANE ENERGIES" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;
    
// load flat plate and prepare the arc crease
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/rectangle_criss_cross.ply");
    MeshTopologySaver plateTopol( plate );
    VectorType plateGeomRef; // plateGeomRef is the geometry without the arc crease in the mesh
    VectorType plateGeomInitial;
    getGeometry(plate,plateGeomRef);
    getGeometry(plate,plateGeomInitial);

    // Also select all the vertices with x > 1.0 in the undeformed geometry
    // We will make a Gravitational Force act on them in order to achieve a flipping effect
    std::vector<int> upperVertices;
    std::vector<int> foldVertices;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);

        if(coords[0] > 1.0){
            upperVertices.push_back(i);
        }

        // we will use the fold vertices to modify the shape functions -> modify the stiffness matrix
        // in order to be able to regularize the mesh better using the Dirichlet minimization
        // In particular, we want a better mesh regularity at the fold
        if(coords[0] == 1.0){
            foldVertices.push_back(i);
        }
    }

    // first, set the bending edge weights for the edge at x = 0.5
    // Next, set the edge weights along the crease to zero
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

        // Set the edge weights for the edge at x = 0.5 to zero -> fold
        if(coords_i[0] == 1.0 && coords_j[0] == 1.0){
            edge_weights[edgeIdx] = 0;
            // Set colors for both vertices of the edge
        plate.set_color(plate.vertex_handle(node_i), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
        plate.set_color(plate.vertex_handle(node_j), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
        }
    }

    // induce the arc crease
    std::vector<int> bdryMaskDirichlet;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);

        if(coords[0] == 1.0){
            bdryMaskDirichlet.push_back(i);
            coords[0] += (0.25 - std::pow(coords[1] - 0.5, 2.0));
        }

        if(coords[0] == 0.0 || coords[0] == 2.0 || coords[1] == 0.0 || coords[1] == 1.0){
            bdryMaskDirichlet.push_back(i);
        }
        setXYZCoord<VectorType,VecType>(plateGeomRef, coords, i);
    }

    setGeometry(plate, plateGeomRef);
    OpenMesh::IO::write_mesh(plate, "plateWithArcCrease.ply");

    int numBdryNodes = bdryMaskDirichlet.size();
    std::cerr << "num of bdry nodes for the mesh regularization with Dirichlet Energy = " << numBdryNodes << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskDirichlet );

    // now, solve the Dirichlet Problem on the reference geometry

    // FIRST CONSIDER SOLUTION OF EULER-LAGRANGE-EQUATION
    std::cerr << "\n\na) OPTIMIZATION BY SOLVING EULER-LAGRANGE EQUATION" << std::endl;
    // assemble and mask stiffness matrix
    std::cerr << "Set up system matrix" << std::endl;
    typename DefaultConfigurator::SparseMatrixType StiffnessMatrix;
    computeStiffnessMatrix<DefaultConfigurator>( plateTopol, plateGeomInitial, StiffnessMatrix );

    // Now, scale the shape functions by multiplying every column and row that belongs to a shape 
    // function from the fold by a factor > 1.0
    // RealType alpha = -1.0;

    applyMaskToMajor<typename DefaultConfigurator::SparseMatrixType>( bdryMaskDirichlet, StiffnessMatrix );

    // set up right hand side and mask
    VectorType rhs( plateGeomRef );
    rhs.setZero();
    for( int i = 0; i < bdryMaskDirichlet.size(); i++ )
        rhs[bdryMaskDirichlet[i]] = plateGeomRef[bdryMaskDirichlet[i]];

    // set up linear system and solve
    std::cerr << "Set up linear system and solve" << std::endl;
    LinearSolver<DefaultConfigurator> directSolver;
    VectorType solution;
    directSolver.prepareSolver( StiffnessMatrix );
    directSolver.backSubstitute( rhs, solution );

    // get final energy value
    computeStiffnessMatrix<DefaultConfigurator>( plateTopol, plateGeomInitial, StiffnessMatrix );
    VectorType temp = StiffnessMatrix * solution;
    std::cerr << "Final Dirichlet energy = " << 0.5 * temp.dot( solution ) << std::endl;

    // saving
    setGeometry( plate, solution );
    getGeometry( plate, plateGeomRef );
    OpenMesh::IO::write_mesh(plate, "solutionDirichlet_EulerLangrange.ply");

    // NOW CONSIDER SOLUTION BY OPTIMIZING DIRICHLET ENERGY
    std::cerr << "\n\nb) OPTIMIZATION BY DIRECT MINIMIZATION" << std::endl;

    VectorType plateGeomDef;
    getGeometry( plate, plateGeomDef );

    OpenMesh::IO::Options opt;
    opt += OpenMesh::IO::Options::VertexColor;

    if(!OpenMesh::IO::write_mesh(plate, "initPlate.ply", opt)){
        std::cerr << "Cannot write mesh to file 'initPlate.ply'" << std::endl;
    }

    // determine boundary mask and deform part of boundary
    std::vector<int> bdryMaskOpt;

    // fix all coordinates of the boundary
    std::vector<int> bdryMaskDirichletDef1;

    // fix only (x,y) coordinates of the boundary
    std::vector<int> bdryMaskDirichletDef2;

    // fix only x coordinates of the boundary
    std::vector<int> bdryMaskDirichletDef3;

    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);

        if( coords[0] <= 1.0 && (coords[1] == 0.0 || coords[1] == 1.0) ){
            bdryMaskOpt.push_back( i );
            bdryMaskDirichletDef1.push_back( i );
            // deform part of boundary
            if(coords[1] == 0.0){
                coords[1] += 0.2;
                coords[2] += 0.15;
            }
            else{
                coords[1] -= 0.2;
                coords[2] += 0.15;
            }
        }
        if(coords[0] >= 1.2){
            bdryMaskDirichletDef2.push_back(i);
        }
        if(coords[0] == 0.0){
            bdryMaskDirichletDef3.push_back(i);
        }

        
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
    }

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "plateWithDeformedBoundary.ply");

    // CONSTRUCTION WORK STARTIN HERE CAUTION

    int numBdryNodesDef = bdryMaskDirichletDef1.size() + bdryMaskDirichletDef2.size() + bdryMaskDirichletDef3.size();
     std::cout<<"Individual sizes: "<<bdryMaskDirichletDef1.size()<<" "<<bdryMaskDirichletDef2.size()<<" "<<bdryMaskDirichletDef3.size()<<std::endl;
    std::cerr << "num of bdry nodes for the mesh regularization with Dirichlet Energy = " << numBdryNodesDef << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskDirichletDef1 );
    
    // want to only fix (x,y) coordinates of the boundary
    std::vector<int> active = (std::vector<int>){1,1,0};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskDirichletDef2 , active);

    std::vector<int> active2 = (std::vector<int>){1,0,0};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskDirichletDef3 , active2);

    // append to the Dirichlet bdry mask
    bdryMaskDirichletDef1.insert(bdryMaskDirichletDef1.end(), bdryMaskDirichletDef2.begin(), bdryMaskDirichletDef2.end());
    bdryMaskDirichletDef1.insert(bdryMaskDirichletDef1.end(), bdryMaskDirichletDef3.begin(), bdryMaskDirichletDef3.end());

    // now, solve the Dirichlet Problem on the reference geometry

    // FIRST CONSIDER SOLUTION OF EULER-LAGRANGE-EQUATION
    std::cerr << "\n\na) OPTIMIZATION BY SOLVING EULER-LAGRANGE EQUATION" << std::endl;

    std::cerr << "Set up system matrix" << std::endl;
    typename DefaultConfigurator::SparseMatrixType StiffnessMatrixDef;

    computeStiffnessMatrix<DefaultConfigurator>( plateTopol, plateGeomInitial, StiffnessMatrixDef );
    // assemble and mask stiffness matrix
    applyMaskToMajor<typename DefaultConfigurator::SparseMatrixType>( bdryMaskDirichletDef1, StiffnessMatrixDef );

    // set up right hand side and mask
    VectorType rhsDef( plateGeomDef );
    rhsDef.setZero();
    for( int i = 0; i < bdryMaskDirichletDef1.size(); i++ )
        rhsDef[bdryMaskDirichletDef1[i]] = plateGeomDef[bdryMaskDirichletDef1[i]];

    // set up linear system and solve
    std::cerr << "Set up linear system and solve" << std::endl;
    LinearSolver<DefaultConfigurator> directSolverDef;
    VectorType solutionDef;
    directSolverDef.prepareSolver( StiffnessMatrixDef );
    directSolverDef.backSubstitute( rhsDef, solutionDef );

    // get final energy value
    computeStiffnessMatrix<DefaultConfigurator>( plateTopol, plateGeomInitial, StiffnessMatrixDef );
    VectorType tempDef = StiffnessMatrixDef * solutionDef;
    std::cerr << "Final Dirichlet energy = " << 0.5 * tempDef.dot( solutionDef ) << std::endl;

    // saving
    setGeometry( plate, solutionDef );
    getGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "solutionDirichlet_EulerLangrange_2.ply");

    std::cerr << "num of bdry nodes for Optimization of SimpleBending + Membrane = " << bdryMaskOpt.size() << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );

    // save initialization
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "initPlate_rectangle2.ply");

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
    // 1000 looks good, but has self intersections. Maybe 300 is better
    RealType factor_gravity = 100.0;

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
    OpenMesh::IO::write_mesh(plate, "bendingFoldSol_withNewton2.ply");
    }
    catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }

    return 0;
}