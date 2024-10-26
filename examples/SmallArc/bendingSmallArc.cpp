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
 * \brief Simulation of a small arc fold
 * \author Johannssen
 *
 * We optimize this energy by direct optimization via gradient descent or BFGS.
 * 
 */

bool isOnUpperLeftDiagonal(RealType coords_x, RealType coords_y){
    bool isOnLine = false;
    if(0.0 <= coords_x && coords_x <= 0.5 && coords_x == 1.0 - coords_y){
        isOnLine = true;
    }
    return isOnLine;
}

bool isOnUpperRightDiagonal(RealType coords_x, RealType coords_y){
    bool isOnLine = false;
    if(0.5 <= coords_x && coords_x <= 1.0 && coords_x == coords_y){
        isOnLine = true;
    }
    return isOnLine;
}

bool isOnLowerLeftDiagonal(RealType coords_x, RealType coords_y){
    bool isOnLine = false;
    if(0.0 <= coords_x && coords_x <= 0.5 && coords_x == coords_y){
        isOnLine = true;
    }
    return isOnLine;
}

bool isOnLowerRightDiagonal(RealType coords_x, RealType coords_y){
    bool isOnLine = false;
    if(0.5 <= coords_x && coords_x <= 1.0 && coords_x == 1.0 - coords_y){
        isOnLine = true;
    }
    return isOnLine;
}

void rotate(RealType &coord_x, RealType &coord_y, RealType angle, RealType rotation_origin_x, RealType rotation_origin_y){
    //First, translate the rotation origin to the true origin
    coord_x -= rotation_origin_x;
    coord_y -= rotation_origin_y;
    //Next, rotate the coordinates around the origin
    RealType new_x = cos(angle)*coord_x - sin(angle)*coord_y;
    RealType new_y = sin(angle)*coord_x + cos(angle)*coord_y;
    //Finally, translate back to the rotation origin
    coord_x = new_x + rotation_origin_x;
    coord_y = new_y + rotation_origin_y;
}

// The arc is the result of translating a diagonal line from the top left corner so that it forms an arc in the upper left quadrant
RealType coordsToAngleNormUpperLeft(RealType coord_x, RealType coord_y){
    RealType rel_x = coord_x;
    RealType rel_y = 1.0 - coord_y;
    RealType norm = std::sqrt(rel_x*rel_x + rel_y*rel_y);
    RealType angle; 
    if(std::abs(norm - 1.0/sqrt(2.0)) <= 1e-10 ){
        angle = M_PI/2.0;
    }
    else{
        angle = std::atan(norm/((1.0/sqrt(2.0)) - norm));
    }
    return angle;
}

void translateDiagToArcUpperLeft(RealType &coords_x, RealType &coords_y){
    RealType coords_old_x = coords_x;
    RealType coords_old_y = coords_y;
    RealType angle = coordsToAngleNormUpperLeft(coords_x, coords_y);
    RealType arc_x = 0.5-cos(angle)*0.5;
    RealType arc_y = 1.0-sin(angle)*0.5;

    // Now, we don't want the arc to be tangential to the square at x=0.0, y=1.0.
    // Otherwise, a triangle in the corner is compressed
    // to remedy this, interpolate between the arc and the diagonal
    // So we don't have a pure arc.
    RealType new_x = arc_x;
    RealType new_y = arc_y;
    if (coords_old_x <= 0.25){
        RealType t = 0.1;
        new_x = (0.25-coords_old_x)*coords_old_x + (1.0-(0.25-coords_old_x))*arc_x;
        new_y = (0.25-coords_old_x)*coords_old_y + (1.0-(0.25-coords_old_x))*arc_y;
    }
    coords_x = new_x;
    coords_y = new_y;
}

void translateDiagToArcLowerLeft(RealType &coords_x, RealType &coords_y){
    RealType pi = M_PI;
    
    // 1. Rotate 90 degrees counterclockwise (instead of clockwise)
    rotate(coords_x, coords_y, -pi/2.0, 0.25, 0.25);

    // Translate the y-coordinate down (by -= 0.5, instead of +=)
    coords_y += 0.5;
    
    // Apply the same arc transformation logic for the upper left
    translateDiagToArcUpperLeft(coords_x, coords_y);

    // Translate back up by 0.5
    coords_y -= 0.5;

    // Rotate back by 90 degrees clockwise
    rotate(coords_x, coords_y, pi/2.0, 0.25, 0.25);
}

void translateDiagToArcLowerRight(RealType &coords_x, RealType &coords_y){
    RealType pi = M_PI;
    
    // 1. Rotate 90 degrees counterclockwise (instead of clockwise)
    rotate(coords_x, coords_y, -pi, 0.75, 0.25);

    // Translate the y-coordinate down (by -= 0.5, instead of +=)
    coords_y += 0.5;
    coords_x -= 0.5;
    
    // Apply the same arc transformation logic for the upper left
    translateDiagToArcUpperLeft(coords_x, coords_y);

    // Translate back up by 0.5
    coords_y -= 0.5;
    coords_x += 0.5;

    // Rotate back by 90 degrees clockwise
    rotate(coords_x, coords_y, pi, 0.75, 0.25);
}

void translateDiagToArcUpperRight(RealType &coords_x, RealType &coords_y){
    RealType pi = M_PI;
    
    // 1. Rotate 90 degrees counterclockwise (instead of clockwise)
    rotate(coords_x, coords_y, -3.0/2.0*pi, 0.75, 0.75);

    // Translate the y-coordinate down (by -= 0.5, instead of +=)
    coords_x -= 0.5;
    
    // Apply the same arc transformation logic for the upper left
    translateDiagToArcUpperLeft(coords_x, coords_y);

    // Translate back up by 0.5
    coords_x += 0.5;

    // Rotate back by 90 degrees clockwise
    rotate(coords_x, coords_y, 3.0/2.0*pi, 0.75, 0.75);
}

int main(int argc, char *argv[])
{
try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMULATION WITH ARC IN UPPER LEFT QUADRANT USING BENDING AND MEMBRAN ENERGIES" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;
    
// load flat plate and prepare the arc crease
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/paperCrissCross_refined.ply");
    MeshTopologySaver plateTopol( plate );
    VectorType plateGeomRef;
    VectorType plateGeomInitial;
    getGeometry(plate,plateGeomRef);
    getGeometry(plate,plateGeomInitial);

    // first, set the bending edge weights for the edge at x = 0.5
    // Next, set the edge weights along the crease to zero
    VectorType edge_weights = VectorType::Ones(plateTopol.getNumEdges());

    plate.request_edge_colors();
    plate.request_vertex_colors();

    for(int edgeIdx = 0; edgeIdx < plateTopol.getNumEdges(); edgeIdx++){
        int node_i = plateTopol.getAdjacentNodeOfEdge(edgeIdx,0);
        int node_j = plateTopol.getAdjacentNodeOfEdge(edgeIdx,1);

        VecType coords_i, coords_j;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords_i, node_i);
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords_j, node_j);

        if(isOnUpperLeftDiagonal(coords_i[0], coords_i[1]) && isOnUpperLeftDiagonal(coords_j[0], coords_j[1])){
            edge_weights[edgeIdx] = 0;
            // Set colors for both vertices of the edge
            plate.set_color(plate.vertex_handle(node_i), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
            plate.set_color(plate.vertex_handle(node_j), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
        }

        if(isOnLowerLeftDiagonal(coords_i[0], coords_i[1]) && isOnLowerLeftDiagonal(coords_j[0], coords_j[1])){
            edge_weights[edgeIdx] = 0;
            // Set colors for both vertices of the edge
            plate.set_color(plate.vertex_handle(node_i), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
            plate.set_color(plate.vertex_handle(node_j), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
        }

        if(isOnLowerRightDiagonal(coords_i[0], coords_i[1]) && isOnLowerRightDiagonal(coords_j[0], coords_j[1])){
            edge_weights[edgeIdx] = 0;
            // Set colors for both vertices of the edge
            plate.set_color(plate.vertex_handle(node_i), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
            plate.set_color(plate.vertex_handle(node_j), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
        }

        if(isOnUpperRightDiagonal(coords_i[0], coords_i[1]) && isOnUpperRightDiagonal(coords_j[0], coords_j[1])){
            edge_weights[edgeIdx] = 0;
            // Set colors for both vertices of the edge
            plate.set_color(plate.vertex_handle(node_i), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
            plate.set_color(plate.vertex_handle(node_j), OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));  // Red
        }
    }

    OpenMesh::IO::Options opt;
    opt += OpenMesh::IO::Options::VertexColor;
    // induce the arc crease
    std::vector<int> bdryMaskDirichlet;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);

        if(coords[0] == 0.0 || coords[0] == 1.0 || coords[1] == 0.0 || coords[1] == 1.0){
            bdryMaskDirichlet.push_back(i);
        }
        else if(isOnUpperLeftDiagonal(coords[0], coords[1])){
            translateDiagToArcUpperLeft(coords[0], coords[1]);
            bdryMaskDirichlet.push_back(i);
        }
        else if(isOnLowerLeftDiagonal(coords[0], coords[1])){
            translateDiagToArcLowerLeft(coords[0], coords[1]);
            bdryMaskDirichlet.push_back(i);
        }
        else if(isOnLowerRightDiagonal(coords[0], coords[1])){
            translateDiagToArcLowerRight(coords[0], coords[1]);
            bdryMaskDirichlet.push_back(i);
        }
        else if(isOnUpperRightDiagonal(coords[0], coords[1])){
            translateDiagToArcUpperRight(coords[0], coords[1]);
            bdryMaskDirichlet.push_back(i);
        }

        setXYZCoord<VectorType,VecType>(plateGeomRef, coords, i);
    }

    //dirichlet
    int numBdryNodesDirichlet = bdryMaskDirichlet.size();
    std::cerr << "num of bdry nodes for the mesh regularization with Dirichlet Energy = " << numBdryNodesDirichlet << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskDirichlet );

    // FIRST CONSIDER SOLUTION OF EULER-LAGRANGE-EQUATION
    std::cerr << "\n\na) OPTIMIZATION BY SOLVING EULER-LAGRANGE EQUATION" << std::endl;
    // assemble and mask stiffness matrix
    std::cerr << "Set up system matrix" << std::endl;
    typename DefaultConfigurator::SparseMatrixType StiffnessMatrix;
    computeStiffnessMatrix<DefaultConfigurator>( plateTopol, plateGeomInitial, StiffnessMatrix );
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
    OpenMesh::IO::write_mesh(plate, "solutionDirichlet_EulerLangrange.ply");

    getGeometry( plate, plateGeomRef );
    VectorType plateGeomDef;
    getGeometry( plate, plateGeomDef );

    std::vector<int> bdryMaskOpt;
    std::vector<int> bdryMaskDirichletDef;
    for(int i = 0; i < plateTopol.getNumVertices(); i++){
        VecType coords;
        getXYZCoord<VectorType, VecType>(plateGeomDef, coords, i);

        const double tolerance = 1e-6;

        if (std::abs(coords[1] - 0.0) < tolerance && std::abs(coords[0] - 0.5) < tolerance) {
            coords[1] += 0.1;
            coords[2] -= 0.1;
            bdryMaskOpt.push_back(i);
            bdryMaskDirichletDef.push_back(i);
        }

        if (std::abs(coords[0] - 0.0) < tolerance && std::abs(coords[1] - 0.5) < tolerance) {
            coords[0] += 0.1;
            coords[2] -= 0.1;
            bdryMaskOpt.push_back(i);
            bdryMaskDirichletDef.push_back(i);
        }

        if (std::abs(coords[1] - 1.0) < tolerance && std::abs(coords[0] - 0.5) < tolerance) {
            coords[1] -= 0.1;
            coords[2] -= 0.1;
            bdryMaskOpt.push_back(i);
            bdryMaskDirichletDef.push_back(i);
        }

        if (std::abs(coords[0] - 1.0) < tolerance && std::abs(coords[1] - 0.5) < tolerance) {
            coords[0] -= 0.1;
            coords[2] -= 0.1;
            bdryMaskOpt.push_back(i);
            bdryMaskDirichletDef.push_back(i);
        }

        if(std::abs(coords[0] - 1.0) <= tolerance &&  (coords[1] <= 0.3 || coords[1] >= 0.7)){
            bdryMaskDirichletDef.push_back(i);
        }
        else if(std::abs(coords[0]) <= tolerance && (coords[1] <= 0.3 || coords[1] >= 0.7)){
            bdryMaskDirichletDef.push_back(i);
        }
        else if(std::abs(coords[1] - 1.0) <= tolerance && (coords[0] <= 0.3 || coords[0] >= 0.7)){
            bdryMaskDirichletDef.push_back(i);
        }
        else if(std::abs(coords[1]) <= tolerance && (coords[0] <= 0.3 || coords[0] >= 0.7)){
            bdryMaskDirichletDef.push_back(i);
        }

        // Set the updated coordinates
        setXYZCoord<VectorType, VecType>(plateGeomDef, coords, i);
    }

    setGeometry( plate, plateGeomDef );

    OpenMesh::IO::write_mesh(plate, "initPlateSmallArcCreaseOpt.ply");

    int numBdryNodesDef = bdryMaskDirichletDef.size();
    std::cerr << "num of bdry nodes for the mesh regularization with Dirichlet Energy = " << numBdryNodesDef << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskDirichletDef );
    
    // now, solve the Dirichlet Problem on the reference geometry

    // FIRST CONSIDER SOLUTION OF EULER-LAGRANGE-EQUATION
    std::cerr << "\n\na) OPTIMIZATION BY SOLVING EULER-LAGRANGE EQUATION" << std::endl;

    std::cerr << "Set up system matrix" << std::endl;
    typename DefaultConfigurator::SparseMatrixType StiffnessMatrixDef;

    computeStiffnessMatrix<DefaultConfigurator>( plateTopol, plateGeomInitial, StiffnessMatrixDef );
    // assemble and mask stiffness matrix
    applyMaskToMajor<typename DefaultConfigurator::SparseMatrixType>( bdryMaskDirichletDef, StiffnessMatrixDef );

    // set up right hand side and mask
    VectorType rhsDef( plateGeomDef );
    rhsDef.setZero();
    for( int i = 0; i < bdryMaskDirichletDef.size(); i++ )
        rhsDef[bdryMaskDirichletDef[i]] = plateGeomDef[bdryMaskDirichletDef[i]];

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

    int numBdryNodesOpt = bdryMaskOpt.size();
    std::cerr << "num of bdry nodes = " << numBdryNodesOpt << std::endl;
    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );

    std::cerr << "\n\nb) OPTIMIZATION BY DIRECT MINIMIZATION" << std::endl;

    std::cout<< "THIS IS THE NUMBER OF EDGES: "<< plateTopol.getNumEdges()<<std::endl;

    SimpleBendingEnergy<DefaultConfigurator> E_bend( plateTopol, plateGeomRef, true , edge_weights);
    SimpleBendingGradientDef<DefaultConfigurator> DE_bend( plateTopol, plateGeomRef , edge_weights);
    SimpleBendingHessianDef<DefaultConfigurator> D2E_bend( plateTopol, plateGeomRef , edge_weights);

    typename DefaultConfigurator::RealType energy;

    NonlinearMembraneEnergy<DefaultConfigurator> E_mem( plateTopol, plateGeomRef, true );
    NonlinearMembraneGradientDef<DefaultConfigurator> DE_mem( plateTopol, plateGeomRef );
    NonlinearMembraneHessianDef<DefaultConfigurator> D2E_mem( plateTopol, plateGeomRef );

    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;

    VectorType factors(2);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;
    //factors[2] = factor_gravity;

    AdditionOp<DefaultConfigurator> E_tot( factors, E_mem, E_bend);//, E_grav);
    AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_mem, DE_bend);//, DE_grav);
    AdditionHessian<DefaultConfigurator> D2E_tot(factors, D2E_mem, D2E_bend);

    // set outer optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setNewtonIterations(1000);
    optPars.setQuietMode( SHOW_ALL );
    VectorType initialization = plateGeomDef;

    
    std::cerr<< "Starting Newton Linesearch..."<<std::endl;
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