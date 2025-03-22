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

int main(int argc, char *argv[])
{

try{
    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE FOLD OPTIMIZATION IN FORM OF A CROSS" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    // load flat plate [0,1]^2
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/testMeshCrissCross.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial;
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );
    getGeometry( plate, plateGeomInitial );

    // determine boundary mask for reference geometry and foldVertices
    // bdryMaskRef_1 fixes x,y,z in the ref geometry
    std::vector<int> bdryMaskRef_1;
    // bdryMaskRef_2 fixes y,z in the ref geometry
    std::vector<int> bdryMaskRef_2;
    // bdryMaskRef_3 fixes x,z in the ref geometry
    std::vector<int> bdryMaskRef_3;
    RealType t_init_x = 0.00;
    RealType t_init_y = 0.00;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
        // first, fix the outer vertices of the square
        if( (std::abs(coords[0] < 1e-6) || std::abs(coords[0] - 1.0) < 1e-6) && (std::abs(coords[1] < 1e-6) || std::abs(coords[1] - 1.0) < 1e-6) ){
            bdryMaskRef_1.push_back( i );
            continue;
        }
        // Now, fix the middle point of the square
        if(std::abs(coords[0] - 0.625) < 1e-6 && std::abs(coords[1] - 0.625) < 1e-6)
        {
            bdryMaskRef_1.push_back( i );
            coords[0] += t_init_x;
            coords[1] += t_init_y;
            setXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
            continue;
        }
        // Now, fix the vertices on the fold lines that intersect the outside of the square
        if(std::abs(coords[0] - 0.625) < 1e-6 && (std::abs(coords[1]) < 1e-6 || std::abs(coords[1] - 1.0) < 1e-6)){
            bdryMaskRef_1.push_back( i );
            coords[0] += t_init_x;
            setXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
            continue;
        }
        if(std::abs(coords[1] - 0.625) < 1e-6 && (std::abs(coords[0]) < 1e-6 || std::abs(coords[0] - 1.0) < 1e-6)){
            bdryMaskRef_1.push_back( i );
            coords[1] += t_init_y;
            setXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
            continue;
        }
        // now, let all other fold vertices move in x-dir for y=0.5 and in y-dir for x=0.5
        if( std::abs(coords[0] - 0.625) < 1e-6 || std::abs(coords[1] - 0.625) < 1e-6 ){
            if(std::abs(coords[0] - 0.625) < 1e-6)
            {
                bdryMaskRef_3.push_back( i );
                coords[0] += t_init_x;
                setXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
                continue;
            }
            if(std::abs(coords[1] - 0.625) < 1e-6){
                bdryMaskRef_2.push_back( i );
                coords[1] += t_init_y;
                setXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
                continue;
            }
        }
        else{
            // vertices at y = 1.0 and y = 0.0 are allowed to move in x-direction
            if( std::abs(coords[1] - 1.0) < 1e-6 || std::abs(coords[1]) < 1e-6 ){
                bdryMaskRef_2.push_back( i );
            }
            // vertices at x = 1.0 and x = 0.0 are allowed to move in y-direction
            if( std::abs(coords[0] - 1.0) < 1e-6 || std::abs(coords[0]) < 1e-6 ){
                bdryMaskRef_3.push_back( i );
            }
        }
        setXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
    }

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef_1 );

    std::vector<int> activeRef_2 = (std::vector<int>){0,1,1};
    std::vector<int> activeRef_3 = (std::vector<int>){1,0,1};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_2 , activeRef_2);
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_3 , activeRef_3);

    std::vector<VectorType> def_geometries;
    std::vector<VectorType> ref_geometries;
    std::vector<VectorType> fold_DOFs;

    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntries(bdryMaskRef_1.begin(), bdryMaskRef_1.end());
    uniqueEntries.insert(bdryMaskRef_2.begin(), bdryMaskRef_2.end());
    uniqueEntries.insert(bdryMaskRef_3.begin(), bdryMaskRef_3.end());

    // Move the unique elements back to bdryMaskRef_1 in order
    bdryMaskRef_1.assign(uniqueEntries.begin(), uniqueEntries.end());
    std::sort(bdryMaskRef_1.begin(), bdryMaskRef_1.end());

    DirichletSmoother<DefaultConfigurator> smoother(plateGeomInitial, bdryMaskRef_1, plateTopol);
    smoother.apply(plateGeomRef,plateGeomRef);

    setGeometry(plate,plateGeomRef);
    OpenMesh::IO::write_mesh(plate,"reference_mesh.ply");

    // determine boundary mask for optimization
    // and deform part of boundary
    std::vector<int> bdryMaskOpt;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomInitial, coords, i);
        // fix all outer vertices and move them slightly towards the middle of the square
        if(std::abs(coords[0] - 1.0) < 1e-6 && std::abs(coords[1] - 1.0) < 1e-6){
            bdryMaskOpt.push_back( i );
            coords[0] -= 0.2;
            coords[1] -= 0.2;
            coords[2] += 0.2;
        }
        else if(std::abs(coords[0]) < 1e-6 && std::abs(coords[1] - 1.0) < 1e-6){
            bdryMaskOpt.push_back( i );
            coords[0] += 0.2;
            coords[1] -= 0.2;
            coords[2] += 0.2;
        }
        else if(std::abs(coords[0] - 1.0) < 1e-6 && std::abs(coords[1]) < 1e-6){
            bdryMaskOpt.push_back( i );
            coords[0] -= 0.2;
            coords[1] += 0.2;
            coords[2] += 0.2;
        }
        else if(std::abs(coords[0]) < 1e-6 && std::abs(coords[1]) < 1e-6){
            bdryMaskOpt.push_back( i );
            coords[0] += 0.2;
            coords[1] += 0.2;
            coords[2] += 0.2;
        }
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
    }

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );

    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "deformed_mesh.ply");

    VectorType edge_weights;
    FoldDofsCross<DefaultConfigurator> foldDofs( plateTopol, plateGeomInitial, plateGeomRef, bdryMaskRef_1);
    foldDofs.getEdgeWeights(edge_weights);

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

    // saving
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "bendingFoldSol_withNewton_cross.ply");
    }
    catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }

    return 0;
}