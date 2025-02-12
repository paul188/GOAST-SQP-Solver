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


int main(int argc, char *argv[])
{
try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "ARC FOLD OPTIMIZATION" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;
    
// load flat plate and prepare the arc crease
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/paperCrissCross.ply");
    MeshTopologySaver plateTopol( plate );
    std::cout<<"num vertices: "<<plateTopol.getNumVertices()<<std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial; // plateGeomRef is the geometry without the arc crease in the mesh
    getGeometry(plate,plateGeomRef);
    getGeometry(plate,plateGeomDef);
    getGeometry(plate,plateGeomInitial);

    // Need two different types of bdrMasks
    /*
    1. bdryMaskRef_i for the reference geometry
    2. bdryMaskOpt_i for the optimization of the deformed geometry
    */

    // Fix all coordinates of the boundary
    std::vector<int> bdryMaskRef_1;
    std::vector<int> bdryMaskOpt_1;
    // Fix only (x,z) coordinates of the boundary
    std::vector<int> bdryMaskRef_2;
    std::vector<int> bdryMaskOpt_2;
    // Fix only y coordinates of the boundary
    std::vector<int> bdryMaskOpt_3;

    RealType t_0 = 1.0; // Initial value of the parameter t

    // Select the vertices with y = 1.0 (topVertices)
    // We will make a Gravitational Force act on them in order to achieve a flipping effect
    std::vector<int> topVertices;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);

        if(std::abs(coords[1] - 1.0) < 1e-6){
            topVertices.push_back(i);
            bdryMaskRef_1.push_back(i);
            continue;
        }

        if(std::abs(coords[1]) < 1e-6)
        {
            bdryMaskRef_1.push_back(i);
            continue;
        }

        if(std::abs(coords[0]) < 1e-6 || std::abs(coords[0] - 1.0) < 1e-6){
            if(std::abs(coords[1] - 0.25) < 1e-6){
                bdryMaskRef_1.push_back(i);
            }
            else{
                bdryMaskRef_2.push_back(i);
            }
            continue;
        }

        if(std::abs(coords[1] - 0.25) < 1e-6){
            bdryMaskRef_1.push_back(i);
            coords[1] += t_0*(0.25 - std::pow(coords[0] - 0.5, 2.0));
            setXYZCoord<VectorType, VecType>(plateGeomRef, coords, i);
            continue;
        }
    }

    // output preliminary meshes
    setGeometry(plate, plateGeomRef);
    OpenMesh::IO::write_mesh(plate, "reference_mesh_before_smoothing.ply");

    // activeRef helper vectors
    std::vector<int> activeRef_2 =  (std::vector<int>){1,0,1};
    std::vector<int> activeRef_3 = (std::vector<int>) {0,1,0};

    // ----------------- Smooth reference geometry -----------------------------

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef_1 );
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_2 , activeRef_2);
    
    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntries(bdryMaskRef_1.begin(), bdryMaskRef_1.end());
    uniqueEntries.insert(bdryMaskRef_2.begin(), bdryMaskRef_2.end());

    // Move the unique elements back to bdryMaskRef in order
    std::vector<int> bdryMaskRef;
    bdryMaskRef.assign(uniqueEntries.begin(), uniqueEntries.end());
    std::sort(bdryMaskRef.begin(), bdryMaskRef.end());

    DirichletSmoother<DefaultConfigurator> smoother(plateGeomInitial, bdryMaskRef, plateTopol);
    smoother.apply(plateGeomRef,plateGeomRef);

    setGeometry(plate,plateGeomRef);
    OpenMesh::IO::write_mesh(plate,"reference_mesh.ply");

    // initialize deformed geometry
    plateGeomDef = plateGeomRef;
    for(int i = 0; i < plateTopol.getNumVertices(); i++){
        VecType coords;
        getXYZCoord<VectorType, VecType>(plateGeomInitial, coords, i);
        if((std::abs(coords[0]) < 1e-6) && (std::abs(coords[1]) <= 0.25))
        {
            if(std::abs(coords[1] - 0.25) < 1e-6 || std::abs(coords[1]) < 1e-6)
            {
                bdryMaskOpt_1.push_back(i);
            }
            else{
                bdryMaskOpt_2.push_back(i);
            }
            coords[0] += 1.0/(2*M_PI);
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>(plateGeomDef, coords, i);
            continue;
        }
        if((std::abs(coords[0] - 1.0) < 1e-6) && (std::abs(coords[1]) <= 0.25))
        {
            if(std::abs(coords[1] - 0.25) < 1e-6 || std::abs(coords[1]) < 1e-6)
            {
                bdryMaskOpt_1.push_back(i);
            }
            else{
                bdryMaskOpt_2.push_back(i);
            }
            coords[0] -= 1.0/(2*M_PI);
            coords[2] += 0.1;
            setXYZCoord<VectorType, VecType>(plateGeomDef, coords, i);
            continue;
        }
    }

    setGeometry(plate,plateGeomDef);
    OpenMesh::IO::write_mesh(plate,"deformed_mesh.ply");

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt_1 );
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskOpt_2 , activeRef_2);

    uniqueEntries.clear();
    uniqueEntries.insert(bdryMaskOpt_1.begin(), bdryMaskOpt_1.end());
    uniqueEntries.insert(bdryMaskOpt_2.begin(), bdryMaskOpt_2.end());

    std::vector<int> bdryMaskOpt;
    bdryMaskOpt.assign(uniqueEntries.begin(), uniqueEntries.end());
    std::sort(bdryMaskOpt.begin(), bdryMaskOpt.end());

    auto foldDofsPtr = std::make_shared<FoldDofsArcLine<DefaultConfigurator>>(plateTopol,plateGeomInitial, plateGeomRef, bdryMaskRef);

    std::vector<int> foldVertices;
    foldDofsPtr -> getFoldVertices(foldVertices);

    std::cout<<"Boundary opt test 1: "<<std::endl;
    for(int i = 0; i < bdryMaskOpt.size(); i++)
    {
        VecType coords;
        if(bdryMaskOpt[i] < plateTopol.getNumVertices())
        {
            getXYZCoord<VectorType, VecType>(plateGeomInitial, coords,bdryMaskOpt[i]);
            std::cout<<"x fixed: "<<i<<": "<<coords[0]<<", "<<coords[1]<<", "<<coords[2]<<std::endl;
        }
        else if(bdryMaskOpt[i] >= plateTopol.getNumVertices() && bdryMaskOpt[i] < 2*plateTopol.getNumVertices())
        {
            getXYZCoord<VectorType, VecType>(plateGeomRef, coords, bdryMaskOpt[i] - plateTopol.getNumVertices());
            std::cout<<"y fixed: "<<i<<": "<<coords[0]<<", "<<coords[1]<<", "<<coords[2]<<std::endl;
        }
        else if(bdryMaskOpt[i] >= 2*plateTopol.getNumVertices())
        {
            getXYZCoord<VectorType, VecType>(plateGeomRef, coords, bdryMaskOpt[i] - 2*plateTopol.getNumVertices());
            std::cout<<"z fixed: "<<i<<": "<<coords[0]<<", "<<coords[1]<<", "<<coords[2]<<std::endl;
        }
    }

    std::cout<<"End boundary vertex test 1"<<std::endl;


    auto DfoldDofsPtr = std::make_shared<FoldDofsArcLineGradient<DefaultConfigurator>>(plateTopol, bdryMaskRef, plateGeomInitial, foldVertices);

    VectorType edge_weights = VectorType::Zero(plateTopol.getNumEdges());
    foldDofsPtr->getEdgeWeights(edge_weights);

    size_t nFoldDOFs = foldDofsPtr->getNumDofs();
    size_t nVertexDOFs = 3*plateTopol.getNumVertices();

    std::cout<<"Number of fold DOFs: "<<nFoldDOFs<<std::endl;
    
    SQPLineSearchParams<DefaultConfigurator> pars;
    CostFunctional<DefaultConfigurator> costFunctional(foldVertices);
    CostFunctionalGradient<DefaultConfigurator> DcostFunctional(plateTopol,foldVertices);
    
    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;
    RealType factor_gravity = 1.0;

    VectorType factors(3);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;
    factors[2] = factor_gravity;

    VectorType mass_distribution = VectorType::Zero(plateTopol.getNumVertices());

    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>(plateGeomRef, coords, i);
        if(std::abs(coords[1] - 1.0) < 1e-4)
        {
            mass_distribution[i] = 1.0;
        }
    }

    VectorType gravity_dir(3);
    gravity_dir << 0.0, 0.0, -1.0;

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