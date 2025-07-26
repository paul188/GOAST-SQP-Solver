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
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
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
        /*
        else if(std::abs(coords[0] - 0.5) < 1e-6 && std::abs(coords[1] - 0.5) < 1e-6)
        {
            coords[2] -= 0.2;
        }
            */
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
    }

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );

    setGeometry(plate,plateGeomDef);
    OpenMesh::IO::write_mesh(plate,"deformed_mesh.ply");

    auto foldDofsPtr = std::make_shared<FoldDofsCross<DefaultConfigurator>>(plateTopol,plateGeomInitial, plateGeomRef, bdryMaskRef_1);
    
    std::vector<int> foldVertices;
    foldDofsPtr -> getFoldVertices(foldVertices);
    std::vector<int> line1Vertices;
    foldDofsPtr -> getLine1Vertices(line1Vertices);
    std::vector<int> line2Vertices;
    foldDofsPtr -> getLine2Vertices(line2Vertices);

    auto DfoldDofsPtr = std::make_shared<FoldDofsCrossGradient<DefaultConfigurator>>(plateTopol, bdryMaskRef_1, plateGeomInitial, foldVertices, line1Vertices, line2Vertices);

    VectorType edge_weights = VectorType::Zero(plateTopol.getNumEdges());
    foldDofsPtr->getEdgeWeights(edge_weights);

    size_t nFoldDOFs = foldDofsPtr->getNumDofs();
    size_t nVertexDOFs = 3*plateTopol.getNumVertices();
    
    SQPLineSearchParams<DefaultConfigurator> pars;
    CostFunctional<DefaultConfigurator> costFunctional(foldVertices);
    CostFunctionalGradient<DefaultConfigurator> DcostFunctional(plateTopol,foldVertices);
    
    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;

    VectorType factors(2);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;

    auto factory = std::make_shared<ElasticEnergyFactory<DefaultConfigurator>>(factors, plateTopol, edge_weights);
    BoundaryDOFS<DefaultConfigurator> boundaryDOFs(bdryMaskOpt, nVertexDOFs, nFoldDOFs);
    // Create the degrees of freedom object
    std::vector<RealType> deviations;
    VectorType vertexDOFs = VectorType::Zero(3*plateTopol.getNumVertices());
    ProblemDOFs<DefaultConfigurator> problemDOFs(VectorType::Zero(nFoldDOFs), vertexDOFs, plateGeomDef, foldDofsPtr, DfoldDofsPtr);
    SQPLineSearchSolver<DefaultConfigurator> solver(pars, costFunctional, DcostFunctional, factory, boundaryDOFs, problemDOFs, 20);
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

    std::ofstream outFile("output_mu_nonmonotone.txt");
    if (!outFile) {
        std::cerr << "Error: Could not open the file for writing.\n";
        return 1;
    }

    // Write the vector to the file
    for (RealType value : deviations) {
        outFile << value << "\n"; // Write each value on a new line
    }

    // Close the file
    outFile.close();
    std::cout<<"File written successfully"<<std::endl;

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}