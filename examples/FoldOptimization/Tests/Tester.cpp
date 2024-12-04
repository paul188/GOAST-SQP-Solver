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

#include "../Headers/SQP/CostFunctional.h"
#include "../Headers/SQP/Merit.h"
#include "../Headers/Utils/SparseMat.h"

#include "../Headers/DOFHandling/BoundaryDOFS.h"
#include "../Headers/DOFHandling/FoldDofs.h"

#include "../Headers/Utils/ObjectFactory.h"
#include "D2E_dt.h"

typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;
typedef DefaultConfigurator::SparseMatrixType MatrixType;
typedef DefaultConfigurator::FullMatrixType FullMatrixType;

int main(int argc, char *argv[])
{

  using MatrixType = DefaultConfigurator::SparseMatrixType;
  using FullMatrixType = DefaultConfigurator::FullMatrixType;

try{
    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE FOLD OPTIMIZATION" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    // load flat plate [0,1]^2
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/testMesh.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial;
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );
    getGeometry( plate, plateGeomInitial );

    // determine boundary mask for reference geometry and foldVertices
    // the part of the reference boundary where everything is fixed
    std::vector<int> bdryMaskRef_1;
    std::vector<int> bdryMaskRef_2;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);
        if( std::abs(coords[1]) < 1e-6 || std::abs(coords[1] - 4.0) < 1e-10 ){
            bdryMaskRef_1.push_back( i );
        }
        // ISSUE HERE WITH VALUES NOT RECOGNIZED
        if(std::abs(coords[1] - 2.3) < 1e-6){
            bdryMaskRef_1.push_back( i );
        }
        if( std::abs(coords[0]) < 1e-6 || std::abs(coords[0] - 1.0) < 6 ){
            bdryMaskRef_2.push_back( i );
        }
    }

    // extend all boundary masks to (x,y,z) coordinates

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef_1 );
    std::vector<int> activeRef_2 = (std::vector<int>){1,0,1};
    extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_2 , activeRef_2);

    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntries(bdryMaskRef_1.begin(), bdryMaskRef_1.end());
    uniqueEntries.insert(bdryMaskRef_2.begin(), bdryMaskRef_2.end());

    // Move the unique elements back to bdryMaskRef_1 in order
    bdryMaskRef_1.assign(uniqueEntries.begin(), uniqueEntries.end());
    std::sort(bdryMaskRef_1.begin(), bdryMaskRef_1.end());

     // determine boundary mask for optimization
     // and deform part of boundary
    std::vector<int> bdryMaskOpt;
    for( int i = 0; i < plateTopol.getNumVertices(); i++ ){
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
        if( coords[1] == 0.0 || coords[1] == 4.0 ){
            bdryMaskOpt.push_back( i );
            if(coords[1] == 0.0){
                coords[1] += 1.5;
            }
            else{
                coords[1] -= 1.5;
            }
        }
        
        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
    }

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskOpt );

    // Test the FoldDofs derivative -> Is successful

    FoldDofsSimpleLine<DefaultConfigurator> foldDofs(plateTopol,plateGeomInitial, plateGeomRef, bdryMaskRef_1);
    std::vector<int> foldVertices;
    foldDofs.getFoldVertices(foldVertices);
    FoldDofsSimpleLineGradient<DefaultConfigurator> foldDofsGradient(plateTopol, bdryMaskRef_1, plateGeomInitial, foldVertices);
    VectorValuedDerivativeTester<DefaultConfigurator> tester(foldDofs, foldDofsGradient, 0.005, 3*plateTopol.getNumVertices());
    VectorType translationDof = VectorType::Zero(1);
    translationDof[0] = -0.5;
    tester.plotAllDirections(translationDof, "foldDofsGradient");

    // Test the CostFunctional derivative -> Is successful
    CostFunctional<DefaultConfigurator> costFunctional(foldVertices);
    CostFunctionalGradient<DefaultConfigurator> DcostFunctional(plateTopol, foldVertices);
    ScalarValuedDerivativeTester<DefaultConfigurator> tester2(costFunctional, DcostFunctional, 0.005, 200);
    tester2.plotAllDirections(plateGeomDef, "costFunctionalGradient");

    // Test second derivative w.r.t. t
    BoundaryDOFS<DefaultConfigurator> boundaryDOFs(bdryMaskOpt, 3*plateTopol.getNumVertices(), 1);
    auto foldDofs2 = std::make_unique<FoldDofsSimpleLine<DefaultConfigurator>>(plateTopol,plateGeomInitial, plateGeomRef, bdryMaskRef_1);
    VectorType edge_weights = VectorType::Zero(plateTopol.getNumEdges());
    foldDofs2 -> getEdgeWeights(edge_weights);
    auto DFoldDofs = std::make_unique<FoldDofsSimpleLineGradient<DefaultConfigurator>>(plateTopol, bdryMaskRef_1, plateGeomInitial, foldVertices);
    MyObjectFactory<DefaultConfigurator> objectFactory;
    EnergyFoldDof<DefaultConfigurator> DE_t(plateTopol, plateGeomDef, edge_weights, objectFactory, std::move(foldDofs2), boundaryDOFs);
    auto foldDofs3 = std::make_unique<FoldDofsSimpleLine<DefaultConfigurator>>(plateTopol,plateGeomInitial, plateGeomRef, bdryMaskRef_1);
    EnergyGradientFoldDof<DefaultConfigurator> D2E_t(plateTopol, plateGeomDef, edge_weights, std::move(foldDofs3), std::move(DFoldDofs), objectFactory, boundaryDOFs);
    VectorValuedDerivativeTester<DefaultConfigurator> tester3(DE_t, D2E_t, 0.005, 3*plateTopol.getNumVertices() - bdryMaskOpt.size());
    VectorType translationDof2 = VectorType::Zero(1);
    translationDof2[0] = 0.5;
    tester3.plotAllDirections(translationDof2, "MixedEnergyTGrad");

    // Test the L1 Merit function: -> works
    L1Merit<DefaultConfigurator> merit(0.5);
    RealType dest = 0;
    VectorType testvec;
    testvec.resize(4);
    testvec << -1,2.4,-3.1,0;
    merit.apply(2.3,testvec,dest);
    std::cout<<"Result of L1Merit: "<<dest<<std::endl;

    MatrixType testMat(4,4);
    testMat.setZero();

    // Test the L1 Merit derivative: -> works
    L1MeritGrad<DefaultConfigurator> Dmerit(0.5, testMat);
    RealType dest2 = 0;
    VectorType testVec2;
    testVec2.resize(4);
    testVec2 << 1,2,3,4;
    VectorType testVec3;
    testVec3.resize(4);
    testVec3 << 5,-2.5,3,1;
    Dmerit.apply(testvec,testVec2,testVec3,dest2);

    std::cout<<dest2<<std::endl;


}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}

    return 0;
}