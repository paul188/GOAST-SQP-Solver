// TO BE DONE

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

#include <goast/SQP/Algorithm/CostFunctional.h>
#include <goast/SQP/Utils/SparseMat.h>
#include <goast/SQP/DOFHandling/FoldDofs.h>

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

  using MatrixType = DefaultConfigurator::SparseMatrixType;
  using FullMatrixType = DefaultConfigurator::FullMatrixType;

try{
    // load flat plate [0,1]^2
    TriMesh plate;
    OpenMesh::IO::read_mesh(plate, "../../data/plate/paperCrissCross.ply");
    MeshTopologySaver plateTopol( plate );
    std::cerr << "num of nodes = " << plateTopol.getNumVertices() << std::endl;
    VectorType plateGeomRef, plateGeomDef, plateGeomInitial;
    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );
    getGeometry( plate, plateGeomInitial );

    // Fix all coordinates of the boundary
    std::vector<int> bdryMaskRef_all;

    RealType t_0 = 0.25; // Initial value of the parameter t

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
            bdryMaskRef_all.push_back(i);
            if(is_near(y,1.0))
            {
                topVertices.push_back(i);
            }
            continue;
        }

        if(is_near(y,0.25)){
            bdryMaskRef_all.push_back(i);
            coords[1] += t_0*(0.25 - std::pow(coords[0] - 0.5, 2.0));
            setXYZCoord<VectorType, VecType>(plateGeomRef, coords, i);
            continue;
        }
    }

    // output preliminary meshes
    setGeometry(plate, plateGeomRef);
    OpenMesh::IO::write_mesh(plate, "reference_mesh_before_smoothing.ply");

    // ----------------- Smooth reference geometry -----------------------------

    extendBoundaryMask( plateTopol.getNumVertices(), bdryMaskRef_all );
    //extendBoundaryMaskPartial( plateTopol.getNumVertices(), bdryMaskRef_xz , activeRef_xz);
    
    // Use an unordered_set to remove duplicates
    std::unordered_set<int> uniqueEntriesRef(bdryMaskRef_all.begin(), bdryMaskRef_all.end());
    //uniqueEntriesRef.insert(bdryMaskRef_xz.begin(), bdryMaskRef_xz.end());

    // Move the unique elements back to bdryMaskRef in order
    std::vector<int> bdryMaskRef;
    bdryMaskRef.assign(uniqueEntriesRef.begin(), uniqueEntriesRef.end());
    std::sort(bdryMaskRef.begin(), bdryMaskRef.end());

    DirichletSmoother<DefaultConfigurator> smoother(plateGeomInitial, bdryMaskRef, plateTopol);
    smoother.apply(plateGeomRef,plateGeomRef);

    setGeometry(plate, plateGeomRef);
    OpenMesh::IO::write_mesh(plate, "test_ref.ply");

    std::vector<int> foldVertices;

    auto foldDofs = FoldDofsArcLine<DefaultConfigurator>(plateTopol,plateGeomInitial, plateGeomRef, bdryMaskRef);
    foldDofs.getFoldVertices(foldVertices);
    auto DfoldDofs = FoldDofsArcLineGradient<DefaultConfigurator>(plateTopol, bdryMaskRef, plateGeomInitial, foldVertices);

    // TEST THE FOLD DOFS ACTION
    VectorType translate_vec = VectorType::Ones(foldDofs.getNumDofs());
    translate_vec[0] = 0.5;
    VectorType testGeom = plateGeomInitial;
    MatrixType testMat;
    foldDofs.apply(translate_vec,testGeom);
    DfoldDofs.apply(translate_vec,testMat);
    setGeometry(plate, testGeom);
    OpenMesh::IO::write_mesh(plate, "deformed_mesh_test_new.ply");
    setGeometry(plate, plateGeomInitial);
    OpenMesh::IO::write_mesh(plate, "deformed_mesh_test_new_2.ply");
    // END TEST

    // Test just applying it:
    MatrixType Dest;
    VectorType Dest_vec;
    foldDofs.apply(VectorType::Zero(foldDofs.getNumDofs()), Dest_vec);
    DfoldDofs.apply(VectorType::Zero(foldDofs.getNumDofs()), Dest);

    VectorValuedDerivativeTester<DefaultConfigurator> tester3(foldDofs, DfoldDofs, 0.005, Dest.rows());
    VectorType translationDof2 = VectorType::Zero(foldDofs.getNumDofs());
    translationDof2[0] = 0.5;
    tester3.plotAllDirections(translationDof2, "deriv_test");


}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}

    return 0;
}