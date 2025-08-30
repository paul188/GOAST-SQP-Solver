#include "goast/SQP/Utils/SparseMat.h"
#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <goast/Developability/Developability.h>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACK

int main()
{
    using VectorType = typename DefaultConfigurator::VectorType;
    using MatrixType = typename DefaultConfigurator::SparseMatrixType;
    using VecType = typename DefaultConfigurator::VecType;
    using RealType = typename DefaultConfigurator::RealType;

    try{

        // ---------------------- TOPOLOGICAL PREPARATIONS --------------------------------------

        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/quadMesh2.ply");
        QuadMeshTopologySaver quadTopol(mesh);
        VectorType quadReferenceGeometry;
        QuadMeshTopologySaver::getGeometry(mesh,quadReferenceGeometry);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        // --------------------- GENERATE BASE CENTROID TOPOLOGY ----------------------------

        TriMesh triMesh;
        OpenMesh::IO::read_mesh(triMesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/quadMesh2_refined.ply");
        MeshTopologySaver triangleTopol(triMesh);
        VectorType triBaseGeometry;
        getGeometry(triMesh, triBaseGeometry);

        double tolerance = 1e-4;
        VectorType quadDeformedGeometry = quadReferenceGeometry;

        std::vector<int> bdryMask;
        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            VecType coords;
            coords[0] = quadReferenceGeometry[3*i];
            coords[1] = quadReferenceGeometry[3*i+1];
            coords[2] = quadReferenceGeometry[3*i+2];
            if(std::abs(coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance){
                bdryMask.push_back(3*i);
                bdryMask.push_back(3*i+1);
                bdryMask.push_back(3*i+2);
            }
            if(coords[0] <= 0.5){
                coords[2] = -coords[0]*sqrt(0.34);
                coords[0] = 0.1 + coords[0] *(0.4/0.5);
            }
            else{
                coords[2] = -(1.0-coords[0])*sqrt(0.34);
                coords[0] = 0.5 + (coords[0]-0.5) *(0.4/0.5);
            }
            quadDeformedGeometry[3*i] = coords[0];
            quadDeformedGeometry[3*i+1] = coords[1];
            quadDeformedGeometry[3*i+2] = coords[2];
        }

        VectorType factors_mem_bend(2);
        factors_mem_bend[0] = 10000.0;
        factors_mem_bend[1] = 1.0;

        VectorType factors_elasticity_dev(2);
        factors_elasticity_dev[0] = 1.0;//1.0; // Elasticity factor
        factors_elasticity_dev[1] = 1.0;//10000.0; // Developability factor
    
        QuadElasticEnergy<DefaultConfigurator> quadElasticEnergy(quadTopol, triangleTopol, triBaseGeometry ,factors_elasticity_dev, factors_mem_bend);
        QuadElasticEnergyGradient<DefaultConfigurator> quadElasticEnergyGrad(quadTopol, triangleTopol, triBaseGeometry ,factors_elasticity_dev, factors_mem_bend);
        QuadElasticEnergyHessian<DefaultConfigurator> quadElasticEnergyHess(quadTopol, triangleTopol, triBaseGeometry ,factors_elasticity_dev, factors_mem_bend);

        // set outer optimization parameters
        OptimizationParameters<DefaultConfigurator> optPars;
        //optPars.setGradientIterations( 2000 );
        //optPars.setBFGSIterations( 2000 );
        optPars.setNewtonIterations( 2000 );
        optPars.setEpsilon( 1e-7);
        optPars.setQuietMode( SHOW_ALL );

        VectorType result = quadDeformedGeometry;
        RealType test_1 = 0.0;

        /*
        GradientDescent<DefaultConfigurator> GD( quadElasticEnergy, quadElasticEnergyGrad, optPars);
        GD.setBoundaryMask( bdryMask );
        GD.solve( quadDeformedGeometry, result );

        QuasiNewtonBFGS<DefaultConfigurator> BFGS( quadElasticEnergy, quadElasticEnergyGrad, optPars);
        BFGS.setBoundaryMask( bdryMask );
        BFGS.solve( result, result );
        */

        LineSearchNewton<DefaultConfigurator> Newton( quadElasticEnergy, quadElasticEnergyGrad, quadElasticEnergyHess, optPars);
        Newton.setBoundaryMask( bdryMask );
        Newton.solve( quadDeformedGeometry, result );

        QuadMeshTopologySaver::setGeometry(mesh, result);
        OpenMesh::IO::write_mesh(mesh, "result_quad_reduced.ply");

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}