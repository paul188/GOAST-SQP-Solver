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

        TriMesh centroid_base_mesh = quadTopol.makeQuadMeshCentroid();
        MeshTopologySaver centroidTopol(centroid_base_mesh);
        VectorType centroidReferenceGeometry;
        getGeometry(centroid_base_mesh, centroidReferenceGeometry);
        OpenMesh::IO::write_mesh(centroid_base_mesh, "centroid_base_coarse.ply");

        VectorType centroidDeformedGeometry = centroidReferenceGeometry;
        double tolerance = 1e-4;

        std::vector<int> bdryMask;
        for(int i = 0; i < centroidTopol.getNumVertices(); i++)
        {
            VecType coords;
            getXYZCoord<VectorType, VecType>(centroidReferenceGeometry,coords, i);
            if(std::abs(coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance){
                bdryMask.push_back(i);
            }
            if(coords[0] <= 0.5){
                coords[2] = -coords[0]*sqrt(0.34);
                coords[0] = 0.1 + coords[0] *(0.4/0.5);
            }
            else{
                coords[2] = -(1.0-coords[0])*sqrt(0.34);
                coords[0] = 0.5 + (coords[0]-0.5) *(0.4/0.5);
            }
            setXYZCoord<VectorType, VecType>(centroidDeformedGeometry,coords, i);
        }

        setGeometry(centroid_base_mesh, centroidDeformedGeometry);
        OpenMesh::IO::write_mesh(centroid_base_mesh, "centroid_base_deformed_0.ply");

        extendBoundaryMask(centroidTopol.getNumVertices(), bdryMask);

        VectorType factors_mem_bend(2);
        factors_mem_bend[0] = 0.0;//10000.0;
        factors_mem_bend[1] = 1.0;

        VectorType factors_elasticity_dev(2);
        factors_elasticity_dev[0] = 0.0;//1.0; // Elasticity factor
        factors_elasticity_dev[1] = 1.0;//10000.0; // Developability factor

        // Create random number generator
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(-0.05, 0.05); // range: [-0.05, 0.05]

        QuadElasticEnergy<DefaultConfigurator> quadElasticEnergy(quadTopol, centroidTopol, quadReferenceGeometry, centroidReferenceGeometry,factors_elasticity_dev, factors_mem_bend);
        QuadElasticEnergyGradient<DefaultConfigurator> quadElasticEnergyGrad(quadTopol, centroidTopol, quadReferenceGeometry, centroidReferenceGeometry, factors_elasticity_dev, factors_mem_bend);

        SimpleBendingEnergy<DefaultConfigurator> E_bend(centroidTopol, centroidReferenceGeometry, true);
        SimpleBendingGradientDef<DefaultConfigurator> DE_bend(centroidTopol, centroidReferenceGeometry);

        NonlinearMembraneEnergy<DefaultConfigurator> E_mem(centroidTopol, centroidReferenceGeometry,true);
        NonlinearMembraneGradientDef<DefaultConfigurator> DE_mem(centroidTopol, centroidReferenceGeometry);

        AdditionOp<DefaultConfigurator> E_tot( factors_mem_bend, E_mem, E_bend);
        AdditionGradient<DefaultConfigurator> DE_tot( factors_mem_bend, DE_mem, DE_bend);
        
        // set outer optimization parameters
        OptimizationParameters<DefaultConfigurator> optPars;
        optPars.setGradientIterations( 200 );
        optPars.setBFGSIterations( 2000 );
        optPars.setEpsilon( 1e-10);
        optPars.setQuietMode( SHOW_ALL );

        VectorType result_centroid = centroidDeformedGeometry;

        GradientDescent<DefaultConfigurator> GD( quadElasticEnergy, quadElasticEnergyGrad, optPars);
        GD.setBoundaryMask( bdryMask );
        GD.solve( centroidDeformedGeometry, result_centroid );

        QuasiNewtonBFGS<DefaultConfigurator> BFGS( quadElasticEnergy, quadElasticEnergyGrad, optPars);
        BFGS.setBoundaryMask( bdryMask );
        BFGS.solve( result_centroid, result_centroid );

        setGeometry(centroid_base_mesh, result_centroid);
        OpenMesh::IO::write_mesh(centroid_base_mesh, "centroid_base_deformed.ply");

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}