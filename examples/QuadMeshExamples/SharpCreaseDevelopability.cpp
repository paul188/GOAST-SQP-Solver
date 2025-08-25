#include <goast/Developability/Developability.h>
#include <goast/Core.h>
#include <random>
#include <iostream>
#include <utility>
#include <goast/Smoothers.h>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACK

int main()
{
    using VectorType = typename DefaultConfigurator::VectorType;
    using MatrixType = typename DefaultConfigurator::SparseMatrixType;
    using VecType = typename DefaultConfigurator::VecType;
    using RealType = typename DefaultConfigurator::RealType;

    try{

        /*
        MyMesh mesh_test;
        OpenMesh::IO::read_mesh(mesh_test,"/home/s24pjoha_hpc/goast_old_old/goast/build/examples/smoothed_quad_newest.ply");
        VectorType testGeometry;
        QuadMeshTopologySaver quadTopolTest(mesh_test);
        QuadMeshTopologySaver::getGeometry(mesh_test, testGeometry);
        ConstraintSqrdReduced<DefaultConfigurator> constraintSqrdReduced_test(quadTopolTest);
        RealType value = 0;
        constraintSqrdReduced_test.apply(testGeometry, value);
        std::cout<<"ConstraintSqrdReduced value: "<<value<<std::endl;

        ConstraintSqrdReducedHessian<DefaultConfigurator> constraintSqrdReducedHess_test(quadTopolTest);
        std::cout<<"test geometry size: "<<testGeometry.size()<<std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        MatrixType hess_test(testGeometry.size(), testGeometry.size());
        constraintSqrdReducedHess_test.apply(testGeometry, hess_test);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "calc_hessian took " << elapsed.count() << " seconds." << std::endl;

        std::cout<<hess_test.rows()<<" x "<<hess_test.cols()<<std::endl;

        for(int i = 0; i < hess_test.outerSize(); i++){
            for(typename MatrixType::InnerIterator it(hess_test,i); it; ++it)
            {
                std::cout<<"Hess entry: ("<<it.row()<<","<<it.col()<<") = "<<it.value()<<std::endl;
            }
        }*/

        // ---------------------- TOPOLOGICAL PREPARATIONS --------------------------------------

        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/basePlateSharpCreaseDevelopability.ply");
        QuadMeshTopologySaver quadTopol(mesh);
        VectorType quadDeformedGeometry, quadReferenceGeometry;
        QuadMeshTopologySaver::getGeometry(mesh,quadDeformedGeometry);
        QuadMeshTopologySaver::getGeometry(mesh,quadReferenceGeometry);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        // Read in the initial deformed quad geometry
        OpenMesh::IO::read_mesh(mesh, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/smoothed_quad_newest.ply");
        QuadMeshTopologySaver::getGeometry(mesh, quadDeformedGeometry);

        VectorType factors(2);
        factors[0] = 0.0; 
        factors[1] = 0.0; 

        int counter_1 = 0;
        int counter_2 = 0;
        std::vector<int> bdryMask;
        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            VecType coords;
            coords[0] = quadReferenceGeometry[3*i];
            coords[1] = quadReferenceGeometry[3*i+1];
            coords[2] = quadReferenceGeometry[3*i+2];

            if((std::abs(coords[0]) < 1e-4) || (std::abs(1.0 - coords[0]) < 1e-3) )
            {
                counter_1++;
                bdryMask.push_back(3*i);
                bdryMask.push_back(3*i + 1);
                bdryMask.push_back(3*i + 2);
            }

            if((std::abs(coords[1]) < 1e-4) || (std::abs(1.0 - coords[1]) < 1e-3) )
            {
                counter_2++;
                bdryMask.push_back(3*i);
                bdryMask.push_back(3*i + 1);
                bdryMask.push_back(3*i + 2);
            }
        }

        /*
        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            VecType coords;
            coords[0] = quadDeformedGeometry[3*i];
            coords[1] = quadDeformedGeometry[3*i+1];
            coords[2] = quadDeformedGeometry[3*i+2];

            coords[1] *= sqrt(1.0 - (0.4*0.4));
            coords[0] *= sqrt(1.0 - (0.4*0.4));

            quadDeformedGeometry[3*i] = coords[0];
            quadDeformedGeometry[3*i + 1] = coords[1];
            quadDeformedGeometry[3*i + 2] = coords[2];
        }*/

        QuadMeshTopologySaver::setGeometry(mesh, quadDeformedGeometry);
        OpenMesh::IO::write_mesh(mesh, "new_2.ply");

        EdgeLengthQuadEnergy<DefaultConfigurator> E_edge(quadTopol, quadReferenceGeometry);
        EdgeLengthQuadEnergyGradient<DefaultConfigurator> DE_edge(quadTopol, quadReferenceGeometry);

        auto constraintSqrdReducedPtr = std::make_shared<ConstraintSqrdReduced<DefaultConfigurator>>(quadTopol);
        auto constraintSqrdReducedGradPtr = std::make_shared<ConstraintSqrdReducedGradient<DefaultConfigurator>>(quadTopol);

        ConstraintSqrdReduced<DefaultConfigurator> constraintSqrdReduced(quadTopol);
        ConstraintSqrdReducedGradient<DefaultConfigurator> constraintSqrdReducedGrad(quadTopol);
        ConstraintSqrdReducedHessian<DefaultConfigurator> constraintSqrdReducedHess(quadTopol);

        // TRY THE COMBINED ENERGY
        VectorType factors_elasticity_dev(2);
        factors_elasticity_dev[0] = 0.0; // Elasticity factor
        factors_elasticity_dev[1] = 10000.0; // Developability factor

        VectorType factors_mem_bend(2);
        factors_mem_bend[0] = 10000.0;
        factors_mem_bend[1] = 1.0;

        QuadMeshTopologySaver::setGeometry(mesh, quadReferenceGeometry);
        OpenMesh::IO::write_mesh(mesh, "quad_base.ply");

        TriMesh centroidMesh = quadTopol.makeQuadMeshCentroid();
        MeshTopologySaver centroidTopol = MeshTopologySaver(centroidMesh);
        VectorType centroidReferenceGeometry;

        auto E_bend = std::make_shared<SimpleBendingEnergy<DefaultConfigurator>>(centroidTopol, centroidReferenceGeometry, true);
        auto E_mem = std::make_shared<NonlinearMembraneEnergy<DefaultConfigurator>>(centroidTopol, centroidReferenceGeometry, true);

        auto DE_bend = std::make_shared<SimpleBendingGradientDef<DefaultConfigurator>>(centroidTopol, centroidReferenceGeometry, true);
        auto DE_mem = std::make_shared<NonlinearMembraneGradientDef<DefaultConfigurator>>(centroidTopol, centroidReferenceGeometry, true);


        getGeometry(centroidMesh, centroidReferenceGeometry);
        CoordConverter<DefaultConfigurator> coordConverter(quadTopol, centroidTopol, quadReferenceGeometry, centroidReferenceGeometry);
        CombinedQuadEnergy<DefaultConfigurator> E_combined(quadTopol, centroidTopol, quadReferenceGeometry, centroidReferenceGeometry, mesh, centroidMesh, factors_elasticity_dev, factors_mem_bend, coordConverter, std::move(E_bend), std::move(E_mem), std::move(constraintSqrdReducedPtr));
        CombinedQuadEnergyGradient<DefaultConfigurator> DE_combined(quadTopol, centroidTopol, quadReferenceGeometry, centroidReferenceGeometry, mesh, centroidMesh, factors_elasticity_dev, factors_mem_bend, coordConverter, std::move(DE_bend), std::move(DE_mem), std::move(constraintSqrdReducedGradPtr));

        VectorType factors_additionop(2);
        factors_additionop[0] = 1.0; // Edge length factor
        factors_additionop[1] = 1.0; // Constraint factor

        AdditionOp<DefaultConfigurator> E_tot( factors_additionop, E_edge, constraintSqrdReduced );
        AdditionGradient<DefaultConfigurator> DE_tot( factors_additionop, DE_edge, constraintSqrdReducedGrad);

        RealType test_energy;
        E_combined.apply(quadDeformedGeometry, test_energy);
        std::cout<<"Test combined energy: "<<test_energy<<std::endl;
        VectorType test_gradient;
        DE_combined.apply(quadDeformedGeometry, test_gradient);
        std::cout<<"Test combined gradient: "<<test_gradient.norm()<<std::endl;
        
        // set outer optimization parameters
        OptimizationParameters<DefaultConfigurator> optPars;
        optPars.setGradientIterations( 200 );
        optPars.setBFGSIterations( 200 );
        optPars.setNewtonIterations( 200 );
        optPars.setEpsilon( 1e-10);
        optPars.setQuietMode( SHOW_ALL );

        VectorType result = quadDeformedGeometry;


        /*
        GradientDescent<DefaultConfigurator> GD( constraintSqrdReduced, DE_combined, optPars);
        GD.setBoundaryMask( bdryMask );
        GD.solve( quadDeformedGeometry, result );

        QuasiNewtonBFGS<DefaultConfigurator> BFGS( E_combined, DE_combined, optPars);
        BFGS.setBoundaryMask( bdryMask );
        BFGS.solve( result, result );
        */

        LineSearchNewton<DefaultConfigurator> newton( constraintSqrdReduced, constraintSqrdReducedGrad, constraintSqrdReducedHess, optPars);
        newton.setBoundaryMask( bdryMask );
        newton.solve( quadDeformedGeometry, result );

        QuadMeshTopologySaver::setGeometry(mesh, result);
        // was centroid_base_deformed_newest.ply for the original without Edge length contributions
        OpenMesh::IO::write_mesh(mesh, "centroid_base_deformed_newest_with_edges.ply");
        OpenMesh::IO::write_mesh(centroidMesh, "/lustre/scratch/data/s24pjoha_hpc-results/thesis_results/developability/result.ply");

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}