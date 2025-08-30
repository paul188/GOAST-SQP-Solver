#include <goast/Developability/Developability.h>
#include <goast/Core.h>
#include <random>
#include <iostream>
#include <utility>
#include <goast/Smoothers.h>
#include <goast/DiscreteShells.h>

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
        VectorType quadBaseGeometry;
        QuadMeshTopologySaver::getGeometry(mesh,quadBaseGeometry);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            std::cout<<"vertex: "<<i<<": "<<quadBaseGeometry[3*i]<<", "<<quadBaseGeometry[3*i+1]<<", "<<quadBaseGeometry[3*i+2]<<std::endl;
        }

        // --------------------- GENERATE BASE CENTROID TOPOLOGY ----------------------------

        TriMesh triMesh;
        OpenMesh::IO::read_mesh(triMesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/quadMesh2_refined.ply");
        MeshTopologySaver triangleTopol(triMesh);
        VectorType triBaseGeometry;
        getGeometry(triMesh, triBaseGeometry);

        for(int i = 0; i < triangleTopol.getNumVertices(); i++)
        {
            VecType coords;
            getXYZCoord<VectorType, VecType>(triBaseGeometry,coords, i);
            std::cout<<"triangle vertex: "<<i<<": "<<coords[0]<<", "<<coords[1]<<", "<<coords[2]<<std::endl;
        }

        //Test: 
        VectorType test_1(3*triangleTopol.getNumVertices());
        QuadTriConverter<DefaultConfigurator> quadTriConverter(quadTopol.getNumVertices());
        quadTriConverter.convertQuadToTri(quadBaseGeometry, test_1);
        std::cout<<"Difference in conversion of base geometries: "<<(test_1 - triBaseGeometry).norm()<<std::endl;
        quadTriConverter.convertTriToQuad(triBaseGeometry, test_1);
        std::cout<<"Difference in conversion of base geometries back: "<<(test_1 - quadBaseGeometry).norm()<<std::endl;

        // ---------------------  GENERATE MAPS BETWEEN QUAD AND CENTROID TOPOLOGIES ----------------------------

        // -------------------------- SETUP BOUNDARY ------------------------------------------

        VarsIdx<DefaultConfigurator> varsIdx(quadTopol);

        double tolerance = 1e-4;

        VectorType quadDeformedGeometry = quadBaseGeometry;

        std::vector<int> bdryMask;

        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            VecType coords;
            coords[0] = quadBaseGeometry[3*i];
            coords[1] = quadBaseGeometry[3*i+1];
            coords[2] = quadBaseGeometry[3*i+2];
            if(std::abs(coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance){
                bdryMask.push_back(varsIdx["vertices"] + 3*i);
                bdryMask.push_back(varsIdx["vertices"] + 3*i + 1);
                bdryMask.push_back(varsIdx["vertices"] + 3*i + 2);
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
            quadDeformedGeometry[3*i + 1] = coords[1];
            quadDeformedGeometry[3*i + 2] = coords[2];
        }

        /*
        for( int i = 0; i < quadTopol.getNumVertices(); i++ ){
            VecType coords;
            coords[0] = quadBaseGeometry[3*i];
            coords[1] = quadBaseGeometry[3*i+1];
            coords[2] = quadBaseGeometry[3*i+2];
            if( coords[0] < 0.01 || coords[0] > 0.99 ){
                bdryMask.push_back( i );
            // deform part of boundary
                if( coords[0] > 0.99 ) {
                    coords[0] -= 0.2;
                }
                if(coords[0] < 0.01){
                    coords[0] += 0.2;
                }
                quadDeformedGeometry[3*i] = coords[0]; 
                quadDeformedGeometry[3*i + 1] = coords[1];
                quadDeformedGeometry[3*i + 2] = coords[2];
            }

            if((coords[1] < 0.01 || coords[1] > 0.99))
            {
                bdryMask.push_back(varsIdx["vertices"] + 3*i);
                bdryMask.push_back(varsIdx["vertices"] + 3*i + 1);
                bdryMask.push_back(varsIdx["vertices"] + 3*i + 2);
            }
        }

        for(int i = 0; i < quadTopol.getNumVertices(); i++){
            VecType coords;
            coords[0] = quadDeformedGeometry[3*i];
            coords[1] = quadDeformedGeometry[3*i+1];
            coords[2] = quadDeformedGeometry[3*i+2];
            coords[2] +=  2*(coords[0] - 0.5)*(coords[0] - 0.5);
            quadDeformedGeometry[3*i] = coords[0]; 
            quadDeformedGeometry[3*i + 1] = coords[1];
            quadDeformedGeometry[3*i + 2] = coords[2];
        }*/


        QuadMeshTopologySaver::setGeometry(mesh, quadDeformedGeometry);
        OpenMesh::IO::write_mesh(mesh, "quad_deformed.ply");

        VectorType factors_mem_bend(2);
        factors_mem_bend[0] = 10000.0;
        factors_mem_bend[1] = 1.0;

        VectorType factors_elasticity_dev(2);
        factors_elasticity_dev[0] = 1.0; // Elasticity factor
        factors_elasticity_dev[1] = 1.0; // Developability factor

        constraint_weights<DefaultConfigurator> weights;
        weights.fair_v = 0.0;//50.0;
        weights.fair_n = 0.0;//10.0;
        weights.fair_r = 0.0;
        weights.ruling_0 = 10.0;
        weights.ruling_1 = 10.0;
        weights.ruling_2 = 10.0;
        weights.dev = 1000.0;

        StripHandler<DefaultConfigurator> stripHandle(quadTopol);
        ConstraintIdx<DefaultConfigurator> constraintIdx(stripHandle, quadTopol);

        MatrixType weights_mat;
        weights.extend_weights(constraintIdx, weights_mat);

        Constraint<DefaultConfigurator> constraint(quadTopol, stripHandle, constraintIdx, varsIdx);
        ConstraintGrad<DefaultConfigurator> Dconstraint(quadTopol, stripHandle, constraintIdx, varsIdx);

        ElasticDevelopabilityCentroidEnergy<DefaultConfigurator> E_tot(quadTopol, triangleTopol, triBaseGeometry, factors_elasticity_dev,factors_mem_bend, constraint, varsIdx, weights_mat);
        ElasticDevelopabilityCentroidGradient<DefaultConfigurator> DE_tot(quadTopol, triangleTopol, triBaseGeometry, factors_elasticity_dev, factors_mem_bend, constraint, Dconstraint, varsIdx, weights_mat);
        //ElasticDevelopabilityCentroidHessian<DefaultConfigurator> D2E_tot(quadTopol, centroidTopol, quadReferenceGeometry, centroidReferenceGeometry, factors_elasticity_dev, factors_mem_bend, constraint, Dconstraint, weights_mat);

        // set outer optimization parameters
        OptimizationParameters<DefaultConfigurator> optPars;
        optPars.setGradientIterations( 200 );
        optPars.setBFGSIterations( 10000 );
        optPars.setQuietMode( SHOW_ALL );
        std::cout<<"Resets: "<<optPars.getBFGSReset()<<std::endl;

        VectorType init = VectorType::Zero(varsIdx["num_dofs"]);
        constraint.initialize_vars(quadDeformedGeometry, init);

        StripHandler<DefaultConfigurator> stripHandler(quadTopol);
        VarsIdx<DefaultConfigurator> variablesIdx(quadTopol);

        // Test:
        RealType test_energy;
        E_tot.apply(init, test_energy);

        VectorType constraint_val;
        constraint.apply(init, constraint_val);

        RealType result_2 = constraint_val.dot(weights_mat*constraint_val);
        std::cout<<"Test energy: "<<test_energy<<"; constraint value: "<<result_2<<std::endl;

        VectorType result = init;

        GradientDescent<DefaultConfigurator> GD( E_tot, DE_tot, optPars);
        GD.setBoundaryMask( bdryMask );
        GD.solve( init, result );

        QuasiNewtonBFGS<DefaultConfigurator> BFGS( E_tot, DE_tot, optPars);
        BFGS.setBoundaryMask( bdryMask );
        BFGS.solve( result, result );

        VectorType result_quad(3*quadTopol.getNumVertices());
        result_quad = result.segment(varsIdx["vertices"], 3*quadTopol.getNumVertices());

        QuadMeshTopologySaver::setGeometry(mesh, result_quad);
        OpenMesh::IO::write_mesh(mesh, "quad_optimized_new.ply");

        std::string normals_file = "/home/s24pjoha_hpc/goast_old_old/goast/examples/QuadMeshExamples/plotting/normals_new.txt";
        export_normals(quadTopol,result, variablesIdx, normals_file);
        std::string rw_rulings_file = "/home/s24pjoha_hpc/goast_old_old/goast/examples/QuadMeshExamples/plotting/rw_rulings_new.txt";
        export_rw_rulings(quadTopol,result, variablesIdx, rw_rulings_file);

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}