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

        // ---------------------  GENERATE MAPS BETWEEN QUAD AND CENTROID TOPOLOGIES ----------------------------
        std::map<int,int> quadToCentroidIdxMap;
        std::vector<int> noQuadVertexIndices;
        bool foundMatch = false;
        for(int i = 0; i < centroidTopol.getNumVertices(); i++){
            VecType coords_tri;
            getXYZCoord<VectorType, VecType>(centroidReferenceGeometry,coords_tri, i);
            for(int j = 0; j < quadTopol.getNumVertices(); j++)
            {
                VecType coords_quad;
                coords_quad[0] = quadReferenceGeometry[3*j];
                coords_quad[1] = quadReferenceGeometry[3*j+1];
                coords_quad[2] = quadReferenceGeometry[3*j+2];

                if((coords_quad - coords_tri).norm() < 1e-4)
                {
                    quadToCentroidIdxMap[j] = i;
                    foundMatch = true;
                    break;
                }
                
            }
            if(!foundMatch)
            {
                noQuadVertexIndices.push_back(i);
            }
            
            // reset foundMatch
            foundMatch = false;
        }

        std::map<int,int> centroidToQuadIdxMap;
        int quadTopol_counter = 0;
        for(int i = 0; i < centroidTopol.getNumVertices(); i++)
        {
            centroidToQuadIdxMap[i] = -1;
        }

        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            centroidToQuadIdxMap[quadToCentroidIdxMap[i]] = i;
        }

        FullMatrixType neighbouringnoQuadVertices(noQuadVertexIndices.size(), 4);
        for(int i = 0; i < noQuadVertexIndices.size(); i++)
        {
            int index = noQuadVertexIndices[i];
            TriMesh::VertexHandle vh = centroid_base_mesh.vertex_handle(index);
            int counter = 0;
            for (TriMesh::VertexVertexIter vv_it = centroid_base_mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
                neighbouringnoQuadVertices(i, counter ) = vv_it->idx();
                counter++;
                if(counter >= 4) break; // only take the first four neighbours
            }

        }
        // -------------------------- SETUP BOUNDARY ------------------------------------------

        VarsIdx<DefaultConfigurator> varsIdx(quadTopol);

        double tolerance = 1e-4;

        VectorType centroidDeformedGeometry = centroidReferenceGeometry;

        std::vector<int> bdryMask;
        for(int i = 0; i < centroidTopol.getNumVertices(); i++)
        {
            VecType coords;
            getXYZCoord<VectorType, VecType>(centroidReferenceGeometry,coords, i);
            if(std::abs(coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance){
                bdryMask.push_back(varsIdx["vertices"] + i);
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

        extendBoundaryMask(centroidTopol.getNumVertices(), bdryMask);

        VectorType quadDeformedGeometry = quadReferenceGeometry;

        // Now, create also the deformed quad geometry
        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            VecType coords;
            coords[0] = quadReferenceGeometry[3*i];
            coords[1] = quadReferenceGeometry[3*i+1];
            coords[2] = quadReferenceGeometry[3*i+2];

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

        QuadMeshTopologySaver::setGeometry(mesh, quadDeformedGeometry);
        OpenMesh::IO::write_mesh(mesh, "quad_deformed.ply");

        // First, output the centroid reference geometry
        setGeometry(centroid_base_mesh, centroidReferenceGeometry);
        OpenMesh::IO::write_mesh(centroid_base_mesh, "centroid_base.ply");
        // Now, output the centroid deformed geometry
        setGeometry(centroid_base_mesh, centroidDeformedGeometry);
        OpenMesh::IO::write_mesh(centroid_base_mesh, "centroid_deformed.ply");

        // get the centroid deformed_2
        //OpenMesh::IO::read_mesh(centroid_base_mesh, "/home/paul_johannssen/Desktop/masterarbeit/goast/goast/build/examples/centroid_deformed_2.ply");
        //setGeometry(centroid_base_mesh, centroidDeformedGeometry);

        ConstraintSqrdReduced<DefaultConfigurator> constraintSqrdReduced(quadTopol);
        RealType squared_dev_energy;
        constraintSqrdReduced.apply(quadDeformedGeometry, squared_dev_energy);
        std::cout<<"Squared developability energy: "<<squared_dev_energy<<std::endl;

        VectorType squared_dev_gradient;
        ConstraintSqrdReducedGradient<DefaultConfigurator> constraintSqrdReducedGradient(quadTopol);
        constraintSqrdReducedGradient.apply(quadDeformedGeometry, squared_dev_gradient);
        std::cout<<"Squared developability gradient: "<<squared_dev_gradient.norm()<<std::endl;

        VectorType factors_mem_bend(2);
        factors_mem_bend[0] = 10000.0;
        factors_mem_bend[1] = 1.0;

        VectorType factors_elasticity_dev(2);
        factors_elasticity_dev[0] = 1.0; // Elasticity factor
        factors_elasticity_dev[1] = 10000.0; // Developability factor

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

        std::cout<<"Centroid topology num vertices before: "<<centroidTopol.getNumVertices()<<std::endl;
        ElasticDevelopabilityCentroidEnergy<DefaultConfigurator> E_tot(quadTopol, centroidTopol, quadReferenceGeometry, centroidReferenceGeometry, factors_elasticity_dev,factors_mem_bend, constraint, varsIdx, weights_mat);
        ElasticDevelopabilityCentroidGradient<DefaultConfigurator> DE_tot(quadTopol, centroidTopol, quadReferenceGeometry, centroidReferenceGeometry, factors_elasticity_dev, factors_mem_bend, constraint, Dconstraint, varsIdx, weights_mat);
        //ElasticDevelopabilityCentroidHessian<DefaultConfigurator> D2E_tot(quadTopol, centroidTopol, quadReferenceGeometry, centroidReferenceGeometry, factors_elasticity_dev, factors_mem_bend, constraint, Dconstraint, weights_mat);

        // set outer optimization parameters
        OptimizationParameters<DefaultConfigurator> optPars;
        optPars.setGradientIterations( 200 );
        optPars.setBFGSIterations( 2000 );
        optPars.setQuietMode( SHOW_ALL );

        VectorType init = VectorType::Zero(varsIdx["num_dofs"]);
        constraint.initialize_vars(quadDeformedGeometry, init);

        StripHandler<DefaultConfigurator> stripHandler(quadTopol);
        VarsIdx<DefaultConfigurator> variablesIdx(quadTopol);

        // Will not use the quad initialization directly, but instead place the centroid reference geometry in the first dofs
        VectorType init_centroid(init.size() + 3*centroidTopol.getNumVertices() - 3*quadTopol.getNumVertices());
        init_centroid.segment(0, variablesIdx["vertices"]) = init.segment(0, variablesIdx["vertices"]);
        init_centroid.segment(variablesIdx["vertices"], 3*centroidTopol.getNumVertices()) = centroidDeformedGeometry;
        init_centroid.segment(variablesIdx["vertices"] + 3*centroidTopol.getNumVertices(), init_centroid.size() - (variablesIdx["vertices"] + 3*centroidTopol.getNumVertices())) = init.segment(variablesIdx["vertices"] + 3*quadTopol.getNumVertices(), init.size() - (variablesIdx["vertices"] + 3*quadTopol.getNumVertices()));

        // Test:
        RealType test_energy;
        E_tot.apply(init_centroid, test_energy);

        VectorType constraint_val;
        constraint.apply(init, constraint_val);

        RealType result_2 = constraint_val.dot(weights_mat*constraint_val);
        std::cout<<"Test energy: "<<test_energy<<"; constraint value: "<<result_2<<std::endl;

        VectorType result_centroid(init_centroid.size());

        GradientDescent<DefaultConfigurator> GD( E_tot, DE_tot, optPars);
        GD.setBoundaryMask( bdryMask );
        GD.solve( init_centroid, result_centroid );

        QuasiNewtonBFGS<DefaultConfigurator> BFGS( E_tot, DE_tot, optPars);
        BFGS.setBoundaryMask( bdryMask );
        BFGS.solve( result_centroid, result_centroid );
        
        VectorType result_centroid_geometry = result_centroid.segment(variablesIdx["vertices"], 3*centroidTopol.getNumVertices());

        // obtain the old geometry and the quad variables
        VectorType quadVariables(variablesIdx["num_dofs"]);
        CoordConverter<DefaultConfigurator> coordConverter(quadTopol, centroidTopol, quadReferenceGeometry, centroidReferenceGeometry);

        quadVariables.segment(0, variablesIdx["vertices"]) = result_centroid.segment(0, variablesIdx["vertices"]);
        VectorType quadGeom_result;
        coordConverter.convertCentroidToQuad(result_centroid_geometry, quadGeom_result);
        quadVariables.segment(variablesIdx["vertices"], 3*quadTopol.getNumVertices()) = quadGeom_result;
        quadVariables.segment(variablesIdx["vertices"] + 3*quadTopol.getNumVertices(), quadVariables.size() - (variablesIdx["vertices"] + 3*quadTopol.getNumVertices())) = result_centroid.segment(variablesIdx["vertices"] + 3*centroidTopol.getNumVertices(), result_centroid.size() - (variablesIdx["vertices"] + 3*centroidTopol.getNumVertices()));

        setGeometry(centroid_base_mesh, result_centroid_geometry);
        OpenMesh::IO::write_mesh(centroid_base_mesh, "centroid_result.ply");

        QuadMeshTopologySaver::setGeometry( mesh, quadGeom_result);
        OpenMesh::IO::write_mesh(mesh, "quad_result.ply");

        std::string normals_file = "/home/s24pjoha_hpc/goast_old_old/goast/examples/QuadMeshExamples/plotting/normals.txt";
        export_normals(quadTopol,quadVariables, variablesIdx, normals_file);
        std::string rw_rulings_file = "/home/s24pjoha_hpc/goast_old_old/goast/examples/QuadMeshExamples/plotting/rw_rulings.txt";
        export_rw_rulings(quadTopol,quadVariables, variablesIdx, rw_rulings_file);

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}