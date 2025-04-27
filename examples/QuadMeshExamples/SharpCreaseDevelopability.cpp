#include "Levenberg_Marquardt.h"
#include "Constraints.h"
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

    try{

        // ---------------------- TOPOLOGICAL PREPARATIONS --------------------------------------

        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "../../data/plate/basePlateSharpCreaseDevelopability_coarse.ply");
        QuadMeshTopologySaver quadTopol(mesh);
        VectorType geom;
        QuadMeshTopologySaver::getGeometry(mesh,geom);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        // --------------------- GENERATE BASE CENTROID TOPOLOGY ----------------------------

        TriMesh centroid_base_mesh = quadTopol.makeQuadMeshCentroid();
        MeshTopologySaver centroidTopol(centroid_base_mesh);
        VectorType centroid_base_geom;
        getGeometry(centroid_base_mesh, centroid_base_geom);

        // ---------------------  GENERATE MAP BETWEEN QUAD AND CENTROID TOPOLOGY ----------------------------
        std::map<int,int> quadToCentroidIdxMap;
        std::vector<int> noQuadVertexIndices;
        bool foundMatch = false;
        for(int i = 0; i < centroidTopol.getNumVertices(); i++){
            VecType coords_tri;
            getXYZCoord<VectorType, VecType>(centroid_base_geom,coords_tri, i);
            for(int j = 0; j < quadTopol.getNumVertices(); j++)
            {
                VecType coords_quad;
                coords_quad[0] = geom[3*j];
                coords_quad[1] = geom[3*j+1];
                coords_quad[2] = geom[3*j+2];

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

        // -------------------------- SETUP BOUNDARY ------------------------------------------

        double tolerance = 1e-4;

        std::vector<int> bdryMask;
        for(int i = 0; i < num_vertices; i++)
        {
            VecType coords;
            coords[0] = geom[3*i];
            coords[1] = geom[3*i+1];
            coords[2] = geom[3*i+2];
            if(std::abs(coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance || std::abs(coords[1]) < tolerance || std::abs(1.0 - coords[1]) < tolerance){
                bdryMask.push_back(i);
            }

            if(std::abs(coords[0]) < tolerance || std::abs(1-coords[0]) < tolerance)
            {
                if(coords[1] <= 0.5)
                {
                    geom[3*i + 2] = 0.4*coords[1];
                }
                else if(coords[1] > 0.5)
                {
                    geom[3*i + 2] = 0.4*(1.0 - coords[1]);
                }
                continue;
            }

            if(std::abs(coords[1]) < tolerance || std::abs(1-coords[1]) < tolerance)
            {
                if(coords[0] <= 0.5)
                {
                    geom[3*i + 2] = 0.4*coords[0];
                }
                else if(coords[0] > 0.5)
                {
                    geom[3*i + 2] = 0.4*(1.0 - coords[0]);
                }
                continue;
            }
        }

        // ----------------- GENERATE DEFORMED CENTROID MESHES --------------------------------

        QuadMeshTopologySaver::setGeometry(mesh, geom);
        QuadMeshTopologySaver new_quad_topol(mesh);
        // Now, we want to refine this into a triangle mesh
        TriMesh centroid_mesh = new_quad_topol.makeQuadMeshCentroid();
        VectorType centroid_geom;
        getGeometry(centroid_mesh, centroid_geom);
        MeshTopologySaver centroid_topol(centroid_mesh);
        OpenMesh::IO::write_mesh(centroid_mesh, "centroid.ply");
        // Now, first recalculate the boundary in the TriMesh format
        std::vector<int> bdryMaskTri;
        for(int i = 0; i < centroid_topol.getNumVertices(); i++)
        {
            VecType coords;
            getXYZCoord<VectorType, VecType>(centroid_base_geom,coords, i);

            if(std::abs(coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance || std::abs(coords[1]) < tolerance || std::abs(1.0 - coords[1]) < tolerance)
            {
                bdryMaskTri.push_back(i);
            }
        }
        extendBoundaryMask(centroid_topol.getNumVertices(), bdryMaskTri);

        // Now, apply the Dirichlet smoothing
        DirichletSmoother<DefaultConfigurator> dirichletSmoother(centroid_base_geom, bdryMaskTri, centroid_topol);
        VectorType smoothedGeom = VectorType::Zero(centroid_topol.getNumVertices() * 3);
        dirichletSmoother.apply(centroid_geom, smoothedGeom);

        setGeometry(centroid_mesh, smoothedGeom);
        OpenMesh::IO::write_mesh(centroid_mesh, "smoothed_centroid.ply");

        // ------------------- TRANSFER SMOOTHED CENTROID GEOMETRY TO QUAD MESH ---------------------

        VectorType smoothedQuadGeom = VectorType::Zero(quadTopol.getNumVertices() * 3);
        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            // get the smoothed centroid coordinates
            VecType centroid_coords;
            getXYZCoord<VectorType, VecType>(smoothedGeom,centroid_coords, quadToCentroidIdxMap[i]);
            smoothedQuadGeom[3*i] = centroid_coords[0];
            smoothedQuadGeom[3*i + 1] = centroid_coords[1];
            smoothedQuadGeom[3*i + 2] = centroid_coords[2];
        }

        QuadMeshTopologySaver::setGeometry(mesh, smoothedQuadGeom);
        OpenMesh::IO::write_mesh(mesh, "smoothed_quad.ply");

        QuadMeshTopologySaver::setGeometry(mesh, geom);

        VectorType bdryCoords(3*bdryMask.size());

        for(int i = 0; i < bdryMask.size(); i++)
        {
            VecType coords;
            coords[0] = geom[3*bdryMask[i]];
            coords[1] = geom[3*bdryMask[i]+1];
            coords[2] = geom[3*bdryMask[i]+2];

            bdryCoords[3*i] = coords[0];
            bdryCoords[3*i+1] = coords[1];
            bdryCoords[3*i+2] = coords[2];
        }

        QuadMeshTopologySaver::setGeometry(mesh, geom);
        OpenMesh::IO::write_mesh(mesh, "prepare_developability.ply");

        std::pair<std::vector<int>, VectorType> bdryData = std::make_pair(bdryMask,bdryCoords);

        StripHandler<DefaultConfigurator> stripHandle(quadTopol);

        ConstraintIdx<DefaultConfigurator> constraintIdx(stripHandle, quadTopol, bdryData);
        VarsIdx<DefaultConfigurator> variablesIdx(quadTopol);

        // --------------------- INITIALIZE CONSTRAINT WEIGHTS --------------------------------

          constraint_weights<DefaultConfigurator> weights;
          weights.fair_v = 50.0;
          weights.fair_n = 10.0;
          weights.fair_r = 0.0;
          weights.ruling_0 = 10.0;
          weights.ruling_1 = 10.0;
          weights.ruling_2 = 10.0;
          weights.bdry_opt = 100000.0;
          weights.dev = 1000.0;

        MatrixType weights_mat;
        weights.extend_weights(constraintIdx, weights_mat);

        //printMatrixToFile(weights_mat, "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/weights.txt");

        // -------------------------- SETUP Constraint and ConstraintGrad ----------------------------

        Constraint<DefaultConfigurator> constraint(quadTopol,stripHandle, constraintIdx, variablesIdx, bdryData);
        ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandle, constraintIdx, variablesIdx, bdryData);

        // -------------------------- INITIALIZE THE VARIABLES --------------------------------------------
        VectorType init = VectorType::Zero(variablesIdx["num_dofs"]);
        constraint.initialize_vars(smoothedQuadGeom, init);


        LevenbergMarquardtParams<DefaultConfigurator> pars;
        LMAlgorithm<DefaultConfigurator> lm(pars, constraint, constraintGrad, variablesIdx["num_dofs"], constraintIdx["num_cons"]);

        //std::string filepath_begin = "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/gauss_image_data_begin.txt";
        //plot_gauss_image(quadTopol,view, filepath_begin);

        VectorType Dest;
        lm.solve(init, Dest, weights_mat);
        auto test_normal = variablesIdx.face_normal(Dest, 0);
        VectorType DestGeom = Dest.segment(variablesIdx["vertices"],3*num_vertices);
        QuadMeshTopologySaver::setGeometry(mesh, DestGeom);
        OpenMesh::IO::write_mesh(mesh, "result_developability.ply");
        std::string normals_file = "/home/paul_johannssen/Desktop/masterarbeit/goast/examples/QuadMeshExamples/plotting/normals.txt";
        export_normals(quadTopol,Dest, variablesIdx, normals_file);
        std::string rw_rulings_file = "/home/paul_johannssen/Desktop/masterarbeit/goast/examples/QuadMeshExamples/plotting/rw_rulings.txt";
        export_rw_rulings(quadTopol,Dest, variablesIdx, rw_rulings_file);

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}