#include "Levenberg_Marquardt.h"
#include "Constraints.h"
#include <goast/Core.h>
#include <random>
#include <iostream>
#include <utility>

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
        OpenMesh::IO::read_mesh(mesh, "../../data/plate/testDevelopabilityCylinder_fine.ply");
        QuadMeshTopologySaver quadTopol(mesh);
        VectorType geom;
        QuadMeshTopologySaver::getGeometry(mesh,geom);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        StripHandler<DefaultConfigurator> stripHandle(quadTopol);

        // -------------------------- SETUP BOUNDARY ------------------------------------------

        double tolerance = 1e-4;

        std::vector<int> bdryMask;
        for(int i = 0; i < num_vertices; i++)
        {
            VecType coords;
            coords[0] = geom[3*i];
            coords[1] = geom[3*i+1];
            coords[2] = geom[3*i+2];
            std::cout<<"Vertex this way: "<<i<<": "<<geom[3*i]<<"; "<<geom[3*i +2]<<std::endl;
            if(std::abs(1.0 + coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance){
                bdryMask.push_back(i);
            }
        }

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

        std::pair<std::vector<int>, VectorType> bdryData = std::make_pair(bdryMask,bdryCoords);

        ConstraintIdx<DefaultConfigurator> constraintIdx(stripHandle, quadTopol, bdryData);
        VarsIdx<DefaultConfigurator> variablesIdx(quadTopol);

        // --------------------- INITIALIZE CONSTRAINT WEIGHTS --------------------------------

          constraint_weights<DefaultConfigurator> weights;
          weights.fair_v = 0.5;//50.0;
          weights.fair_n = 0.5;//10.0;
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
        constraint.initialize_vars(geom, init);


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