#include <goast/Developability/Developability.h>
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
        OpenMesh::IO::read_mesh(mesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/new.ply");
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
            bdryMask.push_back(i);
            //std::cout<<"Vertex this way: "<<i<<": "<<geom[3*i]<<"; "<<geom[3*i +2]<<std::endl;
            //if(std::abs(1.0 + coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance){
                //bdryMask.push_back(i);
            //}
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

        ConstraintIdx<DefaultConfigurator> constraintIdx(stripHandle, quadTopol);
        VarsIdx<DefaultConfigurator> variablesIdx(quadTopol);

        // --------------------- INITIALIZE CONSTRAINT WEIGHTS --------------------------------

          constraint_weights<DefaultConfigurator> weights;
          //weights.fair_v = 50.0;
          //weights.fair_n = 10.0;
          weights.fair_r = 0.0;
          weights.ruling_0 = 10.0;
          weights.ruling_1 = 10.0;
          weights.ruling_2 = 10.0;
          weights.dev = 1000.0;

        MatrixType weights_mat;
        weights.extend_weights(constraintIdx, weights_mat);

        //printMatrixToFile(weights_mat, "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/weights.txt");

        // -------------------------- SETUP Constraint and ConstraintGrad ----------------------------

        Constraint<DefaultConfigurator> constraint(quadTopol,stripHandle, constraintIdx, variablesIdx);
        ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandle, constraintIdx, variablesIdx);
        BoundaryDOFS_quad<DefaultConfigurator> bdryDOFs(bdryMask, num_vertices, variablesIdx);

        // -------------------------- INITIALIZE THE VARIABLES --------------------------------------------
        VectorType init = VectorType::Zero(variablesIdx["num_dofs"]);
        init.segment(variablesIdx["vertices"], 3 * num_vertices) = geom;
        //Initialize all the normals with vectors in the pos z-direction
        for(int i = 0; i < num_faces; i++)
        {
            init[variablesIdx["normals"] + 3*i] = 0.0;
            init[variablesIdx["normals"] + 3*i +1] = 0.0;
            init[variablesIdx["normals"] + 3*i +2] = 1.0;
        }

        LevenbergMarquardtParams<DefaultConfigurator> pars;
        LMAlgorithm<DefaultConfigurator> lm(pars, constraint, constraintGrad, bdryDOFs, variablesIdx["num_dofs"], constraintIdx["num_cons"]);

        //std::string filepath_begin = "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/gauss_image_data_begin.txt";
        //plot_gauss_image(quadTopol,view, filepath_begin);

        VectorType Dest;
        lm.solve(init, Dest, weights_mat);
        auto test_normal = variablesIdx.face_normal(Dest, 0);
        VectorType DestGeom = Dest.segment(variablesIdx["vertices"],3*num_vertices);
        QuadMeshTopologySaver::setGeometry(mesh, DestGeom);
        OpenMesh::IO::write_mesh(mesh, "result_developability.ply");

        std::cout<<"Geometry values: "<<std::endl;
        for(int i = 0; i < num_vertices; i++)
        {
            std::cout<<DestGeom[3*i]<<"; "<<DestGeom[3*i+1]<<"; "<<DestGeom[3*i+2]<<std::endl;
        }
        //std::string normals_file = "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/normals.txt";
        std::string normals_file = "/home/s24pjoha_hpc/goast_old_old/goast/examples/QuadMeshExamples/plotting/normals.txt";
        export_normals(quadTopol,Dest, variablesIdx, normals_file);
        std::string rw_rulings_file = "/home/s24pjoha_hpc/goast_old_old/goast/examples/QuadMeshExamples/plotting/rw_rulings.txt";
        //std::string rw_rulings_file = "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/rw_rulings.txt";
        export_rw_rulings(quadTopol,Dest, variablesIdx, rw_rulings_file);

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}