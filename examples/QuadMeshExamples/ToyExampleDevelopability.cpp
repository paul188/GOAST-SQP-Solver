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

        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "../../data/plate/testPlateDevelopability.ply");

        VectorType geom;

        QuadMeshTopologySaver quadTopol(mesh);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        QuadMeshTopologySaver::getGeometry(mesh,geom);

        constraint_weights<DefaultConfigurator> weights;
        weights.fair_v = 50.0;
        weights.fair_n = 10.0;
        weights.fair_r = 0.0;
        weights.ruling_0 = 10.0;
        weights.ruling_1 = 10.0;
        weights.ruling_2 = 10.0;
        weights.bdry_opt = 100000.0;
        weights.dev = 1000.0;

        VectorType bdryCoords_set(3 * 66);  
        bdryCoords_set <<  
            1.06026, -1.06233 , 0,  
            1.06448, -1.00201, 0.02039,  
            1.06966, -0.94005, 0.040216,  
            1.07484, -0.877396, 0.059336,  
            1.07998, -0.814016, 0.077644,  
            1.08502, -0.749882, 0.095033,  
            1.0899, -0.684976, 0.111397,  
            1.09458, -0.619297, 0.126628,  
            1.099, -0.552853, 0.140623,  
            1.10309, -0.48567, 0.153278,  
            1.10679, -0.41779, 0.164498,  
            1.11005, -0.34927, 0.17419,  
            1.11282, -0.280185, 0.182274,  
            1.11504, -0.210623, 0.188678,  
            1.11667, -0.140687, 0.193344,  
            1.1177, -0.070492, 0.196229,  
            1.11809, -0.000161, 0.197305,  
            1.11784, 0.070181, 0.196562,  
            1.11696, 0.140406, 0.194007,  
            1.11545, 0.210391, 0.189664,  
            1.11336, 0.28002, 0.183573,  
            1.11071, 0.349189, 0.175789,  
            1.10756, 0.417804, 0.166382,  
            1.10395, 0.48579, 0.155432,  
            1.09995, 0.553086, 0.143025,  
            1.09561, 0.619647, 0.12926,  
            1.09099, 0.685446, 0.114237,  
            1.08614, 0.750469, 0.098061,  
            1.08113, 0.814717, 0.080837,  
            1.07602, 0.878205, 0.062672,  
            1.07084, 0.940959, 0.043675,  
            1.06566, 1.00301, 0.023948,  
            1.06145, 1.06338, 0.003637,
            -1.06417, -1.06235, -0.00527, // Now, here startsthe other side of the boundary
            -1.06864, -1.00203, 0.01505,
            -1.07413, -0.940068, 0.034784,
            -1.07962, -0.877408, 0.053812,
            -1.08505, -0.81402, 0.07203,
            -1.09037, -0.749877, 0.089333,
            -1.09553, -0.684962, 0.105614,
            -1.10047, -0.619274, 0.120768,
            -1.10512, -0.552823, 0.13469,
            -1.10943, -0.485636, 0.147281,
            -1.11332, -0.417754, 0.158442,
            -1.11675, -0.349237, 0.168086,
            -1.11965, -0.280158, 0.176132,
            -1.12198, -0.210606, 0.182508,
            -1.12371, -0.140683, 0.187157,
            -1.12478, -0.070503, 0.190037,
            -1.1252, -0.000187, 0.191121,
            -1.12494, 0.07014, 0.190397,
            -1.12402, 0.140354, 0.187873,
            -1.12245, 0.21033, 0.183572,
            -1.12026, 0.279955, 0.177534,
            -1.11748, 0.349125, 0.169812,
            -1.11418, 0.417746, 0.160475,
            -1.11039, 0.485742, 0.149602,
            -1.10618, 0.553052, 0.13728,
            -1.10161, 0.61963, 0.123604,
            -1.09673, 0.685448, 0.108675,
            -1.09162, 0.75049, 0.092595,
            -1.08633, 0.814757, 0.075471,
            -1.08092, 0.878262, 0.057408,
            -1.07544, 0.941031, 0.038512,
            -1.06995, 1.00309, 0.018889,
            -1.06548, 1.06346, -0.001344;


        std::vector<int> bdryMask;
        for(int i = 0; i < num_vertices; i++)
        {
            VecType coords;
            coords[0] = geom[3*i];
            coords[1] = geom[3*i+1];
            coords[2] = geom[3*i+2];
            for(int j = 0; j < (bdryCoords_set.size() / 3); j++){
                VecType bdry_coords;
                bdry_coords[0] = bdryCoords_set[3*j];
                bdry_coords[1] = bdryCoords_set[3*j+1];
                bdry_coords[2] = bdryCoords_set[3*j+2];
                if((coords - bdry_coords).norm() < 0.0005){
                    bdryMask.push_back(i);
                }
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

        LevenbergMarquardtParams<DefaultConfigurator> pars;
        StripHandler<DefaultConfigurator> stripHandle(quadTopol);
        VarsIdx<DefaultConfigurator> vars_idx(quadTopol);
        ConstraintIdx<DefaultConfigurator> cons_idx(stripHandle, quadTopol, bdryData);
        Constraint<DefaultConfigurator> constraint(quadTopol,stripHandle, cons_idx, vars_idx, bdryData);
        ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandle,cons_idx, vars_idx, bdryData);
        
        LMAlgorithm<DefaultConfigurator> lm(pars, constraint, constraintGrad, vars_idx["num_dofs"], cons_idx["num_cons"]);

        // CHOOSE GOOD INITIAL GUESS
        VectorType init = VectorType::Zero(vars_idx["num_dofs"]);
        constraint.initialize_vars(geom, init);

        //::string filepath_begin = "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/gauss_image_data_begin.txt";
        //plot_gauss_image(quadTopol,view, filepath_begin);

        VectorType Dest;
        MatrixType Weights;
        weights.extend_weights(cons_idx, Weights);
        lm.solve(init, Dest, Weights);
        VectorType DestGeom = Dest.segment(vars_idx["vertices"],3*num_vertices);
        QuadMeshTopologySaver::setGeometry(mesh, DestGeom);
        OpenMesh::IO::write_mesh(mesh, "result_developability.ply");
        /*
        std::string filepath_end = "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/gauss_image_data_final.txt";
        plot_gauss_image(quadTopol,constraint, filepath_end);
        std::string normals_file = "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/normals.txt";
        export_normals(quadTopol,constraint, normals_file);
        std::string rw_rulings_file = "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/rw_rulings.txt";
        export_rw_rulings(quadTopol,constraint, rw_rulings_file);
        std::string rulings_file = "/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/rulings.txt";
        export_rulings(quadTopol,constraint, rulings_file);
        */

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}