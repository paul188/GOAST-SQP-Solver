#include "Levenberg_Marquardt.h"
#include "Constraints.h"
#include <goast/Core.h>
#include <random>
#include <iostream>
#include <utility>

int main()
{
    using VectorType = typename DefaultConfigurator::VectorType;
    using MatrixType = typename DefaultConfigurator::SparseMatrixType;
    using VecType = typename DefaultConfigurator::VecType;

    try{

        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "../../data/plate/testDevelopabilityBdry.ply");

        VectorType geom;

        QuadMeshTopologySaver quadTopol(mesh);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        QuadMeshTopologySaver::getGeometry(mesh,geom);

        std::vector<int> bdryMask;
        VectorType bdryCoords;

        // setting the boundary indices
        for(int i = 0; i < num_vertices; i++){
            VecType coords;
            coords[0] = geom[3*i];
            coords[1] = geom[3*i+1];
            coords[2] = geom[3*i+2];
            if(std::abs(coords[0] - 1.0) < 1e-1)
            {
                bdryMask.push_back(i);
            }
        }

        // setting the boundary coordinates
        bdryCoords.resize(3*bdryMask.size());

        for(int i = 0; i < bdryMask.size(); i++){
            VecType coords;
            coords[0] = geom[3*i];
            coords[1] = geom[3*i+1];
            coords[2] = geom[3*i+2];
            if(std::abs(coords[0] - 1.0) < 1e-1)
            {
                bdryCoords[3*i] = coords[0];
                bdryCoords[3*i +1] = coords[1];
                bdryCoords[3*i + 2] = coords[2];
            }
        }

        std::pair<std::vector<int>, VectorType> bdryData = std::make_pair(bdryMask,bdryCoords);

        constraint_weights<DefaultConfigurator> weights;
        weights.fair_v = 0.0;
        weights.fair_n = 0.0;
        weights.fair_r = 0.0;
        weights.ruling_0 = 10.0;
        weights.ruling_1 = 10.0;
        weights.ruling_2 = 10.0;
        weights.bdry_opt = 100000.0;
        weights.dev = 1000.0;

        LevenbergMarquardtParams<DefaultConfigurator> pars;
        StripHandler<DefaultConfigurator> stripHandle(quadTopol);
        Constraint<DefaultConfigurator> constraint(quadTopol,stripHandle, bdryData);
        ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandle, bdryData);
        VectorView<DefaultConfigurator> view(quadTopol);
        ConstraintView<DefaultConfigurator> cons_view(stripHandle, quadTopol, bdryData);

        LMAlgorithm<DefaultConfigurator> lm(pars, constraint, constraintGrad, view._idx["num_dofs"], cons_view._cons_idx["num_cons"]);
        
        // CHOOSE GOOD INITIAL GUESS
        VectorType init = VectorType::Zero(view._idx["num_dofs"]);
        init.segment(view._idx["vertices"],3*num_vertices) = geom;
        init.segment(view._idx["reweighted_vertices"],3*num_vertices) = geom;
        init.segment(view._idx["weights"],num_vertices) = VectorType::Ones(num_vertices);
        
        constraint.initialize_vars(geom, init);
        view.set_vector(init);

        VectorType test_dest;
        constraint.apply(init, test_dest);

        std::cout<<"Test the first constraint vector: "<<std::endl;
        for(int i = 0; i < test_dest.size(); i++){
            std::cout<<"Index: "<<i<<" Test dest: "<<test_dest[i]<<std::endl;
        }

        VectorType Dest;
        MatrixType Weights;
        cons_view.extend_weights(weights, Weights);
        lm.solve(init, Dest, Weights);
        VectorType DestGeom = Dest.segment(view._idx["vertices"],3*num_vertices);
        QuadMeshTopologySaver::setGeometry(mesh, DestGeom);
        OpenMesh::IO::write_mesh(mesh, "result_developability.ply");

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}