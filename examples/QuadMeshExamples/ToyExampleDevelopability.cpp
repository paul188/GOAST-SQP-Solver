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

    try{

        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "../../data/plate/testPlateDevelopability.ply");

        VectorType geom;

        QuadMeshTopologySaver quadTopol(mesh);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        QuadMeshTopologySaver::getGeometry(mesh,geom);

        LevenbergMarquardtParams<DefaultConfigurator> pars;
        StripHandler<DefaultConfigurator> stripHandle(quadTopol);
        Constraint<DefaultConfigurator> constraint(quadTopol,stripHandle);
        ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandle);
        VectorView<DefaultConfigurator> view(quadTopol);
        ConstraintView<DefaultConfigurator> cons_view(stripHandle, quadTopol);

        LMAlgorithm<DefaultConfigurator> lm(pars, constraint, constraintGrad, view._idx["num_dofs"], cons_view._cons_idx["num_cons"]);
        VectorType init = VectorType::Zero(view._idx["num_dofs"]);
        init.segment(view._idx["vertices"],3*num_vertices) = geom;
        init.segment(view._idx["weights"],num_vertices) = VectorType::Ones(num_vertices);
        init.segment(view._idx["normals"], 3*num_faces) = VectorType::Ones(3*num_faces);
        VectorType Dest;
        lm.solve(init, Dest);
        VectorType DestGeom = Dest.segment(view._idx["vertices"],3*num_vertices);
        QuadMeshTopologySaver::setGeometry(mesh, DestGeom);
        OpenMesh::IO::write_mesh(mesh, "result_developability.ply");

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}