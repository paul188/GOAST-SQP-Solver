#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <vector>

#include <goast/QuadMesh/QuadTopology.h>
#include <goast/Core.h>
#include <goast/Developability/Developability.h>

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

int main(int argc, char *argv[])
{

try{

    MyMesh mesh;
    
    OpenMesh::IO::read_mesh(mesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/quadMesh2_coarse_coarse.ply");

    VectorType geom;

    QuadMeshTopologySaver quadTopol(mesh);

    size_t num_vertices = quadTopol.getNumVertices();

    std::cout<<"Number of vertices: "<<num_vertices<<std::endl;

    QuadMeshTopologySaver::getGeometry(mesh,geom);

    VarsIdx<DefaultConfigurator> vars_idx(quadTopol);
    StripHandler<DefaultConfigurator> stripHandle(quadTopol);

    // set the boundary
    std::vector<int> bdryMask;

    for(int i = 0; i < num_vertices; i++){
        VecType coords;
        coords[0] = geom[3*i];
        if(coords[0] == 0.0)
        {
            bdryMask.push_back(i);
        }
    }

    ConstraintIdx<DefaultConfigurator> cons_idx(stripHandle, quadTopol);
    VarsIdx<DefaultConfigurator> variablesIdx(quadTopol);

    Constraint<DefaultConfigurator> constraint(quadTopol,stripHandle,cons_idx, vars_idx);
    ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandle, cons_idx, vars_idx);
    ConstraintHessian<DefaultConfigurator> constraintHess(quadTopol,stripHandle, cons_idx, vars_idx);
    BoundaryDOFS_quad<DefaultConfigurator> bdryDOFs(bdryMask, num_vertices, vars_idx);

    VectorType test_input = VectorType::Ones(vars_idx["num_dofs"]);

    test_input.segment(vars_idx["vertices"],3*num_vertices) = geom;
    
    VectorType test_input_2;
    constraint.initialize_vars(geom,test_input_2);

    VectorType dest;
    MatrixType Dest;
    GenericTensor<DefaultConfigurator::SparseMatrixType> genericTensor;

    constraint.apply(test_input_2,dest);
    constraintGrad.apply(test_input_2,Dest);
    constraintHess.apply(test_input_2,genericTensor);

    for(int i = 0; i < variablesIdx["num_dofs"]; i++)
    {
        for(int j = 0; j < variablesIdx["num_dofs"]; j++)
        {
            if((genericTensor[1]).coeffRef(i,j) != 0.0)
            {
                std::cout<<"Hessian["<<i<<"]["<<j<<"] = "<<(genericTensor[1]).coeffRef(i,j)<<std::endl;
            }
            //std::cout<<"Hessian["<<i<<"]["<<j<<"] = "<<(genericTensor[1]).coeffRef(i,j)<<std::endl;
        }
    }

    printMatrixToFile(Dest.toDense(), "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/deriv_test/constraintGrad.txt");

    std::cout<<"Test the variables idx: "<<std::endl;
    std::cout<<"vertices: "<<variablesIdx["vertices"]<<std::endl;
    std::cout<<"num_dofs: "<<variablesIdx["num_dofs"]<<std::endl;
    std::cout<<"dummy weight: "<<variablesIdx["dummy_weights"]<<std::endl;
    std::cout<<"reweighted vertices: "<<variablesIdx["reweighted_vertices"]<<std::endl;


    //VectorValuedDerivativeTester<DefaultConfigurator> tester(constraint,constraintGrad,0.01,Dest.rows());
    //tester.plotAllDirections(test_input_2,"deriv_test/");

    // Test the hessian of every single constraint
    // first vertex_1
    TensorValuedDerivativeTester<DefaultConfigurator> tester(constraintGrad, constraintHess,0.01, 3*quadTopol.getNumVertices());
    tester.plotAllDirections(test_input_2,"/home/s24pjoha_hpc/goast_old_old/goast/build/examples/deriv_test/");

}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}
return 0;
}