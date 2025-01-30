#include <iostream>
#include <chrono>
#include <ctime>
#include <string>

#include <goast/Quads.h>
#include <goast/Core.h>
#include "Constraints.h"

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

/** 
 * \brief Simulation of a simple fold 
 * \author Johannssen
 *
 * We optimize this energy by direct optimization via gradient descent or BFGS.
 * 
 */

/**/

typedef DefaultConfigurator::SparseMatrixType MatrixType;
typedef DefaultConfigurator::FullMatrixType FullMatrixType;
typedef DefaultConfigurator::VectorType VectorType;

int main(int argc, char *argv[])
{

try{

    MyMesh mesh;
    OpenMesh::IO::read_mesh(mesh, "../../data/plate/quadMesh.ply");

    VectorType geom;

    QuadMeshTopologySaver quadTopol(mesh);

    size_t num_vertices = quadTopol.getNumVertices();

    QuadMeshTopologySaver::getGeometry(mesh,geom);

    VectorView<DefaultConfigurator> view(quadTopol);

    StripHandler<DefaultConfigurator> stripHandle(quadTopol);

    Constraint<DefaultConfigurator> constraint(quadTopol,stripHandle);
    ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandle);

    VectorType test_input = VectorType::Ones(view._idx["num_dofs"]);

    test_input.segment(view._idx["vertices"],3*num_vertices) = geom;

    //printVectorToFile(test_input,"./deriv_test/test_input.txt");
    VectorType dest;
    MatrixType Dest;
    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    constraint.apply(test_input,dest);
    //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    //begin = std::chrono::steady_clock::now();
    constraintGrad.apply(test_input,Dest);
    //end = std::chrono::steady_clock::now();
    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

    VectorValuedDerivativeTester<DefaultConfigurator> tester(constraint,constraintGrad,0.01,Dest.rows());
    tester.plotAllDirections(test_input,"deriv_test_2/");

}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}
return 0;
}