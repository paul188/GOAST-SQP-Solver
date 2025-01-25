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

    auto faceStrips = quadTopol.getFaceStrips();
    auto halfedgeStrips = quadTopol.getHalfEdgeStrips();
    auto vertexStrips = quadTopol.getVertexStrips();

    std::cout<<"Number of face strips: "<<faceStrips.size()<<std::endl;
    std::cout<<"Number of halfedge strips: "<<halfedgeStrips.size()<<std::endl;
    std::cout<<"Number of vertex strips: "<<vertexStrips.size()<<std::endl;
    std::cout<<"unicorn"<<std::endl;

    QuadMeshTopologySaver::getGeometry(mesh,geom);

    VectorView<DefaultConfigurator> view(quadTopol);

    StripHandler<DefaultConfigurator> stripHandle(quadTopol);
    std::cout<<"Test1"<<std::endl;
    std::cout<<stripHandle._num_vertexStrips<<std::endl;
    std::cout<<stripHandle._num_halfedgeStrips<<std::endl;
    std::cout<<stripHandle._num_faceStrips<<std::endl;
    std::cout<<"Test2"<<std::endl;

    Constraint<DefaultConfigurator> constraint(quadTopol,stripHandle);
    ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandle);

    VectorType test_input = VectorType::Ones(view._idx["num_dofs"]);

    test_input.segment(view._idx["vertices"],3*num_vertices) = geom;

    //printVectorToFile(test_input,"./deriv_test/test_input.txt");
    VectorType dest;
    MatrixType Dest;
    constraint.apply(test_input,dest);
    constraintGrad.apply(test_input,Dest);

    VectorValuedDerivativeTester<DefaultConfigurator> tester(constraint,constraintGrad,0.01,Dest.rows());
    tester.plotAllDirections(test_input,"deriv_test_2/");

}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}
return 0;
}