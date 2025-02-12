#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <vector>

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

    for(const auto & vertex: mesh.vertices()){
        VertexHandle vh = mesh.vertex_handle(vertex.idx());
        MyMesh::Point p = mesh.point(vh);
        std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
    }

    VectorType geom;

    QuadMeshTopologySaver quadTopol(mesh);

    size_t num_vertices = quadTopol.getNumVertices();

    std::cout<<"Number of vertices: "<<num_vertices<<std::endl;

    QuadMeshTopologySaver::getGeometry(mesh,geom);

    std::cout<<"Test2"<<std::endl;

    for(int i = 0; i < geom.size(); i++)
    {
        std::cout<<geom[i]<<std::endl;
    }

    std::cout<<"Test2 end"<<std::endl;

    VectorView<DefaultConfigurator> view(quadTopol);
    StripHandler<DefaultConfigurator> stripHandle(quadTopol);

    // set the boundary
    std::vector<int> bdryMask;
    VectorType bdryCoords(3*num_vertices);

    for(int i = 0; i < num_vertices; i++){
        VecType coords;
        coords[0] = geom[3*i];
        coords[1] = geom[3*i+1];
        coords[2] = geom[3*i+2];
        std::cout<<"Coordinates: "<<geom[3*i]<<" "<<geom[3*i+1]<<" "<<geom[3*i+2]<<std::endl;
        if(coords[0] == 0.0)
        {
            std::cout<<"Boundary vertex found"<<std::endl;
            bdryCoords[3*bdryMask.size()] = coords[0];
            bdryCoords[3*bdryMask.size() +1] = coords[1];
            bdryCoords[3*bdryMask.size() + 2] = coords[2];
            bdryMask.push_back(i);
        }
    }

    bdryCoords.conservativeResize(3*bdryMask.size());

    std::pair<std::vector<int>, VectorType> bdryData = std::make_pair(bdryMask,bdryCoords);

    Constraint<DefaultConfigurator> constraint(quadTopol,stripHandle, bdryData);
    ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandle, bdryData);

    VectorType test_input = VectorType::Ones(view._idx["num_dofs"]);

    test_input.segment(view._idx["vertices"],3*num_vertices) = geom;

    //printVectorToFile(test_input,"./deriv_test/test_input.txt");
    VectorType dest;
    MatrixType Dest;
    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    constraint.apply(test_input,dest);
    constraintGrad.apply(test_input,Dest);
    //printSparseMatrix(Dest, "boundary_constraint_matrix.txt");
    //printVectorToFile(dest,"boundary_constraint_vector.txt");

    VectorValuedDerivativeTester<DefaultConfigurator> tester(constraint,constraintGrad,0.01,Dest.rows());
    tester.plotAllDirections(test_input,"deriv_test_2/");

}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}
return 0;
}