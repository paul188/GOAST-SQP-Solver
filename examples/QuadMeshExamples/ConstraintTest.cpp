#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <vector>

#include <goast/QuadMesh/QuadTopology.h>
#include <goast/Core.h>
#include "Constraints.h"

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

int main(int argc, char *argv[])
{

try{

    MyMesh mesh;
    
    OpenMesh::IO::read_mesh(mesh, "../../data/plate/quadMesh.ply");

    VectorType geom;

    QuadMeshTopologySaver quadTopol(mesh);

    size_t num_vertices = quadTopol.getNumVertices();

    QuadMeshTopologySaver::getGeometry(mesh,geom);

    VarsIdx<DefaultConfigurator> vars_idx(quadTopol);
    StripHandler<DefaultConfigurator> stripHandle(quadTopol);

    // set the boundary
    std::vector<int> bdryMask;
    VectorType bdryCoords;

    for(int i = 0; i < num_vertices; i++){
        VecType coords;
        coords[0] = geom[3*i];
        if(coords[0] == 0.0)
        {
            bdryMask.push_back(i);
        }
    }

    bdryCoords.resize(3*bdryMask.size());
    size_t counter = 0;
    for(int i = 0; i < bdryMask.size(); i++){
        VecType coords;
        coords[0] = geom[3*bdryMask[i]];
        coords[1] = geom[3*bdryMask[i]+1];
        coords[2] = geom[3*bdryMask[i]+2];
        if(coords[0] == 0.0)
        {
            bdryCoords[3*counter] = coords[0];
            bdryCoords[3*counter + 1] = coords[1];
            bdryCoords[3*counter + 2] = coords[2];
            counter++;
        }
    }

    std::pair<std::vector<int>, VectorType> bdryData = std::make_pair(bdryMask,bdryCoords);

    ConstraintIdx<DefaultConfigurator> cons_idx(stripHandle, quadTopol, bdryData);
    VarsIdx<DefaultConfigurator> variablesIdx(quadTopol);

    Constraint<DefaultConfigurator> constraint(quadTopol,stripHandle,cons_idx, vars_idx, bdryData);
    ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandle, cons_idx, vars_idx, bdryData);

    VectorType test_input = VectorType::Ones(vars_idx["num_dofs"]);

    test_input.segment(vars_idx["vertices"],3*num_vertices) = geom;
    
    VectorType test_input_2;
    constraint.initialize_vars(geom,test_input_2);

    VectorType dest;
    MatrixType Dest;
    constraint.apply(test_input_2,dest);
    constraintGrad.apply(test_input_2,Dest);

    VectorValuedDerivativeTester<DefaultConfigurator> tester(constraint,constraintGrad,0.01,Dest.rows());
    tester.plotAllDirections(test_input_2,"deriv_test/");

}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}
return 0;
}