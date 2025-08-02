#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <vector>

#include <goast/QuadMesh/QuadTopology.h>
#include <goast/Core.h>
#include <goast/Developability/Developability.h>
#include <random>

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

int main(int argc, char *argv[])
{

try{

    MyMesh mesh;
    
    OpenMesh::IO::read_mesh(mesh, "../../data/plate/quadTest.ply");

    VectorType geom;

    QuadMeshTopologySaver quadTopol(mesh);

    size_t num_vertices = quadTopol.getNumVertices();

    QuadMeshTopologySaver::getGeometry(mesh,geom);

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

    // Create random number generator
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-0.05, 0.05); // range: [-0.05, 0.05]

    /*
    for (int i = 0; i < num_vertices; i++)
    {
        geom[3 * i]     += distribution(generator); // x
        geom[3 * i + 1] += distribution(generator); // y
        geom[3 * i + 2] += distribution(generator); // z
    }*/
    ConstraintSqrdReduced<DefaultConfigurator> constraint(quadTopol);
    ConstraintSqrdReducedGradient<DefaultConfigurator> constraintGrad(quadTopol);

    VectorType test_input = geom;
    
    RealType dest;
    VectorType Dest;

    constraint.apply(test_input,dest);
    constraintGrad.apply(test_input,Dest);

    ScalarValuedDerivativeTester<DefaultConfigurator> tester(constraint,constraintGrad,0.01,50);
    //VectorValuedDerivativeTester<DefaultConfigurator> tester(constraint, constraintGrad,0.01, Dest.rows());
    tester.plotAllDirections(test_input,"deriv_test/");

}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}
return 0;
}