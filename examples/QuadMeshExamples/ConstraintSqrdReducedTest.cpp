#include <goast/Core.h>
#include <goast/Developability/Developability.h>
#include <goast/QuadMesh/QuadTopology.h>
#include <iostream>
#include <ctime>
#include <string>
#include <vector>
#include <random>

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

int main(int argc, char *argv[])
{

try{

    MyMesh mesh;
    
    OpenMesh::IO::read_mesh(mesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/quadTest.ply");

    VectorType geom;

    QuadMeshTopologySaver quadTopol(mesh);

    size_t num_vertices = quadTopol.getNumVertices();

    QuadMeshTopologySaver::getGeometry(mesh,geom);

    VectorType quadGeomRef = geom;

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

    
    for (int i = 0; i < num_vertices; i++)
    {
        geom[3 * i]     += distribution(generator); // x
        geom[3 * i + 1] += distribution(generator); // y
        geom[3 * i + 2] += distribution(generator); // z
    }
    ConstraintSqrdReduced<DefaultConfigurator> constraint(quadTopol);
    ConstraintSqrdReducedGradient<DefaultConfigurator> constraintGrad(quadTopol);
    ConstraintSqrdReducedHessian<DefaultConfigurator> constraintHess(quadTopol);

    //EdgeLengthQuadEnergy<DefaultConfigurator> constraint2(quadTopol, quadGeomRef);
    ///EdgeLengthQuadEnergyGradient<DefaultConfigurator> constraintGrad2(quadTopol, quadGeomRef);
    //EdgeLengthQuadEnergyHessian<DefaultConfigurator> constraintHess2(quadTopol, quadGeomRef);

    //EdgeLengthQuadEnergy2<DefaultConfigurator> constraint(quadTopol, quadGeomRef);
    //EdgeLengthQuadEnergyGradient2<DefaultConfigurator> constraintGrad(quadTopol, quadGeomRef);

    VectorType test_input = geom;
    
    RealType dest;
    VectorType Dest;
    MatrixType Dest_mat;

    constraint.apply(test_input,dest);
    constraintGrad.apply(test_input,Dest);
    constraintHess.apply(test_input,Dest_mat);

    //ScalarValuedDerivativeTester<DefaultConfigurator> tester(constraint,constraintGrad,0.01,50);
    VectorValuedDerivativeTester<DefaultConfigurator> tester(constraintGrad, constraintHess,0.01, Dest_mat.rows());
    tester.plotAllDirections(test_input,"/lustre/scratch/data/s24pjoha_hpc-results/deriv_test_61/");

}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}
return 0;
}