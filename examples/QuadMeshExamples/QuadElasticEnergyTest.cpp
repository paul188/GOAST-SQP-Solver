#include <goast/Core.h>
#include "goast/Core/Auxiliary.h"
#include "goast/Developability/DevelopabilityCentroidEnergy.h"
#include <iostream>
#include <ctime>
#include <string>
#include <vector>

#include <goast/QuadMesh/QuadTopology.h>
#include <goast/Developability/Developability.h>
//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

int main(int argc, char *argv[])
{

try{

    MyMesh mesh;
    
    OpenMesh::IO::read_mesh(mesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/quadMesh2.ply");

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

    TriMesh tri_base_mesh;
    OpenMesh::IO::read_mesh(tri_base_mesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/quadMesh2_refined.ply");
    MeshTopologySaver triangleTopol(tri_base_mesh);
    VectorType triBaseGeometry;
    getGeometry(tri_base_mesh, triBaseGeometry);

    VectorType _factors_elasticity_dev = VectorType::Ones(2);
    _factors_elasticity_dev[0] = 1.0;
    _factors_elasticity_dev[1] = 1.0;

    VectorType _factors_mem_bend = VectorType::Ones(2);
    _factors_mem_bend[0] = 10000.0;
    _factors_mem_bend[1] = 1.0;

    VectorType test_input = quadGeomRef;
    std::cout<<"test input size: "<<test_input.size()<<std::endl;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-0.05, 0.05); // range: [-0.05, 0.05]

    
    for (int i = 0; i < num_vertices; i++)
    {
        test_input[3*i] += distribution(generator); // x
        test_input[3*i + 1] += distribution(generator); // y
        test_input[3*i + 2] += distribution(generator); // z
    }

    RealType dest;
    VectorType Dest;

    QuadElasticEnergy<DefaultConfigurator> quadElasticEnergy(quadTopol, triangleTopol , triBaseGeometry,_factors_elasticity_dev, _factors_mem_bend);
    QuadElasticEnergyGradient<DefaultConfigurator> quadElasticEnergyGrad(quadTopol, triangleTopol, triBaseGeometry, _factors_elasticity_dev, _factors_mem_bend);
    QuadElasticEnergyHessian<DefaultConfigurator> quadElasticEnergyHess(quadTopol, triangleTopol, triBaseGeometry, _factors_elasticity_dev, _factors_mem_bend);

    quadElasticEnergy.apply(test_input, dest);
    quadElasticEnergyGrad.apply(test_input, Dest);

    std::cout<<"Energy: "<<dest<<std::endl;
    std::cout<<"Gradient size: "<<Dest.size()<<std::endl;

    //ScalarValuedDerivativeTester<DefaultConfigurator> tester(quadElasticEnergy,quadElasticEnergyGrad,0.001,50);
    VectorValuedDerivativeTester<DefaultConfigurator> tester(quadElasticEnergyGrad, quadElasticEnergyHess,0.001, Dest.rows());
    tester.plotAllDirections(test_input,"/lustre/scratch/data/s24pjoha_hpc-results/deriv_test_51/");

}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}
return 0;
}