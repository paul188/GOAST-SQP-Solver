#include <goast/Developability/Developability.h>
#include <goast/Core.h>
#include <random>
#include <iostream>
#include <utility>
#include <goast/Smoothers.h>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACK

int main()
{
    using VectorType = typename DefaultConfigurator::VectorType;
    using MatrixType = typename DefaultConfigurator::SparseMatrixType;
    using VecType = typename DefaultConfigurator::VecType;

    try{

        // ---------------------- TOPOLOGICAL PREPARATIONS --------------------------------------

        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "../../data/plate/basePlateSharpCreaseDevelopability.ply");
        QuadMeshTopologySaver quadTopol(mesh);
        VectorType quadDeformedGeometry, quadReferenceGeometry;
        QuadMeshTopologySaver::getGeometry(mesh,quadDeformedGeometry);
        QuadMeshTopologySaver::getGeometry(mesh,quadReferenceGeometry);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        // Read in the initial deformed quad geometry
        OpenMesh::IO::read_mesh(mesh, "/home/s24pjoha_hpc/goast_old_old/goast/build/examples/smoothed_quad_newest.ply");
        QuadMeshTopologySaver::getGeometry(mesh, quadDeformedGeometry);

        VectorType factors(2);
        factors[0] = 0.0; 
        factors[1] = 1.0; 

        int counter_1 = 0;
        int counter_2 = 0;
        std::vector<int> bdryMask;
        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            VecType coords;
            coords[0] = quadReferenceGeometry[3*i];
            coords[1] = quadReferenceGeometry[3*i+1];
            coords[2] = quadReferenceGeometry[3*i+2];

            if((std::abs(coords[0]) < 1e-4) || (std::abs(1.0 - coords[0]) < 1e-3) )
            {
                counter_1++;
                bdryMask.push_back(3*i);
                bdryMask.push_back(3*i + 1);
                bdryMask.push_back(3*i + 2);
            }

            if((std::abs(coords[1]) < 1e-4) || (std::abs(1.0 - coords[1]) < 1e-3) )
            {
                counter_2++;
                bdryMask.push_back(3*i);
                bdryMask.push_back(3*i + 1);
                bdryMask.push_back(3*i + 2);
            }
        }

        /*
        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            VecType coords;
            coords[0] = quadDeformedGeometry[3*i];
            coords[1] = quadDeformedGeometry[3*i+1];
            coords[2] = quadDeformedGeometry[3*i+2];

            coords[1] *= sqrt(1.0 - (0.4*0.4));
            coords[0] *= sqrt(1.0 - (0.4*0.4));

            quadDeformedGeometry[3*i] = coords[0];
            quadDeformedGeometry[3*i + 1] = coords[1];
            quadDeformedGeometry[3*i + 2] = coords[2];
        }*/

        QuadMeshTopologySaver::setGeometry(mesh, quadDeformedGeometry);
        OpenMesh::IO::write_mesh(mesh, "new_2.ply");

        EdgeLengthQuadEnergy<DefaultConfigurator> E_edge(quadTopol, quadReferenceGeometry);
        EdgeLengthQuadEnergyGradient<DefaultConfigurator> DE_edge(quadTopol, quadReferenceGeometry);

        ConstraintSqrdReduced<DefaultConfigurator> constraintSqrdReduced(quadTopol);
        ConstraintSqrdReducedGradient<DefaultConfigurator> constraintSqrdReducedGradient(quadTopol);

        AdditionOp<DefaultConfigurator> E_tot( factors, E_edge ,constraintSqrdReduced);
        AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_edge, constraintSqrdReducedGradient);

        std::cout<<"Hello1"<<std::endl;
        VectorType test = VectorType::Zero(quadDeformedGeometry.size());
        constraintSqrdReducedGradient.apply(quadDeformedGeometry, test);
        std::cout<<"Hello2"<<std::endl;

        if(test.hasNaN())
        {
            std::cout<<"test has NaN!"<<std::endl;
        }

        /*
        for(int i = 0; i < test.size(); i++)
        {
            if(!(test[i]).isFinite())
            {
                std::cout<<"test["<<i<<"] = "<<test[i]<<" is not finite!"<<std::endl;
            }
        }*/
        
        // set outer optimization parameters
        OptimizationParameters<DefaultConfigurator> optPars;
        optPars.setGradientIterations( 200 );
        optPars.setBFGSIterations( 200 );
        optPars.setEpsilon( 1e-10);
        optPars.setQuietMode( SHOW_ALL );

        VectorType result = quadDeformedGeometry;

        GradientDescent<DefaultConfigurator> GD( constraintSqrdReduced, constraintSqrdReducedGradient, optPars);
        GD.setBoundaryMask( bdryMask );
        GD.solve( quadDeformedGeometry, result );

        QuasiNewtonBFGS<DefaultConfigurator> BFGS( constraintSqrdReduced, constraintSqrdReducedGradient, optPars);
        BFGS.setBoundaryMask( bdryMask );
        BFGS.solve( result, result );

        QuadMeshTopologySaver::setGeometry(mesh, result);
        OpenMesh::IO::write_mesh(mesh, "centroid_base_deformed_newest.ply");

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}