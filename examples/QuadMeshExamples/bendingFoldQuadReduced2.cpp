#include "goast/Developability/Constraints.h"
#include "goast/QuadMesh/QuadTopology.h"
#include "goast/SQP/Utils/SparseMat.h"
#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <goast/Developability/Developability.h>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACK

int main()
{
    using VectorType = typename DefaultConfigurator::VectorType;
    using MatrixType = typename DefaultConfigurator::SparseMatrixType;
    using VecType = typename DefaultConfigurator::VecType;
    using RealType = typename DefaultConfigurator::RealType;

    try{

        // ---------------------- TOPOLOGICAL PREPARATIONS --------------------------------------

        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "/home/s24pjoha_hpc/goast_old_old/goast/data/plate/quadMesh2.ply");
        QuadMeshTopologySaver quadTopol(mesh);
        VectorType quadReferenceGeometry;
        QuadMeshTopologySaver::getGeometry(mesh,quadReferenceGeometry);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();
        double tolerance = 1e-4;

        VectorType quadDeformedGeometry = quadReferenceGeometry;

        std::vector<int> bdryMask;

        for(int i = 0; i < num_vertices; i++)
        {
            VecType coords;
            coords[0] = quadReferenceGeometry[3*i];
            coords[1] = quadReferenceGeometry[3*i+1];
            coords[2] = quadReferenceGeometry[3*i+2];
            if(std::abs(coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance || std::abs(coords[1]) < tolerance || std::abs(1.0 - coords[1]) < tolerance){
                bdryMask.push_back(3*i);
                bdryMask.push_back(3*i + 1);
                bdryMask.push_back(3*i + 2);
            }

            if(std::abs(coords[0]) < tolerance || std::abs(1-coords[0]) < tolerance)
            {
                if(coords[1] <= 0.5)
                {
                    coords[2] = 0.4*coords[1];
                }
                else if(coords[1] > 0.5)
                {
                    coords[2] = 0.4*(1.0 - coords[1]);
                }
                quadDeformedGeometry[3*i] = coords[0];
                quadDeformedGeometry[3*i + 1] = coords[1];
                quadDeformedGeometry[3*i + 2] = coords[2];
                continue;
            }

            if(std::abs(coords[1]) < tolerance || std::abs(1-coords[1]) < tolerance)
            {
                if(coords[0] <= 0.5)
                {
                    coords[2] = 0.4*coords[0];
                }
                else if(coords[0] > 0.5)
                {
                    coords[2] = 0.4*(1.0 - coords[0]);
                }
                quadDeformedGeometry[3*i] = coords[0];
                quadDeformedGeometry[3*i + 1] = coords[1];
                quadDeformedGeometry[3*i + 2] = coords[2];
                continue;
            }
        }

        QuadMeshTopologySaver::setGeometry(mesh, quadDeformedGeometry);
        OpenMesh::IO::write_mesh(mesh, "quad_base_deformed_0.ply");

        VectorType factors(2);
        factors[0] = 0.0; 
        factors[1] = 1.0; 

        EdgeLengthQuadEnergy<DefaultConfigurator> E_edge(quadTopol, quadReferenceGeometry);
        EdgeLengthQuadEnergyGradient<DefaultConfigurator> DE_edge(quadTopol, quadReferenceGeometry);

        ConstraintSqrdReduced<DefaultConfigurator> constraintSqrdReduced(quadTopol);
        ConstraintSqrdReducedGradient<DefaultConfigurator> constraintSqrdReducedGradient(quadTopol);

        AdditionOp<DefaultConfigurator> E_tot( factors, E_edge ,constraintSqrdReduced);
        AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_edge, constraintSqrdReducedGradient);
        
        // set outer optimization parameters
        OptimizationParameters<DefaultConfigurator> optPars;
        optPars.setGradientIterations( 200 );
        optPars.setBFGSIterations( 2000 );
        optPars.setEpsilon( 1e-10);
        optPars.setQuietMode( SHOW_ALL );

        VectorType result = quadDeformedGeometry;

        GradientDescent<DefaultConfigurator> GD( E_tot, DE_tot, optPars);
        GD.setBoundaryMask( bdryMask );
        GD.solve( quadDeformedGeometry, result );

        QuasiNewtonBFGS<DefaultConfigurator> BFGS( E_tot, DE_tot, optPars);
        BFGS.setBoundaryMask( bdryMask );
        BFGS.solve( result, result );

        QuadMeshTopologySaver::setGeometry(mesh, result);
        OpenMesh::IO::write_mesh(mesh, "centroid_base_deformed.ply");

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}