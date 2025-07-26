#include <goast/Core.h>
#include <goast/Developability/Developability.h>


int main()
{
    using VectorType = typename DefaultConfigurator::VectorType;
    using MatrixType = typename DefaultConfigurator::SparseMatrixType;
    using VecType = typename DefaultConfigurator::VecType;

    try{

        // ---------------------- TOPOLOGICAL PREPARATIONS --------------------------------------

        MyMesh mesh;
        OpenMesh::IO::read_mesh(mesh, "../../data/plate/testPlateBoundaryDOFsDevelopability.ply");
        QuadMeshTopologySaver quadTopol(mesh);
        VectorType geom;
        QuadMeshTopologySaver::getGeometry(mesh,geom);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        std::vector<int> bdryMask;

        // first, define a boundary mask:
        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            VecType coords;
            coords[0] = geom[3*i];
            coords[1] = geom[3*i+1];
            coords[2] = geom[3*i+2];
            if(std::abs(coords[0]) < 1e-4)
            {
                bdryMask.push_back(i);
            }
        }

        std::cout<<"TEST OUTPUTTING ALL VERTICES: "<<std::endl;
        for(int i = 0; i < quadTopol.getNumVertices(); i++)
        {
            std::cout<<i<<" : "<<geom[3*i]<<"; "<<geom[3*i+1]<<"; "<<geom[3*i+2]<<std::endl;
        }

        std::cout<<"OUTPUTTING BOUNDARY MASK: "<<std::endl;
        for(int i = 0; i < bdryMask.size(); i++)
        {
            std::cout<<bdryMask[i]<<" : "<<geom[3*bdryMask[i]]<<"; "<<geom[3*bdryMask[i]+1]<<"; "<<geom[3*bdryMask[i]+2]<<std::endl;
        }

        VarsIdx<DefaultConfigurator> varsIdx(quadTopol);

        BoundaryDOFS_quad<DefaultConfigurator> bdryDOFs(bdryMask, num_vertices, varsIdx);

        VectorType vars;
        vars.resize(varsIdx["num_dofs"]);
        for(int i = 0; i < vars.size(); i++)
        {
            vars[i] = i;
        }

        VectorType varsReduced;
        bdryDOFs.transformToReducedSpace(vars,varsReduced);

        std::cout<<"TEST VARS REDUCED: "<<std::endl;
        for(int i = 0; i < varsReduced.size(); i++)
        {
            std::cout<<varsReduced[i]<<" ";
        }

        VectorType varsFullAgain;
        bdryDOFs.InverseTransform(varsReduced, varsFullAgain);

        for(int i = 0; i < varsFullAgain.size(); i++)
        {
            std::cout<<varsFullAgain[i]<<" ";
        }

        // ---------------------- WORKS UNTIL HERE --------------------------------------
        StripHandler<DefaultConfigurator> stripHandler(quadTopol);
        ConstraintIdx<DefaultConfigurator> constraintIdx(stripHandler,quadTopol);
        //ConstraintGrad<DefaultConfigurator> constraintGrad(quadTopol,stripHandler,
        MatrixType Id(constraintIdx["num_cons"],varsIdx["num_dofs"]);
        MatrixType Dest;
        Id.setZero();
        for(int i = 0; i  < constraintIdx["num_cons"]; i++)
        {
            Id.coeffRef(i,i) = 1.0;
        }
        bdryDOFs.transformColsToReducedSpace(Id, Dest);
        std::cout<<"Hello1"<<std::endl;
        printSparseMatrix(Dest, "/home/paul_johannssen/Desktop/masterarbeit/goast/build/examples/IdReduced.txt");
        std::cout<<"Hello2"<<std::endl;
    }
    catch(std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        std::cerr << e.what() << std::endl;
        return -1;
    }
    return 0;
}