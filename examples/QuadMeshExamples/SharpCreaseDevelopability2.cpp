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
        VectorType geom;
        QuadMeshTopologySaver::getGeometry(mesh,geom);

        size_t num_vertices = quadTopol.getNumVertices();
        size_t num_faces = quadTopol.getNumFaces();

        // --------------------- GENERATE BASE CENTROID TOPOLOGY ----------------------------

        TriMesh centroid_base_mesh = quadTopol.makeQuadMeshCentroid();
        MeshTopologySaver centroidTopol(centroid_base_mesh);
        VectorType centroid_base_geom;
        getGeometry(centroid_base_mesh, centroid_base_geom);

        // ---------------------  GENERATE MAP BETWEEN QUAD AND CENTROID TOPOLOGY ----------------------------
        std::map<int,int> quadToCentroidIdxMap;
        std::vector<int> noQuadVertexIndices;
        bool foundMatch = false;
        for(int i = 0; i < centroidTopol.getNumVertices(); i++){
            VecType coords_tri;
            getXYZCoord<VectorType, VecType>(centroid_base_geom,coords_tri, i);
            for(int j = 0; j < quadTopol.getNumVertices(); j++)
            {
                VecType coords_quad;
                coords_quad[0] = geom[3*j];
                coords_quad[1] = geom[3*j+1];
                coords_quad[2] = geom[3*j+2];

                if((coords_quad - coords_tri).norm() < 1e-4)
                {
                    quadToCentroidIdxMap[j] = i;
                    foundMatch = true;
                    break;
                }
                
            }
            if(!foundMatch)
            {
                noQuadVertexIndices.push_back(i);
            }
            
            // reset foundMatch
            foundMatch = false;
        }

        // -------------------------- SETUP BOUNDARY ------------------------------------------

        double tolerance = 1e-4;

        std::vector<int> bdryMask;
        for(int i = 0; i < num_vertices; i++)
        {
            VecType coords;
            coords[0] = geom[3*i];
            coords[1] = geom[3*i+1];
            coords[2] = geom[3*i+2];
            if(std::abs(coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance || std::abs(coords[1]) < tolerance || std::abs(1.0 - coords[1]) < tolerance){
                bdryMask.push_back(i);
            }

            if(std::abs(coords[0]) < tolerance || std::abs(1-coords[0]) < tolerance)
            {
                if(coords[1] <= 0.5)
                {
                    geom[3*i + 2] = 0.4*coords[1];
                }
                else if(coords[1] > 0.5)
                {
                    geom[3*i + 2] = 0.4*(1.0 - coords[1]);
                }
                continue;
            }

            if(std::abs(coords[1]) < tolerance || std::abs(1-coords[1]) < tolerance)
            {
                if(coords[0] <= 0.5)
                {
                    geom[3*i + 2] = 0.4*coords[0];
                }
                else if(coords[0] > 0.5)
                {
                    geom[3*i + 2] = 0.4*(1.0 - coords[0]);
                }
                continue;
            }
        }

        // ----------------- GENERATE DEFORMED CENTROID MESHES --------------------------------

        QuadMeshTopologySaver::setGeometry(mesh, geom);
        QuadMeshTopologySaver new_quad_topol(mesh);
        // Now, we want to refine this into a triangle mesh
        TriMesh centroid_mesh = new_quad_topol.makeQuadMeshCentroid();
        VectorType centroid_geom;
        getGeometry(centroid_mesh, centroid_geom);
        MeshTopologySaver centroid_topol(centroid_mesh);
        OpenMesh::IO::write_mesh(centroid_mesh, "centroid.ply");
        // Now, first recalculate the boundary in the TriMesh format
        std::vector<int> bdryMaskTri;
        for(int i = 0; i < centroid_topol.getNumVertices(); i++)
        {
            VecType coords;
            getXYZCoord<VectorType, VecType>(centroid_base_geom,coords, i);

            if(std::abs(coords[0]) < tolerance || std::abs(1.0 - coords[0]) < tolerance || std::abs(coords[1]) < tolerance || std::abs(1.0 - coords[1]) < tolerance)
            {
                bdryMaskTri.push_back(i);
            }
        }

        VectorType smoothedGeom = centroid_geom;

        /*
        for(int i = 0; i < centroidTopol.getNumVertices(); i++)
        {
            VecType coords;
            getXYZCoord<VectorType, VecType>(smoothedGeom,coords, i);

            double x = coords[0];
            double y = coords[1];
            double z = coords[2];

            double eps = 1e-9;
            if (y <= x + 0.5 + eps &&
                y <= 1.5 - x + eps &&
                y >= 0.5 - x - eps &&
                y >= x - 0.5 - eps)
            {
                coords[2] = 0.2; // flatten
                bdryMaskTri.push_back(i);
            }


            if((y <= (x + 0.5)) && ((y <= (1.0 - x))) && (y >= (0.5 - x)) && (y >= (x - 0.5)))
            {
                coords[2] = 0.2;
            }
            setXYZCoord<VectorType, VecType>(smoothedGeom, coords, i);
        }*/

        for(int i = 0; i < centroidTopol.getNumVertices(); i++)
        {
            VecType coords;
            getXYZCoord<VectorType, VecType>(smoothedGeom,coords, i);
            coords[1] *= sqrt(1.0 - (0.4*0.4));
            coords[0] *= sqrt(1.0 - (0.4*0.4));
            setXYZCoord<VectorType, VecType>(smoothedGeom, coords, i);
        }
        /*
        for(int i = 0; i < centroidTopol.getNumVertices(); i++)
        {
            VecType coords;
            getXYZCoord<VectorType, VecType>(smoothedGeom,coords, i);
            
            double x = coords[0];
            double y = coords[1];
            double z = coords[2];
            if(std::abs(x - 0.5) < 0.1 && std::abs(y - 0.5) < 0.1)
            {
                coords[2] = 0.4;
                bdryMaskTri.push_back(i);
            }
            setXYZCoord<VectorType, VecType>(smoothedGeom, coords, i);
        }*/

        extendBoundaryMask(centroid_topol.getNumVertices(), bdryMaskTri);

        setGeometry(centroid_mesh, smoothedGeom);
        OpenMesh::IO::write_mesh(centroid_mesh, "smoothed_centroid_3.ply");

        // Now, apply the Dirichlet smoothing
        DirichletSmoother<DefaultConfigurator> dirichletSmoother(centroid_base_geom, bdryMaskTri, centroid_topol);
        dirichletSmoother.apply(smoothedGeom, smoothedGeom);

        setGeometry(centroid_mesh, smoothedGeom);
        OpenMesh::IO::write_mesh(centroid_mesh, "smoothed_centroid_2.ply");

        ElasticCentroi        

    }catch(std::exception &e){
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    return 0;
}