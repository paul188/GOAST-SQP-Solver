#ifndef QUADMESH_H
#define QUADMESH_H

#include <goast/Core.h>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <vector>
#include <unordered_set>
#include <goast/DiscreteShells.h>
//#include <goast/QuadMesh/ClassicalEnergies.h>

using namespace OpenMesh;

struct QuadTraits : public DefaultTraits
{
  /// Use double precision points
  typedef OpenMesh::Vec3d Point;
  /// Use double precision Normals
  typedef OpenMesh::Vec3d Normal;
  /// Use RGBA Color
  typedef OpenMesh::Vec4f Color;  
  typedef OpenMesh::Vec2d he_weights;
 
  // need status attributes for decimation
  VertexAttributes( OpenMesh::Attributes::Normal );
  EdgeAttributes(  OpenMesh::Attributes::Color );
  FaceAttributes(  OpenMesh::Attributes::Normal );  
};

typedef OpenMesh::PolyMesh_ArrayKernelT<QuadTraits> MyMesh;

class QuadMesh{
    public:
        OpenMesh::VPropHandleT<double> _vertexWeights;
        OpenMesh::FPropHandleT<OpenMesh::Vec3d> _faceNormals;

        QuadMesh(){
                _mesh.add_property(_vertexWeights, "vertex_weights");
                _mesh.property(_vertexWeights).set_persistent(true);
                _mesh.add_property(_faceNormals, "face_normals");
                _mesh.property(_faceNormals).set_persistent(true);
        }
    
        QuadTraits::Point get_edge_midpoint(const HalfedgeHandle &heh){
            
            // Get the first and last vertex of the halfedge
            auto firstvh = _mesh. from_vertex_handle(heh);
            auto lastvh = _mesh. to_vertex_handle(heh);

            QuadTraits::Point p0 = _mesh. point(firstvh);
            QuadTraits::Point p1 = _mesh. point(lastvh);
            return 1/(_mesh. property(_vertexWeights,firstvh) + _mesh. property(_vertexWeights,lastvh))*(_mesh. property(_vertexWeights,firstvh)*p0 + _mesh. property(_vertexWeights,lastvh)*p1);
        }
        // returns one of the average edge vectors of a face
        // the firstEdge bool specifies which average edge Vector
        // the firstEdge respects orientation of the face
        Vec3d get_avg_edge_vec(const FaceHandle &fh, bool firstEdge)
        {
            OpenMesh::ArrayKernel::HalfedgeHandle firstheh;
            OpenMesh::ArrayKernel::HalfedgeHandle secondheh;
            
            if(firstEdge){
                firstheh = _mesh.halfedge_handle(fh);
            }else{
                firstheh = _mesh.next_halfedge_handle(_mesh.halfedge_handle(fh));
            }

            // Now, traverse two halfedges to the opposite halfedge
            secondheh = _mesh.next_halfedge_handle(_mesh.next_halfedge_handle(firstheh));

            //Now, calculate the vector
            auto firstmidpoint = get_edge_midpoint(firstheh);
            auto secondmidpoint = get_edge_midpoint(secondheh);

            return (secondmidpoint - firstmidpoint);
        }

        Vec3d get_face_normal_vec(const FaceHandle &fh)
        {
            return _mesh.property(_faceNormals,fh);
        }
        
        // By choosing which of the two halfedges to take, we can
        // choose the orientation of the ruling vector
        // because the first normal in the cross product is that
        // of the face belonging to the halfedge
        Vec3d get_ruling_vec(const HalfedgeHandle &heh)
        {
            // Check if the halfedge lies on the boundary
            // -> either the halfedge or its opposite halfedge
            // return nothing by default -> undefined
            HalfedgeHandle oppheh = _mesh. opposite_halfedge_handle(heh);
            if(_mesh. is_boundary(heh) || _mesh. is_boundary(oppheh)){
                return Vec3d();
            }

            // Get handles of the adjacent faces
            FaceHandle f1 = _mesh. face_handle(heh);
            FaceHandle f2 = _mesh. face_handle(oppheh);

            Vec3d normal_1 = get_face_normal_vec(f1);
            Vec3d normal_2 = get_face_normal_vec(f2);
            
            return normal_1.cross(normal_2);
        }

        // check if all weights positive
        bool midpoint_weights_valid()
        {
            for(const auto &vh: _mesh.vertices())
            {
                if(_mesh.property(_vertexWeights,vh) <= 0){
                    return false;
                }
            }
            return true;
        }

        // bool firstEdge specifies which ruling vector to take
        Vec3d get_avg_ruling_vec(const FaceHandle &fh, bool firstEdge)
        {

            double tot_faceWeight = 0;

            for(const auto &vh : _mesh. fv_range(fh))
            {
                tot_faceWeight += _mesh.property(_vertexWeights,vh);
            }

            HalfedgeHandle heh1;
            HalfedgeHandle heh2;

            if(firstEdge)
            {
                heh1 = _mesh.halfedge_handle(fh);
            }
            else
            {
                heh1 = _mesh.next_halfedge_handle(_mesh.halfedge_handle(fh));
            }

            Vec3d ruling1 = get_ruling_vec(heh1);
            Vec3d ruling2 = get_ruling_vec(_mesh.next_halfedge_handle(heh1));

            auto vh1 = _mesh.from_vertex_handle(heh1);
            auto vh2 = _mesh.to_vertex_handle(heh1);

            double weight1 = _mesh.property(_vertexWeights,vh1) + _mesh.property(_vertexWeights,vh2);

            heh2 = _mesh.next_halfedge_handle(_mesh.next_halfedge_handle(heh1));

            vh1 = _mesh.from_vertex_handle(heh2);
            vh2 = _mesh.to_vertex_handle(heh2);

            double weight2 = _mesh.property(_vertexWeights,vh1) + _mesh.property(_vertexWeights,vh2);

            return 1/tot_faceWeight*(weight1*ruling1 - weight2*ruling2);
        }

        bool is_face_developable(const FaceHandle &fh)
        {
            // if it is on the boundary, we say true by default
            if(_mesh.is_boundary(fh)){
                return true;
            }

            if(!fh.is_valid()){
                std::cout<<"Invalid face handle"<<std::endl;
            }

            // Get the average ruling vectors
            Vec3d ruling1 = get_avg_ruling_vec(fh, true);
            Vec3d ruling2 = get_avg_ruling_vec(fh, false);

            return ruling1.cross(ruling2).norm() < 1e-6;
        }

        // Check if the entire mesh is developable
        bool is_mesh_developable()
        {
            for(const auto &fh : _mesh.faces())
            {
                if(!is_face_developable(fh)){
                    return false;
                }
            }
            return true;
        }

        bool is_valid_handle(VertexHandle& vh)
        {
            return _mesh.is_valid_handle(vh);
        }

        // Checks if a vertex already exists, if not adds it -> returns the vertex handle
        VertexHandle add_vertex_unique(MyMesh::Point p, double weight)
        {
            // Check if the vertex already exists
            for(const auto &vh : _mesh.vertices())
            {
                if((_mesh.point(vh) - p).norm() < 1e-6){
                    return vh;
                }
            }
            VertexHandle vh = _mesh.add_vertex(p);

            _mesh.property(_vertexWeights,vh) = weight;
            return vh;
        }
    
        // Add a quad face given four vertex handles -> connected in the order
        // v1 -> v2 -> v3 -> v4
        // For every of the eight new halfedges that are created this way,
        // need weight input parameter
        FaceHandle add_quad(std::vector<MyMesh::Point>& quad_points, std::vector<double> weights = {}) {
        if (quad_points.size() != 4) {
            std::cerr << "Error: Face does not have exactly four vertices!" << std::endl;
            return OpenMesh::FaceHandle();
        }

        // Initialize with standard weights
        if(weights.empty()){
            weights = {0.5,0.5,0.5,0.5};
        }
        else{
            if(weights.size()!=4){
                std::cerr << "Error: Need exactly four weights for the vertices!" << std::endl;
                return OpenMesh::FaceHandle();
            }
        }

        // Add the vertices and add the new weights. 
        // keeps already existing weights if the vertex already exists
        VertexHandle vh1 = add_vertex_unique(quad_points[0], weights[0]);
        VertexHandle vh2 = add_vertex_unique(quad_points[1], weights[1]);
        VertexHandle vh3 = add_vertex_unique(quad_points[2], weights[2]);
        VertexHandle vh4 = add_vertex_unique(quad_points[3], weights[3]);

        std::vector<MyMesh::VertexHandle> face_vhandles = {vh1, vh2, vh3, vh4};
        FaceHandle fh = _mesh.add_face(face_vhandles);

        if (!fh.is_valid()) {
            std::cerr << "Error: Could not create a quad face!" << std::endl;
        }

        // update the face normals

        Vec3d avg_edge_1 = get_avg_edge_vec(fh, true);
        Vec3d avg_edge_2 = get_avg_edge_vec(fh, false);

        Vec3d normal = avg_edge_1.cross(avg_edge_2);
        _mesh.property(_faceNormals,fh) = normal.normalized();
        return fh;
    }

    void write_quad_mesh(const std::string &filename) const {
        OpenMesh::IO::write_mesh(_mesh, filename);
    }

    VertexHandle add_vertex_unique_Tri(TriMesh &triMesh, TriMesh::Point p)
    {
        // Check if the vertex already exists
        for(const auto &vh : triMesh.vertices())
        {
            if((triMesh.point(vh) - p).norm() < 1e-6){
                return vh;
            }
        }
        VertexHandle vh = triMesh.add_vertex(p);
        return vh;
    }

    // method to split a quad face into a centroid.
    // Need a start_heh to preserve the orientation of the mesh
    void split_quad(FaceHandle &fh, HalfedgeHandle &start_heh, TriMesh &triMesh){
        
        // This is the outgoing halfedge from vertex vh along the face fh
        VertexHandle v1 = _mesh.from_vertex_handle(start_heh);
        auto heh2 = _mesh.next_halfedge_handle(start_heh);
        VertexHandle v2 = _mesh.from_vertex_handle(heh2);
        VertexHandle v3 = _mesh.to_vertex_handle(heh2);
        VertexHandle v4 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh2));

        VertexHandle v1_tri = add_vertex_unique_Tri(triMesh,_mesh.point(v1));
        VertexHandle v2_tri = add_vertex_unique_Tri(triMesh,_mesh.point(v2));
        VertexHandle v3_tri = add_vertex_unique_Tri(triMesh,_mesh.point(v3));
        VertexHandle v4_tri = add_vertex_unique_Tri(triMesh,_mesh.point(v4));
        VertexHandle v5_tri = add_vertex_unique_Tri(triMesh,0.25*(_mesh.point(v1) + _mesh.point(v2) + _mesh.point(v3) + _mesh.point(v4)));

        FaceHandle trifh_1 = triMesh.add_face(v1_tri,v2_tri,v5_tri);
        FaceHandle trifh_2 = triMesh.add_face(v2_tri,v3_tri,v5_tri);
        FaceHandle trifh_3 = triMesh.add_face(v3_tri,v4_tri,v5_tri);
        FaceHandle trifh_4 = triMesh.add_face(v4_tri,v1_tri,v5_tri);
        
        if(!trifh_1.is_valid() || !trifh_2.is_valid() || !trifh_3.is_valid() || !trifh_4.is_valid()){
            std::cerr << "Error: Could not create a triangle face!" << std::endl;
        }
    }

    TriMesh makeTriMeshCentroid(){
        // The result trimesh
        TriMesh triMesh;

        // Add a property that signals the orientation of the face -> one vertex of the face that lies on the diagonal
        // We consider the halfedge that starts from that vertex
        OpenMesh::FPropHandleT<OpenMesh::HalfedgeHandle> centroid_heh;
        OpenMesh::FPropHandleT<bool> visited;
        _mesh.add_property(centroid_heh, "centroid_heh");
        _mesh.add_property(visited, "visited");

        // We have not visited any faces yet -> means not added to face_stack yet
        // and don't have a diagonal halfedge yet
        for (const auto& fh : _mesh.faces()) {
            _mesh.property(visited, fh) = false;
        }

        // Start with some face handle of the mesh
        FaceHandle init_face;

        for (auto fh : _mesh.faces()) {
            init_face = fh;
            break; // Take the first face
        }

        if (!init_face.is_valid()) {
            std::cerr << "Error: No valid face found!" << std::endl;
        }

        // Stack for DFS traversal
        std::vector<MyMesh::FaceHandle> face_stack;

        // Push the starting face
        face_stack.push_back(init_face);
        _mesh.property(visited, init_face) = true;

        VertexHandle init_vertex;
        HalfedgeHandle init_heh;

        init_heh = _mesh.halfedge_handle(init_face);
        init_vertex = _mesh.from_vertex_handle(init_heh);

        _mesh.property(centroid_heh,init_face) = init_heh;

        // The halfedge handle from the current face going out from a vertex on the diagonal
        HalfedgeHandle start_heh;
        while (!face_stack.empty()) {
            // Pop a face from the stack
            FaceHandle current_face = face_stack.back();
            face_stack.pop_back();

            start_heh = _mesh.property(centroid_heh,current_face);

            // Store all adjacent non triangle faces
            auto current_heh = start_heh;
            int counter = 0;
            do {
                if(! _mesh.is_boundary(current_heh)){
                    FaceHandle adj_face = _mesh.opposite_face_handle(current_heh);
                    if (adj_face.is_valid() && !_mesh.property(visited, adj_face)) {
                        face_stack.push_back(adj_face);
                        // Set the halfedge handle property of the adjacent face
                        HalfedgeHandle centroid_property = _mesh.opposite_halfedge_handle(current_heh);
                        _mesh.property(centroid_heh,adj_face) = centroid_property;
                        _mesh.property(visited, adj_face) = true;
                    }
                }
                current_heh = _mesh.next_halfedge_handle(current_heh);
                counter ++;
            } while(current_heh != start_heh);

            // Split the current face and add to the trimesh
            split_quad(current_face, start_heh, triMesh);
        }
        _mesh.remove_property(centroid_heh);
        _mesh.remove_property(visited);

        return triMesh;
    }
    
    bool is_quad_mesh(){
        for(const auto &fh : _mesh.faces()){
            if(!is_quad_face(fh)){
                return false;
            }
        }
        return true;
    }

    bool is_quad_face(const FaceHandle &fh) const {
        HalfedgeHandle start_heh = _mesh.halfedge_handle(fh);
        if(_mesh.next_halfedge_handle(start_heh) == _mesh.next_halfedge_handle(_mesh.next_halfedge_handle(_mesh.next_halfedge_handle(_mesh.next_halfedge_handle(start_heh))))){
            return true;
        }
        return false;
    }

    void read_overwrite_quad_mesh(const std::string &filename){
        OpenMesh::IO::read_mesh(_mesh, filename);

        // Add the properties
        for(auto vh : _mesh.vertices()){
            _mesh.property(_vertexWeights,vh) = 0.5;
        }

        for(auto fh : _mesh.faces()){
            _mesh.property(_faceNormals,fh) = get_face_normal_vec(fh);
        }

        if(!is_quad_mesh()){
            std::cerr << "Error: Mesh is not a quad mesh!" << std::endl;
        }
    }

    int num_edges(){
        return _mesh.n_edges();
    }

    int num_vertices(){
        return _mesh.n_vertices();
    }

    int num_faces(){
        return _mesh.n_faces();
    }
    

    // -----------------Declarations of all SimpleBending energy objects -----------

    template<typename ConfiguratorType>
    class SimpleBendingEnergy;

    template<typename ConfiguratorType>
    class SimpleBendingGradientDef;

    template<typename ConfiguratorType>
    class SimpleBendingHessianDef;

    template<typename ConfiguratorType>
    class SimpleBendingGradientUndef;

    template<typename ConfiguratorType>
    class SimpleBendingHessianUndef;

    template<typename ConfiguratorType>
    class SimpleBendingHessianMixed;

    // Declarations end

    // -----------------Declarations of all NonlinearMembrane energy objects -----------

    template<typename ConfiguratorType>
    class NonlinearMembraneEnergy;

    template<typename ConfiguratorType>
    class NonlinearMembraneGradientDef;

    template<typename ConfiguratorType>
    class NonlinearMembraneHessianDef;

    template<typename ConfiguratorType>
    class NonlinearMembraneGradientUndef;

    template<typename ConfiguratorType>
    class NonlinearMembraneHessianUndef;

    template<typename ConfiguratorType>
    class NonlinearMembraneHessianMixed;

    // Declarations end

    protected:
        MyMesh _mesh;

};

// -----------------Definitions of all SimpleBending energy objects -----------

template<typename ConfiguratorType>
class QuadMesh::SimpleBendingEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:

    typedef typename ConfiguratorType::RealType   RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType    VecType;
    typedef typename ConfiguratorType::MatType    MatType;

    const ::SimpleBendingEnergy<ConfiguratorType> _E_bend1;
    const ::SimpleBendingEnergy<ConfiguratorType> _E_bend2;

public:
    // only need one TriMesh plate object, since we only extract vertex positions
    // -> the same for both diagonal orientations
    SimpleBendingEnergy(const MeshTopologySaver &plateTopol1,
                        const MeshTopologySaver &plateTopol2, 
                        const VectorType &undefShell,
                        const bool ActiveShellIsDeformed, 
                        const VectorType& Weight) :
                        _E_bend1(plateTopol1, undefShell, ActiveShellIsDeformed, Weight), 
                        _E_bend2(plateTopol2, undefShell, ActiveShellIsDeformed, Weight){}

    SimpleBendingEnergy(const MeshTopologySaver &plateTopol1,
                        const MeshTopologySaver &plateTopol2,
                        const VectorType &undefShell,
                        const bool ActiveShellIsDeformed, 
                        RealType Weight = 1. ) :
                        _E_bend1(plateTopol1, undefShell, ActiveShellIsDeformed, Weight), 
                        _E_bend2(plateTopol2, undefShell, ActiveShellIsDeformed, Weight){}

    void apply( const VectorType& ActiveGeometry, RealType & Dest ) const {
        RealType Dest1, Dest2;
        _E_bend1.apply(ActiveGeometry, Dest1);
        _E_bend2.apply(ActiveGeometry, Dest2);
        Dest = 0.5*(Dest1 + Dest2);
    }
};

template<typename ConfiguratorType>
class QuadMesh::SimpleBendingGradientDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

    typedef typename ConfiguratorType::RealType   RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType    VecType;
    typedef typename ConfiguratorType::MatType    MatType;

    const ::SimpleBendingGradientDef<ConfiguratorType> _DE_bend1;
    const ::SimpleBendingGradientDef<ConfiguratorType> _DE_bend2;

public:
    SimpleBendingGradientDef(const MeshTopologySaver &plateTopol1,
                             const MeshTopologySaver &plateTopol2,
                             const VectorType &undefShell,
                             const VectorType& Weight) :
                             _DE_bend1(plateTopol1, undefShell, Weight),
                             _DE_bend2(plateTopol2, undefShell, Weight){} 

    SimpleBendingGradientDef(const MeshTopologySaver &plateTopol1,
                             const MeshTopologySaver &plateTopol2,
                             const VectorType &undefShell,
                             RealType Weight = 1.) :
                             _DE_bend1(plateTopol1, undefShell, Weight),
                             _DE_bend2(plateTopol2, undefShell, Weight){} 

    void apply(const VectorType& defShell, VectorType& Dest) const {
        VectorType Dest1, Dest2;
        _DE_bend1.apply(defShell, Dest1);
        _DE_bend2.apply(defShell, Dest2);
        Dest = 0.5*(Dest1 + Dest2);
    }
};

template<typename ConfiguratorType>
class QuadMesh::SimpleBendingHessianDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

    protected:
    typedef typename ConfiguratorType::RealType    RealType;
    typedef typename ConfiguratorType::VectorType  VectorType;
    typedef typename ConfiguratorType::SparseMatrixType  MatrixType;

    const ::SimpleBendingHessianDef<ConfiguratorType> _D2E_bend1;
    const ::SimpleBendingHessianDef<ConfiguratorType> _D2E_bend2;

public:
    SimpleBendingHessianDef(const MeshTopologySaver& plateTopol1,
                            const MeshTopologySaver& plateTopol2,
                            const VectorType &undefShell,
                            const VectorType& Weight,
                            int rowOffset = 0, 
                            int colOffset = 0 ) :
                            _D2E_bend1(plateTopol1, undefShell, Weight, rowOffset, colOffset),
                            _D2E_bend2(plateTopol2, undefShell, Weight, rowOffset, colOffset){}


    SimpleBendingHessianDef(const MeshTopologySaver& plateTopol1,
                            const MeshTopologySaver& plateTopol2,
                            const VectorType &undefShell,
                            const RealType Weight = 1,
                            int rowOffset = 0, 
                            int colOffset = 0 ) :
                            _D2E_bend1(plateTopol1, undefShell, Weight, rowOffset, colOffset),
                            _D2E_bend2(plateTopol2, undefShell, Weight, rowOffset, colOffset){}

    void apply( const VectorType& defShell, MatrixType& Dest ) const {    
        MatrixType Dest1, Dest2;
        _D2E_bend1.apply(defShell, Dest1);
        _D2E_bend2.apply(defShell, Dest2);

        Dest = 0.5*(Dest1 + Dest2);
    }

};

template<typename ConfiguratorType>
class QuadMesh::SimpleBendingGradientUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

    typedef typename ConfiguratorType::RealType   RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType    VecType;
    typedef typename ConfiguratorType::MatType    MatType;

const ::SimpleBendingGradientUndef<ConfiguratorType> _DE_bend_undef1;
const ::SimpleBendingGradientUndef<ConfiguratorType> _DE_bend_undef2;

public:
    SimpleBendingGradientUndef( const MeshTopologySaver& plateTopol1,
                                const MeshTopologySaver& plateTopol2,
                                const VectorType& defShell,
                                const VectorType& Weight) :
                                _DE_bend_undef1(plateTopol1, defShell, Weight),
                                _DE_bend_undef2(plateTopol2, defShell, Weight){}

    SimpleBendingGradientUndef( const MeshTopologySaver& plateTopol1,
                                const MeshTopologySaver& plateTopol2,
                                const VectorType& defShell,
                                const RealType Weight = 1. ) :
                                _DE_bend_undef1(plateTopol1, defShell, Weight),
                                _DE_bend_undef2(plateTopol2, defShell, Weight){}

    void apply( const VectorType& undefShell, VectorType& Dest ) const {
        VectorType Dest1, Dest2;
        _DE_bend_undef1.apply(undefShell,Dest1);
        _DE_bend_undef2.apply(undefShell,Dest2);
        Dest = 0.5*(Dest1 + Dest2);
    }

};

template<typename ConfiguratorType>
class QuadMesh::SimpleBendingHessianUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
    typedef typename ConfiguratorType::RealType    RealType;
    typedef typename ConfiguratorType::VectorType  VectorType;
    typedef typename ConfiguratorType::SparseMatrixType  MatrixType;

    const ::SimpleBendingHessianUndef<ConfiguratorType> _D2E_bend_undef1;
    const ::SimpleBendingHessianUndef<ConfiguratorType> _D2E_bend_undef2;

public:
    SimpleBendingHessianUndef(  const MeshTopologySaver& plateTopol1,
                                const MeshTopologySaver& plateTopol2,
                                const VectorType& defShell,
                                const VectorType& Weight,
                                int rowOffset = 0, 
                                int colOffset = 0 ) :
                                _D2E_bend_undef1(plateTopol1, defShell, Weight, rowOffset, colOffset),
                                _D2E_bend_undef2(plateTopol2, defShell, Weight, rowOffset, colOffset){}

    SimpleBendingHessianUndef(  const MeshTopologySaver& plateTopol1,
                                const MeshTopologySaver& plateTopol2,
                                const VectorType& defShell,
                                const RealType Weight,
                                int rowOffset = 0, 
                                int colOffset = 0 ) :
                                _D2E_bend_undef1(plateTopol1, defShell, Weight, rowOffset, colOffset),
                                _D2E_bend_undef2(plateTopol2, defShell, Weight, rowOffset, colOffset){}

    void apply( const VectorType& undefShell, MatrixType& Dest ) const {
        MatrixType Dest1, Dest2;
        _D2E_bend_undef1.apply(undefShell,Dest1);
        _D2E_bend_undef2.apply(undefShell, Dest2);

        Dest = 0.5*(Dest1 + Dest2);
    }

};


template<typename ConfiguratorType>
class QuadMesh::SimpleBendingHessianMixed : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
    typedef typename ConfiguratorType::RealType    RealType;
    typedef typename ConfiguratorType::VectorType  VectorType;
    typedef typename ConfiguratorType::SparseMatrixType  MatrixType;

    const ::SimpleBendingHessianMixed<ConfiguratorType> _D2E_bend_mixed1;
    const ::SimpleBendingHessianMixed<ConfiguratorType> _D2E_bend_mixed2;

public:
    SimpleBendingHessianMixed(  const MeshTopologySaver& plateTopol1,
                                const MeshTopologySaver& plateTopol2,
                                const VectorType& InactiveGeometry,
                                const bool ActiveShellIsDeformed,
                                const bool FirstDerivWRTDef,
                                const VectorType& Weight, 
                                int rowOffset = 0, 
                                int colOffset = 0 ) :
                                _D2E_bend_mixed1(plateTopol1, InactiveGeometry, ActiveShellIsDeformed, FirstDerivWRTDef, Weight, rowOffset, colOffset),
                                _D2E_bend_mixed2(plateTopol2, InactiveGeometry, ActiveShellIsDeformed, FirstDerivWRTDef, Weight, rowOffset, colOffset){}

    SimpleBendingHessianMixed(  const MeshTopologySaver& plateTopol1,
                                const MeshTopologySaver& plateTopol2,
                                const VectorType& InactiveGeometry,
                                const bool ActiveShellIsDeformed,
                                const bool FirstDerivWRTDef,
                                const RealType Weight = 1, 
                                int rowOffset = 0, 
                                int colOffset = 0 ) :
                                _D2E_bend_mixed1(plateTopol1, InactiveGeometry, ActiveShellIsDeformed, FirstDerivWRTDef, Weight, rowOffset, colOffset),
                                _D2E_bend_mixed2(plateTopol2, InactiveGeometry, ActiveShellIsDeformed, FirstDerivWRTDef, Weight, rowOffset, colOffset){}

    void apply( const VectorType& ActiveGeometry, MatrixType& Dest ) const {
        MatrixType Dest1, Dest2;
        _D2E_bend_mixed1.apply(ActiveGeometry, Dest1);
        _D2E_bend_mixed2.apply(ActiveGeometry, Dest2);

        Dest = 0.5*(Dest1 + Dest2);
    }

};

// -----------------Definitions end  ----------------------------------------------

// -----------------Definitions of all NonlinearMembrane energy objects -----------
template<typename ConfiguratorType>
class QuadMesh::NonlinearMembraneEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  const ::NonlinearMembraneEnergy<ConfiguratorType> _E_membrane_1;
  const ::NonlinearMembraneEnergy<ConfiguratorType> _E_membrane_2;

public:

  NonlinearMembraneEnergy( const MeshTopologySaver& topology1,
                           const MeshTopologySaver& topology2,
                           const VectorType& InactiveGeometry,
		                   const bool ActiveShellIsDeformed,
                           RealType Mu = 1., 
                           RealType Lambda = 1. ) :
                            _E_membrane_1(topology1, InactiveGeometry, ActiveShellIsDeformed, Mu, Lambda),
                            _E_membrane_2(topology2, InactiveGeometry, ActiveShellIsDeformed, Mu, Lambda){}

  // energy evaluation
  void apply( const VectorType& ActiveGeometry, RealType & Dest ) const {
    RealType Dest1, Dest2;
    _E_membrane_1.apply(ActiveGeometry, Dest1);
    _E_membrane_2.apply(ActiveGeometry, Dest2);
    Dest = 0.5*(Dest1 + Dest2);
  }
};

template<typename ConfiguratorType>
class QuadMesh::NonlinearMembraneGradientDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
    protected:
	typedef typename ConfiguratorType::RealType   RealType;

	typedef typename ConfiguratorType::VectorType VectorType;

	const ::NonlinearMembraneGradientDef<ConfiguratorType> _DE_membrane_1;
    const ::NonlinearMembraneGradientDef<ConfiguratorType> _DE_membrane_2;

public:
	NonlinearMembraneGradientDef(const MeshTopologySaver& topology1,
                                 const MeshTopologySaver& topology2,
                                 const VectorType& undefShell,
                                 RealType Mu = 1.,
                                 RealType Lambda = 1.) :
                                 _DE_membrane_1(topology1, undefShell, Mu, Lambda),
                                 _DE_membrane_2(topology2, undefShell, Mu, Lambda){}

	void apply(const VectorType& defShell, VectorType& Dest) const {
        VectorType Dest1, Dest2;
        _DE_membrane_1.apply(defShell, Dest1);
        _DE_membrane_2.apply(defShell, Dest2);
        Dest = 0.5*(Dest1 + Dest2);
    }
};

template<typename ConfiguratorType>
class QuadMesh::NonlinearMembraneHessianDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
    protected:
    typedef typename ConfiguratorType::RealType    RealType;
    typedef typename ConfiguratorType::VectorType  VectorType;
    typedef typename ConfiguratorType::SparseMatrixType  MatrixType;

    const ::NonlinearMembraneHessianDef<ConfiguratorType> _D2E_membrane_1;
    const ::NonlinearMembraneHessianDef<ConfiguratorType> _D2E_membrane_2;

public:
    NonlinearMembraneHessianDef(const MeshTopologySaver& topology1,
                                const MeshTopologySaver& topology2,
                                const VectorType& undefShell,
                                const RealType Factor = 1.,
                                int rowOffset = 0,
                                int colOffset = 0,
                                RealType Mu = 1.,
                                RealType Lambda = 1.) :
                                _D2E_membrane_1(topology1, undefShell, Factor, rowOffset, colOffset, Mu, Lambda),
                                _D2E_membrane_2(topology2, undefShell, Factor, rowOffset, colOffset, Mu, Lambda){}

    void apply(const VectorType& defShell, MatrixType& Dest) const {
        MatrixType Dest1, Dest2;
        _D2E_membrane_1.apply(defShell, Dest1);
        _D2E_membrane_2.apply(defShell, Dest2);

        Dest = 0.5*(Dest1);
    }
};

template <typename ConfiguratorType>
class QuadMesh::NonlinearMembraneHessianUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:    
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  
  const ::NonlinearMembraneHessianUndef<ConfiguratorType> _D2E_membrane_undef1;
  const ::NonlinearMembraneHessianUndef<ConfiguratorType> _D2E_membrane_undef2;
  
public:
    NonlinearMembraneHessianUndef( const MeshTopologySaver& topology1,
                                   const MeshTopologySaver& topology2, 
                                   const VectorType& defShell,
                                   RealType Factor = 1.,
                                   int rowOffset = 0, 
                                   int colOffset = 0, 
                                   RealType Mu = 1., 
                                   RealType Lambda = 1. ) :
                                   _D2E_membrane_undef1(topology1, defShell, Factor, rowOffset, colOffset, Mu, Lambda), 
                                   _D2E_membrane_undef2(topology2, defShell, Factor, rowOffset, colOffset, Mu, Lambda){}

    void apply( const VectorType& undefShell, MatrixType& Dest ) const {
        MatrixType Dest1, Dest2;
        _D2E_membrane_undef1.apply(undefShell, Dest1);
        _D2E_membrane_undef2.apply(undefShell, Dest2);

        Dest = 0.5*(Dest1 + Dest2);
    }

};

template <typename ConfiguratorType >
class QuadMesh::NonlinearMembraneHessianMixed : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;

  const ::NonlinearMembraneHessianMixed<ConfiguratorType> _D2E_membrane_mixed1;
  const ::NonlinearMembraneHessianMixed<ConfiguratorType> _D2E_membrane_mixed2;

public:
NonlinearMembraneHessianMixed( const MeshTopologySaver& topology1,
                               const MeshTopologySaver& topology2,
                               const VectorType& InactiveGeometry,
		                       const bool ActiveShellIsDeformed,
			                   const bool FirstDerivWRTDef,
			                   const RealType Factor = 1.,
                               int rowOffset = 0, 
                               int colOffset = 0,
                               RealType Mu = 1., 
                               RealType Lambda = 1. ) :
                               _D2E_membrane_mixed1(topology1, InactiveGeometry, ActiveShellIsDeformed, FirstDerivWRTDef, Factor, rowOffset, colOffset, Mu, Lambda),
                               _D2E_membrane_mixed2(topology2, InactiveGeometry, ActiveShellIsDeformed, FirstDerivWRTDef, Factor, rowOffset, colOffset, Mu, Lambda){}

void apply( const VectorType& ActiveGeometry, MatrixType& Dest ) const {
    MatrixType Dest1, Dest2;
    _D2E_membrane_mixed1.apply(ActiveGeometry, Dest1);
    _D2E_membrane_mixed2.apply(ActiveGeometry, Dest2);

    Dest = 0.5*(Dest1 + Dest2);
}
};
#endif