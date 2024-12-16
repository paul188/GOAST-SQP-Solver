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
 
  // need status attributes for decimation
  VertexAttributes( OpenMesh::Attributes::Normal );
  EdgeAttributes(  OpenMesh::Attributes::Color );
  FaceAttributes(  OpenMesh::Attributes::Normal );  
};

typedef OpenMesh::PolyMesh_ArrayKernelT<QuadTraits> MyMesh;

class QuadMesh{

    friend class Constraint;

    public:
        QuadMesh(){}
    
        // Checks if a vertex already exists, if not adds it -> returns the vertex handle
        VertexHandle add_vertex_unique(MyMesh::Point p)
        {
            // Check if the vertex already exists
            for(const auto &vh : _mesh.vertices())
            {
                if((_mesh.point(vh) - p).norm() < 1e-6){
                    return vh;
                }
            }
            VertexHandle vh = _mesh.add_vertex(p);
            return vh;
        }
    
        // Add a quad face given four vertex handles -> connected in the order
        // v1 -> v2 -> v3 -> v4
        FaceHandle add_quad(std::vector<MyMesh::Point>& quad_points) {
        if (quad_points.size() != 4) {
            std::cerr << "Error: Face does not have exactly four vertices!" << std::endl;
            return OpenMesh::FaceHandle();
        }

        VertexHandle vh1 = add_vertex_unique(quad_points[0]);
        VertexHandle vh2 = add_vertex_unique(quad_points[1]);
        VertexHandle vh3 = add_vertex_unique(quad_points[2]);
        VertexHandle vh4 = add_vertex_unique(quad_points[3]);

        std::vector<MyMesh::VertexHandle> face_vhandles = {vh1, vh2, vh3, vh4};
        FaceHandle fh = _mesh.add_face(face_vhandles);

        if (!fh.is_valid()) {
            std::cerr << "Error: Could not create a quad face!" << std::endl;
        }
        
        return fh;
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


        void write_quad_mesh(const std::string &filename) const {
            OpenMesh::IO::write_mesh(_mesh, filename);
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