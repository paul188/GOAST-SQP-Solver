#ifndef QUADMESH_H
#define QUADMESH_H

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/DefaultTriMesh.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <vector>

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
  VertexAttributes( OpenMesh::Attributes::Status );
  EdgeAttributes( OpenMesh::Attributes::Status );
  FaceAttributes( OpenMesh::Attributes::Status );  
};

typedef OpenMesh::PolyMesh_ArrayKernelT<QuadTraits> MyMesh;

class QuadMesh{
    public:
        OpenMesh::HPropHandleT<OpenMesh::Vec2d> _heWeights;

        QuadMesh(MyMesh &mesh): _mesh(mesh){
            /*
             specify needed properties of the QuadMesh
             1. We want to have weights that specify where our contact element
                points lie, as interpolates
             2. We want to have the contact element point 
                as the oriented sum of the weights of the halfedge
            
            in the standard constructor, we simply take the edge midpoints and all
            weights = 0.5
             as vertices of the contact elements

            -> generally, we need to ensure that the weights are ordered in the same
             manner, i.e. the opposite halfedge should have weights in opposite order

            -> see function bool weights_consistent()
            */

            // Add the property to the mesh
            HPropHandleT<OpenMesh::Vec2d> heWeights;
            _mesh.add_property<OpenMesh::Vec2d>(heWeights, "Contact Element Edge Weights");
            for(const auto &heh : _mesh.halfedges())
            {
                _mesh.property(heWeights, heh) = OpenMesh::Vec2d(0.5,0.5);
            }
            this->_heWeights = heWeights;
        }

        // The weights should always be ordered with the same orientation
        // as the halfedges along the face
        QuadMesh(MyMesh &mesh, const std::vector<OpenMesh::Vec2d> &heWeights):_mesh(mesh)
        {
            this -> _mesh = mesh;

            if(heWeights.size() != _mesh.n_halfedges()){
                std::cerr << "Error: Number of edge weights does not match number of halfedges!" << std::endl;
                return;
            }

              // Add the property to the mesh
            HPropHandleT<OpenMesh::Vec2d> heWeights_temp;
            _mesh.add_property<OpenMesh::Vec2d>(heWeights_temp, "Contact Element Edge Weights");
            for(const auto &heh : _mesh.halfedges())
            {
                _mesh.property(heWeights_temp, heh) = heWeights[heh.idx()];
            }
            this->_heWeights = heWeights_temp;
        }
    
        QuadTraits::Point get_edge_midpoint(const HalfedgeHandle &heh){
            
            // TGet the first and last vertex of the halfedge
            auto firstvh = _mesh.from_vertex_handle(heh);
            auto lastvh = _mesh.to_vertex_handle(heh);

            QuadTraits::Point p0 = _mesh.point(firstvh);
            QuadTraits::Point p1 = _mesh.point(lastvh);
            return 1/(_mesh.property(_heWeights,heh)[0] + _mesh.property(_heWeights,heh)[0])*(_mesh.property(_heWeights,heh)[0]*p0 + _mesh.property(_heWeights,heh)[0]*(p1-p0));
        }
        
        QuadTraits::Point get_edge_midpoint(const EdgeHandle &eh){
            // Obtain the existing property
            auto _heWeights = HProp<OpenMesh::Vec2d>(_mesh,"Contact Element Edge Weights");
            
            //obtain some halfedge handle -> doesnt matter if 0 or 1
            HalfedgeHandle heh = _mesh.halfedge_handle(eh,0);

            // TGet the first and last vertex of the halfedge
            auto firstvh = _mesh.from_vertex_handle(heh);
            auto lastvh = _mesh.to_vertex_handle(heh);

            QuadTraits::Point p0 = _mesh.point(firstvh);
            QuadTraits::Point p1 = _mesh.point(lastvh);
            return 1/(_heWeights[heh][0] + _heWeights[heh][1])*(_heWeights[heh][0]*p0 + _heWeights[heh][1]*(p1-p0));
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
            return _mesh.normal(fh);
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
            HalfedgeHandle oppheh = _mesh.opposite_halfedge_handle(heh);
            if(_mesh.is_boundary(heh) || _mesh.is_boundary(oppheh)){
                return Vec3d();
            }

            // Get handles of the adjacent faces
            FaceHandle f1 = _mesh.face_handle(heh);
            FaceHandle f2 = _mesh.face_handle(oppheh);

            Vec3d normal_1 = get_face_normal_vec(f1);
            Vec3d normal_2 = get_face_normal_vec(f2);
            
            return normal_1.cross(normal_2);
        }

        // Check if the weights of the halfedges are consistent
        // and if they are all > 0
        bool midpoint_weights_valid()
        {

            for (const auto &eh : _mesh.edges())
            {
                auto heh0 = _mesh.halfedge_handle(eh, 0);
                auto heh1 = _mesh.halfedge_handle(eh, 1);

                if(!(_mesh.property(_heWeights,heh0)[0] == _mesh.property(_heWeights,heh1)[1] && _mesh.property(_heWeights,heh0)[1] == _mesh.property(_heWeights,heh1)[0])){
                    std::cerr<<"Weights of opposite halfedges are not consistent!"<<std::endl;
                    return false;
                }

                if(_mesh.property(_heWeights,heh0)[0] <= 0 || _mesh.property(_heWeights,heh0)[1] <= 0){
                    std::cerr<<"Weights of halfedge are not positive!"<<std::endl;
                    return false;
                }

                return true;
            }
        }

        // bool firstEdge specifies which ruling vector to take
        Vec3d get_avg_ruling_vec(const FaceHandle &fh, bool firstEdge)
        {

            double tot_faceWeight = 0;

            for(const auto &heh : _mesh.fh_range(fh))
            {
                tot_faceWeight += (_mesh.property(_heWeights,heh)[0] + _mesh.property(_heWeights,heh)[1]);
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

            double weight1 = _mesh.property(_heWeights,heh1)[0] + _mesh.property(_heWeights,heh1)[1];

            heh2 = _mesh.next_halfedge_handle(_mesh.next_halfedge_handle(heh1));

            double weight2 = _mesh.property(_heWeights,heh2)[0] + _mesh.property(_heWeights,heh2)[1];

            return 1/tot_faceWeight*(weight1*ruling1 - weight2*ruling2);
        }

        bool is_face_developable(const FaceHandle &fh)
        {
            // if it is on the boundary, we say true by default
            if(_mesh.is_boundary(fh)){
                std::cout<<"Trivial"<<std::endl;
                return true;
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
        
        // Add a quad face given four vertex handles
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


            // Check if all faces are quads
        bool is_all_quads() const {
            for (auto fh : _mesh.faces()) {
                if (_mesh.valence(fh) != 4) {
                    std::cerr << "Non-quad face found with valence: " << _mesh.valence(fh) << std::endl;
                    return false;
                }
            }
            return true;
        }

    protected:
        MyMesh _mesh;

};

#endif