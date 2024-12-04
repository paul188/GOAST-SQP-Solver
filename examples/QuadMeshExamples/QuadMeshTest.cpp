#include <goast/Quads.h>
#include <goast/Core.h>

int main(){
    QuadMesh qm;

    //middle row of the mesh
    MyMesh::Point p1 = MyMesh::Point(0.0,0,0);
    MyMesh::Point p2 = MyMesh::Point(1,0,0);
    MyMesh::Point p3 = MyMesh::Point(1,1,0);
    MyMesh::Point p4 = MyMesh::Point(0,1,0);
    MyMesh::Point p5 = MyMesh::Point(2,0,0);
    MyMesh::Point p6 = MyMesh::Point(2,1,0);
    MyMesh::Point p7 = MyMesh::Point(-1,0,0);
    MyMesh::Point p8 = MyMesh::Point(-1,1,0);

    //top row of the mesh
    MyMesh::Point p9 = MyMesh::Point(-1,2,0);
    MyMesh::Point p10 = MyMesh::Point(0,2,0);
    MyMesh::Point p11 = MyMesh::Point(1,2,0);
    MyMesh::Point p12 = MyMesh::Point(2,2,0);

    //bottom row of the mesh
    MyMesh::Point p13 = MyMesh::Point(-1,-1,0);
    MyMesh::Point p14 = MyMesh::Point(0,-1,0);
    MyMesh::Point p15 = MyMesh::Point(1,-1,0);
    MyMesh::Point p16 = MyMesh::Point(2,-1,0);

    //middle row of the mesh
    std::vector<MyMesh::Point> quad_points = {p1,p2,p3,p4};
    std::vector<MyMesh::Point> quad_points2 = {p2,p5,p6,p3};
    std::vector<MyMesh::Point> quad_points3 = {p7,p1,p4,p8};

    //top row of the mesh
    std::vector<MyMesh::Point> quad_points4 = {p8,p4,p10,p9};
    std::vector<MyMesh::Point> quad_points5 = {p4,p3,p11,p10};
    std::vector<MyMesh::Point> quad_points6 = {p3,p6,p12,p11};
    
    //bottom row of the mesh
    std::vector<MyMesh::Point> quad_points7 = {p13,p14,p1,p7};
    std::vector<MyMesh::Point> quad_points8 = {p14,p15,p2,p1};
    std::vector<MyMesh::Point> quad_points9 = {p15,p16,p5,p2};

    FaceHandle fh = qm.add_quad(quad_points);
    FaceHandle fh2 = qm.add_quad(quad_points2);
    FaceHandle fh3 = qm.add_quad(quad_points3);

    FaceHandle fh4 = qm.add_quad(quad_points4);
    FaceHandle fh5 = qm.add_quad(quad_points5);
    FaceHandle fh6 = qm.add_quad(quad_points6);

    FaceHandle fh7 = qm.add_quad(quad_points7);
    FaceHandle fh8 = qm.add_quad(quad_points8);
    FaceHandle fh9 = qm.add_quad(quad_points9);

    std::cout<<qm.is_mesh_developable()<<std::endl;

    qm.write_quad_mesh("quadMesh.ply");

    // Now, test with cube instead -> not developable
    QuadMesh qm2;

    MyMesh::Point q1 = MyMesh::Point(0.0,0,0);
    MyMesh::Point q2 = MyMesh::Point(0,1,0);
    MyMesh::Point q3 = MyMesh::Point(1,1,0);
    MyMesh::Point q4 = MyMesh::Point(1,0,0);
    MyMesh::Point q5 = MyMesh::Point(0.0,0,1);
    MyMesh::Point q6 = MyMesh::Point(0,1,1);
    MyMesh::Point q7 = MyMesh::Point(1,0,1);
    MyMesh::Point q8 = MyMesh::Point(1,1,1);

    std::vector<MyMesh::Point> quad_points_cube1 = {q1,q2,q3,q4};
    std::vector<MyMesh::Point> quad_points_cube2 = {q1,q5,q6,q2};
    std::vector<MyMesh::Point> quad_points_cube3 = {q1,q4,q7,q5};
    std::vector<MyMesh::Point> quad_points_cube4 = {q2,q6,q8,q3};
    std::vector<MyMesh::Point> quad_points_cube5 = {q3,q8,q7,q4};
    std::vector<MyMesh::Point> quad_points_cube6 = {q5,q7,q8,q6};

    qm2.add_quad(quad_points_cube1);
    qm2.add_quad(quad_points_cube2);
    qm2.add_quad(quad_points_cube3);
    qm2.add_quad(quad_points_cube4);
    qm2.add_quad(quad_points_cube5);
    qm2.add_quad(quad_points_cube6);

    std::cout<<qm2.is_mesh_developable()<<std::endl;

    TriMesh triMesh;
    triMesh = qm.makeTriMeshCentroid();

    OpenMesh::IO::write_mesh(triMesh, "triMesh.ply");

    TriMesh triMesh2;
    triMesh2 = qm2.makeTriMeshCentroid();

    OpenMesh::IO::write_mesh(triMesh2, "triMesh2.ply");

    return 0;
}