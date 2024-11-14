#include <goast/Quads.h>

int main(){
    MyMesh mesh;
    QuadMesh qm(mesh);

    //middle row of the mesh
    MyMesh::Point p1 = MyMesh::Point(0,0,0);
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

    std::cout<<"TEST"<<std::endl;
    std::cout<<qm.is_face_developable(fh)<<std::endl;

/*
    std::cout<<qm.is_all_quads()<<std::endl;
    std::cout<<qm.is_mesh_developable()<<std::endl;
    std::cout<<qm.is_face_developable(fh)<<std::endl;
*/

    return 0;
}