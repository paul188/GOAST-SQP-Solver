#pragma once

#include "Constraints.h"
#include <goast/SQP/Utils/SparseMat.h>

// --------------------- Define all the energy classes -----------------------------

template<typename ConfiguratorType>
typename ConfiguratorType::RealType square(typename ConfiguratorType::RealType x){
    return x = x*x;
}

template<typename ConfiguratorType>
typename ConfiguratorType::RealType square_norm(typename ConfiguratorType::VectorType x){
    return x.dot(x);
}

typedef std::vector<std::tuple<int,int,int>> StripType;

template<typename ConfiguratorType>
class Energy_vertex : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;

    QuadMeshTopologySaver &_quadTopol;
    VectorView<ConfiguratorType> _view;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;

public:
    Energy_vertex(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
    {
        VectorView<ConfiguratorType> _view(_quadTopol);
    }

    void apply(const VectorType &vars, RealType &Dest) const override{

        Dest = 0;

        _view.set_vector(vars);

        Dest.setZero();

        // combine the two vertex constraints into one vector
        for(int i = 0; i < _num_vertices; i++){
            Dest += square(view.vertex_weight(i) - (view.dummy_weight(i)*view.dummy_weight(i)) - 1.0);
            Dest += square_norm(view.reweighted_vertex(i) - view.vertex_weight(i)*view.vertex(i)); 
        }
    }
};

template<typename ConfiguratorType>
class Energy_edge : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;

    QuadMeshTopologySaver &_quadTopol;
    VectorView<ConfiguratorType> _view;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;

public:
    Energy_edge(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
    {
        VectorView<ConfiguratorType> _view(_quadTopol);
    }

    void apply(const VectorType &vars, RealType &Dest) const override{

        Dest = 0;

        _view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            // First, obtain first vertex idces of the face
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w_0 = view.vertex_weight(node_0);
            RealType w_1 = view.vertex_weight(node_1);
            RealType w_2 = view.vertex_weight(node_2);
            RealType w_3 = view.vertex_weight(node_3);

            VectorType rw_vertex_0 = view.reweighted_vertex(node_0);
            VectorType rw_vertex_1 = view.reweighted_vertex(node_1);
            VectorType rw_vertex_2 = view.reweighted_vertex(node_2);
            VectorType rw_vertex_3 = view.reweighted_vertex(node_3);

            Dest += square_norm(view.reweighted_edge_1(i) - (w_0 + w_1)*(rw_vertex_2+rw_vertex_3) + (w_2 + w_3)*(rw_vertex_1 + rw_vertex_0));
            Dest += square_norm(view.reweighted_edge_2(i) - (w_1 + w_2)*(rw_vertex_3+rw_vertex_0) + (w_0 + w_3)*(rw_vertex_1 + rw_vertex_2));
        }
    }
};

template<typename ConfiguratorType>
class Energy_normal : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;

    QuadMeshTopologySaver &_quadTopol;
    VectorView<ConfiguratorType> _view;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;

public:
    Energy_normal(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
    {
        VectorView<ConfiguratorType> _view(_quadTopol);
    }

    void apply(const VectorType &vars, RealType &Dest) const override{

        Dest = 0;

        _view.set_vector(vars);

        for(int i = 0; i < _num_faces; i++){
            Dest += square(view.face_normal(i).dot(view.reweighted_edge_1(i)));
            Dest += square(view.face_normal(i).dot(view.reweighted_edge_2(i)));
            Dest += square(view.face_normal(i).dot(view.face_normal(i)) - 1.0);
        }
    }
};

template<typename ConfiguratorType>
class Energy_ruling : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;

    QuadMeshTopologySaver &_quadTopol;
    VectorView<ConfiguratorType> _view;
    SkippingBdryFaceIterator _it_face;
    SkippingBdryHalfEdgeIterator _it_he;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;
public:
    Energy_ruling(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
    {
        _view(_quadTopol);
        _it_he(_quadTopol);
        _it_face(_quadTopol);
    }

    void apply(const VectorType &vars, RealType &Dest) const override{

        _view.set_vector(vars);

        for(_it_he;_it_he.valid();_it_he++){
            int i = _it_he.idx();
            int i_nobdry = _it_he.idx_nobdry();
            int face_1 = _quadTopol.getFaceOfHalfEdge(i,0);
            int face_2 = _quadTopol.getFaceOfHalfEdge(i,1);
            Dest  += square_norm(view.ruling(i_nobdry) - (view.face_normal(face_1).template head<3>()).cross(view.face_normal(face_2).template head<3>())); 
        }

        for(_it_face;it.valid();_it_face++){
            int i = it.value();
            int i_nobdry = it.value_nobdry();
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w0 = view.vertex_weight(node_0);
            RealType w1 = view.vertex_weight(node_1);
            RealType w2 = view.vertex_weight(node_2);
            RealType w3 = view.vertex_weight(node_3);

            VectorType ruling_01 = view.ruling(_it_he.to_nbdry(_quadTopol.getHalfEdgeOfQuad(i,0)));//.getHalfEdgeOfQuad(i,0));
            VectorType ruling_12 = view.ruling(_it_he.to_nbdry(_quadTopol.getHalfEdgeOfQuad(i,1)));
            VectorType ruling_23 = view.ruling(_it_he.to_nbdry(_quadTopol.getHalfEdgeOfQuad(i,2)));
            VectorType ruling_30 = view.ruling(_it_he.to_nbdry(_quadTopol.getHalfEdgeOfQuad(i,3)));

            Dest += square_norm(view.reweighted_ruling_1(i_nobdry) - (w1 + w2)*ruling_12 + (w3 + w0)*(ruling_30));
            Dest += square_norm(view.reweighted_ruling_2(i_nobdry) - (w0 + w1)*ruling_01 + (w2 + w3)*(ruling_23));
        }
    }
};

template<typename ConfiguratorType>
class Energy_dev : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    QuadMeshTopologySaver &_quadTopol;
    VectorView<ConfiguratorType> _view;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;
    SkippingBdryFaceIterator _it_face;
public:
    Energy_dev(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
    {
        _view(_quadTopol);
        _it_face(_quadTopol);
    }

    void apply(const VectorType &vars, RealType &Dest) const override{

        _view.set_vector(vars);

        for(_it_face;it.valid();_it_face++){
            int i = it.value();
            int i_nobdry = it.value_nobdry();
            Dest += square(_view.reweighted_ruling_1(i_nobdry). template<3>()).cross(_view.reweighted_ruling_2(i_nobdry). template<3>());
        }
    }
};

template<typename ConfiguratorType>
class Energy_v_fair : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    QuadMeshTopologySaver &_quadTopol;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;
    StripType _vertex_strips;
public:
    Energy_v_fair(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
    {
        // First, run over all faces to find all face strips
        _vertex_strips = _quadTopol.getVertexStrips();
    }
    void apply(const VectorType &vars, VectorType &Dest) const override{
        Dest = 0;
        VectorView<ConfiguratorType> view(vars, _quadTopol);
        for(int i = 0; i < vertex_strips.size(); i++){
            int node_0 = std::get<0>(vertex_strips[i]);
            int node_1 = std::get<1>(vertex_strips[i]);
            int node_2 = std::get<2>(vertex_strips[i]);
            Dest += square_norm(view.vertex(node_0) - 2*view.vertex(node_1) + view.vertex(node_2));
        }
    }
};

template<typename ConfiguratorType>
class Energy_n_fair : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    QuadMeshTopologySaver &_quadTopol;
    VectorView<ConfiguratorType> _view;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;
    StripType _face_strips;
public:
    Energy_n_fair(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
    {
        // First, run over all faces to find all face strips
        _face_strips = _quadTopol.getFaceStrips();
    }
    void apply(const VectorType &vars, VectorType &Dest) const override{
        Dest = 0;
        _view.set_vector(vars);

        for(int i = 0; i < _face_strips.size(); i++){
            VectorType normal_0 = _view.face_normal(std::get<0>(_face_strips[i]));
            VectorType normal_1 = _view.face_normal(std::get<1>(_face_strips[i]));
            VectorType normal_2 = _view.face_normal(std::get<2>(_face_strips[i]));
            Dest += square_norm(normal_0 - 2*normal_1 + normal_2);
        }
    }
};

template<typename ConfiguratorType>
class Energy_r_fair : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    QuadMeshTopologySaver &_quadTopol;
    VectorView<DefaultConfigurator> _view;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;
    StripType he_strips;
    SkippingBdryHalfEdgeIterator _it_he;
public:
    Energy_r_fair(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
    {
        // First, run over all faces to find all face strips
        he_strips = _quadTopol.getHalfEdgeStrips();
    }
    void apply(const VectorType &vars, RealType &Dest) const override{
        Dest = 0;
        _view.set_vector(vars);

        for(int i = 0; i < he_strips.size(); i++){
            int heh_0 = std::get<0>(he_strips[i]);
            int heh_1 = std::get<1>(he_strips[i]);
            int heh_2 = std::get<2>(he_strips[i]);

            VectorType ruling_0 = _view.ruling(_it_he.to_nbdry(heh_0));
            VectorType ruling_1 = _view.ruling(_it_he.to_nbdry(heh_1));
            VectorType ruling_2 = _view.ruling(_it_he.to_nbdry(heh_2));
            Dest += square_norm(ruling_0 - 2*ruling_1 + ruling_2);
        }
    }
};

// --------------------- Define all the gradients of energies -----------------------------
template<typename ConfiguratorType>
class EnergyGrad_vertex : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    QuadMeshTopologySaver &_quadTopol;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;

public:
    void apply(const VectorType &vars, VectorType &Dest) const override{
        VectorView<ConfiguratorType> view(vars, _quadTopol);

        size_t num_dofs = view._idx["num_dofs"];

        if(Dest.size() != num_dofs){
            Dest.resize(num_dofs);
        }

        for(int i = 0; i < _num_vertices; i++){

            RealType w_i = view.vertex_weight(i);
            RealType omega_i = view.dummy_weight(i);

            // derivatives of the first constraint term
            RealType constr_1 = w_i - (omega_i*omega_i) - 1.0;
            Dest[i + view._idx["weights"]] += 2.0*constr_1;
            Dest[i + view._idx["dummy_weights"]] += -4.0*omega_i*constr_1;

            // derivatives of the second constraint term
            VectorType constr_2 = (view.reweighted_vertex(i) - w_i*view.vertex(i)).template head<3>();
            Dest.segment(i + view._idx["reweighted_vertices"],3) += 2*constr_2;
            Dest.segment(i + view._idx["vertices"],3) += -2*w_i*constr_2;
            Dest[view._idx["weights"] + i] += 2*constr_2.dot(view.vertex(i));
        }
    }
};

template<typename ConfiguratorType>
class EnergyGrad_edge : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    QuadMeshTopologySaver &_quadTopol;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;

public:
    EnergyGrad_edge(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
    {}
    void apply(const VectorType &vars, VectorType &Dest) const override{
        VectorView<ConfiguratorType> view(vars, _quadTopol);
        size_t num_dofs = view._idx["num_dofs"];

        if(Dest.size() != num_dofs){
            Dest.resize(num_dofs);
        }

        // Compute derivatives of the first constraint
        for(int i = 0; i < _num_faces; i++){
            // First, obtain first vertex idces of the face
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            RealType w_0 = view.vertex_weight(node_0);
            RealType w_1 = view.vertex_weight(node_1);
            RealType w_2 = view.vertex_weight(node_2);
            RealType w_3 = view.vertex_weight(node_3);

            VectorType rw_vertex_0 = view.reweighted_vertex(node_0);
            VectorType rw_vertex_1 = view.reweighted_vertex(node_1);
            VectorType rw_vertex_2 = view.reweighted_vertex(node_2);
            VectorType rw_vertex_3 = view.reweighted_vertex(node_3);

            VectorType constr_1 = view.reweighted_edge_1(i) - (w_0 + w_1)*(rw_vertex_2+rw_vertex_3) + (w_2 + w_3)*(rw_vertex_1 + rw_vertex_0);
            VectorType constr_2 = view.reweighted_edge_2(i) - (w_1 + w_2)*(rw_vertex_3+rw_vertex_0) + (w_0 + w_3)*(rw_vertex_1 + rw_vertex_2);

            // First, add all derivatives of the first constraint
            Dest.segment(3*i + view._idx["reweighted_edges_1"],3) += 2*constr_1;
            Dest.segment(3*node_2 + view._idx["reweighted_vertices"]) += -2*(w_0 + w_1)*constr_1;
            Dest.segment(3*node_3 + view._idx["reweighted_vertices"]) += -2*(w_0 + w_1)*constr_1;
            Dest.segment(3*node_0 + view._idx["reweighted_vertices"]) += 2*(w_2 + w_3)*constr_1;
            Dest.segment(3*node_1 + view._idx["reweighted_vertices"]) += 2*(w_2 + w_3)*constr_1;


            // Add all derivatives of the second constraint
            Dest.segment(3*i + view._idx["reweighted_edges_2"],3) += 2*constr_2;
            Dest.segment(3*node_0 + view._idx["reweighted_vertices"],3) += -2*(w_1 + w_2)*constr_2;
            Dest.segment(3*node_3 + view._idx["reweighted_vertices"],3) += -2*(w_1 + w_2)*constr_2;
            Dest.segment(3*node_1 + view._idx["reweighted_vertices"],3) += 2*(w_0 + w_3)*constr_2;
            Dest.segment(3*node_2 + view._idx["reweighted_vertices"],3) += 2*(w_0 + w_3)*constr_2;
        }
    }
};

template<typename ConfiguratorType>
class EnergyGrad_normal : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    QuadMeshTopologySaver &_quadTopol;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;
public:
    EnergyGrad_normal(QuadMeshTopologySaver &quadTopol)
        : _quadTopol(quadTopol),
          _num_vertices(quadTopol.getNumVertices()),
          _num_faces(quadTopol.getNumFaces()),
          _num_edges(quadTopol.getNumEdges()),
          _num_halfedges(quadTopol.getNumHalfEdges())
    {}
    void apply(const VectorType &vars, VectorType &Dest) const override{
        VectorView<ConfiguratorType> view(vars, _quadTopol);
        size_t num_dofs = view._idx["num_dofs"];

        for(int i = 0; i < _num_faces; i++){
            // derivatives of first constraint
            Dest[3*i + view._idx["face_normals"]] += 2*view.reweighted_edge_1(i);
            Dest[3*i + view._idx["reweighted_edges_1"]] += 2*view.face_normal(i); 
            // derivatives of second constraint
            Dest[3*i + view._idx["face_normals"]] += 2*view.reweighted_edge_2(i);
            Dest[3*i + view._idx["reweighted_edges_2"]] += 2*view.face_normal(i);
            // derivatives of third constraint
            Dest[3*i + view._idx["face_normals"]] += 4*view.face_normal(i);
        }
    }
};

template<typename ConfiguratorType>
class EnergyGrad_ruling : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    QuadMeshTopologySaver &_quadTopol;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;
public:
    void apply(const VectorType &vars, VectorType &Dest) const override{
        VectorView<ConfiguratorType> view(vars, _quadTopol);
        size_t num_dofs = view._idx["num_dofs"];
        
        for(SkippingBdryHalfEdgeIterator it(_quadTopol);it.valid();it++){
            int i = it.idx();
            int i_nobdry = _quadTopol.to_nbdry(i);
            int face_1 = _quadTopol.getFaceOfHalfEdge(i,0);
            int face_2 = _quadTopol.getFaceOfHalfEdge(i,1);

            // FIRST, CALCULATE DERIVATES OF CONSTRAINT ZERO

            // First, calculate the constraint vector
            VectorType constr0 = view.ruling(i_nobdry) - (view.face_normal(face_1).template head<3>()).cross(view.face_normal(face_2).template head<3>());
            
            // Calculate the derivatives of the constraint w.r.t. the first normal vector
            FullMatrixType Dconstr0_dn1(3);
            Eigen::Vector3d e = VectorType::Zero(3);
            e[0] = 1.0;
            Dconstr0_dn1.block(0,0,3,1) = - e.cross(view.face_normal(face_2).template head<3>());
            e[0] = 0.0;
            e[1] = 1.0;
            Dconstr0_dn1.block(0,1,3,1) = - e.cross(view.face_normal(face_2).template head<3>());
            e[1] = 0.0;
            e[2] = 1.0;
            Dconstr0_dn1.block(0,2,3,1) = - e.cross(view.face_normal(face_2).template head<3>());

            // Calculate the derivatives of the constraint w.r.t. the second normal vector
            FullMatrixType Dconstr0_dn2(3);
            e[0] = 1.0;
            e[2] = 0.0;
            Dconstr0_dn2.block(0,0,3,1) = view.face_normal(face_1).template head<3>().cross(e);
            e[0] = 0.0;
            e[1] = 1.0;
            Dconstr0_dn2.block(0,1,3,1) = view.face_normal(face_1).template head<3>().cross(e);
            e[1] = 0.0;
            e[2] = 1.0;
            Dconstr0_dn2.block(0,2,3,1) = view.face_normal(face_1).template head<3>().cross(e);

            Dest.segment(3*i_nobdry + view._idx["rulings"],3) += 2*constr0;
            Dest.segment(3*face_1 + view._idx["face_normals"],3) += 2*Dconstr0_dn1*constr0;
            Dest.segment(3*face_2 + view._idx["face_normals"],3) += 2*Dconstr0_dn2*constr0;
        }
        
        // NEXT, CALULATE DERIVATES OF CONSTRAINT ONE
        SkippingBdryHalfEdgeIterator it_he(_quadTopol);

        for(SkippingBdryFaceIterator it(_quadTopol);it.valid();it++){

            int i = it.idx();
            int i_nobdry = it.idx_nobdry();
        
            int node_0 = _quadTopol.getNodeOfQuad(i,0);
            int node_1 = _quadTopol.getNodeOfQuad(i,1);
            int node_2 = _quadTopol.getNodeOfQuad(i,2);
            int node_3 = _quadTopol.getNodeOfQuad(i,3);

            int heh_01 = _quadTopol.getHalfEdgeOfQuad(i,0);
            int heh_12 = _quadTopol.getHalfEdgeOfQuad(i,1);
            int heh_23 = _quadTopol.getHalfEdgeOfQuad(i,2);
            int heh_30 = _quadTopol.getHalfEdgeOfQuad(i,3);

            RealType w0 = view.vertex_weight(node_0);
            RealType w1 = view.vertex_weight(node_1);
            RealType w2 = view.vertex_weight(node_2);
            RealType w3 = view.vertex_weight(node_3);

            VectorType ruling_01 = view.ruling(it_he.to_nbdry(heh_01));
            VectorType ruling_12 = view.ruling(it_he.to_nbdry(heh_12));
            VectorType ruling_23 = view.ruling(it_he.to_nbdry(heh_23));
            VectorType ruling_30 = view.ruling(it_he.to_nbdry(heh_30));

            VectorType constr_1 = view.reweighted_ruling_1(i_nobdry) - (w1 + w2)*ruling_12 + (w3 + w0)*(ruling_30);

            VectorType Dconstr_1_dw0 = ruling_30;
            VectorType Dconstr_1_dw1 = -ruling_12;
            VectorType Dconstr_1_dw2 = -ruling_12;
            VectorType Dconstr_1_dw3 = ruling_30;

            Dest.segment(view._idx["weights"] + node_0,1) += 2*Dconstr_1_dw0.dot(constr_1);
            Dest.segment(view._idx["weights"] + node_1,1) += 2*Dconstr_1_dw1.dot(constr_1);
            Dest.segment(view._idx["weights"] + node_2,1) += 2*Dconstr_1_dw2.dot(constr_1);
            Dest.segment(view._idx["weights"] + node_3,1) += 2*Dconstr_1_dw3.dot(constr_1);

            Dest.segment(3*i_nobdry + view._idx["reweighted_rulings_1"],3) += 2*constr_1;
            Dest.segment(3*it_he.to_nbdry(heh_12) + view._idx["rulings"],3) += -2*(w1+w2)*constr_1;
            Dest.segment(3*it_he.to_nbdry(heh_30) + view._idx["rulings"],3) += 2*(w3+w0)*constr_1;

            // NEXT, CALCULATE DERIVATES OF CONSTRAINT TWO
            VectorType constr_2 = view.reweighted_ruling_2(i_nobdry) - (w0 + w1)*ruling_01 + (w2 + w3)*(ruling_23);

            VectorType Dconstr_2_dw0 = -ruling_01;
            VectorType Dconstr_2_dw1 = -ruling_01;
            VectorType Dconstr_2_dw2 = ruling_23;
            VectorType Dconstr_2_dw3 = ruling_23;

            Dest.segment(view._idx["weights"] + node_0,1) += 2*Dconstr_2_dw0.dot(constr_2);
            Dest.segment(view._idx["weights"] + node_1,1) += 2*Dconstr_2_dw1.dot(constr_2);
            Dest.segment(view._idx["weights"] + node_2,1) += 2*Dconstr_2_dw2.dot(constr_2);
            Dest.segment(view._idx["weights"] + node_3,1) += 2*Dconstr_2_dw3.dot(constr_2);

            Dest.segment(3*i_nobdry + view._idx["reweighted_rulings_2"],3) += 2*constr_2;

            Dest.segment(3*it_he.to_nbdry(heh_01) + view._idx["rulings"],3) += -2*(w0+w1)*constr_2;
            Dest.segment(3*it_he.to_nbdry(heh_23) + view._idx["rulings"],3) += 2*(w2+w3)*constr_2;
        }
    }
};

template<typename ConfiguratorType>
class EnergyGrad_dev : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

    QuadMeshTopologySaver &_quadTopol;
    size_t _num_vertices;
    size_t _num_faces;
    size_t _num_edges;
    size_t _num_halfedges;
public:
    void apply(const VectorType &vars, VectorType &Dest) const override{
        VectorView<ConfiguratorType> view(_quadTopol);
        view.set_vector(vars);

        for(SkippingBdryFaceIterator it(_quadTopol);it.valid();it++){
            int i = it.idx();
            int i_nobdry = it.idx_nobdry();
            
            VectorType constr = view.reweighted_ruling_1(i_nobdry). template<3>().cross(view.reweighted_ruling_2(i_nobdry). template<3>());
            FullMatrixType Dconstr_1(3,3);
            Vec3d e = Vec3d::Zero();
            e[0] = 1.0;
            Dconstr_1.col(0) = e.cross(view.reweighted_ruling_1(i_nobdry). template<3>());
            e[0] = 0.0;
            e[1] = 1.0;
            Dconstr_1.col(1) = e.cross(view.reweighted_ruling_1(i_nobdry). template<3>());
            e[1] = 0.0;
            e[2] = 1.0;
            Dconstr_1.col(2) = e.cross(view.reweighted_ruling_1(i_nobdry). template<3>());

            Dest.segment(3*i_nobdry + view._idx["reweighted_rulings_1"],3) += 2*Dconstr_1*constr;

            e[2] = 0.0;
            e[0] = 1.0;
            FullMatrixType Dconstr_2(3,3);
            Dconstr_2.col(0) = view.reweighted_ruling_2(i_nobdry). template<3>().cross(e);
            e[0] = 0.0;
            e[1] = 1.0;
            Dconstr_2.col(1) = view.reweighted_ruling_2(i_nobdry). template<3>().cross(e);
            e[1] = 0.0;
            e[2] = 1.0;
            Dconstr_2.col(2) = view.reweighted_ruling_2(i_nobdry). template<3>().cross(e);

            Dest.segment(3*i_nobdry + view._idx["reweighted_rulings_2"],3) += 2*Dconstr_2*constr;
        }      
    }
};