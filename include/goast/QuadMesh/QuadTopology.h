#pragma once

//== INCLUDES =================================================================
#include "QuadMesh.h"
#include "QuadMeshIncludes.h"
#include <goast/Core/Auxiliary.h>

#include <queue>

//typedef std::vector<int>      Vec3i;
//typedef std::vector<int>      Vec2i;

//=============================================================================

/**
 * \brief Class to store all topological information of a QuadMesh
 * \author Heeren
 *
 * Class which only stores all important topology information of a given quad mesh,
 * i.e. indices, neighbors and boundary information
 */

using namespace OpenMesh;

class QuadMeshTopologySaver {

protected:  
  typedef MyMesh QuadMeshType;
  typedef OpenMesh::VectorT<signed int,8> Vec8i;
  static const int IndexNotSet = -1;  
  const QuadMeshType& _mesh;

  const int _numOfElems;
  const int _numOfNodes;
  const int _numOfEdges;

  // saves global indices of nodes belonging to each quad
  std::vector< Vec4i > _nodesOfQuads;
  // saves global indices of neighbouring quads for each quad (these are -1 at the boundary )
  std::vector< Vec4i > _neighboursOfQuads;
  // saves global indices of edges belonging to each quad
  std::vector< Vec4i > _edgesOfQuads;
  // saves global indices of neighbouring quads for each edge (these are -1 at the boundary )
  std::vector< Vec2i > _neighboursOfEdges;  
  // saves global indices of nodes belonging to each edge
  std::vector< Vec2i > _nodesOfEdges;
  // saves global indices of nodes opposite of each edge (two of them are -1 for boundary edges)
  std::vector< Vec4i > _oppositeNodesOfEdges;
  // saves the halfedge handle indices belonging to a certain face
  std::vector< Vec4i > _halfedgesOfQuads;
  // saves the face that a halfedge belongs to and the opposite one in that order (one face -1 at the boundary)
  std::vector< Vec2i > _facesOfHalfedges;
  // saves the halfedges that lie at the topological boundary
  std::vector<int> _bdryHalfEdges;
  // saves the faces that lie at the topological boundary
  std::vector<int> _bdryFaces;
   
  // stores unvalid edges and faces
  BitVector _isEdgeValid, _isFaceValid;
  // fixed nodes (e.g. for Dirichlet boundary conditions) 
  BitVector _isFixedNode; // not related to topological boundary nodes!
  int _numOfFixedNodes;

  public:
  QuadMeshTopologySaver( const QuadMeshType& mesh )
      : _mesh( mesh ),
        _numOfElems( mesh.n_faces() ),
        _numOfNodes( mesh.n_vertices() ),
        _numOfEdges( mesh.n_edges() ),
        _nodesOfQuads( _numOfElems ),
        _neighboursOfQuads( _numOfElems ),
        _edgesOfQuads( _numOfElems ),
        _neighboursOfEdges( _numOfEdges ),
        _oppositeNodesOfEdges( _numOfEdges ),
        _halfedgesOfQuads( _numOfElems ),
        _facesOfHalfedges( 2*_numOfEdges ),
        _nodesOfEdges( _numOfEdges ),
        _isEdgeValid( _numOfEdges, true ), 
        _isFaceValid( _numOfElems, true ),
        _isFixedNode( _numOfNodes ),
        _numOfFixedNodes( 0 )
  {    
    // run over all faces
    for ( MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it ){
      // global node indices of current face 
      int locIdx = 0;
      for ( MyMesh::ConstFaceVertexIter cfv_it = mesh.cfv_iter(*f_it); cfv_it.is_valid(); ++cfv_it)
        _nodesOfQuads[f_it->idx()][locIdx++] = cfv_it->idx();
      
      // Ordering consistent with ordering of the halfedges of the face
      locIdx = 0;
      for ( MyMesh::ConstFaceHalfedgeIter cfh_it = mesh.cfh_iter(*f_it); cfh_it.is_valid(); ++cfh_it){
        _neighboursOfQuads[f_it->idx()][locIdx++] = mesh.is_boundary(*cfh_it) ? IndexNotSet : mesh.face_handle( mesh.opposite_halfedge_handle(*cfh_it) ).idx();
      }
  
      //
      locIdx = 0;
      for (MyMesh::ConstFaceEdgeIter cfe_it = mesh.cfe_iter(*f_it); cfe_it.is_valid(); ++cfe_it) {
        _edgesOfQuads[f_it->idx()][locIdx++] = cfe_it->idx();
      }

      locIdx = 0;
      for (MyMesh::ConstFaceHalfedgeIter cfh_it = mesh.cfh_iter(*f_it); cfh_it.is_valid(); ++cfh_it){
        _halfedgesOfQuads[f_it->idx()][locIdx++] = cfh_it->idx();
      }

      for(MyMesh::ConstFaceHalfedgeIter cfh_it = mesh.cfh_iter(*f_it); cfh_it.is_valid(); ++cfh_it){
        if(mesh.is_boundary(*cfh_it) || mesh.is_boundary(mesh.opposite_halfedge_handle(*cfh_it))){
          // one halfedge lies at the bdry, so the face is a bdry face
          _bdryFaces.push_back(f_it->idx());
          break;
        }
      }
    }

    // run over all edges
    for (  MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it ){
       // define two corresponding half edges
       HEH heh0( mesh.halfedge_handle( *e_it, 0) );
       HEH heh1( mesh.halfedge_handle( *e_it, 1) );
       
       // indices of vertices at current edge
       _nodesOfEdges[e_it->idx()][0] = mesh.from_vertex_handle( heh0 ).idx();
       _nodesOfEdges[e_it->idx()][1] = mesh.to_vertex_handle( heh0 ).idx();
        
       // indices of adjacent faces
       _neighboursOfEdges[e_it->idx()][0] = mesh.face_handle( heh0 ).idx();
       _neighboursOfEdges[e_it->idx()][1] = mesh.face_handle( heh1 ).idx();          
     }

    // run over all halfedges
    for ( MyMesh::HalfedgeIter he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it ){
      // define two corresponding half edges
      HEH heh( *he_it );
      
      // indices of vertices at current edge
      _facesOfHalfedges[heh.idx()][0] = mesh.face_handle( heh ).idx();
      _facesOfHalfedges[heh.idx()][1] = mesh.opposite_face_handle( heh ).idx();

      if(mesh.is_boundary(heh) || mesh.is_boundary(mesh.opposite_halfedge_handle(heh))){
        // dont want to push back both halfedge and opposite halfedge as we want to preserve the ordering
        _bdryHalfEdges.push_back(heh.idx());
      }
    }

  }

  typedef std::vector<std::tuple<int,int,int>> StripType;

  // this method works, tested
  StripType getFaceStrips() const
  {
    StripType faceStrips;
    
    for(int i = 0; i < this->getNumFaces(); i++)
    {
      // now, get all halfedges of the face
      int heh_0 = this->getHalfEdgeOfQuad(i,0);
      int heh_1 = this->getHalfEdgeOfQuad(i,1);
      int heh_2 = this->getHalfEdgeOfQuad(i,2);
      int heh_3 = this->getHalfEdgeOfQuad(i,3);

      // face_0 opposite to heh_0, ..., face_3 opposite to heh_3
      int face_0 = getFaceOfHalfEdge(heh_0,1);
      int face_1 = getFaceOfHalfEdge(heh_1,1);
      int face_2 = getFaceOfHalfEdge(heh_2,1);
      int face_3 = getFaceOfHalfEdge(heh_3,1);

      // in one direction, we have found a triple strip
      if((face_0 != -1) && (face_2 != -1))
      {
        faceStrips.push_back(std::make_tuple(face_0,i,face_2));
      }

      // in the other direction, we have found a triple strip
      if((face_1 != -1) && (face_3 != -1))
      {
        faceStrips.push_back(std::make_tuple(face_1,i,face_3));
      }
    }

    return faceStrips;
  }

  // this method works, tested
  StripType getVertexStrips() const
  {
    StripType vertexStrips;

    // iterate over all triple middle vertices
    for(const auto vh_middle : _mesh.vertices()){

      auto adj_vertices = _mesh.vv_range(vh_middle).to_vector();
      auto heh_0 = _mesh.find_halfedge(vh_middle,adj_vertices[0]);
      auto heh_1 = _mesh.find_halfedge(vh_middle,adj_vertices[1]);
      auto heh_2 = _mesh.find_halfedge(vh_middle, adj_vertices[2]);

      int adj_vertices_size = adj_vertices.size();
      std::vector<int> hehs;
      
        if(adj_vertices_size == 1){
          std::cerr<<"Something is topologically wrong, one vertex has only one adjacent vertex"<<std::endl;
        }
        else if(adj_vertices_size == 2){
          // Must both belong to the same face, so no strip
          continue;
        }
        else if(adj_vertices_size == 3){
          // The pattern looks like this:        |
          //                                _____|_____
          // with vertices at each endpoint
          // first obtain the halfedges going out from the middle
          for(int i = 0; i < 3; i++){
            auto he = _mesh.find_halfedge(vh_middle,adj_vertices[i]);
            // Two of those halfedges are adjacent to the boundary -> find out which ones
            // Those are the two opposite ones
            if(this->getFaceOfHalfEdge(he.idx(),0) == -1 || getFaceOfHalfEdge(he.idx(),1) == -1){
              hehs.push_back(he.idx());
            }
          }
          // Now, we have the two opposite halfedges
          // They are all pointing outward from the middle
          auto vh_0 = _mesh.to_vertex_handle(HalfedgeHandle(hehs[0]));
          auto vh_2 = _mesh.to_vertex_handle(HalfedgeHandle(hehs[1]));
          vertexStrips.push_back(std::make_tuple(vh_0.idx(),vh_middle.idx(),vh_2.idx()));
        }
        else if(adj_vertices_size == 4){
          int counter = 0;
          std::vector<int> vh_even;
          std::vector<int> vh_odd;
          // iterate through neighbouring vertices clockwise
          for(const auto& vh_it : _mesh.vv_cw_range(vh_middle)){
            if((counter % 2) == 0){
              vh_even.push_back(vh_it.idx());
            }
            else{
              vh_odd.push_back(vh_it.idx());
            }
            counter ++;
          }
          vertexStrips.push_back(std::make_tuple(vh_even[0],vh_middle.idx(),vh_even[1]));
          vertexStrips.push_back(std::make_tuple(vh_odd[0],vh_middle.idx(),vh_odd[1]));
        }
      }
    return vertexStrips;
  }

// this method works, tested -> only returns halfedge strips that lie in the interior
// since we dont have rulings computed at the boundary, halfedge strips at the boundary cannot contribute to fairness energy
  StripType getHalfEdgeStrips() const
  {
    // First, every HalfEdge strip must also contain a vertex triple.
    // so we can iterate over all vertex triples and check if we can extend them in one of two directions
    // to form a halfedge strip

    auto vertexStrips = getVertexStrips();
    StripType halfedgeStrips;

    for(int i = 0; i < vertexStrips.size(); i++){
      int node_0 = std::get<0>(vertexStrips[i]);
      int node_1 = std::get<1>(vertexStrips[i]);
      int node_2 = std::get<2>(vertexStrips[i]);

      // don't want duplication of tuples. One vertex triple will have a larger index than the other
      for(int j = 0; j < i; j++){
        int node_0_ = std::get<0>(vertexStrips[j]);
        int node_1_ = std::get<1>(vertexStrips[j]);
        int node_2_ = std::get<2>(vertexStrips[j]);

       if(node_0 == node_1_){
          // Still need to check that both vertex strips point in the same direction
          // This can be done by checking that the halfedge from node_0_ to node_0
          // and the halfedge from node_0 to node_1 dont belong to a common face

          auto he_1 = _mesh.find_halfedge(_mesh.vertex_handle(node_0_),_mesh.vertex_handle(node_1_));
          auto he_2 = _mesh.find_halfedge(_mesh.vertex_handle(node_1_),_mesh.vertex_handle(node_1));

          // adjacent faces of halfedge 1
          int face_1_1 = getFaceOfHalfEdge(he_1.idx(),0);
          int face_1_2 = getFaceOfHalfEdge(he_1.idx(),1);

          int face_2_1 = getFaceOfHalfEdge(he_2.idx(),0);
          int face_2_2 = getFaceOfHalfEdge(he_2.idx(),1);

          // if the halfedges have a common face, they dont point in the same direction
          // if they both lie at the boundary, i.e. both have neighbour -1, that is acceptable
          bool common_face = 
          (face_1_1 == face_2_1) || 
          (face_1_1 == face_2_2) || 
          (face_1_2 == face_2_1) || 
          (face_1_2 == face_2_2);

          if(!common_face){
              auto he_3 = _mesh.find_halfedge(_mesh.vertex_handle(node_1),_mesh.vertex_handle(node_2));
              int face_3_1 = getFaceOfHalfEdge(he_3.idx(),0);
              int face_3_2 = getFaceOfHalfEdge(he_3.idx(),1);

              // None of the halfedges adjacent to boundary
              if((face_1_1 != -1) && (face_1_2 != -1) && (face_2_1 != -1) && (face_2_2 != -1) && (face_3_1 != -1) && (face_3_2 != -1))
              {
                halfedgeStrips.push_back(std::make_tuple(he_1.idx(),he_2.idx(),he_3.idx()));
              }
          }
       }
       if(node_2 == node_1_){
          auto he_1 = _mesh.find_halfedge(_mesh.vertex_handle(node_0),_mesh.vertex_handle(node_1));
          auto he_2 = _mesh.find_halfedge(_mesh.vertex_handle(node_1),_mesh.vertex_handle(node_1_));
          
          // adjacent faces of halfedge 1
          int face_1_1 = getFaceOfHalfEdge(he_1.idx(),0);
          int face_1_2 = getFaceOfHalfEdge(he_1.idx(),1);

          int face_2_1 = getFaceOfHalfEdge(he_2.idx(),0);
          int face_2_2 = getFaceOfHalfEdge(he_2.idx(),1);

          // if the halfedges have a common face, they dont point in the same direction
          // if they both lie at the boundary, i.e. both have neighbour -1, that is acceptable
          bool common_face = 
          (face_1_1 == face_2_1) || 
          (face_1_1 == face_2_2) || 
          (face_1_2 == face_2_1) || 
          (face_1_2 == face_2_2);

          if(!common_face)
          {
              auto he_3 = _mesh.find_halfedge(_mesh.vertex_handle(node_1_),_mesh.vertex_handle(node_2_));
              int face_3_1 = getFaceOfHalfEdge(he_3.idx(),0);
              int face_3_2 = getFaceOfHalfEdge(he_3.idx(),1);

              // None of the halfedges adjacent to boundary
              if((face_1_1 != -1) && (face_1_2 != -1) && (face_2_1 != -1) && (face_2_2 != -1) && (face_3_1 != -1) && (face_3_2 != -1))
              {
                halfedgeStrips.push_back(std::make_tuple(he_1.idx(),he_2.idx(),he_3.idx()));
              }
          }
       }
      }
    }
    return halfedgeStrips;
  }

  std::vector<int> getBdryHalfEdges() const{
    return _bdryHalfEdges;
  }

  size_t getNumBdryHalfEdges() const{
    return _bdryHalfEdges.size();
  }

  std::vector<int> getBdryFaces() const{
    return _bdryFaces;
  }

  size_t getNumBdryFaces() const{
    return _bdryFaces.size();
  }

  const QuadMeshType& getGrid() const {
    return _mesh;
  }

  int getNodeOfQuad( int globalFaceIndex, int localNodeIndex ) const {
    return _nodesOfQuads[globalFaceIndex][localNodeIndex];
  }
  
  int getEdgeOfQuad( int globalFaceIndex, int localEdgeIndex ) const {
    return _edgesOfQuads[globalFaceIndex][localEdgeIndex];
  }

  int getHalfEdgeOfQuad( int globalFaceIndex, int localHalfEdgeIndex ) const {
    return _halfedgesOfQuads[globalFaceIndex][localHalfEdgeIndex];
  }

  int getFaceOfHalfEdge( int globalHalfEdgeIndex, int localIndex ) const {
    return _facesOfHalfedges[globalHalfEdgeIndex][localIndex];
  }

  int getNeighbourOfQuad( int globalFaceIndex, int localNodeIndex ) const {
    return _neighboursOfQuads[globalFaceIndex][localNodeIndex];
  }

  int getAdjacentNodeOfEdge( int globalEdgeIndex, int localIndex ) const{
    return _nodesOfEdges[globalEdgeIndex][localIndex];
  }

  std::vector<Vec2i> getAdjacentNodesOfEdges( ) const{
    return _nodesOfEdges;
  }
  
  int getAdjacentQuadOfEdge( int globalEdgeIndex, int localIndex ) const{
    return _neighboursOfEdges[globalEdgeIndex][localIndex];
  }

  int getNumHalfEdges() const {
    return 2*_numOfEdges;
  }

  int getNumEdges() const {
    return _numOfEdges;
  }

  int getNumFaces() const {
    return _numOfElems;
  }

  int getNumVertices() const {
    return _numOfNodes;
  }

  bool isFixedNode( int i ) const {
    return _isFixedNode[i];
  }
  
  void setFixedNode( int i ) {
    if( !_isFixedNode[i] ){
        _isFixedNode[i] = true;
        _numOfFixedNodes++;
    }
  }
    
  void resetFixedNodes() {
    _isFixedNode.setAll( false );
    _numOfFixedNodes = 0;
  }

  void setFixedNodes( const BitVector& mask ) {
    if( mask.size() != _numOfNodes )
      throw BasicException ( "MeshTopologySaver::setFixedNodes(): mask has wrong size!" );  
    _isFixedNode = mask;
    // determine number of fixed nodes
    _numOfFixedNodes = 0;
    for( int i = 0; i < _numOfNodes; i++ )
        if( _isFixedNode[i] )
            _numOfFixedNodes++;
  } 

  const BitVector& getFixedNodes() const {
    return _isFixedNode;
  }
  
  int getNumFixedNodes() const {
    return _numOfFixedNodes;    
  }
  
  bool isEdgeValid( int i ) const {
    return _isEdgeValid[i];
  }
  
  template<typename IntVectorType>
  void addUnvalidEdges( const IntVectorType& unvalidEdges ) {
    for( int i = 0; i < unvalidEdges.size(); i++ )
      _isEdgeValid.set( unvalidEdges[i], false );
  }
  
  void resetUnvalidEdges() {
    _isEdgeValid.setAll( true );
  }
  
  bool isFaceValid( int i ) const {
    return _isFaceValid[i];
  }
  
  template<typename IntVectorType>
  void addUnvalidFaces( const IntVectorType& unvalidFaces ) {
    for( int i = 0; i < unvalidFaces.size(); i++ )
      _isFaceValid.set( unvalidFaces[i], false );
  }
  
  void resetUnvalidFaces() {
    _isFaceValid.setAll( true );
  }

  // only considers one connected component of boundary (containing vertexIndex)!
  void fillPartialBoundaryMask( std::vector<int>& mask, int vertexIndex = -1 ) const {
    mask.clear();
    // get vertex handle at boundary
    MyMesh::VertexHandle vh = MyMesh::InvalidVertexHandle;
    if( vertexIndex < 0 ){      
      for (MyMesh::ConstVertexIter vIter = _mesh.vertices_begin(); vIter != _mesh.vertices_end(); ++vIter)
          if( _mesh.is_boundary ( *vIter ) ){
              vh = *vIter;
              break;
          }      
    }
    else{
      vh = _mesh.vertex_handle(vertexIndex);  
    }
    
    if( vh == MyMesh::InvalidVertexHandle )
        return;
      
    // get hafledge handle at boundary
    MyMesh::ConstVertexIHalfedgeIter hIter ( _mesh, vh );
    for (; !_mesh.is_boundary ( *hIter ); ++hIter ) ;
    
    // run along boundary
    HEH startHEH = *hIter;
    mask.push_back(  _mesh.to_vertex_handle ( startHEH ).idx() );
    HEH bdryHEH = _mesh.next_halfedge_handle ( startHEH );
    while( startHEH != bdryHEH ){
      mask.push_back(  _mesh.to_vertex_handle( bdryHEH ).idx() );
      bdryHEH = _mesh.next_halfedge_handle ( bdryHEH );
    }

  }

  Eigen::MatrixXi getConnectivityMatrix() const {
    Eigen::MatrixXi meshF(getNumFaces(), 4);
    for (int faceIdx = 0; faceIdx < getNumFaces(); faceIdx++) {
      meshF(faceIdx, 0) = getNodeOfQuad(faceIdx, 0);
      meshF(faceIdx, 1) = getNodeOfQuad(faceIdx, 1);
      meshF(faceIdx, 2) = getNodeOfQuad(faceIdx, 2);
      meshF(faceIdx, 3) = getNodeOfQuad(faceIdx, 3); 
    }
    return meshF;
  }


//=================================================================================
//=================================================================================

// geom = (x,y,z) = (x_1, ..., x_n, y_1, ..., y_n, z_1, ..., z_n)
template<typename VectorType>
static void getGeometry( const MyMesh& Mesh, VectorType& geom ){
    
    geom.resize( 3 * Mesh.n_vertices() );
    
    // get pointer to first argument in vertex coordinates list
    typedef MyMesh::Point Point;
    const Point* iter = Mesh.points();
    
    for( uint i = 0; i < Mesh.n_vertices(); i++, iter++ )
        for( int j = 0; j < 3; j++ )
            geom[i*3+j] = (*iter)[j];            
}

// geom = (x,y,z) = (x_1, ..., x_n, y_1, ..., y_n, z_1, ..., z_n)
template<typename VectorType>
static void setGeometry( MyMesh& Mesh, const VectorType& geom ){
    
    if( 3 * Mesh.n_vertices() != (uint)geom.size() )
        throw BasicException( "setGeometry: sizes don't match!" );

    for( MyMesh::VertexIter v_it = Mesh.vertices_begin(); v_it!=Mesh.vertices_end(); ++v_it) {
      MyMesh::Point p;
      for(int j = 0; j < 3; ++j) 
          p[j] = geom[ (v_it->idx())*3 +  j];
      Mesh.set_point( *v_it, p );
    }        
}


};