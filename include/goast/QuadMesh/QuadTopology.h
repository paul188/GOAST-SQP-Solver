#pragma once

//== INCLUDES =================================================================
#include "QuadMesh.h"
#include "OpenMeshIncludes.h"
#include "Auxiliary.h"

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
class QuadMeshTopologySaver {

protected:  
  typedef MyMesh QuadMeshType;
  typedef VectorT<signed int,8> Vec8i;
  static const int IndexNotSet = -1;  
  const QuadMeshType& _mesh;

  const int _numOfElems;
  const int _numOfNodes;
  const int _numOfEdges;

  // saves global indices of nodes belonging to each quad
  std::vector< Vec4i > _nodesOfQuads;
  // saves global indices of neighbouring quads for each quad (these are -1 at the boundary )
  std::vector< Vec4i > _neighboursOfQuads;
  // saves global indices of edges belonging to each triangle
  std::vector< Vec4i > _edgesOfQuads;
  // saves global indices of neighbouring quads for each edge (these are -1 at the boundary )
  std::vector< Vec2i > _neighboursOfEdges;  
  // saves global indices of nodes belonging to each edge
  std::vector< Vec2i > _nodesOfEdges;
  // saves global indices of nodes opposite of each edge (two of them are -1 for boundary edges)
  std::vector< Vec4i > _oppositeNodesOfEdges;
   
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
      
      // indices of neighbours of current face
      // ordering convention: i'th neighboring Quad is opposite of halfedge starting from vertex with local index i
      locIdx = 2;
      for ( MyMesh::ConstFaceHalfedgeIter cfh_it = mesh.cfh_iter(*f_it); cfh_it.is_valid(); ++cfh_it){
        _neighboursOfQuads[f_it->idx()][locIdx%4] = mesh.is_boundary(*cfh_it) ? IndexNotSet : mesh.face_handle( mesh.opposite_halfedge_handle(*cfh_it) ).idx();
        locIdx++;
      }
  
      //
      locIdx = 1;
      for (TriMesh::ConstFaceEdgeIter cfe_it = mesh.cfe_iter(*f_it); cfe_it.is_valid(); ++cfe_it) {
        _edgesOfQuads[f_it->idx()][locIdx] = cfe_it->idx();
        locIdx = (locIdx + 1) % 3;
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

  }

   const QuadMeshType& getGrid() const {
    return _mesh;
  }

  int getNodeOfTriangle( int globalFaceIndex, int localNodeIndex ) const {
    return _nodesOfQuads[globalFaceIndex][localNodeIndex];
  }
  
  int getEdgeOfTriangle( int globalFaceIndex, int localEdgeIndex ) const {
    return _edgesOfQuads[globalFaceIndex][localEdgeIndex];
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
      meshF(faceIdx, 0) = getNodeOfTriangle(faceIdx, 0);
      meshF(faceIdx, 1) = getNodeOfTriangle(faceIdx, 1);
      meshF(faceIdx, 2) = getNodeOfTriangle(faceIdx, 2);
      meshF(faceIdx, 3) = getNodeOfTriangle(faceIdx, 3); 
    }
    return meshF;
  }


//=================================================================================
//=================================================================================

// points of the geometry are vectors in column of the matrix
template<typename FullMatrixType>
void getGeometryMat( const MyMesh& Mesh, FullMatrixType& geom ){
    
    geom.resize( 3, Mesh.n_vertices() );
    
    // get pointer to first argument in vertex coordinates list
    typedef MyMesh::Point Point;
    const Point* iter = Mesh.points();
    
    for( uint i = 0; i < Mesh.n_vertices(); i++, iter++ )
        for( int j = 0; j < 3; j++ )
            geom[j, i] = (*iter)[j];            
}


// geom = (x,y,z) = (x_1, ..., x_n, y_1, ..., y_n, z_1, ..., z_n)
template<typename VectorType>
void getGeometry( const MyMesh& Mesh, VectorType& geom ){
    
    geom.resize( 3 * Mesh.n_vertices() );
    
    // get pointer to first argument in vertex coordinates list
    typedef MyMesh::Point Point;
    const Point* iter = Mesh.points();
    
    for( uint i = 0; i < Mesh.n_vertices(); i++, iter++ )
        for( int j = 0; j < 3; j++ )
            geom[j*Mesh.n_vertices()+i] = (*iter)[j];            
}

// geom = (x,y,z) = (x_1, ..., x_n, y_1, ..., y_n, z_1, ..., z_n)
template<typename VectorType>
void setGeometry( MyMesh& Mesh, const VectorType& geom ){
    
    if( 3 * Mesh.n_vertices() != (uint)geom.size() )
        throw BasicException( "setGeometry: sizes don't match!" );

    for( MyMesh::VertexIter v_it = Mesh.vertices_begin(); v_it!=Mesh.vertices_end(); ++v_it) {
      MyMesh::Point p;
      for(int j = 0; j < 3; ++j) 
          p[j] = geom[ j*Mesh.n_vertices() + v_it->idx() ];
      Mesh.set_point( *v_it, p );
    }        
}


};