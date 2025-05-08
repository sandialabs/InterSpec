#include "InterSpec/RefSpectraModel.h"
#include <Wt/WModelIndex>

#include <algorithm>
#include <iostream>
#include <filesystem>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

namespace fs = std::filesystem;


RefSpectraModel::RefSpectraModel( Wt::WObject *parent )
  : Wt::WAbstractItemModel( parent )
{
}


RefSpectraModel::~RefSpectraModel()
{
}


int RefSpectraModel::columnCount( const Wt::WModelIndex &parent ) const
{
  return 1;  // We only need one column for the file/directory name
}//int columnCount( const Wt::WModelIndex &parent ) const


int RefSpectraModel::rowCount( const Wt::WModelIndex &parent ) const
{
  if( !parent.isValid() )
    return static_cast<int>( m_rootNodes.size() );
  
  Node *node = getNode( parent );
  if( !node )
    return 0;
  
  return static_cast<int>( node->children.size() );
}//int rowCount( const Wt::WModelIndex &parent ) const


Wt::WModelIndex RefSpectraModel::parent( const Wt::WModelIndex &index ) const
{
  if( !index.isValid() )
    return Wt::WModelIndex();
  
  Node *node = getNode( index );
  if( !node || !node->parent )
    return Wt::WModelIndex();
  
  // Find the row of this node in its parent's children
  const auto &siblings = node->parent->parent ? node->parent->parent->children : m_rootNodes;
  auto it = std::find_if( siblings.begin(), siblings.end(),
                         [node]( const auto &n ) { return n.get() == node->parent; } );
  
  if( it == siblings.end() )
    return Wt::WModelIndex();
  
  return createIndex( static_cast<int>(std::distance(begin(siblings), it)), 0, node->parent );
}//Wt::WModelIndex parent( const Wt::WModelIndex &index ) const


Wt::WModelIndex RefSpectraModel::index( int row, int column, const Wt::WModelIndex &parent ) const
{
  if( row < 0 || column != 0 )
    return Wt::WModelIndex();
  
  if( !parent.isValid() )
  {
    if( row >= static_cast<int>(m_rootNodes.size()) )
      return Wt::WModelIndex();
    return createIndex( row, column, m_rootNodes[row].get() );
  }
  
  Node *parentNode = getNode( parent );
  if( !parentNode || row >= static_cast<int>(parentNode->children.size()) )
    return Wt::WModelIndex();
  
  return createIndex( row, column, parentNode->children[row].get() );
}//Wt::WModelIndex index(...) const


boost::any RefSpectraModel::data( const Wt::WModelIndex &index, int role ) const
{
  if( !index.isValid() )
    return boost::any();
  
  Node *node = getNode( index );
  if( !node )
    return boost::any();
  
  if( role == Wt::DisplayRole )
    return boost::any( node->name );
  
  return boost::any();
}//boost::any data(...) const


boost::any RefSpectraModel::headerData( int section, Wt::Orientation orientation, int role ) const
{
  if( orientation == Wt::Orientation::Horizontal && role == Wt::DisplayRole ) {
    return boost::any( Wt::WString::tr("rs-dialog-title") );
  }
  
  return boost::any();
}//boost::any headerData(...) const


void RefSpectraModel::Node::sort_children( Wt::SortOrder order )
{
  std::sort( begin(children), end(children), [order]( const std::unique_ptr<Node> &a, const std::unique_ptr<Node> &b ) {
    return (order == Wt::SortOrder::AscendingOrder) ? (a->name < b->name) : (b->name < a->name);
  } );

  for( const std::unique_ptr<Node> &child : children )
    child->sort_children( order );
}//void sort( Wt::SortOrder order );


void RefSpectraModel::sort( int column, Wt::SortOrder order )
{
  layoutAboutToBeChanged().emit();

  std::sort( m_rootNodes.begin(), m_rootNodes.end(), [order]( const std::unique_ptr<Node> &a, const std::unique_ptr<Node> &b ) {
    return (order == Wt::SortOrder::AscendingOrder) ? (a->name < b->name) : (b->name < a->name);
  } );

  for( const std::unique_ptr<Node> &node : m_rootNodes )
    node->sort_children( order );

  layoutChanged().emit();
}//void sort(...)


Wt::WModelIndex RefSpectraModel::addBaseDirectory( const std::string &path )
{
  try 
  {
    fs::path fsPath( path );
    if( !fs::exists( fsPath ) || !fs::is_directory( fsPath ) )
      return Wt::WModelIndex();

    std::string name = fsPath.filename().string();  
    SpecUtils::ireplace_all( name, "_", " " ); // Replace underscores with spaces
    
    std::unique_ptr<Node> node = std::make_unique<Node>( name, path, true );
    
    //We wont populate child nodes until we actually want to access their data, incase the
    //  directory is some really tons of files.
    //populateNode( node.get() );
    
    // Notify that we've added a new root node
    beginInsertRows( Wt::WModelIndex(), static_cast<int>(m_rootNodes.size()), static_cast<int>(m_rootNodes.size()) );
    m_rootNodes.push_back( std::move(node) );
    endInsertRows();

    return createIndex( static_cast<int>(m_rootNodes.size() - 1), 0, m_rootNodes.back().get() );
  }catch( const std::exception &e )
  {
    std::cerr << "Error adding base directory " << path << ": " << e.what() << std::endl;
  }

  return Wt::WModelIndex();
}//Wt::WModelIndex addBaseDirectory( const std::string &path )


void RefSpectraModel::removeBaseDirectory( const std::string &path )
{
  auto it = std::find_if( m_rootNodes.begin(), m_rootNodes.end(),
                         [&path]( const auto &node ) -> bool { 
    try
    { 
      const std::string normed_node_path = SpecUtils::lexically_normalize_path( node->fullPath );
      const std::string normed_path = SpecUtils::lexically_normalize_path( path );
      return normed_node_path == normed_path; 
    }catch( const std::exception &e )
    {
      std::cerr << "RefSpectraModel::removeBaseDirectory() - Invalid directory: '" << path << "'" << std::endl;
      return node->fullPath == path; 
    }
  } );
  
  if( it != m_rootNodes.end() ) {
    int row = static_cast<int>( std::distance( begin(m_rootNodes), it ) );
    beginRemoveRows( Wt::WModelIndex(), row, row );
    m_rootNodes.erase( it );
    endRemoveRows();
  }
}//void removeBaseDirectory( const std::string &path )


void RefSpectraModel::refresh()
{
  // Store the current paths
  std::vector<std::string> paths;
  for( const auto &node : m_rootNodes ) {
    paths.push_back( node->fullPath );
  }
  
  // Clear the model
  const bool nonEmpty = !m_rootNodes.empty(); 
  if( nonEmpty )
  {
    beginRemoveRows( Wt::WModelIndex(), 0, static_cast<int>(m_rootNodes.size()) - 1 );
    m_rootNodes.clear();
    endRemoveRows();
  }
  
  // Re-add all directories
  for( const auto &path : paths )
    addBaseDirectory( path );
}//void refresh()


std::string RefSpectraModel::getFilePath( const Wt::WModelIndex &index ) const
{
  Node *node = getNode( index );
  return node ? node->fullPath : "";
}//std::string getFilePath( const Wt::WModelIndex &index ) const


std::string RefSpectraModel::getDisplayName( const Wt::WModelIndex &index ) const
{
  Node *node = getNode( index );
  return node ? node->name : "";
}//std::string getDisplayName( const Wt::WModelIndex &index ) const


bool RefSpectraModel::isDirectory( const Wt::WModelIndex &index ) const
{
  Node *node = getNode( index );
  return node ? node->isDirectory : false;
}//bool isDirectory( const Wt::WModelIndex &index ) const


Wt::WModelIndex RefSpectraModel::indexForDisplayName( const std::string &displayName, const Wt::WModelIndex &parent ) const
{
  if( !parent.isValid() )
  {
    for( size_t i = 0; i < m_rootNodes.size(); ++i )
    {
      if( m_rootNodes[i]->name == displayName )
        return createIndex( static_cast<int>(i), 0, m_rootNodes[i].get() );
    }

    return Wt::WModelIndex();
  }//if( !parent.isValid() )

  Node *node = getNode( parent );
  if( !node )
    return Wt::WModelIndex();
  
  for( size_t i = 0; i < node->children.size(); ++i )
  {
    if( node->children[i]->name == displayName )
      return createIndex( static_cast<int>(i), 0, node->children[i].get() );
  }

  return Wt::WModelIndex();
}//Wt::WModelIndex indexForDisplayName(...) const


void RefSpectraModel::populateNode( Node *node ) const
{
  try 
  {
    for( const auto &entry : fs::directory_iterator(node->fullPath) ) {
      const fs::path entryPath = entry.path();
      std::string entryName = entryPath.filename().string();
      const std::string entryPathStr = entryPath.string();
      const bool isDir = entry.is_directory();
      
      if( entryName.empty() || entryName.front() == '.' )
        continue;

      if( !isDir && (SpecUtils::iequals_ascii( entryName, "readme.txt" ) 
                     || SpecUtils::iequals_ascii( entryName, "readme.xml" ) ) )
      {
        continue;
      }

      // If not a directory, remove extension and replace all underscores with spaces.
      if( !isDir )
      {
        const std::string ext = SpecUtils::file_extension( entryName );
        if( !ext.empty() )
          entryName = entryName.substr( 0, entryName.size() - ext.size() );
      }//if( !isDir )

      SpecUtils::ireplace_all( entryName, "_", " " ); // Replace underscores with spaces

      std::unique_ptr<Node> child = std::make_unique<Node>( entryName, entryPathStr, isDir, node );
      
      //We will not populate child directories yet - we'll wait until that data is actually
      //  requested - this saves us recursively iterating over tons of files we wont actually use.
      //if( isDir )
      //  populateNode( child.get() );
      
      node->children.push_back( std::move(child) );
    }
  } catch( const std::exception &e ) {
    std::cerr << "Error reading directory " << node->fullPath << ": " << e.what() << std::endl;
  }
}//void RefSpectraModel::populateNode( Node *node )


RefSpectraModel::Node *RefSpectraModel::getNode( const Wt::WModelIndex &index ) const
{
  Node *node = static_cast<Node *>(index.internalPointer());
  
  // We dont populate directory child nodes until we need the data, so we'll check for this
  //  before returning the node, and populate if necessary.
  //  We only access nodes through this function.
  if( node && node->isDirectory && node->children.empty() )
    populateNode( node );
  
  return node;
}//RefSpectraModel::Node *getNode(...) const


Wt::WModelIndex RefSpectraModel::createIndex( int row, int column, Node *node ) const
{
  return WAbstractItemModel::createIndex( row, column, (void *)node );
}//Wt::WModelIndex createIndex(...) const
