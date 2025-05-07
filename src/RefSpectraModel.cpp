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
}

int RefSpectraModel::rowCount( const Wt::WModelIndex &parent ) const
{
  if( !parent.isValid() ) {
    return m_rootNodes.size();
  }
  
  Node *node = getNode( parent );
  if( !node ) {
    return 0;
  }
  
  return node->children.size();
}

Wt::WModelIndex RefSpectraModel::parent( const Wt::WModelIndex &index ) const
{
  if( !index.isValid() ) {
    return Wt::WModelIndex();
  }
  
  Node *node = getNode( index );
  if( !node || !node->parent ) {
    return Wt::WModelIndex();
  }
  
  // Find the row of this node in its parent's children
  const auto &siblings = node->parent->parent ? node->parent->parent->children : m_rootNodes;
  auto it = std::find_if( siblings.begin(), siblings.end(),
                         [node]( const auto &n ) { return n.get() == node->parent; } );
  
  if( it == siblings.end() ) {
    return Wt::WModelIndex();
  }
  
  return createIndex( std::distance( siblings.begin(), it ), 0, node->parent );
}

Wt::WModelIndex RefSpectraModel::index( int row, int column, const Wt::WModelIndex &parent ) const
{
  if( row < 0 || column != 0 ) {
    return Wt::WModelIndex();
  }
  
  if( !parent.isValid() ) {
    if( row >= static_cast<int>(m_rootNodes.size()) ) {
      return Wt::WModelIndex();
    }
    return createIndex( row, column, m_rootNodes[row].get() );
  }
  
  Node *parentNode = getNode( parent );
  if( !parentNode || row >= static_cast<int>(parentNode->children.size()) ) {
    return Wt::WModelIndex();
  }
  
  return createIndex( row, column, parentNode->children[row].get() );
}

boost::any RefSpectraModel::data( const Wt::WModelIndex &index, int role ) const
{
  if( !index.isValid() ) {
    return boost::any();
  }
  
  Node *node = getNode( index );
  if( !node ) {
    return boost::any();
  }
  
  if( role == Wt::DisplayRole ) {
    return boost::any( node->name );
  }
  
  return boost::any();
}

boost::any RefSpectraModel::headerData( int section, Wt::Orientation orientation, int role ) const
{
  if( orientation == Wt::Orientation::Horizontal && role == Wt::DisplayRole ) {
    return boost::any( Wt::WString::tr("rs-dialog-title") );
  }
  
  return boost::any();
}

void RefSpectraModel::Node::sort_children( Wt::SortOrder order )
{
  std::sort( begin(children), end(children), [order]( const std::unique_ptr<Node> &a, const std::unique_ptr<Node> &b ) {
    return (order == Wt::SortOrder::AscendingOrder) ? (a->name < b->name) : (b->name < a->name);
  } );

  for( const std::unique_ptr<Node> &child : children )
    child->sort_children( order );
}//void sort( Wt::SortOrder order );



void RefSpectraModel::sort(int column, Wt::SortOrder order )
{
  layoutAboutToBeChanged().emit();

  std::sort( m_rootNodes.begin(), m_rootNodes.end(), [order]( const std::unique_ptr<Node> &a, const std::unique_ptr<Node> &b ) {
    return (order == Wt::SortOrder::AscendingOrder) ? (a->name < b->name) : (b->name < a->name);
  } );

  for( const std::unique_ptr<Node> &node : m_rootNodes )
    node->sort_children( order );

  layoutChanged().emit();
}

Wt::WModelIndex RefSpectraModel::addBaseDirectory( const std::string &path )
{
  try 
  {
    fs::path fsPath( path );
    if( !fs::exists( fsPath ) || !fs::is_directory( fsPath ) ) {
      return Wt::WModelIndex();
    }

    std::string name = fsPath.filename().string();  
    SpecUtils::ireplace_all( name, "_", " " ); // Replace underscores with spaces
    
    std::unique_ptr<Node> node = std::make_unique<Node>( name, path, true );
    populateNode( node.get() );
    
    // Notify that we've added a new root node
    beginInsertRows( Wt::WModelIndex(), m_rootNodes.size(), m_rootNodes.size() );
    m_rootNodes.push_back( std::move(node) );
    endInsertRows();

    return createIndex( m_rootNodes.size() - 1, 0, m_rootNodes.back().get() );
  } catch( const std::exception &e ) {
    std::cerr << "Error adding base directory " << path << ": " << e.what() << std::endl;
  }

  return Wt::WModelIndex();
}

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
    int row = std::distance( m_rootNodes.begin(), it );
    beginRemoveRows( Wt::WModelIndex(), row, row );
    m_rootNodes.erase( it );
    endRemoveRows();
  }
}

void RefSpectraModel::refresh()
{
  // Store the current paths
  std::vector<std::string> paths;
  for( const auto &node : m_rootNodes ) {
    paths.push_back( node->fullPath );
  }
  
  // Clear the model
  const bool nonEmpty = !m_rootNodes.empty(); 
  if( nonEmpty ) {
    beginRemoveRows( Wt::WModelIndex(), 0, m_rootNodes.size() - 1 );
    m_rootNodes.clear();
    endRemoveRows();
  }
  
  // Re-add all directories
  for( const auto &path : paths ) {
    addBaseDirectory( path );
  }
}

std::string RefSpectraModel::getFilePath( const Wt::WModelIndex &index ) const
{
  Node *node = getNode( index );
  return node ? node->fullPath : "";
}

std::string RefSpectraModel::getDisplayName( const Wt::WModelIndex &index ) const
{
  Node *node = getNode( index );
  return node ? node->name : "";
}

bool RefSpectraModel::isDirectory( const Wt::WModelIndex &index ) const
{
  Node *node = getNode( index );
  return node ? node->isDirectory : false;
}

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
}

void RefSpectraModel::populateNode( Node *node )
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
      
      if( isDir )
        populateNode( child.get() );
      
      node->children.push_back( std::move(child) );
    }
  } catch( const std::exception &e ) {
    std::cerr << "Error reading directory " << node->fullPath << ": " << e.what() << std::endl;
  }
}

RefSpectraModel::Node* RefSpectraModel::getNode( const Wt::WModelIndex &index ) const
{
  return static_cast<Node*>(index.internalPointer());
}

Wt::WModelIndex RefSpectraModel::createIndex( int row, int column, Node *node ) const
{
  return WAbstractItemModel::createIndex( row, column, (void *)node );
} 