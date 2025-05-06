#ifndef RefSpectraModel_h
#define RefSpectraModel_h

#include <Wt/WModelIndex>
#include <Wt/WAbstractItemModel>

#include <memory>
#include <string>
#include <vector>


#include <boost/any.hpp>


class RefSpectraModel : public Wt::WAbstractItemModel
{
public:
  RefSpectraModel( Wt::WObject *parent = nullptr );
  virtual ~RefSpectraModel();

  // Required WAbstractItemModel overrides
  virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const override;
  virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const override;
  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const override;
  virtual Wt::WModelIndex index( int row, int column, const Wt::WModelIndex &parent = Wt::WModelIndex() ) const override;
  virtual boost::any data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const override;
  virtual boost::any headerData( int section, Wt::Orientation orientation = Wt::Orientation::Horizontal, int role = Wt::DisplayRole ) const override;

  virtual void sort(int column, Wt::SortOrder order = Wt::SortOrder::AscendingOrder );

  // Custom methods for managing the filesystem data

  /** Add a base directory to the model.  Returns the index of the new root node. */
  Wt::WModelIndex addBaseDirectory( const std::string &path );
  /** Remove a base directory from the model. */
  void removeBaseDirectory( const std::string &path );
  void refresh();
  std::string getFilePath( const Wt::WModelIndex &index ) const;
  std::string getDisplayName( const Wt::WModelIndex &index ) const;
  bool isDirectory( const Wt::WModelIndex &index ) const;
private:
  struct Node {
    std::string name;
    std::string fullPath;
    bool isDirectory;
    std::vector<std::unique_ptr<Node>> children;
    Node *parent;
    
    Node( const std::string &n, const std::string &p, bool isDir, Node *par = nullptr )
      : name(n), fullPath(p), isDirectory(isDir), parent(par) {}
    
    void sort_children( Wt::SortOrder order );
  };

  std::vector<std::unique_ptr<Node>> m_rootNodes;
  
  void populateNode( Node *node );
  Node* getNode( const Wt::WModelIndex &index ) const;
  Wt::WModelIndex createIndex( int row, int column, Node *node ) const;
};

#endif // RefSpectraModel_h 