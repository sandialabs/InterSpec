/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov.
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "InterSpec_config.h"

#include <vector>
#include <algorithm>

#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WComboBox>
#include <Wt/WGroupBox>
#include <Wt/WIOService>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WAbstractItemModel>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/MakeDrfChart.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RowStretchTreeView.h"


using namespace std;
using namespace Wt;


class FwhmPeaksModel : public  Wt::WAbstractItemModel
{
public:
  enum class Column : int { Energy, Fwhm, FWhmUncert, UserOrAuto, UseForFit, NumColumn };
  
protected:
  Column m_sort_col;
  SortOrder m_sort_order = AscendingOrder;
  std::vector<MakeFwhmForDrf::TableRow> m_rows;
  
  
public:
  FwhmPeaksModel( Wt::WObject *parent = nullptr )
  : Wt::WAbstractItemModel( parent ),
  m_sort_col( Column::Energy )
  {
  }
  
  virtual ~FwhmPeaksModel() override {}
  
  const std::vector<MakeFwhmForDrf::TableRow> &rowData() const
  {
    return m_rows;
  };
  
  
  void setRowData( const std::vector<MakeFwhmForDrf::TableRow> &data )
  {
    if( !m_rows.empty() )
    {
      beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_rows.size()-1) );
      m_rows.clear();
      endRemoveRows();
    }
    
    if( !data.empty() )
    {
      beginInsertRows( WModelIndex(), 0, static_cast<int>(data.size())-1 );
      
      m_rows = data;
      sortImp( m_sort_col, m_sort_order, m_rows );
      
      endInsertRows();
    }//if( !data.empty() )
  }//void setRowData( const std::vector<MakeFwhmForDrf::TableRow> &data )
  
  
  virtual WModelIndex index( int row, int column,
                                const Wt::WModelIndex &parent = Wt::WModelIndex() ) const override
  {
    if( parent.isValid()
       || (column < 0) || (column >= static_cast<int>(Column::NumColumn))
       || (row < 0) || (row >= static_cast<int>(m_rows.size())) )
      return WModelIndex();
    
    return createIndex( row, column, nullptr );
  }
  
  
  virtual WModelIndex parent( const Wt::WModelIndex &index ) const override
  {
    return WModelIndex();
  }
  
  
  virtual int rowCount( const Wt::WModelIndex & parent = Wt::WModelIndex() ) const override
  {
    return parent.isValid() ? 0 : static_cast<int>(m_rows.size());
  }
  
  
  virtual int columnCount( const Wt::WModelIndex & parent = Wt::WModelIndex() ) const override
  {
    return parent.isValid() ? 0 : static_cast<int>(Column::NumColumn);
  }
  
  
  virtual WFlags<ItemFlag> flags( const WModelIndex &index ) const override
  {
    if( index.column() == static_cast<int>(Column::UseForFit) )
      return WFlags<ItemFlag>(ItemFlag::ItemIsUserCheckable);
    
    return WFlags<ItemFlag>(0);
  }
  
  virtual boost::any headerData( int section, Orientation orientation = Horizontal, int role = DisplayRole) const override
  {
    if( (orientation != Horizontal) || (role != DisplayRole) )
      return boost::any();
    
    switch( Column(section) )
    {
      case Column::Energy:     return WString::tr("Energy (keV)");
      case Column::Fwhm:       return WString::tr("FWHM");
      case Column::FWhmUncert: return WString::tr("fpm-fwhm-uncert");
      case Column::UserOrAuto: return WString::tr("fpm-peak-source");
      case Column::UseForFit:  return WString::tr("Use");
      case Column::NumColumn: break;
    }
    assert( 0 );
    return boost::any();
  }//headerData(...)

  
  virtual boost::any data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const override
  {
    if( !index.isValid()
       || (index.column() < 0) || (index.column() >= static_cast<int>(Column::NumColumn))
       || (index.row() < 0) || (index.row() >= static_cast<int>(m_rows.size())) )
    {
      return boost::any();
    }
    
    if( (Column(index.column()) == Column::UseForFit) && (role != ItemDataRole::CheckStateRole) )
      return boost::any();
    
    if( (Column(index.column()) != Column::UseForFit) && (role != ItemDataRole::DisplayRole) )
      return boost::any();
    
    // Could implement a Wt::ItemDataRole::ToolTipRole
    
    const MakeFwhmForDrf::TableRow &row = m_rows[index.row()];
    
    char buffer[32] = { '\0' };
    switch( Column(index.column()) )
    {
      case Column::Energy:
        snprintf( buffer, sizeof(buffer), "%.2f", row.m_peak->mean() );
        break;
        
      case Column::Fwhm:
        snprintf( buffer, sizeof(buffer), "%.2f", row.m_peak->fwhm() );
        break;
      
      case Column::FWhmUncert:
        snprintf( buffer, sizeof(buffer), "%.2f", 2.35482*row.m_peak->uncertainty(PeakDef::CoefficientType::Sigma) );
        break;
        
      case Column::UserOrAuto:
        return WString::tr( row.m_is_user_peak ? "fpm-user" : "fpm-auto-fit" );
        break;
        
      case Column::UseForFit:
        return boost::any(row.m_use_for_fit);
        
      case Column::NumColumn:
        assert( 0 );
        break;
    }//switch( Column(index.column()) )
    
    return boost::any( WString(buffer) );
  }//data(...)
  
  
  void set_use( const shared_ptr<const PeakDef> &p, const bool use_peak )
  {
    for( size_t i = 0; i < m_rows.size(); ++i )
    {
      if( m_rows[i].m_peak == p )
      {
        if( m_rows[i].m_use_for_fit != use_peak )
        {
          const WModelIndex index = createIndex(static_cast<int>(i),
                                                static_cast<int>(Column::UseForFit), nullptr );
          m_rows[i].m_use_for_fit = use_peak;
          dataChanged().emit( index, index );
        }
        
        return;
      }//if( m_rows[i].m_peak == p )
    }//for( size_t i = 0; i < m_rows.size(); ++i )
    
    assert( 0 );
  }//void set_use( distances[index].second, const bool use_peak )
  
  
  virtual bool setData( const WModelIndex &index, const boost::any &value, int role = EditRole ) override
  {
    if( index.parent().isValid() || (index.row() < 0)
       || (index.row() >= static_cast<int>(m_rows.size()))
       || (index.column() != static_cast<int>(Column::UseForFit))
       || value.empty()
       || (value.type() != typeid(bool)) )
    {
      assert( 0 );
      return false;
    }
    
    m_rows[index.row()].m_use_for_fit = boost::any_cast<bool>(value);
    dataChanged().emit( index, index );
    
    return true;
  }//bool setData( const WModelIndex &index, const boost::any &value, int role = EditRole )
  
  
  static void sortImp( const Column column, const SortOrder order,
                      vector<MakeFwhmForDrf::TableRow> &rows )
  {
    const bool littleToBig = (order == AscendingOrder);
    
    std::sort( begin(rows), end(rows),
              [column, littleToBig]( const MakeFwhmForDrf::TableRow &lhs,
                                    const MakeFwhmForDrf::TableRow &rhs ) -> bool {
      switch( Column(column) )
      {
        case Column::Energy:
          return littleToBig ? lhs.m_peak->mean() < rhs.m_peak->mean() : lhs.m_peak->mean() > rhs.m_peak->mean();
          
        case Column::Fwhm:
          return littleToBig ? lhs.m_peak->fwhm() < rhs.m_peak->fwhm() : lhs.m_peak->fwhm() > rhs.m_peak->fwhm();
        
        case Column::FWhmUncert:
        {
          const double l = lhs.m_peak->uncertainty(PeakDef::CoefficientType::Sigma);
          const double r = rhs.m_peak->uncertainty(PeakDef::CoefficientType::Sigma);
          return littleToBig ? l < r : l > r;
        }
          
        case Column::UserOrAuto:
          return littleToBig ? lhs.m_is_user_peak < rhs.m_is_user_peak : lhs.m_is_user_peak > rhs.m_is_user_peak;
          
        case Column::UseForFit:
          return littleToBig ? lhs.m_use_for_fit < rhs.m_use_for_fit : lhs.m_use_for_fit > rhs.m_use_for_fit;
          
        case Column::NumColumn:
          assert( 0 );
          break;
      }//switch( Column(index.column()) )
      
      return false;
    } );
  }//void sortImp( Column column, SortOrder order )
  
  
  virtual void sort( int column, SortOrder order = AscendingOrder ) override
  {
    layoutAboutToBeChanged().emit();
    
    switch( Column(column) )
    {
      case Column::Energy:     case Column::Fwhm:
      case Column::FWhmUncert: case Column::UserOrAuto:
      case Column::UseForFit:
        break;
        
      case Column::NumColumn:
      default:
        throw runtime_error( "FwhmPeaksModel: invalid sort order" );
    }//switch( Column(column) )
    
    m_sort_order = order;
    m_sort_col = Column(column);
    
    sortImp( m_sort_col, m_sort_order, m_rows );
    
    layoutChanged().emit();
  }//sort(...)
  
  
  void set_peaks( const vector<shared_ptr<const PeakDef>> &user_peaks,
                 const vector<shared_ptr<const PeakDef>> &auto_fit_peaks )
  {
    if( !m_rows.empty() )
    {
      beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_rows.size()-1) );
      m_rows.clear();
      endRemoveRows();
    }
    
    if( user_peaks.empty() && auto_fit_peaks.empty() )
      return;
    
    beginInsertRows( WModelIndex(), 0, static_cast<int>(user_peaks.size() + auto_fit_peaks.size())-1 );
    
    for( const auto &p : user_peaks )
    {
      MakeFwhmForDrf::TableRow r;
      r.m_is_user_peak = true;
      r.m_use_for_fit = p->useForDrfFwhmFit();
      r.m_peak = p;
      m_rows.push_back( std::move(r) );
    }//for( const auto &p : user_peaks )
    
    for( const auto &p : auto_fit_peaks )
    {
      MakeFwhmForDrf::TableRow r;
      r.m_is_user_peak = false;
      r.m_use_for_fit = p->useForDrfFwhmFit();
      r.m_peak = p;
      m_rows.push_back( std::move(r) );
    }//for( const auto &p : user_peaks )
    
    sortImp( m_sort_col, m_sort_order, m_rows );
    
    endInsertRows();
  }//void setPeaks(...)
  
  
  static vector<shared_ptr<const PeakDef>> peaks_for_use( const vector<MakeFwhmForDrf::TableRow> &rows )
  {
    vector<shared_ptr<const PeakDef>> answer;
    for( const auto &row : rows )
    {
      if( row.m_use_for_fit )
        answer.push_back( row.m_peak );
    }
    
    return answer;
  }//vector<shared_ptr<const PeakDef>> peaks_to_use() const
  
  
  vector<shared_ptr<const PeakDef>> peaks_to_use() const
  {
    return peaks_for_use( m_rows );
  }//vector<shared_ptr<const PeakDef>> peaks_to_use() const
};//class FwhmPeaksModel
 

MakeFwhmForDrfWindow::MakeFwhmForDrfWindow( const bool use_auto_fit_peaks_too )
 : AuxWindow( WString::tr("window-title-add-fwhm"),
             (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
              | AuxWindowProperties::SetCloseable
              | AuxWindowProperties::DisableCollapse
              | AuxWindowProperties::EnableResize
              | AuxWindowProperties::IsModal) ),
  m_tool( nullptr )
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  if( !viewer )
    return;
  
  viewer->useMessageResourceBundle( "MakeFwhmForDrf" );
  
  AuxWindow *window = this;
  
  const int ww = viewer->renderedWidth();
  const int wh = viewer->renderedHeight();
  if( ww > 100 && wh > 100 )
  {
    const int width = std::min( (8*ww)/9, 750 );
    const int height = std::min( 700, ((wh < 420) ? wh : (19*wh)/20 ) );
    
    window->resizeWindow( width, height );
    window->setMinimumSize( std::min(width,640), std::min(height,480) );
  }//if( ww > 100 && wh > 100 )
  
  shared_ptr<SpecMeas> foreground = viewer->measurment(SpecUtils::SpectrumType::Foreground );
  shared_ptr<DetectorPeakResponse> drf = foreground ? foreground->detector() : nullptr;
  
  m_tool = new MakeFwhmForDrf( use_auto_fit_peaks_too, viewer, drf );
    
  window->stretcher()->addWidget( m_tool, 0, 0 );
  window->stretcher()->setContentsMargins( 0, 0, 0, 0 );
  
  AuxWindow::addHelpInFooter( window->footer(), "add-fwhm-to-drf" );
    
  WPushButton *closeButton = window->addCloseButtonToFooter( WString::tr("Cancel") );
  closeButton->clicked().connect( window, &AuxWindow::hide );
    
  WPushButton *saveAs = new WPushButton( WString::tr("mffdw-use-fwhm-btn"), window->footer() );
  saveAs->clicked().connect( m_tool, &MakeFwhmForDrf::setToDrf );
  m_tool->validationChanged().connect( boost::bind( &WPushButton::setEnabled, saveAs,
                                                                   boost::placeholders::_1 ) );
  // Maybe
  //m_tool->validationChanged().connect( std::bind([m_tool,saveAs](){
  //  saveAs->setEnabled( m_tool->isValidFwhm() );
  //}) );
  
  saveAs->disable();
  window->show();
    
  window->resizeToFitOnScreen();
  window->centerWindow();
  window->rejectWhenEscapePressed( false );
}//MakeFwhmForDrfWindow constructor

  
MakeFwhmForDrf *MakeFwhmForDrfWindow::tool()
{
  return m_tool;
}
  

MakeFwhmForDrf::MakeFwhmForDrf( const bool auto_fit_peaks,
                               InterSpec *viewer,
               std::shared_ptr<DetectorPeakResponse> drf,
               Wt::WContainerWidget *parent )
 : WContainerWidget( parent ),
  m_interspec( viewer ),
  m_refit_scheduled( false ),
  m_undo_redo_scheduled( false ),
  m_orig_drf( drf ),
  m_user_peaks{},
  m_auto_fit_peaks{},
  m_chart( nullptr ),
  m_fwhmEqnType( nullptr ),
  m_sqrtEqnOrder( nullptr ),
  m_parEdits{},
  m_error( nullptr ),
  m_equation( nullptr ),
  m_table( nullptr ),
  m_model( nullptr ),
  m_validationChanged( this ),
  m_updatedDrf( this )
{
  wApp->useStyleSheet( "InterSpec_resources/MakeFwhmForDrf.css" );
  
  addStyleClass( "MakeFwhmForDrf" );
  
  m_interspec->useMessageResourceBundle( "MakeFwhmForDrf" );
    
  // Using a Wt layout since the chart requires this
  WGridLayout *layout = new WGridLayout( this );
  layout->setContentsMargins( 0, 0, 0, 0 );
    
  // TODO: make chart be full energy range of spectrum
  //       Maybe have a parameter be "fixed"...
  m_chart = new MakeDrfChart();
  m_chart->showEfficiencyPoints( false );
  DrfChartHolder *chartholder = new DrfChartHolder( m_chart, nullptr );
  layout->addWidget( chartholder, layout->rowCount(), 0 );
  chartholder->setHeight( 250 );
  //layout->setRowResizable( 0, true, WLength(250, WLength::Pixel) );
  //layout->setRowStretch( 1, 1 );
    
  m_error = new WText( "", Wt::XHTMLText );
  layout->addWidget( m_error, 1, 0 );
  m_error->setInline( false );
  m_error->addStyleClass( "ErrTxt" );
  m_error->hide();
  m_equation = new WText( "", Wt::XHTMLText );
  layout->addWidget( m_equation, 2, 0 );
  m_equation->setInline( false );
  m_equation->addStyleClass( "EqnTxt" );
  m_equation->hide();
  
  WContainerWidget *lowerDiv = new WContainerWidget();
  layout->addWidget( lowerDiv, 3, 0 );
  layout->setRowStretch( 3, 1 );
    
  WGridLayout *lowerLayout = new WGridLayout( lowerDiv );
  //lowerLayout->setVerticalSpacing( 0 );
  //lowerLayout->setHorizontalSpacing( 0 );
  lowerLayout->setContentsMargins( 0, 0, 0, 0 );
    
  WContainerWidget *optionsDiv = new WContainerWidget();
  lowerLayout->addWidget( optionsDiv, 0, 0 );
  optionsDiv->addStyleClass( "Options" );
    
  WGroupBox *box = new WGroupBox( WString::tr("mffd-eqn-type-label"), optionsDiv );
  box->addStyleClass( "OptDiv" );
  m_fwhmEqnType = new WComboBox( box );
  
  for( auto fcnfrm = DetectorPeakResponse::ResolutionFnctForm(0);
      fcnfrm != DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm;
      fcnfrm = DetectorPeakResponse::ResolutionFnctForm( static_cast<int>(fcnfrm) + 1 ) )
  {
    switch( fcnfrm )
    {
      case DetectorPeakResponse::kGadrasResolutionFcn:
        m_fwhmEqnType->addItem( "GADRAS" );
        break;
          
      case DetectorPeakResponse::kSqrtPolynomial:
        m_fwhmEqnType->addItem( "sqrt(A0 + A1*E + ...)" );
        break;
          
      case DetectorPeakResponse::kSqrtEnergyPlusInverse:
        m_fwhmEqnType->addItem( "sqrt(A0 + A1*E + A2/E)" );
        break;
        
      case DetectorPeakResponse::kConstantPlusSqrtEnergy:
        m_fwhmEqnType->addItem( "A0 + A1*sqrt(E)" );
        break;
          
      case DetectorPeakResponse::kNumResolutionFnctForm:
        break;
    }//switch( fcnfrm )
  }//for( loop over functional forms )
    
  m_fwhmEqnType->setCurrentIndex( DetectorPeakResponse::kSqrtPolynomial );
  m_fwhmEqnType->activated().connect( this, &MakeFwhmForDrf::handleFwhmEqnTypeChange );
  
  box = new WGroupBox( WString::tr("mffd-num-terms-label"), optionsDiv );
  box->addStyleClass( "OptDiv" );
  m_sqrtEqnOrder = new WComboBox( box );
  m_sqrtEqnOrder->addItem( "1" );
  m_sqrtEqnOrder->addItem( "2" );
  m_sqrtEqnOrder->addItem( "3" );
  m_sqrtEqnOrder->addItem( "4" );
  m_sqrtEqnOrder->addItem( "5" );
  m_sqrtEqnOrder->setCurrentIndex( 1 );
  m_sqrtEqnOrder->activated().connect( this, &MakeFwhmForDrf::handleSqrtEqnOrderChange );
  
  WGroupBox *parametersDiv = new WGroupBox( WString::tr("mffd-par-vals-label") );
  //WContainerWidget *parametersDiv = new WContainerWidget();
  lowerLayout->addWidget( parametersDiv, 1, 0, Wt::AlignmentFlag::AlignTop );
  lowerLayout->setRowStretch( 1, 1 );
    
  parametersDiv->addStyleClass( "Parameters" );
    
  for( int i = 0; i < m_sqrtEqnOrder->count(); ++i )
  {
    WContainerWidget *parDiv = new WContainerWidget( parametersDiv );
    parDiv->addStyleClass( "ParDiv" );
    WLabel *label = new WLabel( "A" + std::to_string(i), parDiv );
    label->setInline( false );
    NativeFloatSpinBox *sb = new NativeFloatSpinBox( parDiv );
    label->setBuddy( sb );
    sb->setText( "" );
    sb->setSpinnerHidden( true );
    sb->setWidth( 75 );
    sb->valueChanged().connect( boost::bind(&MakeFwhmForDrf::coefficientManuallyChanged, this, i) );
    sb->disable();
    parDiv->setHidden( (i > m_sqrtEqnOrder->currentIndex()) );
    m_parEdits.push_back( sb );
  }//for( int i = 0; i < m_sqrtEqnOrder->count(); ++i )
    
  m_model = new FwhmPeaksModel( this );
  m_table = new RowStretchTreeView();
  m_table->setRootIsDecorated( false ); //makes the tree look like a table! :)
  m_table->addStyleClass( "PeakTable" );
  lowerLayout->addWidget( m_table, 0, 1, 2, 1 );
  m_table->setModel( m_model );
  lowerLayout->setColumnStretch( 1, 1 );
  
  m_table->setColumnWidth( static_cast<int>(FwhmPeaksModel::Column::Energy), 115 );
  m_table->setColumnWidth( static_cast<int>(FwhmPeaksModel::Column::Fwhm), 85 );
  m_table->setColumnWidth( static_cast<int>(FwhmPeaksModel::Column::FWhmUncert), 125 );
  m_table->setColumnWidth( static_cast<int>(FwhmPeaksModel::Column::UserOrAuto), 105 );
  m_table->setColumnWidth( static_cast<int>(FwhmPeaksModel::Column::UseForFit), 60 );
  
    
  m_model->dataChanged().connect( this, &MakeFwhmForDrf::handleTableDataChange );
  m_model->rowsInserted().connect( this, &MakeFwhmForDrf::handleTableDataChange );
  m_model->rowsRemoved().connect( this, &MakeFwhmForDrf::handleTableDataChange );
  m_model->layoutChanged().connect( this, &MakeFwhmForDrf::handleTableDataChange );
    
  if( auto_fit_peaks )
  {
    startAutomatedPeakSearch();
  }else
  {
    const vector<shared_ptr<const PeakDef>> user_peaks = get_user_peaks();
    const auto dummy_auto_fit_peaks = make_shared<vector<shared_ptr<const PeakDef>>>();
    setPeaksFromAutoSearch( user_peaks, dummy_auto_fit_peaks );
  }//if( auto_fit_peaks )
    
  scheduleUndoRedoStep();  //Fills in initial m_current_state
    
  scheduleRefit();
}//MakeFwhmForDrf( constructor )


MakeFwhmForDrf::~MakeFwhmForDrf()
{
  
}//~MakeFwhmForDrf()


void MakeFwhmForDrf::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  if( m_refit_scheduled )
    doRefitWork();
    
  if( m_undo_redo_scheduled )
    doAddUndoRedoStep();
  
  m_refit_scheduled = false;
  m_undo_redo_scheduled = false;
  
  WContainerWidget::render( flags );
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


vector<shared_ptr<const PeakDef>> MakeFwhmForDrf::get_user_peaks()
{
  vector<shared_ptr<const PeakDef>> user_peaks;
  
  PeakModel *peakModel = m_interspec->peakModel();
  assert( peakModel );
  if( !peakModel )
    return user_peaks; //shouldnt happen
  
  const auto originalPeaks = peakModel->peakVec();
  std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > startingPeaks = peakModel->peaks();
  if( startingPeaks && !startingPeaks->empty() )
    user_peaks.insert( end(user_peaks), begin(*startingPeaks), end(*startingPeaks) );
  
  return user_peaks;
}//get_user_peaks() const


void MakeFwhmForDrf::startAutomatedPeakSearch()
{
  auto dataPtr = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  assert( dataPtr );
  if( !dataPtr )
  {
    passMessage( WString::tr("mffd-err-no-foreground"), 1 );
    return;
  }//if( !dataPtr )
    
  vector<shared_ptr<const PeakDef>> user_peaks = get_user_peaks();
  
  const bool isHPGe = PeakFitUtils::is_likely_high_res( m_interspec );
  
  //The results of the peak search will be placed into the vector pointed to by searchresults
  auto searchresults = std::make_shared< vector<std::shared_ptr<const PeakDef> > >();
    
  //I think WApplication::bind should protect against calling back after this widget is deleted
  boost::function<void(void)> callback
                            = wApp->bind( boost::bind( &MakeFwhmForDrf::setPeaksFromAutoSearch,
                                          this, user_peaks, searchresults ) );
    
  
  Wt::WServer *server = Wt::WServer::instance();
  assert( server );
  if( !server )
    return;
  
  const string seshid = wApp->sessionId();
  const shared_ptr<const DetectorPeakResponse> drf = m_orig_drf;
    
  server->ioService().boost::asio::io_service::post( std::bind( [=](){
    const bool singleThread = false;
    auto existingPeaks = make_shared<deque<shared_ptr<const PeakDef>>>();
    existingPeaks->insert( end(*existingPeaks), begin(user_peaks), end(user_peaks) );
    
    *searchresults = ExperimentalAutomatedPeakSearch::search_for_peaks( dataPtr, drf, existingPeaks, singleThread, isHPGe );
    
    Wt::WServer *server = Wt::WServer::instance();
    if( server )
      server->post( seshid, callback );
  } ) );
}//void startAutomatedPeakSearch();


void MakeFwhmForDrf::scheduleRefit()
{
  m_refit_scheduled = true;
  scheduleRender();
}//void scheduleRefit()


void MakeFwhmForDrf::doRefitWork()
{
  m_parameters.clear();
  m_uncertainties.clear();
  
  try
  {
    vector<shared_ptr<const PeakDef>> peaks = m_model->peaks_to_use();
    auto peaks_deque = make_shared<deque<shared_ptr<const PeakDef>>>();
    peaks_deque->insert( end(*peaks_deque), begin(peaks), end(peaks) );
    
    int sqrtEqnOrder = -1;
    
    const int fwhm_type_int = m_fwhmEqnType->currentIndex();
    
    if( (fwhm_type_int < 0)
       || (fwhm_type_int >= DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm) )
      throw runtime_error( "Invalid function type" );
    
    const auto fwhm_type = DetectorPeakResponse::ResolutionFnctForm( fwhm_type_int );
    
    string eqn;
    switch( fwhm_type )
    {
      case DetectorPeakResponse::kGadrasResolutionFcn:
      case DetectorPeakResponse::kConstantPlusSqrtEnergy:
        eqn = "FWHM(x): ";
        break;
        
      case DetectorPeakResponse::kSqrtEnergyPlusInverse:
        //sqrt(pars[0] + pars[1]*energy + pars[2]/energy);
        eqn = "FWHM(x) = sqrt( ";
        break;
        
      case DetectorPeakResponse::kSqrtPolynomial:
        //FWHM = sqrt( Sum_i{A_i*pow(x/1000,i)} );
        eqn = "FWHM(x) = sqrt( ";
        sqrtEqnOrder = m_sqrtEqnOrder->currentIndex() + 1;
        break;
        
      case DetectorPeakResponse::kNumResolutionFnctForm:
        assert( 0 );
        break;
    }//switch( fwhm_type )
    
    auto meas = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    
    vector<float> result, uncerts;
    const double chi2 = MakeDrfFit::performResolutionFit( peaks_deque, fwhm_type,
                                                         sqrtEqnOrder, result, uncerts );
    
    m_parameters = result;
    m_uncertainties = uncerts;
    
    for( size_t i = 0; (i < m_parEdits.size()) && (i < result.size()); ++i )
    {
      m_parEdits[i]->setValue( result[i] );
      m_parEdits[i]->enable();
      assert( !m_parEdits[i]->isHidden() );
      m_parEdits[i]->setHidden(false);
      if( i < uncerts.size() )
        m_parEdits[i]->setToolTip( WString::tr("mffd-tt-par").arg(SpecUtils::printCompact(uncerts[i], 5)) );
      else
        m_parEdits[i]->setToolTip( "" );
      
      const string val_str = SpecUtils::printCompact(result[i], 6);
      
      switch( fwhm_type )
      {
        case DetectorPeakResponse::kGadrasResolutionFcn:
          eqn += (i ? ", " : "") + string("A") + std::to_string(i) + "=" + val_str;
          break;
          
        case DetectorPeakResponse::kSqrtEnergyPlusInverse:
          //sqrt(pars[0] + pars[1]*energy + pars[2]/energy);
          if( i == 0 )
            eqn += val_str;
          else if( i == 1 )
            eqn += " + " + val_str + "*x";
          else if( i == 2 )
            eqn += " + " + val_str + "/x";
          break;
          
        case DetectorPeakResponse::kConstantPlusSqrtEnergy:
          //pars[0] + pars[1]*sqrt(energy);
          if( i == 0 )
            eqn += val_str;
          else if( i == 1 )
            eqn += " + " + val_str + "*sqrt(x)";
          break;
          
        case DetectorPeakResponse::kSqrtPolynomial:
          //FWHM = sqrt( Sum_i{A_i*pow(x/1000,i)} );
          if( i == 0 )
            eqn += val_str;
          else if( i == 1 )
            eqn += " + " + val_str + "*x/1000";
          else
            eqn += " + " + val_str + "*(x/1000)^" + std::to_string(i);
          break;
          
        case DetectorPeakResponse::kNumResolutionFnctForm:
          assert( 0 );
          break;
      }//switch( fwhm_type )
    }//for( size_t i = 0; (i < m_parEdits.size()) && (i < result.size()); ++i )
    
    if( fwhm_type != DetectorPeakResponse::kGadrasResolutionFcn )
      eqn += " )";
    
    m_error->hide();
    m_error->setText( "" );
    m_equation->show();
    m_equation->setText( eqn );
    
    m_validationChanged.emit(true);
  }catch( std::exception &e )
  {
    m_error->show();
    m_equation->hide();
    m_equation->setText( "" );
    
    for( NativeFloatSpinBox *sb : m_parEdits )
    {
      sb->disable();
      sb->setText( "" );
    }
    
    m_error->setText( WString::tr("mffd-err-fail-fit").arg(e.what()) );
    m_validationChanged.emit(false);
  }//try / catch
  
  setEquationToChart();
}//void doRefitWork()


void MakeFwhmForDrf::setEquationToChart()
{
  auto meas = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  vector<shared_ptr<const PeakDef>> peaks = m_model->peaks_to_use();
  
  const int fwhm_type_int = m_fwhmEqnType->currentIndex();
  
  if( (fwhm_type_int < 0)
     || (fwhm_type_int >= DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm) )
    throw runtime_error( "Invalid function type" );
  
  const auto fwhm_type = DetectorPeakResponse::ResolutionFnctForm( fwhm_type_int );
  
  MakeDrfChart::FwhmCoefType chart_fwhm_type = MakeDrfChart::FwhmCoefType::Gadras;
  switch( fwhm_type )
  {
    case DetectorPeakResponse::kGadrasResolutionFcn:
      chart_fwhm_type = MakeDrfChart::FwhmCoefType::Gadras;
      break;
    case DetectorPeakResponse::kSqrtPolynomial:
      chart_fwhm_type = MakeDrfChart::FwhmCoefType::SqrtEqn;
      break;
    case DetectorPeakResponse::kSqrtEnergyPlusInverse:
      chart_fwhm_type = MakeDrfChart::FwhmCoefType::SqrtEnergyPlusInverse;
      break;
    case DetectorPeakResponse::kConstantPlusSqrtEnergy:
      chart_fwhm_type = MakeDrfChart::FwhmCoefType::ConstantPlusSqrtEnergy;
      break;
      
    case DetectorPeakResponse::kNumResolutionFnctForm:
      assert( 0 );
      break;
  }//switch( fwhm_type )
  
  float min_energy = -999, max_energy = -999;
  vector<MakeDrfChart::DataPoint> datapoints;
  for( const auto &peak : peaks )
  {
    MakeDrfChart::DataPoint point;
    point.energy = static_cast<float>( peak->mean() );
    point.livetime = meas ? meas->live_time() : 1.0f;
    point.peak_area = peak->peakArea();
    point.peak_area_uncertainty = 2.35482f*peak->peakAreaUncert();
    point.peak_fwhm = peak->fwhm();
    point.peak_fwhm_uncertainty = 2.35482f*peak->sigmaUncert();
    point.source_count_rate = 0;
    point.source_count_rate_uncertainty = 0;
    point.distance = 1.0*PhysicalUnits::meter;
    point.source_information = "";
    point.peak_color = peak->lineColor();
    point.background_peak_area = 0;
    point.background_peak_live_time = 0;
    
    min_energy = (min_energy < 0.0) ? point.energy : std::min(point.energy, min_energy);
    max_energy = (max_energy < 0.0) ? point.energy : std::max(point.energy, max_energy);
    
    datapoints.push_back( std::move(point) );
  }//for( const auto &peak : peaks )
  
  const float det_diameter = m_orig_drf ? m_orig_drf->detectorDiameter()
                                        : static_cast<float>(2.54f*PhysicalUnits::cm);
  
  m_chart->setDataPoints( datapoints, det_diameter, min_energy, max_energy );
  m_chart->setFwhmCoefficients( m_parameters, m_uncertainties, chart_fwhm_type, MakeDrfChart::EqnEnergyUnits::keV );
}//void setEquationToChart()


void MakeFwhmForDrf::setPeaksFromAutoSearch( vector<shared_ptr<const PeakDef>> user_peaks,
                             shared_ptr<vector<shared_ptr<const PeakDef>>> auto_search_peaks )
{
  assert( auto_search_peaks );
  
  // `auto_search_peaks` will contain both the original users peaks, as well as the auto-fit
  //  peaks, but we want to keep them separate to indicate to the user, so we'll just remove
  //  the peaks in `auto_search_peaks` that overlap with the user peaks.  Not perfect, but
  //  good enough probably.
  vector<shared_ptr<const PeakDef>> filtered_auto_search;
  if( auto_search_peaks )
  {
    for( const auto &p : *auto_search_peaks )
    {
      bool overlaps_user = false;
      const double mean = p->mean();
      for( size_t i = 0; !overlaps_user && (i < user_peaks.size()); ++i )
        overlaps_user = (fabs(user_peaks[i]->mean() - mean) < 0.5*user_peaks[i]->fwhm());
      if( !overlaps_user )
        filtered_auto_search.push_back( p );
    }//for( const auto &p : *auto_search_peaks )
  }//if( auto_search_peaks )
  
  
  m_model->set_peaks( user_peaks, filtered_auto_search ); //should trigger a re-fit
  
  
  // If more than a handful of peaks, fit to a function, then go through and uncheck outliers
  vector<shared_ptr<const PeakDef>> initial_fit_peaks = m_model->peaks_to_use();
  if( initial_fit_peaks.size() > 9 )  //9 is arbitrary
  {
    try
    {
      auto peaks_deque = make_shared<deque<shared_ptr<const PeakDef>>>();
      peaks_deque->insert( end(*peaks_deque), begin(initial_fit_peaks), end(initial_fit_peaks) );
      
      int sqrt_eqn_order = 3;
      const auto fwhm_type = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
      
      auto meas = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
      
      vector<float> result, uncerts;
      const double chi2 = MakeDrfFit::performResolutionFit( peaks_deque, fwhm_type,
                                                            sqrt_eqn_order, result, uncerts );
      
      vector<pair<double,shared_ptr<const PeakDef>>> distances;
      for( const auto &p : initial_fit_peaks )
      {
        double pred_fwhm = 0.0;
        const double mean = 0.001 * p->mean();
        for( size_t i = 0; i < result.size(); ++i )
          pred_fwhm += result[i] * std::pow( mean, 1.0*i );
        pred_fwhm = sqrt( pred_fwhm );
        const double frac_diff = fabs( p->fwhm() - pred_fwhm ) / p->fwhm();
        if( !IsNan(frac_diff) && !IsInf(frac_diff) )
          distances.emplace_back( frac_diff, p );
      }//for( const auto &p : initial_fit_peaks )
      
      std::sort( begin(distances), end(distances),
        []( const pair<double,shared_ptr<const PeakDef>> &lhs,
            const pair<double,shared_ptr<const PeakDef>> &rhs) -> bool {
          return lhs.first > rhs.first;
        } );

      // Limit to un-selecting max of 20% of peaks (arbitrarily chosen), if the deviate
      // more than 17.5% from the fit (again, arbitrarily chosen).
      const size_t max_remove = static_cast<size_t>( std::ceil( 0.2*distances.size() ) );
      for( size_t index = 0; (index < max_remove) && (index < distances.size()); ++index )
      {
        if( distances[index].first < 0.175 ) //0.175 chosen arbitrarily
          break;
        m_model->set_use( distances[index].second, false );
      }
    }catch( std::exception & )
    {
      // Oh, well, wasnt meant to be
    }//try / catch
  }//if( a good number of peaks )
  
  
  m_current_state = currentState();
  wApp->triggerUpdate();
}//setPeaksFromAutoSearch(...)


void MakeFwhmForDrf::handleTableDataChange()
{
  scheduleUndoRedoStep();
  scheduleRefit();
}//void handleTableDataChange();


void MakeFwhmForDrf::handleFwhmEqnTypeChange()
{
  size_t num_pars = 0;
  const int fwhm_index = std::max(0,m_fwhmEqnType->currentIndex());
  const auto fwhm_type = DetectorPeakResponse::ResolutionFnctForm( fwhm_index );
  bool showSqrtEqn = false;
  switch( fwhm_type )
  {
    case DetectorPeakResponse::kGadrasResolutionFcn:
      num_pars = 3;
      break;
      
    case DetectorPeakResponse::kSqrtPolynomial:
      showSqrtEqn = true;
      num_pars = static_cast<int>( std::max(0,m_sqrtEqnOrder->currentIndex()) + 1);
      break;
      
    case DetectorPeakResponse::kSqrtEnergyPlusInverse:
      num_pars = 3;
      break;
      
    case DetectorPeakResponse::kConstantPlusSqrtEnergy:
      num_pars = 2;
      break;
      
    case DetectorPeakResponse::kNumResolutionFnctForm:
      assert( 0 );
      return;
      break;
  }//switch( fwhm_type )
  
  if( m_sqrtEqnOrder->parent() )
    m_sqrtEqnOrder->parent()->setHidden( !showSqrtEqn );
  
  for( size_t i = 0; i < num_pars; ++i )
  {
    if( m_parEdits[i]->parent() )
      m_parEdits[i]->parent()->show();
  }
  
  for( size_t i = num_pars; i < m_parEdits.size(); ++i )
  {
    m_parEdits[i]->setValue( 0.0f );
    
    if( m_parEdits[i]->parent() )
      m_parEdits[i]->parent()->hide();
  }
  
  scheduleUndoRedoStep();
  scheduleRefit();
}//void handleFwhmEqnTypeChange();


void MakeFwhmForDrf::handleSqrtEqnOrderChange()
{
  handleFwhmEqnTypeChange();
  scheduleUndoRedoStep();
}//void handleSqrtEqnOrderChange()


void MakeFwhmForDrf::coefficientManuallyChanged( const int coef_num )
{
  assert( (coef_num >= 0) && (coef_num < static_cast<int>(m_parEdits.size())) );
  
  if( (coef_num < 0) || (coef_num >= static_cast<int>(m_parEdits.size())) )
  {
    return;
  }
  
  const float value = m_parEdits[coef_num]->value();
  assert( coef_num < m_parameters.size() );
  if( coef_num >= m_parameters.size() )
    m_parameters.resize(coef_num + 1, 0.0f );
  m_parameters[coef_num] = value;
  m_uncertainties.clear();
  
  setEquationToChart();
  scheduleUndoRedoStep();
}//void coefficientManuallyChanged( const int coef_num )


Wt::Signal<bool> &MakeFwhmForDrf::validationChanged()
{
  return m_validationChanged;
}


bool MakeFwhmForDrf::isValidFwhm() const
{
  return !m_parameters.empty();
}//bool isValidFwhm() const


void MakeFwhmForDrf::setToDrf()
{
  assert( !m_parameters.empty() );
  if( m_parameters.empty() )
    return;
  
  // Set current FWHM to DRF and broadcast out to the rest of the app
  shared_ptr<SpecMeas> foreground = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  shared_ptr<DetectorPeakResponse> prev_det = foreground ? foreground->detector() : nullptr;
  
  shared_ptr<DetectorPeakResponse> new_det;
  if( prev_det )
  {
    new_det = make_shared<DetectorPeakResponse>( *prev_det );
    new_det->setParentHashValue( prev_det->hashValue() );
  }else
  {
    new_det = make_shared<DetectorPeakResponse>( "Flat Response", "FWHM info only" );
    new_det->setIntrinsicEfficiencyFormula( "1.0", 2.54*PhysicalUnits::cm, PhysicalUnits::keV,
                                            0.0f, 0.0f, DetectorPeakResponse::EffGeometryType::FarField );
  }
  
  const auto fwhm_type = DetectorPeakResponse::ResolutionFnctForm( std::max(0,m_fwhmEqnType->currentIndex()) );
  new_det->setFwhmCoefficients( m_parameters, fwhm_type );
  m_interspec->detectorChanged().emit(new_det);
  
  UndoRedoManager *undoManager = m_interspec->undoRedoManager();
  if( undoManager && undoManager->canAddUndoRedoNow() )
  {
    auto undo = [prev_det](){
      InterSpec::instance()->detectorChanged().emit(prev_det);
    };
    
    auto redo = [new_det](){
      InterSpec::instance()->detectorChanged().emit(new_det);
    };
    
    undoManager->addUndoRedoStep( undo, redo, "Change FWHM for DRF" );
  }//if( m_parent && undoManager )
  
  m_updatedDrf.emit(new_det);
}//void setToDrf()


Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &MakeFwhmForDrf::updatedDrf()
{
  return m_updatedDrf;
}//Wt::Signal<> &updatedDrf()


std::shared_ptr<MakeFwhmForDrf::ToolState> MakeFwhmForDrf::currentState() const
{
  auto answer = make_shared<MakeFwhmForDrf::ToolState>();
  
  answer->m_fwhm_index = m_fwhmEqnType->currentIndex();
  answer->m_sqrt_eqn_index = m_sqrtEqnOrder->currentIndex();
  answer->m_rows = m_model->rowData();
  answer->m_parameters = m_parameters;
  answer->m_uncertainties = m_uncertainties;
  answer->m_orig_drf = m_orig_drf;
  
  return answer;
}//std::shared_ptr<ToolState> currentState() const


void MakeFwhmForDrf::setState( shared_ptr<const MakeFwhmForDrf::ToolState> state )
{
  assert( state );
  if( !state )
    return;
  
  m_current_state = state;
  
  assert( (state->m_fwhm_index >= 0) && (state->m_fwhm_index < m_fwhmEqnType->count()) );
  if( (state->m_fwhm_index >= 0) && (state->m_fwhm_index < m_fwhmEqnType->count()) )
    m_fwhmEqnType->setCurrentIndex( state->m_fwhm_index );
  
  assert( (state->m_sqrt_eqn_index >= 0) && (state->m_sqrt_eqn_index < m_sqrtEqnOrder->count()) );
  if( (state->m_sqrt_eqn_index >= 0) && (state->m_sqrt_eqn_index < m_sqrtEqnOrder->count()) )
    m_sqrtEqnOrder->setCurrentIndex( state->m_sqrt_eqn_index );
  
  handleFwhmEqnTypeChange();
  
  m_model->setRowData( state->m_rows );
  
  if( state->m_uncertainties.empty() && !state->m_parameters.empty() )
  {
    m_parameters = state->m_parameters;
    m_uncertainties = state->m_uncertainties;
    for( size_t i = 0; (i < m_parEdits.size()) && (i < m_parameters.size()); ++i )
    {
      m_parEdits[i]->setValue( m_parameters[i] );
      m_parEdits[i]->setToolTip( "" );
    }
  }//if( the user manually set parameter values )
  
  m_orig_drf = state->m_orig_drf;
  
  scheduleRefit(); //It should already be scheduled, but just to be explicit
}//void setState( std::shared_ptr<const ToolState> state )


void MakeFwhmForDrf::doAddUndoRedoStep()
{
  UndoRedoManager *undoManager = m_interspec->undoRedoManager();
  if( !undoManager || !undoManager->canAddUndoRedoNow() )
    return;
  
  shared_ptr<MakeFwhmForDrf::ToolState> state = currentState();
  assert( state );
  if( !state )
    return;
  
  if( !m_current_state )
  {
    m_current_state = state;
    return;
  }
  
  if( (*state) == (*m_current_state) )
    return;
  
  shared_ptr<const MakeFwhmForDrf::ToolState> prev_state = m_current_state;
  m_current_state = state;
  
  auto undo = [prev_state](){
    InterSpec *viewer = InterSpec::instance();
    MakeFwhmForDrfWindow *window = viewer ? viewer->fwhmFromForegroundWindow(false) : nullptr;
    if( window )
      window->tool()->setState( prev_state );
  };
    
  auto redo = [this,state](){
    InterSpec *viewer = InterSpec::instance();
    MakeFwhmForDrfWindow *window = viewer ? viewer->fwhmFromForegroundWindow(false) : nullptr;
    if( window )
      window->tool()->setState( state );
  };
    
  undoManager->addUndoRedoStep( undo, redo, "FWHM from spectrum tool change" );
}//void doAddUndoRedoStep()


void MakeFwhmForDrf::scheduleUndoRedoStep()
{
  UndoRedoManager *undoManager = m_interspec->undoRedoManager();
  if( undoManager && undoManager->canAddUndoRedoNow() )
  {
    m_undo_redo_scheduled = true;
    scheduleRender();
  }//if( undoManager && undoManager->canAddUndoRedoNow() )
}//void scheduleUndoRedoStep()


bool MakeFwhmForDrf::ToolState::operator==( const MakeFwhmForDrf::ToolState &rhs ) const
{
  if( (m_fwhm_index != rhs.m_fwhm_index)
     || (m_sqrt_eqn_index != rhs.m_sqrt_eqn_index)
     || (m_parameters != rhs.m_parameters)
     || (m_uncertainties != rhs.m_uncertainties)
     || (m_rows.size() != rhs.m_rows.size()) )
  {
    return false;
  }
     
  vector<MakeFwhmForDrf::TableRow> lhs_row = m_rows;
  vector<MakeFwhmForDrf::TableRow> rhs_row = rhs.m_rows;
  FwhmPeaksModel::sortImp( FwhmPeaksModel::Column::Energy, Wt::SortOrder::AscendingOrder, lhs_row );
  FwhmPeaksModel::sortImp( FwhmPeaksModel::Column::Energy, Wt::SortOrder::AscendingOrder, rhs_row );
  
  assert( lhs_row.size() == rhs_row.size() );
  for( size_t i = 0; (i < lhs_row.size()) && (i < rhs_row.size()); ++i )
  {
    if( (lhs_row[i].m_is_user_peak != rhs_row[i].m_is_user_peak)
       || (lhs_row[i].m_use_for_fit != rhs_row[i].m_use_for_fit)
      || (lhs_row[i].m_peak != rhs_row[i].m_peak) )
    {
      return false;
    }
  }
  
  return true;
}//ToolState::operator==
