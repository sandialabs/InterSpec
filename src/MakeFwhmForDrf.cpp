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
#include <Wt/WIOService>
#include <Wt/WTableView>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WAbstractItemModel>

#include "SpecUtils/SpecFile.h"

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


using namespace std;
using namespace Wt;


class FwhmPeaksModel : public  Wt::WAbstractItemModel
{
public:
  enum class Column : int { Energy, Fwhm, FWhmUncert, UserOrAuto, UseForFit, NumColumn };
  
  struct Row
  {
    bool m_is_user_peak;
    bool m_use_for_fit;
    std::shared_ptr<const PeakDef> m_peak;
  };//struct Row
  
protected:
  Column m_sort_col;
  SortOrder m_sort_order = AscendingOrder;
  std::vector<Row> m_rows;
  
  
public:
  FwhmPeaksModel( Wt::WObject *parent = nullptr )
  : Wt::WAbstractItemModel( parent ),
  m_sort_col( Column::Energy )
  {
  }
  
  virtual ~FwhmPeaksModel() override {}
  
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
      case Column::Energy:     return WString("Energy (keV)");
      case Column::Fwhm:       return WString("FWHM");
      case Column::FWhmUncert: return WString("FWHM Uncert");
      case Column::UserOrAuto: return WString("Peak Source");
      case Column::UseForFit:  return WString("Use");
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
    
    const Row &row = m_rows[index.row()];
    
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
        snprintf( buffer, sizeof(buffer), "%s", (row.m_is_user_peak ? "User" : "Auto-fit") );
        break;
        
      case Column::UseForFit:
        return boost::any(row.m_use_for_fit);
        
      case Column::NumColumn:
        assert( 0 );
        break;
    }//switch( Column(index.column()) )
    
    return boost::any( WString(buffer) );
  }//data(...)
  
  
  virtual bool setData( const WModelIndex &index, const boost::any &value, int role = EditRole ) override
  {
    if( index.parent().isValid() || (index.row() < 0)
       || (index.row() >= static_cast<int>(m_rows.size()))
       || (index.column() != static_cast<int>(Column::UseForFit))
       || value.empty()
       || (value.type() != typeid(bool)) )
    {
      assert( 0 );
      return;
    }
    
    m_rows[index.row()].m_use_for_fit = boost::any_cast<bool>(value);
    dataChanged().emit( index, index );
  }//bool setData( const WModelIndex &index, const boost::any &value, int role = EditRole )
  
  
  void sortImp( Column column, SortOrder order )
  {
    const bool littleToBig = (order == AscendingOrder);
    
    std::sort( begin(m_rows), end(m_rows), [column, littleToBig]( const Row &lhs, const Row &rhs ) -> bool {
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
    
    sortImp( m_sort_col, m_sort_order );
    
    layoutChanged().emit();
  }//sort(...)
  
  
  void set_peaks( const vector<shared_ptr<const PeakDef>> &user_peaks,
                 const vector<shared_ptr<const PeakDef>> &auto_fit_peaks )
  {
    if( m_rows.empty() )
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
      Row r;
      r.m_is_user_peak = true;
      r.m_use_for_fit = p->useForDrfFwhmFit();
      r.m_peak = p;
      m_rows.push_back( std::move(r) );
    }//for( const auto &p : user_peaks )
    
    for( const auto &p : auto_fit_peaks )
    {
      Row r;
      r.m_is_user_peak = false;
      r.m_use_for_fit = p->useForDrfFwhmFit();
      r.m_peak = p;
      m_rows.push_back( std::move(r) );
    }//for( const auto &p : user_peaks )
    
    sortImp( m_sort_col, m_sort_order );
    
    endInsertRows();
  }//void setPeaks(...)
  
  vector<shared_ptr<const PeakDef>> peaks_to_use() const
  {
    vector<shared_ptr<const PeakDef>> answer;
    for( const auto &row : m_rows )
    {
      if( row.m_use_for_fit )
        answer.push_back( row.m_peak );
    }
    
    return answer;
  }//vector<shared_ptr<const PeakDef>> peaks_to_use() const
};//class FwhmPeaksModel
 


pair<AuxWindow *,MakeFwhmForDrf *> MakeFwhmForDrf::makeAddFwhmToDrfWindow()
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  if( !viewer )
    return pair<AuxWindow *,MakeFwhmForDrf *>{nullptr, nullptr};
  
  AuxWindow *window = new AuxWindow( "Add FWHM to Detector Response Function",
                                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                                     | AuxWindowProperties::SetCloseable
                                     | AuxWindowProperties::DisableCollapse
                                     | AuxWindowProperties::EnableResize
                                     | AuxWindowProperties::IsModal) );
  
  const int ww = viewer->renderedWidth();
  const int wh = viewer->renderedHeight();
  if( ww > 100 && wh > 100 )
  {
    const int width = std::min( 3*ww/4, 900 );
    const int height = ((wh < 420) ? wh : (19*wh)/20 );
    
    window->resizeWindow( width, height );
    window->setMinimumSize( std::min(width,640), std::min(height,480) );
  }//if( ww > 100 && wh > 100 )
  
  shared_ptr<const SpecMeas> foreground = viewer->measurment(SpecUtils::SpectrumType::Foreground );
  shared_ptr<const DetectorPeakResponse> drf = foreground ? foreground->detector() : nullptr;
  
  MakeFwhmForDrf *makeFwhmWidget = new MakeFwhmForDrf( viewer, drf );
    
  window->stretcher()->addWidget( makeFwhmWidget, 0, 0 );
  window->stretcher()->setContentsMargins( 0, 0, 0, 0 );
  
  AuxWindow::addHelpInFooter( window->footer(), "add-fwhm-to-drf" );
    
  WPushButton *closeButton = window->addCloseButtonToFooter( "Close" );
  closeButton->clicked().connect( window, &AuxWindow::hide );
    
  WPushButton *saveAs = new WPushButton( "Use", window->footer() );
  saveAs->clicked().connect( makeFwhmWidget, &MakeFwhmForDrf::setToDrf );
  makeFwhmWidget->validationChanged().connect( boost::bind( &WPushButton::setEnabled, saveAs,
                                                                   boost::placeholders::_1 ) );
  // Maybe 
  //makeFwhmWidget->validationChanged().connect( std::bind([makeFwhmWidget,saveAs](){
  //  saveAs->setEnabled( makeFwhmWidget->isValidFwhm() );
  //}) );
  
  saveAs->disable();
    
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    
  window->show();
    
  window->resizeToFitOnScreen();
  window->centerWindow();
  window->rejectWhenEscapePressed( false );
  
  //PeakModel *peakModel = viewer->peakModel();
  //shared_ptr<const deque< PeakModel::PeakShrdPtr > > peakModel->peaks();
  
  return pair<AuxWindow *,MakeFwhmForDrf *>{ window, makeFwhmWidget };
}//AuxWindow *makeDrfWindow(...)


MakeFwhmForDrf::MakeFwhmForDrf( InterSpec *viewer,
               std::shared_ptr<const DetectorPeakResponse> drf,
               Wt::WContainerWidget *parent )
 : WContainerWidget( parent ),
  m_interspec( viewer ),
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
  
  // The chart requires using a Wt layout, so we'll go with it
  // I think this should be made four rows.
  //  First is the chart
  //  Second is equation txt
  //  Third is error txt
  //  The Fourth row has all the options and values on left, and then the table on the right
    
  WGridLayout *layout = new WGridLayout( this );
  //layout->setVerticalSpacing( 0 );
  //layout->setHorizontalSpacing( 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
    
  m_chart = new MakeDrfChart();
  m_chart->showEfficiencyPoints( false );
  DrfChartHolder *chartholder = new DrfChartHolder( m_chart, nullptr );
  layout->addWidget( chartholder, layout->rowCount(), 0 );
  chartholder->setHeight( 250 );
  //layout->setRowResizable( 0, true, WLength(250, WLength::Pixel) );
  //layout->setRowStretch( 0, 2 );
    
    
  WContainerWidget *optionsDiv = new WContainerWidget();
  layout->addWidget( optionsDiv, layout->rowCount(), 0 );
  optionsDiv->addStyleClass( "Options" );
    
  WContainerWidget *optDiv = new WContainerWidget( optionsDiv );
  optDiv->addStyleClass( "OptDiv" );
  WLabel *label = new WLabel( "FWHM Equation Type:", optDiv );
  m_fwhmEqnType = new WComboBox( optDiv );
  
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
          
      case DetectorPeakResponse::kNumResolutionFnctForm:
        break;
    }//switch( fcnfrm )
  }//for( loop over functional forms )
    
  m_fwhmEqnType->setCurrentIndex( DetectorPeakResponse::kSqrtPolynomial );
  m_fwhmEqnType->activated().connect( this, &MakeFwhmForDrf::handleFwhmEqnTypeChange );
  
  optDiv = new WContainerWidget( optionsDiv );
  optDiv->addStyleClass( "OptDiv" );
  label = new WLabel( "Num Terms:", optDiv );
  m_sqrtEqnOrder = new WComboBox( optDiv );
  m_sqrtEqnOrder->addItem( "1" );
  m_sqrtEqnOrder->addItem( "2" );
  m_sqrtEqnOrder->addItem( "3" );
  m_sqrtEqnOrder->addItem( "4" );
  m_sqrtEqnOrder->addItem( "5" );
  m_sqrtEqnOrder->setCurrentIndex( 1 );
  m_sqrtEqnOrder->activated().connect( this, &MakeFwhmForDrf::handleSqrtEqnOrderChange );
  
  WContainerWidget *parametersDiv = new WContainerWidget();
  layout->addWidget( parametersDiv, layout->rowCount(), 0 );
  parametersDiv->addStyleClass( "Parameters" );
    
  for( int i = 0; i < m_sqrtEqnOrder->count(); ++i )
  {
    WContainerWidget *parDiv = new WContainerWidget( parametersDiv );
    parDiv->addStyleClass( "ParDiv" );
    label = new WLabel( "A" + std::to_string(i), parDiv );
    NativeFloatSpinBox *sb = new NativeFloatSpinBox( parDiv );
    sb->setText( "" );
    sb->disable();
    sb->valueChanged().connect( boost::bind(&MakeFwhmForDrf::coefficientManuallyChanged, this, i) );
    if( i >= m_sqrtEqnOrder->currentIndex() )
      parDiv->hide();
    m_parEdits.push_back( sb );
  }//for( int i = 0; i < m_sqrtEqnOrder->count(); ++i )
    
  m_error = new WText( "", Wt::XHTMLText, this );
  layout->addWidget( m_error, layout->rowCount(), 0 );
  m_error->setInline( false );
  m_error->addStyleClass( "ErrTxt" );
  m_error->hide();
  m_equation = new WText( "", Wt::XHTMLText, this );
  layout->addWidget( m_equation, layout->rowCount(), 0 );
  m_equation->setInline( false );
  m_equation->addStyleClass( "EqnTxt" );
  m_equation->hide();
    
  m_model = new FwhmPeaksModel( this );
  m_table = new WTableView();
  layout->addWidget( m_table, layout->rowCount(), 0 );
  //layout->setRowStretch( 4, 3 );
  m_table->setModel( m_model );
  
  m_model->dataChanged().connect( this, &MakeFwhmForDrf::refit );
  m_model->rowsInserted().connect( this, &MakeFwhmForDrf::refit );
  m_model->rowsRemoved().connect( this, &MakeFwhmForDrf::refit );
  m_model->layoutChanged().connect( this, &MakeFwhmForDrf::refit );
   
  startAutomatedPeakSearch();
    
  refit();
}//MakeFwhmForDrf( constructor )


MakeFwhmForDrf::~MakeFwhmForDrf()
{
  
}//~MakeFwhmForDrf()


void MakeFwhmForDrf::startAutomatedPeakSearch()
{
  PeakModel *peakModel = m_interspec->peakModel();
  assert( peakModel );
  if( !peakModel )
    return; //shouldnt happen
  
  auto dataPtr = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  assert( dataPtr );
  if( !dataPtr )
  {
    passMessage( "There is no data to search for peaks on.", 1 );
    return;
  }//if( !dataPtr )
    
  const auto originalPeaks = peakModel->peakVec();
  std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > startingPeaks = peakModel->peaks();
  vector<shared_ptr<const PeakDef>> user_peaks;
  if( startingPeaks && !startingPeaks->empty() )
    user_peaks.insert( end(user_peaks), begin(*startingPeaks), end(*startingPeaks) );
    
  
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
    
    *searchresults = ExperimentalAutomatedPeakSearch::search_for_peaks( dataPtr, drf, existingPeaks, singleThread );
    
    Wt::WServer *server = Wt::WServer::instance();
    if( server )
      server->post( seshid, callback );
  } ) );
}//void startAutomatedPeakSearch();


void MakeFwhmForDrf::refit()
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
    const bool highres = PeakFitUtils::is_high_res(meas);
    
    vector<float> result, uncerts;
    const double chi2 = MakeDrfFit::performResolutionFit( peaks_deque, fwhm_type, highres,
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
        m_parEdits[i]->setToolTip( "With uncertainty: " + PhysicalUnits::printCompact(uncerts[i], 5) );
      else
        m_parEdits[i]->setToolTip( "" );
      
      const string val_str = PhysicalUnits::printCompact(result[i], 6);
      
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
    
    string msg = "Failed to fit FWHM equation: " + string( e.what() );
    m_error->setText( msg );
    m_validationChanged.emit(false);
  }//try / catch
  
  setEquationToChart();
}//void refit()


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


void MakeFwhmForDrf::setPeaksFromAutoSearch( std::vector<std::shared_ptr<const PeakDef>> user_peaks,
                             std::shared_ptr<std::vector<std::shared_ptr<const PeakDef>>> auto_search_peaks )
{
  assert( auto_search_peaks );
  
  m_model->set_peaks( user_peaks, *auto_search_peaks ); //should trigger a re-fit
  
  wApp->triggerUpdate();
}//setPeaksFromAutoSearch(...)


void MakeFwhmForDrf::handleFwhmEqnTypeChange()
{
  size_t num_pars = 0;
  const auto fwhm_type = DetectorPeakResponse::ResolutionFnctForm( std::max(0,m_fwhmEqnType->currentIndex()) );
  switch( fwhm_type )
  {
    case DetectorPeakResponse::kGadrasResolutionFcn:
      num_pars = 3;
      break;
      
    case DetectorPeakResponse::kSqrtPolynomial:
      if( m_sqrtEqnOrder->parent() )
        m_sqrtEqnOrder->parent()->setHidden( false );
      num_pars = static_cast<int>( std::max(0,m_sqrtEqnOrder->currentIndex()) + 1);
      break;
      
    case DetectorPeakResponse::kSqrtEnergyPlusInverse:
      num_pars = 3;
      break;
      
    case DetectorPeakResponse::kNumResolutionFnctForm:
      assert( 0 );
      return;
      break;
  }//switch( fwhm_type )
  
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
  
  refit();
}//void handleFwhmEqnTypeChange();


void MakeFwhmForDrf::handleSqrtEqnOrderChange()
{
  handleFwhmEqnTypeChange();
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
    new_det = make_shared<DetectorPeakResponse>( "Flat Detector", "FWHM info only" );
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
