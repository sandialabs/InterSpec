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

#include <map>
#include <string>
#include <vector>
#include <memory>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <exception>
#include <sys/stat.h>

#include <boost/any.hpp>

#include <Wt/WPen>
#include <Wt/WText>
#include <Wt/WDate>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WRectF>
#include <Wt/WLength>
#include <Wt/WSlider>
#include <Wt/WAnchor>
#include <Wt/WDialog>
#include <Wt/WServer>
#include <Wt/WPainter>
#include <Wt/WSpinBox>
#include <Wt/WDateEdit>
#include <Wt/WCheckBox>
#include <Wt/WGroupBox>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WResource>
#include <Wt/WTabWidget>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WPaintDevice>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WStandardItemModel>
#include <Wt/Chart/WCartesianChart>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayChainChart.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DecayActivityDiv.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecaySelectNuclideDiv.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

using namespace Wt;
using namespace std;
using namespace PhysicalUnits;


#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif


namespace
{
  string remove_space_copy( string input )
  {
    input.erase( std::remove_if(input.begin(), input.end(), [](char a){return isspace(a);}), input.end() );
    return input;
  }


}//namespace

class DecayActivityModel : public Wt::WStandardItemModel
{
#if( WT_VERSION > 0x3040000 )
  // HACK:
  // Wt 3.7.1 has an issue for really, really small values, it will actually draw the line at
  //  the top of the chart, when it should be at the bottom.
  //  So we will do this hack that sets the value to zero, for anything effectively at zero
  //  compared to what we are showing.
  double m_max_activity;
#endif
  
public:
  DecayActivityModel( Wt::WObject *parent = 0 )
    : Wt::WStandardItemModel( parent ),
      m_max_activity( 0.0 )
  {
  }
  
  virtual ~DecayActivityModel()
  {
  }
  
  virtual boost::any data( const Wt::WModelIndex &index,
                          int role = Wt::DisplayRole ) const
  {
    if( index.column() == 0 )
      return WStandardItemModel::data( index, role );
    
    boost::any show = headerData( index.column(), Wt::Horizontal, Wt::UserRole );
    try
    {
      const bool doShow = boost::any_cast<bool>( show );
      if( !doShow )
        return boost::any();
    }catch( std::exception &e )
    {
      cerr << "DecayActivityModel::data(...) unexpectedly caught: " << e.what()
           << endl << index.row() << ", " << index.column()
           << endl;
    }//try / catch
    
#if( WT_VERSION > 0x3040000 )
    if( role == Wt::DisplayRole )
    {
      boost::any val = WStandardItemModel::data( index, role );
      if( val.empty() )
        return val;
      
      const double valdbl = asNumber(val);
      if( IsNan(valdbl) )
        return val;
      
      // 0.00001 chosen arbitrarily, but values this small wont show up as non-zero on the chart
      if( valdbl < 0.00001*m_max_activity )
        return boost::any( 0.0 );
      
      return val;
    }//if( role == Wt::DisplayRole )
#endif
    
    
    return WStandardItemModel::data( index, role );
  }//data(...)
  
  void setNuclide( int colum, const std::string nuc )
  {
    setHeaderData( colum, Wt::Horizontal, boost::any(nuc), Wt::UserRole+1 );
  }
  
  std::string nuclide( int colum )
  {
    if( colum <= 0 || colum >= columnCount() )
      return "";
    boost::any val = headerData( colum, Wt::Horizontal, Wt::UserRole+1 );
    try
    {
      return boost::any_cast<std::string>( val );
    }catch(...)
    {}
    
    return "";
  }
  
  bool showSeries( int colum ) const
  {
    if( colum <= 0 || colum >= columnCount() )
      return true;
    
    boost::any val = headerData( colum, Wt::Horizontal, Wt::UserRole );
    try
    {
      return boost::any_cast<bool>( val );
    }catch(...)
    {}
    
    return true;
  }//bool showSeries( int colum ) const

#if( WT_VERSION > 0x3040000 )
  void updateMaxActivity()
  {
    m_max_activity = 0.0;
    
    for( int row = 0; row < rowCount(); ++row )
    {
      for( int col = 1; col < columnCount(); ++col )
      {
        if( !showSeries( col ) )
          continue;
        
        boost::any val = WStandardItemModel::data( index(row, col) );
        if( val.empty() )
          continue;
        
        const double valdbl = asNumber(val);
        if( IsNan(valdbl) )
          continue;
        
        m_max_activity = std::max( m_max_activity, valdbl );
      }
    }
  }//void updateMaxActivity()
#endif
  
  void setShowSeries( int colum, bool show )
  {
    const bool oldval = showSeries( colum );
    if( oldval == show )
      return;
    
    if( colum < 1 || colum >= columnCount() )
      return;
    
    setHeaderData( colum, Wt::Horizontal, boost::any(show), Wt::UserRole );
    
#if( WT_VERSION > 0x3040000 )
    updateMaxActivity();
#endif
    
    if( rowCount() )
    {
      WModelIndex left = index( 0, colum );
      WModelIndex right = index( rowCount()-1, colum );
      dataChanged().emit( left, right );
    }//if( rowCount() )
  }//bool setShowSeries( int colum, bool show ) const
};//class DecayActivityModel


class DecayActivityChart : public Wt::Chart::WCartesianChart
{
protected:
  double m_widthInPixels, m_heightInPixels;
  std::string m_xAxisTitle;
  bool m_phone;
  std::shared_ptr<const ColorTheme> m_colorTheme;  //TODO: actually use this information
  WBrush m_chartMarginBrush;
  WPen m_textPen;
  
public:
  DecayActivityChart(  Wt::WContainerWidget *parent = NULL  )
  : Wt::Chart::WCartesianChart( parent ),
    m_widthInPixels(250),
    m_heightInPixels(250),
    m_phone( false )
  {
    setStyleClass( "DecayActivityChart" );
    //setPreferredMethod( WPaintedWidget::Method::InlineSvgVml );
    setSeriesSelectionEnabled(true);
    //Signal<const WDataSeries *, WPointF>& seriesSelected() { return seriesSelected_; }
    //setFollowCurve(const WDataSeries *series);
  }//DecayActivityChart( constructor )
  
  virtual ~DecayActivityChart()
  {
  }
  
  std::shared_ptr<const ColorTheme> colorTheme()
  {
    return m_colorTheme;
  }
  
  void setIsPhone()
  {
    m_phone = true;
  }
  
  void setColorTheme( std::shared_ptr<const ColorTheme> theme )
  {
    m_colorTheme = theme;
    if( !m_colorTheme )
      return;  //Should reset colors to default, but whatever for now since we only ever call this funciton once
    
    if( theme->theme_name.toUTF8().find("Default") != string::npos )
    {
    }else
    {
      //Need to customize 'WColor seriesColor( int i )' funciton to use theme->referenceLineColor colors
      //theme->referenceLineColor
      //theme->spectrumChartBackground;
      //theme->spectrumAxisLines;
      
      if( theme->spectrumChartText.isDefault() )
      {
        m_textPen = WPen( WColor(GlobalColor::black) );
      }else
      {
        m_textPen = WPen(theme->spectrumChartText);
      }
      setTextPen( m_textPen );
      
      if( theme->spectrumAxisLines.isDefault() )
      {
        axis(Chart::XAxis).setPen( WPen(GlobalColor::black) );
        axis(Chart::YAxis).setPen( WPen(GlobalColor::black) );
      }else
      {
        axis(Chart::XAxis).setPen( WPen(theme->spectrumAxisLines) );
        axis(Chart::YAxis).setPen( WPen(theme->spectrumAxisLines) );
      }
      
      if( theme->spectrumChartMargins.isDefault() )
        m_chartMarginBrush = WBrush();
      else
        m_chartMarginBrush = WBrush(theme->spectrumChartMargins);
      
      if( theme->spectrumChartBackground.isDefault() )
        setBackground( Wt::NoBrush );
      else
        setBackground( WBrush(theme->spectrumChartBackground) );
    }
  }//setColorTheme()
  
  
  void setXAxisTitle( const std::string &xtitle )
  {
    m_xAxisTitle = xtitle;
    if( !m_phone )
      axis(Chart::XAxis).setTitle( xtitle );
  }
  
  int paintedWidth() const { return static_cast<int>(m_widthInPixels); }
  int paintedHeight() const { return static_cast<int>(m_heightInPixels); }
  
  void showGridLines( const bool draw = true )
  {
    Chart::WCartesianChart::axis(Chart::XAxis).setGridLinesEnabled( draw );
    Chart::WCartesianChart::axis(Chart::YAxis).setGridLinesEnabled( draw );
  }//void showGridLines( const bool draw )
  
  virtual void paintEvent( Wt::WPaintDevice *paintDevice )
  {
    m_widthInPixels  = paintDevice->width().toPixels();
    m_heightInPixels = paintDevice->height().toPixels();
    WCartesianChart::paintEvent( paintDevice );
  }
  
  virtual void paint( Wt::WPainter& painter,
                     const Wt::WRectF& rectangle = Wt::WRectF() ) const
  {
    Chart::WCartesianChart::paint( painter, rectangle );
    
    painter.save();
    //TODO 20110325: for chrome the below seems to have no effect, however
    //               it works on FireFox
    painter.rotate( -90 );
    painter.setPen( m_textPen );
    painter.setBrush( WBrush() );
    const std::string ytitle = axis(Chart::YAxis).title().toUTF8();
    if( ytitle != "" )
    {
      const double height = painter.window().height();
      painter.drawText( -0.45*height, 0.0, 0.2, 0.1,
                       AlignCenter, axis(Chart::YAxis).title() );
    }//if( ytitle != "" )
    
    painter.restore();
    
    const double height = painter.window().height();
    const double width = painter.window().width();
    
    if( m_phone )
    {
      //We will make the x-axis title "compact".  This has the downside that
      //  since we arent manually drawing x-axis labels (Wt is doing it), we
      //  dont know if we
      const double w = 7*m_xAxisTitle.size();
      const double x = width - w - 5 - plotAreaPadding(Wt::Right);
      const double y = height-15;
      WBrush fillbrush = background();
      if( fillbrush.color().isDefault() )
        fillbrush = WBrush( WColor(120,120,120) );
      
      painter.fillRect( x-5, y-2.5, w+12, 17, fillbrush );
      painter.setPen( m_textPen );
      painter.drawText( x, y, w, 12, AlignMiddle, m_xAxisTitle );
    }//if( m_phone )
    
    
    //Color theme support code not complete or tested!
    //  Also the MakeDrfChart essentially copies this code, so if you improve here, update there too
    //  (or refactor into common base class)
    auto plotArea = [&]() -> WRectF {
      int w, h;
      if( rectangle.isNull() || rectangle.isEmpty() )
      {
        w = static_cast<int>( width );
        h = static_cast<int>( height );
      }else
      {
        w = static_cast<int>( rectangle.width() );
        h = static_cast<int>( rectangle.height() );
      }
      
      const int padLeft = plotAreaPadding(Left);
      const int padRight = plotAreaPadding(Right);
      const int padTop = plotAreaPadding(Top);
      const int padBottom = plotAreaPadding(Bottom);
      
      WRectF area;
      if( orientation() == Vertical )
        area = WRectF( padLeft, padTop, std::max(10, w - padLeft - padRight), std::max(10, h - padTop - padBottom) );
      else
        area = WRectF( padTop, padRight, std::max(10, w - padTop - padBottom), std::max(10, h - padRight - padLeft) );
      
      return area;
    };//plotArea
    
    auto paintMargins = [&](){
      auto area = plotArea();
      const double rx = area.x();
      const double ry = area.y();
      const double rw = area.width();
      const double rh = area.height();
      if( ry > 0 )  //Top strip
        painter.fillRect( hv(WRectF(rx,0,rw,ry)), m_chartMarginBrush );
      if( (ry+rh) < rectangle.height() ) //Bottom strip
        painter.fillRect( hv(WRectF(rx,ry+rh,rw,rectangle.height()-ry)), m_chartMarginBrush );
      if( rx > 0 )
        painter.fillRect( hv(WRectF(0,0,rx,rectangle.height())), m_chartMarginBrush );
      if( (rx+rw) < rectangle.width() )
        painter.fillRect( hv(WRectF(rx+rw,0,rectangle.width()-rw-rx,rectangle.height())), m_chartMarginBrush );
    };
    
    if( background().style() != NoBrush && m_chartMarginBrush.style() != NoBrush )
    {
      paintMargins();
      painter.fillRect( hv(plotArea()), background() );
    }else if( background().style() != NoBrush )
    {
      painter.fillRect( hv(rectangle), background() );
    }else if( m_chartMarginBrush.style() != NoBrush )
    {
      paintMargins();
    }
  }//paint( ... )
};//class DecayActivityChart




class DateLengthCalculator : public WContainerWidget
{
  DecayActivityDiv *m_activityDiv;
  WGridLayout *m_layout;
  WDateEdit *m_begindate;
  WDateEdit *m_enddate;
  WLineEdit *m_duration;
  //WPushButton *m_updateParent;
  WText *m_error;
  WContainerWidget *m_info;
  
  
  public:
  DateLengthCalculator( DecayActivityDiv *activityDiv, WContainerWidget *parent )
   : WContainerWidget( parent ),
     m_activityDiv( activityDiv ),
     m_layout( NULL ),
     m_begindate( NULL ),
     m_enddate( NULL ),
     m_duration( NULL ),
     //m_updateParent( nullptr ),
     m_error( NULL ),
     m_info( NULL )
  {
    assert( m_activityDiv );
    
    int ndaysBack = 364;
    const double currentDur = activityDiv->timeToDisplayTill();
    if( currentDur > 0 && currentDur < 100*PhysicalUnits::year )
      ndaysBack = (int)floor(0.5 + currentDur / PhysicalUnits::day);
    
    
    m_layout = new WGridLayout( this );
    
    WLabel *label = new WLabel( "Initial Date" );
    m_layout->addWidget( label, 0, 0, AlignMiddle );
    
    m_begindate = new WDateEdit();
    label->setBuddy( m_begindate );
    m_begindate->setFormat( "MM/dd/yyyy" );
    m_begindate->setPlaceholderText( "mm/dd/yyy" );
    try
    {
      m_begindate->setDate(Wt::WDate::currentServerDate().addDays(-ndaysBack) );
    }catch(...)
    {
    }
    
    m_layout->addWidget( m_begindate, 0, 1, 1, 2, AlignLeft );
    
    label = new WLabel( "End Date" );
    m_layout->addWidget( label, 1, 0, AlignMiddle );
    m_enddate = new WDateEdit();
    label->setBuddy( m_enddate );
    m_enddate->setFormat( "MM/dd/yyyy" );
    m_enddate->setPlaceholderText( "mm/dd/yyy" );
    // \TODO: consider setting to current foreground spectrum date
    m_enddate->setDate( Wt::WDate::currentServerDate() );
    m_layout->addWidget( m_enddate, 1, 1, 1, 2, AlignLeft );
    
    m_begindate->changed().connect( this, &DateLengthCalculator::dateChanged );
    m_enddate->changed().connect( this, &DateLengthCalculator::dateChanged );
    
    label = new WLabel( "Time Span" );
    m_layout->addWidget( label, 2, 0, AlignMiddle );
    
    m_duration = new WLineEdit();
    
    m_duration->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_duration->setAttributeValue( "autocorrect", "off" );
    m_duration->setAttributeValue( "spellcheck", "off" );
#endif
    label->setBuddy( m_duration );
    
    WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalPosOrNegRegex(), m_duration );
    validator->setFlags(Wt::MatchCaseInsensitive);
    m_duration->setValidator(validator);
    m_duration->changed().connect( this, &DateLengthCalculator::durationChanged );
    m_duration->enterPressed().connect( this, &DateLengthCalculator::durationChanged );
    m_layout->addWidget( m_duration, 2, 1, AlignLeft );
    
    //m_updateParent = new WPushButton( "Update Chart");
    //m_updateParent->clicked().connect( this, &DateLengthCalculator::pushCurrentToParent );
    //m_layout->addWidget( m_updateParent, 2, 2, AlignCenter );
    m_layout->addWidget( new WContainerWidget(), 2, 2, AlignCenter );
    
    m_error = new WText( "" );
    m_error->setInline( false );
    m_error->hide();
    m_layout->addWidget( m_error, 3, 0, 1, 3 );
    
    m_info = new WContainerWidget( this );
    m_info->hide();
    m_layout->addWidget( m_info, 4, 0, 1, 3 );
    
    m_layout->setRowStretch( 3, 1 );
    m_layout->setRowStretch( 4, 1 );
    
    m_layout->setColumnStretch( 2, 1 );
  }//DateLengthCalculator
  
  /** Returns 0.0 on error.
   Return a negative time period if we want to calculate an initial activity.
   */
  double getValidatedTimeSpan()
  {
    double halfLife = 0.0;
    const SandiaDecay::Nuclide *nuc = NULL;
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    const std::vector<DecayActivityDiv::Nuclide> &nucs = m_activityDiv->m_nuclides;
    if( !nucs.empty() )
      nuc = db->nuclide( nucs[0].z, nucs[0].a, nucs[0].iso );
    if( nuc )
      halfLife = nuc->halfLife;
    
    //Let be optimistic things will be okay
    m_error->hide();
    m_error->setText("");
    //if( m_updateParent )
    //  m_updateParent->enable();
    
    try
    {
      const std::string durtxt = m_duration->text().toUTF8();
      const bool validBeginDate = (m_begindate->validate() == Wt::WValidator::Valid);
      const bool validEndDate = (m_enddate->validate() == Wt::WValidator::Valid);
      bool duration_valid = false;
      double duration = 0.0;
      try
      {
        duration = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( durtxt, halfLife );
//        if( duration < 1*PhysicalUnits::second || duration > 100*PhysicalUnits::year )
  //        throw runtime_error( "Time span must be greater than one second and less than 100 years" );
        duration_valid = true;
      }catch(...)
      {
        if( !durtxt.empty() )
          m_duration->setText( "" );
      }
      
      
      if( duration_valid && validBeginDate && validEndDate )
      {
        // Should make sure duration corresponds to date range...
        return duration;
      }else if( !duration_valid )
      {
        if( !validBeginDate )
          throw runtime_error( "Begin date is invalid" );
        if( !validEndDate )
          throw runtime_error( "End date is invalid" );
        
        //if( m_begindate->date() >= m_enddate->date() )
        //  throw runtime_error( "End date must be later than start date" );
        
        const int days = m_begindate->date().daysTo(m_enddate->date());
        duration = 24.0*3600.0*days*PhysicalUnits::second;
        const string datestr = PhysicalUnitsLocalized::printToBestTimeUnits( duration, 2 );
        m_duration->setText( datestr );
        
        return duration;
      }else if( !validBeginDate && !validEndDate )
      {
        return duration;
      }else if( !validBeginDate )
      {
        const int ndays = (int)floor( 0.5 + duration/PhysicalUnits::day );
        m_begindate->setDate( m_enddate->date().addDays( -ndays ) );
        return duration;
      }else if( !validEndDate )
      {
        const int ndays = (int)floor( 0.5 + duration/PhysicalUnits::day );
        m_enddate->setDate( m_begindate->date().addDays( ndays ) );
        return duration;
      }
    }catch( std::exception &e )
    {
      m_error->show();
      m_error->setText( e.what() );
      m_duration->setText( "" );
      //if( m_updateParent )
      //  m_updateParent->disable();
      m_info->hide();
    }//try / catch
    
    return 0.0;
  }//double getValidatedTimeSpan()
  
  
  void pushCurrentToParent()
  {
    if( getValidatedTimeSpan() == 0.0 )
      return;
    
    string txt = m_duration->text().toUTF8();
    if( txt.length() && (txt[0] == '-') )
      txt = txt.substr(1);
    
    m_activityDiv->m_displayTimeLength->setText( txt );
    m_activityDiv->refreshDecayDisplay( false );
  }//void pushCurrentToParent()
  
  
  void dateChanged()
  {
    m_duration->setText( "" );
    
    if( getValidatedTimeSpan() == 0.0 )
      return;
    
    int days = m_begindate->date().daysTo(m_enddate->date());
    //days = abs(days);
    string datestr = PhysicalUnitsLocalized::printToBestTimeUnits( 24.0*3600.0*days*PhysicalUnits::second, 2 );
    m_duration->setText( datestr );
    
    pushCurrentToParent();
    updateInfo();
  }//void dateChanged()
  
  
  void durationChanged()
  {
    // We get here when the user edits the duration text
    const string durtxt = m_duration->text().toUTF8();
    if( !durtxt.empty() && (m_duration->validate() == WValidator::State::Valid) )
    {
      try
      {
        const double duration = PhysicalUnitsLocalized::stringToTimeDuration( durtxt );
        
        double mintime = 1.0 * SandiaDecay::second;
        const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
        for( const auto &n : m_activityDiv->m_nuclides )
        {
          const auto nuc = db->nuclide( n.z, n.a, n.iso );
          if( nuc )
            mintime = std::min( mintime, nuc->halfLife );
        }
        
        if( (fabs(duration) < 0.01*mintime) && m_activityDiv && m_activityDiv->m_displayTimeLength )
        {
          m_duration->setText( m_activityDiv->m_displayTimeLength->text() );
          updateInfo();
          return;
        }
        
        int64_t ndays = std::floor(fabs(duration) / PhysicalUnits::day);
        double remanderSeconds = fabs(duration) - (ndays * PhysicalUnits::day);
        if( fabs(PhysicalUnits::day - remanderSeconds) <= 1.0 )
        {
          remanderSeconds = PhysicalUnits::day - remanderSeconds;
          ndays += 1.0;
        }
        
        // \TODO: do the same thing as days, but for years to try and set an initial date
        if( duration < 0.0 )
          ndays = -ndays;

        if( (abs(ndays) < 1) || (remanderSeconds > 1.0))
        {
          m_begindate->setText( "" );
        }else if( m_enddate->validate() == WValidator::State::Valid )
        {
          m_begindate->setDate( m_enddate->date().addDays( -static_cast<int>(ndays) ) );
        }else if( m_begindate->validate() == WValidator::State::Valid )
        {
          m_enddate->setDate( m_begindate->date().addDays( static_cast<int>(ndays) ) );
        }else
        {
          m_begindate->setText( "" );
          m_enddate->setText( "" );
        }
      }catch(...)
      {
        m_begindate->setText( "" );
      }
    }//if( we have valid duration string )
    
    if( getValidatedTimeSpan() == 0.0 )
      return;
    
    pushCurrentToParent();
    updateInfo();
  }//void durationChanged()
  
  void setTimeRangeTxt( std::string txt )
  {
    const string oldtxt = m_duration->text().toUTF8();
    if( txt != oldtxt )
    {
      m_duration->setText( txt );
      if( oldtxt.empty() || (oldtxt[0]!= '-') || m_begindate->text().empty() )
        m_begindate->setText( "" );
      else
        m_enddate->setText( "" );
    }
    
    updateInfo(); //the
  }//void setTimeRangeTxt( std::string txt )
  
  
  void updateInfo()
  {
    m_info->clear();
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    
    const double timeSpan = getValidatedTimeSpan();
    const double absTimeSpan = fabs( timeSpan );
    
    if( timeSpan == 0.0 )
    {
      m_info->hide();
      return;
    }
    
    m_info->show();
    
    
    SandiaDecay::NuclideMixture *mix = m_activityDiv->m_currentMixture;
    
    int nSigFigsActivity = 0;
    
    // If we have a negative time span, we need to create our own mixture
    std::unique_ptr<SandiaDecay::NuclideMixture> local_mix;
    if( timeSpan < 0.0 )
    {
      local_mix.reset( new SandiaDecay::NuclideMixture() );
      mix = local_mix.get();
      // We'll add the nuclides in the loop over nuclides below
    }//if( timeSpan < 0.0 )
    
    if( !mix )
      return;
    
    //const string durationTxt = m_duration->text().toUTF8();
    const string durationTxt = PhysicalUnitsLocalized::printToBestTimeUnits( absTimeSpan );
    
    const std::vector<DecayActivityDiv::Nuclide> &nucs = m_activityDiv->m_nuclides;
    if( nucs.empty() )
      return;
    
    
    for( size_t i = 0; i < nucs.size(); ++i )
    {
      const DecayActivityDiv::Nuclide &nucinfo = nucs[i];
      const SandiaDecay::Nuclide *nuc = db->nuclide( nucinfo.z, nucinfo.a, nucinfo.iso );
      if( !nuc )
        continue;
      
      double entered_activity = nucinfo.activity;
      const bool useCurie = nucinfo.useCurie;
      //const bool useCurie = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      
      if( timeSpan < 0.0 )
      {
        // TODO: estimate how many half-lives is reasonable to calculate this initial activity for.
        const double decrease_factor = std::exp( timeSpan * nuc->decayConstant() );
        double initial_activity = entered_activity / decrease_factor;
        
        if( IsNan(initial_activity)
            || IsInf(initial_activity)
            || (initial_activity < std::numeric_limits<float>::epsilon()) )
        {
          m_info->hide();
          m_error->show();
          string errmsg = "Duration of " + printToBestTimeUnits(absTimeSpan) + " is to long for "
                          + nuc->symbol + " (half life " + printToBestTimeUnits(nuc->halfLife) + ")";
          m_error->setText( errmsg );
          m_duration->setText( "" );
          
          return;
        }//if( probably too many half-lives to back calculate or zero initial activity )
        
        assert( local_mix && (mix == local_mix.get()) );
        
        if( nucinfo.age > DBL_EPSILON )
          local_mix->addAgedNuclideByActivity( nuc, initial_activity, nucinfo.age );
        else
          local_mix->addNuclideByActivity( nuc, initial_activity );
      }//if( timeSpan < 0.0 )
      
      string actTxt = PhysicalUnits::printToBestActivityUnits( entered_activity, 3, useCurie );
      if (!nucinfo.activityStr.empty())
      {
        assert((entered_activity < 1.0E-6)
            || (fabs(entered_activity - PhysicalUnits::stringToActivity(nucinfo.activityStr)) < 0.001 * entered_activity));
       
        int nThisSigFig = 0;
        for (const char c : nucinfo.activityStr)
          nThisSigFig += ((c >= '0') && (c <= '9')); //(std::isdigit(c) != 0)
        nSigFigsActivity = std::max(nThisSigFig, nSigFigsActivity);

        actTxt = nucinfo.activityStr;
      }//if (!nucinfo.activityStr.empty())


      WTable *nuctbl = new WTable( m_info );
      nuctbl->setMargin( 10, Wt::Side::Top );
      
      WTableCell *cell = nullptr;
      
      const bool showNucHeader = (nucs.size() > 1);
      const int rowOffset = showNucHeader ? 1 : 0;
      
      if( showNucHeader )
      {
        nuctbl->setHeaderCount( 1 );
        cell = nuctbl->elementAt(0, 0);
        new WLabel( "For " + nuc->symbol + ":", cell );
        cell->setColumnSpan( 2 );
        cell->setAttributeValue( "style", "text-align: left;" + cell->attributeValue("style") );
      }//if( showNucHeader )
      
      cell = nuctbl->elementAt(1 + rowOffset, 0);
      cell->setVerticalAlignment( AlignmentFlag::AlignMiddle );
      new WLabel( ((timeSpan < 0.0) ? "Final Activity&nbsp;" : "Initial Activity&nbsp;"), cell );
      
      cell = nuctbl->elementAt(1 + rowOffset, 1);
      WLineEdit *activityEdit = new WLineEdit(cell);
      activityEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
      activityEdit->setAttributeValue( "autocorrect", "off" );
      activityEdit->setAttributeValue( "spellcheck", "off" );
#endif
      
      WRegExpValidator *actvalidator = new WRegExpValidator( PhysicalUnits::sm_activityRegex, activityEdit );
      actvalidator->setFlags(Wt::MatchCaseInsensitive);
      activityEdit->setValidator( actvalidator );
      activityEdit->setTextSize( 10 );
      activityEdit->setText( actTxt );
      
      auto doActivityUpdate = [this, activityEdit, nucinfo, actTxt](){
        if( activityEdit->validate() != WValidator::State::Valid )
        {
          activityEdit->setText( actTxt );
          return;
        }
        
        try
        {
          const string input_act_str = activityEdit->text().toUTF8();
          const double activity = PhysicalUnits::stringToActivity(input_act_str);

          if( activity <= 0.0 )
            throw runtime_error( "Activity must be greater than zero" );
          
          //if( timeSpan < 0.0 )
          //  activity = activity / std::exp( timeSpan * nuc->decayConstant() );
          
          if( IsNan(activity) || IsInf(activity) )
            throw runtime_error( "Invalid time-span or activity" );
          
          bool updated = false;
          for( DecayActivityDiv::Nuclide &n : m_activityDiv->m_nuclides )
          {
            if( (n.z == nucinfo.z) && (n.a == nucinfo.a) && (n.iso == nucinfo.iso)
                && fabs(n.age - nucinfo.age) < 1.0 )
            {
              updated = true;
              n.activity = activity;
              n.activityStr = input_act_str;
              n.updateTxt();
              break;
            }
          }//for( loop over nuclides )
          
          if( !updated )
            throw runtime_error( "Couldnt find source to update" );
        }catch(...)
        {
          activityEdit->setText( actTxt );
          return;
        }// try / catch
        
        m_activityDiv->refreshDecayDisplay( false );
        updateInfo();
      };//doActivityUpdate lambda
      
      activityEdit->changed().connect( std::bind( doActivityUpdate ) );
      activityEdit->enterPressed().connect( std::bind( doActivityUpdate ) );
      
      
      const string ageTxt = PhysicalUnitsLocalized::printToBestTimeUnits( nucinfo.age );
      cell = nuctbl->elementAt(2 + rowOffset, 0);
      cell->setVerticalAlignment( AlignmentFlag::AlignMiddle );
      new WLabel( "Initial Age&nbsp;", cell );
      
      cell = nuctbl->elementAt(2 + rowOffset, 1);
      WLineEdit *ageEdit = new WLineEdit(cell);
      ageEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
      ageEdit->setAttributeValue( "autocorrect", "off" );
      ageEdit->setAttributeValue( "spellcheck", "off" );
#endif

      WRegExpValidator *agevalidator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationRegex(), ageEdit );
      agevalidator->setFlags(Wt::MatchCaseInsensitive);
      ageEdit->setValidator( agevalidator );
      ageEdit->setTextSize( 10 );
      ageEdit->setText( ageTxt );
      
      auto doAgeUpdate = [this, ageEdit, ageTxt, nucinfo](){
        if( ageEdit->validate() != WValidator::State::Valid )
        {
          ageEdit->setText( ageTxt );
          return;
        }
        
        try
        {
          double age = PhysicalUnitsLocalized::stringToTimeDuration( ageEdit->text().toUTF8() );
          if( age < 0.0 )
            throw runtime_error( "Initial age must be zero or more" );
          
          bool updated = false;
          for( DecayActivityDiv::Nuclide &n : m_activityDiv->m_nuclides )
          {
            if( (n.z == nucinfo.z) && (n.a == nucinfo.a) && (n.iso == nucinfo.iso)
               && fabs(n.age - nucinfo.age) < 1.0 )
            {
              updated = true;
              n.age = age;
              n.updateTxt();
              break;
            }
          }//for( loop over nuclides )
          
          if( !updated )
            throw runtime_error( "Couldnt find source to update" );
        }catch(...)
        {
          ageEdit->setText( ageTxt );
          return;
        }// try / catch
        
        m_activityDiv->refreshDecayDisplay( false );
        updateInfo();
      };//doActivityUpdate lambda
      
      ageEdit->changed().connect( std::bind( doAgeUpdate ) );
      ageEdit->enterPressed().connect( std::bind( doAgeUpdate ) );
    }//for( size_t i = 0; i < nucs.size(); ++i )
    
    auto checkUseCuries = [&nucs,db]( const SandiaDecay::NuclideActivityPair &nap ) -> bool {
      
      for( size_t nucn = 0; nucn < nucs.size(); ++nucn )
      {
        const vector<const SandiaDecay::Nuclide *> parents = nap.nuclide->forebearers();
        const SandiaDecay::Nuclide *nuc = db->nuclide( nucs[nucn].z, nucs[nucn].a, nucs[nucn].iso );
        if( std::find(parents.begin(),parents.end(),nuc) != parents.end() )
          return nucs[nucn].useCurie;
      }//for( size_t nucn = 0; nucn < nucs.size(); ++nucn )
      
      return true;
    };//checkUseCuries
    
    
    if( timeSpan < 0.0 )
    {
      WText *line = new WText( "&nbsp;", m_info );
      line->setInline( false );
      line = new WText( "Before " + durationTxt + ":", m_info );
      line->setInline( false );
      
      WTable *beforeTable = new WTable( m_info );
      beforeTable->addStyleClass( "BeforeDecayNucInfoTable DecayedNucInfoTable" );
      beforeTable->setHeaderCount( 1, Wt::Horizontal );
      beforeTable->elementAt(0, 0)->addWidget( new WLabel("Nuclide") );
      beforeTable->elementAt(0, 1)->addWidget( new WLabel("Activity") );
      beforeTable->elementAt(0, 2)->addWidget( new WLabel("Mass") );
      
      const vector<SandiaDecay::NuclideActivityPair> activities = mix->activity( 0.0 );
      const vector<SandiaDecay::NuclideNumAtomsPair> numatoms = mix->numAtoms( 0.0 );
      
      for( size_t actnum = 0, row_num = 0; actnum < activities.size(); ++actnum )
      {
        const SandiaDecay::NuclideActivityPair &nap = activities[activities.size()-1-actnum];
        if( !nap.nuclide )
          continue;
        
        bool is_orig = false;
        const int nInitialNuc = mix->numInitialNuclides();
        for( int i = 0; !is_orig && (i < nInitialNuc); ++i )
        {
          auto n = mix->initialNuclide(i);
          is_orig = (n && (n == nap.nuclide));
        }
          
        if( !is_orig )
          continue;
        
        row_num += 1;
        
        //Figure out if user inputed activity in Ci or Bq for this nuclide
        const bool useCuries = checkUseCuries( nap );
        
        beforeTable->elementAt(static_cast<int>(row_num),0)->addWidget( new WText(nap.nuclide->symbol) );
        
        const string actstr = PhysicalUnits::printToBestActivityUnits( nap.activity, 3, useCuries, SandiaDecay::becquerel );
        beforeTable->elementAt(static_cast<int>(row_num),1)->addWidget( new WText(actstr) );
        
        const SandiaDecay::NuclideNumAtomsPair *natomp = NULL;
        for( size_t a = 0; !natomp && a < numatoms.size(); ++a )
        if( numatoms[a].nuclide == nap.nuclide )
          natomp = &(numatoms[a]);
        
        if( natomp )
        {
          const double num_atoms_in_gram = nap.nuclide->atomsPerGram();
          const double ngrams = natomp->numAtoms / num_atoms_in_gram;
          const string massstr = PhysicalUnits::printToBestMassUnits( ngrams*PhysicalUnits::gram );
          beforeTable->elementAt(static_cast<int>(row_num),2)->addWidget( new WText(massstr) );
        }
      }//for( size_t i = 0; i < activities.size(); ++i )
    }//if( timeSpan < 0.0 )
    
    
    {
      WText *line = new WText( "&nbsp;", m_info );
      line->setInline( false );
      line = new WText( "After " + durationTxt + ":", m_info );
      line->setInline( false );
    }
    
    WTable *infotable = new WTable( m_info );
    infotable->addStyleClass( "DecayedNucInfoTable" );
    infotable->setHeaderCount( 1, Wt::Horizontal );
    infotable->elementAt(0, 0)->addWidget( new WLabel("Nuclide") );
    infotable->elementAt(0, 1)->addWidget( new WLabel("Activity") );
    infotable->elementAt(0, 2)->addWidget( new WLabel("Mass") );
    
    //infotable->elementAt(0, 3)->addWidget( new WLabel("Gammas/s") );  //could have mouseover showing intensity of each energy
    //infotable->elementAt(0, 4)->addWidget( new WLabel("Mass Frac") );
    //infotable->elementAt(0, 5)->addWidget( new WLabel("Act. Frac") );
    
    const vector<SandiaDecay::NuclideActivityPair> activities = mix->activity( absTimeSpan );
    const vector<SandiaDecay::NuclideNumAtomsPair> numatoms = mix->numAtoms( absTimeSpan );
        
    // What we actually want is number of sig-figs, but since the print function
    //  take number of digits after decimal, we'll just assume a single significant 
    //  digit before the decimal, and then print `nSigFigsActivity` digits after
    //  (i.e., one more sig-fig than input)
    const int maxDecimal = std::max(2, nSigFigsActivity);
    //activityEdit->tex
    
    for( int actnum = 0, row_num = 1; actnum < static_cast<int>(activities.size()); ++actnum )
    {
      const SandiaDecay::NuclideActivityPair &nap = activities[activities.size()-1-actnum];
      if( !nap.nuclide )
        continue;
      
      //Figure out if user inputed activity in Ci or Bq for this nuclide
      const bool useCuries = checkUseCuries( nap );
      
      infotable->elementAt(row_num,0)->addWidget( new WText(nap.nuclide->symbol) );
      
      string actstr = PhysicalUnits::printToBestActivityUnits( nap.activity, maxDecimal, useCuries );
      if( IsInf(nap.nuclide->halfLife) )
        actstr = "stable";
      infotable->elementAt(row_num,1)->addWidget( new WText(actstr) );
      
      
      const SandiaDecay::NuclideNumAtomsPair *natomp = NULL;
      for( size_t a = 0; !natomp && a < numatoms.size(); ++a )
        if( numatoms[a].nuclide == nap.nuclide )
          natomp = &(numatoms[a]);
      
      if( natomp )
      {
        const double num_atoms_in_gram = nap.nuclide->atomsPerGram();
        const double ngrams = natomp->numAtoms / num_atoms_in_gram;
        const string massstr = PhysicalUnits::printToBestMassUnits( ngrams*PhysicalUnits::gram, maxDecimal);
        infotable->elementAt(static_cast<int>(row_num),2)->addWidget( new WText(massstr) );
      }
      
      row_num += 1;
    }//for( size_t i = 0; i < activities.size(); ++i )
  }//void updateInfo()
  
};//class DateLengthCalculator

namespace
{
  enum YAxisType
  {
    ActivityAxis,
    GammasAxis,
    BetasAxis,
    AlphasAxis,
    NumYAxisType
  };//enum YAxisType
  
  
  //We will use the HSV space to select sequential colors for the data series,
  //  that are reasonably distinct, and then translate back into RGB.
  //  Bassically we'll hold the saturation and value constant, and increment
  //  the hue by one over the golden ratio (and then modulo this by 1.0).
  //see:
  //  http://martin.ankerl.com/2009/12/09/how-to-create-random-colors-programmatically/
  void hsv_to_rgb( const float h, const float s, const float v,
                  uint8_t &r_ans, uint8_t &g_ans, uint8_t &b_ans )
  {
    // h, s, and v should all vary between [0 ,1)
    unsigned int h_i = static_cast<unsigned int>( h * 6.0f );
    float f = h*6.0f - h_i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - f*s);
    float t = v * (1.0f - (1.0f - f) * s);
    
    float r, g, b;
    
    if( h_i == 0 )    { r = v; g = t; b = p; }
    else if( h_i==1 ) { r = q; g = v; b = p; }
    else if( h_i==2 ) { r = p; g = v; b = t; }
    else if( h_i==3 ) { r = p; g = q; b = v; }
    else if( h_i==4 ) { r = t; g = p; b = v; }
    else              { r = v; g = p; b = q; }
    
    r_ans = static_cast<uint8_t>( r*256 );
    g_ans = static_cast<uint8_t>( g*256 );
    b_ans = static_cast<uint8_t>( b*256 );
  }//void hsv_to_rgb( uint8_t &h_r, uint8_t &s_g, uint8_t &v_b )
  
  WColor seriesColor( int i )
  {
    uint8_t r, g, b;
    
    const float golden_ratio_conjugate = 0.618033988749895f;
    float hue = fabs(0.9f + i*golden_ratio_conjugate);
    hue -= floor( hue );
    const float saturation = 0.88f;
    const float value = 0.98f;
    hsv_to_rgb( hue, saturation, value, r, g, b );
    return WColor( r, g, b );
  }//WColor seriesColor( int i )
  
  //ChartAndLegendHolder: a helper class so we can display a legend over the
  //  chart using a WContainerWidget for the legend.  This is kind a non-ideal
  //  solution, but I couldnt find any other ones that worked super easy.
  class ChartAndLegendHolder : public WContainerWidget
  {
  public:
    ChartAndLegendHolder( DecayActivityChart *chart, WContainerWidget *legend,
                          const bool isPhone, WContainerWidget *parent = 0 )
    : WContainerWidget( parent ),
      m_chart( chart ),
      m_legend( legend )
    {
      const string chartref = "$('#" + m_chart->id() + "').get(0)";
      const string legref = "$('#" + legend->id() + "')";
      
      addWidget( m_chart );
      addWidget( m_legend );
      
      const string extraheight = isPhone ? "20" : "100";
      
      const string js =
      """function(self, w, h,layout) {"
      ""  "if (!self.wtWidth || self.wtWidth!=w "
      ""      "|| !self.wtHeight || self.wtHeight!=h) {"
      ""    "self.wtWidth=w; self.wtHeight=h;"
      ""    "self.style.height=h + 'px';"
      ""  "}"
      ""  "var child = $('#" + m_chart->id() + "').get(0);"
      ""  "if(child && child.wtResize){"
      ""     "child.wtResize(child,w,h,layout);"
      ""  "}else{ console.log('ChartAndLegendHolder: couldnt resize chart, please look into.'); }"
      ""  + legref + ".css('max-height',(h-" + extraheight + ") +'px');"
      "" "}";
      setJavaScriptMember( WT_RESIZE_JS, js );
    }//ChartAndLegendHolder constructor
    
    
    DecayActivityChart *m_chart;
    WContainerWidget *m_legend;
  };//class ChartAndLegendHolder
  
  
  class DecayCsvResource : public Wt::WResource
  {
    //ToDo: Decouple from DecayCsvResource class and put in own source/header,
    //      or maybe with CsvDownloadGui
  public:
    DecayCsvResource( DecayActivityDiv *display )
      : WResource( display ),
        m_display( display ),
        m_numRows( 400 ),
        m_timeSpan( 0 ),
        m_activity( true ),
        m_xrays( false ),
        m_gammas( false ),
        m_alphas( false ),
        m_betas( false ),
        m_sessionId( wApp->sessionId() ),
        m_app( wApp )
    {
      if( m_display )
        m_timeSpan = m_display->timeToDisplayTill();
      
      //When a download is complete, we will delete the dialog, but jic something
      //  wacky (and unexpected) is going on, we'll wrap it so that it wont
      //  get called if the DecayActivityDiv has been deleted
      //  (I'm a little adverse to relying on a 'modal' user effect to enforce
      //   object lifetimes)
      m_cleanup = wApp->bind( boost::bind( &DecayActivityDiv::deleteCsvDownloadGuiTriggerUpdate, m_display ) );
    }
    
    virtual ~DecayCsvResource()
    {
      beingDeleted();
    }
    
    void setTimeSpan( const double span )    { m_timeSpan = span; }
    void setGiveActivities( const bool give ){ m_activity = give; }
    void setGiveXrays( const bool give )     { m_xrays    = give; }
    void setGiveGammas( const bool give )    { m_gammas   = give; }
    void setGiveAlphas( const bool give )    { m_alphas   = give; }
    void setGiveBetas( const bool give )     { m_betas    = give; }
    void setNumberRows( const int nrow ){
      if( nrow < 1 || nrow > 10001 )
        throw runtime_error( "DecayCsvResource(): Invalid number of points" );
      m_numRows  = static_cast<size_t>(nrow);
    }
    
  private:
    virtual void handleRequest( const Wt::Http::Request &request,
                                Wt::Http::Response &response )
    {
      WApplication::UpdateLock lock( m_app );
      
      if( !lock )
      {
        log("error") << "Failed to WApplication::UpdateLock in DecayCsvResource.";
        
        response.out() << "Error grabbing application lock to form DecayCsvResource resource; please report to InterSpec@sandia.gov.";
        response.setStatus(500);
        assert( 0 );
        
        return;
      }//if( !lock )
      
      
      const string eol_char = "\r\n"; //for windows - could potentially cosutomize this for the users operating system
      
      try
      {
        SandiaDecay::NuclideMixture *mixture = m_display->m_currentMixture;
        if( !mixture || !m_display )
          throw runtime_error( "No data avalaiable" );
    
        const int nparents = mixture->numInitialNuclides();
        const int nchildnucs = mixture->numSolutionNuclides();
        if( nparents < 1 || nchildnucs < 1 )
          throw runtime_error( "No nuclides selected" );
        
        double timeSpan = m_timeSpan;
        if( timeSpan <= 0.0 )
        {
          for( int i = 0; i < nparents; ++i )
          {
            auto nuc = mixture->initialNuclide(i);
            if( nuc && !IsInf(nuc->halfLife) )
              timeSpan = std::max( timeSpan, 5.0*nuc->halfLife );
          }
        }
      
        if( timeSpan <= 0 )
          throw runtime_error( "Ivalid time range selcted" );  //shouldnt ever happed
        
        const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      
        string name;
        bool use_curies = false;
        for( const DecayActivityDiv::Nuclide &n : m_display->m_nuclides )
        {
          const SandiaDecay::Nuclide *nuc = db->nuclide( n.z, n.a, n.iso );
          if( nuc && !IsInf(nuc->halfLife) )  //dont print out activity for stable nuclides.
          {
            use_curies = (use_curies || n.useCurie);
            name += (!name.empty() ? "_" : "") + nuc->symbol;
            if( n.age > std::numeric_limits<float>::epsilon() )
              name += remove_space_copy( PhysicalUnits::printToBestTimeUnits(n.age) );
          }
        }
        
        const double act_unit = use_curies ? PhysicalUnits::curie : PhysicalUnits::becquerel;
        const string act_unit_str = use_curies ? "ci" : "bq";
        
        const auto bestTimeUnit = PhysicalUnits::bestTimeUnit( timeSpan );
        
        const double time_unit = bestTimeUnit.second;
        const string time_unit_str = bestTimeUnit.first;
        
        const string timestr = PhysicalUnits::printToBestTimeUnits(timeSpan);
        
        name += "_over_" + remove_space_copy(timestr) + ".csv";
  
        suggestFileName( name, WResource::Attachment );
        response.setMimeType( "text/csv" );
      
        //const bool mutliple_types = (int(m_activity) + int(m_xrays) + int(m_gammas) + int(m_alphas) + int(m_betas)) > 1;
        const auto order = SandiaDecay::NuclideMixture::OrderByEnergy;
        
        const vector<SandiaDecay::NuclideTimeEvolution> &decayEvolutions = mixture->decayedToNuclidesEvolutions();
        
        response.out() << "Time (" << time_unit_str << ")";
        if( m_activity )
        {
          for( const auto &evo : decayEvolutions )
          {
            if( evo.nuclide && !IsInf(evo.nuclide->halfLife) )  //dont print out activity for stable nuclides.
              response.out() << "," << evo.nuclide->symbol << " (" << act_unit_str << ")";
          }
        }//if( m_activity )
        
        if( m_xrays )
        {
          for( auto &x : mixture->xrays(timeSpan/10.0, order) )
            response.out() << "," << x.energy << " keV (xray/s)";
        }//if( m_xrays )
        
        if( m_gammas )
        {
          for( auto &x : mixture->gammas(timeSpan/10.0, order, true) )
            response.out() << "," << x.energy << " keV (gamma/s)";
        }//if( m_gammas )
        
        if( m_alphas )
        {
          for( auto &x : mixture->alphas(timeSpan/10.0, order) )
            response.out() << "," << x.energy << " keV (alpha/s)";
        }//if( m_alphas )
        
        if( m_betas )
        {
          for( auto &x : mixture->betas(timeSpan/10.0, order) )
            response.out() << "," <<  x.energy << " keV (b-/s)";
          for( auto &x : mixture->betaPlusses(timeSpan/10.0, order) )
            response.out() << "," << x.energy << " keV (b+/s)";
        }//if( m_betas )
        
        response.out() << eol_char;
      
        const double numerator = (m_numRows > 1 ? (m_numRows-1.0) : 1.0);
        
        for( size_t row = 0; row < m_numRows; ++row )
        {
          // If you are only requesting one row, make it the full time span time, not t=0.
          const double time_now = m_timeSpan * ((m_numRows > 1) ? row : 1u) / numerator;
          
          response.out() << (time_now/time_unit);
          
          if( m_activity )
          {
            for( auto &x : mixture->activity(time_now) )
              if( x.nuclide && !IsInf(x.nuclide->halfLife) )
              response.out() << "," << (x.activity/act_unit);
          }//if( m_activity )
          
          if( m_xrays )
          {
            for( auto &x : mixture->xrays(time_now, order) )
              response.out() << "," << x.numPerSecond;
          }//if( m_xrays )
          
          if( m_gammas )
          {
            for( auto &x : mixture->gammas(time_now, order, true) )
              response.out() << "," << x.numPerSecond;
          }//if( m_gammas )
          
          if( m_alphas )
          {
            for( auto &x : mixture->alphas(time_now, order) )
              response.out() << "," << x.numPerSecond;
          }//if( m_alphas )
          
          if( m_betas )
          {
            for( auto &x : mixture->betas(time_now, order) )
              response.out() << "," << x.numPerSecond;
            for( auto &x : mixture->betaPlusses(time_now, order) )
              response.out() << "," << x.numPerSecond;
          }//if( m_betas )
        
          response.out() << eol_char;
        }//for( int row = 0; row < nrow; ++row )
      }catch( std::exception &e )
      {
        suggestFileName( "no_data.csv", WResource::Attachment );
        response.out() << eol_char << e.what() << eol_char;
      }
      
      WServer::instance()->post( m_sessionId, m_cleanup );
    }//handleRequest(...)
    
    DecayActivityDiv *m_display;
    size_t m_numRows;
    double m_timeSpan;
    bool m_activity;
    bool m_xrays;
    bool m_gammas;
    bool m_alphas;
    bool m_betas;
    std::string m_sessionId;
    Wt::WApplication *m_app; //it looks like WApplication::instance() will be valid in handleRequest, but JIC
    boost::function<void()> m_cleanup;
  };//class DecayCsvResource : public Wt::WResource
}//namespace


class CsvDownloadGui : public AuxWindow
{
  DecayActivityDiv *m_parent;
  DecayCsvResource *m_csvResouce;
  WLineEdit *m_ageEdit;
  
public:
  
  CsvDownloadGui( DecayActivityDiv *parent )
  : AuxWindow( "CSV Export", (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal) | AuxWindowProperties::TabletNotFullScreen | AuxWindowProperties::DisableCollapse | AuxWindowProperties::SetCloseable) ),
  m_parent( parent ),
  m_csvResouce( nullptr ),
  m_ageEdit( nullptr )
  {
    assert( m_parent );

    setResizable( false );
    finished().connect( m_parent, &DecayActivityDiv::deleteCsvDownloadGui );
    
    m_csvResouce = new DecayCsvResource( parent );
    
    const double displayTime = parent->timeToDisplayTill();
    const string timeSpanStr = PhysicalUnitsLocalized::printToBestTimeUnits( displayTime );
    
    m_csvResouce->setTimeSpan( displayTime );
    m_csvResouce->setNumberRows( 400 );
    m_csvResouce->setGiveActivities( true );
    m_csvResouce->setGiveXrays( false );
    m_csvResouce->setGiveGammas( false );
    m_csvResouce->setGiveAlphas( false );
    m_csvResouce->setGiveBetas( false );
    
    WTable *table = new WTable( contents() );
    
    WTableCell *el = table->elementAt(0,0);
    table->setMargin( 8, Wt::Top | Wt::Left | Wt::Right );
    WLabel *label = new WLabel( "Time Span:", el );
    label->addStyleClass( "DecayChartTimeSpanLabel" );
    el->setVerticalAlignment( Wt::AlignMiddle );
    label->setMargin( 2, Wt::Right );
    
    el = table->elementAt(0,1);
    m_ageEdit = new WLineEdit( timeSpanStr, el );
    m_ageEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_ageEdit->setAttributeValue( "autocorrect", "off" );
    m_ageEdit->setAttributeValue( "spellcheck", "off" );
#endif
    label->setToolTip( "Enter time period you would like decay data for."
                       " Ex. '5.2 y', '52.3 s', '00:01:2.1', or '3.2d 15h'"
                       " would al be valid time periods."
                       " You can also enter time periods such as '3.2hl' or"
                       " '3.2 halflife' for multiples of the first parent"
                       " nuclides half life.");
    label->setBuddy( m_ageEdit );
    WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), this );
    validator->setFlags( Wt::MatchCaseInsensitive );
    m_ageEdit->setValidator(validator);
    m_ageEdit->changed().connect( this, &CsvDownloadGui::updateTimeSpan );
    m_ageEdit->blurred().connect( this, &CsvDownloadGui::updateTimeSpan );
    m_ageEdit->enterPressed().connect( this, &CsvDownloadGui::updateTimeSpan );
    
    el = table->elementAt(0,2);
    label = new WLabel( "Num Rows", el );
    el->setVerticalAlignment( Wt::AlignMiddle );
    label->setMargin( 20, Wt::Left );
    label->setMargin( 2, Wt::Right );
    el = table->elementAt(0,3);
    WSpinBox *numPointsSB = new WSpinBox( el );
    label->setBuddy( numPointsSB );
    label->setToolTip( "Number of time points (rows) to include data for." );
    numPointsSB->setValue( 400 );
    numPointsSB->setMinimum(1);
    numPointsSB->setMaximum( 10000 );
    numPointsSB->valueChanged().connect( m_csvResouce, &DecayCsvResource::setNumberRows );
    
    table = new WTable( contents() );
    table->setMargin( 8 );
    el = table->elementAt(0,0);
    new WLabel( "Include:", el );
    
    el = table->elementAt(0,1);
    WCheckBox *cb = new WCheckBox( "act.", el );
    cb->setMargin( 5, Wt::Left );
    cb->setChecked( true );
    cb->addStyleClass( "CbNoLineBreak" );
    cb->checked().connect( boost::bind( &DecayCsvResource::setGiveActivities, m_csvResouce, true ) );
    cb->unChecked().connect( boost::bind( &DecayCsvResource::setGiveActivities, m_csvResouce, false ) );
    
    el = table->elementAt(0,2);
    cb = new WCheckBox( "xrays", el );
    cb->setMargin( 5, Wt::Left );
    cb->addStyleClass( "CbNoLineBreak" );
    cb->checked().connect( boost::bind( &DecayCsvResource::setGiveXrays, m_csvResouce, true ) );
    cb->unChecked().connect( boost::bind( &DecayCsvResource::setGiveXrays, m_csvResouce, false ) );
    
    el = table->elementAt(0,3);
    cb = new WCheckBox( "gammas", el );
    cb->setMargin( 5, Wt::Left );
    cb->addStyleClass( "CbNoLineBreak" );
    cb->checked().connect( boost::bind( &DecayCsvResource::setGiveGammas, m_csvResouce, true ) );
    cb->unChecked().connect( boost::bind( &DecayCsvResource::setGiveGammas, m_csvResouce, false ) );
    
    el = table->elementAt(0,4);
    cb = new WCheckBox( "alphas", el );
    cb->setMargin( 5, Wt::Left );
    cb->addStyleClass( "CbNoLineBreak" );
    cb->checked().connect( boost::bind( &DecayCsvResource::setGiveAlphas, m_csvResouce, true ) );
    cb->unChecked().connect( boost::bind( &DecayCsvResource::setGiveAlphas, m_csvResouce, false ) );
    
    el = table->elementAt(0,5);
    cb = new WCheckBox( "betas", el );
    cb->setMargin( 5, Wt::Left );
    cb->addStyleClass( "CbNoLineBreak" );
    cb->checked().connect( boost::bind( &DecayCsvResource::setGiveBetas, m_csvResouce, true ) );
    cb->unChecked().connect( boost::bind( &DecayCsvResource::setGiveBetas, m_csvResouce, false ) );
    
    //table = new WTable( footer() );
    //table->setWidth( WLength(100,WLength::Percentage) );
    //footer()->setOverflow( WContainerWidget::OverflowHidden );
    
    //el = table->elementAt( 0, 0 );
    //WPushButton *button = new WPushButton( "Cancel", el );
    //button->setAttributeValue( "style", "margin:auto;display:block;" + button->attributeValue("style") );
    //button->clicked().connect( m_parent, &DecayActivityDiv::deleteCsvDownloadGui );
    
    WAnchor *csvancor = new WAnchor( footer() );
    csvancor->setLink( WLink(m_csvResouce) );
    csvancor->setTarget( TargetNewWindow );
    csvancor->setText( "Export CSV" );
    csvancor->setPadding( 8, Wt::Top );
    csvancor->setPadding( 2, Wt::Bottom );
    csvancor->setStyleClass( "ExportCsvLink" );
    footer()->setStyleClass( "modal-footer" );
    footer()->setContentAlignment( Wt::AlignCenter );
    
    //button = new WPushButton( "Export CSV", el );
    //button->setAttributeValue( "style", "margin:auto;display:block;" + button->attributeValue("style") );
    //button->setLink( WLink(m_csvResouce) );
    //button->setLinkTarget( Wt::TargetNewWindow );
        
#if( ANDROID )
    // Using hacked saving to temporary file in Android, instead of via network download of file.
    csvancor->clicked().connect( std::bind([this](){
      android_download_workaround(m_csvResouce, "nuc_decay.csv");
    }) );
#endif //ANDROID

    
    show();
    centerWindow();
  }//CsvDownloadGui
  
  void updateTimeSpan()
  {
    double timespan = 0.0;
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    
    try
    {
      const string input = m_ageEdit->valueText().toUTF8();
      
      const SandiaDecay::Nuclide *nuc = NULL;
      if( m_parent->m_nuclides.size() > 0 )
        nuc = db->nuclide( m_parent->m_nuclides[0].z, m_parent->m_nuclides[0].a, m_parent->m_nuclides[0].iso );
      
      if( nuc )
        timespan = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( input, nuc->halfLife );
      else
        timespan = PhysicalUnitsLocalized::stringToTimeDuration( input );
    }catch(...)
    {
      timespan = m_parent->timeToDisplayTill();
      m_ageEdit->setText( PhysicalUnitsLocalized::printToBestTimeUnits(timespan) );
    }
    
    m_csvResouce->setTimeSpan( timespan );
  }//void updateTimeSpan()
  
};//class CsvDownloadGui

DecayActivityDiv::DecayActivityDiv( InterSpec *viewer, Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_viewer( viewer ),
  m_nuclides( 0 ),
  m_displayTimeLength( NULL ),
  m_displayActivityUnitsCombo( NULL ),
  m_displayActivityUnitsLabel( NULL ),
  m_logYScale( NULL ),
  m_photopeakLogYScale( NULL ),
  m_showGridLines( NULL ),
  m_yAxisType( NULL ),
  m_parentNuclidesDiv( NULL ),
  m_nuclidesAddedDiv( NULL ),
  m_createNewNuclideButton( NULL ),
  m_clearNuclidesButton( NULL ),
  m_nuclideSelectDialog( NULL ),
  m_nuclideSelect( NULL ),
  m_chartTabWidget( NULL ),
  m_decayChart( NULL ),
  m_decayModel( NULL ),
  m_moreInfoDialog( NULL ),
  m_decayLegend( NULL ),
  m_calc( NULL ),
  m_csvDownloadDialog( NULL ),
  m_currentTimeUnits( -1.0 ),
  m_currentTimeRange( -1.0 ),
  m_currentNumXPoints( 250 ),
  m_currentMixture( new SandiaDecay::NuclideMixture() )
{
  addStyleClass( "DecayActivityDiv" );
  
  wApp->useStyleSheet( "InterSpec_resources/NuclideDecayInfo.css");
  
  if( m_viewer )
    m_viewer->useMessageResourceBundle( "DecayActivity" );
  
  init();
}//DecayActivityDiv constructor


void DecayActivityDiv::globalKeyPressed( const Wt::WKeyEvent &e )
{
  if( (e.key() == Wt::Key_A)
      && ((e.modifiers() & Wt::ShiftModifier)
          || (e.modifiers() & Wt::AltModifier)
          || (e.modifiers() & Wt::ControlModifier)
          || (e.modifiers() & Wt::MetaModifier))
     )
  {
    m_nuclideSelectDialog->show();
    m_nuclideSelectDialog->centerWindow();
    m_nuclideSelect->setNuclideSearchToFocus();
  }else if( (e.key() == Wt::Key_Left) || (e.key() == Wt::Key_Right) )
  {
    const int num = m_chartTabWidget->count();
    const int current = m_chartTabWidget->currentIndex();
    
    int next = current + ((e.key()==Wt::Key_Left) ? -1 : 1);
    if( next < 0 )
      next = num - 1;
    next = (next % num);
    m_chartTabWidget->setCurrentIndex( next );
  }//if / else determine key pressed
}//void DecayActivityDiv::globalKeyPressed( const Wt::WKeyEvent &e )


void DecayActivityDiv::init()
{
  //First delete any widgets that may be owned my *this
  WContainerWidget::clear();

//  wApp->globalKeyPressed().connect( this, &DecayActivityDiv::globalKeyPressed );
  wApp->globalKeyWentDown().connect( this, &DecayActivityDiv::globalKeyPressed );
  
  const bool isPhone = m_viewer ? m_viewer->isPhone() : false;
  
  //Initialize all member widgets; will assign parentage later on
//  m_nuclides( 0 ),
  m_displayTimeLength      = new WLineEdit( "" );
  m_displayTimeLength->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_displayTimeLength->setAttributeValue( "autocorrect", "off" );
  m_displayTimeLength->setAttributeValue( "spellcheck", "off" );
#endif
  
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), m_displayTimeLength );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_displayTimeLength->setValidator(validator);

  m_displayActivityUnitsCombo      = new WComboBox();
  m_displayActivityUnitsLabel      = new WLabel( WString::tr(isPhone ? "dad-units-label-phone" : "dad-units-label") );
  m_logYScale                      = new WCheckBox( WString::tr("dad-logy-cb") );
  m_showGridLines                  = new WCheckBox( WString::tr("dad-grid-cb") );
  m_yAxisType                      = new WComboBox();
  m_parentNuclidesDiv              = new WContainerWidget();
  m_nuclidesAddedDiv               = new WContainerWidget();
  m_createNewNuclideButton         = new WPushButton( WString::tr("dad-add-nucs") );
  m_clearNuclidesButton            = new WPushButton( WString::tr(isPhone ? "Clear" : "dad-remove-all")  );
  m_nuclideSelectDialog            = new AuxWindow( WString::tr("dad-sel-nuc-window-title"),
                                      (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen) | AuxWindowProperties::DisableCollapse) );
  
  m_nuclideSelect                  = new DecaySelectNuclide( isPhone, nullptr, m_nuclideSelectDialog );
  m_decayLegend                    = new WContainerWidget();

  m_chartTabWidget                 = new WTabWidget();
  m_decayChart                     = new DecayActivityChart();
  m_decayModel                     = new DecayActivityModel( this );

  //Set elements names so we can style with css
  m_parentNuclidesDiv->addStyleClass( "m_parentNuclidesDiv" );
  m_nuclidesAddedDiv->addStyleClass( "m_nuclidesAddedDiv" );

  m_createNewNuclideButton->setIcon( "InterSpec_resources/images/plus_min_white.svg" );
  m_createNewNuclideButton->addStyleClass( "m_createNewNuclideButton" );
  
  m_clearNuclidesButton->setIcon( "InterSpec_resources/images/remove_all.png" );
  m_clearNuclidesButton->addStyleClass( "m_clearNuclidesButton" );
    
  m_displayActivityUnitsCombo->addStyleClass( isPhone ? "DisplayActivityUnitsComboPhone" : "DisplayActivityUnitsCombo" );
  m_displayActivityUnitsLabel->addStyleClass( "DisplayActivityUnitsLabel" );
  m_logYScale->addStyleClass( "DecayChartLogYScale" );
  m_yAxisType->addStyleClass( "DecayChartYaxisType" );
  m_showGridLines->addStyleClass( "DecayChartShowGridLines" );
  m_decayChart->addStyleClass( "DecayChart" );
  m_chartTabWidget->addStyleClass( "DecayChartTabWidget" );
  m_chartTabWidget->contentsStack()->addStyleClass( "ChartTabWidgetStack" );
  
  m_decayLegend->addStyleClass( "DecayLegend" );

  m_displayTimeLength->addStyleClass( "m_displayTimeLength" );

  WContainerWidget *displayOptionsDiv = initDisplayOptionWidgets();

  //Time to get to work!
  initCharts();

  //setup m_parentNuclidesDiv
  m_parentNuclidesDiv->setInline( false );
  m_parentNuclidesDiv->addWidget( m_createNewNuclideButton );
  m_parentNuclidesDiv->addWidget( m_nuclidesAddedDiv );
  m_parentNuclidesDiv->addWidget( m_clearNuclidesButton );

  m_nuclidesAddedDiv->setInline( true );
  m_nuclideSelectDialog->contents()->addWidget( m_nuclideSelect );
  if( m_viewer && (m_viewer->renderedHeight() > 100) )
    m_nuclideSelectDialog->setMaximumSize( WLength::Auto, 0.9*m_viewer->renderedHeight() );
  
  m_nuclideSelectDialog->rejectWhenEscapePressed();
  m_createNewNuclideButton->clicked().connect( boost::bind( &WDialog::show,
                                                      m_nuclideSelectDialog ) );
  m_createNewNuclideButton->clicked().connect( m_nuclideSelect,
                                   &DecaySelectNuclide::setNuclideSearchToFocus );
  
  m_nuclideSelect->done().connect( m_nuclideSelectDialog, &AuxWindow::hide );
  m_nuclideSelectDialog->finished().connect( this, &DecayActivityDiv::nuclideSelectDialogDone );
  m_nuclideSelectDialog->centerWindow();
  
  m_nuclideSelectDialog->hide();
  
  m_nuclideSelect->selected().connect( boost::bind( &DecayActivityDiv::addTheNuclide, this,
                                                   boost::placeholders::_1 ) );
  m_clearNuclidesButton->clicked().connect( this, &DecayActivityDiv::clearAllNuclides );
  m_decayChart->clicked().connect( this, &DecayActivityDiv::decayChartClicked );

  WGridLayout *decayLayout = new WGridLayout();
  WContainerWidget *decayDiv = new WContainerWidget();
  decayDiv->addStyleClass( "DecayActivityChartDiv" ); /* TODO: Hacked to avoid dealing with color theme for the moment */
  decayDiv->addStyleClass( "decayDiv" );
  decayDiv->setLayout( decayLayout );
  
  ChartAndLegendHolder *chartHolder
                     = new ChartAndLegendHolder( m_decayChart, m_decayLegend, isPhone );
  
  decayLayout->addWidget( chartHolder, 0, 0, 1, 1 );
  decayLayout->addWidget( displayOptionsDiv, 1, 0, 1, 1 );
  decayLayout->setRowStretch( 0, 5 );

  //////////////////
  // Decay Chain
  m_decayChainChart = new DecayChainChart();
  WGridLayout *decayChainLayout = new WGridLayout();
  WContainerWidget *decayChainDiv = new WContainerWidget();
  decayChainDiv->setMinimumSize( 300, 350 );
  decayChainDiv->addStyleClass( "decayChainDiv" ); //ToDo: Currently decayChainDiv is hacked to avoid dealing with color theme ish
  decayChainDiv->setOverflow( WContainerWidget::OverflowAuto );
  decayChainDiv->setLayout( decayChainLayout );
  
  //decayChainLayout->addWidget(elementTemplate, 0, 0, 1, 1, AlignCenter);
  decayChainLayout->addWidget( m_decayChainChart, 0, 0, 1, 1 );
  decayChainLayout->setRowStretch( 0, 1 );

  //////////////////
  // tab widget
  const WTabWidget::LoadPolicy loadPolicy = WTabWidget::LazyLoading; //WTabWidget::PreLoading;
  m_chartTabWidget->addTab( decayDiv, WString::tr("dad-mi-act-chart"), loadPolicy );

  m_chartTabWidget->addTab( decayChainDiv, WString::tr("dad-mi-decay-chain"), loadPolicy );
  m_chartTabWidget->currentChanged().connect( m_decayChainChart,
                                    &DecayChainChart::deleteMoreInfoDialog );
  m_decayChainChart->nuclideChanged().connect( this,
                            &DecayActivityDiv::manageActiveDecayChainNucStyling );


  m_chartTabWidget->currentChanged().connect( this,
                                        &DecayActivityDiv::deleteMoreInfoDialog );
  m_chartTabWidget->currentChanged().connect( this,
                            &DecayActivityDiv::manageActiveDecayChainNucStyling );

  m_calc = new DateLengthCalculator( this, NULL );
  m_chartTabWidget->addTab( m_calc, WString::tr("dad-mi-calc"), loadPolicy );
  
  
  WGridLayout *layout = new WGridLayout();
  layout->setContentsMargins( 3, 3, 3, 3 );
  layout->addWidget(  m_parentNuclidesDiv, 0, 0, 1, 1 );
  layout->addWidget(  m_chartTabWidget,    1, 0, 1, 1 );
  layout->setRowStretch( 1, 5 );
  WContainerWidget::setLayout( layout );
  
  if( m_viewer )
    m_viewer->colorThemeChanged().connect( this, &DecayActivityDiv::colorThemeChanged );
}//void DecayActivityDiv::init()



void DecayActivityDiv::handleAppUrl( std::string path, std::string query_str )
{
  auto isNucKey = []( const std::string &key ) -> bool {
    return (key == "iso") || (key == "nuc") || (key == "isotope") || (key == "nuclide");
  };
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  
  clearAllNuclides();

  //setDecayChartTimeRange( double dt );
  if( !query_str.empty() && (query_str[0] == '?') )
    query_str = query_str.substr(1);
  
  vector<string> query_args;
  SpecUtils::split( query_args, query_str, "&" );
  
  size_t num_nucs = 0, last_nuc_index = 0;
  vector<pair<string,string>> field_values;
  for( size_t i = 0; i < query_args.size(); ++i )
  {
    string comp = query_args[i];
    SpecUtils::trim( comp );
    if( comp.empty() )
      continue;
    
    auto pos = comp.find('=');
    if( pos == string::npos )
      pos = comp.find(':');
    
    string key = (pos == string::npos) ? comp : comp.substr(0,pos);
    string value = (pos == string::npos) ? string("") : comp.substr(pos+1);
    
    SpecUtils::trim(key);
    SpecUtils::trim(value);
    
    //key = Wt::Utils::urlDecode( key );
    //value = Wt::Utils::urlDecode( value );
    
    SpecUtils::to_lower_ascii(key);
    
    field_values.push_back( {key,value} );
    
    if( isNucKey(key) )
    {
      num_nucs += 1;
      last_nuc_index = i;
    }
  }//for( string comp : query_args )
  
  // If only a single nuclide, the ordering or URL arguments wont matter, but if there are multiple
  //  nuclides, then it will matter (e.g., age, activity, etc, must come AFTER the nuclide)
  if( (num_nucs == 1) && (last_nuc_index != 0) )
    std::swap( field_values[0], field_values[last_nuc_index] );
  
  string dispTimeLen;
  vector<string> error_messages;
  vector<const SandiaDecay::Nuclide *> nuclides;
  
  // Loop over and check for activity units to use, as well as time span to display of
  for( size_t i = 0; i < field_values.size(); ++i )
  {
    const string &key = field_values[i].first;
    string &value = field_values[i].second;
    
    if( key == "actunits" )
    {
      SpecUtils::to_lower_ascii( value );
      
      if( (value == "becquerel") || (value == "bq") )
        useBq = true;
      else if( (value == "curie") || (value == "ci") )
        useBq = false;
      else
        error_messages.push_back( "Unrecognized activity units '" + value + "', must be one of"
                                  " 'becquerel', 'bq', 'curie', 'ci'" );
    }else if( (key == "time") || (key == "timespan") )
    {
      dispTimeLen = value;
    }else if( (key == "grid") || (key == "gridlines") )
    {
      const bool show = ((value == "1") || (value == "true") || (value == "yes"));
      m_showGridLines->setChecked( show );
      m_decayChart->showGridLines( show );
    }else if( key == "charttype" )
    {
      if( (value == "act") || (value == "activity") )
        m_yAxisType->setCurrentIndex( YAxisType::ActivityAxis );
      else if( value == "gamma" )
        m_yAxisType->setCurrentIndex( YAxisType::GammasAxis );
      else if( value == "beta" )
        m_yAxisType->setCurrentIndex( YAxisType::BetasAxis );
      else if( value == "alpha" )
        m_yAxisType->setCurrentIndex( YAxisType::AlphasAxis );
      else
        error_messages.push_back( "Unrecognized chart type '" + value
                                  + "', must be Activity, Gamma, Beta, or Alpha." );
    }else if( key == "logy" )
    {
      const bool logy = ((value == "1") || (value == "true") || (value == "yes"));
      m_logYScale->setChecked( logy );
      updateYScale();
    }else if( !isNucKey(key) && (key != "age") && (key != "initialage")
             && (key != "act") && (key != "activity") )
    {
      error_messages.push_back( "Unrecognized query key '" + key + "'." );
    }
  }//for( size_t i = 0; i < field_values.size(); ++i )
  
  
  // Loop over and add all nuclides
  for( size_t i = 0; i < field_values.size(); ++i )
  {
    const string &key = field_values[i].first;
    string &value = field_values[i].second;
    
    if( isNucKey(key) )
    {
      const SandiaDecay::Nuclide *nuc = db->nuclide(value);
      if( !nuc )
      {
        error_messages.push_back( "'" + value + "' not a valid nuclide." );
        continue;
      }
      
      if( nuc->isStable() )
      {
        error_messages.push_back( "'" + value + "' is a stable nuclide." );
        continue;
      }
            
      string activity_str, age_str;
      for( size_t j = i + 1; j < field_values.size(); ++j )
      {
        if( isNucKey(field_values[j].first) )
          break;
        
        if( (field_values[j].first == "age") || (field_values[j].first == "initialage") )
        {
          age_str = field_values[j].second;
        }else if( (field_values[j].first == "act")|| (field_values[j].first == "activity") )
        {
          activity_str = field_values[j].second;
        }
      }//for( size_t j = i + 1; j < field_values.size(); ++j )
      
      
      double act = 1.0*PhysicalUnits::mCi, age = 0.0;
      
      try
      {
        const double actFromStr = PhysicalUnits::stringToActivity( activity_str );
        if( actFromStr > 0.0 )
          act = actFromStr;
      }catch(...)
      {
      }
      
      try
      {
        const double ageFromStr = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( age_str, nuc->halfLife );
        if( ageFromStr >= 0.0 )
          age = ageFromStr;
      }catch(...)
      {
      }//try / catch for age
      
      nuclides.push_back( nuc );
      addNuclide( nuc->atomicNumber, nuc->massNumber, nuc->isomerNumber, act, !useBq, age, activity_str);
    }
  }//for( size_t i = 0; i < field_values.size(); ++i )
  
  
  if( SpecUtils::iequals_ascii(path, "chart") )
  {
    m_chartTabWidget->setCurrentIndex( 0 );
  }else if( SpecUtils::iequals_ascii(path, "chain") )
  {
    m_chartTabWidget->setCurrentIndex( 1 );
  }else if( SpecUtils::iequals_ascii(path, "calc") )
  {
    m_chartTabWidget->setCurrentIndex( 2 );
  }else if( !path.empty() )
  {
    error_messages.push_back( "Unrecognized component '" + path + "' to show on decay window." );
  }
  
  m_displayTimeLength->setText( dispTimeLen );
  
  refreshDecayDisplay( true );
  
  // Display up to the first 5 (arbitrary) error messages
  for( size_t i = 0; (i < 5) && (i < error_messages.size()); ++i )
    passMessage( error_messages[i], WarningWidget::WarningMsgHigh );
}//void handleAppUrl( std::string path, std::string query_str )


std::string DecayActivityDiv::encodeStateToUrl()
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  string path, query_str;
  
  switch( m_chartTabWidget->currentIndex() )
  {
    case 0: path += "chart"; break;
    case 1: path += "chain"; break;
    case 2: path += "calc"; break;
    default: assert( 0 ); break;
  }//
  
  query_str = "grid=" + string(m_showGridLines->isChecked() ? "1" : "0");
  query_str += "&logy=" + string(m_logYScale->isChecked() ? "1" : "0");
  
  const bool usebq = SpecUtils::icontains( m_displayActivityUnitsCombo->currentText().toUTF8(), "bq" );
  query_str += "&actunits=" + string(usebq ? "bq" : "ci");
  
  
  switch( static_cast<YAxisType>( m_yAxisType->currentIndex() ) )
  {
    case ActivityAxis: query_str += "&charttype=activity"; break;
    case GammasAxis:   query_str += "&charttype=gamma";    break;
    case BetasAxis:    query_str += "&charttype=beta";     break;
    case AlphasAxis:   query_str += "&charttype=alpha";    break;
    case NumYAxisType: break;
  }//switch( y-axis type )
  
  
  const double timespan = timeToDisplayTill();
  if( timespan > FLT_EPSILON )
    query_str += "&timespan=" + PhysicalUnits::printToBestTimeUnits(timespan,6);
  
  
  for( const Nuclide &nucinfo : m_nuclides )
  {
    const SandiaDecay::Nuclide *nuc = db->nuclide( nucinfo.z, nucinfo.a, nucinfo.iso );
    if( !nuc )
      continue;  //shouldnt ever happen
    query_str += "&nuc="  + nuc->symbol;
    query_str += "&act="  + PhysicalUnits::printToBestActivityUnits(nucinfo.activity,4); // TODO: sanitize and use nucinfo.activityStr
    if( nucinfo.age > FLT_EPSILON )
      query_str += "&initialage="  + PhysicalUnits::printToBestTimeUnits(nucinfo.age,6);
  }//for( const Nuclide &nucinfo : m_nuclides )
  
  return path + "?" + query_str;
}//std::string encodeStateToUrl()


void DecayActivityDiv::setGridLineStatus()
{
  m_decayChart->showGridLines( m_showGridLines->isChecked() );
}//void DecayActivityDiv::setGridLineStatus()


Wt::WContainerWidget *DecayActivityDiv::initDisplayOptionWidgets()
{
  const bool isPhone = m_viewer ? m_viewer->isPhone() : false;
  
  m_showGridLines->setUnChecked();
  m_logYScale->setUnChecked();
  updateYScale();
  m_logYScale->checked().connect( this, &DecayActivityDiv::updateYScale );
  m_logYScale->unChecked().connect( this, &DecayActivityDiv::updateYScale );
  m_showGridLines->checked().connect(this, &DecayActivityDiv::setGridLineStatus );
  m_showGridLines->unChecked().connect(this, &DecayActivityDiv::setGridLineStatus );

  for( const auto &nuclide : PhysicalUnits::sm_activityUnitNameValues )
  {
    m_displayActivityUnitsCombo->addItem( nuclide.first );
  }//for( each defined unit of activity )

  m_displayActivityUnitsCombo->addItem( WString::tr("dad-arbitrary") );

  const bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  m_displayActivityUnitsCombo->setCurrentIndex( (useBq ? 2 : 7) ); //MBq : mCi


  WContainerWidget *displayOptionsDiv = new WContainerWidget();
  displayOptionsDiv->addStyleClass( "displayOptionsDiv" );
  WContainerWidget *displOptUpper = new WContainerWidget( displayOptionsDiv );
  
  WLabel *endLabel = new WLabel( WString::tr("dad-time-span-label") );
  endLabel->setBuddy( m_displayTimeLength );
  displOptUpper->addWidget( endLabel );
  displOptUpper->addWidget( m_displayTimeLength );
  m_displayTimeLength->changed().connect( boost::bind( &DecayActivityDiv::refreshDecayDisplay, this, true ) );
  m_displayTimeLength->enterPressed().connect( boost::bind( &DecayActivityDiv::refreshDecayDisplay, this, true ) );
  
  
  const bool showToolTips = m_viewer ? UserPreferences::preferenceValue<bool>( "ShowTooltips", m_viewer ) : false;
  if( m_viewer )
    HelpSystem::attachToolTipOn( m_displayTimeLength, WString::tr("dad-tt-age"), showToolTips );
  
  WLabel *yaxisTypeLabel = new WLabel( WString::tr(isPhone ? "dad-yaxis-label-phone" : "dad-yaxis-label"), displOptUpper );
  yaxisTypeLabel->addStyleClass( "DecayChartYaxisTypeLabel" );
  
  displOptUpper->addWidget( m_yAxisType );
  for( YAxisType y = YAxisType(0); y < NumYAxisType; y = YAxisType(int(y)+1) )
  {
    switch( y )
    {
      case ActivityAxis: m_yAxisType->addItem( WString::tr("Activity") );   break;
      case GammasAxis:   m_yAxisType->addItem( WString::tr("dad-gamma-rate") ); break;
      case BetasAxis:    m_yAxisType->addItem( WString::tr("dad-beta-rate") );  break;
      case AlphasAxis:   m_yAxisType->addItem( WString::tr("dad-alpha-rate") ); break;
      case NumYAxisType:                                       break;
    }//switch( y )
  }//for( Y-Axis type )
  
  //If phone, skip creating displOptLower to make things more compact
  WContainerWidget *displOptLower = nullptr;
  if( m_viewer->isPhone() )
  {
    displOptLower = displOptUpper;
  }else
  {
    displOptLower = new WContainerWidget( displayOptionsDiv );
    displOptLower->setMargin( 5, Wt::Top );
  }
  
  
  displOptLower->addWidget( m_displayActivityUnitsLabel );
  displOptLower->addWidget( m_displayActivityUnitsCombo );
  
  displOptLower->addWidget( m_logYScale );
  displOptLower->addWidget( m_showGridLines );
  
  if( m_viewer->isPhone() )
  {
    m_logYScale->hide();
    m_showGridLines->hide();
  }
  
  m_yAxisType->setCurrentIndex( ActivityAxis );
  m_yAxisType->activated().connect( boost::bind( &DecayActivityDiv::refreshDecayDisplay, this, true ) );
  
  if( !m_viewer->isPhone() )
  {
    WPushButton *csvButton = new WPushButton( displOptLower );
    csvButton->setIcon( "InterSpec_resources/images/download_small.svg" );
    csvButton->setText( WString("{1}...").arg( WString::tr("CSV") ) );
    csvButton->setStyleClass( "LinkBtn DownloadBtn" );
    csvButton->clicked().connect( this, &DecayActivityDiv::createCsvDownloadGui );
  }

  m_displayActivityUnitsCombo->changed().connect( boost::bind( &DecayActivityDiv::refreshDecayDisplay, this, true ) );

  return displayOptionsDiv;
}//Wt::WContainerWidget *initDisplayOptionWidgets()



void DecayActivityDiv::createCsvDownloadGui()
{
  deleteCsvDownloadGui();
  m_csvDownloadDialog = new CsvDownloadGui( this );
}//void createCsvDownloadGui()


void DecayActivityDiv::deleteCsvDownloadGui()
{
  if( m_csvDownloadDialog )
  {
    m_csvDownloadDialog->finished().setBlocked(true);
    AuxWindow::deleteAuxWindow( m_csvDownloadDialog );
    m_csvDownloadDialog = nullptr;
  }
}//void deleteCsvDownloadGui()


void DecayActivityDiv::deleteCsvDownloadGuiTriggerUpdate()
{
  deleteCsvDownloadGui();
  auto app = wApp;
  if( app )
    app->triggerUpdate();
}//void deleteCsvDownloadGuiTriggerUpdate()


void DecayActivityDiv::initCharts()
{
  //Init basic features of the chart, first the decay chart
  m_decayChart->setMinimumSize( WLength(250), WLength(100) );
//  m_decayChart->resize( WLength(800), WLength(500) );

  m_decayChart->setModel( m_decayModel );
  m_decayChart->setXSeriesColumn( 0 );
//  m_decayChart->axis(Chart::XAxis).setScale( Chart::DateTimeScale );
  m_decayChart->setType( Chart::ScatterPlot );
  m_decayChart->axis(Chart::XAxis).setScale( Chart::LinearScale );

  //TODO 20110325: The either of next line doesnt seem to have any effect
  m_decayChart->axis(Chart::XAxis).setLabelFormat( "%.3g" );

  //  decayLayout->addWidget( m_decayChart, 0, 0, 1, 1 );
  if( m_viewer && m_viewer->isPhone() )
  {
    m_decayChart->setIsPhone();
    m_decayChart->setPlotAreaPadding( 75, Left );
    m_decayChart->setPlotAreaPadding( 25, Bottom );
    m_decayChart->setPlotAreaPadding( 20, Right );
    m_decayChart->setPlotAreaPadding( 0, Top );
  }else
  {
    m_decayChart->setPlotAreaPadding( 80, Left );
    m_decayChart->setPlotAreaPadding( 50, Bottom );
    m_decayChart->setPlotAreaPadding( 20, Right );
    m_decayChart->setPlotAreaPadding( 5, Top );
  }
  
  if( m_viewer && m_viewer->getColorTheme() )
    m_decayChart->setColorTheme( m_viewer->getColorTheme() );

  //do not display the legend to the right of the chat
  m_decayChart->setLegendEnabled( false);
  //m_decayChart->axis(Chart::XAxis).setTitleOffset( 0.0 );
  
//  m_decayChart->setMargin(0);
//  m_decayChart->initLayout();

  m_decayChart->setToolTip( WString::tr("dad-click-for-more-info") );
//  m_decayChart->mouseMoved().connect( this, &DecayActivityDiv::updateMouseOver );
}//void initCharts()



double DecayActivityDiv::timeToDisplayTill()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  std::string txt = m_displayTimeLength->text().narrow();
  
  const SandiaDecay::Nuclide *nuc = NULL;
  if( m_nuclides.size() > 0 )
    nuc = db->nuclide( m_nuclides[0].z, m_nuclides[0].a, m_nuclides[0].iso );
  
  try
  {
    if( nuc )
      return PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( txt, nuc->halfLife, SandiaDecay::second);
    else
      return PhysicalUnitsLocalized::stringToTimeDuration( txt );
  }catch( std::exception & )
  {
    const double fracT0 = 0.1;
    const double finalTime = findTimeForActivityFrac( m_currentMixture, fracT0 );
    
    txt = PhysicalUnitsLocalized::printToBestTimeUnits( finalTime );
    m_displayTimeLength->setText( txt );
    return finalTime;
  }
}//double DecayActivityDiv::timeToDisplayTill()



void DecayActivityDiv::nuclideSelectDialogDone()
{
//  m_nuclideSelectDialog->accept();
  for( Nuclide &nuc : m_nuclides )
  {
    const string styleclass = nuc.display->styleClass().toUTF8();
    
    if( styleclass.find("EditingNuclide") != string::npos )
      nuc.display->removeStyleClass( "EditingNuclide" );
  }//for( Nuclide &nuc, m_nuclides )
  
  m_nuclideSelect->setAddButtonToAdd();
}//nuclideSelectDialogDone()


void DecayActivityDiv::sourceNuclideDoubleClicked( Wt::WContainerWidget *w )
{
  for( Nuclide &nuc : m_nuclides )
  {
    if( nuc.display != w )
      continue;
    
    nuc.display->addStyleClass( "EditingNuclide" );
    m_nuclideSelectDialog->show();
    m_nuclideSelect->setNuclideSearchToFocus();
    m_nuclideSelect->setCurrentInfo( nuc.a, nuc.z, nuc.iso,
                                     nuc.age, nuc.activity, nuc.useCurie, nuc.activityStr );
    m_nuclideSelect->setAddButtonToAccept();
    m_nuclideSelectDialog->centerWindow();
    m_nuclideSelectDialog->resizeToFitOnScreen();
    break;
  }//for( Nuclide &nuc : m_nuclides )
}//void sourceNuclideDoubleClicked( Wt::WContainerWidget *w );


void DecayActivityDiv::addTheNuclide( const NuclideSelectedInfo &n )
{
  addNuclide( n.z, n.a, n.metasable, n.activity, 
              n.useCurie, n.initialAge, n.activityStr);
}//void addTheNuclide( const NuclideSelectedInfo &nuc )


void DecayActivityDiv::Nuclide::updateTxt()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Element *element = db->element( z );
  const string elementName = (element ? element->symbol : string(""));
  
  stringstream label;
  label << "<sup><font size=\"1.5\">" << a;
  
  if( iso )
    label << "m";
  
  if( iso > 1 )
    label << iso;
  
  label << "</font></sup>" << elementName;
  label << " " << PhysicalUnits::printToBestActivityUnits(activity, 2, useCurie );
  
  if( age > 0.0 )
  {
    label << " " << PhysicalUnitsLocalized::printToBestTimeUnits(age);
  }//if( age > 0.0 )
  
  txt->setText( label.str() );
}//void DecayActivityDiv::Nuclide::updateTxt()


void DecayActivityDiv::addNuclide( const int z, const int a, const int iso,
                                 const double activity, const bool useCurie,
                                 const double age, const std::string &activityStr )
{
  //See if we are actually editing a Nuclide
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  WContainerWidget *editedNuclide = (WContainerWidget *)0;
  for( Nuclide &nuc : m_nuclides )
  {
    const string styleclass = nuc.display->styleClass().toUTF8();
    if( styleclass.find("EditingNuclide") != string::npos )
    {
      editedNuclide = nuc.display;
      break;
    }//if( styleclass has "EditingNuclide"  in it )
  }//for( Nuclide &nuc : m_nuclides )
  
  
  Nuclide nuclide;
  nuclide.a         = a;
  nuclide.z         = z;
  nuclide.activity  = activity;
  nuclide.useCurie  = useCurie;
  nuclide.activityStr = activityStr;
  nuclide.iso       = iso;
  nuclide.age       = age;
  nuclide.display   = new WContainerWidget();
  nuclide.display->setInline( true );
  nuclide.display->addStyleClass( "Nuclide" );
  nuclide.txt = new WText( "", XHTMLText, nuclide.display );

  if( !editedNuclide )
  {
    m_nuclidesAddedDiv->addWidget( nuclide.display );
  }else
  {
    m_nuclidesAddedDiv->insertBefore( nuclide.display, editedNuclide );
    removeNuclide( editedNuclide );
  }

  const SandiaDecay::Nuclide *nuc = db->nuclide( z, a );
  if( nuc )
  {
    nuclide.display->clicked().connect(
          boost::bind( &DecayChainChart::setNuclide,
                       m_decayChainChart, nuc, useCurie,
                       DecayChainChart::DecayChainType::DecayFrom) );
  }//if( nuc )
  
  nuclide.display->doubleClicked().connect(
                    boost::bind( &DecayActivityDiv::sourceNuclideDoubleClicked,
                                 this, nuclide.display) );
  nuclide.display->setToolTip( WString::tr("dad-tt-double-click") );
  
  nuclide.updateTxt();

//  WText *closeIcon = new WText( "" );
//  WText *closeIcon = new WText();
//  WApplication::instance()->theme()->apply( this, closeIcon,
//                                           MenuItemCloseRole );
  WPushButton *closeIcon = new WPushButton();
  closeIcon->addStyleClass( "mycloseicon" );
  
  closeIcon->clicked().connect( boost::bind( &DecayActivityDiv::removeNuclide,
                                             this, nuclide.display ) );
  nuclide.display->addWidget( closeIcon );

  m_nuclides.push_back( nuclide );

  setTimeLimitToDisplay();

  const auto nucptr = db->nuclide(z, a, iso);
  m_decayChainChart->setNuclide( nucptr, useCurie, DecayChainChart::DecayChainType::DecayFrom );

  refreshDecayDisplay( true );
}//void addNuclide(..)


void DecayActivityDiv::removeNuclide( Wt::WContainerWidget *frame )
{
  int to_be_removed = -1;
  for( size_t index =0; index < m_nuclides.size(); ++index )
  {
    if( m_nuclides[index].display == frame )
    {
      to_be_removed = static_cast<int>( index );
    }//if( we found it )
  }//for( loop over indeces )

  if( to_be_removed < 0 ) //should never happen
  {
    cerr << "removeNuclide(...): Could not find the desired nuclide :(" << endl;
    return;
  }//if( to_be_removed < 0 )

  if( to_be_removed == (m_nuclides.size() - 1) )
  {
    if (m_nuclides.size() == 1)
      m_decayChainChart->setNuclide( nullptr, true, DecayChainChart::DecayChainType::DecayThrough );
    else
    {
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      const DecayActivityDiv::Nuclide nuc = m_nuclides[m_nuclides.size() - 2];
      const SandiaDecay::Nuclide * const nucptr = db->nuclide(nuc.z, nuc.a, nuc.iso);
      m_decayChainChart->setNuclide( nucptr, nuc.useCurie, DecayChainChart::DecayChainType::DecayThrough );
    }
  }

  delete m_nuclides[to_be_removed].display;
  m_nuclides.erase( m_nuclides.begin() + to_be_removed );

  setTimeLimitToDisplay();

  refreshDecayDisplay( true );
}//void DecayActivityDiv::removeNuclide( Wt::WContainerWidget *frame )


void DecayActivityDiv::clearAllNuclides()
{
  for( const Nuclide &nuclide : m_nuclides )
    delete nuclide.display;
  m_nuclides.clear();

  setTimeLimitToDisplay();
  
  m_decayChainChart->setNuclide( nullptr, true, DecayChainChart::DecayChainType::DecayFrom );

  deleteMoreInfoDialog();
  
  if( m_decayChainChart )
    m_decayChainChart->deleteMoreInfoDialog();

  refreshDecayDisplay( true );
}//void clearAllNuclides();





void DecayActivityDiv::updateInitialMixture() const
{
  using SandiaDecay::NuclideActivityPair;
  using SandiaDecay::NuclideNumAtomsPair;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();

  m_currentMixture->clear();

  for( const Nuclide &nuc : m_nuclides )
  {
    const SandiaDecay::Nuclide *dbnuclide = db->nuclide( nuc.z, nuc.a,
                                                                  nuc.iso );
    const double activity = nuc.activity;
    if( nuc.age > DBL_EPSILON )
      m_currentMixture->addAgedNuclideByActivity( dbnuclide, activity, nuc.age );
    else
      m_currentMixture->addNuclideByActivity( dbnuclide, activity );
  }//for( each m_nuclides )
}//void updateInitialMixture()


double DecayActivityDiv::findTimeForActivityFrac(
                                    const SandiaDecay::NuclideMixture *mixture,
                                    const double fracOfT0,
                                    const double startTime,
                                    const double endTime )
{
  //TODO 20110325: This function should be broken up, it is a little long
  //

  const double epsilon = fracOfT0 * 0.01;

  if( (fracOfT0 <= 0.0) || (fracOfT0>1.0) || (startTime > endTime) )
  {
    stringstream msg;
    msg << "findTimeForActivityFrac(...): Invalid fracOfT0=" << fracOfT0
        << ", or time range (" << startTime << ", " << endTime << ")";
    throw std::runtime_error( msg.str() );
  }//if( invalid input )

  const std::vector<SandiaDecay::NuclideTimeEvolution> &
                           evolutions = mixture->decayedToNuclidesEvolutions();

  double maxHalfLife = 0.0;

  for( size_t i = 0; i < evolutions.size(); ++i )
  {
    const double halfLife = evolutions[i].nuclide->halfLife;
    if( !IsInf(halfLife) && !IsNan(halfLife) )
      maxHalfLife = max( maxHalfLife, halfLife );
  }//for( loop over nuclides )


  const double initialActivty = mixture->totalActivity(0.0);

  if( initialActivty == 0.0 )
    return 0.0;

  const bool ranged = ( (startTime >= 0.0) && (endTime > 0.0) );
  const double dt = ( ranged  ? (endTime-startTime)/10.0 : maxHalfLife/2.0 );
  const double startT = ( ranged ? startTime : 0.0 );

  int iteration = 0;
  double frac = 1.0, time = startT, lastTime = startT;

  while( frac > fracOfT0 && (iteration++ < 5000000) )
  {
    const double activity = mixture->totalActivity( time  );
    frac = activity / initialActivty;

    if( abs( frac - fracOfT0 ) < epsilon )
      return time;

    if( frac <= fracOfT0 )
    {
      if( time == startT )//shouldnt ever happen
      {
        stringstream msg;
        msg << "DecayActivityDiv::findTimeForActivityFrac(...): Error, the start"
            << " time (=" << startTime << ") already has an  activity of "
            << frac << " times original activity";
        cerr << msg.str() << endl;    // throw std::runtime_error( msg );
        return time; //not throwing; maybe the caller to this function got lucky
      }//if( time == startT )

      //need to recursively search between the last time range and this one
      return findTimeForActivityFrac( mixture, fracOfT0, lastTime, time );
    }//if( frac < fracOfT0 )

    lastTime = time;
    time += dt;
  }//while( frac > fracOfT0 )

  //I doubt we will ever make it here, but just incase put a message in the
  //  apache log file, and return an invalid time
  stringstream msg;
  msg << "DecayActivityDiv::findTimeForActivityFrac(...):"
      << " couldnt find the time to end showing graph; searched "
     << (startT+iteration*dt)/second << " seconds at a dt=" << dt/second;

  cerr << msg.str() << endl;

  return -1.0;
}//double findTimeForActivityFrac(...)




void DecayActivityDiv::updateYScale()
{
  const Chart::AxisScale scale = m_logYScale->isChecked() ? Chart::LogScale : Chart::LinearScale;
  m_decayChart->axis(Chart::YAxis).setScale( scale );
    
  // If we dont do a full refresh, for some reason sometimes the lines wont show up
  refreshDecayDisplay( true );
}//void updateYScale()


void DecayActivityDiv::displayMoreInfoPopup( const double time )
{
  WContainerWidget *summary = isotopesSummary( time );

  summary->setMaximumSize( WLength::Auto, WLength( 0.8*m_viewer->renderedHeight() ,WLength::Pixel) );
  
  WString title = WString::tr("dad-summary-at-time").arg( PhysicalUnitsLocalized::printToBestTimeUnits(time) );

  if( m_moreInfoDialog )
  {
    m_moreInfoDialog->contents()->clear();
    m_moreInfoDialog->setWindowTitle( title );
  }else
  {
    m_moreInfoDialog = new AuxWindow( title );
    m_moreInfoDialog->setClosable( true );
    m_moreInfoDialog->disableCollapse();
    m_moreInfoDialog->rejectWhenEscapePressed();
    m_moreInfoDialog->finished().connect(
                  boost::bind( &DecayActivityDiv::deleteMoreInfoDialog, this ) );
      
    WPushButton *ok = m_moreInfoDialog->addCloseButtonToFooter();
    ok->clicked().connect(boost::bind( &DecayActivityDiv::deleteMoreInfoDialog, this ) );
  }//if( m_moreInfoDialog ) / else

  m_moreInfoDialog->contents()->addWidget( summary );
  
  //WGridLayout *layout = new WGridLayout();
  //container->setLayout( layout );
  //layout->addWidget( summary, 0, 0 );
  
  //layout->setRowStretch( 0, 10 );

  //I dont really care about centering the dialog, but without doing the below
  //  the scrolling of the contents can be buggy
//  const string jsthis = "$('#" + dialog->id() + "')";
//  stringstream moveJs;
//  moveJs << "var el = " << dialog->jsRef() << ";"
//         << "var ws = " << wApp->javaScriptClass() << ".WT.windowSize();"
//         << "el.style.top = Math.round( ( ws.y - " << jsthis << ".height() ) / 2 )+'px';"
//         << "el.style.left = Math.round( ( ws.x - " << jsthis << ".width() ) / 2 )+'px';";
//  doJavaScript( moveJs.str() );
  
  m_moreInfoDialog->setHidden( false );
  m_moreInfoDialog->resizeToFitOnScreen();
  m_moreInfoDialog->centerWindow();
}//void displayMoreInfoPopup( const double time )


void DecayActivityDiv::checkTimeRangeValid()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  std::string txt = m_displayTimeLength->text().narrow();
  
  const SandiaDecay::Nuclide *nuc = NULL;
  if( m_nuclides.size() > 0 )
    nuc = db->nuclide( m_nuclides[0].z, m_nuclides[0].a, m_nuclides[0].iso );
  
  double timelen = 0.0;
  
  try
  {
    if( nuc )
      timelen = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( txt, nuc->halfLife, SandiaDecay::second);
    else
      timelen = PhysicalUnitsLocalized::stringToTimeDuration( txt );
  }catch( std::exception & )
  {
    timelen = findTimeForActivityFrac( m_currentMixture, 0.1 );
    txt = PhysicalUnitsLocalized::printToBestTimeUnits( timelen );
    m_displayTimeLength->setText( txt );
  }

  
  double minhl = DBL_MAX;
  for( const Nuclide &nuc : m_nuclides )
  {
    const SandiaDecay::Nuclide *dbnuclide = db->nuclide( nuc.z, nuc.a,
                                                                 nuc.iso );
    minhl = std::min( minhl, dbnuclide->halfLife );
  }//for( const Nuclide &nuc : m_nuclides )
  
  minhl = std::max( minhl, 1.0E-6*SandiaDecay::second );
  minhl = std::min( minhl, 1000.0*SandiaDecay::second );
  if( timelen < 0.001*minhl )
  {
    timelen = 0.001*minhl;
    txt = PhysicalUnitsLocalized::printToBestTimeUnits( timelen );
    m_displayTimeLength->setText( txt );
  }//if( nunits*unitpair.second < 0.001*minhl )
}//void DecayActivityDiv::checkTimeRangeValid()


void DecayActivityDiv::setTimeLimitToDisplay()
{
  //NOTE 20110906: below TODO from 20110505 is possibly irreevant
  //TODO 20110505: we have the current population cached in m_currentPopulation
  //               however this may not be up to date if the user added/removed
  //               nuclides.  Myabe shold have a seperate function to update
  //               the m_current... varables, instead of doing it in
  //               refreshDecayDisplay(), so that we we could use these
  //               variables in this funcion, or similar function.

  updateInitialMixture();

  const double fracT0 = 0.1;
  const double finalTime = findTimeForActivityFrac( m_currentMixture, fracT0 );
  
  string txt = PhysicalUnitsLocalized::printToBestTimeUnits( finalTime );
  m_displayTimeLength->setText( txt );
}//void DecayActivityDiv::setTimeLimitToDisplay()


void DecayActivityDiv::setDecayChartTimeRange( double finalTime )
{
  string txt = PhysicalUnitsLocalized::printToBestTimeUnits( finalTime );
  m_displayTimeLength->setText( txt );
}//void setDecayChartTimeRange()


void DecayActivityDiv::colorThemeChanged()
{
  m_decayChainChart->colorThemeChanged();
}//void colorThemeChanged();


void DecayActivityDiv::deleteMoreInfoDialog()
{
  if( m_moreInfoDialog )
  {
    delete m_moreInfoDialog;
    m_moreInfoDialog = 0;
  }
}//void deleteMoreInfoDialog()


void DecayActivityDiv::manageActiveDecayChainNucStyling()
{
  const int tab = m_chartTabWidget->currentIndex();
  
  for( const Nuclide &nuc : m_nuclides )
  {
    if( nuc.display && nuc.display->hasStyleClass( "ActiveDecayChainNuc" ) )
      nuc.display->removeStyleClass( "ActiveDecayChainNuc" );
  }//for( const Nuclide &nuc : m_nuclides )
  
  const SandiaDecay::Nuclide *nuc = m_decayChainChart->nuclide();
  if( !nuc || (tab != (m_chartTabWidget->count()-1)) )
    return;
  
  for( const Nuclide &n : m_nuclides )
  {
    if( n.a == nuc->massNumber
        && n.z == nuc->atomicNumber && n.iso == nuc->isomerNumber )
    {
      n.display->addStyleClass( "ActiveDecayChainNuc" );
    }//if( this is a matching nuclide )
  }//for( const Nuclide &n : m_nuclides )
}//void manageActiveDecayChainNucStyling()


void DecayActivityDiv::decayChartClicked( const WMouseEvent& event )
{
  if( (m_currentTimeUnits<0.0)
     || (m_currentTimeRange<0.0) )
    return;

  const WMouseEvent::Coordinates coords = event.widget();
  const WPointF devicePoint( coords.x, coords.y );
  const WPointF cartesianPoint = m_decayChart->mapFromDevice( devicePoint );

  double mouseTime = 0.0;

  if( m_decayChart->type() != Chart::ScatterPlot )
  {
    //xBin goes from 0 to sm_numTimeAxisDivsions, which maps
    // from time 0 to m_currentTimeRange
    
    const double xBin = cartesianPoint.x();
    mouseTime = m_currentTimeRange * xBin / m_currentNumXPoints;
  }else
  {
    const double xPosition = cartesianPoint.x();
    mouseTime = xPosition * m_currentTimeUnits;
  }//if m_decayChart->mapFromDevice() gives the bin number / else actual X value

  displayMoreInfoPopup( mouseTime );
}//void DecayActivityDiv::decayChartClicked( const WMouseEvent& event )



WContainerWidget *DecayActivityDiv::isotopesSummary( const double time ) const
{
  WContainerWidget *cont = new WContainerWidget();
  cont->addStyleClass( "isotopesSummary" );

  if( !m_currentMixture )
    return cont;

  //string mixtureInfo = m_currentMixture->info(time);
  
  
   //Or the equivalent, but not compiling code is:
  string mixtureInfo = [this, time]() -> std::string {
    using namespace SandiaDecay;
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    
      if( !m_currentMixture->numInitialNuclides() )
        return "";
      
      const vector<NuclideActivityPair> activities = m_currentMixture->activity( time );
      const vector<EnergyRatePair> gammasAbundances = m_currentMixture->gammas( time, NuclideMixture::OrderByAbundance, true );
      const vector<EnergyRatePair> alphasAbundances = m_currentMixture->alphas( time );
      const vector<EnergyRatePair> betasAbundances = m_currentMixture->betas( time );
      const vector<EnergyRatePair> betaPlussesAbundances = m_currentMixture->betaPlusses( time );
      
      stringstream infostrm;
    
    auto use_curry = [this,db]( const SandiaDecay::Nuclide * const initial_nuc ) -> bool {
      for( const Nuclide &nuc : m_nuclides )
      {
        const SandiaDecay::Nuclide * const nuclide = db->nuclide( nuc.z, nuc.a, nuc.iso );
        const vector<const SandiaDecay::Nuclide *> children = nuclide->descendants();
        if( std::find(begin(children),end(children),initial_nuc) != end(children) )
          return nuc.useCurie;
      }
      
      return !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    };
    
      infostrm << "<div>Starting from:</div><div style=\"margin-left: 20px; max-width: 60ex;\">";
    
      {//begin codeblock to put starting nuclides information into the stream
        for( int i = 0; i < m_currentMixture->numInitialNuclides(); ++i )
        {
          const SandiaDecay::Nuclide *initial_nuc = m_currentMixture->initialNuclide(i);
          const bool useCurie = use_curry( initial_nuc );
          
          const double initial_act = m_currentMixture->initialActivity(i);
          const string initial_act_str = PhysicalUnits::printToBestActivityUnits(initial_act,2,useCurie,SandiaDecay::becquerel);
          infostrm << (i ? "," : "") << initial_act_str << " " << initial_nuc->symbol;
        }//for( loop over orignal nuclides )
      }//end codeblock to put starting nuclides information into the stream
    
    infostrm << "</div>";
    
      
    infostrm << "<div style=\"margin-top: 10px;\">The following nuclides are present at "
      << PhysicalUnitsLocalized::printToBestTimeUnits(time/SandiaDecay::second,4)
      << " :</div>\n"
      << "<div style=\"margin-left: 20px; max-width: 60ex;\">";
    
    int nparents = 0;
    for( size_t i = 0; i < activities.size(); ++i )
    {
      const NuclideActivityPair &pair = activities[i];
      
      if( IsInf(pair.nuclide->halfLife) || IsNan(pair.nuclide->halfLife) )
        continue;
      
      const bool useCurie = use_curry( pair.nuclide );
      const string act_str = PhysicalUnits::printToBestActivityUnits(pair.activity,2,useCurie,SandiaDecay::becquerel);
      infostrm << (nparents?", ":"") << act_str << " " << pair.nuclide->symbol;
      ++nparents;
    }//for( loop over activities )
    
    infostrm << "</div>\n";
    
      if( gammasAbundances.size() )
      {
        infostrm << "<div style=\"margin-top: 10px;\">Gammas Present:</div>\n"
                 << "<table style=\"margin-left: 20px;\">\n";
        for( size_t i = 0; i < gammasAbundances.size(); ++i )
          infostrm << "<tr><td>" << gammasAbundances[i].energy/SandiaDecay::keV << " keV</td><td>&nbsp;-&nbsp;</td><td>"
          << gammasAbundances[i].numPerSecond << " per second</td></tr>\n";
        infostrm << "</table>\n";
      }//if( gammasAbundances.size() )
      
      if( alphasAbundances.size() )
      {
        infostrm << "<div style=\"margin-top: 10px;\">Alphas Present:</div>\n"
        << "<table style=\"margin-left: 20px;\">\n";
        for( size_t i = 0; i < alphasAbundances.size(); ++i )
          infostrm << "<tr><td>" << alphasAbundances[i].energy/SandiaDecay::keV << " keV</td><td>&nbsp;-&nbsp;</td><td>"
          << alphasAbundances[i].numPerSecond << " per second</td></tr>\n";
        infostrm << "</table>\n";
      }//if( alphasAbundances.size() )
      
      if( betasAbundances.size() )
      {
        infostrm << "<div style=\"margin-top: 10px;\">Betas Present:</div>\n"
        << "<table style=\"margin-left: 20px;\">\n";
        for( size_t i = 0; i < betasAbundances.size(); ++i )
          infostrm << "<tr><td>" << betasAbundances[i].energy/SandiaDecay::keV << " keV</td><td>&nbsp;-&nbsp;</td><td>"
          << betasAbundances[i].numPerSecond << " per second</td></tr>\n";
        infostrm << "</table>\n";
      }//if( alphasAbundances.size() )
      
      
      if( betaPlussesAbundances.size() )
      {
        infostrm << "<div style=\"margin-top: 10px;\">Beta+'s Present:</div>"
        << "<table style=\"margin-left: 20px;\">\n";
        for( size_t i = 0; i < betaPlussesAbundances.size(); ++i )
          infostrm << "<tr><td>" << betaPlussesAbundances[i].energy/SandiaDecay::keV << " keV</td><td>&nbsp;-&nbsp;</td><td>"
          << betaPlussesAbundances[i].numPerSecond << " per second</td></tr>\n";
        infostrm << "</table>\n";
      }//if( betaPlussesAbundances.size() )
      
      return infostrm.str();
    }();
  
   
  mixtureInfo = "<div class=\"NuclideInfo\">\n"
                + mixtureInfo + "\n</div>";
  new WText( mixtureInfo, TextFormat::XHTMLText, cont );

  new WText( "<br /><br /><H3>Further information about the nuclides:</H3>",
             TextFormat::XHTMLText, cont );

  const int nSolutionNuclides = m_currentMixture->numSolutionNuclides();

  for( int i = 0; i < nSolutionNuclides; ++i )
  {
    const SandiaDecay::Nuclide *nuclide = m_currentMixture->solutionNuclide(i);
    WContainerWidget *nucInfo = nuclideInformation( nuclide );
    cont->addWidget( nucInfo );
  }//for( loop over widgets )


  return cont;
}//WContainerWidget *isotopesSummary( const double time ) const


WContainerWidget *DecayActivityDiv::nuclideInformation(
                                    const SandiaDecay::Nuclide *nuclide ) const
{
  WContainerWidget *cont = new WContainerWidget();
  cont->addStyleClass( "nuclideInformation" );

  if( !nuclide )
    return cont;
  
  /*
  stringstream textstr;
  textstr << "<pre class=\"NuclideInfo\">\n";

  if( nuclide )
    textstr << SandiaDecay::human_str_summary(*nuclide) << "\n";

  textstr << "</pre>\n";

  new WText( textstr.str(), XHTMLUnsafeText, cont );
*/
  
  
  stringstream ostr;
  ostr << "<pre class=\"NuclideInfo\">\n";
  ostr << nuclide->symbol << " Atomic Number " << nuclide->atomicNumber <<", Atomic Mass "
  << nuclide->massNumber << ", Isomer Number " << nuclide->isomerNumber << " "
  << nuclide->atomicMass << " AMU, HalfLife=" << PhysicalUnitsLocalized::printToBestTimeUnits(nuclide->halfLife);
  const size_t nParents = nuclide->decaysFromParents.size();
  if( nParents )
    ostr << "\n  Parent";
  if( nParents > 1 )
    ostr << "s";
  if( nParents )
    ostr << ": ";
  
  for( size_t i = 0; i < nParents; ++i )
  {
    if( i )
      ostr << ", ";
    if( nuclide->decaysFromParents[i] && nuclide->decaysFromParents[i]->parent )
      ostr << nuclide->decaysFromParents[i]->parent->symbol;
  }//for( loop over parents )
  
  const size_t nChilds = nuclide->decaysToChildren.size();
  
  for( size_t i = 0; i < nChilds; ++i )
  {
    const auto child = nuclide->decaysToChildren[i];
    //ostr << "\n    " << human_str_summary(*nuclide->decaysToChildren[i]);
    ostr << "Transition ";
    if( child->parent )
      ostr << child->parent->symbol;
    ostr << "&rarr;";
    if( child->child )
      ostr << child->child->symbol;
    else
      ostr << "various";
    
    ostr << ": mode=" << child->mode << " branchRatio=" << child->branchRatio;
    if( child->products.size() )
      ostr << "; Products:";
    for( size_t i = 0; i < child->products.size(); ++i )
    {
      const auto &product = child->products[i];
      ostr << "\n    ";
      
      switch( product.type )
      {
        case SandiaDecay::BetaParticle:            ostr << "beta:"; break;
        case SandiaDecay::GammaParticle:           ostr << "gamma:"; break;
        case SandiaDecay::AlphaParticle:           ostr << "alpha:"; break;
        case SandiaDecay::PositronParticle:        ostr << "positron:"; break;
        case SandiaDecay::CaptureElectronParticle: ostr << "electronCapture:"; break;
        case SandiaDecay::XrayParticle:            ostr << "xray:"; break;
      };//switch( obj.type )
      
      ostr << " energy=" << product.energy << "keV intensity=" << product.intensity;
      switch( product.type )
      {
        case SandiaDecay::GammaParticle:
          break;
          
        case SandiaDecay::AlphaParticle:
          ostr << " hinderence=" << product.hindrance;
          break;
          
        case SandiaDecay::BetaParticle: case SandiaDecay::PositronParticle:
        case SandiaDecay::CaptureElectronParticle:
          ostr << " forbiddenness=" << product.forbiddenness << " logFT=" << product.logFT;
          break;
          
        case SandiaDecay::XrayParticle:
          break;
      };//switch( obj.type )
    }//for( size_t i = 0; i < child->products.size(); ++i )
    
    ostr << "\n";
  }//for( size_t i = 0; i < nChilds; ++i )
  
  ostr << "</pre>\n";
  new WText( ostr.str(), XHTMLUnsafeText, cont );
  
  return cont;
}//WContainerWidget *nuclideInformation( const Nuclide &nuclide ) const;


void DecayActivityDiv::refreshDecayDisplay( const bool update_calc )
{
  checkTimeRangeValid();
  
  bool showsum = true;
  const int noldcol = m_decayModel->columnCount();
  try
  {
    if( noldcol > 2 )
      showsum = boost::any_cast<bool>(m_decayModel->headerData(noldcol - 1,
                                                Wt::Horizontal, Wt::UserRole) );
  }catch( std::exception &e )
  {
    cerr << "refreshDecayDisplay(): unexpectedly caught" << e.what() << endl;
  }
  
  //Grab if we should show the nuclides before clearing the model.
  //  Note that the very first and very last column no not contain nuclide info
  map<string,bool> oldShowNuclide;
  for( int column = 0; column < m_decayModel->columnCount(); ++column )
  {
    const string nuc = m_decayModel->nuclide(column);
    //cout << "Nuc:'" << nuc << "'=" << m_decayModel->showSeries( column ) << endl;
    if( nuc != "" )
      oldShowNuclide[m_decayModel->nuclide(column)] = m_decayModel->showSeries( column );
  }
  
  m_currentNumXPoints = m_decayChart->paintedWidth() / 2;
  
  //Remove the series from last time, and reset the model
#if( WT_VERSION >= 0x3030800 )
  const std::vector<Wt::Chart::WDataSeries *> &series = m_decayChart->series();
  for( size_t i = 0; i < series.size(); ++i )
    m_decayChart->removeSeries( series[i]->modelColumn() );
#else
  const std::vector<Wt::Chart::WDataSeries> &series = m_decayChart->series();
  for( size_t i = 0; i < series.size(); ++i )
    m_decayChart->removeSeries( series[i].modelColumn() );
#endif
  
  m_decayModel->clear();
  m_currentTimeUnits  = -1.0;
  m_currentTimeRange  = -1.0;

  updateInitialMixture();
  
  if( m_nuclides.empty() )
  {
    m_decayLegend->clear();
    return;
  }//if( m_nuclides.empty() )

  const int activityIndex = m_displayActivityUnitsCombo->currentIndex();

  double actunit;
  string unitStr;

  if( activityIndex < static_cast<int>( sm_activityUnitNameValues.size() ) )
  {
    const UnitNameValuePair act_name_unit
                               = sm_activityUnitNameValues.at( activityIndex );
    unitStr = act_name_unit.first;
    actunit = act_name_unit.second;
  }else
  {
    unitStr = "I0";
    actunit = m_currentMixture->totalActivity(0.0);
  }//if( user selected a standard unit ) / else

  const double maxDiplayTime = timeToDisplayTill();
//  const double endActivity = m_currentMixture->totalActivity( maxDiplayTime );
  
  const YAxisType yaxis = YAxisType( m_yAxisType->currentIndex() );

  switch( yaxis )
  {
    case ActivityAxis:
      m_displayActivityUnitsCombo->show();
      m_displayActivityUnitsLabel->show();
    break;
    case GammasAxis: case BetasAxis: case AlphasAxis:
      m_displayActivityUnitsCombo->hide();
      m_displayActivityUnitsLabel->hide();
    break;
    case NumYAxisType:
    break;
  }//switch( yaxis )
  
  
  if( maxDiplayTime <= 0.0 )
    return;
  
  const int nRows = m_currentNumXPoints;
  const vector<SandiaDecay::NuclideTimeEvolution> &origEvolutions
                              = m_currentMixture->decayedToNuclidesEvolutions();
  vector<SandiaDecay::NuclideTimeEvolution> evolutions;
  for( const SandiaDecay::NuclideTimeEvolution &evo : origEvolutions )
  {
    if( evo.nuclide && !IsInf(evo.nuclide->halfLife) && !IsNan( evo.nuclide->halfLife ) )
      evolutions.push_back( evo );
  }//for( const SandiaDecay::NuclideTimeEvolution &evo : origEvolutions )
  
  const int nElements = static_cast<int>( evolutions.size() );

  const double dt = maxDiplayTime / nRows;
  m_decayModel->insertColumns( 0, nElements + 2 );
  m_decayModel->insertRows( 0, nRows );

  double maxActivity = 0.0, minActivity = DBL_MAX;
  vector<double> totalActivities( nRows, 0.0 );

  const UnitNameValuePair xUnitsPair = bestTimeUnit( maxDiplayTime );
  const double tunit = xUnitsPair.second;
  //m_decayChart->axis(Chart::XAxis).setTitle( xUnitsPair.first + " after t\u2080" );
  string xtitle = xUnitsPair.first;
  if( xtitle.size() && islower((int)xtitle[0]) )
    xtitle[0] = (char)toupper( (int)xtitle[0] );
  m_decayChart->setXAxisTitle( xtitle );
  
  for( int row = 0; row < nRows; ++row )
  {
    const double row_time = (row * dt);

    //We will set the x-axis data as a formatted string since
    //  WAxis::setLabelFormat( "%.3g" ); doesnt seem to work
    stringstream labelText;
    labelText << fixed << setprecision(3) << (row * dt)/tunit;
    const WString label( labelText.str() );
    m_decayModel->setData( row, 0, boost::any( label ) );

    for( int elN = 0; elN < nElements; ++elN )
    {
      const int column = elN + 1;
      double yval = 0;
      double activity = evolutions[elN].activity(row_time);

      if( yaxis == ActivityAxis )
      {
        yval = activity / actunit;
      }else
      {
        SandiaDecay::ProductType particletype = SandiaDecay::GammaParticle;
        
        switch( yaxis )
        {
          case ActivityAxis:  case NumYAxisType:                         break;
          case GammasAxis:    particletype = SandiaDecay::GammaParticle; break;
          case BetasAxis:     particletype = SandiaDecay::BetaParticle;  break;
          case AlphasAxis:    particletype = SandiaDecay::AlphaParticle; break;
        }//switch( yaxis )
        
        activity /= SandiaDecay::becquerel;
        for( const SandiaDecay::Transition *trans : evolutions[elN].nuclide->decaysToChildren )
        {
          for( const SandiaDecay::RadParticle &particle : trans->products )
          {
            if( particle.type == particletype )
              yval += activity * trans->branchRatio * particle.intensity;
          }
        }//for( const SandiaDecay::Transition *trans : decays )
      }//if( yaxis == ActivityAxis ) / else

      if( IsInf(yval) || IsNan(yval) )
        continue;

      totalActivities[row] += yval;
      maxActivity = std::max( maxActivity, yval );
      minActivity = std::min( minActivity, yval );
      
#if( WT_VERSION > 0x3040000 )
      // TODO: In Wt 3.7.1, having really small values causes the line to be drawn to really go whack; I'm not sure how to fix this, or where the problem actually comes in.  Maybe the thing to do is just create a D3.js version of this chart
      //if( yval < 10*FLT_EPSILON )
      //{
      //  yval = 0.0;
      //  cout << "[" << row << "," << evolutions[elN].nuclide->symbol << "," << elN << "," << labelText.str() << "]=" << yval << endl;
      //}
#endif
//      if( activity >= (0.00001*endActivity) )
      
      m_decayModel->setData( row, column, boost::any( yval ) );
      
      //WString tt = "My Tool Tip";
      //m_decayModel->setData( row, column, boost::any( tt ), ToolTipRole );
    }//for( loop over nuclides to add )
  }//for( loop over time points to add )

  //Now put in the sum of all the activities
  const WString dateHeader = WString::tr("dad-hdr-date").arg(xUnitsPair.first);
  m_decayModel->setHeaderData( 0, boost::any( dateHeader ) );

  for( int column = 1; column <= nElements; ++column )
  {
    const WString name = evolutions[column-1].nuclide->symbol;
    m_decayModel->setNuclide( column, evolutions[column-1].nuclide->symbol ); //duplicating lots here
    m_decayModel->setHeaderData( column, boost::any( name ) );
    m_decayModel->setHeaderData( column, Wt::Horizontal, true, Wt::UserRole );
  }//for( int column = 0; column < nElements; ++column )

  //Add columns to the chart, and set their header data
  for( int row = 0; row < nRows; ++row )
  {
    maxActivity = max( maxActivity, totalActivities[row] );
    minActivity = std::min( minActivity, totalActivities[row] );
    const double data = totalActivities[row];
    m_decayModel->setData( row, nElements+1, boost::any(data) );
  }//for( loop over rows )
  
  switch( yaxis )
  {
    case ActivityAxis:
      m_decayModel->setHeaderData( nElements+1,
                                   boost::any( WString::tr("dad-hdr-total-act") ) );
      m_decayChart->axis(Chart::YAxis).setTitle( unitStr );
    break;
    case GammasAxis:
      m_decayModel->setHeaderData( nElements+1,
                                   boost::any( WString::tr("dad-hdr-total-gammas") ) );
      m_decayChart->axis(Chart::YAxis).setTitle( WString::tr("dad-yaxis-title-total-gammas") );
    break;
    case BetasAxis:
      m_decayModel->setHeaderData( nElements+1,
                                   boost::any( WString::tr("dad-hdr-total-betas") ) );
      m_decayChart->axis(Chart::YAxis).setTitle( WString::tr("dad-yaxis-title-total-betas") );
    break;
    case AlphasAxis:
      m_decayModel->setHeaderData( nElements+1,
                                   boost::any( WString::tr("dad-hdr-total-alphas") ) );
      m_decayChart->axis(Chart::YAxis).setTitle( WString::tr("dad-yaxis-title-total-alphas") );
    break;
    case NumYAxisType:
    break;
  }//switch( yaxis )
  
  
  m_decayModel->setHeaderData( nElements+1, Wt::Horizontal,
                                 boost::any(showsum), Wt::UserRole );
  std::vector<string> nuclidesset;
  for( int column = 1; column <= nElements; ++column )
  {
    //nuclide also containted in m_decayModel->headerData(column, Wt::Horizontal, Wt::UserRole);
    const string name = evolutions[column-1].nuclide->symbol;
    nuclidesset.push_back( name );
  }
  
  //Now check if we have the same nuclides as before, and if so, set the
  //  visability properties for each of the nuclides charts to be the same
  //  as before.
  bool sameNuclides = (nuclidesset.size() == oldShowNuclide.size());
  for( const string &n : nuclidesset )
    sameNuclides = (sameNuclides && (oldShowNuclide.count(n)>0));
  
  if( sameNuclides )
  {
    for( int column = 0; column < m_decayModel->columnCount(); ++column )
    {
      string nuc = m_decayModel->nuclide(column);
      if( nuc != "" )
        m_decayModel->setShowSeries( column, oldShowNuclide[nuc] );
    }
  }//if( sameNuclides )
  
  
  addDecaySeries();
  
  /*
  {//begin section to add interactive areas
    for( WAbstractArea *area : m_decayChart->areas() )
    {
      m_decayChart->removeArea( area );
      delete area;
    }
    
    const int nColumns = m_decayModel->columnCount();
    for( int column = nColumns-1; column > 0; --column )
    {
    }
  }//end section to add interactive areas
  */
  

  //fill in the 'cache' variables for updateMouseOver(...)
  m_currentTimeUnits  = tunit;
  m_currentTimeRange  = maxDiplayTime;
  
  if( m_calc && update_calc )
    m_calc->setTimeRangeTxt( m_displayTimeLength->text().toUTF8() );
  
#if( WT_VERSION > 0x3040000 )
  m_decayModel->updateMaxActivity();
#endif
  
  updateYAxisRange();
}//void refreshDecayDisplay()


void DecayActivityDiv::updateYAxisRange()
{
  //The next line of letting Wt set y-axis limits doesnt always work super great, so we'll do it
  //  manually
  // m_decayChart->axis(Chart::YAxis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
  
  double miny = std::numeric_limits<double>::max(), maxy = 0.0;
  
  for( int row = 0; row < m_decayModel->rowCount(); ++row )
  {
    for( int col = 1; col < m_decayModel->columnCount(); ++col )
    {
      if( !m_decayModel->showSeries( col ) )
        continue;
      
      boost::any data = m_decayModel->data( m_decayModel->index(row, col) );
      const double yval = asNumber(data);
      if( !IsNan(yval) )
      {
        miny = std::min( miny, yval );
        maxy = std::max( maxy, yval );
      }//if( we are showing this value )
    }//for( loop over columns )
  }//for( loop over rows )

  if( maxy <= 0.0 )
    maxy = 1.0;
  
  if( miny > maxy )
    miny = 0.0;
  
  double min_disp_act = 0.0, max_disp_act = 1.1*maxy;
  switch( m_decayChart->axis(Chart::YAxis).scale() )
  {
    case Wt::Chart::LogScale:
      if( miny < 0.000001*maxy )
        min_disp_act = 0.00001*maxy;
      else
        min_disp_act = 0.75*miny;
      max_disp_act = 1.5*maxy;
      break;
    
    case Wt::Chart::LinearScale:
      //If it wont change the dynamic range of the chart much anyways, might as well
      //  anchor the y-axis to zero
      break;
      
    case Wt::Chart::CategoryScale:
    case Wt::Chart::DateScale:
    case Wt::Chart::DateTimeScale:
      assert( 0 );
      break;
  }//switch( m_decayChart->axis(Chart::YAxis).scale() )
  
  
  m_decayChart->axis(Chart::YAxis).setRange( min_disp_act, max_disp_act );
}//void updateYAxisRange();


void DecayActivityDiv::userSetShowSeries( int series, bool show )
{
  m_decayModel->setShowSeries( series, show );
  updateYAxisRange();
}//void userSetShowSeries()



void DecayActivityDiv::addDecaySeries()
{
#if( WT_VERSION >= 0x3030800 )
  const std::vector<Wt::Chart::WDataSeries *> &series = m_decayChart->series();
  for( size_t i = 0; i < series.size(); ++i )
    m_decayChart->removeSeries( series[i]->modelColumn() );
#else
  const std::vector<Wt::Chart::WDataSeries> &series = m_decayChart->series();
  for( size_t i = 0; i < series.size(); ++i )
    m_decayChart->removeSeries( series[i].modelColumn() );
#endif
  
  const int nColumns = m_decayModel->columnCount();
  
  for( int column = 1; column < nColumns; ++column )
  {
    Chart::WDataSeries series(column, Chart::CurveSeries, Chart::YAxis);
    WPen pen;
    
    auto theme = m_decayChart->colorTheme();
    if( theme && theme->theme_name.toUTF8().find("Default")!=string::npos )
      theme.reset();
    
    if( column != (nColumns-1) )
    {
      if( theme && (theme->referenceLineColor.size()>6) )
      {
        WColor color = theme->referenceLineColor[(column-1) % theme->referenceLineColor.size()];
        pen.setColor( color );
      }else
      {
        pen.setColor( seriesColor( column ) );
      }
    }else
    {
      if( theme && !theme->foregroundLine.isDefault() )
        pen.setColor( theme->foregroundLine );
      pen.setStyle( DashLine );
    }
    
    series.setPen( pen );
    
    m_decayChart->addSeries( series );
  }//for( int column = 1; column < nElements; ++column )
  
  m_decayLegend->clear();
  
  
  for( int column = nColumns-1; column > 0; --column )
  {
    WWidget *item = m_decayChart->createLegendItemWidget( column );
    WContainerWidget *ww = dynamic_cast<WContainerWidget *>( item );
    
    if( nColumns == 3 )  //there is exactly one nuclide
    {
      if( column == (nColumns-1) )  //we're on the Total Activity Entry
      {
        delete item;
        m_decayModel->setShowSeries( column, false );
        continue;
      }
      m_decayModel->setShowSeries( column, true );
    }//if( nColumns == 3 )
    
    if( ww )
    {
      ww->setInline( false );
      WCheckBox *cb = new WCheckBox();
      ww->insertWidget( 0, cb );
      ww->addStyleClass( "LegendEntry" );
      m_decayLegend->addWidget( ww );
      const bool check = m_decayModel->showSeries( column );
      cb->setChecked( check );
      cb->addStyleClass( "ShowSeriesCb" );
      cb->setDisabled( (nColumns == 3) );
      cb->checked().connect( boost::bind( &DecayActivityDiv::userSetShowSeries, this, column, true ) );
      cb->unChecked().connect( boost::bind( &DecayActivityDiv::userSetShowSeries, this, column, false ) );
    }else
    {
      delete item;
      cerr << "DecayActivityDiv::addDecaySeries(): legend items are no longer"
           << " WContainerWidget, please update the code"<< endl;
    }//if( ww )
  }
  
}//void addDecaySeries()


void DecayActivityDiv::updateMouseOver( const Wt::WMouseEvent& event )
{
  if( (m_currentTimeUnits<0.0)
      || (m_currentTimeRange<0.0) )
    return;

  const WMouseEvent::Coordinates coords = event.widget();
  const WPointF devicePoint( coords.x, coords.y );
  const WPointF cartesianPoint = m_decayChart->mapFromDevice( devicePoint );

  double mouseTime = 0.0;

  if( m_decayChart->type() != Chart::ScatterPlot )
  {
  //xBin goes from 0 to sm_numTimeAxisDivsions, which maps
  // from time 0 to m_currentTimeRange
    const double xBin = cartesianPoint.x();
    mouseTime = m_currentTimeRange * xBin / m_currentNumXPoints;
  }else
  {
    const double xPosition = cartesianPoint.x();
    mouseTime = xPosition * m_currentTimeUnits;
  }//if m_decayChart->mapFromDevice() gives us the bin / else actual X value


  const UnitNameValuePair unitpair = bestTimeUnit( mouseTime );
  const double tUnit = unitpair.second;
  const string tUnitStr = unitpair.first;

  const int activityIndex = m_displayActivityUnitsCombo->currentIndex();

  double actunit;
  string actUnitStr;

  if( activityIndex < static_cast<int>( sm_activityUnitNameValues.size() ) )
  {
    const UnitNameValuePair &act_name_unit
                                = sm_activityUnitNameValues.at( activityIndex );
    actUnitStr = act_name_unit.first;
    actunit    = act_name_unit.second;
  }else
  {
    actUnitStr = "I0";
    actunit = m_currentMixture->totalActivity(0.0);
  }//if( user selected a standard unit ) /else


  
  
  stringstream tip;
  tip << WString::tr("dad-chart-click-more-info")
          .arg( SpecUtils::printCompact(mouseTime/tUnit, 4) )
          .arg( tUnitStr ).toUTF8();

  double totalActivity = 0.0;

  const std::vector<SandiaDecay::NuclideTimeEvolution> &evolutions
                              = m_currentMixture->decayedToNuclidesEvolutions();


  const size_t nChildren = evolutions.size();
  for( size_t i = 0; i < nChildren; ++i )
  {
    const size_t index = nChildren - i - 1;
    const SandiaDecay::Nuclide *el = evolutions[index].nuclide;
    const double activity = evolutions[index].activity( mouseTime );
    totalActivity += activity;
    tip << "\n  " << el->symbol << " " << activity/actunit << " " << actUnitStr;
  }//for( loop over population children )

  tip << "\n  " << WString::tr("dad-hdr-total-act").toUTF8() << ": " << totalActivity/actunit << " " << actUnitStr;

  m_decayChart->setToolTip( tip.str() );
}//void updateMouseOver()

DecayActivityDiv::~DecayActivityDiv()
{
  if( m_currentMixture )
    delete m_currentMixture;

  if( m_nuclideSelectDialog )
    delete m_nuclideSelectDialog;
  
  if( m_csvDownloadDialog )
    deleteCsvDownloadGui();
}//DecayActivityDiv destructor
