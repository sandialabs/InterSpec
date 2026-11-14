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
#include <cmath>
#include <cfloat>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <exception>

#include <Wt/WAny.h>
#include <Wt/WPen.h>
#include <Wt/WDate.h>
#include <Wt/WText.h>
#include <Wt/WEvent.h>
#include <Wt/WLabel.h>
#include <Wt/WRectF.h>
#include <Wt/WTable.h>
#include <Wt/WAnchor.h>
#include <Wt/WLength.h>
#include <Wt/WPointF.h>
#include <Wt/WServer.h>
#include <Wt/WPainter.h>
#include <Wt/WSpinBox.h>
#include <Wt/WCheckBox.h>
#include <Wt/WComboBox.h>
#include <Wt/WDateEdit.h>
#include <Wt/WLineEdit.h>
#include <Wt/WResource.h>
#include <Wt/WTabWidget.h>
#include <Wt/WGridLayout.h>
#include <Wt/WPushButton.h>
#include <Wt/Http/Request.h>
#include <Wt/WPaintDevice.h>
#include <Wt/Http/Response.h>
#include <Wt/WStackedWidget.h>
#include <Wt/WContainerWidget.h>
#include <Wt/WRegExpValidator.h>
#include <Wt/WStandardItemModel.h>
#include <Wt/Chart/WCartesianChart.h>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/WidgetUtils.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/DecayChainChart.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DecayActivityDiv.h"
#include "InterSpec/DecayActivityCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DecayBatchCalcWidget.h"
#include "InterSpec/DecaySelectNuclideDiv.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

using namespace Wt;
using namespace std;
using namespace PhysicalUnits;


namespace
{
  string remove_space_copy( string input )
  {
    input.erase( std::remove_if(input.begin(), input.end(), [](char a){return isspace(a);}), input.end() );
    return input;
  }

  // Chart sample count bounds: the painted width is not known before the first paint, and is zero
  //  for a hidden chart.
  const int sm_minChartPoints = 150;
  const int sm_maxChartPoints = 600;

}//namespace

class DecayActivityModel : public Wt::WStandardItemModel
{
  // HACK:
  // Wt has had an issue for really, really small values, where it draws the line at the top of the
  //  chart, when it should be at the bottom.  So we set anything effectively zero, compared to
  //  what we are showing, to zero.
  double m_max_activity;
  
public:
  DecayActivityModel()
    : Wt::WStandardItemModel(),
      m_max_activity( 0.0 )
  {
  }
  
  virtual ~DecayActivityModel()
  {
  }
  
  virtual Wt::cpp17::any data( const Wt::WModelIndex &index,
                          Wt::ItemDataRole role = Wt::ItemDataRole::Display ) const
  {
    if( index.column() == 0 )
      return WStandardItemModel::data( index, role );
    
    Wt::cpp17::any show = headerData( index.column(), Wt::Orientation::Horizontal, Wt::ItemDataRole::User );
    try
    {
      const bool doShow = Wt::cpp17::any_cast<bool>( show );
      if( !doShow )
        return Wt::cpp17::any();
    }catch( std::exception &e )
    {
      cerr << "DecayActivityModel::data(...) unexpectedly caught: " << e.what()
           << endl << index.row() << ", " << index.column()
           << endl;
    }//try / catch
    
    if( role == Wt::ItemDataRole::Display )
    {
      Wt::cpp17::any val = WStandardItemModel::data( index, role );
      if( !val.has_value() )
        return val;
      
      const double valdbl = asNumber(val);
      if( IsNan(valdbl) )
        return val;
      
      // 0.00001 chosen arbitrarily, but values this small wont show up as non-zero on the chart
      if( valdbl < 0.00001*m_max_activity )
        return Wt::cpp17::any( 0.0 );
      
      return val;
    }//if( role == Wt::ItemDataRole::Display )
    
    
    return WStandardItemModel::data( index, role );
  }//data(...)
  
  void setNuclide( int colum, const std::string nuc )
  {
    setHeaderData( colum, Wt::Orientation::Horizontal, Wt::cpp17::any(nuc), Wt::ItemDataRole( Wt::ItemDataRole::User + 1 ) );
  }
  
  std::string nuclide( int colum )
  {
    if( colum <= 0 || colum >= columnCount() )
      return "";
    Wt::cpp17::any val = headerData( colum, Wt::Orientation::Horizontal, Wt::ItemDataRole::User+1 );
    try
    {
      return Wt::cpp17::any_cast<std::string>( val );
    }catch(...)
    {}
    
    return "";
  }
  
  bool showSeries( int colum ) const
  {
    if( colum <= 0 || colum >= columnCount() )
      return true;
    
    Wt::cpp17::any val = headerData( colum, Wt::Orientation::Horizontal, Wt::ItemDataRole::User );
    try
    {
      return Wt::cpp17::any_cast<bool>( val );
    }catch(...)
    {}
    
    return true;
  }//bool showSeries( int colum ) const

  void updateMaxActivity()
  {
    m_max_activity = 0.0;
    
    for( int row = 0; row < rowCount(); ++row )
    {
      for( int col = 1; col < columnCount(); ++col )
      {
        if( !showSeries( col ) )
          continue;
        
        Wt::cpp17::any val = WStandardItemModel::data( index(row, col) );
        if( !val.has_value() )
          continue;
        
        const double valdbl = asNumber(val);
        if( IsNan(valdbl) )
          continue;
        
        m_max_activity = std::max( m_max_activity, valdbl );
      }
    }
  }//void updateMaxActivity()
  
  void setShowSeries( int colum, bool show )
  {
    const bool oldval = showSeries( colum );
    if( oldval == show )
      return;
    
    if( colum < 1 || colum >= columnCount() )
      return;
    
    setHeaderData( colum, Wt::Orientation::Horizontal, Wt::cpp17::any(show), Wt::ItemDataRole::User );
    
    updateMaxActivity();
    
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
  std::shared_ptr<const ColorTheme> m_colorTheme;
  WBrush m_chartMarginBrush;
  WPen m_textPen;
  
public:
  DecayActivityChart()
  : Wt::Chart::WCartesianChart(),
    m_widthInPixels(250),
    m_heightInPixels(250),
    m_phone( false ),
    m_textPen( WColor(Wt::StandardColor::Black) )
  {
    setStyleClass( "DecayActivityChart" );
    setSeriesSelectionEnabled(true);
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
  
  /** Applies the spectrum-chart colors of a theme; a null or "Default" theme restores the defaults. */
  void setColorTheme( std::shared_ptr<const ColorTheme> theme )
  {
    m_colorTheme = theme;

    const bool useTheme = theme && (theme->theme_name.toUTF8().find("Default") == string::npos);
    auto themeColor = [useTheme, &theme]( WColor ColorTheme::*member ) -> WColor {
      return useTheme ? ((*theme).*member) : WColor();
    };

    const WColor textColor = themeColor( &ColorTheme::spectrumChartText );
    m_textPen = WPen( textColor.isDefault() ? WColor(Wt::StandardColor::Black) : textColor );
    setTextPen( m_textPen );

    const WColor axisColor = themeColor( &ColorTheme::spectrumAxisLines );
    const WPen axisPen( axisColor.isDefault() ? WColor(Wt::StandardColor::Black) : axisColor );
    axis(Chart::Axis::X).setPen( axisPen );
    axis(Chart::Axis::Y).setPen( axisPen );

    const WColor marginColor = themeColor( &ColorTheme::spectrumChartMargins );
    m_chartMarginBrush = marginColor.isDefault() ? WBrush() : WBrush( marginColor );

    // Without a theme, Wts default white; a theme with no chart background gets none.
    const WColor backgroundColor = themeColor( &ColorTheme::spectrumChartBackground );
    if( !backgroundColor.isDefault() )
      setBackground( WBrush(backgroundColor) );
    else
      setBackground( useTheme ? WBrush() : WBrush( WColor(Wt::StandardColor::White) ) );
  }//setColorTheme()
  
  
  void setXAxisTitle( const std::string &xtitle )
  {
    m_xAxisTitle = xtitle;
    if( !m_phone )
      axis(Chart::Axis::X).setTitle( xtitle );
  }
  
  int paintedWidth() const { return static_cast<int>(m_widthInPixels); }
  int paintedHeight() const { return static_cast<int>(m_heightInPixels); }
  
  void showGridLines( const bool draw = true )
  {
    Chart::WCartesianChart::axis(Chart::Axis::X).setGridLinesEnabled( draw );
    Chart::WCartesianChart::axis(Chart::Axis::Y).setGridLinesEnabled( draw );
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
    const WRectF rect = (rectangle.isNull() || rectangle.isEmpty()) ? painter.window() : rectangle;
    const double width = rect.width();
    const double height = rect.height();

    // Fill the margins around the plot area first; the base class then fills the plot area itself
    //  with `background()`, and draws on top.
    //  (MakeDrfChart has a copy of this margin painting.)
    if( m_chartMarginBrush.style() != Wt::BrushStyle::None )
    {
      const double padLeft = plotAreaPadding(Wt::Side::Left);
      const double padRight = plotAreaPadding(Wt::Side::Right);
      const double padTop = plotAreaPadding(Wt::Side::Top);
      const double padBottom = plotAreaPadding(Wt::Side::Bottom);

      const bool vertical = (orientation() == Wt::Orientation::Vertical);
      const double w = vertical ? width : height, h = vertical ? height : width;
      const double rx = vertical ? padLeft : padTop;
      const double ry = vertical ? padTop : padRight;
      const double rw = std::max( 10.0, w - (vertical ? (padLeft + padRight) : (padTop + padBottom)) );
      const double rh = std::max( 10.0, h - (vertical ? (padTop + padBottom) : (padRight + padLeft)) );

      if( ry > 0 )  //Top strip
        painter.fillRect( hv(WRectF(rx, 0, rw, ry)), m_chartMarginBrush );
      if( (ry + rh) < h ) //Bottom strip
        painter.fillRect( hv(WRectF(rx, ry + rh, rw, h - ry - rh)), m_chartMarginBrush );
      if( rx > 0 ) //Left strip
        painter.fillRect( hv(WRectF(0, 0, rx, h)), m_chartMarginBrush );
      if( (rx + rw) < w ) //Right strip
        painter.fillRect( hv(WRectF(rx + rw, 0, w - rw - rx, h)), m_chartMarginBrush );
    }//if( we have a margin color )

    Chart::WCartesianChart::paint( painter, rectangle );
    
    painter.save();
    //TODO 20110325: for chrome the below seems to have no effect, however
    //               it works on FireFox
    painter.rotate( -90 );
    painter.setPen( m_textPen );
    painter.setBrush( WBrush() );
    if( !axis(Chart::Axis::Y).title().empty() )
      painter.drawText( -0.45*height, 0.0, 0.2, 0.1,
                       Wt::AlignmentFlag::Center, axis(Chart::Axis::Y).title() );
    painter.restore();
    
    if( m_phone )
    {
      //We will make the x-axis title "compact".  This has the downside that
      //  since we arent manually drawing x-axis labels (Wt is doing it), we
      //  dont know if we
      const double w = 7*m_xAxisTitle.size();
      const double x = width - w - 5 - plotAreaPadding(Wt::Side::Right);
      const double y = height-15;
      WBrush fillbrush = background();
      if( fillbrush.color().isDefault() )
        fillbrush = WBrush( WColor(120,120,120) );
      
      painter.fillRect( x-5, y-2.5, w+12, 17, fillbrush );
      painter.setPen( m_textPen );
      painter.drawText( x, y, w, 12, Wt::AlignmentFlag::Middle, m_xAxisTitle );
    }//if( m_phone )
  }//paint( ... )
};//class DecayActivityChart




class DateLengthCalculator : public WContainerWidget
{
  /** Time spans longer than this are not put into the date fields (a WDate can not represent
   them, and the day count would overflow an int).
   */
  static constexpr double sm_maxDateSpanDays = 1.0E6;

  DecayActivityDiv *m_activityDiv;
  WGridLayout *m_layout;
  WDateEdit *m_begindate;
  WDateEdit *m_enddate;
  WLineEdit *m_duration;
  WText *m_error;

  /** Holds #m_inputs and #m_results; hidden while there is no valid time span. */
  WContainerWidget *m_info;

  /** Activity and age inputs for each nuclide.  Only rebuilt when the set of nuclides changes
   elsewhere, so they are not destroyed while the user is editing, or tabbing between, them.
   */
  WContainerWidget *m_inputs;

  /** The "Initial Activity" / "Final Activity" label of each nuclide in #m_inputs. */
  std::vector<WLabel *> m_activityLabels;

  /** Whether #m_inputs needs rebuilding. */
  bool m_inputsStale;

  /** Tables of the activities and masses before/after the time span. */
  WContainerWidget *m_results;
  
public:
  DateLengthCalculator( DecayActivityDiv *activityDiv )
   : WContainerWidget(),
     m_activityDiv( activityDiv ),
     m_layout( nullptr ),
     m_begindate( nullptr ),
     m_enddate( nullptr ),
     m_duration( nullptr ),
     m_error( nullptr ),
     m_info( nullptr ),
     m_inputs( nullptr ),
     m_activityLabels{},
     m_inputsStale( true ),
     m_results( nullptr )
  {
    assert( m_activityDiv );
    
    int ndaysBack = 364;
    const double currentDur = activityDiv->timeToDisplayTill();
    if( currentDur > 0 && currentDur < 100*PhysicalUnits::year )
      ndaysBack = (int)floor(0.5 + currentDur / PhysicalUnits::day);
    
    m_layout = setLayout( make_unique<WGridLayout>() );

    // The date format is localizable; the placeholder shows it in lower case (e.g., "mm/dd/yyyy").
    const WString dateFormat = WString::tr( "dad-calc-date-format" );
    const WString datePlaceholder = WString::fromUTF8( SpecUtils::to_lower_ascii_copy( dateFormat.toUTF8() ) );

    WLabel *label = m_layout->addWidget( make_unique<WLabel>( WString::tr("dad-calc-initial-date") ),
                                         0, 0, Wt::AlignmentFlag::Middle );

    m_begindate = m_layout->addWidget( make_unique<WDateEdit>(), 0, 1, 1, 2, Wt::AlignmentFlag::Left );
    label->setBuddy( m_begindate );
    m_begindate->setFormat( dateFormat );
    m_begindate->setPlaceholderText( datePlaceholder );
    m_begindate->setDate( Wt::WDate::currentServerDate().addDays(-ndaysBack) );

    label = m_layout->addWidget( make_unique<WLabel>( WString::tr("dad-calc-end-date") ),
                                 1, 0, Wt::AlignmentFlag::Middle );
    m_enddate = m_layout->addWidget( make_unique<WDateEdit>(), 1, 1, 1, 2, Wt::AlignmentFlag::Left );
    label->setBuddy( m_enddate );
    m_enddate->setFormat( dateFormat );
    m_enddate->setPlaceholderText( datePlaceholder );
    // \TODO: consider setting to current foreground spectrum date
    m_enddate->setDate( Wt::WDate::currentServerDate() );

    m_begindate->changed().connect( this, &DateLengthCalculator::dateChanged );
    m_enddate->changed().connect( this, &DateLengthCalculator::dateChanged );

    label = m_layout->addWidget( make_unique<WLabel>( WString::tr("dad-calc-time-span") ),
                                 2, 0, Wt::AlignmentFlag::Middle );

    m_duration = m_layout->addWidget( make_unique<WLineEdit>(), 2, 1, Wt::AlignmentFlag::Left );

#if( BUILD_AS_OSX_APP || IOS )
    m_duration->setAttributeValue( "autocorrect", "off" );
    m_duration->setAttributeValue( "spellcheck", "off" );
#endif
    label->setBuddy( m_duration );

    auto validator = make_shared<WRegExpValidator>( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalPosOrNegRegex() );
    validator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
    m_duration->setValidator( validator );
    m_duration->changed().connect( this, &DateLengthCalculator::durationChanged );
    m_duration->enterPressed().connect( this, &DateLengthCalculator::durationChanged );
    
    WPushButton *batchButton = m_layout->addWidget( make_unique<WPushButton>( WString::tr("dad-batch-decay") ),
                                                    2, 2, Wt::AlignmentFlag::Right );
    batchButton->setStyleClass( "LinkBtn" );
    batchButton->clicked().connect( [this](){ m_activityDiv->openBatchDecayTool(); } );

    m_error = m_layout->addWidget( make_unique<WText>(), 3, 0, 1, 3 );
    m_error->setInline( false );
    m_error->hide();

    m_info = m_layout->addWidget( make_unique<WContainerWidget>(), 4, 0, 1, 3 );
    m_inputs = m_info->addNew<WContainerWidget>();
    m_results = m_info->addNew<WContainerWidget>();
    m_results->addStyleClass( "SelectableTxt" );
    m_info->hide();
    
    m_layout->setRowStretch( 3, 1 );
    m_layout->setRowStretch( 4, 1 );
    
    m_layout->setColumnStretch( 2, 1 );
  }//DateLengthCalculator
  
  
  /** Call when the source nuclides change, other than through this calculator. */
  void setInputsStale()
  {
    m_inputsStale = true;
  }
  
  
  /** Updates the results, and the inputs if the nuclides changed. */
  void update()
  {
    if( m_inputsStale )
      rebuildInputs();
    updateResults();
  }//void update()
  
  
  /** Sets the time span from the chart (always positive); keeps the direction the user chose. */
  void setTimeRangeTxt( const std::string &txt )
  {
    const string oldtxt = m_duration->text().toUTF8();
    if( txt != oldtxt )
    {
      m_duration->setText( WString::fromUTF8(txt) );
      if( oldtxt.empty() || (oldtxt[0] != '-') || m_begindate->text().empty() )
        m_begindate->setText( "" );
      else
        m_enddate->setText( "" );
    }
  }//void setTimeRangeTxt( std::string txt )
  
  
private:
  void showError( const WString &msg )
  {
    m_error->setText( msg );
    m_error->show();
    m_duration->setText( "" );
    m_info->hide();
  }//void showError( const WString &msg )
  
  
  /** Returns 0.0 on error.
   Return a negative time period if we want to calculate an initial activity.
   */
  double getValidatedTimeSpan()
  {
    //Let be optimistic things will be okay
    m_error->hide();
    m_error->setText( "" );
    
    const std::string durtxt = m_duration->text().toUTF8();
    const bool validBeginDate = (m_begindate->validate() == Wt::ValidationState::Valid);
    const bool validEndDate = (m_enddate->validate() == Wt::ValidationState::Valid);
    bool duration_valid = false;
    double duration = 0.0;
    try
    {
      duration = m_activityDiv->parseTimeSpan( durtxt );
      duration_valid = true;
    }catch( std::exception & )
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
      {
        showError( WString::tr("dad-calc-err-begin-date") );
        return 0.0;
      }
      
      if( !validEndDate )
      {
        showError( WString::tr("dad-calc-err-end-date") );
        return 0.0;
      }
      
      const int days = m_begindate->date().daysTo(m_enddate->date());
      duration = 24.0*3600.0*days*PhysicalUnits::second;
      m_duration->setText( PhysicalUnitsLocalized::printToBestTimeUnits( duration, 2 ) );
      
      return duration;
    }else if( (!validBeginDate && !validEndDate)
              || (fabs(duration/PhysicalUnits::day) > sm_maxDateSpanDays) )
    {
      return duration;
    }else if( !validBeginDate )
    {
      const int ndays = (int)floor( 0.5 + duration/PhysicalUnits::day );
      m_begindate->setDate( m_enddate->date().addDays( -ndays ) );
      return duration;
    }
    
    assert( !validEndDate );
    const int ndays = (int)floor( 0.5 + duration/PhysicalUnits::day );
    m_enddate->setDate( m_begindate->date().addDays( ndays ) );
    return duration;
  }//double getValidatedTimeSpan()
  
  
  /** Puts the (absolute) time span into the chart, which redraws, and then updates us. */
  void pushCurrentToParent()
  {
    string txt = m_duration->text().toUTF8();
    if( txt.length() && (txt[0] == '-') )
      txt = txt.substr(1);
    
    m_activityDiv->m_displayTimeLength->setText( WString::fromUTF8(txt) );
    m_activityDiv->refreshDecayDisplay( false );
  }//void pushCurrentToParent()
  
  
  void dateChanged()
  {
    m_duration->setText( "" );
    
    if( getValidatedTimeSpan() == 0.0 )
    {
      updateResults();
      return;
    }
    
    const int days = m_begindate->date().daysTo(m_enddate->date());
    const string datestr = PhysicalUnitsLocalized::printToBestTimeUnits( 24.0*3600.0*days*PhysicalUnits::second, 2 );
    m_duration->setText( WString::fromUTF8(datestr) );
    
    pushCurrentToParent();
  }//void dateChanged()
  
  
  void durationChanged()
  {
    // We get here when the user edits the duration text
    const string durtxt = m_duration->text().toUTF8();
    if( !durtxt.empty() && (m_duration->validate() == WValidator::State::Valid) )
    {
      try
      {
        const double duration = m_activityDiv->parseTimeSpan( durtxt );
        
        double mintime = 1.0 * SandiaDecay::second;
        for( const DecayActivityDiv::Nuclide &n : m_activityDiv->m_nuclides )
        {
          if( n.nuclide )
            mintime = std::min( mintime, n.nuclide->halfLife );
        }
        
        if( fabs(duration) < 0.01*mintime )
        {
          m_duration->setText( m_activityDiv->m_displayTimeLength->text() );
          updateResults();
          return;
        }
        
        int64_t ndays = std::floor(fabs(duration) / PhysicalUnits::day);
        double remanderSeconds = fabs(duration) - (ndays * PhysicalUnits::day);
        if( fabs(PhysicalUnits::day - remanderSeconds) <= 1.0 )
        {
          remanderSeconds = PhysicalUnits::day - remanderSeconds;
          ndays += 1;
        }
        
        // \TODO: do the same thing as days, but for years to try and set an initial date
        if( duration < 0.0 )
          ndays = -ndays;

        if( (std::abs(ndays) < 1) || (remanderSeconds > 1.0) || (std::abs(ndays) > sm_maxDateSpanDays) )
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
      }catch( std::exception & )
      {
        m_begindate->setText( "" );
      }
    }//if( we have valid duration string )
    
    if( getValidatedTimeSpan() == 0.0 )
    {
      updateResults();
      return;
    }
    
    pushCurrentToParent();
  }//void durationChanged()
  
  
  /** Adds an activity or age input, for the nuclide at `index`, to a table cell. */
  WLineEdit *addNuclideInput( WTableCell *cell, const WString &regex, const string &txt )
  {
    WLineEdit *edit = cell->addNew<WLineEdit>();
#if( BUILD_AS_OSX_APP || IOS )
    edit->setAttributeValue( "autocorrect", "off" );
    edit->setAttributeValue( "spellcheck", "off" );
#endif
    auto validator = make_shared<WRegExpValidator>( regex );
    validator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
    edit->setValidator( validator );
    edit->setTextSize( 10 );
    edit->setText( WString::fromUTF8(txt) );
    
    return edit;
  }//addNuclideInput(...)
  
  
  void rebuildInputs()
  {
    m_inputsStale = false;
    m_inputs->clear();
    m_activityLabels.clear();
    
    const std::vector<DecayActivityDiv::Nuclide> &nucs = m_activityDiv->m_nuclides;
    
    for( size_t index = 0; index < nucs.size(); ++index )
    {
      const DecayActivityDiv::Nuclide &nucinfo = nucs[index];
      if( !nucinfo.nuclide )
        continue;
      
      WTable *nuctbl = m_inputs->addNew<WTable>();
      nuctbl->setMargin( 10, Wt::Side::Top );

      const bool showNucHeader = (nucs.size() > 1);
      const int rowOffset = showNucHeader ? 1 : 0;

      if( showNucHeader )
      {
        nuctbl->setHeaderCount( 1 );
        WTableCell *cell = nuctbl->elementAt(0, 0);
        cell->addNew<WLabel>( WString::tr("dad-calc-for-nuc").arg(nucinfo.nuclide->symbol) );
        cell->setColumnSpan( 2 );
        cell->setAttributeValue( "style", "text-align: left;" + cell->attributeValue("style") );
      }//if( showNucHeader )

      WTableCell *cell = nuctbl->elementAt(rowOffset, 0);
      cell->setVerticalAlignment( Wt::AlignmentFlag::Middle );
      WLabel *activityLabel = cell->addNew<WLabel>( WString::tr("dad-calc-initial-act") );
      activityLabel->setMargin( 4, Wt::Side::Right );
      m_activityLabels.push_back( activityLabel );

      WLineEdit *activityEdit = addNuclideInput( nuctbl->elementAt(rowOffset, 1),
                                                 PhysicalUnits::sm_activityRegex, nucinfo.activityStr );
      activityLabel->setBuddy( activityEdit );
      
      // The edit is addressed by index (not a copy of the nuclide), as the nuclides change.
      auto doActivityUpdate = [this, index, activityEdit](){
        std::vector<DecayActivityDiv::Nuclide> &nucs = m_activityDiv->m_nuclides;
        if( index >= nucs.size() )
          return;
        
        DecayActivityDiv::Nuclide &nuc = nucs[index];
        const string input_act_str = activityEdit->text().toUTF8();
        
        double activity = -1.0;
        try
        {
          if( activityEdit->validate() == WValidator::State::Valid )
            activity = PhysicalUnits::stringToActivity( input_act_str );
        }catch( std::exception & )
        {
        }
        
        if( (activity <= 0.0) || IsNan(activity) || IsInf(activity) )
        {
          activityEdit->setText( WString::fromUTF8(nuc.activityStr) );
          return;
        }
        
        if( (activity == nuc.activity) && (input_act_str == nuc.activityStr) )
          return;
        
        nuc.activity = activity;
        nuc.activityStr = input_act_str;
        nuc.useCurie = DecayActivityCalc::activityStringIsCurie( input_act_str ).value_or( nuc.useCurie );
        nuc.updateTxt();
        
        m_activityDiv->refreshDecayDisplay( false );
      };//doActivityUpdate lambda
      
      activityEdit->changed().connect( activityEdit, doActivityUpdate );
      activityEdit->enterPressed().connect( activityEdit, doActivityUpdate );
      
      
      cell = nuctbl->elementAt(1 + rowOffset, 0);
      cell->setVerticalAlignment( Wt::AlignmentFlag::Middle );
      WLabel *ageLabel = cell->addNew<WLabel>( WString::tr("dad-calc-initial-age") );
      ageLabel->setMargin( 4, Wt::Side::Right );
      
      const string ageTxt = PhysicalUnitsLocalized::printToBestTimeUnits( nucinfo.age );
      WLineEdit *ageEdit = addNuclideInput( nuctbl->elementAt(1 + rowOffset, 1),
                                            PhysicalUnitsLocalized::timeDurationRegex(), ageTxt );
      ageLabel->setBuddy( ageEdit );
      
      auto doAgeUpdate = [this, index, ageEdit](){
        std::vector<DecayActivityDiv::Nuclide> &nucs = m_activityDiv->m_nuclides;
        if( index >= nucs.size() )
          return;
        
        DecayActivityDiv::Nuclide &nuc = nucs[index];
        
        double age = -1.0;
        try
        {
          if( ageEdit->validate() == WValidator::State::Valid )
            age = PhysicalUnitsLocalized::stringToTimeDuration( ageEdit->text().toUTF8() );
        }catch( std::exception & )
        {
        }
        
        if( (age < 0.0) || IsNan(age) || IsInf(age) )
        {
          ageEdit->setText( WString::fromUTF8( PhysicalUnitsLocalized::printToBestTimeUnits(nuc.age) ) );
          return;
        }
        
        if( age == nuc.age )
          return;
        
        nuc.age = age;
        nuc.updateTxt();
        
        m_activityDiv->refreshDecayDisplay( false );
      };//doAgeUpdate lambda
      
      ageEdit->changed().connect( ageEdit, doAgeUpdate );
      ageEdit->enterPressed().connect( ageEdit, doAgeUpdate );
    }//for( size_t index = 0; index < nucs.size(); ++index )
  }//void rebuildInputs()
  
  
  /** Rebuilds the result tables. */
  void updateResults()
  {
    m_results->clear();
    
    const double timeSpan = getValidatedTimeSpan();
    const double absTimeSpan = fabs( timeSpan );
    
    if( timeSpan == 0.0 )
    {
      m_info->hide();
      return;
    }
    
    m_info->show();
    
    for( WLabel *label : m_activityLabels )
      label->setText( WString::tr( (timeSpan < 0.0) ? "dad-calc-final-act" : "dad-calc-initial-act" ) );
    
    const std::vector<DecayActivityDiv::Nuclide> &nucs = m_activityDiv->m_nuclides;
    if( nucs.empty() )
      return;
    
    // For a negative time span the entered activities are at the end, so we make our own mixture,
    //  with the activities they must have started at.
    SandiaDecay::NuclideMixture *mix = m_activityDiv->m_currentMixture.get();
    std::unique_ptr<SandiaDecay::NuclideMixture> local_mix;
    if( timeSpan < 0.0 )
    {
      local_mix = std::make_unique<SandiaDecay::NuclideMixture>();
      mix = local_mix.get();
    }
    
    if( !mix )
      return;
    
    int nSigFigsActivity = 0;
    
    for( const DecayActivityDiv::Nuclide &nucinfo : nucs )
    {
      const SandiaDecay::Nuclide * const nuc = nucinfo.nuclide;
      if( !nuc )
        continue;
      
      int nThisSigFig = 0;
      for( const char c : nucinfo.activityStr )
        nThisSigFig += ((c >= '0') && (c <= '9'));
      nSigFigsActivity = std::max( nThisSigFig, nSigFigsActivity );
      
      if( timeSpan < 0.0 )
      {
        // TODO: estimate how many half-lives is reasonable to calculate this initial activity for.
        const double decrease_factor = std::exp( timeSpan * nuc->decayConstant() );
        const double initial_activity = nucinfo.activity / decrease_factor;
        
        if( IsNan(initial_activity)
            || IsInf(initial_activity)
            || (initial_activity < std::numeric_limits<float>::epsilon()) )
        {
          showError( WString::tr("dad-calc-err-too-long")
                      .arg( PhysicalUnitsLocalized::printToBestTimeUnits(absTimeSpan) )
                      .arg( nuc->symbol )
                      .arg( PhysicalUnitsLocalized::printToBestTimeUnits(nuc->halfLife) ) );
          return;
        }//if( probably too many half-lives to back calculate or zero initial activity )
        
        if( nucinfo.age > DBL_EPSILON )
          local_mix->addAgedNuclideByActivity( nuc, initial_activity, nucinfo.age );
        else
          local_mix->addNuclideByActivity( nuc, initial_activity );
      }//if( timeSpan < 0.0 )
    }//for( const DecayActivityDiv::Nuclide &nucinfo : nucs )
    
    const string durationTxt = PhysicalUnitsLocalized::printToBestTimeUnits( absTimeSpan );
    
    // What we actually want is number of sig-figs, but since the print function
    //  take number of digits after decimal, we'll just assume a single significant
    //  digit before the decimal, and then print `nSigFigsActivity` digits after
    //  (i.e., one more sig-fig than input)
    const int maxDecimal = std::max(2, nSigFigsActivity);
    
    // Adds a table of the nuclides activities and masses at `time`; if `onlyInitial`, only those
    //  nuclides the mixture started with.
    auto addTable = [&]( const WString &heading, const double time, const bool onlyInitial,
                         const int nDecimal, const string &styleClass ){
      WText *line = m_results->addNew<WText>( "&nbsp;" );
      line->setInline( false );
      line = m_results->addNew<WText>( heading );
      line->setInline( false );
      
      WTable *table = m_results->addNew<WTable>();
      table->addStyleClass( styleClass );
      table->setHeaderCount( 1, Wt::Orientation::Horizontal );
      table->elementAt(0, 0)->addNew<WLabel>( WString::tr("Nuclide") );
      table->elementAt(0, 1)->addNew<WLabel>( WString::tr("Activity") );
      table->elementAt(0, 2)->addNew<WLabel>( WString::tr("Mass") );
      
      const vector<SandiaDecay::NuclideActivityPair> activities = mix->activity( time );
      const vector<SandiaDecay::NuclideNumAtomsPair> numatoms = mix->numAtoms( time );
      
      int row = 1;
      for( size_t i = 0; i < activities.size(); ++i )
      {
        const SandiaDecay::NuclideActivityPair &nap = activities[activities.size() - 1 - i];
        if( !nap.nuclide )
          continue;
        
        if( onlyInitial )
        {
          bool is_orig = false;
          for( int j = 0; !is_orig && (j < mix->numInitialNuclides()); ++j )
            is_orig = (mix->initialNuclide(j) == nap.nuclide);
          if( !is_orig )
            continue;
        }//if( onlyInitial )
        
        table->elementAt(row, 0)->addNew<WText>( nap.nuclide->symbol );
        
        const bool useCuries = m_activityDiv->useCurieFor( nap.nuclide );
        const WString actstr = IsInf(nap.nuclide->halfLife)
                ? WString::tr("stable")
                : WString::fromUTF8( PhysicalUnits::printToBestActivityUnits( nap.activity, nDecimal, useCuries ) );
        table->elementAt(row, 1)->addNew<WText>( actstr );
        
        for( const SandiaDecay::NuclideNumAtomsPair &natoms : numatoms )
        {
          if( natoms.nuclide != nap.nuclide )
            continue;
          
          const double ngrams = natoms.numAtoms / nap.nuclide->atomsPerGram();
          const string massstr = PhysicalUnits::printToBestMassUnits( ngrams*PhysicalUnits::gram, nDecimal );
          table->elementAt(row, 2)->addNew<WText>( massstr );
          break;
        }//for( loop over number of atoms )
        
        row += 1;
      }//for( loop over activities )
    };//addTable lambda
    
    if( timeSpan < 0.0 )
      addTable( WString::tr("dad-calc-before").arg(durationTxt), 0.0, true, 3,
                "BeforeDecayNucInfoTable DecayedNucInfoTable" );
    
    addTable( WString::tr("dad-calc-after").arg(durationTxt), absTimeSpan, false, maxDecimal,
              "DecayedNucInfoTable" );
  }//void updateResults()
};//class DateLengthCalculator

namespace
{
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
                          const bool isPhone )
    : WContainerWidget(),
      m_chart( chart ),
      m_legend( legend )
    {
      m_chart = addWidget( unique_ptr<DecayActivityChart>(chart) );
      m_legend = addWidget( unique_ptr<WContainerWidget>(legend) );

      const string extraheight = isPhone ? "20" : "100";

      const string js =
      """function(self, w, h,layout) {"
      ""  "if (!self.wtWidth || self.wtWidth!=w "
      ""      "|| !self.wtHeight || self.wtHeight!=h) {"
      ""    "self.wtWidth=w; self.wtHeight=h;"
      ""    "self.style.height=h + 'px';"
      ""  "}"
      ""  "var child = document.getElementById('" + m_chart->id() + "');"
      ""  "if(child && child.wtResize){"
      ""     "child.wtResize(child,w,h,layout);"
      ""  "}else{ console.log('ChartAndLegendHolder: couldnt resize chart, please look into.'); }"
      ""  "var leg=document.getElementById('" + m_legend->id() + "');"
      ""  "if(leg)leg.style.maxHeight=(h-" + extraheight + ")+'px';"
      "" "}";
      setJavaScriptMember( WT_RESIZE_JS, js );
    }//ChartAndLegendHolder constructor
    
    
    DecayActivityChart *m_chart;
    WContainerWidget *m_legend;
  };//class ChartAndLegendHolder
}//namespace


// Not in the anonymous namespace, as DecayActivityDiv befriends it.
class DecayCsvResource : public Wt::WResource
{
  //ToDo: Decouple from DecayCsvResource class and put in own source/header,
  //      or maybe with CsvDownloadGui
public:
  DecayCsvResource( DecayActivityDiv *display )
    : WResource(),
      m_display( display ),
      m_displayHandle( display ),
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
    // Have Wt hold the session lock across handleRequest instead of dropping it for us.
    //  By default (`takesUpdateLock() == false`) `WResource::handle` calls `handler->unlock()`
    //  before invoking us, so the session thread can be tearing this resource down - inside
    //  `~WResource`'s `beingDeleted()`, which waits for our use-count to drop - while we are
    //  blocked acquiring the very lock it holds.  That is a deadlock, and it is why the
    //  `WApplication::UpdateLock` below had to exist at all.  Costs nothing: that lock already
    //  spanned the whole of handleRequest, so this only moves the acquisition earlier.
    setTakesUpdateLock( true );

    if( m_display )
      m_timeSpan = m_display->timeToDisplayTill();
    
    //When a download is complete, we will delete the dialog, but jic something wacky (and
    //  unexpected) is going on, we'll wrap it so that it wont get called if the DecayActivityDiv
    //  has been deleted.
    //  (I'm a little adverse to relying on a 'modal' user effect to enforce object lifetimes)
    //
    //  This used to be `wApp->bind(...)`, which Wt 4 removed.  It cannot capture `this` or
    //  `m_display`: the task is posted from `handleRequest` (which runs on a request/IO thread)
    //  and `WServer::post` only guarantees the *session* is still alive when it runs - by then
    //  both the resource (a child of the transient CSV dialog) and the DecayActivityDiv (inside a
    //  transient DecayWindow) may be gone.  A `WidgetHandle` is just an inert DOM-id string, so it
    //  is safe to copy across the thread boundary, unlike an observing_ptr.
    const WidgetUtils::WidgetHandle handle = m_displayHandle;
    m_cleanup = [handle](){
      DecayActivityDiv * const div = handle.resolve_as<DecayActivityDiv>();
      if( div )
        div->deleteCsvDownloadGuiTriggerUpdate();
    };
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
  // WSpinBox::valueChanged() passes on what the user typed, so clamp rather than throw.
  void setNumberRows( const int nrow ){ m_numRows = static_cast<size_t>( std::clamp(nrow, 1, 10000) ); }
  
private:
  virtual void handleRequest( const Wt::Http::Request &request,
                              Wt::Http::Response &response )
  {
    // With setTakesUpdateLock(true) (see ctor) Wt already holds the session lock for us here, so
    //  this is a same-thread no-op (WApplication::UpdateLock detects that case and returns
    //  immediately).  Kept so the null-session path below still has something to test.
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
      // Re-resolve rather than trusting `m_display`: this runs on a request thread, and while the
      //  UpdateLock we hold keeps the session thread out for the duration, the widget could
      //  already have been destroyed before we got here.  `findById` is legal under the lock.
      DecayActivityDiv * const display = m_displayHandle.resolve_as<DecayActivityDiv>();
      if( !display )
        throw runtime_error( "No data avalaiable" );

      SandiaDecay::NuclideMixture *mixture = display->m_currentMixture.get();
      if( !mixture )
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
      
      string name;
      bool use_curies = false;
      for( const DecayActivityDiv::Nuclide &n : display->m_nuclides )
      {
        const SandiaDecay::Nuclide * const nuc = n.nuclide;
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

      suggestFileName( name, ContentDisposition::Attachment );
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
        const double time_now = timeSpan * ((m_numRows > 1) ? row : 1u) / numerator;
        
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
      suggestFileName( "no_data.csv", ContentDisposition::Attachment );
      response.out() << eol_char << e.what() << eol_char;
    }
    
    WServer::instance()->post( m_sessionId, m_cleanup );
  }//handleRequest(...)
  
  /** Only valid during construction (which happens on the session thread) - use #m_displayHandle
   from `handleRequest` and from anything posted back to the session. */
  DecayActivityDiv *m_display;
  /** Id-based handle to #m_display, for re-resolving it off the session thread; see #m_cleanup. */
  const WidgetUtils::WidgetHandle m_displayHandle;
  size_t m_numRows;
  double m_timeSpan;
  bool m_activity;
  bool m_xrays;
  bool m_gammas;
  bool m_alphas;
  bool m_betas;
  std::string m_sessionId;
  Wt::WApplication *m_app; //it looks like WApplication::instance() will be valid in handleRequest, but JIC
  std::function<void()> m_cleanup;
};//class DecayCsvResource : public Wt::WResource


class CsvDownloadGui : public AuxWindow
{
  friend class AuxWindow;

  DecayActivityDiv *m_parent;
  DecayCsvResource *m_csvResouce;
  WLineEdit *m_ageEdit;

protected:
  // Constructor is protected; use AuxWindow::make<CsvDownloadGui>() to create.
  CsvDownloadGui( DecayActivityDiv *parent )
  : AuxWindow( WString::tr("dad-csv-window-title"), (AuxWindowProperties::IsModal
                              | AuxWindowProperties::TabletNotFullScreen
                              | AuxWindowProperties::DisableCollapse
                              | AuxWindowProperties::SetCloseable) ),
  m_parent( parent ),
  m_csvResouce( nullptr ),
  m_ageEdit( nullptr )
  {
    assert( m_parent );

    setResizable( false );
    finished().connect( m_parent, &DecayActivityDiv::deleteCsvDownloadGui );
    
    m_csvResouce = addChild( make_unique<DecayCsvResource>( parent ) );
    
    const double displayTime = parent->timeToDisplayTill();
    const string timeSpanStr = PhysicalUnitsLocalized::printToBestTimeUnits( displayTime );
    
    m_csvResouce->setTimeSpan( displayTime );
    m_csvResouce->setNumberRows( 400 );
    m_csvResouce->setGiveActivities( true );
    m_csvResouce->setGiveXrays( false );
    m_csvResouce->setGiveGammas( false );
    m_csvResouce->setGiveAlphas( false );
    m_csvResouce->setGiveBetas( false );
    
    WTable *table = contents()->addNew<WTable>();

    WTableCell *el = table->elementAt(0,0);
    table->setMargin( 8, Wt::Side::Top | Wt::Side::Left | Wt::Side::Right );
    WLabel *label = el->addNew<WLabel>( WString::tr("dad-time-span-label") );
    label->addStyleClass( "DecayChartTimeSpanLabel" );
    el->setVerticalAlignment( Wt::AlignmentFlag::Middle );
    label->setMargin( 2, Wt::Side::Right );

    el = table->elementAt(0,1);
    m_ageEdit = el->addNew<WLineEdit>( timeSpanStr );
#if( BUILD_AS_OSX_APP || IOS )
    m_ageEdit->setAttributeValue( "autocorrect", "off" );
    m_ageEdit->setAttributeValue( "spellcheck", "off" );
#endif
    label->setToolTip( WString::tr("dad-csv-tt-time-span") );
    label->setBuddy( m_ageEdit );
    auto validator = make_shared<WRegExpValidator>( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex() );
    validator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
    m_ageEdit->setValidator( validator );
    m_ageEdit->changed().connect( this, &CsvDownloadGui::updateTimeSpan );
    m_ageEdit->blurred().connect( this, &CsvDownloadGui::updateTimeSpan );
    m_ageEdit->enterPressed().connect( this, &CsvDownloadGui::updateTimeSpan );
    
    el = table->elementAt(0,2);
    label = el->addNew<WLabel>( WString::tr("dad-csv-num-rows") );
    el->setVerticalAlignment( Wt::AlignmentFlag::Middle );
    label->setMargin( 20, Wt::Side::Left );
    label->setMargin( 2, Wt::Side::Right );
    el = table->elementAt(0,3);
    WSpinBox *numPointsSB = el->addNew<WSpinBox>();
    label->setBuddy( numPointsSB );
    label->setToolTip( WString::tr("dad-csv-tt-num-rows") );
    numPointsSB->setValue( 400 );
    numPointsSB->setMinimum(1);
    numPointsSB->setMaximum( 10000 );
    numPointsSB->valueChanged().connect( m_csvResouce, &DecayCsvResource::setNumberRows );
    
    table = contents()->addNew<WTable>();
    table->setMargin( 8 );
    el = table->elementAt(0,0);
    el->addNew<WLabel>( WString::tr("dad-csv-include") );

    el = table->elementAt(0,1);
    WCheckBox *cb = el->addNew<WCheckBox>( WString::tr("dad-csv-cb-act") );
    cb->setMargin( 5, Wt::Side::Left );
    cb->setChecked( true );
    cb->addStyleClass( "CbNoLineBreak" );
    cb->checked().connect( this, [this](){ m_csvResouce->setGiveActivities( true ); } );
    cb->unChecked().connect( this, [this](){ m_csvResouce->setGiveActivities( false ); } );

    el = table->elementAt(0,2);
    cb = el->addNew<WCheckBox>( WString::tr("dad-csv-cb-xrays") );
    cb->setMargin( 5, Wt::Side::Left );
    cb->addStyleClass( "CbNoLineBreak" );
    cb->checked().connect( this, [this](){ m_csvResouce->setGiveXrays( true ); } );
    cb->unChecked().connect( this, [this](){ m_csvResouce->setGiveXrays( false ); } );

    el = table->elementAt(0,3);
    cb = el->addNew<WCheckBox>( WString::tr("dad-csv-cb-gammas") );
    cb->setMargin( 5, Wt::Side::Left );
    cb->addStyleClass( "CbNoLineBreak" );
    cb->checked().connect( this, [this](){ m_csvResouce->setGiveGammas( true ); } );
    cb->unChecked().connect( this, [this](){ m_csvResouce->setGiveGammas( false ); } );

    el = table->elementAt(0,4);
    cb = el->addNew<WCheckBox>( WString::tr("dad-csv-cb-alphas") );
    cb->setMargin( 5, Wt::Side::Left );
    cb->addStyleClass( "CbNoLineBreak" );
    cb->checked().connect( this, [this](){ m_csvResouce->setGiveAlphas( true ); } );
    cb->unChecked().connect( this, [this](){ m_csvResouce->setGiveAlphas( false ); } );

    el = table->elementAt(0,5);
    cb = el->addNew<WCheckBox>( WString::tr("dad-csv-cb-betas") );
    cb->setMargin( 5, Wt::Side::Left );
    cb->addStyleClass( "CbNoLineBreak" );
    cb->checked().connect( this, [this](){ m_csvResouce->setGiveBetas( true ); } );
    cb->unChecked().connect( this, [this](){ m_csvResouce->setGiveBetas( false ); } );
    
    WLink csvLink( WLink::Type::Url, m_csvResouce->url() );
    csvLink.setTarget( Wt::LinkTarget::NewWindow );
    WAnchor *csvancor = footer()->addNew<WAnchor>();
    csvancor->setLink( csvLink );
    csvancor->setText( WString::tr("dad-csv-export") );
    csvancor->setPadding( 8, Wt::Side::Top );
    csvancor->setPadding( 2, Wt::Side::Bottom );
    csvancor->setStyleClass( "ExportCsvLink" );
    footer()->setStyleClass( "modal-footer" );
    footer()->setContentAlignment( Wt::AlignmentFlag::Center );
    
    show();
    centerWindow();
  }//CsvDownloadGui

public:
  void updateTimeSpan()
  {
    double timespan = -1.0;
    try
    {
      timespan = m_parent->parseTimeSpan( m_ageEdit->valueText().toUTF8() );
    }catch( std::exception & )
    {
    }
    
    if( (timespan <= 0.0) || IsInf(timespan) || IsNan(timespan) )
    {
      timespan = m_parent->timeToDisplayTill();
      m_ageEdit->setText( WString::fromUTF8( PhysicalUnitsLocalized::printToBestTimeUnits(timespan) ) );
    }
    
    m_csvResouce->setTimeSpan( timespan );
  }//void updateTimeSpan()
  
};//class CsvDownloadGui

namespace
{
  /** Name of the selected display unit: an entry of `PhysicalUnits::sm_activityUnitNameValues`,
   or `DecayActivityCalc::sm_arbitraryDisplayUnit` (the last combo entry).
   */
  string displayUnitName( const WComboBox *combo )
  {
    const int index = combo->currentIndex();
    if( (index >= 0) && (index < static_cast<int>(PhysicalUnits::sm_activityUnitNameValues.size())) )
      return PhysicalUnits::sm_activityUnitNameValues[index].first;
    return DecayActivityCalc::sm_arbitraryDisplayUnit;
  }//displayUnitName(...)


  /** Selects a display unit by exact name; an empty name selects MBq or mCi. */
  void setDisplayUnit( WComboBox *combo, string unit, const bool useBq )
  {
    if( unit.empty() )
      unit = useBq ? "MBq" : "mCi";

    const PhysicalUnits::UnitNameValuePairV &units = PhysicalUnits::sm_activityUnitNameValues;
    int index = static_cast<int>( units.size() ); //The "Arbitrary" entry
    for( size_t i = 0; i < units.size(); ++i )
    {
      if( units[i].first == unit )
        index = static_cast<int>( i );
    }

    combo->setCurrentIndex( index );
  }//setDisplayUnit(...)
}//namespace


DecayActivityDiv::DecayActivityDiv( InterSpec *viewer )
: WContainerWidget(),
  m_viewer( viewer ),
  m_nuclides{},
  m_displayTimeLength( nullptr ),
  m_displayActivityUnitsCombo( nullptr ),
  m_displayActivityUnitsLabel( nullptr ),
  m_logYScale( nullptr ),
  m_showGridLines( nullptr ),
  m_yAxisType( nullptr ),
  m_parentNuclidesDiv( nullptr ),
  m_nuclidesAddedDiv( nullptr ),
  m_noNuclidesTxt( nullptr ),
  m_createNewNuclideButton( nullptr ),
  m_clearNuclidesButton( nullptr ),
  m_nuclideSelectDialog( nullptr ),
  m_nuclideSelect( nullptr ),
  m_editingNuclide( nullptr ),
  m_chartTabWidget( nullptr ),
  m_decayChart( nullptr ),
  m_decayModel( nullptr ),
  m_decayChainChart( nullptr ),
  m_moreInfoDialog( nullptr ),
  m_decayLegend( nullptr ),
  m_calc( nullptr ),
  m_calcNeedsUpdate( true ),
  m_csvDownloadDialog( nullptr ),
  m_currentTimeUnits( -1.0 ),
  m_currentMixture( std::make_unique<SandiaDecay::NuclideMixture>() )
{
  addStyleClass( "DecayActivityDiv" );
  
  wApp->useStyleSheet( "InterSpec_resources/NuclideDecayInfo.css");
  
  if( m_viewer )
    m_viewer->useMessageResourceBundle( "DecayActivity" );
  
  init();
}//DecayActivityDiv constructor


void DecayActivityDiv::init()
{
  //First delete any widgets that may be owned my *this
  WContainerWidget::clear();

  const bool isPhone = m_viewer ? m_viewer->isPhone() : false;
  
  //Initialize all member widgets; will assign parentage later on
  m_displayTimeLength      = new WLineEdit( "" );
#if( BUILD_AS_OSX_APP || IOS )
  m_displayTimeLength->setAttributeValue( "autocorrect", "off" );
  m_displayTimeLength->setAttributeValue( "spellcheck", "off" );
#endif
  
  auto validator = make_shared<WRegExpValidator>( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex() );
  validator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
  m_displayTimeLength->setValidator( validator );

  m_displayActivityUnitsCombo      = new WComboBox();
  m_displayActivityUnitsLabel      = new WLabel( WString::tr(isPhone ? "dad-units-label-phone" : "dad-units-label") );
  m_logYScale                      = new WCheckBox( WString::tr("dad-logy-cb") );
  m_showGridLines                  = new WCheckBox( WString::tr("dad-grid-cb") );
  m_yAxisType                      = new WComboBox();
  m_parentNuclidesDiv              = new WContainerWidget();
  m_nuclidesAddedDiv               = new WContainerWidget();
  m_createNewNuclideButton         = new WPushButton();
  m_clearNuclidesButton            = new WPushButton();
  m_nuclideSelectDialog            = AuxWindow::make( WString::tr("dad-sel-nuc-window-title"),
                                    (AuxWindowProperties::TabletNotFullScreen
                                     | AuxWindowProperties::DisableCollapse
                                     | AuxWindowProperties::StartHidden) );
  
  m_nuclideSelect                  = new DecaySelectNuclide( isPhone, m_nuclideSelectDialog.get() );
  m_decayLegend                    = new WContainerWidget();

  m_chartTabWidget                 = new WTabWidget();
  m_decayChart                     = new DecayActivityChart();
  m_decayModel                     = make_shared<DecayActivityModel>();

  //Set elements names so we can style with css
  m_parentNuclidesDiv->addStyleClass( "m_parentNuclidesDiv" );
  m_nuclidesAddedDiv->addStyleClass( "m_nuclidesAddedDiv" );

  m_createNewNuclideButton->setIcon( "InterSpec_resources/images/plus_min_black.svg" );
  m_createNewNuclideButton->addStyleClass( "m_createNewNuclideButton" );
  
  m_clearNuclidesButton->setIcon( "InterSpec_resources/images/remove_all_black.png" );
  m_clearNuclidesButton->addStyleClass( "m_clearNuclidesButton" );

  // The buttons show only their icon now, so the label each one used to carry becomes its tooltip.
  m_createNewNuclideButton->setToolTip( WString::tr("dad-add-nucs") );
  m_clearNuclidesButton->setToolTip( WString::tr("dad-remove-all") );
    
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

  //setup m_parentNuclidesDiv: the nuclides sit in their own outlined area, and the two actions
  //  (add / remove-all) are grouped together at the right of it.
  m_parentNuclidesDiv->setInline( false );
  m_parentNuclidesDiv->addWidget( unique_ptr<WContainerWidget>(m_nuclidesAddedDiv) );

  WContainerWidget *nucActions = m_parentNuclidesDiv->addNew<WContainerWidget>();
  nucActions->addStyleClass( "NucActions" );
  nucActions->addWidget( unique_ptr<WPushButton>(m_createNewNuclideButton) );
  nucActions->addWidget( unique_ptr<WPushButton>(m_clearNuclidesButton) );

  m_nuclidesAddedDiv->setInline( false );

  // Shown only while no nuclide has been added, so the empty area says what it is for.
  m_noNuclidesTxt = m_nuclidesAddedDiv->addNew<WText>( WString::tr("dad-no-nuclides") );
  m_noNuclidesTxt->addStyleClass( "NoNuclidesTxt" );
  m_nuclideSelectDialog->contents()->addWidget( unique_ptr<DecaySelectNuclide>(m_nuclideSelect) );
  if( m_viewer && (m_viewer->renderedHeight() > 100) )
    m_nuclideSelectDialog->setMaximumSize( WLength::Auto, 0.9*m_viewer->renderedHeight() );
  
  m_nuclideSelectDialog->rejectWhenEscapePressed();
  m_createNewNuclideButton->clicked().connect( this, &DecayActivityDiv::showAddNuclideDialog );

  m_nuclideSelect->done().connect( m_nuclideSelectDialog.get(), &AuxWindow::hide );
  m_nuclideSelectDialog->finished().connect( this, &DecayActivityDiv::nuclideSelectDialogDone );
  // Dialog starts hidden via AuxWindowProperties::StartHidden — no explicit hide() here.
  
  m_nuclideSelect->selected().connect( this, [this]( const NuclideSelectedInfo &info ){ addTheNuclide( info ); } );
  m_clearNuclidesButton->clicked().connect( this, &DecayActivityDiv::confirmClearAllNuclides );
  m_decayChart->clicked().connect( this, &DecayActivityDiv::decayChartClicked );

  auto decayDivOwner = make_unique<WContainerWidget>();
  WContainerWidget *decayDiv = decayDivOwner.get();
  decayDiv->addStyleClass( "DecayActivityChartDiv" ); /* TODO: Hacked to avoid dealing with color theme for the moment */
  decayDiv->addStyleClass( "decayDiv" );
  WGridLayout *decayLayout = decayDiv->setLayout( make_unique<WGridLayout>() );

  decayLayout->addWidget( make_unique<ChartAndLegendHolder>( m_decayChart, m_decayLegend, isPhone ), 0, 0, 1, 1 );
  decayLayout->addWidget( unique_ptr<WContainerWidget>(displayOptionsDiv), 1, 0, 1, 1 );
  decayLayout->setRowStretch( 0, 5 );

  //////////////////
  // Decay Chain
  m_decayChainChart = new DecayChainChart();
  auto decayChainDivOwner = make_unique<WContainerWidget>();
  WContainerWidget *decayChainDiv = decayChainDivOwner.get();
  decayChainDiv->setMinimumSize( 300, 350 );
  decayChainDiv->addStyleClass( "decayChainDiv" ); //ToDo: Currently decayChainDiv is hacked to avoid dealing with color theme ish
  decayChainDiv->setOverflow( Wt::Overflow::Auto );
  WGridLayout *decayChainLayout = decayChainDiv->setLayout( make_unique<WGridLayout>() );

  decayChainLayout->addWidget( unique_ptr<DecayChainChart>(m_decayChainChart), 0, 0, 1, 1 );
  decayChainLayout->setRowStretch( 0, 1 );

  //////////////////
  // tab widget; the tab order must match DecayActivityCalc::Tab
  const Wt::ContentLoading loadPolicy = Wt::ContentLoading::Lazy;
  m_chartTabWidget->addTab( std::move(decayDivOwner), WString::tr("dad-mi-act-chart"), loadPolicy );
  m_chartTabWidget->addTab( std::move(decayChainDivOwner), WString::tr("dad-mi-decay-chain"), loadPolicy );
  {
    auto calcOwner = make_unique<DateLengthCalculator>( this );
    m_calc = calcOwner.get();
    m_chartTabWidget->addTab( std::move(calcOwner), WString::tr("dad-mi-calc"), loadPolicy );
  }
  assert( m_chartTabWidget->count() == 3 );

  m_chartTabWidget->currentChanged().connect( m_decayChainChart, &DecayChainChart::deleteMoreInfoDialog );
  m_chartTabWidget->currentChanged().connect( this, &DecayActivityDiv::deleteMoreInfoDialog );
  m_chartTabWidget->currentChanged().connect( this, &DecayActivityDiv::manageActiveDecayChainNucStyling );
  m_chartTabWidget->currentChanged().connect( this, &DecayActivityDiv::chartTabChanged );
  m_decayChainChart->nuclideChanged().connect( this, &DecayActivityDiv::manageActiveDecayChainNucStyling );

  WGridLayout *layout = setLayout( make_unique<WGridLayout>() );
  layout->setContentsMargins( 3, 3, 3, 3 );
  layout->addWidget( unique_ptr<WContainerWidget>(m_parentNuclidesDiv), 0, 0, 1, 1 );
  layout->addWidget( unique_ptr<WTabWidget>(m_chartTabWidget), 1, 0, 1, 1 );
  layout->setRowStretch( 1, 5 );
  
  if( m_viewer )
    m_viewer->colorThemeChanged().connect( this, &DecayActivityDiv::colorThemeChanged );
}//void DecayActivityDiv::init()


void DecayActivityDiv::handleAppUrl( std::string path, std::string query_str )
{
  InterSpec * const viewer = m_viewer ? m_viewer : InterSpec::instance();
  const bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", viewer );
  
  vector<WString> errors;
  DecayActivityCalc::ToolState state;
  try
  {
    state = DecayActivityCalc::fromAppUrl( path, query_str, useBq, errors );
  }catch( std::exception &e )
  {
    passMessage( WString::fromUTF8(e.what()), WarningWidget::WarningMsgHigh );
    return;
  }
  
  deleteMoreInfoDialog();
  m_decayChainChart->deleteMoreInfoDialog();
  
  clearNuclidesNoRefresh();
  for( const DecayActivityCalc::SourceNuclide &src : state.nuclides )
    addNuclideNoRefresh( src );
  updateNuclideAreaState();
  if( m_calc )
    m_calc->setInputsStale();
  
  const Nuclide * const lastNuc = m_nuclides.empty() ? nullptr : &m_nuclides.back();
  m_decayChainChart->setNuclide( (lastNuc ? lastNuc->nuclide : nullptr),
                                 (lastNuc ? lastNuc->useCurie : !useBq),
                                 DecayChainChart::DecayChainType::DecayFrom );
  
  m_showGridLines->setChecked( state.gridLines );
  setGridLineStatus();
  m_logYScale->setChecked( state.logY );
  applyYScale();
  m_yAxisType->setCurrentIndex( static_cast<int>(state.chartType) );
  setDisplayUnit( m_displayActivityUnitsCombo, state.displayUnit, useBq );
  
  if( state.timeSpan.empty() )
    setTimeLimitToDisplay();
  else
    m_displayTimeLength->setText( WString::fromUTF8(state.timeSpan) );
  
  refreshDecayDisplay( true );
  
  // After the refresh, so the calculator (if shown) uses the new nuclides
  m_chartTabWidget->setCurrentIndex( static_cast<int>(state.tab) );
  manageActiveDecayChainNucStyling();
  
  // Display up to the first 5 (arbitrary) error messages
  for( size_t i = 0; (i < 5) && (i < errors.size()); ++i )
    passMessage( errors[i], WarningWidget::WarningMsgHigh );
}//void handleAppUrl( std::string path, std::string query_str )


std::string DecayActivityDiv::encodeStateToUrl()
{
  DecayActivityCalc::ToolState state;
  
  const int tab = m_chartTabWidget->currentIndex();
  assert( (tab >= 0) && (tab <= static_cast<int>(DecayActivityCalc::Tab::Calc)) );
  if( (tab >= 0) && (tab <= static_cast<int>(DecayActivityCalc::Tab::Calc)) )
    state.tab = static_cast<DecayActivityCalc::Tab>( tab );
  
  state.gridLines = m_showGridLines->isChecked();
  state.logY = m_logYScale->isChecked();
  
  const int chartType = m_yAxisType->currentIndex();
  if( (chartType >= 0) && (chartType <= static_cast<int>(DecayActivityCalc::ChartType::Alpha)) )
    state.chartType = static_cast<DecayActivityCalc::ChartType>( chartType );
  
  state.displayUnit = displayUnitName( m_displayActivityUnitsCombo );
  
  // Written in English units (not the users locale), so the URL works anywhere.
  const double timespan = timeToDisplayTill();
  if( timespan > FLT_EPSILON )
    state.timeSpan = DecayActivityCalc::printTimeForUrl( timespan );
  
  for( const Nuclide &nuc : m_nuclides )
    state.nuclides.push_back( nuc );
  
  return DecayActivityCalc::toAppUrl( state );
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
  applyYScale();
  m_logYScale->checked().connect( this, &DecayActivityDiv::updateYScale );
  m_logYScale->unChecked().connect( this, &DecayActivityDiv::updateYScale );
  m_showGridLines->checked().connect(this, &DecayActivityDiv::setGridLineStatus );
  m_showGridLines->unChecked().connect(this, &DecayActivityDiv::setGridLineStatus );

  for( const PhysicalUnits::UnitNameValuePair &unit : PhysicalUnits::sm_activityUnitNameValues )
    m_displayActivityUnitsCombo->addItem( unit.first );
  m_displayActivityUnitsCombo->addItem( WString::tr("dad-arbitrary") );

  const bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", m_viewer ? m_viewer : InterSpec::instance() );
  setDisplayUnit( m_displayActivityUnitsCombo, "", useBq );

  WContainerWidget *displayOptionsDiv = new WContainerWidget();
  displayOptionsDiv->addStyleClass( "displayOptionsDiv" );
  WContainerWidget *displOptUpper = displayOptionsDiv->addNew<WContainerWidget>();

  WLabel *endLabel = displOptUpper->addNew<WLabel>( WString::tr("dad-time-span-label") );
  endLabel->setBuddy( m_displayTimeLength );
  displOptUpper->addWidget( unique_ptr<WLineEdit>(m_displayTimeLength) );
  m_displayTimeLength->changed().connect( this, [this](){ refreshDecayDisplay( true ); } );
  m_displayTimeLength->enterPressed().connect( this, [this](){ refreshDecayDisplay( true ); } );
  
  const bool showToolTips = m_viewer ? UserPreferences::preferenceValue<bool>( "ShowTooltips", m_viewer ) : false;
  if( m_viewer )
    HelpSystem::attachToolTipOn( m_displayTimeLength, WString::tr("dad-tt-age"), showToolTips );
  
  WLabel *yaxisTypeLabel = displOptUpper->addNew<WLabel>( WString::tr(isPhone ? "dad-yaxis-label-phone" : "dad-yaxis-label") );
  yaxisTypeLabel->addStyleClass( "DecayChartYaxisTypeLabel" );

  // Item order must match DecayActivityCalc::ChartType
  displOptUpper->addWidget( unique_ptr<WComboBox>(m_yAxisType) );
  m_yAxisType->addItem( WString::tr("Activity") );
  m_yAxisType->addItem( WString::tr("dad-gamma-rate") );
  m_yAxisType->addItem( WString::tr("dad-beta-rate") );
  m_yAxisType->addItem( WString::tr("dad-alpha-rate") );
  
  //If phone, skip creating displOptLower to make things more compact
  WContainerWidget *displOptLower = nullptr;
  if( isPhone )
  {
    displOptLower = displOptUpper;
  }else
  {
    displOptLower = displayOptionsDiv->addNew<WContainerWidget>();
    displOptLower->setMargin( 5, Wt::Side::Top );
    displOptLower->addStyleClass( "LowerOptions" );
  }

  displOptLower->addWidget( unique_ptr<WLabel>(m_displayActivityUnitsLabel) );
  displOptLower->addWidget( unique_ptr<WComboBox>(m_displayActivityUnitsCombo) );

  displOptLower->addWidget( unique_ptr<WCheckBox>(m_logYScale) );
  displOptLower->addWidget( unique_ptr<WCheckBox>(m_showGridLines) );
  
  if( isPhone )
  {
    m_logYScale->hide();
    m_showGridLines->hide();
  }
  
  m_yAxisType->setCurrentIndex( static_cast<int>(DecayActivityCalc::ChartType::Activity) );
  m_yAxisType->activated().connect( this, [this]( int ){ refreshDecayDisplay( true ); } );
  
  if( !isPhone )
  {
    WContainerWidget *spacer = displOptLower->addNew<WContainerWidget>();
    spacer->addStyleClass( "LowerOptionsSpacer" );

    WPushButton *batchButton = displOptLower->addNew<WPushButton>( WString::tr("dad-batch-decay") );
    batchButton->setStyleClass( "LinkBtn" );
    batchButton->clicked().connect( this, &DecayActivityDiv::openBatchDecayTool );

    WPushButton *csvButton = displOptLower->addNew<WPushButton>();
    csvButton->setIcon( "InterSpec_resources/images/download_small.svg" );
    csvButton->setText( WString("{1}...").arg( WString::tr("CSV") ) );
    csvButton->setStyleClass( "LinkBtn DownloadBtn" );
    csvButton->clicked().connect( this, &DecayActivityDiv::createCsvDownloadGui );
  }

  m_displayActivityUnitsCombo->changed().connect( this, [this](){ refreshDecayDisplay( true ); } );

  return displayOptionsDiv;
}//Wt::WContainerWidget *initDisplayOptionWidgets()


void DecayActivityDiv::createCsvDownloadGui()
{
  deleteCsvDownloadGui();
  m_csvDownloadDialog = AuxWindow::make<CsvDownloadGui>( this );
}//void createCsvDownloadGui()


void DecayActivityDiv::openBatchDecayTool()
{
  InterSpec *viewer = m_viewer ? m_viewer : InterSpec::instance();
  if( !viewer )
    return;

  DecayBatchCalcWindow *win = viewer->createDecayBatchCalcWindow();
  if( !win )
    return;

  win->clearNuclides();

  for( const Nuclide &nuc : m_nuclides )
  {
    const SandiaDecay::Nuclide * const n = nuc.nuclide;
    if( n )
      win->addNuclide( n->atomicNumber, n->massNumber, n->isomerNumber, nuc.activity, nuc.useCurie,
                       nuc.age, nuc.activityStr );
  }//for( const Nuclide &nuc : m_nuclides )
}//void openBatchDecayTool()


void DecayActivityDiv::deleteCsvDownloadGui()
{
  if( m_csvDownloadDialog )
    AuxWindow::deleteAuxWindow( m_csvDownloadDialog.get() );
  assert( !m_csvDownloadDialog );
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

  m_decayChart->setModel( m_decayModel );
  m_decayChart->setXSeriesColumn( 0 );
  m_decayChart->setType( Chart::ChartType::Scatter );
  m_decayChart->axis(Chart::Axis::X).setScale( Chart::AxisScale::Linear );

  //TODO 20110325: The either of next line doesnt seem to have any effect
  m_decayChart->axis(Chart::Axis::X).setLabelFormat( "%.3g" );

  if( m_viewer && m_viewer->isPhone() )
  {
    m_decayChart->setIsPhone();
    m_decayChart->setPlotAreaPadding( 75, Wt::Side::Left );
    m_decayChart->setPlotAreaPadding( 25, Wt::Side::Bottom );
    m_decayChart->setPlotAreaPadding( 20, Wt::Side::Right );
    m_decayChart->setPlotAreaPadding( 0, Wt::Side::Top );
  }else
  {
    m_decayChart->setPlotAreaPadding( 80, Wt::Side::Left );
    m_decayChart->setPlotAreaPadding( 50, Wt::Side::Bottom );
    m_decayChart->setPlotAreaPadding( 20, Wt::Side::Right );
    m_decayChart->setPlotAreaPadding( 5, Wt::Side::Top );
  }
  
  if( m_viewer )
    m_decayChart->setColorTheme( m_viewer->getColorTheme() );

  //do not display the legend to the right of the chat
  m_decayChart->setLegendEnabled( false);
  
  m_decayChart->setToolTip( WString::tr("dad-click-for-more-info") );
}//void initCharts()


double DecayActivityDiv::firstNuclideHalfLife() const
{
  if( m_nuclides.empty() || !m_nuclides.front().nuclide )
    return -1.0;
  return m_nuclides.front().nuclide->halfLife;
}//double firstNuclideHalfLife() const


double DecayActivityDiv::parseTimeSpan( const std::string &txt ) const
{
  return PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( txt, firstNuclideHalfLife() );
}//double parseTimeSpan( const std::string &txt ) const


bool DecayActivityDiv::useCurieFor( const SandiaDecay::Nuclide *nuc ) const
{
  for( const Nuclide &src : m_nuclides )
  {
    if( !src.nuclide || !nuc )
      continue;
    
    if( src.nuclide == nuc )
      return src.useCurie;
    
    const vector<const SandiaDecay::Nuclide *> descendants = src.nuclide->descendants();
    if( std::find( begin(descendants), end(descendants), nuc ) != end(descendants) )
      return src.useCurie;
  }//for( const Nuclide &src : m_nuclides )
  
  InterSpec * const viewer = m_viewer ? m_viewer : InterSpec::instance();
  return !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", viewer );
}//bool useCurieFor( const SandiaDecay::Nuclide *nuc ) const


double DecayActivityDiv::timeToDisplayTill()
{
  try
  {
    return parseTimeSpan( m_displayTimeLength->text().toUTF8() );
  }catch( std::exception & )
  {
  }
  
  const double finalTime = DecayActivityCalc::findTimeForActivityFrac( *m_currentMixture, 0.1 );
  m_displayTimeLength->setText( WString::fromUTF8( PhysicalUnitsLocalized::printToBestTimeUnits(finalTime) ) );
  return finalTime;
}//double DecayActivityDiv::timeToDisplayTill()


void DecayActivityDiv::showAddNuclideDialog()
{
  clearEditingNuclide();
  m_nuclideSelectDialog->show();
  m_nuclideSelectDialog->centerWindow();
  m_nuclideSelect->setNuclideSearchToFocus();
}//void showAddNuclideDialog()


void DecayActivityDiv::clearEditingNuclide()
{
  if( m_editingNuclide )
    m_editingNuclide->removeStyleClass( "EditingNuclide" );
  m_editingNuclide = nullptr;
  
  if( m_nuclideSelect )
    m_nuclideSelect->setAddButtonToAdd();
}//void clearEditingNuclide()


void DecayActivityDiv::nuclideSelectDialogDone()
{
  clearEditingNuclide();
}//nuclideSelectDialogDone()


void DecayActivityDiv::sourceNuclideDoubleClicked( Wt::WContainerWidget *w )
{
  for( const Nuclide &nuc : m_nuclides )
  {
    if( (nuc.display != w) || !nuc.nuclide )
      continue;
    
    if( m_editingNuclide && (m_editingNuclide.get() != w) )
      m_editingNuclide->removeStyleClass( "EditingNuclide" );
    m_editingNuclide = w;
    w->addStyleClass( "EditingNuclide" );
    
    m_nuclideSelectDialog->show();
    m_nuclideSelect->setNuclideSearchToFocus();
    m_nuclideSelect->setCurrentInfo( nuc.nuclide->massNumber, nuc.nuclide->atomicNumber,
                                     nuc.nuclide->isomerNumber, nuc.age, nuc.activity,
                                     nuc.useCurie, nuc.activityStr );
    m_nuclideSelect->setAddButtonToAccept();
    m_nuclideSelectDialog->centerWindow();
    m_nuclideSelectDialog->resizeToFitOnScreen();
    break;
  }//for( const Nuclide &nuc : m_nuclides )
}//void sourceNuclideDoubleClicked( Wt::WContainerWidget *w );


void DecayActivityDiv::addTheNuclide( const NuclideSelectedInfo &n )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const nuc = db ? db->nuclide( n.z, n.a, n.metasable ) : nullptr;
  if( !nuc )
    return;
  
  int replaceIndex = -1;
  for( size_t i = 0; m_editingNuclide && (i < m_nuclides.size()); ++i )
  {
    if( m_nuclides[i].display == m_editingNuclide.get() )
      replaceIndex = static_cast<int>( i );
  }
  
  // Editing a nuclides activity or age keeps the users time span; a different nuclide does not.
  const bool sameNuclide = (replaceIndex >= 0) && (m_nuclides[replaceIndex].nuclide == nuc);
  
  clearEditingNuclide();
  
  DecayActivityCalc::SourceNuclide src;
  src.nuclide = nuc;
  src.activity = n.activity;
  src.activityStr = n.activityStr;
  src.useCurie = n.useCurie;
  src.age = n.initialAge;
  addNuclideNoRefresh( src, replaceIndex );
  
  m_decayChainChart->setNuclide( nuc, n.useCurie, DecayChainChart::DecayChainType::DecayFrom );
  
  if( sameNuclide )
  {
    manageActiveDecayChainNucStyling();
    if( m_calc )
      m_calc->setInputsStale();
    refreshDecayDisplay( true );
  }else
  {
    nuclidesChanged();
  }
}//void addTheNuclide( const NuclideSelectedInfo &nuc )


void DecayActivityDiv::Nuclide::updateTxt()
{
  if( !txt || !nuclide )
    return;
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Element *element = db ? db->element( nuclide->atomicNumber ) : nullptr;
  const string elementName = (element ? element->symbol : string(""));
  
  stringstream label;
  label << "<sup><font size=\"1.5\">" << nuclide->massNumber;
  
  if( nuclide->isomerNumber )
    label << "m";
  
  if( nuclide->isomerNumber > 1 )
    label << nuclide->isomerNumber;
  
  label << "</font></sup>" << elementName;
  label << " " << PhysicalUnits::printToBestActivityUnits(activity, 2, useCurie );
  
  if( age > 0.0 )
    label << " " << PhysicalUnitsLocalized::printToBestTimeUnits(age);
  
  txt->setText( WString::fromUTF8(label.str()) );
}//void DecayActivityDiv::Nuclide::updateTxt()


void DecayActivityDiv::addNuclide( const int z, const int a, const int iso,
                                   const double activity, const bool useCurie,
                                   const double age, const std::string &activityStr )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const nuc = db ? db->nuclide( z, a, iso ) : nullptr;
  if( !nuc )
  {
    cerr << "DecayActivityDiv::addNuclide: invalid nuclide z=" << z << ", a=" << a
         << ", iso=" << iso << endl;
    return;
  }
  
  DecayActivityCalc::SourceNuclide src;
  src.nuclide = nuc;
  src.activity = activity;
  src.activityStr = activityStr;
  src.useCurie = useCurie;
  src.age = age;
  addNuclideNoRefresh( src );
  
  m_decayChainChart->setNuclide( nuc, useCurie, DecayChainChart::DecayChainType::DecayFrom );
  
  nuclidesChanged();
}//void addNuclide(..)


void DecayActivityDiv::addNuclideNoRefresh( DecayActivityCalc::SourceNuclide src, const int replaceIndex )
{
  assert( src.nuclide );
  if( !src.nuclide )
    return;
  
  DecayActivityCalc::normalizeActivityString( src );
  
  Nuclide nuclide;
  static_cast<DecayActivityCalc::SourceNuclide &>( nuclide ) = src;
  
  auto displayOwner = make_unique<WContainerWidget>();
  WContainerWidget * const disp = displayOwner.get();
  nuclide.display = disp;
  disp->setInline( true );
  disp->addStyleClass( "Nuclide" );
  nuclide.txt = disp->addNew<WText>( "", Wt::TextFormat::XHTML );
  nuclide.updateTxt();
  
  // A single click shows this nuclides decay chain (looked up when clicked, as its units may change).
  disp->clicked().connect( this, [this, disp](){
    for( const Nuclide &n : m_nuclides )
    {
      if( n.display == disp )
        m_decayChainChart->setNuclide( n.nuclide, n.useCurie, DecayChainChart::DecayChainType::DecayFrom );
    }
  } );
  disp->doubleClicked().connect( this, [this, disp](){ sourceNuclideDoubleClicked( disp ); } );
  disp->setToolTip( WString::tr("dad-tt-double-click") );
  
  WPushButton *editIcon = disp->addNew<WPushButton>();
  editIcon->addStyleClass( "NucEditIcon" );
  editIcon->setIcon( "InterSpec_resources/images/edit_pencil.svg" );
  editIcon->setToolTip( WString::tr("dad-tt-edit-nuclide") );
  editIcon->clicked().preventPropagation();
  editIcon->clicked().connect( this, [this, disp](){ sourceNuclideDoubleClicked( disp ); } );
  
  WPushButton *closeIcon = disp->addNew<WPushButton>();
  closeIcon->addStyleClass( "mycloseicon" );
  closeIcon->setToolTip( WString::tr("dad-tt-remove-nuclide") );
  closeIcon->clicked().preventPropagation();
  closeIcon->clicked().connect( this, [this, disp](){ removeNuclide( disp ); } );
  
  if( (replaceIndex >= 0) && (replaceIndex < static_cast<int>(m_nuclides.size())) )
  {
    // Put the new pill where the old one was, and keep the nuclide order (which nuclide is first
    //  matters for "half-life" time spans, and the URL).
    WContainerWidget * const old = m_nuclides[replaceIndex].display;
    m_nuclidesAddedDiv->insertBefore( std::move(displayOwner), old );
    { std::unique_ptr<WWidget> removed = m_nuclidesAddedDiv->removeWidget( old ); }
    m_nuclides[replaceIndex] = nuclide;
  }else
  {
    m_nuclidesAddedDiv->addWidget( std::move(displayOwner) );
    m_nuclides.push_back( nuclide );
  }
}//void addNuclideNoRefresh(...)


void DecayActivityDiv::nuclidesChanged()
{
  updateNuclideAreaState();
  manageActiveDecayChainNucStyling(); // New pills need the highlight, if on the Decay Chain tab
  
  if( m_calc )
    m_calc->setInputsStale();
  
  setTimeLimitToDisplay();
  refreshDecayDisplay( true );
}//void nuclidesChanged()


void DecayActivityDiv::removeNuclide( Wt::WContainerWidget *frame )
{
  const vector<Nuclide>::iterator pos = std::find_if( begin(m_nuclides), end(m_nuclides),
                                          [frame]( const Nuclide &n ){ return n.display == frame; } );
  if( pos == end(m_nuclides) ) //should never happen
  {
    cerr << "removeNuclide(...): Could not find the desired nuclide :(" << endl;
    return;
  }//if( pos == end(m_nuclides) )
  
  if( m_editingNuclide.get() == frame )
    clearEditingNuclide();
  
  const SandiaDecay::Nuclide * const removedNuc = pos->nuclide;
  { std::unique_ptr<WWidget> removed = m_nuclidesAddedDiv->removeWidget( pos->display ); }
  m_nuclides.erase( pos );
  
  // If the decay chain was showing this nuclide, show the last remaining one instead.
  const bool stillHaveNuc = std::any_of( begin(m_nuclides), end(m_nuclides),
                                  [removedNuc]( const Nuclide &n ){ return n.nuclide == removedNuc; } );
  if( !stillHaveNuc && (m_decayChainChart->nuclide() == removedNuc) )
  {
    const Nuclide * const last = m_nuclides.empty() ? nullptr : &m_nuclides.back();
    m_decayChainChart->setNuclide( (last ? last->nuclide : nullptr), (last ? last->useCurie : true),
                                   DecayChainChart::DecayChainType::DecayFrom );
  }
  
  nuclidesChanged();
}//void DecayActivityDiv::removeNuclide( Wt::WContainerWidget *frame )


void DecayActivityDiv::confirmClearAllNuclides()
{
  // Losing a single nuclide is cheap to redo; losing a whole set is not, and nothing in this tool
  //  is undoable, so anything more than one gets a confirmation.
  if( m_nuclides.size() < 2 )
  {
    clearAllNuclides();
    return;
  }

  SimpleDialog *dialog = SimpleDialog::make( WString::tr("dad-remove-all"),
                                             WString::tr("dad-remove-all-confirm") );
  WPushButton *yes = dialog->addButton( WString::tr("Yes"), WidgetUtils::ButtonRole::Affirm );
  yes->clicked().connect( this, &DecayActivityDiv::clearAllNuclides );
  dialog->addButton( WString::tr("No"), WidgetUtils::ButtonRole::Dismiss )->setFocus();
}//void DecayActivityDiv::confirmClearAllNuclides()


void DecayActivityDiv::updateNuclideAreaState()
{
  const bool haveNuclides = !m_nuclides.empty();
  if( m_noNuclidesTxt )
    m_noNuclidesTxt->setHidden( haveNuclides );
  if( m_clearNuclidesButton )
    m_clearNuclidesButton->setDisabled( !haveNuclides );
}//void DecayActivityDiv::updateNuclideAreaState()


void DecayActivityDiv::clearNuclidesNoRefresh()
{
  clearEditingNuclide();
  
  for( const Nuclide &nuclide : m_nuclides )
  {
    std::unique_ptr<WWidget> removed = m_nuclidesAddedDiv->removeWidget( nuclide.display );
  }
  m_nuclides.clear();
}//void clearNuclidesNoRefresh()


void DecayActivityDiv::clearAllNuclides()
{
  clearNuclidesNoRefresh();
  
  m_decayChainChart->setNuclide( nullptr, true, DecayChainChart::DecayChainType::DecayFrom );

  deleteMoreInfoDialog();
  m_decayChainChart->deleteMoreInfoDialog();

  nuclidesChanged();
}//void clearAllNuclides();


void DecayActivityDiv::updateInitialMixture() const
{
  m_currentMixture->clear();

  for( const Nuclide &nuc : m_nuclides )
  {
    if( !nuc.nuclide )
      continue;
    
    if( nuc.age > DBL_EPSILON )
      m_currentMixture->addAgedNuclideByActivity( nuc.nuclide, nuc.activity, nuc.age );
    else
      m_currentMixture->addNuclideByActivity( nuc.nuclide, nuc.activity );
  }//for( each m_nuclides )
}//void updateInitialMixture()


void DecayActivityDiv::applyYScale()
{
  const Chart::AxisScale scale = m_logYScale->isChecked() ? Chart::AxisScale::Log : Chart::AxisScale::Linear;
  m_decayChart->axis(Chart::Axis::Y).setScale( scale );
}//void applyYScale()


void DecayActivityDiv::updateYScale()
{
  applyYScale();
  
  // If we dont do a full refresh, for some reason sometimes the lines wont show up
  refreshDecayDisplay( true );
}//void updateYScale()


void DecayActivityDiv::displayMoreInfoPopup( const double time )
{
  WContainerWidget *summary = isotopesSummary( time );

  if( m_viewer && (m_viewer->renderedHeight() > 100) )
    summary->setMaximumSize( WLength::Auto, WLength( 0.8*m_viewer->renderedHeight(), WLength::Unit::Pixel) );
  
  WString title = WString::tr("dad-summary-at-time").arg( PhysicalUnitsLocalized::printToBestTimeUnits(time) );

  if( m_moreInfoDialog )
  {
    m_moreInfoDialog->contents()->clear();
    m_moreInfoDialog->setWindowTitle( title );
  }else
  {
    m_moreInfoDialog = AuxWindow::make( title );
    m_moreInfoDialog->setClosable( true );
    m_moreInfoDialog->disableCollapse();
    m_moreInfoDialog->rejectWhenEscapePressed();
    m_moreInfoDialog->finished().connect( this, [this](){ deleteMoreInfoDialog(); } );

    WPushButton *ok = m_moreInfoDialog->addCloseButtonToFooter();
    ok->clicked().connect( this, [this](){ deleteMoreInfoDialog(); } );
  }//if( m_moreInfoDialog ) / else

  m_moreInfoDialog->contents()->addWidget( unique_ptr<WWidget>(summary) );
  
  m_moreInfoDialog->setHidden( false );
  m_moreInfoDialog->resizeToFitOnScreen();
  m_moreInfoDialog->centerWindow();
}//void displayMoreInfoPopup( const double time )


void DecayActivityDiv::checkTimeRangeValid()
{
  double timelen = 0.0;
  
  try
  {
    timelen = parseTimeSpan( m_displayTimeLength->text().toUTF8() );
  }catch( std::exception & )
  {
    timelen = DecayActivityCalc::findTimeForActivityFrac( *m_currentMixture, 0.1 );
    m_displayTimeLength->setText( WString::fromUTF8( PhysicalUnitsLocalized::printToBestTimeUnits(timelen) ) );
  }
  
  double minhl = DBL_MAX;
  for( const Nuclide &nuc : m_nuclides )
  {
    if( nuc.nuclide )
      minhl = std::min( minhl, nuc.nuclide->halfLife );
  }//for( const Nuclide &nuc : m_nuclides )
  
  minhl = std::max( minhl, 1.0E-6*SandiaDecay::second );
  minhl = std::min( minhl, 1000.0*SandiaDecay::second );
  if( timelen < 0.001*minhl )
  {
    timelen = 0.001*minhl;
    m_displayTimeLength->setText( WString::fromUTF8( PhysicalUnitsLocalized::printToBestTimeUnits(timelen) ) );
  }//if( nunits*unitpair.second < 0.001*minhl )
}//void DecayActivityDiv::checkTimeRangeValid()


void DecayActivityDiv::setTimeLimitToDisplay()
{
  updateInitialMixture();

  const double finalTime = DecayActivityCalc::findTimeForActivityFrac( *m_currentMixture, 0.1 );
  m_displayTimeLength->setText( WString::fromUTF8( PhysicalUnitsLocalized::printToBestTimeUnits(finalTime) ) );
}//void DecayActivityDiv::setTimeLimitToDisplay()


void DecayActivityDiv::setDecayChartTimeRange( double finalTime )
{
  m_displayTimeLength->setText( WString::fromUTF8( PhysicalUnitsLocalized::printToBestTimeUnits(finalTime) ) );

  // Redraw with the new range, as editing the field does (addNuclide draws with its own range).
  refreshDecayDisplay( true );
}//void setDecayChartTimeRange()


void DecayActivityDiv::colorThemeChanged( std::shared_ptr<const ColorTheme> theme )
{
  m_decayChainChart->colorThemeChanged();
  
  // The series pens come from the theme, so re-add them too.
  m_decayChart->setColorTheme( theme );
  addDecaySeries();
}//void colorThemeChanged( std::shared_ptr<const ColorTheme> theme )


void DecayActivityDiv::deleteMoreInfoDialog()
{
  if( m_moreInfoDialog )
    AuxWindow::deleteAuxWindow( m_moreInfoDialog.get() );
  assert( !m_moreInfoDialog );
}//void deleteMoreInfoDialog()


void DecayActivityDiv::manageActiveDecayChainNucStyling()
{
  for( const Nuclide &nuc : m_nuclides )
  {
    if( nuc.display && nuc.display->hasStyleClass( "ActiveDecayChainNuc" ) )
      nuc.display->removeStyleClass( "ActiveDecayChainNuc" );
  }//for( const Nuclide &nuc : m_nuclides )
  
  const SandiaDecay::Nuclide * const nuc = m_decayChainChart->nuclide();
  const int chainTab = static_cast<int>( DecayActivityCalc::Tab::Chain );
  if( !nuc || (m_chartTabWidget->currentIndex() != chainTab) )
    return;
  
  for( const Nuclide &n : m_nuclides )
  {
    if( (n.nuclide == nuc) && n.display )
      n.display->addStyleClass( "ActiveDecayChainNuc" );
  }//for( const Nuclide &n : m_nuclides )
}//void manageActiveDecayChainNucStyling()


void DecayActivityDiv::decayChartClicked( const WMouseEvent& event )
{
  if( m_currentTimeUnits <= 0.0 )
    return;

  // The chart is a scatter chart, with the x-axis in units of m_currentTimeUnits
  const Coordinates coords = event.widget();
  const WPointF cartesianPoint = m_decayChart->mapFromDevice( WPointF(coords.x, coords.y) );
  
  displayMoreInfoPopup( cartesianPoint.x() * m_currentTimeUnits );
}//void DecayActivityDiv::decayChartClicked( const WMouseEvent& event )


WContainerWidget *DecayActivityDiv::isotopesSummary( const double time ) const
{
  using namespace SandiaDecay;
  
  WContainerWidget *cont = new WContainerWidget();
  cont->addStyleClass( "isotopesSummary" );
  cont->addStyleClass( "SelectableTxt" );

  if( !m_currentMixture || !m_currentMixture->numInitialNuclides() )
    return cont;
  
  const vector<NuclideActivityPair> activities = m_currentMixture->activity( time );
  
  stringstream infostrm;
  infostrm << "<div class=\"NuclideInfo\">\n"
           << "<div>" << WString::tr("dad-summary-starting-from").toUTF8() << "</div>"
           << "<div style=\"margin-left: 20px; max-width: 60ex;\">";
  
  for( int i = 0; i < m_currentMixture->numInitialNuclides(); ++i )
  {
    const SandiaDecay::Nuclide *initial_nuc = m_currentMixture->initialNuclide(i);
    const double initial_act = m_currentMixture->initialActivity(i);
    infostrm << (i ? ", " : "")
             << PhysicalUnits::printToBestActivityUnits(initial_act, 2, useCurieFor(initial_nuc), SandiaDecay::becquerel)
             << " " << initial_nuc->symbol;
  }//for( loop over orignal nuclides )
  
  infostrm << "</div>";
  
  infostrm << "<div style=\"margin-top: 10px;\">"
           << WString::tr("dad-summary-nucs-present")
                .arg( PhysicalUnitsLocalized::printToBestTimeUnits(time/SandiaDecay::second,4) ).toUTF8()
           << "</div>\n"
           << "<div style=\"margin-left: 20px; max-width: 60ex;\">";
  
  int nparents = 0;
  for( const NuclideActivityPair &pair : activities )
  {
    if( IsInf(pair.nuclide->halfLife) || IsNan(pair.nuclide->halfLife) )
      continue;
    
    const bool useCurie = useCurieFor( pair.nuclide );
    infostrm << (nparents ? ", " : "")
             << PhysicalUnits::printToBestActivityUnits(pair.activity, 2, useCurie, SandiaDecay::becquerel)
             << " " << pair.nuclide->symbol;
    ++nparents;
  }//for( loop over activities )
  
  infostrm << "</div>\n";
  
  // Lines with a zero rate (e.g., from a progeny not yet built in) are left out.
  const string perSecond = WString::tr("dad-summary-per-second").toUTF8();
  auto addParticles = [&infostrm, &perSecond]( const char *headingKey, const vector<EnergyRatePair> &particles ){
    if( std::none_of( begin(particles), end(particles), []( const EnergyRatePair &p ){ return p.numPerSecond > 0.0; } ) )
      return;
    
    infostrm << "<div style=\"margin-top: 10px;\">" << WString::tr(headingKey).toUTF8() << "</div>\n"
             << "<table style=\"margin-left: 20px;\">\n";
    for( const EnergyRatePair &p : particles )
    {
      if( p.numPerSecond > 0.0 )
        infostrm << "<tr><td>" << p.energy/SandiaDecay::keV << " keV</td><td>&nbsp;-&nbsp;</td><td>"
                 << p.numPerSecond << " " << perSecond << "</td></tr>\n";
    }
    infostrm << "</table>\n";
  };//addParticles
  
  addParticles( "dad-summary-gammas", m_currentMixture->gammas( time, NuclideMixture::OrderByAbundance, true ) );
  addParticles( "dad-summary-alphas", m_currentMixture->alphas( time ) );
  addParticles( "dad-summary-betas", m_currentMixture->betas( time ) );
  addParticles( "dad-summary-beta-plus", m_currentMixture->betaPlusses( time ) );
  
  infostrm << "\n</div>";
  cont->addNew<WText>( WString::fromUTF8(infostrm.str()), TextFormat::XHTML );

  cont->addNew<WText>( "<br /><br /><H3>" + WString::tr("dad-summary-further-info").toUTF8() + "</H3>",
                       TextFormat::XHTML );

  const int nSolutionNuclides = m_currentMixture->numSolutionNuclides();
  for( int i = 0; i < nSolutionNuclides; ++i )
  {
    const SandiaDecay::Nuclide *nuclide = m_currentMixture->solutionNuclide(i);
    cont->addWidget( unique_ptr<WContainerWidget>( nuclideInformation(nuclide) ) );
  }//for( loop over widgets )

  return cont;
}//WContainerWidget *isotopesSummary( const double time ) const


WContainerWidget *DecayActivityDiv::nuclideInformation( const SandiaDecay::Nuclide *nuclide ) const
{
  WContainerWidget *cont = new WContainerWidget();
  cont->addStyleClass( "nuclideInformation" );

  if( !nuclide )
    return cont;
  
  // A dump of the SandiaDecay data, so field names are left as-is (not translated).
  stringstream ostr;
  ostr << "<pre class=\"NuclideInfo\">\n";
  ostr << nuclide->symbol << " Atomic Number " << nuclide->atomicNumber <<", Atomic Mass "
  << nuclide->massNumber << ", Isomer Number " << nuclide->isomerNumber << " "
  << nuclide->atomicMass << " AMU, HalfLife=" << PhysicalUnitsLocalized::printToBestTimeUnits(nuclide->halfLife);
  const size_t nParents = nuclide->decaysFromParents.size();
  if( nParents )
    ostr << "\n  Parent" << ((nParents > 1) ? "s" : "") << ": ";
  
  for( size_t i = 0; i < nParents; ++i )
  {
    if( i )
      ostr << ", ";
    if( nuclide->decaysFromParents[i] && nuclide->decaysFromParents[i]->parent )
      ostr << nuclide->decaysFromParents[i]->parent->symbol;
  }//for( loop over parents )
  
  for( const SandiaDecay::Transition * const child : nuclide->decaysToChildren )
  {
    ostr << "\nTransition ";
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
    for( const SandiaDecay::RadParticle &product : child->products )
    {
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
    }//for( const SandiaDecay::RadParticle &product : child->products )
  }//for( loop over decays to children )
  
  ostr << "\n</pre>\n";
  cont->addNew<WText>( ostr.str(), Wt::TextFormat::UnsafeXHTML );

  return cont;
}//WContainerWidget *nuclideInformation( const Nuclide &nuclide ) const;


void DecayActivityDiv::refreshDecayDisplay( const bool syncCalcTimeSpan )
{
  refreshDecayChart();
  
  if( m_calc && syncCalcTimeSpan )
    m_calc->setTimeRangeTxt( m_displayTimeLength->text().toUTF8() );
  
  updateCalculator();
}//void refreshDecayDisplay( const bool syncCalcTimeSpan )


void DecayActivityDiv::updateCalculator()
{
  if( !m_calc )
    return;
  
  // The calculator results can be a lot of widgets, so only build them when they are seen.
  const int calcTab = static_cast<int>( DecayActivityCalc::Tab::Calc );
  m_calcNeedsUpdate = (m_chartTabWidget->currentIndex() != calcTab);
  if( !m_calcNeedsUpdate )
    m_calc->update();
}//void updateCalculator()


void DecayActivityDiv::chartTabChanged( int index )
{
  if( m_calcNeedsUpdate && (index == static_cast<int>(DecayActivityCalc::Tab::Calc)) )
    updateCalculator();
}//void chartTabChanged( int index )


void DecayActivityDiv::refreshDecayChart()
{
  checkTimeRangeValid();
  
  // Keep whether the user chose to show the total, unless the previous chart had a single nuclide,
  //  which always hides its total (it is the same as the nuclide).
  bool showsum = true;
  const int noldcol = m_decayModel->columnCount();
  if( noldcol > 3 )
    showsum = m_decayModel->showSeries( noldcol - 1 );
  
  //Grab if we should show the nuclides before clearing the model.
  //  Note that the very first and very last column no not contain nuclide info
  map<string,bool> oldShowNuclide;
  for( int column = 0; column < m_decayModel->columnCount(); ++column )
  {
    const string nuc = m_decayModel->nuclide(column);
    if( nuc != "" )
      oldShowNuclide[nuc] = m_decayModel->showSeries( column );
  }
  
  //Remove the series from last time, and reset the model
  {
    std::vector<Wt::Chart::WDataSeries *> seriesToRemove;
    for( const auto &s : m_decayChart->series() )
      seriesToRemove.push_back( s.get() );
    for( Wt::Chart::WDataSeries *s : seriesToRemove )
      m_decayChart->removeSeries( s );
  }
  
  m_decayModel->clear();
  m_currentTimeUnits  = -1.0;

  updateInitialMixture();
  
  if( m_nuclides.empty() )
  {
    m_decayLegend->clear();
    return;
  }//if( m_nuclides.empty() )

  const int activityIndex = m_displayActivityUnitsCombo->currentIndex();

  double actunit;
  string unitStr;

  if( (activityIndex >= 0) && (activityIndex < static_cast<int>( sm_activityUnitNameValues.size() )) )
  {
    const UnitNameValuePair &act_name_unit = sm_activityUnitNameValues.at( activityIndex );
    unitStr = act_name_unit.first;
    actunit = act_name_unit.second;
  }else
  {
    unitStr = "I0";
    actunit = m_currentMixture->totalActivity(0.0);
  }//if( user selected a standard unit ) / else

  const double maxDiplayTime = timeToDisplayTill();
  
  const int yaxisIndex = m_yAxisType->currentIndex();
  const DecayActivityCalc::ChartType yaxis
          = ((yaxisIndex >= 0) && (yaxisIndex <= static_cast<int>(DecayActivityCalc::ChartType::Alpha)))
            ? static_cast<DecayActivityCalc::ChartType>( yaxisIndex )
            : DecayActivityCalc::ChartType::Activity;
  
  const bool showActUnits = (yaxis == DecayActivityCalc::ChartType::Activity);
  m_displayActivityUnitsCombo->setHidden( !showActUnits );
  m_displayActivityUnitsLabel->setHidden( !showActUnits );
  
  if( maxDiplayTime <= 0.0 )
    return;
  
  SandiaDecay::ProductType particletype = SandiaDecay::GammaParticle;
  switch( yaxis )
  {
    case DecayActivityCalc::ChartType::Activity:                                             break;
    case DecayActivityCalc::ChartType::Gamma: particletype = SandiaDecay::GammaParticle; break;
    case DecayActivityCalc::ChartType::Beta:  particletype = SandiaDecay::BetaParticle;  break;
    case DecayActivityCalc::ChartType::Alpha: particletype = SandiaDecay::AlphaParticle; break;
  }//switch( yaxis )
  
  // Two pixels per point, within limits.
  const int nRows = std::clamp( m_decayChart->paintedWidth() / 2, sm_minChartPoints, sm_maxChartPoints );
  
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

  vector<double> totalActivities( nRows, 0.0 );

  const UnitNameValuePair xUnitsPair = bestTimeUnit( maxDiplayTime );
  const double tunit = xUnitsPair.second;
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
    m_decayModel->setData( row, 0, Wt::cpp17::any( label ) );

    for( int elN = 0; elN < nElements; ++elN )
    {
      const int column = elN + 1;
      const double activity = evolutions[elN].activity(row_time);
      const double yval = showActUnits
                    ? (activity / actunit)
                    : DecayActivityCalc::particleRate( *evolutions[elN].nuclide, activity, particletype );

      if( IsInf(yval) || IsNan(yval) )
        continue;

      totalActivities[row] += yval;
      m_decayModel->setData( row, column, Wt::cpp17::any( yval ) );
    }//for( loop over nuclides to add )
  }//for( loop over time points to add )

  //Now put in the sum of all the activities
  const WString dateHeader = WString::tr("dad-hdr-date").arg(xUnitsPair.first);
  m_decayModel->setHeaderData( 0, Wt::cpp17::any( dateHeader ) );

  for( int column = 1; column <= nElements; ++column )
  {
    const string &symbol = evolutions[column-1].nuclide->symbol;
    m_decayModel->setNuclide( column, symbol );
    m_decayModel->setHeaderData( column, Wt::cpp17::any( WString::fromUTF8(symbol) ) );
    m_decayModel->setHeaderData( column, Wt::Orientation::Horizontal, true, Wt::ItemDataRole::User );
  }//for( int column = 0; column < nElements; ++column )

  for( int row = 0; row < nRows; ++row )
    m_decayModel->setData( row, nElements+1, Wt::cpp17::any(totalActivities[row]) );
  
  switch( yaxis )
  {
    case DecayActivityCalc::ChartType::Activity:
      m_decayModel->setHeaderData( nElements+1, Wt::cpp17::any( WString::tr("dad-hdr-total-act") ) );
      m_decayChart->axis(Chart::Axis::Y).setTitle( unitStr );
    break;
    case DecayActivityCalc::ChartType::Gamma:
      m_decayModel->setHeaderData( nElements+1, Wt::cpp17::any( WString::tr("dad-hdr-total-gammas") ) );
      m_decayChart->axis(Chart::Axis::Y).setTitle( WString::tr("dad-yaxis-title-total-gammas") );
    break;
    case DecayActivityCalc::ChartType::Beta:
      m_decayModel->setHeaderData( nElements+1, Wt::cpp17::any( WString::tr("dad-hdr-total-betas") ) );
      m_decayChart->axis(Chart::Axis::Y).setTitle( WString::tr("dad-yaxis-title-total-betas") );
    break;
    case DecayActivityCalc::ChartType::Alpha:
      m_decayModel->setHeaderData( nElements+1, Wt::cpp17::any( WString::tr("dad-hdr-total-alphas") ) );
      m_decayChart->axis(Chart::Axis::Y).setTitle( WString::tr("dad-yaxis-title-total-alphas") );
    break;
  }//switch( yaxis )
  
  m_decayModel->setHeaderData( nElements+1, Wt::Orientation::Horizontal,
                                 Wt::cpp17::any(showsum), Wt::ItemDataRole::User );
  
  //Now check if we have the same nuclides as before, and if so, set the
  //  visability properties for each of the nuclides charts to be the same
  //  as before.
  bool sameNuclides = (static_cast<size_t>(nElements) == oldShowNuclide.size());
  for( int column = 1; sameNuclides && (column <= nElements); ++column )
    sameNuclides = (oldShowNuclide.count( evolutions[column-1].nuclide->symbol ) > 0);
  
  if( sameNuclides )
  {
    for( int column = 0; column < m_decayModel->columnCount(); ++column )
    {
      const string nuc = m_decayModel->nuclide(column);
      if( nuc != "" )
        m_decayModel->setShowSeries( column, oldShowNuclide[nuc] );
    }
  }//if( sameNuclides )
  
  addDecaySeries();
  
  m_currentTimeUnits = tunit;
  
  m_decayModel->updateMaxActivity();
  
  updateYAxisRange();
}//void refreshDecayChart()


void DecayActivityDiv::updateYAxisRange()
{
  //The next line of letting Wt set y-axis limits doesnt always work super great, so we'll do it
  //  manually
  // m_decayChart->axis(Chart::Axis::Y).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
  
  double miny = std::numeric_limits<double>::max(), maxy = 0.0;
  
  for( int row = 0; row < m_decayModel->rowCount(); ++row )
  {
    for( int col = 1; col < m_decayModel->columnCount(); ++col )
    {
      if( !m_decayModel->showSeries( col ) )
        continue;
      
      Wt::cpp17::any data = m_decayModel->data( m_decayModel->index(row, col) );
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
  switch( m_decayChart->axis(Chart::Axis::Y).scale() )
  {
    case Wt::Chart::AxisScale::Log:
      if( miny < 0.000001*maxy )
        min_disp_act = 0.00001*maxy;
      else
        min_disp_act = 0.75*miny;
      max_disp_act = 1.5*maxy;
      break;

    case Wt::Chart::AxisScale::Linear:
      //If it wont change the dynamic range of the chart much anyways, might as well
      //  anchor the y-axis to zero
      break;

    case Wt::Chart::AxisScale::Discrete:
    case Wt::Chart::AxisScale::Date:
    case Wt::Chart::AxisScale::DateTime:
      assert( 0 );
      break;
  }//switch( m_decayChart->axis(Chart::Axis::Y).scale() )
  
  m_decayChart->axis(Chart::Axis::Y).setRange( min_disp_act, max_disp_act );
}//void updateYAxisRange();


void DecayActivityDiv::userSetShowSeries( int series, bool show )
{
  m_decayModel->setShowSeries( series, show );
  updateYAxisRange();
}//void userSetShowSeries()


void DecayActivityDiv::addDecaySeries()
{
  {
    std::vector<Wt::Chart::WDataSeries *> seriesToRemove;
    for( const auto &s : m_decayChart->series() )
      seriesToRemove.push_back( s.get() );
    for( Wt::Chart::WDataSeries *s : seriesToRemove )
      m_decayChart->removeSeries( s );
  }
  
  const int nColumns = m_decayModel->columnCount();
  
  std::shared_ptr<const ColorTheme> theme = m_decayChart->colorTheme();
  if( theme && theme->theme_name.toUTF8().find("Default")!=string::npos )
    theme.reset();
  
  for( int column = 1; column < nColumns; ++column )
  {
    auto seriesOwned = make_unique<Chart::WDataSeries>( column, Chart::SeriesType::Curve, Chart::Axis::Y );
    WPen pen;

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
      pen.setStyle( Wt::PenStyle::DashLine );
    }

    seriesOwned->setPen( pen );

    m_decayChart->addSeries( std::move(seriesOwned) );
  }//for( int column = 1; column < nElements; ++column )
  
  m_decayLegend->clear();
  
  for( int column = nColumns-1; column > 0; --column )
  {
    unique_ptr<WWidget> itemOwned = m_decayChart->createLegendItemWidget( column );
    WContainerWidget *ww = dynamic_cast<WContainerWidget *>( itemOwned.get() );

    if( nColumns == 3 )  //there is exactly one nuclide
    {
      if( column == (nColumns-1) )  //we're on the Total Activity Entry
      {
        // itemOwned goes out of scope and is deleted
        m_decayModel->setShowSeries( column, false );
        continue;
      }
      m_decayModel->setShowSeries( column, true );
    }//if( nColumns == 3 )

    if( ww )
    {
      ww->setInline( false );
      WCheckBox *cb = ww->insertWidget( 0, make_unique<WCheckBox>() );
      ww->addStyleClass( "LegendEntry" );
      m_decayLegend->addWidget( std::move(itemOwned) );
      const bool check = m_decayModel->showSeries( column );
      cb->setChecked( check );
      cb->addStyleClass( "ShowSeriesCb" );
      cb->setDisabled( (nColumns == 3) );
      cb->checked().connect( this, [this, column](){ userSetShowSeries( column, true ); } );
      cb->unChecked().connect( this, [this, column](){ userSetShowSeries( column, false ); } );
    }else
    {
      // itemOwned goes out of scope and is deleted
      cerr << "DecayActivityDiv::addDecaySeries(): legend items are no longer"
           << " WContainerWidget, please update the code"<< endl;
    }//if( ww )
  }
}//void addDecaySeries()


DecayActivityDiv::~DecayActivityDiv()
{
  if( m_nuclideSelectDialog )
    AuxWindow::deleteAuxWindow( m_nuclideSelectDialog.get() );
  
  deleteCsvDownloadGui();

  // Owned by wApp (AuxWindow::make), so it does not go away with us; its Close button is a tracked
  //  connection to this widget, so once we are gone the window is up and un-closeable.
  deleteMoreInfoDialog();
}//DecayActivityDiv destructor
