/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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
#include <Wt/WRectArea>
#include <Wt/WTemplate>
#include <Wt/WCheckBox>
#include <Wt/WGroupBox>
#include <Wt/WCalendar>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WResource>
#include <Wt/WTabWidget>
#include <Wt/WPopupMenu>
#include <Wt/WGridLayout>
#include <Wt/WDatePicker>
#include <Wt/WPushButton>
#include <Wt/WPaintDevice>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WDoubleSpinBox>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WStandardItemModel>
#include <Wt/Chart/WCartesianChart>

#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "sandia_decay/SandiaDecay.h"
#include "InterSpec/DecayActivityDiv.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DecayChainChart.h"
#if( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK )
#include "InterSpec/ChartToImageResource.h"
#endif
#include "InterSpec/DecaySelectNuclideDiv.h"
#include "InterSpec/MassAttenuationTool.h"


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


}//namespace

class SrbActivityModel : public Wt::WStandardItemModel
{
  
public:
  SrbActivityModel( Wt::WObject *parent = 0 )
    : Wt::WStandardItemModel( parent )
  {
  }
  
  virtual ~SrbActivityModel()
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
      cerr << "SrbActivityModel::data(...) unexpectedly caught: " << e.what()
           << endl << index.row() << ", " << index.column()
           << endl;
    }//try / catch
    
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

  
  void setShowSeries( int colum, bool show )
  {
    const bool oldval = showSeries( colum );
    if( oldval == show )
      return;
    
    if( colum < 1 || colum >= columnCount() )
      return;
    
    setHeaderData( colum, Wt::Horizontal, boost::any(show), Wt::UserRole );
    
    if( rowCount() )
    {
      WModelIndex left = index( 0, colum );
      WModelIndex right = index( rowCount()-1, colum );
      dataChanged().emit( left, right );
    }//if( rowCount() )
  }//bool setShowSeries( int colum, bool show ) const
};//class SrbActivityModel

class SrbActivityChart : public Wt::Chart::WCartesianChart
{
protected:
  double m_widthInPixels, m_heightInPixels;
  
public:
  SrbActivityChart(  Wt::WContainerWidget *parent = NULL  )
  : Wt::Chart::WCartesianChart( parent ),
    m_widthInPixels(250),
    m_heightInPixels(250)
  {
    setStyleClass( "SrbActivityChart" );
  }//SrbActivityChart( constructor )
  
  virtual ~SrbActivityChart()
  {
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
    painter.setPen( WPen() );
    painter.setBrush( WBrush() );
    const std::string ytitle = axis(Chart::YAxis).title().toUTF8();
    if( ytitle != "" )
    {
      const double height = painter.window().height();
      painter.drawText( -0.45*height, 0.0, 0.2, 0.1,
                       AlignCenter, axis(Chart::YAxis).title() );
    }//if( ytitle != "" )
    
    painter.restore();
  }//paint( ... )
};//class SrbActivityChart




class DateLengthCalculator : public WContainerWidget
{
  DecayActivityDiv *m_activityDiv;
  WGridLayout *m_layout;
  WDateEdit *m_begindate;
  WDateEdit *m_enddate;
  WLineEdit *m_duration;
  WPushButton *m_updateParent;
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
     m_updateParent( NULL ),
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
    m_begindate->setFormat( "MM/dd/yyyy" );
    m_begindate->setPlaceholderText( "mm/dd/yyy" );
    m_begindate->setDate(Wt::WDate::currentServerDate().addDays(-ndaysBack) );
    m_layout->addWidget( m_begindate, 0, 1, 1, 2, AlignLeft );
    
    label = new WLabel( "End Date" );
    m_layout->addWidget( label, 1, 0, AlignMiddle );
    m_enddate = new WDateEdit();
    m_enddate->setFormat( "MM/dd/yyyy" );
    m_enddate->setPlaceholderText( "mm/dd/yyy" );
    m_enddate->setDate( Wt::WDate::currentServerDate() );
    m_layout->addWidget( m_enddate, 1, 1, 1, 2, AlignLeft );
    
    m_begindate->changed().connect( this, &DateLengthCalculator::dateChanged );
    m_enddate->changed().connect( this, &DateLengthCalculator::dateChanged );
    
    label = new WLabel( "Time Span" );
    m_layout->addWidget( label, 2, 0, AlignMiddle );
    
    m_duration = new WLineEdit();
    WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex, m_duration );
    validator->setFlags(Wt::MatchCaseInsensitive);
    m_duration->setValidator(validator);
    m_duration->changed().connect( this, &DateLengthCalculator::durationChanged );
    m_duration->enterPressed().connect( this, &DateLengthCalculator::durationChanged );
    m_layout->addWidget( m_duration, 2, 1 );
    
    m_updateParent = new WPushButton( "Update Chart");
    m_updateParent->clicked().connect( this, &DateLengthCalculator::pushCurrentToParent );
    
    m_layout->addWidget( m_updateParent, 2, 2, AlignCenter );
    
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
  
  
  double getValidatedTimeSpan()
  {
    double halfLife = 0.0;
    const SandiaDecay::Nuclide *nuc = NULL;
    const std::vector<DecayActivityDiv::Nuclide> &nucs = m_activityDiv->m_nuclides;
    if( !nucs.empty() )
      nuc = m_activityDiv->m_nuclideDB->nuclide( nucs[0].z, nucs[0].a, nucs[0].iso );
    if( nuc )
      halfLife = nuc->halfLife;
    
    //Let be optimistic things will be okay
    m_error->hide();
    m_error->setText("");
    m_updateParent->enable();
    
    try
    {
      const std::string durtxt = m_duration->text().narrow();
      const bool validBeginDate = (m_begindate->validate() == Wt::WValidator::Valid);
      const bool validEndDate = (m_enddate->validate() == Wt::WValidator::Valid);
      bool duration_valid = false;
      double duration = 0.0;
      try
      {
        duration = PhysicalUnits::stringToTimeDurationPossibleHalfLife( durtxt, halfLife );
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
        return duration;
      }else if( !duration_valid )
      {
        if( !validBeginDate )
          throw runtime_error( "Begin date is invalid" );
        if( !validEndDate )
          throw runtime_error( "End date is invalid" );
        if( m_begindate->date() >= m_enddate->date() )
          throw runtime_error( "End date must be later than start date" );
        
        const int days = m_begindate->date().daysTo(m_enddate->date());
        duration = 24*3600*days*PhysicalUnits::second;
        const string datestr = PhysicalUnits::printToBestTimeUnits( duration, 2 );
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
      m_updateParent->disable();
      m_info->hide();
      
      return 0.0;
    }//try / catch
    
    return 0.0; //avoid warning
  }//double getValidatedTimeSpan()
  
  
  void pushCurrentToParent()
  {
    if( getValidatedTimeSpan() == 0.0 )
      return;
    
    const string txt = m_duration->text().narrow();
    m_activityDiv->m_displayTimeLength->setText( txt );
    m_activityDiv->refreshDecayDisplay();
    m_activityDiv->updatePhotopeakSliderEndDateText();
    m_activityDiv->m_chartTabWidget->setCurrentIndex( 0 );
  }//void pushCurrentToParent()
  
  
  void dateChanged()
  {
    m_duration->setText( "" );
    
    if( getValidatedTimeSpan() == 0.0 )
      return;
    
    int days = m_begindate->date().daysTo(m_enddate->date());
    days = abs(days);
    string datestr = PhysicalUnits::printToBestTimeUnits( 24*3600*days*PhysicalUnits::second, 2 );
    m_duration->setText( datestr );
    
    updateInfo();
  }//void dateChanged()
  
  
  void durationChanged()
  {
    if( getValidatedTimeSpan() == 0.0 )
      return;
  
    updateInfo();
  }//void durationChanged()
  
  void setTimeRangeTxt( std::string txt )
  {
    if( txt != m_duration->text().narrow() )
    {
      m_duration->setText( txt );
      m_begindate->setText( "" );
    }
    
    updateInfo(); //the
  }//void setTimeRangeTxt( std::string txt )
  
  
  void updateInfo()
  {
    m_info->clear();
    
    const double timeSpan = getValidatedTimeSpan();
    
    if( timeSpan == 0.0 )
    {
      m_info->hide();
      return;
    }
    
    m_info->show();
    
    
    SandiaDecay::NuclideMixture *mix = m_activityDiv->m_currentMixture;
    if( !mix )
      return;
    
    const string txt = m_duration->text().narrow();
    
    const std::vector<DecayActivityDiv::Nuclide> &nucs = m_activityDiv->m_nuclides;
    if( nucs.empty() )
      return;

    for( size_t i = 0; i < nucs.size(); ++i )
    {
      const SandiaDecay::Nuclide *nuc = m_activityDiv->m_nuclideDB->nuclide( nucs[i].z, nucs[i].a, nucs[i].iso );
      if( !nuc )
        continue;
      
      string info = nuc->symbol + " had an initial activity of "
                    + PhysicalUnits::printToBestActivityUnits( nucs[i].activity, 2, nucs[i].useCurrie );
      
      if( m_begindate->validate() == Wt::WValidator::Valid && !m_begindate->text().empty() )
        info += " on " + m_begindate->text().narrow();
      
      if( nucs[i].age > 0 )
        info += " with an initial age of " + PhysicalUnits::printToBestTimeUnits(nucs[i].age);
      
      WText *txt = new WText( info, m_info );
      txt->setInline( false );
    }
    
    {
      WText *line = new WText( "&nbsp;", m_info );
      line->setInline( false );
      line = new WText( "After " + txt + ":", m_info );
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
    
    const vector<SandiaDecay::NuclideActivityPair> activities = mix->activity( timeSpan );
    const vector<SandiaDecay::NuclideNumAtomsPair> numatoms = mix->numAtoms( timeSpan );
    
    for( size_t actnum = 0; actnum < activities.size(); ++actnum )
    {
      const SandiaDecay::NuclideActivityPair &nap = activities[activities.size()-1-actnum];
      if( !nap.nuclide )
        continue;
      
      //Figure out if user inputed activity in Ci or Bq for this nuclide
      bool useCurries = true;
      for( size_t nucn = 0; nucn < nucs.size(); ++nucn )
      {
        const vector<const SandiaDecay::Nuclide *> parents = nap.nuclide->forebearers();
        const SandiaDecay::Nuclide *nuc = m_activityDiv->m_nuclideDB->nuclide( nucs[nucn].z, nucs[nucn].a, nucs[nucn].iso );
        if( std::find(parents.begin(),parents.end(),nuc) != parents.end() )
        {
          useCurries = nucs[nucn].useCurrie;
          break;
        }
      }//for( size_t nucn = 0; nucn < nucs.size(); ++nucn )
      
      infotable->elementAt(1+actnum,0)->addWidget( new WText(nap.nuclide->symbol) );
      
      string actstr = PhysicalUnits::printToBestActivityUnits( nap.activity, 2, useCurries, SandiaDecay::becquerel );
      if( IsInf(nap.nuclide->halfLife) )
        actstr = "stable";
      infotable->elementAt(1+actnum,1)->addWidget( new WText(actstr) );
      
      
      const SandiaDecay::NuclideNumAtomsPair *natomp = NULL;
      for( size_t a = 0; !natomp && a < numatoms.size(); ++a )
        if( numatoms[a].nuclide == nap.nuclide )
          natomp = &(numatoms[a]);
      
      if( natomp )
      {
        const double num_atoms_in_gram = nap.nuclide->atomsPerGram();
        const double ngrams = natomp->numAtoms / num_atoms_in_gram;
        const string massstr = PhysicalUnits::printToBestMassUnits( ngrams*PhysicalUnits::gram );
        infotable->elementAt(1+actnum,2)->addWidget( new WText(massstr) );
      }
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
    ChartAndLegendHolder( SrbActivityChart *chart, WContainerWidget *legend,
                          WContainerWidget *parent = 0 )
    : WContainerWidget( parent ),
      m_chart( chart ),
      m_legend( legend )
    {
      const string chartref = "$('#" + m_chart->id() + "').get(0)";
      const string legref = "$('#" + legend->id() + "')";
      
      addWidget( m_chart );
      addWidget( m_legend );
      
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
      ""  + legref + ".css('max-height',(h-100) +'px');"
      "" "}";
      setJavaScriptMember( "wtResize", js );
    }//ChartAndLegendHolder constructor
    
    
    SrbActivityChart *m_chart;
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
        m_sessionId( wApp->sessionId() )
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
        
        const SandiaDecay::SandiaDecayDataBase * const db = m_display->m_nuclideDB;
      
        string name;
        bool use_curries = false;
        for( const DecayActivityDiv::Nuclide &n : m_display->m_nuclides )
        {
          const SandiaDecay::Nuclide *nuc = db->nuclide( n.z, n.a, n.iso );
          if( nuc && !IsInf(nuc->halfLife) )  //dont print out activity for stable nuclides.
          {
            use_curries = (use_curries || n.useCurrie);
            name += (!name.empty() ? "_" : "") + nuc->symbol;
            if( n.age > std::numeric_limits<float>::epsilon() )
              name += remove_space_copy( PhysicalUnits::printToBestTimeUnits(n.age) );
          }
        }
        
        const double act_unit = use_curries ? PhysicalUnits::curie : PhysicalUnits::becquerel;
        const string act_unit_str = use_curries ? "ci" : "bq";
        
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
          const double time_now = m_timeSpan * row / numerator;
          
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
              response.out() <<  x.numPerSecond;
            for( auto &x : mixture->betaPlusses(time_now, order) )
              response.out() << x.numPerSecond;
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
  : AuxWindow( "CSV Export", true, false ),
  m_parent( parent ),
  m_csvResouce( nullptr ),
  m_ageEdit( nullptr )
  {
    assert( m_parent );
    
    disableCollapse();
    setClosable( true );
    setResizable( false );
    finished().connect( m_parent, &DecayActivityDiv::deleteCsvDownloadGui );
    
    m_csvResouce = new DecayCsvResource( parent );
    
    const double displayTime = parent->timeToDisplayTill();
    const string timeSpanStr = PhysicalUnits::printToBestTimeUnits( displayTime );
    
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
    WLabel *label = new WLabel( "Time Span", el );
    el->setVerticalAlignment( Wt::AlignMiddle );
    label->setMargin( 2, Wt::Right );
    
    el = table->elementAt(0,1);
    m_ageEdit = new WLineEdit( timeSpanStr, el );
    label->setToolTip( "Enter time period you would like decay data for."
                       " Ex. '5.2 y', '52.3 s', '00:01:2.1', or '3.2d 15h'"
                       " would al be valid time periods."
                       " You can also enter time periods such as '3.2hl' or"
                       " '3.2 halflife' for multiples of the first parent"
                       " nuclides half life.");
    label->setBuddy( m_ageEdit );
    WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex, this );
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
    label = new WLabel( "Include:", el );
    
    el = table->elementAt(0,1);
    WCheckBox *cb = new WCheckBox( "act.", el );
    cb->setMargin( 5, Wt::Left );
    cb->setChecked( true );
    cb->checked().connect( boost::bind( &DecayCsvResource::setGiveActivities, m_csvResouce, true ) );
    cb->unChecked().connect( boost::bind( &DecayCsvResource::setGiveActivities, m_csvResouce, false ) );
    
    el = table->elementAt(0,2);
    cb = new WCheckBox( "xrays", el );
    cb->setMargin( 5, Wt::Left );
    cb->checked().connect( boost::bind( &DecayCsvResource::setGiveXrays, m_csvResouce, true ) );
    cb->unChecked().connect( boost::bind( &DecayCsvResource::setGiveXrays, m_csvResouce, false ) );
    
    el = table->elementAt(0,3);
    cb = new WCheckBox( "gammas", el );
    cb->setMargin( 5, Wt::Left );
    cb->checked().connect( boost::bind( &DecayCsvResource::setGiveGammas, m_csvResouce, true ) );
    cb->unChecked().connect( boost::bind( &DecayCsvResource::setGiveGammas, m_csvResouce, false ) );
    
    el = table->elementAt(0,4);
    cb = new WCheckBox( "alphas", el );
    cb->setMargin( 5, Wt::Left );
    cb->checked().connect( boost::bind( &DecayCsvResource::setGiveAlphas, m_csvResouce, true ) );
    cb->unChecked().connect( boost::bind( &DecayCsvResource::setGiveAlphas, m_csvResouce, false ) );
    
    el = table->elementAt(0,5);
    cb = new WCheckBox( "betas", el );
    cb->setMargin( 5, Wt::Left );
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
    
    
    show();
    centerWindow();
  }//CsvDownloadGui
  
  void updateTimeSpan()
  {
    double timespan = 0.0;
    try
    {
      const string input = m_ageEdit->valueText().toUTF8();
      
      const SandiaDecay::Nuclide *nuc = NULL;
      if( m_parent->m_nuclides.size() > 0 )
        nuc = m_parent->m_nuclideDB->nuclide( m_parent->m_nuclides[0].z, m_parent->m_nuclides[0].a, m_parent->m_nuclides[0].iso );
      
      if( nuc )
        timespan = PhysicalUnits::stringToTimeDurationPossibleHalfLife( input, nuc->halfLife );
      else
        timespan = PhysicalUnits::stringToTimeDuration( input );
    }catch(...)
    {
      timespan = m_parent->timeToDisplayTill();
      m_ageEdit->setText( PhysicalUnits::printToBestTimeUnits(timespan) );
    }
    
    m_csvResouce->setTimeSpan( timespan );
  }//void updateTimeSpan()
  
};//class CsvDownloadGui

DecayActivityDiv::DecayActivityDiv( InterSpec *viewer, Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_viewer( viewer ),
  m_nuclides( 0 ),
#if( DECAY_CHART_USE_X_AXIS_DATES )
  m_initialTimePicker( NULL ),
#endif
#if( DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
  m_fracActivityToEndAtEdit( NULL ),
#else
  m_displayTimeLength( NULL ),
#endif
  m_displayActivityUnitsCombo( NULL ),
  m_displayActivityUnitsLabel( NULL ),
  m_logYScale( NULL ),
  m_photopeakLogYScale( NULL ),
  m_showGridLines( NULL ),
#if( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK )
  m_pdfAnchor( NULL ),
#endif
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
  m_photoPeakChart( NULL ),
  m_photoPeakModel( NULL ),
  m_photopeakAgeSlider( NULL ),
  m_sliderCurrentAgeText( NULL ),
  m_sliderEndAgeText( NULL ),
  m_photoPeakYScaleFixed( NULL ),
  m_photoPeakShieldingZ( NULL ),
  m_photoPeakShieldingAD( NULL ),
  m_moreInfoDialog( NULL ),
  m_decayLegend( NULL ),
#if( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK )
  m_pdfResource( NULL ),
#endif
  m_calc( NULL ),
  m_csvDownloadDialog( NULL ),
  m_currentTimeUnits( -1.0 ),
  m_currentTimeRange( -1.0 ),
  m_currentNumXPoints( 250 ),
  m_width( -1 ),
  m_height( -1 ),
  m_nuclideDB( DecayDataBaseServer::database() ),
  m_currentMixture( new SandiaDecay::NuclideMixture() )
{
  setInline( false );
  addStyleClass( "DecayActivityDiv" );
  setLayoutSizeAware( true );
  
  wApp->useStyleSheet( "InterSpec_resources/NuclideDecayInfo.css");
  
  init();
}//DecayActivityDiv constructor


void DecayActivityDiv::layoutSizeChanged( int width, int height )
{
  m_width = width;
  m_height = height;
}

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
  
  //Initialize all member widgets; will assign parentage later on
//  m_nuclides( 0 ),
#if( DECAY_CHART_USE_X_AXIS_DATES )
  m_initialTimePicker              = new WDatePicker();
#endif
#if( DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
  m_fracActivityToEndAtEdit        = new WDoubleSpinBox();
#else
  m_displayTimeLength      = new WLineEdit( "" );
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex, m_displayTimeLength );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_displayTimeLength->setValidator(validator);
#endif
  m_displayActivityUnitsCombo      = new WComboBox();
  m_displayActivityUnitsLabel      = new WLabel( "Display Units" );
  m_logYScale                      = new WCheckBox( "Log-Y Scale" );
  m_showGridLines                  = new WCheckBox( "Grid Lines" );
#if( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK )
  m_pdfAnchor                      = new WAnchor();
#endif
  m_yAxisType                      = new WComboBox();
  m_parentNuclidesDiv              = new WContainerWidget();
  m_nuclidesAddedDiv               = new WContainerWidget();
  m_createNewNuclideButton         = new WPushButton( "Add Nuclide..." );
  m_clearNuclidesButton            = new WPushButton( "Remove All" );
  m_nuclideSelectDialog            = new AuxWindow( "Select Nuclide To Add" );
  
  m_nuclideSelect                  = new DecaySelectNuclide( m_nuclideDB, NULL, m_nuclideSelectDialog);
  m_decayLegend                    = new WContainerWidget();

  m_chartTabWidget                 = new WTabWidget();
  m_decayChart                     = new SrbActivityChart();
  m_decayModel                     = new SrbActivityModel( this );

  m_photoPeakChart                 = new Chart::WCartesianChart();
  m_photoPeakModel                 = new WStandardItemModel( this );

  m_photopeakAgeSlider             = new WSlider();
  m_sliderCurrentAgeText           = new WText( "", XHTMLUnsafeText );
  m_sliderEndAgeText               = new WText( "", XHTMLUnsafeText );
  m_photoPeakYScaleFixed           = new WCheckBox( "Fix Y-Range to Maximum" );

  m_photoPeakShieldingZ            = new WDoubleSpinBox();
  m_photoPeakShieldingAD           = new WDoubleSpinBox();


#if( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK )
  m_pdfResource                    = new ChartToImageResource( m_decayChart );
#endif

  //Sete elements names so we can style with css
  m_parentNuclidesDiv->addStyleClass( "m_parentNuclidesDiv" );
  m_nuclidesAddedDiv->addStyleClass( "m_nuclidesAddedDiv" );

  m_createNewNuclideButton->addStyleClass( "AddIcon" );
  m_createNewNuclideButton->addStyleClass( "m_createNewNuclideButton" );
  
  m_clearNuclidesButton->addStyleClass( "CrossIcon" );
  m_clearNuclidesButton->addStyleClass( "m_clearNuclidesButton" );
    
  m_displayActivityUnitsCombo->addStyleClass( "m_displayActivityUnitsCombo" );
  m_displayActivityUnitsLabel->addStyleClass( "m_displayActivityUnitsLabel" );
  m_logYScale->addStyleClass( "m_logYScale" );
  m_yAxisType->addStyleClass( "m_yAxisType" );
  m_showGridLines->addStyleClass( "m_showGridLines" );
  m_decayChart->addStyleClass( "m_decayChart" );
  m_chartTabWidget->addStyleClass( "m_chartTabWidget" );
  m_photoPeakYScaleFixed->addStyleClass( "m_photoPeakYScaleFixed" );
  m_sliderCurrentAgeText->addStyleClass( "m_sliderCurrentAgeText" );
  m_photoPeakShieldingZ->addStyleClass( "m_photoPeakShieldingZ" );
  m_photoPeakShieldingAD->addStyleClass( "m_photoPeakShieldingAD" );
  m_decayLegend->addStyleClass( "DecayLegend" );

#if( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK )
  m_pdfAnchor->addStyleClass( "m_pdfAnchor" );
#endif

#if( DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
  m_fracActivityToEndAtEdit->addStyleClass( "m_fracActivityToEndAtEdit" );
#else
  m_displayTimeLength->addStyleClass( "m_displayTimeLength" );
#endif

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
  m_nuclideSelectDialog->rejectWhenEscapePressed();
  m_nuclideSelectDialog->setModal( false );
  m_nuclideSelectDialog->disableCollapse();
//  m_nuclideSelectDialog->setClosable( true );
  if( m_height > 0 )
    m_nuclideSelectDialog->setMaximumSize( WLength::Auto, int(1.25*m_height) );
  m_createNewNuclideButton->clicked().connect( boost::bind( &WDialog::show,
                                                      m_nuclideSelectDialog ) );
  m_createNewNuclideButton->clicked().connect( m_nuclideSelect,
                                   &DecaySelectNuclide::setNuclideSearchToFocus );
  
  m_nuclideSelect->done().connect( m_nuclideSelectDialog, &AuxWindow::hide );
//  m_nuclideSelect->done().connect( this, &DecayActivityDiv::nuclideSelectDialogDone );
  m_nuclideSelectDialog->finished().connect( this,
                                    &DecayActivityDiv::nuclideSelectDialogDone );
  m_nuclideSelectDialog->centerWindow();
  
  m_nuclideSelectDialog->hide();
  
//  m_nuclideSelect->selected().connect(
//                      boost::bind( &DecayActivityDiv::addNuclide,
//                                   this, _1, _2, _3, _4, _5, _6 ) );
  m_nuclideSelect->selected().connect(
                                   boost::bind( &DecayActivityDiv::addTheNuclide,
                                                this, _1 ) );
  
  m_clearNuclidesButton->clicked().connect( this,
                                            &DecayActivityDiv::clearAllNuclides );


  m_decayChart->clicked().connect( this, &DecayActivityDiv::decayChartClicked );


  WGridLayout *decayLayout = new WGridLayout();
  WContainerWidget *decayDiv = new WContainerWidget();
  decayDiv->addStyleClass( "SrbActivityChartDiv" );
  decayDiv->addStyleClass( "decayDiv" );
  decayDiv->setLayout( decayLayout );
  
//  decayLayout->addWidget( m_decayChart, 0, 0, 1, 1 );
  ChartAndLegendHolder *chartHolder
                     = new ChartAndLegendHolder( m_decayChart, m_decayLegend );
  

  decayLayout->addWidget( chartHolder, 0, 0, 1, 1 );
  
  
  decayLayout->addWidget( displayOptionsDiv, 1, 0, 1, 1 );
  decayLayout->setRowStretch( 0, 5 );

#if( ADD_DECAY_CHAIN_CHART )
  //////////////////
  // Decay Chain
  m_decayChainChart = new DecayChainChart();
  WGridLayout *decayChainLayout = new WGridLayout();
  WContainerWidget *decayChainDiv = new WContainerWidget();
  decayChainDiv->setMinimumSize( 300, 350 );
  decayChainDiv->addStyleClass( "decayChainDiv" );
  decayChainDiv->setOverflow( WContainerWidget::OverflowAuto );
  decayChainDiv->setLayout( decayChainLayout );
  
  //decayChainLayout->addWidget(elementTemplate, 0, 0, 1, 1, AlignCenter);
  decayChainLayout->addWidget( m_decayChainChart, 0, 0, 1, 1 );
  decayChainLayout->setRowStretch( 0, 1 );
#endif //ADD_DECAY_CHAIN_CHART

  //////////////////
  // tab widget
  const WTabWidget::LoadPolicy loadPolicy = WTabWidget::LazyLoading; //WTabWidget::PreLoading;
  m_chartTabWidget->addTab( decayDiv, "Activity Chart", loadPolicy );
    //  m_chartTabWidget->addTab( photopeakDiv, "Photo Peaks", loadPolicy );
#if( ADD_DECAY_CHAIN_CHART )
  m_chartTabWidget->addTab( decayChainDiv, "Decay Chain", loadPolicy );
  m_chartTabWidget->currentChanged().connect( m_decayChainChart,
                                    &DecayChainChart::deleteMoreInfoDialog );
  m_decayChainChart->nuclideChanged().connect( this,
                            &DecayActivityDiv::manageActiveDecayChainNucStyling );
#endif

  m_chartTabWidget->currentChanged().connect( this,
                                        &DecayActivityDiv::deleteMoreInfoDialog );
  m_chartTabWidget->currentChanged().connect( this,
                            &DecayActivityDiv::manageActiveDecayChainNucStyling );

  m_calc = new DateLengthCalculator( this, NULL );
  m_chartTabWidget->addTab( m_calc, "Calculator", loadPolicy );
  
  
  WGridLayout *layout = new WGridLayout();
  layout->setContentsMargins( 3, 3, 3, 3 );
  layout->addWidget(  m_parentNuclidesDiv, 0, 0, 1, 1 );
  layout->addWidget(  m_chartTabWidget,    1, 0, 1, 1 );
  layout->setRowStretch( 1, 5 );
  WContainerWidget::setLayout( layout );
}//void DecayActivityDiv::init()


void DecayActivityDiv::showDecayTab()
{
  m_chartTabWidget->setCurrentIndex( 0 );
}//void showDecayTab()



void DecayActivityDiv::showPhotopeakTab()
{
  m_chartTabWidget->setCurrentIndex( 1 );
}//void showPhotopeakTab()



bool DecayActivityDiv::dirExists( const std::string &dir )
{
  struct stat st;
  const int err = stat( dir.c_str(), &st );
  if( err < 0 ) return false;
  return (st.st_nlink > 1);
}//bool dirExists( std::string file )

void DecayActivityDiv::setGridLineStatus()
{
  m_decayChart->showGridLines( m_showGridLines->isChecked() );
}//void DecayActivityDiv::setGridLineStatus()



Wt::WContainerWidget *DecayActivityDiv::initDisplayOptionWidgets()
{
#if( DECAY_CHART_USE_X_AXIS_DATES )
  m_initialTimePicker->setDate( WDate::currentDate() );
#endif

#if( DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
  m_fracActivityToEndAtEdit->setRange( 1.0E-12, 1.0 );
  m_fracActivityToEndAtEdit->setSingleStep( 0.01 );
  m_fracActivityToEndAtEdit->setText( "0.01" );
#else
#endif


  m_showGridLines->setUnChecked();
  m_logYScale->setUnChecked();
  updateYScale();
  m_logYScale->changed().connect( this, &DecayActivityDiv::updateYScale );
  m_showGridLines->changed().connect(this, &DecayActivityDiv::setGridLineStatus );

  for( const auto &nuclide : PhysicalUnits::sm_activityUnitNameValues )
  {
    m_displayActivityUnitsCombo->addItem( nuclide.first );
  }//for( each defined unit of activity )

  m_displayActivityUnitsCombo->addItem( "Arbitrary" );

  m_displayActivityUnitsCombo->setCurrentIndex( 7 );

#if( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK )
  m_pdfAnchor->setResource( m_pdfResource );
  m_pdfAnchor->setTarget( TargetNewWindow );
  string anchorText = string("Export ") + m_pdfResource->imageType();
  m_pdfAnchor->setText( anchorText );
#endif

#if( DECAY_CHART_USE_X_AXIS_DATES )
    WLabel *beginLabel = new WLabel( "Start Date: " );
//  beginLabel->setBuddy( m_initialTimePicker );
#endif


  WContainerWidget *displayOptionsDiv = new WContainerWidget();
  displayOptionsDiv->addStyleClass( "displayOptionsDiv" );
  WContainerWidget *displOptUpper = new WContainerWidget( displayOptionsDiv );
  
#if( DECAY_CHART_USE_X_AXIS_DATES )
  displOptUpper->addWidget( beginLabel );
  displOptUpper->addWidget( m_initialTimePicker );
#endif


#if( DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
  WLabel *endLabel = new WLabel( "Ending fraction of initial activity:" );
  endLabel->setBuddy( m_fracActivityToEndAtEdit );
  displOptUpper->addWidget( endLabel );
  displOptUpper->addWidget( m_fracActivityToEndAtEdit );
  m_fracActivityToEndAtEdit->changed().connect( this,
                                        &DecayActivityDiv::refreshDecayDisplay );
#else
  WLabel *endLabel = new WLabel( "Time Span" );
  endLabel->setBuddy( m_displayTimeLength );
  displOptUpper->addWidget( endLabel );
  displOptUpper->addWidget( m_displayTimeLength );
  m_displayTimeLength->changed().connect( this,
                                        &DecayActivityDiv::refreshDecayDisplay );
  m_displayTimeLength->changed().connect( this,
                            &DecayActivityDiv::updatePhotopeakSliderEndDateText );
  m_displayTimeLength->enterPressed().connect( this, &DecayActivityDiv::refreshDecayDisplay );
  m_displayTimeLength->enterPressed().connect( this, &DecayActivityDiv::updatePhotopeakSliderEndDateText );
  
  const char *tooltip = "<div>Age can be specified using a combination of time units, "
  "similar to '<b>5.3y 8d 22m</b>' or in half lives like "
  "'<b>2.5 HL</b>'.</div>"
  "<div>"
  "Acceptible time units: <b>year</b>, <b>yr</b>, <b>y</b>, <b>day</b>, <b>d</b>, <b>hrs</b>, <b>hour</b>, <b>h</b>, <b>minute</b>, "
  "<b>min</b>, <b>m</b>, <b>second</b>, <b>s</b>, <b>ms</b>, <b>microseconds</b>, <b>us</b>, <b>nanoseconds</b>, <b>ns</b>, or "
  "you can specify time period by <b>hh:mm:ss</b>. Half life units can be "
  "specified using <b>hl</b>, <b>halflife</b>, <b>halflives</b>, <b>half-life</b>, <b>half-lives</b>, "
  "<b>half lives</b>, or <b>half life</b>."
  "</div>"
  "<div>"
  "Half life units or time periods can not be mixed with "
  "other units, and if multiple nuclides are specified the first one is assumed. When multiple time periods are "
  "specified, they are summed, e.x. '1y6months 3m' is interpreted as "
  "18 months and 3 minutes"
  "</div>";
//  m_displayTimeLength->setToolTip( tooltip );
  
  const bool showToolTipInstantly = m_viewer ? InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer ) : false;
  if( m_viewer )
    HelpSystem::attachToolTipOn( m_displayTimeLength, tooltip, showToolTipInstantly );
#endif
  
  WLabel *yaxisTypeLabel = new WLabel( "Y-Axis Type:", displOptUpper );
  yaxisTypeLabel->addStyleClass( "yaxisTypeLabel" );
  
  displOptUpper->addWidget( m_yAxisType );
  for( YAxisType y = YAxisType(0); y < NumYAxisType; y = YAxisType(int(y)+1) )
  {
    switch( y )
    {
      case ActivityAxis: m_yAxisType->addItem( "Activity" );   break;
      case GammasAxis:   m_yAxisType->addItem( "Gamma Rate" ); break;
      case BetasAxis:    m_yAxisType->addItem( "Beta Rate" );  break;
      case AlphasAxis:   m_yAxisType->addItem( "Alpha Rate" ); break;
      case NumYAxisType:                                       break;
    }//switch( y )
  }//for( Y-Axis type )
  
  WContainerWidget *displOptLower = new WContainerWidget( displayOptionsDiv );
  displOptLower->setMargin( 5, Wt::Top );
  displOptLower->addWidget( m_displayActivityUnitsLabel );
  displOptLower->addWidget( m_displayActivityUnitsCombo );
  
  displOptLower->addWidget( m_logYScale );
  displOptLower->addWidget( m_showGridLines );
  
  m_yAxisType->setCurrentIndex( ActivityAxis );
  m_yAxisType->activated().connect( this,&DecayActivityDiv::refreshDecayDisplay );
  
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
  if( app->viewer()->isSupportFile() )
  {
    WText *csvbutton = new WText( "csv&hellip;" );
    csvbutton->addStyleClass( "csvAnchor" );
    displOptLower->addWidget( csvbutton );
  
    csvbutton->clicked().connect( this, &DecayActivityDiv::createCsvDownloadGui );
  }

#if( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK )
  displOptLower->addWidget( m_pdfAnchor );
#endif

#if( DECAY_CHART_USE_X_AXIS_DATES )
  m_initialTimePicker->changed().connect( this,
                                        &DecayActivityDiv::refreshDecayDisplay );
  m_initialTimePicker->changed().connect(
                                   boost::bind( &WDatePicker::setPopupVisible,
                                                m_initialTimePicker, false ) );
#endif


  m_displayActivityUnitsCombo->changed().connect( this,
                                        &DecayActivityDiv::refreshDecayDisplay );

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
#if( DECAY_CHART_USE_X_AXIS_DATES )
  m_decayChart->setType( Chart::CategoryChart );
  m_decayChart->axis(Chart::XAxis).setScale( Chart::CategoryScale );
  m_decayChart->axis(Chart::XAxis).setLabelFormat( "MMM dd yyyy hh:mm" );
#endif

  m_decayChart->setPlotAreaPadding( 80, Left );
  m_decayChart->setPlotAreaPadding( 50, Bottom );
  m_decayChart->setPlotAreaPadding( 20, Right );
  m_decayChart->setPlotAreaPadding( 5, Top );
  //do not display the legend to the right of the chat
  m_decayChart->setLegendEnabled( false);
//  m_decayChart->setMargin(0);
//  m_decayChart->initLayout();

  m_decayChart->setToolTip( "Click for more information" );
//  m_decayChart->mouseMoved().connect( this, &DecayActivityDiv::updateMouseOver );

  //Now lets init the phoptopeak chart
  m_photoPeakChart->setModel( m_photoPeakModel );
  m_photoPeakChart->setXSeriesColumn( 0 );
  m_photoPeakChart->setMinimumSize( WLength(250), WLength(100) );
//  m_photoPeakChart->setType( Chart::CategoryChart );
  m_photoPeakChart->setType( Chart::ScatterPlot );
  m_photoPeakChart->axis(Chart::XAxis).setScale( Chart::LinearScale );
  m_photoPeakChart->setPlotAreaPadding( 80, Left );
  m_photoPeakChart->setPlotAreaPadding( 50, Bottom );

  m_photoPeakChart->setPlotAreaPadding( 20, Right );
  m_photoPeakChart->setPlotAreaPadding( 5, Top );
  m_photoPeakChart->setLegendEnabled( false );
  m_photoPeakChart->axis(Chart::XAxis).setTitle( "Photopeak Energy (keV)" );
  m_photoPeakChart->axis(Chart::XAxis).setRange( 0.0, 3000.0 );
  m_photoPeakChart->initLayout();
}//void initCharts()



double DecayActivityDiv::timeToDisplayTill()
{
#if( DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
  const double dispTo = m_fracActivityToEndAtEdit->value();
  return findTimeForActivityFrac( m_currentMixture, dispTo );
#else
  
  std::string txt = m_displayTimeLength->text().narrow();
  
  const SandiaDecay::Nuclide *nuc = NULL;
  if( m_nuclides.size() > 0 )
    nuc = m_nuclideDB->nuclide( m_nuclides[0].z, m_nuclides[0].a, m_nuclides[0].iso );
  
  try
  {
    if( nuc )
      return PhysicalUnits::stringToTimeDurationPossibleHalfLife( txt, nuc->halfLife, SandiaDecay::second);
    else
      return PhysicalUnits::stringToTimeDuration( txt );
  }catch( std::exception & )
  {
    const double fracT0 = 0.1;
    const double finalTime = findTimeForActivityFrac( m_currentMixture, fracT0 );
    
    txt = PhysicalUnits::printToBestTimeUnits( finalTime );
    m_displayTimeLength->setText( txt );
    return finalTime;
  }
#endif
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
                                     nuc.age, nuc.activity, nuc.useCurrie );
    m_nuclideSelect->setAddButtonToAccept();
    break;
  }//for( Nuclide &nuc : m_nuclides )
}//void sourceNuclideDoubleClicked( Wt::WContainerWidget *w );


void DecayActivityDiv::addTheNuclide( const NuclideSelectedInfo &n )
{
  addNuclide( n.z, n.a, n.metasable, n.activity, n.useCurrie, n.initialAge );
}//void addTheNuclide( const NuclideSelectedInfo &nuc )

void DecayActivityDiv::addNuclide( const int z, const int a, const int iso,
                                 const double activity, const bool useCurrie,
                                 const double age )
{
  //See if we are actually editing a Nuclide
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
  nuclide.useCurrie = useCurrie;
  nuclide.iso       = iso;
  nuclide.age       = age;
  nuclide.display   = new WContainerWidget();
  nuclide.display->setInline( true );
  nuclide.display->addStyleClass( "Nuclide" );

  if( !editedNuclide )
  {
    m_nuclidesAddedDiv->addWidget( nuclide.display );
  }else
  {
    m_nuclidesAddedDiv->insertBefore( nuclide.display, editedNuclide );
    removeNuclide( editedNuclide );
  }
  
  const SandiaDecay::Element *element = m_nuclideDB->element( z );

  string name = (element ? element->symbol : string(""));

#if( ADD_DECAY_CHAIN_CHART )
  const SandiaDecay::Nuclide *nuc = m_nuclideDB->nuclide( z, a );
  if( nuc )
    nuclide.display->clicked().connect(
         boost::bind(&DecayChainChart::setNuclide, m_decayChainChart, nuc) );
#endif
  
  nuclide.display->doubleClicked().connect(
                    boost::bind( &DecayActivityDiv::sourceNuclideDoubleClicked,
                                 this, nuclide.display) );
  nuclide.display->setToolTip( "Double click to edit this source nuclide. "
                               "Single click to make the decay chain to display"
                               " this nuclide" );

  stringstream label;
  label << "<sup><font size=\"1.5\">" << a;
  if( nuclide.iso ) label << "m";
  if( nuclide.iso > 1 ) label << nuclide.iso;
  label << "</font></sup>" << name;
  label << " " << fixed << setprecision(2) << PhysicalUnits::printToBestActivityUnits(activity, 2, useCurrie );

  if( age > 0.0 )
  {
    //const PhysicalUnits::UnitNameValuePair units = PhysicalUnits::bestTimeUnitShortHtml( age );
    //label << " " << fixed << setprecision(2) <<  age/units.second << units.first;
    label << " " << PhysicalUnits::printToBestTimeUnits(age);
  }//if( age > 0.0 )

  new WText( label.str(), XHTMLUnsafeText, nuclide.display );

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

#if( !DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
  setTimeLimitToDisplay();
#endif

#if( ADD_DECAY_CHAIN_CHART )
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  m_decayChainChart->setNuclide(db->nuclide(z, a, iso));
#endif

  refreshDecayDisplay();
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

#if( ADD_DECAY_CHAIN_CHART )
  if (to_be_removed == m_nuclides.size() - 1)
  {
    if (m_nuclides.size() == 1)
      m_decayChainChart->setNuclide(NULL);
    else
    {
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      const DecayActivityDiv::Nuclide nuc = m_nuclides[m_nuclides.size() - 2];
      m_decayChainChart->setNuclide(db->nuclide(nuc.z, nuc.a, nuc.iso));
    }
  }
#endif

  delete m_nuclides[to_be_removed].display;
  m_nuclides.erase( m_nuclides.begin() + to_be_removed );

#if( !DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
  setTimeLimitToDisplay();
#endif

  refreshDecayDisplay();
}//void DecayActivityDiv::removeNuclide( Wt::WContainerWidget *frame )


void DecayActivityDiv::clearAllNuclides()
{
  for( const Nuclide &nuclide : m_nuclides )
    delete nuclide.display;
  m_nuclides.clear();

#if( !DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
  setTimeLimitToDisplay();
#endif
    
#if( ADD_DECAY_CHAIN_CHART )
  m_decayChainChart->setNuclide(NULL);
#endif

  deleteMoreInfoDialog();
  
#if( ADD_DECAY_CHAIN_CHART )
  if( m_decayChainChart )
    m_decayChainChart->deleteMoreInfoDialog();
#endif

  refreshDecayDisplay();
}//void clearAllNuclides();





void DecayActivityDiv::updateInitialMixture() const
{
  using SandiaDecay::NuclideActivityPair;
  using SandiaDecay::NuclideNumAtomsPair;

  m_currentMixture->clear();

  for( const Nuclide &nuc : m_nuclides )
  {
    const SandiaDecay::Nuclide *dbnuclide = m_nuclideDB->nuclide( nuc.z, nuc.a,
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
  if( m_logYScale->isChecked() )
    m_decayChart->axis(Chart::YAxis).setScale( Chart::LogScale );
  else
    m_decayChart->axis(Chart::YAxis).setScale( Chart::LinearScale );
}//void updateYScale()

void DecayActivityDiv::setPhotoPeakChartLogY( bool logy )
{
  if( logy )
  {
    if( m_photoPeakChart->axis(Chart::YAxis).minimum() < 0.0001 )
      m_photoPeakChart->axis(Chart::YAxis).setMinimum( 0.0001 );
    m_photoPeakChart->axis(Chart::YAxis).setScale( Chart::LogScale );
  }else
  {
    m_photoPeakChart->axis(Chart::YAxis).setScale( Chart::LinearScale );
  }

  if( m_photopeakLogYScale->isChecked() != logy )
    m_photopeakLogYScale->setChecked( logy );
}//void setPhotoPeakChartLogY( bool logy )


void DecayActivityDiv::displayMoreInfoPopup( const double time )
{
  WContainerWidget *summary = isotopesSummary( time );

  string title = "Summary at t=";
  title += PhysicalUnits::printToBestTimeUnits( time );

  if( m_moreInfoDialog )
  {
    m_moreInfoDialog->contents()->clear();
    m_moreInfoDialog->setWindowTitle( title );
  }else
  {
    m_moreInfoDialog = new AuxWindow( title );
    m_moreInfoDialog->setModal( false );
//    m_moreInfoDialog->setResizable( true );
//    m_moreInfoDialog->setMaximumSize( WLength(90.0,WLength::Percentage),
//                                      WLength(90.0,WLength::Percentage) );
    m_moreInfoDialog->setClosable( true );
    m_moreInfoDialog->rejectWhenEscapePressed();
    m_moreInfoDialog->finished().connect(
                  boost::bind( &DecayActivityDiv::deleteMoreInfoDialog, this ) );
      
    WPushButton *ok = m_moreInfoDialog->addCloseButtonToFooter();
    ok->clicked().connect(boost::bind( &DecayActivityDiv::deleteMoreInfoDialog, this ) );

  }//if( m_moreInfoDialog ) / else

  m_moreInfoDialog->resizeToFitOnScreen();
  WGridLayout *layout = new WGridLayout();
  m_moreInfoDialog->contents()->setLayout( layout );
  layout->addWidget( summary, 0, 0 );
  
  layout->setRowStretch( 0, 10 );

  //I dont really care about centering the dialog, but without doing the bellow
  //  the scrolling of the contents can be buggy
//  const string jsthis = "$('#" + dialog->id() + "')";
//  stringstream moveJs;
//  moveJs << "var el = " << dialog->jsRef() << ";"
//         << "var ws = " << wApp->javaScriptClass() << ".WT.windowSize();"
//         << "el.style.top = Math.round( ( ws.y - " << jsthis << ".height() ) / 2 )+'px';"
//         << "el.style.left = Math.round( ( ws.x - " << jsthis << ".width() ) / 2 )+'px';";
//  doJavaScript( moveJs.str() );
  
  m_moreInfoDialog->setHidden( false );
}//void displayMoreInfoPopup( const double time )


void DecayActivityDiv::photopeakDisplayMoreInfo()
{
  const double time = photopeakSliderTime();
  displayMoreInfoPopup( time );
}//void photopeakDisplayMoreInfo()


#if( !DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
void DecayActivityDiv::checkTimeRangeValid()
{
  std::string txt = m_displayTimeLength->text().narrow();
  
  const SandiaDecay::Nuclide *nuc = NULL;
  if( m_nuclides.size() > 0 )
    nuc = m_nuclideDB->nuclide( m_nuclides[0].z, m_nuclides[0].a, m_nuclides[0].iso );
  
  double timelen = 0.0;
  
  try
  {
    if( nuc )
      timelen = PhysicalUnits::stringToTimeDurationPossibleHalfLife( txt, nuc->halfLife, SandiaDecay::second);
    else
      timelen = PhysicalUnits::stringToTimeDuration( txt );
  }catch( std::exception & )
  {
    timelen = findTimeForActivityFrac( m_currentMixture, 0.1 );
    txt = PhysicalUnits::printToBestTimeUnits( timelen );
    m_displayTimeLength->setText( txt );
  }

  
  double minhl = DBL_MAX;
  for( const Nuclide &nuc : m_nuclides )
  {
    const SandiaDecay::Nuclide *dbnuclide = m_nuclideDB->nuclide( nuc.z, nuc.a,
                                                                 nuc.iso );
    minhl = std::min( minhl, dbnuclide->halfLife );
  }//for( const Nuclide &nuc : m_nuclides )
  
  minhl = std::max( minhl, 1.0E-6*SandiaDecay::second );
  minhl = std::min( minhl, 1000.0*SandiaDecay::second );
  if( timelen < 0.001*minhl )
  {
    timelen = 0.001*minhl;
    txt = PhysicalUnits::printToBestTimeUnits( timelen );
    m_displayTimeLength->setText( txt );
  }//if( nunits*unitpair.second < 0.001*minhl )
}//void DecayActivityDiv::checkTimeRangeValid()

void DecayActivityDiv::updatePhotopeakSliderEndDateText()
{
  using namespace PhysicalUnits;
  checkTimeRangeValid();
  
  m_sliderEndAgeText->setText( m_displayTimeLength->text() );
}//void updatePhotopeakSliderEndDateText()
#endif


#if( !DECAY_CHART_LIMIT_TIME_BY_ACTIVITY )
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
  
  string txt = PhysicalUnits::printToBestTimeUnits( finalTime );
  m_displayTimeLength->setText( txt );

  updatePhotopeakSliderEndDateText();
}//void DecayActivityDiv::setTimeLimitToDisplay()
#endif


void DecayActivityDiv::setDecayChartTimeRange( double finalTime )
{
  string txt = PhysicalUnits::printToBestTimeUnits( finalTime );
  m_displayTimeLength->setText( txt );
  
  updatePhotopeakSliderEndDateText();
}//void setDecayChartTimeRange()



double DecayActivityDiv::attentuationCoeff( const double energy )
{
  if( (m_photoPeakShieldingZ->validate() != WValidator::Valid)
      || (m_photoPeakShieldingAD->validate() != WValidator::Valid) )
    return 0.0;

  const int atomicNumber = static_cast<int>(m_photoPeakShieldingZ->value() + 0.5f);
  const float arealDensity = static_cast<float>(m_photoPeakShieldingAD->value());

  if( (atomicNumber<=0) || (arealDensity<=0.0f) )
    return 1.0;

  if( atomicNumber>98 )
    return 0.0;

  const double xsenergy = (energy/SandiaDecay::keV)
                              * PhysicalUnits::keV;
  const double mu = MassAttenuation::massAttenuationCoeficient( atomicNumber,
                                        static_cast<float>(xsenergy) );

  const static double gPerCm2 = PhysicalUnits::g / PhysicalUnits::cm2;
  return exp( -mu * (arealDensity*gPerCm2) );
}//double attentuationCoeff( const double energy )


void DecayActivityDiv::setPhotopeakXScaleRange()
{
  const int nrow = m_photoPeakModel->rowCount();
  if( nrow < 2 )
  {
    m_photoPeakChart->axis(Chart::XAxis).setRange( 0.0, 3000.0 );
    return;
  }//if( nrow < 2 )

  const WModelIndex index = m_photoPeakModel->index( nrow-1, 0 );
  boost::any last_data = m_photoPeakModel->data( index );
  double energy = 3000.0;
  try{ energy = boost::any_cast<double>( last_data ); }catch(...){}

  //Display to at leat 3 MeV to be consistent with converntion, but if isotope
  //  has heigher energy gamma lines, display those too, up to 10 MeV.
  energy = max( energy, 3000.0 );
  energy = min( energy, 10000.0 );

  m_photoPeakChart->axis(Chart::XAxis).setRange( 0.0, energy );
}//void DecayActivityDiv::setPhotopeakXScaleRange()



void DecayActivityDiv::setPhotopeakYScaleRange()
{
  double minValue = DBL_MAX, maxValue = -1.0;

  if( m_photoPeakYScaleFixed->isChecked() )
  {
    const double maxDiplayTime = timeToDisplayTill();
    const double dt = maxDiplayTime / 100.0;

    for( double age = 0; age <= maxDiplayTime; age += dt )
    {
      const vector<SandiaDecay::EnergyRatePair> gammas
           = m_currentMixture->gammas( age,
                                     SandiaDecay::NuclideMixture::OrderByAbundance, true);

      if( gammas.size() )
      {
        const double maxAtten = attentuationCoeff( gammas[0].energy );
        const double minAtten = attentuationCoeff( gammas.back().energy );
        maxValue = max( maxValue, maxAtten*gammas[0].numPerSecond );
        minValue = min( minValue, minAtten*gammas.back().numPerSecond );
      }//if( gammas.size() )
    }//for( loop over times to find max value )
  }else
  {
    const double age = photopeakSliderTime();
    const vector<SandiaDecay::EnergyRatePair> gammas
         = m_currentMixture->gammas( age,
                                   SandiaDecay::NuclideMixture::OrderByAbundance, true );

    if( gammas.size() )
    {
      const double maxAtten = attentuationCoeff( gammas[0].energy );
      const double minAtten = attentuationCoeff( gammas.back().energy );
      maxValue = max( maxValue, maxAtten*gammas[0].numPerSecond );
      minValue = 0.1*min( minValue, minAtten*gammas.back().numPerSecond );
    }//if( gammas.size() )

//    m_photoPeakChart->axis(Chart::YAxis).setAutoLimits( Chart::MaximumValue );
//    m_photoPeakChart->axis(Chart::YAxis).setAutoLimits( Chart::MinimumValue );
  }//if( fixed y axis ) / else

  if( minValue==maxValue )
  {
    if( minValue == 0.0 )
    {
      minValue = 0.0;
      maxValue = 1.0;
    }else
    {
      minValue = 0.1 * minValue;
      maxValue = 1.1 * maxValue;
    }//
  }//if( fabs(minValue-maxValue) < 0.00000001 )

  if( (m_photoPeakChart->axis(Chart::YAxis).scale() != Chart::LinearScale)
      && (minValue <= 0.00001) )
    minValue = 0.00001;

  if( m_photoPeakChart->axis(Chart::YAxis).scale() == Chart::LinearScale )
    minValue = 1.0;

  if( maxValue > 0.0 )
    m_photoPeakChart->axis(Chart::YAxis).setMaximum( maxValue );

  if( (minValue > 0.0) && (minValue < DBL_MAX) )
    m_photoPeakChart->axis(Chart::YAxis).setMinimum( minValue );
  else
    m_photoPeakChart->axis(Chart::YAxis).setAutoLimits( Chart::MinimumValue );
}//void setPhotopeakYScaleRange()



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
#if( ADD_DECAY_CHAIN_CHART )
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
#endif
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

  string mixtureInfo = m_currentMixture->info( time );
  mixtureInfo = "<pre class=\"NuclideInfo\">\n"
                + mixtureInfo + "\n</pre>";
  new WText( mixtureInfo, XHTMLUnsafeText, cont );

  new WText( "<br><br><H3>Further information about the nuclides:</H3>",
             XHTMLUnsafeText, cont );

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

  stringstream textstr;
  textstr << "<pre class=\"NuclideInfo\">\n";

  if( nuclide )
    textstr << SandiaDecay::human_str_summary(*nuclide) << "\n";

  textstr << "</pre>\n";

  new WText( textstr.str(), XHTMLUnsafeText, cont );

  return cont;
}//WContainerWidget *nuclideInformation( const Nuclide &nuclide ) const;


double DecayActivityDiv::photopeakSliderTime()
{
  const double minSlider = static_cast<double>(m_photopeakAgeSlider->minimum());
  const double maxSlider = static_cast<double>(m_photopeakAgeSlider->maximum());
  const double sliderValue = static_cast<double>(m_photopeakAgeSlider->value());
  const double sliderFraction = (sliderValue-minSlider) / (maxSlider-minSlider);
  const double maxDiplayTime = timeToDisplayTill();
  return sliderFraction * maxDiplayTime;
}//double DecayActivityDiv::photopeakSliderTime()


void DecayActivityDiv::setPhopeakSliderTime( double time )
{
  time = max( time, 0.0 );

  {//begin codeblock to set display time range to make sure resolution of slider
   // doesnt significantly effect desired display range
    //m_displayTimeLength->setValue( 20 );
  }//end codeblock to set display time range

  const double maxDisplayTime = timeToDisplayTill();
  const double minSlider = static_cast<double>(m_photopeakAgeSlider->minimum());
  const double maxSlider = static_cast<double>(m_photopeakAgeSlider->maximum());
  const double sliderRange = maxSlider - minSlider;

  const double fractionSlider = time / maxDisplayTime;
  const double slidVal =  minSlider + fractionSlider*sliderRange;
  m_photopeakAgeSlider->setValue( static_cast<int>( floor(0.5 + slidVal)  ));

  refreshPhotopeakDisplay();
}//void setPhopeakSliderTime( double time )


void DecayActivityDiv::refreshPhotopeakDisplay()
{
  m_photoPeakModel->clear();

  if( m_nuclides.empty() )
    return;

  const double age = photopeakSliderTime();
  m_sliderCurrentAgeText->setText( PhysicalUnits::printToBestTimeUnits(age) );


  const vector<SandiaDecay::EnergyRatePair> gammas
                        = m_currentMixture->gammas( age,
                             SandiaDecay::NuclideMixture::OrderByEnergy, true );
  const int nGamma = static_cast<int>( gammas.size() );

  m_photoPeakModel->insertRows( 0, nGamma );
  m_photoPeakModel->insertColumns( 0, 2 );

  for( int i = 0; i < nGamma; ++i )
  {
    const SandiaDecay::EnergyRatePair &gamma = gammas[i];
    const double energy_keV = gamma.energy/SandiaDecay::keV;
    const double attenuation = attentuationCoeff( gamma.energy );
    const double numPerSecond = attenuation * gamma.numPerSecond;

    m_photoPeakModel->setData( i, 0, boost::any(energy_keV));
    m_photoPeakModel->setData( i, 1, boost::any(numPerSecond) );
  }//for( loop over gammas, i )

  Chart::WDataSeries series(1, Chart::BarSeries, Chart::YAxis);
  WPen photopeakPen;
  photopeakPen.setWidth(2);
  series.setPen( photopeakPen );
  m_photoPeakChart->addSeries( series );
  setPhotopeakXScaleRange();
  setPhotopeakYScaleRange();
}//void refreshPhotopeakDisplay()



void DecayActivityDiv::refreshDecayDisplay()
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
    cout << "Nuc:'" << nuc << "'=" << m_decayModel->showSeries( column ) << endl;
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

  //TODO: call refreshPhotopeakDisplay() seperately from refreshDecayDisplay()
  refreshPhotopeakDisplay();

  if( m_nuclides.empty() )
  {
    m_decayLegend->clear();
    return;
  }//if( m_nuclides.empty() )

#if( DECAY_CHART_USE_X_AXIS_DATES )
  const WDateTime t0( m_initialTimePicker->date() );
#endif

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
  m_decayChart->axis(Chart::XAxis).setTitle( xUnitsPair.first + " after t0" );
  
  for( int row = 0; row < nRows; ++row )
  {
    const double row_time = (row * dt);

#if( DECAY_CHART_USE_X_AXIS_DATES )
    try
    {
      const int dt_seconds = floor( (row * dt)/second + 0.5 );
      const WDateTime x = t0.addSecs( dt_seconds );
      if( !x.isValid() )
        throw std::runtime_error( "" );
      m_decayModel->setData( row, 0, boost::any( x ) );
      //m_decayModel->setData( row, 0, boost::any( x ), ToolTipRole );
    }catch(...)
    {
      cerr << "refreshDecayDisplay(): Warning, couldnt add " << dt_seconds
           << " secs to " << t0.toString() << " or issue with m_decayModel"
           << endl;
      continue;
    }//try / catch
#else
    //We will set the x-axis data as a formatted string since
    //  WAxis::setLabelFormat( "%.3g" ); doesnt seem to work
    stringstream labelText;
    labelText << fixed << setprecision(3) << (row * dt)/tunit;
    const WString label( labelText.str() );
    m_decayModel->setData( row, 0, boost::any( label ) );
#endif

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

      if( yval==0.0 || IsInf(yval) || IsNan(yval) )
        continue;

      totalActivities[row] += yval;
      maxActivity = std::max( maxActivity, yval );
      minActivity = std::min( minActivity, yval );
      
//      if( activity >= (0.00001*endActivity) )
      m_decayModel->setData( row, column, boost::any( yval ) );
      
      //WString tt = "My Tool Tip";
      //m_decayModel->setData( row, column, boost::any( tt ), ToolTipRole );
    }//for( loop over nuclides to add )
  }//for( loop over time points to add )

  //Now put in the sum of all the activities
  const WString dateHeader = "Date (" + xUnitsPair.first + ")";
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
                                   boost::any( WString( "Total Activity" ) ) );
      m_decayChart->axis(Chart::YAxis).setTitle( unitStr );
    break;
    case GammasAxis:
      m_decayModel->setHeaderData( nElements+1,
                                   boost::any( WString( "Total Gammas" ) ) );
      m_decayChart->axis(Chart::YAxis).setTitle( "Gammas/Second" );
    break;
    case BetasAxis:
      m_decayModel->setHeaderData( nElements+1,
                                   boost::any( WString( "Total Betas" ) ) );
      m_decayChart->axis(Chart::YAxis).setTitle( "Betas/Second" );
    break;
    case AlphasAxis:
      m_decayModel->setHeaderData( nElements+1,
                                   boost::any( WString( "Total Alphas" ) ) );
      m_decayChart->axis(Chart::YAxis).setTitle( "Alphas/Second" );
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
  
  
//  double miny = 0.1*endActivity;
//  double miny = 0.9*minActivity;
//  const double maxy = 1.1*maxActivity;
  
  //If it wont change the dynamic range of the chart much anyways, might as well
  //  anchor the y-axis to zero
//  if( miny < 0.05*maxy )
//    miny = 0.0; //never read
//  m_decayChart->axis(Chart::YAxis).setRange( miny, maxy );
  m_decayChart->axis(Chart::YAxis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );

  //fill in the 'cache' variables for updateMouseOver(...)
  m_currentTimeUnits  = tunit;
  m_currentTimeRange  = maxDiplayTime;
  
  if( m_calc )
    m_calc->setTimeRangeTxt( m_displayTimeLength->text().narrow() );
}//void refreshDecayDisplay()



void DecayActivityDiv::userSetShowSeries( int series, bool show )
{
  m_decayModel->setShowSeries( series, show );
  
  
//  double miny = 0.9*minActivity/actunit;
//  const double maxy = 1.1*maxActivity/actunit;
  //If it wont change the dynamic range of the chart much anyways, might as well
  //  anchor the y-axis to zero
//  if( miny < 0.05*maxy )
//    miny = 0.0;
//  m_decayChart->axis(Chart::YAxis).setRange( miny, maxy );
  m_decayChart->axis(Chart::YAxis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
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
    if( column != (nColumns-1) )
      pen.setColor( seriesColor( column ) );
    else
      pen.setStyle( DashLine );
    
    series.setPen( pen );
    
    m_decayChart->addSeries( series );
  }//for( int column = 1; column < nElements; ++column )
  
  m_decayLegend->clear();
  
  
  
  cout << "nColumns=" << nColumns <<endl;
  
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
  tip << "Click for more information - "
      << fixed << setprecision(4) << mouseTime/tUnit
      << " " << tUnitStr <<  " after T0 we have:";

  double totalActivity = 0.0;

  const std::vector<SandiaDecay::NuclideTimeEvolution> &evolutions
                              = m_currentMixture->decayedToNuclidesEvolutions();


  const size_t nDaughter = evolutions.size();
  for( size_t i = 0; i < nDaughter; ++i )
  {
    const size_t index = nDaughter - i - 1;
    const SandiaDecay::Nuclide *el = evolutions[index].nuclide;
    const double activity = evolutions[index].activity( mouseTime );
    totalActivity += activity;
    tip << "\n  " << el->symbol << " " << activity/actunit << " " << actUnitStr;
  }//for( loop over population daughters )

  tip << "\n  Total Activity: " << totalActivity/actunit << " " << actUnitStr;

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
