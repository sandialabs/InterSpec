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


#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>
#include <Wt/Chart/WCartesianChart>

#include "InterSpec/DrfChart.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace Wt;


DrfChart::DrfChart( WContainerWidget *parent )
: Wt::Chart::WCartesianChart( parent )
{
  setBackground(Wt::WColor(220, 220, 220));
  setXSeriesColumn(0);
  setType(Wt::Chart::ScatterPlot);
  
  //Before, m_efficiencyModel was actually a memory leak, because it didnt have
  //  a parent (m_chart in this case), so everytime you created a new one, the
  //  old one would be leaked since there was no longer a refernce to it anywhere
  m_efficiencyModel = new WStandardItemModel( this );
  setModel( m_efficiencyModel );
  
  m_chartEnergyLineColor = WColor("#B02B2C"); //Ruby on Rails red
  m_chartFwhmLineColor = WColor("#3F4C6B");  //Mozilla blue
  
  //We should check the color theme for colors
  InterSpec *viewer = InterSpec::instance();
  
  if( viewer )
    viewer->colorThemeChanged().connect( this, &DrfChart::handleColorThemeChange );
  
  auto theme = viewer ? viewer->getColorTheme() : nullptr;
  handleColorThemeChange( theme );
  
  const bool is_phone = (!viewer || viewer->isPhone() || (viewer->renderedHeight() < 500 ));
  
  //setAutoLayoutEnabled(); //Leaves a lot of room at the top, but maybe because font-metrics not available?
  setPlotAreaPadding(5, Wt::Top);
  
  if( is_phone )
  {
    setPlotAreaPadding(20, Wt::Bottom);
    setPlotAreaPadding(50, Wt::Right | Wt::Left);
    
    WFont labelFont = axis(Wt::Chart::XAxis).labelFont();
    labelFont.setSize( WFont::Size::XXSmall );
    axis(Wt::Chart::XAxis).setLabelFont( labelFont );
    axis(Wt::Chart::Y1Axis).setLabelFont( labelFont );
    axis(Wt::Chart::Y2Axis).setLabelFont( labelFont );
    
    WFont titleFont = axis(Wt::Chart::XAxis).titleFont();
    titleFont.setSize( WFont::Size::XSmall );
    axis(Wt::Chart::XAxis).setTitleFont( titleFont );
    axis(Wt::Chart::Y1Axis).setTitleFont( titleFont );
    axis(Wt::Chart::Y2Axis).setTitleFont( titleFont );
    
    // titleOffset doesnt seem to have an effect, so leaving commented out (would be nice to make title closer to labels)
    // axis(Wt::Chart::Y1Axis).setTitleOffset( 1 );
    // axis(Wt::Chart::Y2Axis).setTitleOffset( 1 );
    
    setMinimumSize(WLength(200), WLength(100));
  }else
  {
    setPlotAreaPadding(55, Wt::Right | Wt::Left);
    setPlotAreaPadding(60, Wt::Bottom);
    axis(Wt::Chart::XAxis).setTitle("Energy (keV)");
    
    setMinimumSize(WLength(350), WLength(200));
  }//if( is_phone ) / else
  
  axis(Wt::Chart::Y1Axis).setVisible(true);
  axis(Wt::Chart::Y1Axis).setTitle("Efficiency");
  
  axis(Wt::Chart::Y1Axis).setTitleOrientation( Wt::Vertical );
  axis(Wt::Chart::Y2Axis).setTitleOrientation( Wt::Vertical );
  
  setLegendEnabled( true );
  
  setLegendLocation(Wt::Chart::LegendInside, Wt::Top, Wt::AlignRight);
  
#if( WT_VERSION >= 0x3040000 )  //I'm not too sure when these features became available, but I dont think in Wt 3.3.4
  setCurveManipulationEnabled( true );
  
  Wt::Chart::WheelActions wheelAction;
  //InteractiveAction: ZoomY, ZoomXY, ZoomMatching, PanX, PanY, PanMatching
  wheelAction[KeyboardModifier::NoModifier] = Chart::InteractiveAction::ZoomX;
  
  setWheelActions(wheelAction);
  
  setPanEnabled(true);
  setZoomEnabled(true);
  
  //setRubberBandEffectEnabled(true);
#endif
}//DrfChart constructor


void DrfChart::handleColorThemeChange( std::shared_ptr<const ColorTheme> theme )
{
  if( !theme )
    return;
  
  if( !theme->foregroundLine.isDefault() )
    m_chartEnergyLineColor = theme->foregroundLine;
  else
    m_chartEnergyLineColor = WColor("#B02B2C"); //Ruby on Rails red
  
  
  if( !theme->backgroundLine.isDefault() )
    m_chartFwhmLineColor = theme->backgroundLine;
  else
    m_chartFwhmLineColor = WColor("#3F4C6B");  //Mozilla blue
  
  const WColor txtColor = theme->spectrumChartText.isDefault()
  ? WColor(GlobalColor::black)
  : theme->spectrumChartText;
  
  WPen txtpen( txtColor );
  setTextPen( txtpen );
  axis(Chart::XAxis).setTextPen( txtpen );
  axis(Chart::YAxis).setTextPen( txtpen );
  axis(Chart::Y2Axis).setTextPen( txtpen );
  
  setCrosshairColor( txtColor );
  
  if( theme->spectrumChartBackground.isDefault() )
    setBackground( Wt::NoBrush );
  else
    setBackground( WBrush(theme->spectrumChartBackground) );
  
  
  //From what I can tell, we cant change the legend text color easily, so
  //  we'll just cheat and back the legend background different enough from
  //  black so we can always read the text.  Lame, but whatever.
  setLegendStyle( legendFont(), legendBorder(), WBrush(Wt::WColor(220, 220, 220, 120)) );
  
  if( (theme->spectrumChartMargins.isDefault() && !theme->spectrumChartBackground.isDefault()) )
  {
    //theme->spectrumChartBackground
  }else if( !theme->spectrumChartMargins.isDefault() )
  {
    //theme->spectrumChartMargins
  }
  
  const WColor axisLineColor = theme->spectrumAxisLines.isDefault()
  ? WColor(GlobalColor::black)
  : theme->spectrumAxisLines;
  
  WPen defpen = axis(Chart::XAxis).pen();
  defpen.setColor( axisLineColor );
  axis(Chart::XAxis).setPen( defpen );
  axis(Chart::Y1Axis).setPen( defpen );
  axis(Chart::Y2Axis).setPen( defpen );
}//handleColorThemeChange(...)


void DrfChart::updateChart( std::shared_ptr<const DetectorPeakResponse> det )
{
  m_detector = det;
  
  // clear series if any
  removeSeries(1);
  removeSeries(2);
  
  const Wt::Chart::SeriesType seriesType = Wt::Chart::SeriesType::CurveSeries;
  
  const bool hasEfficiency = !!m_detector && m_detector->isValid();
  
  if( hasEfficiency )
  {
    int nValidPoints = 0, nEffPoints = 0;
    try
    {
      const bool hasResloution = m_detector->hasResolutionInfo();
      
      float minEnergy = 0.0f, maxEnergy = 3000.0f;
      
      // If DRF has a defined energy range, of at least 100 keV, use it.
      if( m_detector && m_detector->isValid()
         && (m_detector->upperEnergy() > (m_detector->lowerEnergy() + 100.0)) )
      {
        minEnergy = m_detector->lowerEnergy();
        maxEnergy = m_detector->upperEnergy();
      }
      
      
      // TODO: pick this better using the chart width or whatever
      int numEnergyPoints = floor(maxEnergy - minEnergy) / 2.5; //max of 8k points, but usually 12.5
      if( numEnergyPoints > 4500 ) //4500 chosen arbitrarily
        numEnergyPoints = 4500;
      
      m_efficiencyModel->clear();
      m_efficiencyModel->insertRows( 0, numEnergyPoints );
      m_efficiencyModel->insertColumns( 0, hasResloution? 3 : 2 );
      float energy = 0.0f, efficiency = 0.0f;
      for( int row = 0; row < numEnergyPoints; ++row )
      {
        energy = minEnergy + (float(row)/float(numEnergyPoints)) * (maxEnergy-minEnergy);
        efficiency = static_cast<float>( m_detector->intrinsicEfficiency( energy ) );
        m_efficiencyModel->setData(row, 0, energy );
        
        
        //Skip any points outside where we would expect.
        if( IsNan(efficiency) || IsInf(efficiency) || efficiency < 0.0f )
        {
          m_efficiencyModel->setData(row, 1, boost::any() );
        }else
        {
          ++nValidPoints;
          m_efficiencyModel->setData(row, 1, efficiency );
        }
        
        if( hasResloution )
        {
          const float fwhm = m_detector->peakResolutionFWHM(energy);
          if( !IsNan(fwhm) && !IsInf(fwhm) && fwhm>=0.0f && fwhm<9999.9f )
          {
            ++nEffPoints;
            m_efficiencyModel->setData(row, 2, fwhm);
          }
        }
      }//for( int row = 0; row < numEnergyPoints; ++row )
      
      m_efficiencyModel->setHeaderData(0, Wt::WString("Energy"));
      m_efficiencyModel->setHeaderData(1, Wt::WString("Efficiency"));
      if( m_detector->hasResolutionInfo() )
        m_efficiencyModel->setHeaderData(2, Wt::WString("FWHM"));
      
      setXSeriesColumn(0);  //Not having this line after creating a new model was why the x-axis was the wrong scale
      Wt::Chart::WDataSeries s1(1, seriesType, Wt::Chart::Y1Axis);
      s1.setPen( WPen(m_chartEnergyLineColor) );
      addSeries(s1);
      
#if( WT_VERSION >= 0x3040000 )
      setFollowCurve( 1 );
#endif
      
      if( nValidPoints == 0 )
      {
        axis(Wt::Chart::Y1Axis).setRange( 0.0, 1.0 );
      }else
      {
        //m_chart->axis(Wt::Chart::Y1Axis).setRoundLimits(<#WFlags<Wt::Chart::AxisValue> locations#>)
        axis(Wt::Chart::Y1Axis).setAutoLimits( Chart::MinimumValue | Chart::MaximumValue );
      }
      
      if( nEffPoints )
      { //only if there is resolution FWHM
        Wt::Chart::WDataSeries s2(2, seriesType, Wt::Chart::Y2Axis);
        s2.setPen( WPen(m_chartFwhmLineColor) );
        addSeries(s2);
        axis(Wt::Chart::Y2Axis).setTitle("FWHM");
        axis(Wt::Chart::Y2Axis).setVisible(true);
      }else
      { //no ResolutionInfo
        axis(Wt::Chart::Y2Axis).setVisible(false);
        axis(Wt::Chart::Y2Axis).setTitle("");
      }//no ResolutionInfo
    }catch( std::exception &e )
    {
      cerr << "DrfSelect::updateChart()\n\tCaught: " << e.what() << endl;
    }//try / catch
  }//if( hasEfficiency )
} //DrfChart::updateChart()

