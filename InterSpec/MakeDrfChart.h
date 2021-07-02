#ifndef MakeDrfChart_h
#define MakeDrfChart_h
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

#include <string>
#include <vector>

#include <Wt/WColor>
#include <Wt/Chart/WCartesianChart>


class MakeDrfChart : public Wt::Chart::WCartesianChart
{
public:
  
  /** Units efficiency and FWHM equation coefficients are given in. */
  enum class EqnEnergyUnits{ keV, MeV };
  
  /** Fwhm equation type; equivalent to DetectorPeakResponse::ResolutionFnctForm
   but redifined here to not have to include that header.
   */
  enum class FwhmCoefType{ Gadras, SqrtEqn };
  
  /** Information from the calibration data for a specific peak.  Used to draw
     data points on the chart.
   */
  struct DataPoint
  {
    /** Energy, in keV, of peak. */
    float energy;
    
    /** Live time, in seconds, of spectrum peak belongs to. */
    float livetime;
    
    /** Fit counts in the peak; does not include continuum area. */
    float peak_area;
    
    /** Uncertainty in peak counts; e.g. usually ~sqrt(counts_in_peak).
        If you background subtracted, make sure to include background error to.
     */
    float peak_area_uncertainty;
    
    /** Full width at half maximum (eg 2.35*gaus_sigma), in keV, of the peak. */
    float peak_fwhm;
    
    /** Uncertainty, in keV, of FWHM.  Zero or less means dont use. */
    float peak_fwhm_uncertainty;
    
    /** The gammas emmitted into 4pi by the source at the energy (or at least
        that contribute to the peak).
     */
    float source_count_rate;
    
    /** Source assay uncert, in CPS. Zero or less means dont use. */
    float source_count_rate_uncertainty;
    
    /** Distance from source to detector face. */
    float distance;
    
    /** Tooltip for data point about the source, ex. "22 uCi Cs137 @ 1.2m" */
    std::string source_information;
    
    /** Color of the peak for this data point. */
    Wt::WColor peak_color;
    
    /** If zero or negative, then background subtraction was not done.
       Otherwise the background peak area is indicated here.
     */
    float background_peak_area;
    
    /** Must be >0 for background peak subtraction to be done. */
    float background_peak_live_time;
  };//struct DataPoint
  
public:
  MakeDrfChart( Wt::WContainerWidget *parent = nullptr );
  virtual ~MakeDrfChart();
  virtual void paint( Wt::WPainter &painter,
                     const Wt::WRectF &rectangle = Wt::WRectF() ) const;
  virtual void paintEvent( Wt::WPaintDevice *paintDevice );
  
  virtual void layoutSizeChanged( int width, int height );
  
  /** Set the coefficients to draw the FWHM response with.
      Setting with zero coefficients will remove line from chart.
  */
  void setFwhmCoefficients( const std::vector<float> &coeffs,
                            const std::vector<float> &uncerts,
                            const FwhmCoefType eqnType,
                            const EqnEnergyUnits units  );
  
  
  /** Set the coefficients to draw the detector efficiency response with.
      Setting with zero coefficients will remove line from chart.
   */
  void setEfficiencyCoefficients( const std::vector<float> &coeffs,
                                  const std::vector<float> &uncerts,
                                  const EqnEnergyUnits units );
  
  /**
   */
  void setDataPoints( const std::vector<DataPoint> &datapoints,
                      const float det_diameter,
                      const float det_low_energy, const float det_up_energy );
  
  /** Returns the current data points. */
  const std::vector<DataPoint> &currentDataPoints() const;
  
  
  /** Sets whther the FWHM data points and equation should be shown. */
  void showFwhmPoints( const bool show );
  
  /** Set the displayed x-range, in keV. */
  void setXRange( double lower, double upper );
  
  /** Emitted when the x-range changes due to data; not emitted when setXRange
     called.
   */
  Wt::Signal<double,double> &xRangeChanged();
  
protected:
  void updateYAxisRange();
  void updateDataToModel();
  void updateEqnEnergyToModel();
  void updateEffEquationToModel();
  void updateFwhmEquationToModel();
  
  void updateColorTheme( std::shared_ptr<const ColorTheme> theme );

  float m_det_diameter;
  float m_det_lower_energy;
  float m_det_upper_energy;
  std::vector<DataPoint> m_datapoints;  //probably dont need this, but whtaever for now
  
  /** Energy units equation for FWHM is given in.  */
  EqnEnergyUnits m_fwhmEnergyUnits;
  
  /** We are currently using a GADRAS style FWHM equation with exactly three
      coefficients, so this vector should have either zero or three entries.
   */
  std::vector<float> m_fwhmCoefs;
  
  std::vector<float> m_fwhmCoefUncerts;
  
  /** Energy units equation for detector efficiency is given in.  */
  EqnEnergyUnits m_efficiencyEnergyUnits;
  
  /** Coefficients for intrinsic efficiency formula.
      E.x., exp( a + log(energy)^1 + log(energy)^2 + ...
   */
  std::vector<float> m_efficiencyCoefs;
  
  std::vector<float> m_efficiencyCoefUncerts;
  
  FwhmCoefType m_fwhmEqnType;
  
  Wt::Signal<double,double> m_xRangeChanged;
  
  Wt::WBrush m_chartMarginBrush;
  Wt::WPen m_textPen;
};//class WCartesianChart


#endif  //MakeDrfChart_h
