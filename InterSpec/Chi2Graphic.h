#ifndef Chi2Graphic_h
#define Chi2Graphic_h

#include "InterSpec_config.h"

#include <map>
#include <vector>

#include <Wt/WContainerWidget>

#include "InterSpec/PeakDef.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcManual.h"

//Forward declarations
namespace Wt
{
  class WColor;
  class WCssTextRule;
}//namespace Wt



class Chi2Graphic : public Wt::WContainerWidget
{
public:
  
  struct PeakFitInfo
  {
    /** Energy of the gamma associated with the peak. */
    double energy;
    
    /** `(observed_counts - expected_counts) / observed_uncertainty` */
    double numSigmaOff;
    
    /** `observed_counts / expected_counts` */
    double observedOverExpected;
    double observedOverExpectedUncert;
    
    /** `peak.lineColor()` */
    Wt::WColor peakColor;
    
    Wt::WString nuclidename;
  };//struct PeakFitInfo
  
public:
  Chi2Graphic( Wt::WContainerWidget *parent = 0 );
  virtual ~Chi2Graphic();
  
  /** Enumation for what to display. */
  enum class ChiDispType
  {
    /** The dots showing how dar away each point is away from the model answer.. */
    Pull,
    
    /** Dots showing the multiple each peak is from the models prediction. */
    Scale
  };//enum class ChiDispType
  
  /** Sets how the data should be displayed. */
  void setChartType( const ChiDispType type );
  
  void setData( const int ndof, const std::vector<PeakFitInfo> &used_points );
  
  void setXAxisTitle( const Wt::WString &title );
  void setYAxisTitle( const Wt::WString &title );
  
  /** Set the chart content margins (e.g. how many pixels inside the <svg /> element the titles or axises should be drawn).
   */
  void setContentMargins( int top, int right, int bottom, int left );
  
protected:
  void defineJavaScript();
  
  void setDataToClient();
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  /** The javascript variable name used to refer to the RelEffChart object.
   Currently is `jsRef() + ".chart"`.
   */
  const std::string m_jsgraph;
  
  Wt::WString m_xAxisTitle;
  Wt::WString m_yAxisTitle;
  
  ChiDispType m_displayType;
  
  /** Margins between the edge of the SVG, and where the plot contents start getting drawn, in pixels.
   It would be reasonable to instead implement this as CSS margin on the <div /> object the SVG is in,
   but leaving for the moment.
   */
  int m_topMargin;
  int m_rightMargin;
  int m_bottomMargin;
  int m_leftMargin;
  
  /** The distance between the axis title (if present), and the axis numbers.  Defaults to -3 to make compact. */
  int m_titlePadding;
  
  /** JS calls requested before the widget has been rendered, so wouldnt have
     ended up doing anything are saved here, and then executed once the widget
     is rendered.
     Note that not all calls to the D3 chart before Wt's rendering need to go
     here as they will be options set to the D3 chart during first rendering.
   */
  std::vector<std::string> m_pendingJs;

  enum Chi2RenderActions
  {
    SetData = 0x01
  };
  
  Wt::WFlags<Chi2RenderActions> m_renderFlags;
  
  int m_ndof;
  std::vector<PeakFitInfo> m_used_points;
};//class Chi2Graphic


#endif //Chi2Graphic_h
