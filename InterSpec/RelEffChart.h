#ifndef RelEffChart_h
#define RelEffChart_h

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


class RelEffChart : public Wt::WContainerWidget
{
public:
  RelEffChart( Wt::WContainerWidget *parent = 0 );
  virtual ~RelEffChart();
  
  /** Set data from an "auto" relative efficiency fit. */
  void setData( const double live_time,
               const std::vector<PeakDef> &fit_peaks,
               const std::vector<RelActCalcAuto::NuclideRelAct> &rel_acts,
               const std::string &jsRelEffEqn );
  
  /** Set data from an "manual" relative efficiency fit. */
  void setData( const std::vector<RelActCalcManual::GenericPeakInfo> &peaks,
               const std::map<std::string,std::pair<double,std::string>> &relActsColors,
               std::string relEffEqn );
  
  void setLineColor( const Wt::WColor &color );
  void setTextColor( const Wt::WColor &color );
  void setAxisLineColor( const Wt::WColor &color );
  void setChartBackgroundColor( const Wt::WColor &color );
  void setDefaultMarkerColor( const Wt::WColor &color );
  
protected:
  void defineJavaScript();
  
  void setCssRules();
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  /** The javascript variable name used to refer to the RelEffChart object.
   Currently is `jsRef() + ".chart"`.
   */
  const std::string m_jsgraph;
  
  std::map<std::string,Wt::WCssTextRule *> m_cssRules;
  
  /** JS calls requested before the widget has been rendered, so wouldnt have
     ended up doing anything are saved here, and then executed once the widget
     is rendered.
     Note that not all calls to the D3 chart before Wt's rendering need to go
     here as they will be options set to the D3 chart during first rendering.
   */
  std::vector<std::string> m_pendingJs;

};//class RelEffChart


#endif //RelEffChart_h
