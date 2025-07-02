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
  
  struct ReCurveInfo
  {
    double live_time = 0.0;
    std::vector<PeakDef> fit_peaks;
    std::vector<RelActCalcAuto::NuclideRelAct> rel_acts;
    std::string js_rel_eff_eqn;
    std::string js_rel_eff_uncert_eqn;
    Wt::WString re_curve_name;
    Wt::WString re_curve_eqn_txt;
  };//struct ReCurveInfo
  
  /** Set data from an "auto" relative efficiency fit. */
  void setData( const ReCurveInfo &info );
  
  /** Set data from multiple "auto" relative efficiency fits. */
  void setData( const std::vector<ReCurveInfo> &infoSets );
  
  /** Set data from an "manual" relative efficiency fit. */
  void setData( const std::vector<RelActCalcManual::GenericPeakInfo> &peaks,
               const std::map<std::string,std::pair<double,std::string>> &nuc_to_act_and_color,
               std::string relEffEqn,
               const Wt::WString &chi2_title_str,
               const std::string &jsRelEffEqnUncert );
  
  /** Generate JSON from ReCurveInfo objects */
  static std::string jsonForData(const std::vector<ReCurveInfo> &infoSets);
  
  void setLineColor( const Wt::WColor &color );
  void setTextColor( const Wt::WColor &color );
  void setAxisLineColor( const Wt::WColor &color );
  void setChartBackgroundColor( const Wt::WColor &color );
  void setDefaultMarkerColor( const Wt::WColor &color );
  
  /** Set custom colors for datasets based on ColorTheme */
  void setRelEffCurveColors();
  
  void setXAxisTitle( const std::string &title );
  void setYAxisTitle( const std::string &title );
  
  /** Set the chart content margins (e.g. how many pixels inside the <svg /> element the titles or axises should be drawn).
   */
  void setContentMargins( int top, int right, int bottom, int left );
  
protected:
  void defineJavaScript();
  
  void setCssRules();
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
private:
  /** Helper struct for dataset information */
  struct RelEffChartDataset
  {
    std::vector<RelActCalcManual::GenericPeakInfo> peaks;
    std::map<std::string, std::pair<double, std::string>> relActsColors;
    std::string relEffEqn;
    Wt::WString chi2_title_str;
    std::string relEffEqnUncert;
  };
  
  /** Generate JSON for multiple datasets */
  static std::string jsonForData(const std::vector<RelEffChartDataset> &datasets);
  
  /** Generate JSON for a single dataset */
  static std::string jsonForDataset(const RelEffChartDataset &dataset, bool isFirstDataset);
  
  /** Send JavaScript to set the datasets data */
  void sendDataToJavaScript(const std::string &jsonData);
  
  /** Private implementation for multiple datasets */
  void setData(const std::vector<RelEffChartDataset> &datasets);
  
  /** The javascript variable name used to refer to the RelEffChart object.
   Currently is `jsRef() + ".chart"`.
   */
  const std::string m_jsgraph;
  
  std::string m_xAxisTitle;
  std::string m_yAxisTitle;
  
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
