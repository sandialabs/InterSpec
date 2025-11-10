#ifndef ShieldingSourceFitPlot_h
#define ShieldingSourceFitPlot_h

#include "InterSpec_config.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Wt/WContainerWidget>
#include <Wt/WString>

//Forward declarations
namespace Wt
{
  class WColor;
  class WCssTextRule;
}//namespace Wt

namespace ShieldingSourceFitCalc
{
  struct ModelFitResults;
}//namespace ShieldingSourceFitCalc


class ShieldingSourceFitPlot : public Wt::WContainerWidget
{
public:
  ShieldingSourceFitPlot( Wt::WContainerWidget *parent = 0 );
  virtual ~ShieldingSourceFitPlot();

  /** Set data from activity/shielding fit results. */
  void setData( const ShieldingSourceFitCalc::ModelFitResults &results );

  /** Set whether to show Chi or Mult mode. */
  void setShowChi( bool show_chi );

  /** Generate JSON string from ModelFitResults for passing to JavaScript or static HTML. */
  static std::string jsonForData( const ShieldingSourceFitCalc::ModelFitResults &results );

  /** Set the x-axis title. */
  void setXAxisTitle( const Wt::WString &title );

  /** Set the y-axis title for Chi mode. */
  void setYAxisTitleChi( const Wt::WString &title );

  /** Set the y-axis title for Mult mode. */
  void setYAxisTitleMult( const Wt::WString &title );

  /** Set the chart content margins (e.g. how many pixels inside the <svg /> element the titles or axises should be drawn).
   */
  void setContentMargins( int top, int right, int bottom, int left );

  /** Signal emitted when user toggles between Chi and Mult display modes.
      Argument is the new state of showChi (true = Chi mode, false = Mult mode).
   */
  Wt::Signal<bool> &displayModeChanged();

  /** Signal emitted when user clicks on a data point.
      Argument is the energy of the clicked point.
   */
  Wt::Signal<double> &dataPointClicked();

protected:
  void defineJavaScript();

  void setCssRules();

  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );

  virtual void refresh();

  /** Generate JSON object with localized strings for JavaScript. */
  std::string localizedStringsJson();

  /** Send data JSON to JavaScript chart. */
  void sendDataToJavaScript( const std::string &jsonData );

  /** Callback for when display mode is changed in JavaScript */
  void handleDisplayModeChanged( bool showChi );

  /** Callback for when a data point is clicked in JavaScript */
  void handleDataPointClicked( double energy );

private:
  /** The javascript variable name used to refer to the ShieldingSourceFitPlot object.
   Currently is `jsRef() + ".chart"`.
   */
  const std::string m_jsgraph;

  Wt::WString m_xAxisTitle;
  Wt::WString m_yAxisTitleChi;
  Wt::WString m_yAxisTitleMult;

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

  /** Whether to show Chi mode (true) or Mult mode (false). */
  bool m_showChi;

  std::map<std::string,Wt::WCssTextRule *> m_cssRules;

  /** JS calls requested before the widget has been rendered, so wouldnt have
     ended up doing anything are saved here, and then executed once the widget
     is rendered.
     Note that not all calls to the D3 chart before Wt's rendering need to go
     here as they will be options set to the D3 chart during first rendering.
   */
  std::vector<std::string> m_pendingJs;

  /** JSignal for display mode changes from JavaScript. */
  std::unique_ptr<Wt::JSignal<bool>> m_displayModeChangedJS;

  /** JSignal for data point clicks from JavaScript. */
  std::unique_ptr<Wt::JSignal<double>> m_dataPointClickedJS;

  /** Wt::Signal for display mode changes. */
  Wt::Signal<bool> m_displayModeChanged;

  /** Wt::Signal for data point clicks. */
  Wt::Signal<double> m_dataPointClicked;

};//class ShieldingSourceFitPlot


#endif //ShieldingSourceFitPlot_h
