#ifndef RelActAutoGui_h
#define RelActAutoGui_h
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

#include <atomic>
#include <memory>
#include <utility>

#include <Wt/WContainerWidget>

#include "InterSpec/RelActCalcAuto.h"


class PeakDef;
class PeakModel;
class AuxWindow;
class InterSpec;
class RelEffChart;
class PopupDivMenu;
class PopupDivMenuItem;
class RelActTxtResults;
class RelEffShieldWidget;
class DetectorPeakResponse;
class D3SpectrumDisplayDiv;
class RelActAutoGuiNuclide;
class RelActAutoGuiRelEffOptions;

namespace SpecUtils
{
  class Measurement;
  enum class SpectrumType : int;
}//namespace SpecUtils

namespace RelActCalcAuto
{
  struct RelActAutoSolution;
}//namespace RelActCalcAuto

namespace Wt
{
  class WMenu;
  class WText;
  class WCheckBox;
  class WComboBox;
  class WInPlaceEdit;
}//namespace Wt

namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml



class RelActAutoGui : public Wt::WContainerWidget
{
public:
  enum class AddUncert : int
  {
    StatOnly, OneHundrethPercent, OneTenthPercent, OnePercent, FivePercent, TenPercent, 
    TwentyFivePercent, FiftyPercent, SeventyFivePercent, OneHundredPercent, NumAddUncert
  };//enum class AddUncert
  
  static const char *to_str( const AddUncert val );
  
public:
  RelActAutoGui( InterSpec *viewer, Wt::WContainerWidget *parent = nullptr );
  
  ~RelActAutoGui();
  
  static std::pair<RelActAutoGui *,AuxWindow *> createWindow( InterSpec *viewer  );
  
  void updateDuringRenderForSpectrumChange();
  void updateSpectrumToDefaultEnergyRange();
  void updateDuringRenderForNuclideChange();
  void updateDuringRenderForRefGammaLineChange();
  void updateDuringRenderForEnergyRangeChange();
  void updateDuringRenderForFreePeakChange();
  void startUpdatingCalculation();
  void updateFromCalc( std::shared_ptr<RelActCalcAuto::RelActAutoSolution> answer,
                      std::shared_ptr<std::atomic_bool> cancel_flag );
  void handleCalcException( std::shared_ptr<std::string> message,
                           std::shared_ptr<std::atomic_bool> cancel_flag );
  
  
  void handleDisplayedSpectrumChange( SpecUtils::SpectrumType );
  void handlePresetChange();
  void handleRelEffEqnFormChanged();
  void handleRelEffEqnOrderChanged();
  void handleFwhmFormChanged();
  void handleFwhmEstimationMethodChanged();
  void handleFitEnergyCalChanged();
  void handleBackgroundSubtractChanged();
  void handleSameAgeChanged();
  void handlePuByCorrelationChanged();
  void handleSkewTypeChanged();
  void handleNucDataSrcChanged();
  void handleAddNuclide();
  void handleAddEnergy();
  void handleClearAllEnergyRanges();
  void handleShowFreePeaks();
  void handleHideFreePeaks();
  void handleAddFreePeak( const double energy,
                          const bool constrain_fwhm,
                          const bool apply_energy_cal );
  void handleRemoveEnergy( Wt::WWidget *w );
  void handleRemoveNuclide( Wt::WWidget *w );
  void handleRemoveFreePeak( Wt::WWidget *w );
  void handleRemovePartOfEnergyRange( Wt::WWidget *energy_range,
                                      double lower_energy,
                                     double upper_energy );
  void handleSplitEnergyRange( Wt::WWidget *energy_range, const double energy );
  
  void handleConvertEnergyRangeToIndividuals( Wt::WWidget *energy_range );
  

  /** Called when a nuclide is added or removed (or changed from one to another) */
  void handleNuclidesChanged();

  /** Called when a nuclide's fit age is changed. */
  void handleNuclideFitAgeChanged( RelActAutoGuiNuclide *nuc, bool fit_age );

  void handleNuclideAgeChanged( RelActAutoGuiNuclide *nuc );
  
  /** Called when energy ranges are added, removed, or edited. */
  void handleEnergyRangeChange();
  
  /** Called when free peaks are added, removed, or edited. */
  void handleFreePeakChange();
  
  void handleAdditionalUncertChanged();
  
  void setOptionsForNoSolution();
  void setOptionsForValidSolution();
  void makeZeroAmplitudeRoisToChart();
  
  /** Checks if the m_presets is in a "custom" state, and if not, puts it there
   
   @param force_create if true, will create a new state, even if in a user created/modified state,
          if false, will only create a new state if not already in a custom state.
   
   TODO: save the previous user to allow going back to them
   */
  void checkIfInUserConfigOrCreateOne( const bool force_create );
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags ) override;
  
  ::rapidxml::xml_node<char> *serialize( ::rapidxml::xml_node<char> *parent ) const;
  void deSerialize( const rapidxml::xml_node<char> *base_node );
  
  std::unique_ptr<rapidxml::xml_document<char>> guiStateToXml() const;
  void setGuiStateFromXml( const rapidxml::xml_document<char> *doc );
  
  void setCalcOptionsGui( const RelActCalcAuto::Options &options );
  

  RelActCalcAuto::Options getCalcOptions() const;
  std::vector<RelActCalcAuto::NucInputInfo> getNucInputInfo( const int rel_eff_curve_index ) const;
  std::vector<RelActCalcAuto::RelEffCurveInput::ActRatioConstraint> getActRatioConstraints( const int rel_eff_curve_index ) const;
  std::vector<RelActCalcAuto::RelEffCurveInput::MassFractionConstraint> getMassFractionConstraints( const int rel_eff_curve_index ) const;

  std::vector<RelActCalcAuto::RoiRange> getRoiRanges() const;
  std::vector<RelActCalcAuto::FloatingPeak> getFloatingPeaks() const;
  
  
  std::shared_ptr<const RelActCalcAuto::RelActAutoSolution> getCurrentSolution() const;
  
protected:
  void handleRoiDrag( double new_roi_lower_energy,
                     double new_roi_upper_energy,
                     double new_roi_lower_px,
                     double new_roi_upper_px,
                     const double original_roi_lower_energy,
                     const bool is_final_range );
  
  void handleCreateRoiDrag( const double lower_energy,
                            const double upper_energy,
                            const int num_peaks_to_force,
                           const bool is_final_range );
  
  void handleShiftDrag( const double lower_energy, const double upper_energy );
  void handleDoubleLeftClick( const double energy, const double counts );
  void handleRightClick( const double energy, const double counts,
                        const int page_x_px, const int page_y_px );
  
  void handleToggleForceFullRange( Wt::WWidget *w );
  
  /** Combines the two RelActAutoGuiEnergyRange's into a single range, and returns the resulting RelActAutoGuiEnergyRange.
   
   Will only return nullptr if one/both of the inputs is not a RelActAutoGuiEnergyRange, or not in #m_energy_ranges.
   */
  Wt::WWidget *handleCombineRoi( Wt::WWidget *left_roi, Wt::WWidget *right_roi );
  
  void removeAllEnergyRanges();
  void startApplyFitEnergyCalToSpecFile();
  void applyFitEnergyCalToSpecFile();
  void handleShowRefLines( const bool show );
  void setPeaksToForeground();
  
  void handleRelEffModelOptionsChanged();
  void handleDetectorChange();

  void handleAddRelEffCurve();

  void handleDelRelEffCurve( RelActAutoGuiRelEffOptions *curve );
  void handleRelEffCurveNameChanged( RelActAutoGuiRelEffOptions *curve, const Wt::WString &name );

  protected:
  RelActAutoGuiRelEffOptions *getRelEffCurveOptions( const int index );
  const RelActAutoGuiRelEffOptions *getRelEffCurveOptions( const int index ) const;

  /** Calculation has been started. */
  Wt::Signal<> &calculationStarted();
  
  /** Calculation finished and the solution is valid. */
  Wt::Signal<> &calculationSuccessful();
  
  /** Calculation failed, and m_solution is nullptr. */
  Wt::Signal<> &calculationFailed();
  
  /** Calculation is finished. */
  Wt::Signal< std::shared_ptr<const RelActCalcAuto::RelActAutoSolution> > &solutionUpdated();
  
  void addDownloadAndUploadLinks( Wt::WContainerWidget *parent );
  void handleRequestToUploadXmlConfig();
protected:
  
  enum RenderActions
  {
    UpdateSpectra         = 0x01,
    UpdateNuclidesPresent = 0x02,
    UpdateEnergyRanges    = 0x04,
    UpdateCalculations    = 0x08,
    ChartToDefaultRange   = 0x10,
    UpdateFreePeaks       = 0x20,
    UpdateFitEnergyCal    = 0x40,
    UpdateRefGammaLines   = 0x80
  };//enum D3RenderActions
  
  Wt::WFlags<RenderActions> m_render_flags;
  
  std::string m_default_par_sets_dir;
  std::string m_user_par_sets_dir;
  
  InterSpec *m_interspec;
  std::shared_ptr<const SpecUtils::Measurement> m_foreground;
  std::shared_ptr<const SpecUtils::Measurement> m_background;
  double m_background_sf;
  
  
  Wt::WText *m_status_indicator;
  
  D3SpectrumDisplayDiv *m_spectrum;
  PeakModel *m_peak_model;
  RelEffChart *m_rel_eff_chart;
  RelActTxtResults *m_txt_results;
  Wt::WMenu *m_upper_menu;
  
  /** The paths to XML files representing user-selectable states; both ones that come with
   InterSpec, and ones in the users data directory..  Initialized during widget construction,
   and the very first entry will be blank.
   */
  std::vector<std::string> m_preset_paths;
  
  /** A selecter for user states (i.e. states from XML) - currently the "preset" name is out
   of date and a new name for gui state needs to be decided on.
   */
  Wt::WComboBox *m_presets;
  
  /** A flag to indicate if we are currently loading state from XML; it is set to true in
   #deSerialize, and then false at the end of render.
   #checkIfInUserConfigOrCreateOne will check this flag, and if true, it will not create
   a new config.
   */
  bool m_loading_preset;
  
  /** The currently selected preset index; we keep this around to compare the new index
   to when #m_preset gets changed.
   */
  int m_current_preset_index;
  
  /** When the user selects a new state, we'll save thier old ones in this variable, so
   they can go back to it. The map index corresponds to the #m_preset index for the state.
   */
  std::map<int,std::unique_ptr<rapidxml::xml_document<char>>> m_previous_presets;

  /** The place to indicate errors in calculation, when calc is not successful. */
  Wt::WText *m_error_msg;
  /* The place to give a summary of the fit, so like "chi2 = 123", or "Chi2 = 123, Uranium Enrichment 23.2%" */
  Wt::WText *m_fit_chi2_msg;

  Wt::WMenu *m_rel_eff_opts_menu;
  Wt::WStackedWidget *m_rel_eff_opts_stack;

  Wt::WComboBox *m_fwhm_eqn_form;
  Wt::WComboBox *m_fwhm_estimation_method;
  
  Wt::WCheckBox *m_fit_energy_cal;
  Wt::WCheckBox *m_background_subtract;
  Wt::WCheckBox *m_same_z_age;
  
  
  Wt::WComboBox *m_skew_type;
  Wt::WComboBox *m_add_uncert;
  
  // Wt::WComboBox *m_u_pu_data_source;
  PopupDivMenu *m_more_options_menu;
  PopupDivMenuItem *m_apply_energy_cal_item;
  PopupDivMenuItem *m_show_ref_lines_item;
  PopupDivMenuItem *m_hide_ref_lines_item;
  PopupDivMenuItem *m_set_peaks_foreground;
  
  /** If the user wants to show reference gamma lines, we'll use a #ReferencePhotopeakDisplay
   widget to calculate them and load them to m_spectrum; this is primarily for code re-use
   until I bother to refactor #ReferencePhotopeakDisplay to more easily generate reference
   lines.
   This #ReferencePhotopeakDisplay wont be inserted into the Wt widget hierarchy.
   */
  std::unique_ptr<ReferencePhotopeakDisplay> m_photopeak_widget;

  
  Wt::WPushButton *m_clear_energy_ranges;
  Wt::WPushButton *m_show_free_peak;
  Wt::WContainerWidget *m_free_peaks_container;
  
  Wt::WContainerWidget *m_nuclides;
  Wt::WContainerWidget *m_energy_ranges;
  Wt::WContainerWidget *m_free_peaks;
  
  // For the future
  // - Free peaks
  // - Load peaks to foreground
  // - Additional uncertainty
  // - Export to XML setup
  // - HTML Report
  
  
  bool m_is_calculating;
  std::shared_ptr<std::atomic_bool> m_cancel_calc;
  std::shared_ptr<RelActCalcAuto::RelActAutoSolution> m_solution;
  
  /** A good amount of calculation time is spent determining all the "auto-searched" peaks
   and consequently the FWHM DRF, so we will cache these so we only compute them
   */
  std::shared_ptr<const DetectorPeakResponse> m_cached_drf;
  std::vector<std::shared_ptr<const PeakDef>> m_cached_all_peaks;
  
  /** Emitted when a new calculation is started, even if there was already one running. */
  Wt::Signal<> m_calc_started;
  
  /** Calculation finished and the solution is valid. */
  Wt::Signal<> m_calc_successful;
  
  /** Emitted when calculation failed (m_solution is nullptr), or failed to find solution. */
  Wt::Signal<> m_calc_failed;
  
  /** Emitted whenever m_solution is set. */
  Wt::Signal< std::shared_ptr<const RelActCalcAuto::RelActAutoSolution> > m_solution_updated;
  
  Wt::WResource *m_html_download_rsc;
  
  Wt::WResource *m_xml_download_rsc;
};//class RelActAutoGui


#endif //RelActAutoGui_h
