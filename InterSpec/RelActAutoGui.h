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
class RelActTxtResults;
class DetectorPeakResponse;
class D3SpectrumDisplayDiv;

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
}//namespace Wt

namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml



class RelActAutoGui : public Wt::WContainerWidget
{
public:
  RelActAutoGui( InterSpec *viewer, Wt::WContainerWidget *parent = nullptr );
  
  static std::pair<RelActAutoGui *,AuxWindow *> createWindow( InterSpec *viewer  );
  
  void updateDuringRenderForSpectrumChange();
  void updateSpectrumToDefaultEnergyRange();
  void updateDuringRenderForNuclideChange();
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
  void handleFitEnergyCalChanged();
  void handleBackgroundSubtractChanged();
  void handleSameAgeChanged();
  void handleUPuByCorrelationChanged();
  void handleNucDataSrcChanged();
  void handleAddNuclide();
  void handleAddEnergy();
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
  
  /** Called when a nuclides information (such as age) is changed. */
  void handleNuclidesInfoEdited();
  
  /** Called when energy ranges are added, removed, or edited. */
  void handleEnergyRangeChange();
  
  /** Called when free peaks are added, removed, or edited. */
  void handleFreePeakChange();
  
  void setOptionsForNoSolution();
  void setOptionsForValidSolution();
  void makeZeroAmplitudeRoisToChart();
  
  /** Checks if the m_presets is in a "custom" state, and if not, puts it there
   
   TODO: save the previous user to allow going back to them
   */
  void checkIfInUserConfigOrCreateOne();
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags ) override;
  
  ::rapidxml::xml_node<char> *serialize( ::rapidxml::xml_node<char> *parent ) const;
  void deSerialize( const rapidxml::xml_node<char> *base_node );
  
  std::unique_ptr<rapidxml::xml_document<char>> guiStateToXml() const;
  void setGuiStateFromXml( const rapidxml::xml_document<char> *doc );
  
  void setCalcOptionsGui( const RelActCalcAuto::Options &options );
  
  RelActCalcAuto::Options getCalcOptions() const;
  std::vector<RelActCalcAuto::NucInputInfo> getNucInputInfo() const;
  std::vector<RelActCalcAuto::RoiRange> getRoiRanges() const;
  std::vector<RelActCalcAuto::FloatingPeak> getFloatingPeaks() const;
  
  
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
  
  void handleRightClick( const double energy, const double counts,
                        const int page_x_px, const int page_y_px );
  
  void handleToggleForceFullRange( Wt::WWidget *w );
  void handleCombineRoi( Wt::WWidget *left_roi, Wt::WWidget *right_roi );
  
  
protected:
  
  enum RenderActions
  {
    UpdateSpectra         = 0x01,
    UpdateNuclidesPresent = 0x02,
    UpdateEnergyRanges    = 0x04,
    UpdateCalculations    = 0x08,
    ChartToDefaultRange   = 0x10,
    UpdateFreePeaks       = 0x20,
    UpdateFitEnergyCal    = 0x40
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
  
  Wt::WText *m_error_msg;
  
  Wt::WComboBox *m_rel_eff_eqn_form;
  Wt::WComboBox *m_rel_eff_eqn_order;
  
  Wt::WComboBox *m_fwhm_eqn_form;
  
  Wt::WCheckBox *m_fit_energy_cal;
  Wt::WCheckBox *m_background_subtract;
  
  Wt::WCheckBox *m_same_z_age;
  Wt::WCheckBox *m_u_pu_by_correlation;
  
  // Wt::WComboBox *m_u_pu_data_source;
  
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
};//class RelActAutoGui


#endif //RelActAutoGui_h
