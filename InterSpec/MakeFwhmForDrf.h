#ifndef MakeFwhmForDrf_h
#define MakeFwhmForDrf_h
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

#include <memory>
#include <vector>
#include <utility>

#include <Wt/WContainerWidget>

class PeakDef;
class AuxWindow;
class InterSpec;
class MakeDrfChart;
class FwhmPeaksModel;
class NativeFloatSpinBox;
class DetectorPeakResponse;

namespace Wt
{
  class WText;
  class WComboBox;
  class WTableView;
}//namespace Wt

/* Adds FWHM functional information to the DRF passed in; resulting
 DRF can be set to the current one for the application
 */
class MakeFwhmForDrf : public Wt::WContainerWidget
{
public:
  MakeFwhmForDrf( InterSpec *viewer,
                 std::shared_ptr<const DetectorPeakResponse> drf,
                 Wt::WContainerWidget *parent = nullptr );
  
  virtual ~MakeFwhmForDrf() override;

  Wt::Signal<bool> &validationChanged();
  bool isValidFwhm() const;
  
  void setToDrf();
  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &updatedDrf();
  
  
public: //Some stuff for undo/redo support
  struct TableRow
  {
    bool m_is_user_peak;
    bool m_use_for_fit;
    std::shared_ptr<const PeakDef> m_peak;
  };//struct TableRow
  
  
  struct ToolState
  {
    int m_fwhm_index;
    int m_sqrt_eqn_index;
    std::vector<MakeFwhmForDrf::TableRow> m_rows;
    std::vector<float> m_parameters, m_uncertainties;
    
    bool operator==( const ToolState &rhs ) const;
  };//struct ToolState

  std::shared_ptr<ToolState> currentState() const;
  void setState( std::shared_ptr<const ToolState> state );
  
protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags ) override;
  
  void handleTableDataChange();
  void handleFwhmEqnTypeChange();
  void handleSqrtEqnOrderChange();
  
  void coefficientManuallyChanged( const int coef_num );
  
  void scheduleRefit();
  void doRefitWork();
  void setEquationToChart();
  
  void startAutomatedPeakSearch();
  void setPeaksFromAutoSearch( std::vector<std::shared_ptr<const PeakDef>> user_peaks,
                               std::shared_ptr<std::vector<std::shared_ptr<const PeakDef>>> auto_search_peaks );
  
  void doAddUndoRedoStep();
  void scheduleUndoRedoStep();
  
protected:
  InterSpec *m_interspec;
  
  bool m_refit_scheduled;
  bool m_undo_redo_scheduled;
  
  std::shared_ptr<const DetectorPeakResponse> m_orig_drf;
  std::vector<std::shared_ptr<const PeakDef>> m_user_peaks;
  std::vector<std::shared_ptr<const PeakDef>> m_auto_fit_peaks;
  
  MakeDrfChart *m_chart;
  
  Wt::WComboBox *m_fwhmEqnType;
  Wt::WComboBox *m_sqrtEqnOrder;
  std::vector<NativeFloatSpinBox *> m_parEdits;
  std::vector<float> m_parameters;
  std::vector<float> m_uncertainties;
  
  Wt::WText *m_error;
  Wt::WText *m_equation;
  
  Wt::WTableView *m_table;
  FwhmPeaksModel *m_model;
  
  Wt::Signal<bool> m_validationChanged;
  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> m_updatedDrf;
  
  std::shared_ptr<const ToolState> m_current_state;
};//class MakeFwhmForDrf


/** Create a AuxWindow with the MakeFwhmForDrf widget as the primary content.
    Returns pointer to the created AuxWindow, but is will already be shown,
    and have the signals to delete it when closed hooked up, so you probably
    wont need the window.
 */
class MakeFwhmForDrfWindow : public AuxWindow
{
public:
  MakeFwhmForDrfWindow();
  
  MakeFwhmForDrf *tool();
  
protected:
  MakeFwhmForDrf *m_tool;
};//class MakeFwhmForDrfWindow

#endif //MakeFwhmForDrf_h

