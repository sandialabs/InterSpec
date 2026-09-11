#ifndef EccUncertOptions_h
#define EccUncertOptions_h
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

#include <Wt/WSignal.h>
#include <Wt/WContainerWidget.h>

class DetectorEfficiencyUncert;

namespace Wt
{
  class WText;
  class WTable;
  class WLabel;
  class WComboBox;
  class WCheckBox;
  class WLineEdit;
}//namespace Wt

/** Reusable UI for choosing how ISOCS .ecc file uncertainties are imported.

 Given the raw per-energy baseline (correlated %err) and convergence
 (uncorrelated %cnvrg) fractional uncertainties parsed from an .ecc file, this
 widget lets the user:
   - toggle importing the uncertainties at all,
   - choose the baseline correlation mode across energy (fully correlated -
     the default, matching GENIE; Gaussian in log-energy with a user-set
     correlation length; or uncorrelated), and
   - preview the resulting correlation for a few representative energy pairs.

 The build logic lives here so both interactive .ecc import paths (the
 drag-drop dialog and the "DRF Select" upload tab) share one implementation.
 */
class EccUncertOptions : public Wt::WContainerWidget
{
public:
  /** @param energies Node energies (keV) parsed from the .ecc file.
      @param baselineFrac Correlated fractional 1-sigma uncertainties (%err/100).
      @param convergenceFrac Uncorrelated fractional 1-sigma uncertainties
             (%cnvrg/100); may be empty.
   */
  EccUncertOptions( const std::vector<float> &energies,
                    const std::vector<float> &baselineFrac,
                    const std::vector<float> &convergenceFrac );

  virtual ~EccUncertOptions();

  /** Whether the user has opted to import the uncertainties. */
  bool importUncertainties() const;

  /** The correlation length (ln-energy units) implied by the current UI state:
     <= 0 for uncorrelated, the parsed field value for Gaussian, or the
     fully-correlated sentinel otherwise.
   */
  double effectiveCorrLength() const;

  /** Builds the efficiency uncertainty for the current UI selection; returns
     nullptr when import is unchecked or there are < 2 nodes.
   */
  std::shared_ptr<DetectorEfficiencyUncert> buildUncert() const;

  /** Emitted when the import toggle, correlation mode, or length changes. */
  Wt::Signal<> &changed();

protected:
  void handleModeChanged();
  void handleImportToggled();
  void rebuildExampleTable();

  const std::vector<float> m_energies;
  const std::vector<float> m_baselineFrac;
  const std::vector<float> m_convergenceFrac;

  Wt::WCheckBox *m_import;
  Wt::WComboBox *m_mode;
  Wt::WLabel *m_corrLenLabel;
  Wt::WLineEdit *m_corrLen;
  Wt::WText *m_exampleTitle;
  Wt::WTable *m_exampleTable;

  Wt::Signal<> m_changed;

  /** Combo indices for #m_mode. */
  enum Mode
  {
    FullyCorrelated = 0,
    Gaussian = 1,
    Uncorrelated = 2
  };//enum Mode
};//class EccUncertOptions

#endif //EccUncertOptions_h
