#ifndef PeakEdit_h
#define PeakEdit_h
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

#include <Wt/WString>
#include <Wt/WModelIndex>
#include <Wt/WContainerWidget>

#include "InterSpec/PeakDef.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/ReactionGamma.h"

class PeakEdit;
class PeakModel;
class InterSpec;
class ColorSelect;
class IsotopeNameFilterModel;
class PeakIsotopeNameFilterModel;

namespace Wt
{
  class WText;
  class WTable;
  class WComboBox;
  class WCheckBox;
  class WLineEdit;
  class WPushButton;
  class WSuggestionPopup;
}//namespace Wt

class PeakEditWindow : public AuxWindow
{
public:
  PeakEditWindow( const double energy,
                  PeakModel *peakmodel,
                  InterSpec *viewer );
  PeakEdit *peakEditor();
  Wt::Signal<> &editingDone();
protected:
  PeakEdit *m_edit;
};//class PeakEditWindow


class PeakEdit : public Wt::WContainerWidget
{
public:
  //PeakPars gives the order of fields in the edit widget.
  //  Currently extends the PeakDef::CoefficientType enum, in a kinda hacky way;
  //  there are some static asserts that try to ensure consistentcy incase
  //  PeakDef::CoefficientType is changed, hopefully to compile,
  //  PeakEdit::PeakPars will have to be altered as well.
  enum PeakPars
  {
    Mean,
    Sigma,
    SigmaDrfPredicted,
    GaussAmplitude,
    LandauAmplitude,
    LandauMode,
    LandauSigma,
    Chi2DOF,
    RangeStartEnergy,
    RangeEndEnergy,
    OffsetPolynomial0,
    OffsetPolynomial1,
    OffsetPolynomial2,
    OffsetPolynomial3,
    OffsetPolynomial4,
    PeakColor,
    NumPeakPars
  };//enum PeakPars
  
  /** Convert from PeakPars to PeakDef::CoefficientType.
   
   Throws exception if does not directly map; e.g., `type` is larger than #PeakPars::RangeStartEnergy or is equal
   to #PeakPars::SigmaDrfPredicted.
   */
  static PeakDef::CoefficientType row_to_peak_coef_type( const PeakPars type );
  
public:
  PeakEdit( const double energy,
            PeakModel *peakModel,
            InterSpec *viewer,
            AuxWindow* aux = 0 );
  virtual ~PeakEdit();
  
  //changePeak(...): only changes the peak if the energy is within the ROI of
  //  a valid peak
  void changePeak( const double energy );

  double currentPeakEnergy() const;
  
  Wt::Signal<> &done();

  //isEditingValidPeak(): is the editor actually editing a valid peak, or blank
  bool isEditingValidPeak() const;
  
  static const char *rowLabel( const PeakPars par );
  
protected:
  void init();
  void refreshPeakInfo();
  void peakModelRowsRemoved( Wt::WModelIndex, int, int );
  void peakModelRowsAdded( Wt::WModelIndex, int, int );
  void checkIfDirty( PeakPars type, bool uncert );
  void checkIfColorDirty();
  void checkIfUserLabelDirty();
  void peakTypeChanged();
  void contnuumTypeChanged();
  void skewTypeChanged();
  void validateMeanOrRoiChange();
  void changeToNextPeakInRoi( const bool back );
  
  //fitTypeChanged(...): updates the peak in PeakModel, as well as
  //  m_currentPeak, but not m_originalPeak
  void fitTypeChanged( PeakPars t );
  
  void setAmplitudeForDataDefinedPeak();
  
  void isotopeChanged();
  void transitionChanged();
  bool nuclideInfoIsDirty() const;
  //setNuclideFields(...): populates m_nuclide and m_photoPeakEnergy.  If a
  //  transition and particle_index is specified, then the correct
  //  m_photoPeakEnergy index will be chosen, otherwise one will be guessed.
  void setNuclideFields( const SandiaDecay::Nuclide *nuclide,
                         const double mean, const double sigma,
                         const SandiaDecay::Transition *transition,
                         const int particle_index,
                         const PeakDef::SourceGammaType sourceGammaType );
  void setNuclideFields( const SandiaDecay::Element *el, const float mean );
  void setNuclideFields( const ReactionGamma::Reaction *rctn, const float mean );
  
  /** Returns if there is more than one color associated with the
      nuclide/element/reaction of the current peak.  Also sets the label
      for #m_applyColorForAllNuc appriately.
   */
  bool checkNuclideForDiffColors();
  
  void refit();
  void apply();  //throws if an input is invalid
  void accept();
  void cancel();
  void deletePeak();
  
  bool isDirty() const;
  
  void updateDrfFwhmTxt();
protected:
  double m_energy;
  PeakModel *m_peakModel;
  InterSpec *m_viewer;
  
  Wt::WModelIndex m_peakIndex;
  PeakDef m_originalPeak;
  PeakDef m_currentPeak;
  bool m_blockInfoRefresh;
  
  Wt::WTable *m_valueTable;
  bool m_valIsDirty[NumPeakPars];
  bool m_uncertIsDirty[NumPeakPars];
  
  //m_valTxts, and m_uncertTxts could be illiminated if the changed() signal is
  //  used instead of the blurred() and enterPressed() signals
  Wt::WString m_valTxts[NumPeakPars];
  Wt::WString m_uncertTxts[NumPeakPars];
  
  Wt::WLineEdit *m_values[NumPeakPars];
  Wt::WLineEdit *m_uncertainties[NumPeakPars];
  Wt::WCheckBox *m_fitFors[NumPeakPars];

  Wt::WLineEdit *m_nuclide;
  Wt::WSuggestionPopup *m_suggestions;
  
  //I have temporarily switched from the preffered PeakIsotopeNameFilterModel
  //  to IsotopeNameFilterModel, to allow entry of reactions and xrays.
//  PeakIsotopeNameFilterModel *m_filterModel;
  IsotopeNameFilterModel *m_filterModel;
  Wt::WComboBox *m_photoPeakEnergy;

  Wt::WLineEdit *m_userLabel;
  ColorSelect *m_color;
  Wt::WCheckBox *m_applyColorForAllNuc;
  Wt::WComboBox *m_peakType;
  Wt::WComboBox *m_continuumType;
  Wt::WComboBox *m_skewType;
  
  Wt::WPushButton *m_apply;
  Wt::WPushButton *m_accept;
  Wt::WPushButton *m_cancel;
  Wt::WPushButton *m_refit;
  
  Wt::WContainerWidget *m_otherPeaksDiv;
  Wt::WText *m_otherPeakTxt;
  Wt::WPushButton *m_prevPeakInRoi;
  Wt::WPushButton *m_nextPeakInRoi;
  
  Wt::WText *m_drfFwhm;
  
  Wt::Signal<> m_doneSignal;
  Wt::WContainerWidget *m_footer;
  AuxWindow *m_aux;
};//class PeakEdit



#endif  //PeakEdit_h
