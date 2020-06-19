#ifndef Recalibrator_h
#define Recalibrator_h
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

#include <map>
#include <vector>
#include <time.h>
#include <utility>

#include <Wt/WLayout>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"


//Forward Declarations
class PeakDef;
class Recalibrator;
class InterSpec;
class NativeFloatSpinBox;
class SpectraFileModel;
#if ( USE_SPECTRUM_CHART_D3 )
class D3SpectrumDisplayDiv;
#else
class SpectrumDisplayDiv;
#endif
class RowStretchTreeView;
namespace SpecUtils{ enum class EnergyCalType : int; }
namespace SpecUtils{ enum class SpectrumType : int; }

namespace Wt
{
  class WText;
  class WTextArea;
  class WLineEdit;
  class WPushButton;
  class WLayout;
}//namespace Wt


class DeviationPairDisplay;
class DevPair : public Wt::WContainerWidget
{
protected:
  DevPair( Wt::WContainerWidget *parent = 0 );
  void setDevPair( const std::pair<float,float> &d );
  std::pair<float,float> devPair() const;
  void visuallyIndicateChanged();
  //Wt::WDoubleSpinBox *m_energy, *m_offset;
  NativeFloatSpinBox *m_energy, *m_offset;
  Wt::WContainerWidget *m_delete;
  friend class DeviationPairDisplay;
};//class DevPair
  
class DeviationPairDisplay : public Wt::WContainerWidget
{
public:
  DeviationPairDisplay( Wt::WContainerWidget *parent = 0 );
  void setDeviationPairs( std::vector< std::pair<float,float> > d );
  std::vector< std::pair<float,float> > deviationPairs() const;
  void removeDevPair( DevPair *devpair );
  DevPair *newDevPair( const bool emitChangedNow );
  Wt::Signal<> &changed(); //emits when dev pair is added, deleted, or changed
  
  void setInvalidValues();
  void setValidValues();
  
protected:
  void sortDisplayOrder( const bool indicateVisually );
  
  void emitChanged();
    
  Wt::Signal<> m_changed;
  Wt::WContainerWidget *m_pairs;
};//class DeviationPairDisplay


//PreserveCalibWindow: propogates calibration from one SpecMeas to
//  another.
class PreserveCalibWindow : public AuxWindow
{
public:
  PreserveCalibWindow( std::shared_ptr<SpecMeas> newmeas,
                              const SpecUtils::SpectrumType newtype,
                              std::shared_ptr<SpecMeas> oldmeas,
                              const SpecUtils::SpectrumType oldtype,
                              Recalibrator *calibrator );
  
  static bool candidate( std::shared_ptr<SpecMeas> newmeas,
                         std::shared_ptr<SpecMeas> oldmeas );
  
  virtual ~PreserveCalibWindow();
  
  void doRecalibration();
  
protected:
  Recalibrator *m_calibrator;
  std::shared_ptr<SpecMeas> m_newmeas;
  std::vector<float> m_coeffs;
  SpecUtils::EnergyCalType m_type;
  std::vector< std::pair<float,float> > m_devPairs;
  const SpecUtils::SpectrumType m_newtype, m_oldtype;
};//class PreserveCalibWindow


class Recalibrator : public Wt::WContainerWidget
{
#if( !SpecUtils_REBIN_FILES_TO_SINGLE_BINNING )
#error Recalibrator is currently only setup to work properly when SpecUtils_REBIN_FILES_TO_SINGLE_BINNING is defined
#endif

public:
  enum LayoutStyle{ kWide, kTall };
  enum RecalAction
  {
    RevertRecal,  //reverts calibration back to last time refreshRecalibrator()
                  //  was called.
    ApplyRecal    //Applies the currently displayed calibration, but does not
                  //  reset the revert to point
  };//enum RecalAction
  
public:
  Recalibrator(
#if ( USE_SPECTRUM_CHART_D3 )
                D3SpectrumDisplayDiv *displayDiv,
#else
                SpectrumDisplayDiv *displayDiv,
#endif
                InterSpec *hostViewer,
                PeakModel *peakModel);
  virtual ~Recalibrator();
  
  //refreshRecalibrator(): fills out GUI from values in the data structures,
  //  makes the backup copies of the data so you can undo, and clears out
  //  the information about the last calibration (m_lastGraphicalRecal,...)
  void refreshRecalibrator();
  
  /** Sets the displayed coefficients labels (keV, keV/chnl, etc.) based on
   m_coeffEquationType.
   */
  void updateCoefLabels();
  
  void engageRecalibration( RecalAction action );
  
  void userChangedDeviationPairs();
  
  void specTypeCheckedCallback( const SpecUtils::SpectrumType type, const bool checked );
  
  
  void handleGraphicalRecalRequest( double xstart, double xfinish );
  void deleteGraphicalRecalConfirmWindow();
  
  //currently old_eqn_type and and the new eqn type can only be Polyniomial or
  //  LowerChannelEdge, or an exception will be thrown.
  //Call this function after recalibrating the passed in meas, so it already
  //  has the new calibration the peaks will be moved to.
  static void shiftPeaksForEnergyCalibration(
                                      PeakModel *peakmodel,
                                      const std::vector<float> &new_pars,
                                      const std::vector< std::pair<float,float> > &new_devpairs,
                                      SpecUtils::EnergyCalType new_eqn_type,
                                      std::shared_ptr<SpecMeas> meas,
                                      const SpecUtils::SpectrumType spectype,
                                      std::vector<float> old_pars,
                                      const std::vector< std::pair<float,float> > &old_devpairs,
                                      SpecUtils::EnergyCalType old_eqn_type );

  static bool checkFullRangeFractionCoeffsValid( const std::vector<float> &pars,
                          const std::vector< std::pair<float,float> > &devpairs,
                                                size_t nbin );
  
  //Note that to accomidate the "tall" and "wide" layout, some fairly
  //  signinifacnt amount of haggling with m_layout had to be done, so things
  //  with the layout are a bit awkward.  Aldso note that the m_isotopeView
  //  couldnt span rows all the way down to the bottom of m_layout or else
  //  the rendering was bad
  //XXX - if we are coming back from a "Wide" view, then I cant figure
  //  out how to place m_fitCoefButton back into m_layout reliably - so I just
  //  leave it in a kinda awkward spot for now
  void setWideLayout();
  void setTallLayout(AuxWindow* parent);

  LayoutStyle currentLayoutStyle() const;
    
protected:
  /** The number of energy calibration coefficients to display. */
  static const size_t sm_numCoefs = 4;
  
  void initWidgets( LayoutStyle style, AuxWindow* parent = NULL);

  void recalibrateByPeaks();
    
  void createMultifileFitter();

  /** Function called when user requests that a lower channel energy calibration
      be converted to polynomial.
   */
  void startConvertToPolynomial();
  
  /** Function that performs the work of converting lower channel energy to
     polynomial calibraiton.
   
   ToDo: Currently discards any non-linearities in the spectrum by just using a
         linear calibration that matches the old energy range.  Could try to
         fit a higher order polynomial, maybe with deviation pairs, or could
         provide an option that would allow the user to choose if they want this
         or to destructivly rebin counts to match a new polynomial...
   ToDo: Currently applies changes to the entire files of foreground,
         background, and secondary spectrum.  Should provide more fine-grained
         controls for this.
   */
  void finishConvertToPolynomial( const bool followThrough );
  
  /** Sets the "Fit Coeffs" button as enabled/disabled based on:
    - Foreground is displayed and has "Apply to" checkbox checked.
    - Peaks to use as calibration peaks are present.
    - At least one term is checked for "fit"
   */
  void checkIfCanFitCoefs();
  
  /** Enables/disables the "Multi. Files..." button, based on if multiple files are loaded, or multiple spectra from
   current file have peaks.
   
   Since not all loaded files may be in memory, doesnt check for peaks except for in current foreground specral file,
   and instead assumes if more than one file is loaded, than button should be enabled.  Also doesnt (currently) look
   that defined peaks are selected to be used for calibration.  E.g., could be tightened up a okay amount.
   */
  void checkIfCanMultiFileCal();
  
  LayoutStyle m_currentLayout;
  
  WContainerWidget *footer;
#if ( USE_SPECTRUM_CHART_D3 )
  D3SpectrumDisplayDiv *m_spectrumDisplayDiv;
#else
  SpectrumDisplayDiv *m_spectrumDisplayDiv;
#endif
  InterSpec *m_hostViewer;
  PeakModel *m_peakModel;

  RowStretchTreeView *m_isotopeView;
  Wt::WCheckBox *m_fitFor[sm_numCoefs];
  Wt::WPushButton *m_fitCoefButton;

  Wt::WPushButton *m_multiFileButton;
  
  class GraphicalRecalConfirm;
  friend class GraphicalRecalConfirm;
  
  class MultiFileCalibFit;
  friend class MultiFileCalibFit;
  
  class MultiFileModel;
  friend class MultiFileModel;
  
  friend class PreserveCalibWindow;
  
  GraphicalRecalConfirm *m_graphicalRecal;
  

  //Fields which apply to manual and peak calibration
  std::vector< Wt::WCheckBox * > m_detectors;
  
  Wt::WText *m_applyToLabel;
  
  // Which datafields (foreground, background, secondary) to apply it to
  Wt::WCheckBox *m_applyTo[3];
  
  Wt::WText *m_convertToPolynomialLabel;
  Wt::WPushButton *m_convertToPolynomial;
  /** Workaround for apparent bug in (at least) Wt 3.3.4, where #m_convertToPolynomialLabel never actually gets hidden
   when #WText::hide is called on it - very bizare, so instead will hide this container, which does actually work.
   */
  Wt::WContainerWidget *m_polyConverDiv;
  
  AuxWindow *m_convertToPolyDialog;
  
  Wt::WPushButton *m_revert;

  Wt::WText *m_acceptText;

  NativeFloatSpinBox *m_coefficientDisplay[sm_numCoefs];
  
  // The boxes for the manual recalibration coefficients
  //Wt::WDoubleSpinBox *m_offset, *m_linear, *m_quadratic;
  
  // holders for the exponents
  Wt::WLabel *m_exponentLabels[sm_numCoefs];
  
  SpecUtils::EnergyCalType m_coeffEquationType;

  //The following will be set to a negative number when they arent defined
//  double m_offsetUncert, m_linearUncert, m_quadraticUncert;
  double m_uncertainties[sm_numCoefs];

  DeviationPairDisplay *m_devPairs;
  
  Wt::WGridLayout* m_layout;
  time_t m_lastGraphicalRecal;
  int m_lastGraphicalRecalType;
  float m_lastGraphicalRecalEnergy;
  
  //setWasGraphicalRecal( int type ): 'type' is of type
  //  GraphicalRecalConfirm::RecalTypes,
  void setWasGraphicalRecal( int type, float energy );
  
  
  struct CalibrationInformation
  {
    void reset();
    CalibrationInformation();
    SpecUtils::EnergyCalType type;
    std::vector<float> coefficients;
    std::vector< std::pair<float,float> > deviationpairs;
    
    std::set<int> sample_numbers;
    std::set<int> detectors_numbers;
    
    //could also hold original SpecMeas::SampleNumsToPeakMap ...
  };//struct CalibrationInformation
  
  CalibrationInformation m_originalCal[3];
  
public:
  struct RecalPeakInfo
  {
    double peakMean;
    double peakMeanUncert;
    double peakMeanBinNumber;
    double photopeakEnergy;
  };//struct RecalPeakInfo

protected:

  
  class GraphicalRecalConfirm : public AuxWindow
  {
  public:
    enum RecalTypes
    {
      kOffset, kLinear,
      /*kQuadratic,*/
      kDeviation, NumRecalTypes
    };//enum RecalTypes
    
    GraphicalRecalConfirm( double lowe, double highe, Recalibrator *cal,
                      time_t lastCal, RecalTypes lastType, float lastEnergy );
    virtual ~GraphicalRecalConfirm();
    void apply();
    void setEnergies( double xstart, double xfinish );

  protected:
    Recalibrator *m_calibrator;
    Wt::WButtonGroup *m_typeButtons;
    Wt::WCheckBox *m_foregroundOnly;
    Wt::WDoubleSpinBox *m_startE, *m_finalE;
    Wt::WCheckBox *m_preserveLastCal;
    const RecalTypes m_lastType;
    const float m_lastEnergy;
  };//class GraphicalRecalConfirm
  
  //MultiFileCalibFit and MultiFileModel are to allow the user to calibrate
  //  a number of spectrum files simultaneously using peaks from each file.
  //  As of 20131213 this capability is only roughed in and not tested very well
  //XXX - applying the calibration parameters does not shift peaks!
  class MultiFileCalibFit : public AuxWindow
  {
  public:
    MultiFileCalibFit( Recalibrator *cal );
    virtual ~MultiFileCalibFit();
    
    void doFit();
    void handleFinish( Wt::WDialog::DialogCode result );
  protected:
    void updateCoefDisplay();
    
    Recalibrator *m_calibrator;
    MultiFileModel *m_model;
    Wt::WCheckBox *m_fitFor[sm_numCoefs];
    Wt::WLineEdit *m_coefvals[sm_numCoefs];
    Wt::WPushButton *m_use, *m_cancel, *m_fit;
    Wt::WTextArea *m_fitSumary;
    SpecUtils::EnergyCalType m_eqnType;
    double m_calVal[sm_numCoefs], m_calUncert[sm_numCoefs];
  };//class MultiFileCalibFit
  
  class MultiFileModel : public  Wt::WAbstractItemModel
  {
  public:
    MultiFileModel( SpectraFileModel *fileModel, Wt::WObject *parent = 0 );
    virtual ~MultiFileModel();
    
    virtual Wt::WModelIndex index( int row, int column,
                     const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
    
    virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;
    virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
    virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
    
    virtual boost::any data( const Wt::WModelIndex &index,
                             int role = Wt::DisplayRole ) const;
    virtual bool setData( const Wt::WModelIndex &index,
                          const boost::any &value, int role = Wt::EditRole );
    virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const;
    virtual boost::any headerData( int section,
                                   Wt::Orientation orientation = Wt::Horizontal,
                                   int role = Wt::DisplayRole) const;
    void refreshData();
  protected:
    std::vector<std::vector< std::pair<bool,std::shared_ptr<const PeakDef> > > > m_peaks;
    Recalibrator *m_calibrator;
    SpectraFileModel *m_fileModel;
    
    friend class MultiFileCalibFit;
  };//class MultiFileModel
  
};//class Recalibrator

#endif


