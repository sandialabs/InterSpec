#ifndef Recalibrator_h
#define Recalibrator_h
/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>
#include <Wt/WLayout>

#include "Minuit2/FCNBase.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/RowStretchTreeView.h"

//Forward Declarations
class PeakDef;
class Recalibrator;
class InterSpec;
class SpectraFileModel;
#if ( USE_SPECTRUM_CHART_D3 )
class D3SpectrumDisplayDiv;
#else
class SpectrumDisplayDiv;
#endif

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
  Wt::WDoubleSpinBox *m_energy, *m_offset;
  Wt::WText *m_delete;
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
protected:
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
                              const SpectrumType newtype,
                              std::shared_ptr<SpecMeas> oldmeas,
                              const SpectrumType oldtype,
                              Recalibrator *calibrator );
  
  static bool candidate( std::shared_ptr<SpecMeas> newmeas,
                         std::shared_ptr<SpecMeas> oldmeas );
  
  virtual ~PreserveCalibWindow();
  
  void doRecalibration();
  
protected:
  Recalibrator *m_calibrator;
  std::shared_ptr<SpecMeas> m_newmeas;
  std::vector<float> m_coeffs;
  Measurement::EquationType m_type;
  std::vector< std::pair<float,float> > m_devPairs;
  const SpectrumType m_newtype, m_oldtype;
};//class PreserveCalibWindow


class Recalibrator : public Wt::WContainerWidget
{
#if( !REBIN_FILES_TO_SINGLE_BINNING )
#error Recalibrator is currently only setup to work properly when REBIN_FILES_TO_SINGLE_BINNING is defined
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

  //setRecalibratorLayout(...): Needed to set the layout before it toggles
  //  wide/tall display layout.
  //Note that this really doesnt seem like the best solution.  It would be nice
  //  if we didnt have to externally set the layout for this widget, and instead
  //  everything was self consistent
  void setRecalibratorLayout( Wt::WGridLayout *layout );
  
  //refreshRecalibrator(): fills out GUI from values in the data structures,
  //  makes the backup copies of the data so you can undo, and clears out
  //  the information about the last calibration (m_lastGraphicalRecal,...)
  void refreshRecalibrator();
  
  void engageRecalibration( RecalAction action );
  
  void specTypeCheckedCallback( const SpectrumType type, const bool checked );
  
  
  void handleGraphicalRecalRequest( double xstart, double xfinish );
  void deleteGraphicalRecalConfirmWindow();
  
  WContainerWidget* getFooter();
  
  //currently old_eqn_type and and the new eqn type can only be Polyniomial or
  //  LowerChannelEdge, or an exception will be thrown.
  //Call this function after recalibrating the passed in meas, so it already
  //  has the new calibration the peaks will be moved to.
  static void shiftPeaksForEnergyCalibration(
                                      PeakModel *peakmodel,
                                      const std::vector<float> &new_pars,
                                      const std::vector< std::pair<float,float> > &new_devpairs,
                                      Measurement::EquationType new_eqn_type,
                                      std::shared_ptr<SpecMeas> meas,
                                      const SpectrumType spectype,
                                      std::vector<float> old_pars,
                                      const std::vector< std::pair<float,float> > &old_devpairs,
                                      Measurement::EquationType old_eqn_type );

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
  void initWidgets( LayoutStyle style, AuxWindow* parent = NULL);

  void recalibrateByPeaks();
    
  void createMultifileFitter();

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
  Wt::WCheckBox *m_fitFor[3];
  Wt::WPushButton *m_fitCoefButton;

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
  
  // Which datafields to apply it to
  Wt::WCheckBox *m_applyTo[3];
  
  // The buttons that get pushed beep boop
  Wt::WPushButton *m_revert, *m_closeButton;

  Wt::WText *m_acceptText;
  Wt::WContainerWidget *m_acceptButtonDiv;

  Wt::WDoubleSpinBox *m_coefficientDisplay[3];
  
  // The boxes for the manual recalibration coefficients
  //Wt::WDoubleSpinBox *m_offset, *m_linear, *m_quadratic;
  
  // holders for the exponents
  Wt::WLabel *m_exponentLabels[3];
  
  // and the ints for the scientific notation
  int m_coeffExps[3];

  Measurement::EquationType m_coeffEquationType;

  //The following will be set to a negative number when they arent defined
//  double m_offsetUncert, m_linearUncert, m_quadraticUncert;
  double m_uncertainties[3];

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
    Measurement::EquationType type;
    std::vector<float> coefficients;
    std::vector< std::pair<float,float> > deviationpairs;
    
    std::set<int> sample_numbers;
    std::set<int> detectors_numbers;
    
    //could also hold original SpecMeas::SampleNumsToPeakMap ...
  };//struct CalibrationInformation
  
  CalibrationInformation m_originalCal[3];
  
protected:


  struct RecalPeakInfo
  {
    double peakMean;
    double peakMeanUncert;
    double peakMeanBinNumber;
    double photopeakEnergy;
  };//struct RecalPeakInfo

  class PolyCalibCoefMinFcn
      : public ROOT::Minuit2::FCNBase
  {
  public:
    PolyCalibCoefMinFcn( const std::vector<RecalPeakInfo> &peakInfo,
                         const size_t nbin,
                         Measurement::EquationType eqnType,
                         const std::vector< std::pair<float,float> > &devpair );
    virtual ~PolyCalibCoefMinFcn();
    virtual double Up() const;
    virtual double operator()( const std::vector<double> &x ) const;

  protected:
    size_t m_nbin;
    Measurement::EquationType m_eqnType;
    std::vector<RecalPeakInfo> m_peakInfo;
    std::vector< std::pair<float,float> > m_devpair;
  };//class PolyCalibCoefMinFcn

  
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
    Wt::WCheckBox *m_fitFor[3];
    Wt::WLineEdit *m_coefvals[3];
    Wt::WPushButton *m_use, *m_cancel, *m_fit;
    Wt::WTextArea *m_fitSumary;
    Measurement::EquationType m_eqnType;
    double m_calVal[3], m_calUncert[3];
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


