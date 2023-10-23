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

#include <vector>

#include <Wt/WGridLayout>
#include <Wt/WPushButton>

#include "SpecUtils/SpecFile.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/MakeFwhmForDrf.h"


using namespace std;
using namespace Wt;



AuxWindow *MakeFwhmForDrf::makeAddFwhmToDrfWindow()
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  if( !viewer )
    return;
  
  AuxWindow *window = new AuxWindow( "Add FWHM to Detector Response Function",
                                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                                     | AuxWindowProperties::SetCloseable
                                     | AuxWindowProperties::DisableCollapse
                                     | AuxWindowProperties::EnableResize
                                     | AuxWindowProperties::IsModal) );
  
  const int ww = viewer->renderedWidth();
  const int wh = viewer->renderedHeight();
  if( ww > 100 && wh > 100 )
  {
    const int width = std::min( 3*ww/4, 900 );
    const int height = ((wh < 420) ? wh : (19*wh)/20 );
    
    window->resizeWindow( width, height );
    window->setMinimumSize( std::min(width,640), std::min(height,480) );
  }//if( ww > 100 && wh > 100 )
  
  shared_ptr<const SpecMeas> foreground = viewer->measurment(SpecUtils::SpectrumType::Foreground );
  shared_ptr<const DetectorPeakResponse> drf = foreground ? foreground->detector() : nullptr;
  
  MakeFwhmForDrf *makeFwhmWidget = new MakeFwhmForDrf( viewer, drf );
    
  window->stretcher()->addWidget( makeFwhmWidget, 0, 0 );
  window->stretcher()->setContentsMargins( 0, 0, 0, 0 );
  
  AuxWindow::addHelpInFooter( window->footer(), "make-drf" );
    
  WPushButton *closeButton = window->addCloseButtonToFooter( "Close" );
  closeButton->clicked().connect( window, &AuxWindow::hide );
    
  WPushButton *saveAs = new WPushButton( "Use", window->footer() );
  saveAs->clicked().connect( makeFwhmWidget, &MakeFwhmForDrf::setToDrf );
  makeFwhmWidget->validationChanged().connect( boost::bind( &WPushButton::setEnabled, saveAs,
                                                                   boost::placeholders::_1 ) );
  // Maybe 
  //makeFwhmWidget->validationChanged().connect( std::bind([makeFwhmWidget,saveAs](){
  //  saveAs->setEnabled( makeFwhmWidget->isValidFwhm() );
  //}) );
  
  saveAs->disable();
    
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    
  window->show();
    
  window->resizeToFitOnScreen();
  window->centerWindow();
  window->rejectWhenEscapePressed( false );
  
  //PeakModel *peakModel = viewer->peakModel();
  //shared_ptr<const deque< PeakModel::PeakShrdPtr > > peakModel->peaks();
  
  return window;
}//AuxWindow *makeDrfWindow(...)


MakeFwhmForDrf::MakeFwhmForDrf( InterSpec *viewer,
               std::shared_ptr<const DetectorPeakResponse> drf,
               Wt::WContainerWidget *parent )
 : WContainerWidget( parent )

InterSpec *m_interspec;
std::shared_ptr<const DetectorPeakResponse> m_orig_drf;
std::vector<std::shared_ptr<const PeakDef>> m_user_peaks;
std::vector<std::shared_ptr<const PeakDef>> m_auto_fit_peaks;

Wt::WComboBox *m_fwhmEqnType;
Wt::WComboBox *m_sqrtEqnOrder;
std::vector<NativeFloatSpinBox *> m_parEdits;

Wt::Signal<bool> m_validationChanged;

{
  
  /*
  
  we will call
  double MakeDrfFit::performResolutionFit( std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
                                          const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                                          const bool highResolution,
                                          const int sqrtEqnOrder,
                                          std::vector<float> &result,
                                          std::vector<float> &uncerts );
WE will show a
  MakeDrfChart(...)
  
We will also call some variant of
   = ExperimentalAutomatedPeakSearch::search_for_peaks( data, drf, existingPeaks, singleThread );
  */
}//MakeFwhmForDrf( constructor )


MakeFwhmForDrf::~MakeFwhmForDrf()
{
  
}//~MakeFwhmForDrf()


Wt::Signal<bool> &MakeFwhmForDrf::validationChanged()
{
  return m_validationChanged;
}


bool MakeFwhmForDrf::isValidFwhm() const
{
  // Returns if current functinoal form is valid
}//bool isValidFwhm() const


void MakeFwhmForDrf::setToDrf()
{
  // Set current FWHM to DRF and broadcast out to the rest of the app
}//void setToDrf()
