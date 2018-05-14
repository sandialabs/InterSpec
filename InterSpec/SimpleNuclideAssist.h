#ifndef SimpleNuclideAssist_h
#define SimpleNuclideAssist_h

#include "InterSpec_config.h"

#include <string>

#include <Wt/WContainerWidget>

#include "InterSpec/PopupDiv.h"

class SimpleNuclide;
class SpectrumViewer;
class SimpleNuclideAssist;

namespace Wt
{
  class WText;
}

class SimpleNuclideAssistPopup : public PopupDivMenu
{
public:
  SimpleNuclideAssistPopup( const float energy, SpectrumViewer *viewer,
                            int pagex, int pagey );
  virtual ~SimpleNuclideAssistPopup();
  
  bool isValid() const;
  void handleAboutToHide();
  
protected:
  PopupDivMenuItem *m_item;
  SimpleNuclideAssist *m_assist;
};//class SimpleNuclideAssistPopup



class SimpleNuclideAssist : public Wt::WContainerWidget
{
/*
 TODO 20150708 (sorted roughly by priority):
  -Finish filling out information in the XML!
  -Add in parsing <someshielding> and <moreshielding> from XML and use them
  -Improve the window width for the ROI, especially for low resolution detectors
  -Add in 'equivalent' nuclides into the more information popup
  -Make sure all nuclides/reactions in XML are available in Reference Photopeaks
  -Finish the styling of everything
  -Improve the message to users about what this widget actually is
  -Test on mobile
  -Add in a link to a help page describing this widget
  -Maybe add in a close button/box
  -Check on XML parse speed and maybe make it so its only parsed once, or try 
   memory mapping the XML file
  -Maybe make energy range changeable - maybe
*/
public:
  SimpleNuclideAssist( const float energy,
                       SpectrumViewer *viewer,
                       Wt::WContainerWidget *parent = 0 );
  
  virtual ~SimpleNuclideAssist();
  
  //setEnergy(): convience function that determines an aproximately appropriate
  //  window around `energy` using the detector response function, or fit peaks,
  //  but enforces at least 7 px on either side of energy.
  void setEnergy( const float energy );
  
  
  void setEnergyRange( const float lowerEnergy, const float upperEnergy );
  
  void updateSourceDisplayed( SimpleNuclide *widget );
  
  void setShielding( const std::string &name, const std::string &thickness );
  
protected:
  SpectrumViewer *m_viewer;
  std::string m_refPhotopeakIntialStateXml;
  Wt::WContainerWidget *m_sources;
  Wt::WContainerWidget *m_messages;
  Wt::WContainerWidget *m_legend;
  Wt::WText *m_energyRangeTxt;
  Wt::WText *m_contaminantLegTxt;
  Wt::WContainerWidget *m_contaminantLegImg;
  
  bool m_initialized;
  float m_lowerEnergy;
  float m_upperEnergy;
  
  std::vector<size_t> m_regionids;

public:
  static const std::string sm_dataFileName;
};//class SimpleNuclideAssist

#endif  //SimpleNuclideAssist_h
