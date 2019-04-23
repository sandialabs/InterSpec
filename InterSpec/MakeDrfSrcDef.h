#ifndef MakeDrfSrcDef_h
#define MakeDrfSrcDef_h
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

#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>

#include <Wt/WContainerWidget>

class MaterialDB;
class ShieldingSelect;

namespace Wt
{
  class WText;
  class WLabel;
  class WTable;
  class WLineEdit;
  class WCheckBox;
  class WDateEdit;
  class WComboBox;
  class WDoubleSpinBox;
  class WSuggestionPopup;
}//namespace Wt

namespace SandiaDecay
{
  struct Nuclide;
}//namespace SandiaDecay

/* Widget to allow entering:
   - Nuclide
   - Activity
   - Date of activity
   - Date of measurment being used for (can likely be auto-filed out)
   - Potentially age when activity was determined
   - Activity uncertainty
   - Distance to source
   - Shielding
 */
class MakeDrfSrcDef : public Wt::WContainerWidget
{
public:
  MakeDrfSrcDef( const SandiaDecay::Nuclide *nuc,
                 const boost::posix_time::ptime &measDate,
                 MaterialDB *materialDB,
                 Wt::WSuggestionPopup *materialSuggest,
                 Wt::WContainerWidget *parent = nullptr );
  
  virtual ~MakeDrfSrcDef();
  
  
protected:
  /** Creates the widget, but doesnt fill out any of the information. */
  void create();
  
  void handleNuclideChanged();
  
  void useAgeInfoUserToggled();
  void useShiledingInfoUserToggled();
  
  void handleUserChangedActivity();
  double enteredActivity();
  
  Wt::WTable *m_table;
  
  /** The name of the nuclide this source represents. */
  const SandiaDecay::Nuclide *m_nuclide;
  
  /** Pointer to shielding material database, garunteed non-null.
      Not owned by this object.
   */
  MaterialDB *m_materialDB;
  
  /** The material suggestion widget shared everywher within InterSpec.
      Not owned by this object.
   */
  Wt::WSuggestionPopup *m_materialSuggest;
  
  /** The user input of the nuclide. */
  Wt::WLineEdit *m_nuclideEdit;
  
  Wt::WLineEdit *m_distanceEdit;
  
  
  /** */
  Wt::WLineEdit *m_activityEdit;
  
  Wt::WComboBox *m_activityUnits;
  
  Wt::WDoubleSpinBox *m_activityUncertainty;
  
  /** Whether or not age/decay info should be entered (default: no) */
  Wt::WCheckBox *m_useAgeInfo;
  
  /** The date entered activity was determined. */
  Wt::WDateEdit *m_activityDate;
  
  /** The date the spectrum being used for creating the DRF was taken on. */
  Wt::WDateEdit *m_drfMeasurementDate;
  
  /** The computed source activity when measurement was taken. */
  Wt::WText *m_sourceActivityAtMeasurement;
  
  /** The source age at the time of the measurement. */
  Wt::WText *m_sourceAgeAtMeasurement;
  
  /** The date the source was created on, if applicable */
  //For the moment, to save complexity, we wont enable having sources aged
  //  relative to their assay date.
  //Wt::WDateEdit *m_sourceCreationDate;
  
  /** Whether or not to have shielding for the source (default: no) */
  Wt::WCheckBox *m_useShielding;
  
  /** The widget to allow selecting shielding, is selected to do so. */
  ShieldingSelect *m_shieldingSelect;
  
};//MakeDrfSrcDef

#endif //MakeDrfSrcDef_h
