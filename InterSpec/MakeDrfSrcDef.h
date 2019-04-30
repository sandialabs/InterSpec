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
  
  /** Returns the nuclide this MakeDrfSrcDef is for.
   
   Will be nullptr if that is what was passed in.
   */
  const SandiaDecay::Nuclide *nuclide() const;
  
  /** Returns the nuclide activity, at the time the characterization measurment
     was taken.
   
     Throws exception if there are errors in any of the user input fileds.
   */
  double activityAtSpectrumTime() const;
  
  /** Returns the nuclide age (in units of PhysicalUnits), at the time the
   characterization measurment was taken.  Returns the age you should assume the
   source is at the time of the characterization, not necassarily the actual
   age.
   
   Throws exception if there are errors in any of the user input fileds.
   */
  double ageAtSpectrumTime() const;
  
  /** Returns nullptr if the user has not selected a shielding, other wise
     returns the ShieldingSelect.
   */
  ShieldingSelect *shielding();
  
  
protected:
  /** Creates the widget, but doesnt fill out any of the information. */
  void create();
  
  void setNuclide( const SandiaDecay::Nuclide *nuc );
  
  void updateAgedText();
  
  void useAgeInfoUserToggled();
  
  void useShieldingInfoUserToggled();
  
  void handleUserChangedActivity();
  
  void handleUserChangedAgeAtAssay();
  
  void handleEnteredDatesUpdated();
  
  /** Returns user entered activity.  Throws exception if invalid. */
  double enteredActivity() const;
  
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
  
  /** Display of the nuclide. */
  Wt::WText *m_nuclideLabel;
  
  Wt::WLineEdit *m_distanceEdit;
  
  
  /** */
  Wt::WLineEdit *m_activityEdit;
  
  Wt::WComboBox *m_activityUnits;
  
  Wt::WDoubleSpinBox *m_activityUncertainty;
  
  /** Whether or not age/decay info should be entered (default: no) */
  Wt::WCheckBox *m_useAgeInfo;
  
  /** The date entered activity was determined. */
  Wt::WDateEdit *m_assayDate;
  
  /** The date the spectrum being used for creating the DRF was taken on. */
  Wt::WDateEdit *m_drfMeasurementDate;
  
  /** The source age and activity at the time of the measurement. */
  Wt::WText *m_sourceInfoAtMeasurement;
  
  /** The source age at assay, if applicable */
  Wt::WLineEdit *m_sourceAgeAtAssay;
  
  /** Whether or not to have shielding for the source (default: no) */
  Wt::WCheckBox *m_useShielding;
  
  /** The widget to allow selecting shielding, is selected to do so. */
  ShieldingSelect *m_shieldingSelect;
  
};//MakeDrfSrcDef

#endif //MakeDrfSrcDef_h
