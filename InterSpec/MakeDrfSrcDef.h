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
  
  /** Signal emitted when widget is upadated by the user. */
  Wt::Signal<> &updated();
  
  /** Returns user entered distance.
   
     Throws exception if input field has invalid text.
   */
  double distance() const;
  
  /** Returns the nuclide this MakeDrfSrcDef is for.
   
   Will be nullptr if that is what was passed in.
   */
  const SandiaDecay::Nuclide *nuclide() const;
  
  /** Returns the nuclide activity, at the time the characterization measurment
     was taken.
   
     Throws exception if there are errors in any of the user input fileds.
   */
  double activityAtSpectrumTime() const;
  
  /** Returns the fractional (so 0.0 to 1.0) uncertainty on the activity.
   
   Throws exception if there are errors in the user input field.
   */
  double fractionalActivityUncertainty() const;
  
  /** Returns the nuclide age (in units of PhysicalUnits), at the time the
   characterization measurment was taken.  Returns the age you should assume the
   source is at the time of the characterization, not necassarily the actual
   age.
   
   Throws exception if there are errors in any of the user input fileds.
   */
  double ageAtSpectrumTime() const;
  
  /** Returns nullptr if the user has not selected a shielding, currently
   a fixed geometry, other wise returns the ShieldingSelect.
   */
  ShieldingSelect *shielding();
  
  /** Set the distance, maybe from a hint in the spectrum file.
   Does not cause the updated() signal to be emitted.
   */
  void setDistance( const double dist );
  
  /** Set the activity, maybe from a hint in the spectrum file.
   Does not cause the updated() signal to be emitted.
   */
  void setActivity( const double act );
  
  /** Sets the activity and assay date, enabling aging if it isnt already.
   Does not cause the updated() signal to be emitted.
   */
  void setAssayInfo( const double activity, const boost::posix_time::ptime &assay_date );
  
  /** Sets the age at assay date, enabling aging if it isnt already.
   Does not cause the updated() signal to be emitted.
   //Function was implemented but not tested, so leaving commented out until its needed
   */
  //void setAgeAtAssay( const double age );
  
  
  /** Sets the age at the time of measurement.
     If there is already a measurement and assay date, and 'age' is larger than
     their difference, the age at assay field is modified.  If this is not the
     case, the assay date is set equal to measurement date, and age at assay
     field is set to the value passed in.
     If a negative value is passed in, the setNuclide(...) is called, effectively
     reseting things.
   
     Does not cause the updated() signal to be emitted.
   */
  void setAgeAtMeas( const double age );
  
  
  /** Sets the shielding to be generic material with given atomoic number and
     areal density.
     Causes shielding widget to be shown.
   */
  void setShielding( const float atomic_number, const float areal_density );
  
  
  /** Sets geometry type of DRF being made (e.g., hide distance and shielding), or not (default).
   The int passed in is defined by DetectorPeakResponse::EffGeometryType.
   */
  void setIsEffGeometryType( const int drf_eff_geom_type );
  
  
  /** Returns a GADRAS style source string.
      May include source age, if it will matter.
      Will not include source distance.
   
      Examples are: "133Ba,10uCi"
                    "232U,10uC{26,10}"
                    "232U,10uC{26,10} Age=20y"  //The Age is a InterSpec extension and not followed by GADRAS.
   
      Throughs exception if any input widgets are invalid.
   */
  std::string toGadrasLikeSourceString() const;
  
protected:
  /** Creates the widget, but doesnt fill out any of the information. */
  void create();
  
  void setNuclide( const SandiaDecay::Nuclide *nuc );
  
  void updateAgedText();
  
  void useAgeInfoUserToggled();
  
  void useShieldingInfoUserToggled();
  
  void handleUserChangedShielding();
  
  void handleUserChangedDistance();
  
  void handleUserChangedActivity();
  
  void handleUserChangedActivityUncertainty();
  
  void handleUserChangedAgeAtAssay();
  
  void handleEnteredDatesUpdated();
  
  void validateDateFields();
  
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
  
  
  Wt::WLabel *m_distanceLabel;
  
  /** The distance of source to face of detector.
   
   When geometry is fixed, this will be hidden.
   */
  Wt::WLineEdit *m_distanceEdit;
  
  Wt::WLabel *m_activityLabel;
  
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
  
  Wt::Signal<> m_updated;
};//MakeDrfSrcDef

#endif //MakeDrfSrcDef_h
