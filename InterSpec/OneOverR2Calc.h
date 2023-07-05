#ifndef OneOverR2Calc_h
#define OneOverR2Calc_h
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

#include "InterSpec/AuxWindow.h"

namespace Wt
{
  class WText;
  class WLineEdit;
  class WComboBox;
  class WDoubleSpinBox;
}//namespace Wt


class OneOverR2Calc : public AuxWindow
{
public:
  OneOverR2Calc();
  virtual ~OneOverR2Calc();

  void doCalc();

  
  /** Handles receiving a "deep-link" url starting with "interspec://1overr2/...".
   
   Example URIs:
   - "interspec://1overr2?near=3.5&far=1.1&back=2&dist=1.9&power=1.75"
   
   @param query_str The query portion of the URI.  So for example, if the URI has a value of
          "interspec://1overr2?near=3.5&far=1.1&b...", then this string would be "near=3.5&far=1.1&b...".
          This string is is in standard URL format of "key1=value1&key2=value2&..." with ordering not mattering.
          Capitalization is not important.
          Assumes the string passed in has alaready been url-decoded.
          If not a valid query_str, throws exception.
   */
  void handleAppUrl( std::string query_str );
  
  /** Encodes current tool state to app-url format.  Returned string does not include the
   "interspec://" protocol, or "1overr2" path; so will look something like "near=3.5&far=1.1&back=2&dit=1.9&power=1.75",
   and it will not be url-encoded.
   */
  std::string encodeStateToUrl() const;
  
protected:
  void powerLawSelected();
  
protected:
  Wt::WDoubleSpinBox *m_nearMeasurement;
  Wt::WDoubleSpinBox *m_farMeasurement;
  Wt::WDoubleSpinBox *m_backgroundMeasurment;
  Wt::WDoubleSpinBox *m_distance;

  Wt::WComboBox *m_powerLawSelect;
  
  Wt::WLineEdit *m_answer;
  Wt::WText *m_message;
  
  /** For tracking of undo/redo, we will store all all the values, and power law index (as a float). */
  std::array<float,5> m_prevValues;
};//class OneOverR2Calc


#endif //OneOverR2Calc
