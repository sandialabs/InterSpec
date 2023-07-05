#ifndef UnitsConverterTool_h
#define UnitsConverterTool_h
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
  class WDoubleSpinBox;
}//namespace Wt


class UnitsConverterTool : public AuxWindow
{
public:
  UnitsConverterTool();
  virtual ~UnitsConverterTool();
  
  /** Performs the actual conversion.
   
   Throws exception if a sucessful conversion couldnt be done.
   */
  static std::string convert( std::string input );
  
  /** Handles receiving a "deep-link" url starting with "interspec://unit?input=1.2m".
   
   Example URIs:
   - "interspec://unit?input=1.2m"
   
   @param query_str The query portion of the URI.  So for example, if the URI has a value of
          "interspec://unit?input=1.2m", then this string would be "input=1.2m".
          Capitalization is not important.
          Assumes the string passed in has already been url-decoded.
          If not a valid query_str, throws exception.
   */
  void handleAppUrl( std::string query_str );
  
  /** Encodes current tool state to app-url format.  Returned string does not include the
   "interspec://" protocol, or "unit" authority; so will look something like "input=1.2m",
   The path part of the URI specifies tab the tool is on.
   and it will not be url-encoded.
   */
  std::string encodeStateToUrl() const;
  
protected:
  void convert();
  
  Wt::WLineEdit *m_input;
  Wt::WLineEdit *m_output;
  Wt::WText *m_message;
  Wt::WString m_prevInput;  // For tracking undo/redo
  std::string m_prevAnswer; // For tracking undo/redo
};//class UnitsConverterTool


#endif //UnitsConverterTool_h
