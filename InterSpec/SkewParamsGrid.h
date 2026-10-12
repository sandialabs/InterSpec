#ifndef SkewParamsGrid_h
#define SkewParamsGrid_h
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
#include <optional>

#include <Wt/WSignal.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/PeakDef.h"

class NativeFloatSpinBox;

namespace Wt
{
  class WText;
  class WCheckBox;
}//namespace Wt


/** A grid to view or edit the parameters of a peak skew type: a row for each parameter, with its name,
 its value (at the low and high energy of the spectrum, for parameters that may vary with energy - see
 `PeakDef::is_energy_dependent`), and optionally a checkbox of whether to fit it, and its fit result.

 Used by the peak-fit preferences, the fit-skew tool, and the isotopics-by-nuclides tool.
 */
class SkewParamsGrid : public Wt::WContainerWidget
{
public:
  /** @param show_fit_column Add a checkbox, for each parameter, of whether it should be fit.
      @param allow_blank Allow a value to be left blank (meaning no value is given).
   */
  SkewParamsGrid( const bool show_fit_column, const bool allow_blank );

  /** Rebuilds the rows for `type`, with each value at its default starting value (see
   `PeakDef::skew_parameter_range`), and each fit checkbox per `PeakDef::skew_parameter_fit_by_default`.
   The widget is hidden for skew types without parameters.
   */
  void setSkewType( const PeakDef::SkewType type );

  PeakDef::SkewType skewType() const;

  /** Sets a parameter's value; `upper` is only used by energy-dependent parameters.  A value that is not
   given is left blank if blanks are allowed; otherwise an upper value that is not given is set to `lower`.
   */
  void setValue( const size_t index, const std::optional<double> &lower, const std::optional<double> &upper );

  /** The parameter's value; nullopt if it is blank, or not a parameter of the skew type.  The upper value
   is nullopt for parameters that are not energy dependent.
   */
  std::optional<double> lowerValue( const size_t index ) const;
  std::optional<double> upperValue( const size_t index ) const;

  bool isFit( const size_t index ) const;
  void setFit( const size_t index, const bool fit );

  /** Shows each parameter's fit result in an extra column (requires the fit column), with an optional
   tooltip for each; an empty `results` removes the column.
   */
  void setResults( const std::vector<Wt::WString> &results, const std::vector<Wt::WString> &tooltips = {} );

  void setEditable( const bool editable );

  /** Emitted when the user changes a value or fit checkbox (not for programmatic changes). */
  Wt::Signal<> &userChanged();

protected:
  const bool m_show_fit;
  const bool m_allow_blank;
  bool m_editable;
  PeakDef::SkewType m_skew_type;

  NativeFloatSpinBox *m_lower[6];
  NativeFloatSpinBox *m_upper[6];
  Wt::WCheckBox *m_fit[6];
  Wt::WText *m_result[6];
  Wt::WText *m_result_header;

  Wt::Signal<> m_user_changed;
};//class SkewParamsGrid

#endif //SkewParamsGrid_h
