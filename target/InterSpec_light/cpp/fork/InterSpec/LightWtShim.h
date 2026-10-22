#ifndef LightWtShim_h
#define LightWtShim_h
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

/* Minimal stand-ins for the two Wt types the forked peak-fitting code uses, so the
 forked sources only need their #include lines changed.
 */

#include <string>

namespace Wt
{
  enum NoFlagsTag { None = 0 };

  /** Bit-flag set over a plain enum; `test()` is a bitwise AND, like Wt's. */
  template<class Enum>
  class WFlags
  {
  public:
    WFlags() : m_flags( 0 ) {}
    WFlags( NoFlagsTag ) : m_flags( 0 ) {}
    WFlags( const Enum flag ) : m_flags( static_cast<int>(flag) ) {}
    explicit WFlags( const int flags ) : m_flags( flags ) {}

    bool test( const Enum flag ) const { return (m_flags & static_cast<int>(flag)) != 0; }
    WFlags &operator|=( const Enum flag ) { m_flags |= static_cast<int>(flag); return *this; }
    WFlags &operator|=( const WFlags &rhs ) { m_flags |= rhs.m_flags; return *this; }
    WFlags operator|( const Enum flag ) const { WFlags r( *this ); r |= flag; return r; }
    WFlags &clear( const Enum flag ) { m_flags &= ~static_cast<int>(flag); return *this; }
    int value() const { return m_flags; }

  private:
    int m_flags;
  };//class WFlags


  /** A CSS color string; empty means "not set". */
  class WColor
  {
  public:
    WColor() = default;
    explicit WColor( const std::string &css ) : m_css( css ) {}

    bool isDefault() const { return m_css.empty(); }
    const std::string &cssText( const bool /*withAlpha*/ = false ) const { return m_css; }

    bool operator==( const WColor &rhs ) const { return m_css == rhs.m_css; }
    bool operator!=( const WColor &rhs ) const { return m_css != rhs.m_css; }

  private:
    std::string m_css;
  };//class WColor
}//namespace Wt

#endif //LightWtShim_h
