#ifndef InterSpec_config_h
#define InterSpec_config_h
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
/* Standalone configuration for the InterSpec-light hard fork of the InterSpec
 peak-fitting code.  Only the macros the forked sources actually use are defined.
 */

#include <cmath>
#include <iostream>

#define PERFORM_DEVELOPER_CHECKS 0
#define USE_REL_ACT_TOOL 1
#define SpecUtils_ENABLE_D3_CHART 1
#define INTERSPEC_CERES_JET_NUMTRAITS_HAS_NONFINITE 0

#ifndef IsInf
#define IsInf(x) (std::isinf)(x)
#endif

#ifndef IsNan
#define IsNan(x) (std::isnan)(x)
#endif

#define passMessage(message,priority) \
  { std::cerr << "passMessage: " << message << std::endl; }

#define InterSpec_API

#endif // InterSpec_config_h
