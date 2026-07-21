#ifndef LlmJsSandboxBridge_h
#define LlmJsSandboxBridge_h
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

#if( USE_LLM_INTERFACE )

#include <map>
#include <string>
#include <variant>
#include <functional>

#include <Wt/WSignal>
#include <Wt/WContainerWidget>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

/** Bridge between the C++ `run_javascript` LLM tool and the browser-side
 sandboxed iframe defined in `InterSpec_resources/LlmJsSandbox.js`.

 One bridge instance lives as a hidden child container of `LlmToolGui`;
 it owns the JSignal that delivers iframe results back to C++ and tracks
 in-flight request IDs so that asynchronous callbacks can be matched up
 with their original requests.

 The iframe itself is created (and re-created after a timeout-kill) on
 the JavaScript side; this class never directly manipulates DOM beyond
 hosting the container element the JS code attaches the iframe to.
 */
class LlmJsSandboxBridge : public Wt::WContainerWidget
{
public:
  /** Callback fired on the GUI thread when a sandboxed JS run completes
   (successfully, with an error, or via timeout).  The variant holds a
   JSON result on success or an error string on hard dispatch failure. */
  using AsyncCallback
    = std::function<void( std::variant<nlohmann::json, std::string> )>;

  LlmJsSandboxBridge( Wt::WContainerWidget *parent = nullptr );
  ~LlmJsSandboxBridge();

  /** Run a sandboxed JavaScript snippet.

   `params` is the raw tool-parameter object - expected keys are:
     - `code` (string, required): JS source.
     - `data` (any JSON, optional): exposed inside the sandbox as `data`.
     - `timeoutMs` (integer, optional, default 5000, max 30000).

   The callback is invoked exactly once on the GUI thread.
   */
  void run( const nlohmann::json &params, AsyncCallback callback );

private:
  /** Slot connected to the `jsSandboxResult` JSignal.  Parses the JSON
   payload, looks up the matching pending callback, and dispatches. */
  void handleJsResult( std::string payload );

  /** Wt-side signal whose JS counterpart is `Wt.emit(..., 'jsSandboxResult', payload)`. */
  Wt::JSignal<std::string> m_resultSignal;

  /** Pending request IDs awaiting callbacks.  Stale IDs (e.g. results
   that arrive after the JS-side timeout already fired) are dropped. */
  std::map<std::string, AsyncCallback> m_pending;

  /** Monotonic counter used to keep request IDs unique within a session. */
  uint64_t m_nextRequestId;
};

#endif // USE_LLM_INTERFACE
#endif // LlmJsSandboxBridge_h
