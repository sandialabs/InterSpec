#ifndef QuickJsSandbox_h
#define QuickJsSandbox_h
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

#include <string>
#include <variant>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"


/** Runs the assistant's `run_javascript` code with QuickJS (quickjs-ng), for where there is no
 browser to run it in a sandboxed iframe (the stdio MCP server).

 Gives the same results as the browser sandbox (`InterSpec_resources/LlmJsSandbox.js` and
 `LlmJsSandboxBridge`): the same wrapping of the code, `data` argument, console capture, and result
 serialization, run by the same JavaScript.  Each call gets a fresh runtime, with limits on memory,
 stack, and time; only the core language is available - no modules, files, network, or timers.
 */
namespace QuickJsSandbox
{
  /** Memory and stack limits of each run. */
  constexpr size_t sm_memory_limit_bytes = 64 * 1024 * 1024;
  constexpr size_t sm_max_stack_bytes = 256 * 1024;  // well within the 512 kB of a macOS secondary thread

  /** Runs the tool call; `params` has `code` (required), `data` (optional, any JSON), and
   `timeoutMs` (optional; default 5000, at most 30000), as for the browser sandbox.

   Blocks for up to the timeout, so should be called on a worker thread; it is independent of Wt.
   Returns the tool result (which reports errors in the code itself), or an error string if the
   call could not be made at all.
   */
  std::variant<nlohmann::json, std::string> run( const nlohmann::json &params );
}//namespace QuickJsSandbox

#endif //QuickJsSandbox_h
