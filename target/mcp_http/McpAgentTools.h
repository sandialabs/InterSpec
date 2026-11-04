#ifndef McpAgentTools_h
#define McpAgentTools_h
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
#include <vector>
#include <cstdint>
#include <iosfwd>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

static_assert( BUILD_AS_MCP_HTTP_SERVER, "McpAgentTools is only for the MCP HTTP server build (target/mcp_http)" );

class InterSpec;

/** The MCP tools, and supporting functions, offered only when InterSpec runs as a multi-session
 agent server (see `McpSessionManager`): session management, file exchange, and exporting the
 analysis.  `LlmMcpResource` dispatches them; the tools shared with the
 in-app assistant are in `LlmTools::ToolRegistry`.
 */
namespace McpAgentTools
{
  struct ToolDef
  {
    std::string name;
    std::string description;
    nlohmann::json schema;

    /** True if the tool acts on a session (and so needs the `userSession` argument). */
    bool needs_session = false;

    /** True if the tool runs inside the session (on its thread); else it is handled by the resource. */
    bool runs_in_session = false;
  };//struct ToolDef

  /** The agent-server-only tools. */
  const std::vector<ToolDef> &toolDefinitions();

  /** The definition of `name`, or nullptr if it is not an agent-server tool. */
  const ToolDef *toolDefinition( const std::string &name );

  /** The `instructions` returned by MCP `initialize`, from "mcp_agent_instructions.md" in the LLM
   config directory, or else the static data directory.  Empty if neither has it.
   */
  std::string serverInstructions();

  /** The client skill: "skills/interspec-mcp-agent/SKILL.md", from the LLM config directory or
   the static data directory.
   */
  std::string usageGuide();

  /** Writes `in` to a new file in the session's `input/` directory, named from `filename`
   (sanitized, and made unique).  Enforces `MaxUploadBytes` and `MaxSessionBytes`.

   @param num_bytes Set to the number of bytes written.
   @returns the path relative to the session directory, e.g. "input/spectrum.n42".
   @throws std::exception if a limit was exceeded, or the file could not be written.
   */
  std::string saveInputFile( const std::string &handle, const std::string &filename,
                             std::istream &in, uint64_t &num_bytes );

  /** Writes `data` to a new file in the session's `output/` directory, as for `saveInputFile`.

   @returns the path relative to the session directory, e.g. "output/spectrum.png".
   */
  std::string saveOutputFile( const std::string &handle, const std::string &filename,
                              const std::string &data );

  /** Lists the files of the session's `input/` and `output/` directories:
   `{input: [{path, bytes}], output: [{path, bytes}]}`.
   */
  nlohmann::json listSessionFiles( const std::string &handle );

  /** The `export_analysis` tool: writes the requested items to the session's `output/` directory.

   Must be called from the session's thread.  Each item is one of
   - "n42": the foreground as an N42 file with all the analysis (peaks, the activity/shielding
     model, isotopics configuration, the assistant conversation, ...), and the background, if it is
     from a different file; a later session can resume from these with `load_spectrum_file`.
   - "peaks_csv": the fitted peaks as CSV.
   - "spectrum_html": a standalone, interactive HTML plot of the displayed spectra and peaks.

   @returns `{files: [{path, item, bytes}]}`, with paths relative to the session directory.
   */
  nlohmann::json exportAnalysis( InterSpec *viewer, const std::string &handle,
                                 const nlohmann::json &args );
}//namespace McpAgentTools

#endif //McpAgentTools_h
