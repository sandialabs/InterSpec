#ifndef McpSessionManager_h
#define McpSessionManager_h
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

#include <map>
#include <set>
#include <mutex>
#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <optional>
#include <functional>

#include "InterSpec/LlmConfig.h"

// Compiled into the MCP HTTP server (target/mcp_http), and into the unit tests of its helpers.
static_assert( USE_MCP_SERVER, "McpSessionManager is only for builds with the MCP server" );

class CdpBrowser;
class InterSpec;

/** Owns the analysis sessions of InterSpec built as the multi-session MCP HTTP server
 (`BUILD_AS_MCP_HTTP_SERVER`; see target/mcp_http).

 Each session is an ordinary Wt session of InterSpec, driven by its own tab of one headless
 Chromium (a Wt session only lives while a browser keeps it alive, and timers, server-push, and
 chart rendering all need one).  A session is identified to MCP clients by its "handle", which is
 also the Wt external token the tab was opened with (see `InterSpecServer::add_allowed_session_token`);
 the internal Wt session id is never handed out.

 Files are exchanged through a per-session directory, `<ExchangeDirectory>/<YYYY-MM-DD>_<handle>/`,
 holding `input/`, `output/`, `logs/`, and a marker file; clients may use the directory directly
 (a host mount), or through an upload URL for the session and a download URL per output file - each
 a random token, separate from the handle, and revoked when the session ends.  Every path a session
 uses is resolved inside its own directory (`resolveSessionPath`).

 Sessions end on `endSession` (the `close_session` tool), the idle timeout, the maximum lifetime,
 their app going away, or the browser dying; their directory is then deleted unless
 `KeepSessionFiles` is set.  At startup, session directories left behind by an earlier run are
 removed the same way.

 Thread safety: all public functions may be called from any thread.
 */
class McpSessionManager : public std::enable_shared_from_this<McpSessionManager>
{
public:
  using AgentSettings = LlmConfig::McpServer::AgentServer;

  enum class State
  {
    Starting,  ///< Directories made and tab requested; waiting for the InterSpecApp to start.
    Live,      ///< The app has started, and can be used.
    Ended      ///< Ended (see `SessionInfo::end_reason`); kept a while to give clients a clear error.
  };

  struct SessionInfo
  {
    std::string handle;
    State state = State::Starting;
    std::string end_reason;

    /** Directory name relative to the exchange root, e.g. "2026-10-17_Xy3...". */
    std::string dir_name;

    /** Absolute, canonical path of the session directory. */
    std::string dir;

    std::string upload_token;

    /** Empty until the app has started. */
    std::string wt_session_id;

    std::chrono::system_clock::time_point created;
    std::chrono::system_clock::time_point expires;
  };//struct SessionInfo

  /** What a download or upload token refers to. */
  struct FileToken
  {
    std::string handle;
    bool is_upload = false;
    /** For downloads, the path relative to the session directory (e.g. "output/x.n42"). */
    std::string relpath;
  };//struct FileToken

  /** Creates the manager, sweeps session directories left by an earlier run, and starts the reaper.
   The browser is started lazily, by the first `createSession`.

   @param settings  The `<AgentServer>` settings.
   @param app_host  Host the browser reaches the server at, e.g. "127.0.0.1".
   @param app_port  Port the server listens on.

   @throws std::exception if the exchange directory can not be created or written to.
   */
  static void start( const AgentSettings &settings, const std::string &app_host, const int app_port );

  /** Ends all sessions, and stops the browser.  Called as the server shuts down. */
  static void stop();

  /** The manager, or nullptr if not running as an agent server. */
  static std::shared_ptr<McpSessionManager> instance();

  ~McpSessionManager();

  const AgentSettings &settings() const;

  /** Creates a session: its directories and upload token, authorizes its handle as a Wt external
   token, and asks the browser to open a tab for it (starting the browser if needed).

   The session is `State::Starting` until its app calls `appStarted`; use `addSessionListener` to
   be told when that, or a failure, happens.

   @returns the new session's handle.
   @throws std::exception if `MaxSessions` are already open, or setup failed.
   */
  std::string createSession();

  /** Information on a session, including recently ended ones; nullopt if the handle is unknown. */
  std::optional<SessionInfo> session( const std::string &handle ) const;

  /** Ends a session: asks its app to quit, closes its tab, revokes its tokens, and (once the app
   has gone, or a short grace period passed) deletes its directory unless `KeepSessionFiles`.
   Idempotent.
   */
  void endSession( const std::string &handle, const std::string &reason );

  /** Records client activity, resetting the idle timeout. */
  void touch( const std::string &handle );

  /** Bracket a tool call (or assistant turn) running in the session: a session with calls in
   progress is not idle, however long they take (the maximum lifetime still applies).
   */
  void callStarted( const std::string &handle );
  void callFinished( const std::string &handle );

  /** Calls `fn` (once, on an arbitrary thread) the next time the session's state changes; or right
   away if it is no longer `State::Starting`.  Used to wait for a new session to come up.
   */
  void addSessionListener( const std::string &handle, std::function<void()> fn );

  /** Called by the InterSpecApp of a managed session once its widgets are set up; only the first
   app to start is the session's.
   */
  void appStarted( const std::string &handle, const std::string &wt_session_id );

  /** Called from the destructor of the InterSpecApp of a managed session; ignored if
   `wt_session_id` is not the session's app.
   */
  void appEnded( const std::string &handle, const std::string &wt_session_id );

  /** True if `handle` is (or was) a session of this manager. */
  bool isManagedHandle( const std::string &handle ) const;

  /** Resolves `relpath` (e.g. "input/spec.n42") inside the session's directory.

   The result is the canonical absolute path, with symbolic links resolved; it need not exist.

   @throws std::runtime_error, with a deliberately generic message, if `relpath` is empty or
           absolute, or resolves outside the session directory, or the session is not live.
   */
  std::string resolveSessionPath( const std::string &handle, const std::string &relpath ) const;

  /** Resolves `relpath` inside `base_dir` (which must be canonical): the canonical absolute path,
   with symbolic links resolved (it need not exist).  The building block of `resolveSessionPath`.

   @throws std::runtime_error, with a deliberately generic message, if `relpath` is empty or
           absolute, or resolves outside `base_dir`.
   */
  static std::string resolvePathWithin( const std::string &base_dir, const std::string &relpath );

  /** Returns a name safe to create in a session directory: the final path component only, with
   characters outside [A-Za-z0-9._-] replaced by '_', no leading '.', and at most 128 characters.
   Returns "file" if nothing usable is left.
   */
  static std::string sanitizeFilename( const std::string &name );

  /** Picks a name, based on `filename` (sanitized), that does not yet exist in `directory`,
   appending "-1", "-2", ... before the extension as needed.
   */
  static std::string uniqueFilename( const std::string &directory, const std::string &filename );

  /** Total size of the files in the session's directory. */
  uint64_t sessionBytesUsed( const std::string &handle ) const;

  /** Registers a download token for a file of the session (`relpath` relative to the session
   directory, e.g. "output/x.n42"), returning the token; an existing token for the file is reused.
   */
  std::string addDownload( const std::string &handle, const std::string &relpath );

  /** What `token` refers to, if it is a current upload or download token. */
  std::optional<FileToken> lookupFileToken( const std::string &token ) const;

  /** Session counts and browser state, for the unauthenticated health check (no identifiers). */
  struct HealthInfo
  {
    size_t live_sessions = 0;
    size_t starting_sessions = 0;
    size_t max_sessions = 0;
    bool browser_running = false;
  };
  HealthInfo health() const;

  /** The handle of the session the current WApplication belongs to, or empty if the current
   application is not a managed session (or there is none).  Must be called from a session thread.
   */
  static std::string currentHandle();

private:
  McpSessionManager( const AgentSettings &settings, const std::string &app_host, const int app_port );

  /** The URL a new session's tab is opened at. */
  std::string sessionAppUrl( const std::string &handle );

  struct Session
  {
    SessionInfo info;
    std::chrono::steady_clock::time_point created;
    std::chrono::steady_clock::time_point last_activity;
    std::chrono::steady_clock::time_point ended;
    std::string cdp_target_id;
    bool app_ended = false;
    bool dir_removed = false;
    int active_calls = 0;
    /** The `m_browser_generation` of the browser the session's tab is in. */
    uint64_t browser_generation = 0;
    std::set<std::string> download_tokens;
    std::vector<std::function<void()>> listeners;
  };//struct Session

  /** Sweeps directories left behind by an earlier run (only marker-bearing session directories). */
  void sweepOrphanedDirectories();

  /** Deletes the session directory, if allowed; caller must hold `m_mutex`. */
  void removeSessionDirectory( Session &session );

  /** Periodic check for timed out sessions, sessions that never started, and a dead browser. */
  void reap();
  void scheduleReaper();

  /** Called by the browser when a tab goes away on its own (renderer crash, etc). */
  void handleTargetGone( const std::string &target_id );

  /** Called by the browser of `browser_generation` when its process exits; ends its sessions. */
  void handleBrowserExited( uint64_t browser_generation );

  /** Ends a session; caller must hold `m_mutex`.  Returns the listeners to call, once unlocked. */
  std::vector<std::function<void()>> endSessionLocked( Session &session, const std::string &reason );

  const AgentSettings m_settings;
  const std::string m_app_host;
  const int m_app_port;

  /** Numbers the per-session host names (see `sessionAppUrl`); protected by `m_mutex`. */
  uint64_t m_num_sessions_created = 0;

  /** Canonical absolute path of the exchange directory. */
  std::string m_exchange_root;

  mutable std::mutex m_mutex;
  std::map<std::string,Session> m_sessions;              // keyed by handle
  std::map<std::string,FileToken> m_file_tokens;         // keyed by token
  std::shared_ptr<CdpBrowser> m_browser;

  /** Incremented for each browser started; protected by `m_mutex`. */
  uint64_t m_browser_generation = 0;

  /** Set false by `stop()`; protected by `m_mutex`. */
  bool m_running;
};//class McpSessionManager

#endif //McpSessionManager_h
