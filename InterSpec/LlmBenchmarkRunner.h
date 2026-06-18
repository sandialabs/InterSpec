#ifndef LlmBenchmarkRunner_h
#define LlmBenchmarkRunner_h
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

#include <set>
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include <optional>
#include <functional>

#include <Wt/WSignal>

#include "SpecUtils/SpecFile.h"

// Forward declarations
class InterSpec;
class LlmToolGui;
class LlmInterface;
class LlmConfig;
struct LlmInteraction;


/** Data structures for benchmark problems and results. */

struct BenchmarkCriterion
{
  int points = 0;
  std::string description;
};//struct BenchmarkCriterion


struct BenchmarkJudgement
{
  enum class Type
  {
    LlmJudge,       // Use an LLM to judge the answer
    ContainsString   // Simple string containment check
  };

  Type type = Type::LlmJudge;
  bool caseSensitive = true;  // Only for ContainsString
  std::string expectedAnswer;
  std::vector<BenchmarkCriterion> rubric;
};//struct BenchmarkJudgement


struct BenchmarkQuestion
{
  int part = 1;
  bool skip = false;  // If true, skip this question (e.g., requires capabilities we don't have)
  std::string skipReason;  // Why the question is skipped
  std::string prompt;
  std::string imagePath;  // Optional image file (resolved relative to benchmark XML)
  std::optional<BenchmarkJudgement> judgement;  // Absent means ungraded
};//struct BenchmarkQuestion


struct BenchmarkSpectrumFile
{
  SpecUtils::SpectrumType specType = SpecUtils::SpectrumType::Foreground;
  std::string filePath;  // Resolved relative to benchmark XML
};//struct BenchmarkSpectrumFile


struct BenchmarkSequenceStep
{
  std::vector<BenchmarkSpectrumFile> spectra;
  std::string prompt;
};//struct BenchmarkSequenceStep


struct BenchmarkProblem
{
  std::string id;
  std::string category;
  std::string difficulty;
  std::string context;  // Optional context prepended to first question

  // Simple case: single spectrum load for all questions in this problem
  std::vector<BenchmarkSpectrumFile> spectrumFiles;

  // Multi-spectrum case: sequential steps building up conversation context
  std::vector<BenchmarkSequenceStep> sequenceSteps;

  // Judged questions (asked after any sequence steps complete)
  std::vector<BenchmarkQuestion> questions;
};//struct BenchmarkProblem


struct BenchmarkQuestionResult
{
  std::string problemId;
  int questionPart = 0;
  std::string prompt;
  std::string llmAnswer;
  std::string expectedAnswer;
  bool graded = false;    // false if no judgement was specified
  bool correct = false;
  int score = 0;
  int maxScore = 0;
  std::string judgementReason;
  std::optional<size_t> promptTokens;
  std::optional<size_t> completionTokens;
  std::chrono::milliseconds duration{0};
  bool hadError = false;
  std::string errorMessage;
};//struct BenchmarkQuestionResult


struct BenchmarkResults
{
  std::string benchmarkName;
  std::string modelName;
  std::string benchmarkFilePath;
  std::chrono::system_clock::time_point startTime;
  std::chrono::system_clock::time_point endTime;
  std::vector<BenchmarkQuestionResult> questionResults;

  // Summary stats (computed from questionResults)
  int totalQuestions = 0;
  int gradedQuestions = 0;
  int correctCount = 0;
  int errorCount = 0;
  int totalScore = 0;
  int maxTotalScore = 0;
};//struct BenchmarkResults


/** Drives LLM benchmark evaluation through the LlmToolGui.

 This class parses a benchmark XML file, loops over problems, loads spectra,
 sends questions through the LLM chat GUI (as if the user typed them),
 extracts answers, judges them against expected answers, and produces a report.

 The runner is owned by LlmToolGui and uses an async state-machine pattern
 driven by Wt signals — no new threads are needed.
 */
class LlmBenchmarkRunner
{
public:
  /** Construct a benchmark runner.
   @param viewer The InterSpec instance
   @param toolGui The LlmToolGui that owns this runner
   */
  LlmBenchmarkRunner( InterSpec *viewer, LlmToolGui *toolGui );

  ~LlmBenchmarkRunner();

  /** Parse and start running a benchmark file.
   @param xmlFilePath Absolute path to the benchmark XML file
   */
  void startBenchmark( const std::string &xmlFilePath );

  /** Check if a benchmark is currently running. */
  bool isRunning() const;

  /** Cancel a running benchmark. */
  void cancelBenchmark();

  /** Signal emitted when the benchmark finishes (or is cancelled). */
  Wt::Signal<> &benchmarkFinished();

  /** Parse a benchmark XML file into problems.
   @param xmlPath Absolute path to the benchmark XML file
   @return Parsed problems
   @throws std::runtime_error on parse errors
   */
  static std::vector<BenchmarkProblem> parseXml( const std::string &xmlPath );

  /** Scan a directory for benchmark XML files matching *_llm_benchmark.xml.
   @param directory The directory to scan
   @return Vector of (display_name, absolute_path) pairs
   */
  static std::vector<std::pair<std::string,std::string>> findBenchmarkFiles( const std::string &directory );

private:
  enum class State
  {
    Idle,
    LoadingSpectrum,
    LoadingSequenceStep,
    WaitingForResponse,
    WaitingForSequenceStepResponse,
    Judging,
    Finished
  };

  InterSpec *m_viewer;
  LlmToolGui *m_toolGui;

  State m_state;
  std::vector<BenchmarkProblem> m_problems;
  BenchmarkResults m_results;
  std::string m_benchmarkBaseDir;  // Directory of the benchmark XML file

  size_t m_currentProblem;
  size_t m_currentQuestion;
  size_t m_currentSequenceStep;

  // Timing for current question
  std::chrono::system_clock::time_point m_questionStartTime;

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  std::unique_ptr<std::ostream> m_log_file;
#endif

  // Judging interface (separate from main, with tool calls blocked)
  std::shared_ptr<LlmInterface> m_judgeInterface;

  // Signal connections we manage
  boost::signals2::connection m_conversationFinishedConn;
  boost::signals2::connection m_responseErrorConn;
  boost::signals2::connection m_spectrumChangedConn;  // For sequence step spectrum loads
  boost::signals2::connection m_judgeConversationFinishedConn;
  boost::signals2::connection m_judgeResponseErrorConn;

  Wt::Signal<> m_benchmarkFinished;

  // Whether handleSpectrumChanged is currently disconnected (during SpectrumSequence)
  bool m_spectrumChangeDisconnected;

  /** Advance the state machine to the next step. */
  void advanceToNextStep();

  /** Load spectrum files for the current problem. */
  void loadCurrentSpectrum();

  /** Load spectrum for the current sequence step. */
  void loadSequenceStep();

  /** Called when spectrum loading finishes (for simple problems). */
  void handleSpectrumLoaded( SpecUtils::SpectrumType specType,
                             std::shared_ptr<SpecMeas> meas,
                             std::set<int> samples,
                             std::vector<std::string> detectors );

  /** Send the current judged question to the LLM. */
  void sendCurrentQuestion();

  /** Send the current sequence step prompt to the LLM. */
  void sendSequenceStepPrompt();

  /** Called when the LLM finishes responding. */
  void handleConversationFinished();

  /** Called when the LLM encounters an error. */
  void handleResponseError();

  /** Extract the LLM's answer text and dispatch judging. */
  void extractAnswerAndJudge();

  /** Judge via simple string containment. */
  void judgeWithStringMatch( const std::string &llmAnswer,
                             const BenchmarkJudgement &judgement );

  /** Judge via a separate LLM call. */
  void judgeWithLlm( const std::string &llmAnswer,
                     const BenchmarkQuestion &question );

  /** Called when the judge LLM finishes. */
  void handleJudgeConversationFinished();

  /** Called when the judge LLM encounters an error. */
  void handleJudgeResponseError();

  /** Record a result and advance. */
  void recordResult( BenchmarkQuestionResult result );

  /** Generate the final report (console + JSON). */
  void generateReport();

  /** Show results summary in a dialog and offer save. */
  void showResultsDialog();

  /** Disconnect LlmToolGui::handleSpectrumChanged for SpectrumSequence. */
  void disconnectSpectrumChangedHandler();

  /** Reconnect LlmToolGui::handleSpectrumChanged after SpectrumSequence. */
  void reconnectSpectrumChangedHandler();

  /** Create the judge LLM interface. */
  void createJudgeInterface();

  /** Get the current problem, or nullptr if done. */
  const BenchmarkProblem *currentProblem() const;

  /** Get the current question, or nullptr if done. */
  const BenchmarkQuestion *currentQuestion() const;

  /** Log a message with [LlmBenchmark] prefix. */
  void log( const std::string &message ) const;

  /** Log an error with [LlmBenchmark] prefix. */
  void logError( const std::string &message ) const;

  /** Resolve a relative path against the benchmark base directory. */
  std::string resolvePath( const std::string &relativePath ) const;

  /** Compute summary statistics from questionResults. */
  void computeSummaryStats();
};//class LlmBenchmarkRunner


#endif // USE_LLM_INTERFACE
#endif // LlmBenchmarkRunner_h
