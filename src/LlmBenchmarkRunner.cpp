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
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include <Wt/WText>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include <rapidxml/rapidxml.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmToolGui.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/LlmBenchmarkRunner.h"
#include "InterSpec/LlmConversationHistory.h"

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

using namespace std;


LlmBenchmarkRunner::LlmBenchmarkRunner( InterSpec *viewer, LlmToolGui *toolGui )
  : m_viewer( viewer ),
    m_toolGui( toolGui ),
    m_state( State::Idle ),
    m_currentProblem( 0 ),
    m_currentQuestion( 0 ),
    m_currentSequenceStep( 0 ),
    m_spectrumChangeDisconnected( false )
{
  assert( m_viewer );
  assert( m_toolGui );
}


LlmBenchmarkRunner::~LlmBenchmarkRunner()
{
  // Disconnect any active signal connections
  m_conversationFinishedConn.disconnect();
  m_responseErrorConn.disconnect();
  m_spectrumChangedConn.disconnect();
  m_judgeConversationFinishedConn.disconnect();
  m_judgeResponseErrorConn.disconnect();
}


bool LlmBenchmarkRunner::isRunning() const
{
  return (m_state != State::Idle) && (m_state != State::Finished);
}


void LlmBenchmarkRunner::cancelBenchmark()
{
  if( !isRunning() )
    return;

  log( "Benchmark cancelled by user" );

  // Reconnect spectrum changed handler if we disconnected it
  reconnectSpectrumChangedHandler();

  m_state = State::Finished;
  computeSummaryStats();

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  m_log_file.reset();
#endif

  m_benchmarkFinished.emit();
}


Wt::Signal<> &LlmBenchmarkRunner::benchmarkFinished()
{
  return m_benchmarkFinished;
}


void LlmBenchmarkRunner::log( const string &message ) const
{
  cout << "[LlmBenchmark] " << message << endl;

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  if( m_log_file )
  {
    (*m_log_file) << message << endl;
    m_log_file->flush();
  }
#endif
}


void LlmBenchmarkRunner::logError( const string &message ) const
{
  cerr << "[LlmBenchmark] ERROR: " << message << endl;

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  if( m_log_file )
  {
    (*m_log_file) << "[LlmBenchmark] ERROR: " << message << endl;
    m_log_file->flush();
  }
#endif
}


string LlmBenchmarkRunner::resolvePath( const string &relativePath ) const
{
  if( relativePath.empty() )
    return string();

  // If already absolute, return as-is
  if( !relativePath.empty() && (relativePath[0] == '/') )
    return relativePath;

#ifdef _WIN32
  if( relativePath.size() >= 2 && relativePath[1] == ':' )
    return relativePath;
#endif

  return SpecUtils::append_path( m_benchmarkBaseDir, relativePath );
}


const BenchmarkProblem *LlmBenchmarkRunner::currentProblem() const
{
  if( m_currentProblem < m_problems.size() )
    return &m_problems[m_currentProblem];
  return nullptr;
}


const BenchmarkQuestion *LlmBenchmarkRunner::currentQuestion() const
{
  const BenchmarkProblem *problem = currentProblem();
  if( problem && (m_currentQuestion < problem->questions.size()) )
    return &problem->questions[m_currentQuestion];
  return nullptr;
}


// ---- XML Parsing ----

vector<BenchmarkProblem> LlmBenchmarkRunner::parseXml( const string &xmlPath )
{
  vector<BenchmarkProblem> problems;

  const string xmlDir = SpecUtils::parent_path( xmlPath );

  // Read file contents
  vector<char> data;
  SpecUtils::load_file_data( xmlPath.c_str(), data );
  data.push_back( '\0' );

  rapidxml::xml_document<char> doc;
  const int flags = rapidxml::parse_trim_whitespace;
  doc.parse<flags>( data.data() );

  const rapidxml::xml_node<char> *root = XML_FIRST_NODE( &doc, "LlmBenchmark" );
  if( !root )
    throw runtime_error( "LlmBenchmarkRunner::parseXml: Missing <LlmBenchmark> root element in " + xmlPath );

  for( const rapidxml::xml_node<char> *problemNode = XML_FIRST_NODE( root, "Problem" );
       problemNode;
       problemNode = XML_NEXT_TWIN( problemNode ) )
  {
    BenchmarkProblem problem;

    const rapidxml::xml_attribute<char> *idAttr = XML_FIRST_ATTRIB( problemNode, "id" );
    if( idAttr )
      problem.id = SpecUtils::xml_value_str( idAttr );

    // Check for problem-level skip attribute (applies to all questions in this problem)
    bool problemSkip = false;
    string problemSkipReason;
    const rapidxml::xml_attribute<char> *problemSkipAttr = XML_FIRST_ATTRIB( problemNode, "skipQuestion" );
    if( problemSkipAttr )
    {
      const string skipStr = SpecUtils::xml_value_str( problemSkipAttr );
      problemSkip = SpecUtils::iequals_ascii( skipStr, "true" ) || (skipStr == "1");
    }
    const rapidxml::xml_attribute<char> *problemSkipReasonAttr = XML_FIRST_ATTRIB( problemNode, "skipReason" );
    if( problemSkipReasonAttr )
      problemSkipReason = SpecUtils::xml_value_str( problemSkipReasonAttr );

    const rapidxml::xml_node<char> *catNode = XML_FIRST_NODE( problemNode, "Category" );
    if( catNode )
      problem.category = SpecUtils::xml_value_str( catNode );

    const rapidxml::xml_node<char> *diffNode = XML_FIRST_NODE( problemNode, "Difficulty" );
    if( diffNode )
      problem.difficulty = SpecUtils::xml_value_str( diffNode );

    const rapidxml::xml_node<char> *ctxNode = XML_FIRST_NODE( problemNode, "Context" );
    if( ctxNode )
      problem.context = SpecUtils::xml_value_str( ctxNode );

    // Parse SpectrumFile elements (simple case)
    for( const rapidxml::xml_node<char> *sfNode = XML_FIRST_NODE( problemNode, "SpectrumFile" );
         sfNode;
         sfNode = XML_NEXT_TWIN( sfNode ) )
    {
      BenchmarkSpectrumFile sf;
      sf.filePath = SpecUtils::xml_value_str( sfNode );
      SpecUtils::trim( sf.filePath );

      const rapidxml::xml_attribute<char> *typeAttr = XML_FIRST_ATTRIB( sfNode, "type" );
      const string typeStr = typeAttr ? SpecUtils::xml_value_str( typeAttr ) : string( "Foreground" );

      if( SpecUtils::iequals_ascii( typeStr, "Background" ) )
        sf.specType = SpecUtils::SpectrumType::Background;
      else if( SpecUtils::iequals_ascii( typeStr, "Secondary" ) || SpecUtils::iequals_ascii( typeStr, "SecondForeground" ) )
        sf.specType = SpecUtils::SpectrumType::SecondForeground;
      else
        sf.specType = SpecUtils::SpectrumType::Foreground;

      if( !sf.filePath.empty() )
        problem.spectrumFiles.push_back( sf );
    }//for each SpectrumFile

    // Parse SpectrumSequence element (multi-spectrum case)
    const rapidxml::xml_node<char> *seqNode = XML_FIRST_NODE( problemNode, "SpectrumSequence" );
    if( seqNode )
    {
      for( const rapidxml::xml_node<char> *stepNode = XML_FIRST_NODE( seqNode, "Step" );
           stepNode;
           stepNode = XML_NEXT_TWIN( stepNode ) )
      {
        BenchmarkSequenceStep step;

        for( const rapidxml::xml_node<char> *sfNode = XML_FIRST_NODE( stepNode, "SpectrumFile" );
             sfNode;
             sfNode = XML_NEXT_TWIN( sfNode ) )
        {
          BenchmarkSpectrumFile sf;
          sf.filePath = SpecUtils::xml_value_str( sfNode );
          SpecUtils::trim( sf.filePath );

          const rapidxml::xml_attribute<char> *typeAttr = XML_FIRST_ATTRIB( sfNode, "type" );
          const string typeStr = typeAttr ? SpecUtils::xml_value_str( typeAttr ) : string( "Foreground" );

          if( SpecUtils::iequals_ascii( typeStr, "Background" ) )
            sf.specType = SpecUtils::SpectrumType::Background;
          else
            sf.specType = SpecUtils::SpectrumType::Foreground;

          if( !sf.filePath.empty() )
            step.spectra.push_back( sf );
        }//for each SpectrumFile in Step

        const rapidxml::xml_node<char> *promptNode = XML_FIRST_NODE( stepNode, "Prompt" );
        if( promptNode )
          step.prompt = SpecUtils::xml_value_str( promptNode );

        problem.sequenceSteps.push_back( step );
      }//for each Step
    }//if( seqNode )

    // Parse Question elements
    for( const rapidxml::xml_node<char> *qNode = XML_FIRST_NODE( problemNode, "Question" );
         qNode;
         qNode = XML_NEXT_TWIN( qNode ) )
    {
      BenchmarkQuestion question;

      const rapidxml::xml_attribute<char> *partAttr = XML_FIRST_ATTRIB( qNode, "part" );
      if( partAttr )
      {
        const string partStr = SpecUtils::xml_value_str( partAttr );
        question.part = std::stoi( partStr );
      }

      const rapidxml::xml_attribute<char> *skipAttr = XML_FIRST_ATTRIB( qNode, "skipQuestion" );
      if( skipAttr )
      {
        const string skipStr = SpecUtils::xml_value_str( skipAttr );
        question.skip = SpecUtils::iequals_ascii( skipStr, "true" ) || (skipStr == "1");
      }

      const rapidxml::xml_attribute<char> *skipReasonAttr = XML_FIRST_ATTRIB( qNode, "skipReason" );
      if( skipReasonAttr )
        question.skipReason = SpecUtils::xml_value_str( skipReasonAttr );

      const rapidxml::xml_node<char> *promptNode = XML_FIRST_NODE( qNode, "Prompt" );
      if( promptNode )
        question.prompt = SpecUtils::xml_value_str( promptNode );

      const rapidxml::xml_node<char> *imageNode = XML_FIRST_NODE( qNode, "Image" );
      if( imageNode )
        question.imagePath = SpecUtils::xml_value_str( imageNode );

      const rapidxml::xml_node<char> *judgNode = XML_FIRST_NODE( qNode, "Judgement" );
      if( judgNode )
      {
        BenchmarkJudgement judgement;

        const rapidxml::xml_attribute<char> *typeAttr = XML_FIRST_ATTRIB( judgNode, "type" );
        const string typeStr = typeAttr ? SpecUtils::xml_value_str( typeAttr ) : string( "LlmJudge" );

        if( SpecUtils::iequals_ascii( typeStr, "ContainsString" ) )
          judgement.type = BenchmarkJudgement::Type::ContainsString;
        else
          judgement.type = BenchmarkJudgement::Type::LlmJudge;

        const rapidxml::xml_attribute<char> *csAttr = XML_FIRST_ATTRIB( judgNode, "caseSensitive" );
        if( csAttr )
        {
          const string csStr = SpecUtils::xml_value_str( csAttr );
          judgement.caseSensitive = !(SpecUtils::iequals_ascii( csStr, "false" ) || csStr == "0");
        }

        const rapidxml::xml_node<char> *expNode = XML_FIRST_NODE( judgNode, "ExpectedAnswer" );
        if( expNode )
          judgement.expectedAnswer = SpecUtils::xml_value_str( expNode );

        const rapidxml::xml_node<char> *rubricNode = XML_FIRST_NODE( judgNode, "Rubric" );
        if( rubricNode )
        {
          for( const rapidxml::xml_node<char> *critNode = XML_FIRST_NODE( rubricNode, "Criterion" );
               critNode;
               critNode = XML_NEXT_TWIN( critNode ) )
          {
            BenchmarkCriterion criterion;
            criterion.description = SpecUtils::xml_value_str( critNode );

            const rapidxml::xml_attribute<char> *ptsAttr = XML_FIRST_ATTRIB( critNode, "points" );
            if( ptsAttr )
              criterion.points = std::stoi( SpecUtils::xml_value_str( ptsAttr ) );

            judgement.rubric.push_back( criterion );
          }//for each Criterion
        }//if( rubricNode )

        question.judgement = judgement;
      }//if( judgNode )

      // Inherit skip from problem level if not set on question
      if( problemSkip && !question.skip )
      {
        question.skip = true;
        if( question.skipReason.empty() )
          question.skipReason = problemSkipReason;
      }

      problem.questions.push_back( question );
    }//for each Question

    if( !problem.questions.empty() || !problem.sequenceSteps.empty() )
      problems.push_back( std::move( problem ) );
  }//for each Problem

  return problems;
}//parseXml


vector<pair<string,string>> LlmBenchmarkRunner::findBenchmarkFiles( const string &directory )
{
  vector<pair<string,string>> result;

  const vector<string> files = SpecUtils::recursive_ls( directory, "_llm_benchmark.xml" );

  for( const string &filePath : files )
  {
    const string filename = SpecUtils::filename( filePath );
    // Strip the _llm_benchmark.xml suffix for display name
    string displayName = filename;
    const size_t suffixPos = displayName.rfind( "_llm_benchmark.xml" );
    if( suffixPos != string::npos )
      displayName = displayName.substr( 0, suffixPos );

    // Replace underscores with spaces for readability
    std::replace( displayName.begin(), displayName.end(), '_', ' ' );

    result.push_back( make_pair( displayName, filePath ) );
  }

  // Sort by display name
  std::sort( result.begin(), result.end(),
    []( const pair<string,string> &a, const pair<string,string> &b ) {
      return a.first < b.first;
    }
  );

  return result;
}//findBenchmarkFiles


// ---- Benchmark Execution ----

void LlmBenchmarkRunner::startBenchmark( const string &xmlFilePath )
{
  if( isRunning() )
  {
    logError( "Cannot start benchmark - already running" );
    return;
  }

  log( "Starting benchmark from: " + xmlFilePath );

  m_benchmarkBaseDir = SpecUtils::parent_path( xmlFilePath );

  try
  {
    m_problems = parseXml( xmlFilePath );
  }catch( const exception &e )
  {
    logError( "Failed to parse benchmark XML: " + string( e.what() ) );
    return;
  }

  log( "Parsed " + to_string( m_problems.size() ) + " problems" );

  // Initialize results
  m_results = BenchmarkResults();
  m_results.benchmarkFilePath = xmlFilePath;
  m_results.startTime = chrono::system_clock::now();

  // Get benchmark name from XML (look for Name element)
  // We already parsed it, but the name isn't in BenchmarkProblem - parse again quickly
  {
    vector<char> data;
    SpecUtils::load_file_data( xmlFilePath.c_str(), data );
    data.push_back( '\0' );
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_trim_whitespace>( data.data() );
    const rapidxml::xml_node<char> *root = XML_FIRST_NODE( &doc, "LlmBenchmark" );
    if( root )
    {
      const rapidxml::xml_node<char> *nameNode = XML_FIRST_NODE( root, "Name" );
      if( nameNode )
        m_results.benchmarkName = SpecUtils::xml_value_str( nameNode );
    }
  }

  if( m_results.benchmarkName.empty() )
    m_results.benchmarkName = SpecUtils::filename( xmlFilePath );

  // Get model name from config
  LlmInterface *llm = m_toolGui ? m_toolGui->llmInterface() : nullptr;
  if( llm && llm->config() )
    m_results.modelName = llm->config()->llmApi.model();

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  const string session = wApp ? wApp->sessionId() : string("null");
  string quiz_progress_filename = "quiz_progress_" + session + "_" + SpecUtils::filename(m_benchmarkBaseDir) + ".log";

  m_log_file = std::make_unique<std::ofstream>( quiz_progress_filename.c_str(), ios::binary | ios::out );
#endif

  log( "Benchmark: " + m_results.benchmarkName + ", Model: " + m_results.modelName );

  // Reset state
  m_currentProblem = 0;
  m_currentQuestion = 0;
  m_currentSequenceStep = 0;
  m_state = State::Idle;

  // Clear reference photopeak lines so they dont interfere with benchmark spectra
  ReferencePhotopeakDisplay *refLines = m_viewer->referenceLinesWidget();
  if( refLines )
    refLines->clearAllLines();

  // Connect to the main LLM interface signals
  if( llm )
  {
    m_conversationFinishedConn = llm->conversationFinished().connect(
      this, &LlmBenchmarkRunner::handleConversationFinished );
    m_responseErrorConn = llm->responseError().connect(
      this, &LlmBenchmarkRunner::handleResponseError );
  }

  // Create judge interface
  createJudgeInterface();

  // Start the first step
  advanceToNextStep();
}


void LlmBenchmarkRunner::createJudgeInterface()
{
  LlmInterface *llm = m_toolGui ? m_toolGui->llmInterface() : nullptr;
  if( !llm || !llm->config() )
    return;

  // Use main config (or JudgeApi config if available in future)
  shared_ptr<const LlmConfig> config = llm->config();

  m_judgeInterface = make_shared<LlmInterface>( m_viewer, config );
  m_judgeInterface->setBlockToolCalls( true );

  m_judgeConversationFinishedConn = m_judgeInterface->conversationFinished().connect(
    this, &LlmBenchmarkRunner::handleJudgeConversationFinished );
  m_judgeResponseErrorConn = m_judgeInterface->responseError().connect(
    this, &LlmBenchmarkRunner::handleJudgeResponseError );

  log( "Created judge LLM interface (tool calls blocked)" );
}


void LlmBenchmarkRunner::advanceToNextStep()
{
  // Check if we're done with all problems
  if( m_currentProblem >= m_problems.size() )
  {
    log( "All problems complete" );
    m_state = State::Finished;
    m_results.endTime = chrono::system_clock::now();
    computeSummaryStats();
    generateReport();
    showResultsDialog();

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  m_log_file.reset();
#endif

    m_benchmarkFinished.emit();
    return;
  }

  const BenchmarkProblem &problem = m_problems[m_currentProblem];

  log( "Problem " + to_string( m_currentProblem + 1 ) + "/" + to_string( m_problems.size() )
       + " [" + problem.id + "] category=" + problem.category
       + " difficulty=" + problem.difficulty );

  // Check if this problem has a SpectrumSequence with remaining steps
  if( !problem.sequenceSteps.empty() && (m_currentSequenceStep < problem.sequenceSteps.size()) )
  {
    loadSequenceStep();
    return;
  }

  // If we just finished a sequence, proceed to questions
  if( !problem.sequenceSteps.empty() && (m_currentSequenceStep >= problem.sequenceSteps.size()) )
  {
    // Reconnect spectrum changed handler after sequence
    reconnectSpectrumChangedHandler();

    if( m_currentQuestion < problem.questions.size() )
    {
      sendCurrentQuestion();
      return;
    }
    // else fall through to advance to next problem
  }

  // Check if we need to load spectrum for this problem (first time entering this problem)
  if( (m_currentQuestion == 0) && (m_currentSequenceStep == 0) )
  {
    if( !problem.spectrumFiles.empty() )
    {
      loadCurrentSpectrum();
      return;
    }

    // No spectrum files — this is a knowledge-only question.
    //  Clear all displayed spectra so the LLM isnt distracted by irrelevant data
    //  from a previous problem, and start a fresh conversation.
    log( "No spectrum files for this problem - clearing spectra for knowledge-only question" );

    const set<int> empty_samples;
    m_viewer->setSpectrum( nullptr, empty_samples, SpecUtils::SpectrumType::Foreground );
    m_viewer->setSpectrum( nullptr, empty_samples, SpecUtils::SpectrumType::SecondForeground );
    m_viewer->setSpectrum( nullptr, empty_samples, SpecUtils::SpectrumType::Background );

    // Explicitly clear conversation history; we cant rely on the spectrum-changed
    //  signal flow because InterSpec::setSpectrum skips the history swap when
    //  meas is nullptr.
    m_toolGui->clearConversationHistory();
  }//if( first time entering problem, no spectra )

  // Check if there are more questions
  if( m_currentQuestion < problem.questions.size() )
  {
    sendCurrentQuestion();
    return;
  }

  // Done with this problem - advance to next
  m_currentProblem++;
  m_currentQuestion = 0;
  m_currentSequenceStep = 0;
  advanceToNextStep();
}


void LlmBenchmarkRunner::loadCurrentSpectrum()
{
  const BenchmarkProblem *problem = currentProblem();
  if( !problem )
    return;

  m_state = State::LoadingSpectrum;

  // Save the conversation history from the previous problem for later reference,
  //  then let the normal spectrum-change flow clear it so the new problem starts fresh.
  {
    auto prevHistory = m_toolGui->getConversationHistory();
    if( prevHistory && !prevHistory->empty() )
    {
      log( "Saved " + to_string( prevHistory->size() )
           + " conversation(s) from previous problem" );
      // The conversations are stored in m_results.questionResults[].llmAnswer already,
      //  so we don't need to store the full history separately.
    }

    // Reconnect the LlmToolGui spectrum-changed handler so the normal setSpectrum()
    //  flow will clear the conversation history when the foreground changes.
    m_toolGui->reconnectSpectrumChangedForBenchmark();
  }

  // Load each spectrum file.
  //  Note: loadFromFileSystem() calls setSpectrum() synchronously, which fires
  //  displayedSpectrumChanged and then swaps the LLM conversation history.
  //  For the foreground load, this clears the conversation (new file has no history).
  //  We do NOT send the question from a signal handler — we send it after all files load.
  SpecMeasManager *fileMgr = m_viewer->fileManager();
  for( const BenchmarkSpectrumFile &sf : problem->spectrumFiles )
  {
    const string fullPath = resolvePath( sf.filePath );
    log( "Loading spectrum: " + fullPath + " as " +
         (sf.specType == SpecUtils::SpectrumType::Foreground ? "Foreground" :
          sf.specType == SpecUtils::SpectrumType::Background ? "Background" : "Secondary") );

    if( !fileMgr->loadFromFileSystem( fullPath, sf.specType, SpecUtils::ParserType::Auto, false ) )
    {
      logError( "Failed to load spectrum file: " + fullPath );

      // Record error for all questions in this problem
      for( const BenchmarkQuestion &q : problem->questions )
      {
        BenchmarkQuestionResult result;
        result.problemId = problem->id;
        result.questionPart = q.part;
        result.prompt = q.prompt;
        result.hadError = true;
        result.errorMessage = "Failed to load spectrum file: " + sf.filePath;
        if( q.judgement.has_value() )
        {
          result.graded = true;
          result.expectedAnswer = q.judgement->expectedAnswer;
        }
        m_results.questionResults.push_back( result );
      }

      // Skip to next problem
      m_currentProblem++;
      m_currentQuestion = 0;
      m_currentSequenceStep = 0;
      m_state = State::Idle;
      advanceToNextStep();
      return;
    }
  }//for( each spectrum file )

  // Now disconnect the spectrum-changed handler so that subsequent tool-triggered
  //  spectrum changes (e.g., from sub-agents) don't cancel in-flight requests.
  m_toolGui->disconnectSpectrumChangedForBenchmark();

  // Verify the conversation history is empty (the foreground load should have cleared it)
  {
    auto history = m_toolGui->getConversationHistory();
    const size_t numConvos = (history ? history->size() : 0);
    if( numConvos > 0 )
      log( "Warning: conversation history has " + to_string( numConvos )
           + " entries after loading new problem spectra (expected 0)" );
    else
      log( "Conversation history is empty - fresh start for new problem" );
  }

  // All spectrum files loaded — proceed to ask the first question
  m_state = State::Idle;
  sendCurrentQuestion();
}


void LlmBenchmarkRunner::handleSpectrumLoaded( SpecUtils::SpectrumType /*specType*/,
                                                shared_ptr<SpecMeas> /*meas*/,
                                                set<int> /*samples*/,
                                                vector<string> /*detectors*/ )
{
  if( m_state != State::LoadingSpectrum && m_state != State::LoadingSequenceStep )
    return;

  // Disconnect this one-shot handler
  m_spectrumChangedConn.disconnect();

  log( "Spectrum loaded successfully" );

  if( m_state == State::LoadingSequenceStep )
  {
    sendSequenceStepPrompt();
    return;
  }

  // For normal loading, give a brief moment for the UI to settle, then send question
  // In Wt, we can just proceed directly since everything is on the same event loop
  m_state = State::Idle;
  sendCurrentQuestion();
}


void LlmBenchmarkRunner::loadSequenceStep()
{
  const BenchmarkProblem *problem = currentProblem();
  if( !problem || m_currentSequenceStep >= problem->sequenceSteps.size() )
    return;

  const BenchmarkSequenceStep &step = problem->sequenceSteps[m_currentSequenceStep];

  log( "Sequence step " + to_string( m_currentSequenceStep + 1 )
       + "/" + to_string( problem->sequenceSteps.size() ) );

  // Disconnect spectrum changed handler so conversation is preserved across loads
  disconnectSpectrumChangedHandler();

  m_state = State::LoadingSequenceStep;

  // Connect our own handler for this load
  m_spectrumChangedConn = m_viewer->displayedSpectrumChanged().connect(
    this, &LlmBenchmarkRunner::handleSpectrumLoaded );

  SpecMeasManager *fileMgr = m_viewer->fileManager();
  for( const BenchmarkSpectrumFile &sf : step.spectra )
  {
    const string fullPath = resolvePath( sf.filePath );
    log( "  Loading: " + fullPath );

    if( !fileMgr->loadFromFileSystem( fullPath, sf.specType, SpecUtils::ParserType::Auto, false ) )
    {
      logError( "Failed to load sequence step spectrum: " + fullPath );
      m_spectrumChangedConn.disconnect();
      reconnectSpectrumChangedHandler();

      // Skip entire problem on sequence step failure
      m_currentProblem++;
      m_currentQuestion = 0;
      m_currentSequenceStep = 0;
      m_state = State::Idle;
      advanceToNextStep();
      return;
    }
  }
}


void LlmBenchmarkRunner::sendSequenceStepPrompt()
{
  const BenchmarkProblem *problem = currentProblem();
  if( !problem || m_currentSequenceStep >= problem->sequenceSteps.size() )
    return;

  const BenchmarkSequenceStep &step = problem->sequenceSteps[m_currentSequenceStep];

  m_state = State::WaitingForSequenceStepResponse;

  // Build the prompt, prepending a notice that the spectrum has changed.
  //  The LLM needs to know that previous peaks, analysis results, and tool outputs
  //  no longer apply to the newly loaded spectrum.
  string prompt;

  // Describe what spectra were just loaded
  string specDesc;
  for( const BenchmarkSpectrumFile &sf : step.spectra )
  {
    if( !specDesc.empty() )
      specDesc += "; ";
    const string typeStr = (sf.specType == SpecUtils::SpectrumType::Foreground) ? "Foreground"
      : (sf.specType == SpecUtils::SpectrumType::Background) ? "Background" : "Secondary";
    specDesc += typeStr + ": " + SpecUtils::filename( sf.filePath );
  }

  if( m_currentSequenceStep == 0 )
  {
    prompt = "[New spectrum loaded: " + specDesc + "]\n\n";
  }
  else
  {
    prompt = "[New spectrum loaded: " + specDesc + ". "
      "The previous spectrum has been replaced. "
      "All previously fit peaks and analysis results no longer apply to this spectrum "
      "- you must start your analysis fresh.]\n\n";
  }

  prompt += step.prompt;

  log( "Sending sequence step prompt: "
       + (prompt.size() > 100 ? prompt.substr( 0, 100 ) + "..." : prompt) );

  m_toolGui->submitMessageAsUser( prompt );
}


void LlmBenchmarkRunner::sendCurrentQuestion()
{
  const BenchmarkQuestion *question = currentQuestion();
  const BenchmarkProblem *problem = currentProblem();
  if( !question || !problem )
    return;

  // Handle skipped questions
  if( question->skip )
  {
    const string reason = question->skipReason.empty()
      ? "Question marked as skip" : question->skipReason;
    log( "Skipping question part " + to_string( question->part ) + ": " + reason );

    BenchmarkQuestionResult result;
    result.problemId = problem->id;
    result.questionPart = question->part;
    result.prompt = question->prompt;
    result.llmAnswer = "[SKIPPED] " + reason;
    result.graded = false;
    result.duration = chrono::milliseconds( 0 );
    if( question->judgement.has_value() )
    {
      result.expectedAnswer = question->judgement->expectedAnswer;
    }
    m_results.questionResults.push_back( result );

    m_currentQuestion++;
    m_state = State::Idle;
    advanceToNextStep();
    return;
  }

  m_state = State::WaitingForResponse;
  m_questionStartTime = chrono::system_clock::now();

  string prompt = question->prompt;

  // Prepend context to first question of a problem
  if( (m_currentQuestion == 0) && !problem->context.empty() && problem->sequenceSteps.empty() )
  {
    prompt = problem->context + "\n\n" + prompt;
  }

  log( "Sending question part " + to_string( question->part ) + ": "
       + (prompt.size() > 100 ? prompt.substr( 0, 100 ) + "..." : prompt) );

  // TODO: handle question->imagePath for image questions (sendUserMessageWithImage)
  m_toolGui->submitMessageAsUser( prompt );
}


void LlmBenchmarkRunner::handleConversationFinished()
{
  if( m_state == State::WaitingForSequenceStepResponse )
  {
    log( "Sequence step response received" );

    // Advance to next sequence step
    m_currentSequenceStep++;
    m_state = State::Idle;
    advanceToNextStep();
    return;
  }

  if( m_state == State::WaitingForResponse )
  {
    log( "Question response received" );
    extractAnswerAndJudge();
    return;
  }

  // Ignore conversationFinished signals when not waiting for a response
  // (could be from judging or other activity)
}


void LlmBenchmarkRunner::handleResponseError()
{
  if( m_state != State::WaitingForResponse && m_state != State::WaitingForSequenceStepResponse )
    return;

  const BenchmarkQuestion *question = currentQuestion();
  const BenchmarkProblem *problem = currentProblem();

  logError( "LLM response error for problem " + (problem ? problem->id : "?")
            + " question " + to_string( question ? question->part : -1 ) );

  if( m_state == State::WaitingForSequenceStepResponse )
  {
    // Skip entire problem on sequence step error
    reconnectSpectrumChangedHandler();
    m_currentProblem++;
    m_currentQuestion = 0;
    m_currentSequenceStep = 0;
    m_state = State::Idle;
    advanceToNextStep();
    return;
  }

  // Record error result for the current question
  if( question && problem )
  {
    BenchmarkQuestionResult result;
    result.problemId = problem->id;
    result.questionPart = question->part;
    result.prompt = question->prompt;
    result.hadError = true;
    result.errorMessage = "LLM response error";
    const auto now = chrono::system_clock::now();
    result.duration = chrono::duration_cast<chrono::milliseconds>( now - m_questionStartTime );
    if( question->judgement.has_value() )
    {
      result.graded = true;
      result.expectedAnswer = question->judgement->expectedAnswer;
    }
    m_results.questionResults.push_back( result );
  }

  // Advance to next question
  m_currentQuestion++;
  m_state = State::Idle;
  advanceToNextStep();
}


void LlmBenchmarkRunner::extractAnswerAndJudge()
{
  const BenchmarkQuestion *question = currentQuestion();
  const BenchmarkProblem *problem = currentProblem();
  if( !question || !problem )
    return;

  // Extract the LLM's answer from the conversation history
  string llmAnswer;
  LlmInterface *llm = m_toolGui ? m_toolGui->llmInterface() : nullptr;
  if( llm )
  {
    shared_ptr<LlmConversationHistory> history = llm->getHistory();
    if( history )
    {
      const vector<shared_ptr<LlmInteraction>> &conversations = history->getConversations();
      if( !conversations.empty() )
      {
        const shared_ptr<LlmInteraction> &lastConvo = conversations.back();
        // Walk responses backward to find the last FinalLlmResponse
        for( int i = static_cast<int>( lastConvo->responses.size() ) - 1; i >= 0; --i )
        {
          const shared_ptr<LlmInteractionTurn> &turn = lastConvo->responses[i];
          if( turn && turn->type() == LlmInteractionTurn::Type::FinalLlmResponse )
          {
            const LlmInteractionFinalResponse *finalResp =
              dynamic_cast<const LlmInteractionFinalResponse *>( turn.get() );
            if( finalResp )
            {
              llmAnswer = finalResp->content();
              break;
            }
          }
        }
      }
    }
  }

  const auto now = chrono::system_clock::now();
  const chrono::milliseconds duration = chrono::duration_cast<chrono::milliseconds>( now - m_questionStartTime );

  log( "Extracted answer (" + to_string( llmAnswer.size() ) + " chars, "
       + to_string( duration.count() ) + " ms): "
       + (llmAnswer.size() > 200 ? llmAnswer.substr( 0, 200 ) + "..." : llmAnswer) );

  if( !question->judgement.has_value() )
  {
    // Ungraded question - just record the answer
    BenchmarkQuestionResult result;
    result.problemId = problem->id;
    result.questionPart = question->part;
    result.prompt = question->prompt;
    result.llmAnswer = llmAnswer;
    result.graded = false;
    result.duration = duration;
    recordResult( std::move( result ) );
    return;
  }

  const BenchmarkJudgement &judgement = question->judgement.value();

  if( judgement.type == BenchmarkJudgement::Type::ContainsString )
  {
    judgeWithStringMatch( llmAnswer, judgement );
  }
  else
  {
    judgeWithLlm( llmAnswer, *question );
  }
}


void LlmBenchmarkRunner::judgeWithStringMatch( const string &llmAnswer,
                                                const BenchmarkJudgement &judgement )
{
  const BenchmarkQuestion *question = currentQuestion();
  const BenchmarkProblem *problem = currentProblem();
  if( !question || !problem )
    return;

  bool found = false;
  if( judgement.caseSensitive )
  {
    found = (llmAnswer.find( judgement.expectedAnswer ) != string::npos);
  }
  else
  {
    found = SpecUtils::icontains( llmAnswer, judgement.expectedAnswer );
  }

  const auto now = chrono::system_clock::now();

  BenchmarkQuestionResult result;
  result.problemId = problem->id;
  result.questionPart = question->part;
  result.prompt = question->prompt;
  result.llmAnswer = llmAnswer;
  result.expectedAnswer = judgement.expectedAnswer;
  result.graded = true;
  result.correct = found;
  result.score = found ? 1 : 0;
  result.maxScore = 1;
  result.judgementReason = found
    ? "Answer contains expected string: \"" + judgement.expectedAnswer + "\""
    : "Answer does not contain expected string: \"" + judgement.expectedAnswer + "\"";
  result.duration = chrono::duration_cast<chrono::milliseconds>( now - m_questionStartTime );

  log( "ContainsString judge: " + string( found ? "PASS" : "FAIL" )
       + " (looking for \"" + judgement.expectedAnswer + "\")" );

  recordResult( std::move( result ) );
}


void LlmBenchmarkRunner::judgeWithLlm( const string &llmAnswer,
                                        const BenchmarkQuestion &question )
{
  if( !m_judgeInterface )
  {
    logError( "No judge interface available - falling back to ungraded" );

    const BenchmarkProblem *problem = currentProblem();
    BenchmarkQuestionResult result;
    result.problemId = problem ? problem->id : "";
    result.questionPart = question.part;
    result.prompt = question.prompt;
    result.llmAnswer = llmAnswer;
    result.graded = false;
    result.judgementReason = "No judge interface available";
    const auto now = chrono::system_clock::now();
    result.duration = chrono::duration_cast<chrono::milliseconds>( now - m_questionStartTime );
    recordResult( std::move( result ) );
    return;
  }

  m_state = State::Judging;

  const BenchmarkJudgement &judgement = question.judgement.value();

  // Build the judge prompt
  string judgePrompt = "You are evaluating an LLM's answer to a spectral analysis question.\n\n";
  judgePrompt += "**Question asked:**\n" + question.prompt + "\n\n";
  judgePrompt += "**Expected answer:**\n" + judgement.expectedAnswer + "\n\n";
  judgePrompt += "**Actual answer from LLM:**\n" + llmAnswer + "\n\n";

  if( !judgement.rubric.empty() )
  {
    judgePrompt += "**Scoring rubric:**\n";
    int maxPoints = 0;
    for( const BenchmarkCriterion &c : judgement.rubric )
    {
      judgePrompt += "- " + to_string( c.points ) + " points: " + c.description + "\n";
      if( c.points > 0 )
        maxPoints += c.points;
    }
    judgePrompt += "\nMax possible score: " + to_string( maxPoints ) + "\n\n";
  }

  judgePrompt += "Respond with ONLY a JSON object (no markdown code fences, no other text):\n";
  judgePrompt += "{\"correct\": true/false, \"score\": <points_earned>, \"maxScore\": <max_possible>, \"reason\": \"<brief explanation>\"}\n";
  judgePrompt += "\nIf there is no rubric, use score 1 for correct and 0 for incorrect, with maxScore 1.\n";
  judgePrompt += "Be generous with partial credit — the answer doesn't need to be word-for-word, just scientifically correct and addressing the key points.";

  log( "Sending to judge LLM (" + to_string( judgePrompt.size() ) + " chars)" );

  // Clear judge history for a fresh conversation
  shared_ptr<LlmConversationHistory> judgeHistory = m_judgeInterface->getHistory();
  if( judgeHistory )
    judgeHistory->clear();

  m_judgeInterface->sendUserMessage( judgePrompt );
}


void LlmBenchmarkRunner::handleJudgeConversationFinished()
{
  if( m_state != State::Judging )
    return;

  const BenchmarkQuestion *question = currentQuestion();
  const BenchmarkProblem *problem = currentProblem();

  // Extract judge's response
  string judgeResponse;
  shared_ptr<LlmConversationHistory> judgeHistory = m_judgeInterface->getHistory();
  if( judgeHistory )
  {
    const vector<shared_ptr<LlmInteraction>> &conversations = judgeHistory->getConversations();
    if( !conversations.empty() )
    {
      const shared_ptr<LlmInteraction> &lastConvo = conversations.back();
      for( int i = static_cast<int>( lastConvo->responses.size() ) - 1; i >= 0; --i )
      {
        const shared_ptr<LlmInteractionTurn> &turn = lastConvo->responses[i];
        if( turn && turn->type() == LlmInteractionTurn::Type::FinalLlmResponse )
        {
          const LlmInteractionFinalResponse *finalResp =
            dynamic_cast<const LlmInteractionFinalResponse *>( turn.get() );
          if( finalResp )
          {
            judgeResponse = finalResp->content();
            break;
          }
        }
      }
    }
  }

  log( "Judge response: " + (judgeResponse.size() > 300 ? judgeResponse.substr( 0, 300 ) + "..." : judgeResponse) );

  BenchmarkQuestionResult result;
  result.problemId = problem ? problem->id : "";
  result.questionPart = question ? question->part : 0;
  result.prompt = question ? question->prompt : "";
  result.expectedAnswer = (question && question->judgement.has_value())
    ? question->judgement->expectedAnswer : "";
  result.graded = true;

  const auto now = chrono::system_clock::now();
  result.duration = chrono::duration_cast<chrono::milliseconds>( now - m_questionStartTime );

  // Extract the LLM answer from the main conversation (same as extractAnswerAndJudge did)
  LlmInterface *llm = m_toolGui ? m_toolGui->llmInterface() : nullptr;
  if( llm )
  {
    shared_ptr<LlmConversationHistory> history = llm->getHistory();
    if( history )
    {
      const vector<shared_ptr<LlmInteraction>> &conversations = history->getConversations();
      if( !conversations.empty() )
      {
        const shared_ptr<LlmInteraction> &lastConvo = conversations.back();
        for( int i = static_cast<int>( lastConvo->responses.size() ) - 1; i >= 0; --i )
        {
          const shared_ptr<LlmInteractionTurn> &turn = lastConvo->responses[i];
          if( turn && turn->type() == LlmInteractionTurn::Type::FinalLlmResponse )
          {
            const LlmInteractionFinalResponse *finalResp =
              dynamic_cast<const LlmInteractionFinalResponse *>( turn.get() );
            if( finalResp )
            {
              result.llmAnswer = finalResp->content();
              break;
            }
          }
        }
      }
    }
  }

  // Parse judge JSON response
  try
  {
    // Strip markdown code fences if present
    string jsonStr = judgeResponse;
    {
      const size_t startFence = jsonStr.find( "```" );
      if( startFence != string::npos )
      {
        size_t contentStart = jsonStr.find( '\n', startFence );
        if( contentStart != string::npos )
          contentStart++;
        else
          contentStart = startFence + 3;

        const size_t endFence = jsonStr.find( "```", contentStart );
        if( endFence != string::npos )
          jsonStr = jsonStr.substr( contentStart, endFence - contentStart );
        else
          jsonStr = jsonStr.substr( contentStart );
      }
    }

    SpecUtils::trim( jsonStr );

    const nlohmann::json j = nlohmann::json::parse( jsonStr );
    result.correct = j.value( "correct", false );
    result.score = j.value( "score", result.correct ? 1 : 0 );
    result.maxScore = j.value( "maxScore", 1 );
    result.judgementReason = j.value( "reason", "" );

    log( "LLM Judge: " + string( result.correct ? "PASS" : "FAIL" )
         + " score=" + to_string( result.score ) + "/" + to_string( result.maxScore )
         + " reason: " + result.judgementReason );
  }catch( const exception &e )
  {
    logError( "Failed to parse judge response as JSON: " + string( e.what() ) );
    result.correct = false;
    result.score = 0;
    result.maxScore = 1;
    result.judgementReason = "Failed to parse judge response: " + judgeResponse;
  }

  recordResult( std::move( result ) );
}


void LlmBenchmarkRunner::handleJudgeResponseError()
{
  if( m_state != State::Judging )
    return;

  logError( "Judge LLM response error" );

  const BenchmarkQuestion *question = currentQuestion();
  const BenchmarkProblem *problem = currentProblem();

  BenchmarkQuestionResult result;
  result.problemId = problem ? problem->id : "";
  result.questionPart = question ? question->part : 0;
  result.prompt = question ? question->prompt : "";
  result.graded = true;
  result.correct = false;
  result.hadError = true;
  result.errorMessage = "Judge LLM response error";
  const auto now = chrono::system_clock::now();
  result.duration = chrono::duration_cast<chrono::milliseconds>( now - m_questionStartTime );

  recordResult( std::move( result ) );
}


void LlmBenchmarkRunner::recordResult( BenchmarkQuestionResult result )
{
  log( "Recorded result for " + result.problemId + " Q" + to_string( result.questionPart )
       + ": " + (result.hadError ? "ERROR" : (result.graded ? (result.correct ? "PASS" : "FAIL") : "UNGRADED")) );

#if( PERFORM_DEVELOPER_CHECKS && BUILD_AS_LOCAL_SERVER )
  const BenchmarkProblem *prev_problem = currentProblem();
  assert( prev_problem );

  const string session = wApp ? wApp->sessionId() : string("null");
  string quiz_progress_filename = "quiz_progress_" + session + "_" + SpecUtils::filename(m_benchmarkBaseDir) + ".log";
  std::stringstream output;
  output
  << "-----------------------------------------------------------------------------------" << endl
  << "Problem " << m_currentQuestion << " of " << (m_problems.size() - 1) << ":" << endl
  << "\tproblemId: " << result.problemId << endl
  //std::string prompt;
  << "\tquestionPart: " << result.questionPart << endl
  << "\tllmAnswer: " << result.llmAnswer << endl
  << "\texpectedAnswer: " << result.expectedAnswer << endl
  << "\thadError: " << result.hadError << (result.hadError ? string(" errorm msg: " + result.errorMessage) : string("")) << endl
  << "\tgraded: " << result.graded << endl
  << "\tcorrect: " << result.correct << endl
  << "\tscore: " << result.score << " out of " << result.maxScore << endl
  << "\tjudgementReason: " << result.judgementReason << endl
  << "\tduration: " << (result.duration.count() / 1000) << " seconds" << endl
  << "\tcompletionTokens: " << (result.completionTokens.has_value() ? std::to_string(result.completionTokens.value()) : string("n/a")) << endl
  << endl;

  log( output.str() );
#endif


  m_results.questionResults.push_back( std::move( result ) );

  // Advance to next question
  m_currentQuestion++;

  const BenchmarkProblem *problem = currentProblem();
  if( problem && m_currentQuestion >= problem->questions.size() )
  {
    // Done with this problem
    m_currentProblem++;
    m_currentQuestion = 0;
    m_currentSequenceStep = 0;
  }

  m_state = State::Idle;
  advanceToNextStep();
}


void LlmBenchmarkRunner::disconnectSpectrumChangedHandler()
{
  if( m_spectrumChangeDisconnected )
    return;

  m_toolGui->disconnectSpectrumChangedForBenchmark();
  m_spectrumChangeDisconnected = true;
  log( "Disconnected spectrum changed handler (SpectrumSequence)" );
}


void LlmBenchmarkRunner::reconnectSpectrumChangedHandler()
{
  if( !m_spectrumChangeDisconnected )
    return;

  m_toolGui->reconnectSpectrumChangedForBenchmark();
  m_spectrumChangeDisconnected = false;
  log( "Reconnected spectrum changed handler" );
}


void LlmBenchmarkRunner::computeSummaryStats()
{
  m_results.totalQuestions = 0;
  m_results.gradedQuestions = 0;
  m_results.correctCount = 0;
  m_results.errorCount = 0;
  m_results.totalScore = 0;
  m_results.maxTotalScore = 0;

  for( const BenchmarkQuestionResult &r : m_results.questionResults )
  {
    m_results.totalQuestions++;
    if( r.hadError )
      m_results.errorCount++;
    if( r.graded )
    {
      m_results.gradedQuestions++;
      if( r.correct )
        m_results.correctCount++;
      m_results.totalScore += r.score;
      m_results.maxTotalScore += r.maxScore;
    }
  }
}


void LlmBenchmarkRunner::generateReport()
{
  computeSummaryStats();

  const auto totalDuration = chrono::duration_cast<chrono::seconds>(
    m_results.endTime - m_results.startTime );

  log( "=== LLM Benchmark Results ===" );
  log( "Benchmark: " + m_results.benchmarkName );
  log( "Model: " + m_results.modelName );
  log( "Total time: " + to_string( totalDuration.count() ) + "s" );
  log( "Questions: " + to_string( m_results.totalQuestions )
       + " total, " + to_string( m_results.gradedQuestions ) + " graded" );

  if( m_results.gradedQuestions > 0 )
  {
    log( "Correct: " + to_string( m_results.correctCount )
         + "/" + to_string( m_results.gradedQuestions )
         + " (" + to_string( m_results.gradedQuestions > 0
                             ? (100 * m_results.correctCount / m_results.gradedQuestions) : 0 ) + "%)" );
  }

  if( m_results.maxTotalScore > 0 )
  {
    log( "Score: " + to_string( m_results.totalScore )
         + "/" + to_string( m_results.maxTotalScore )
         + " (" + to_string( 100 * m_results.totalScore / m_results.maxTotalScore ) + "%)" );
  }

  if( m_results.errorCount > 0 )
    log( "Errors: " + to_string( m_results.errorCount ) );

  log( "" );
  log( "Problem Results:" );

  for( const BenchmarkQuestionResult &r : m_results.questionResults )
  {
    const bool isSkipped = (r.llmAnswer.substr( 0, 9 ) == "[SKIPPED]");
    string status;
    if( isSkipped )
      status = "SKIP";
    else if( r.hadError )
      status = "ERROR";
    else if( !r.graded )
      status = "UNGRADED";
    else if( r.correct )
      status = "PASS";
    else
      status = "FAIL";

    string line = "  [" + status + "] " + r.problemId + " Q" + to_string( r.questionPart );

    if( r.graded && !r.hadError )
      line += ": " + to_string( r.score ) + "/" + to_string( r.maxScore ) + " pts";

    if( !isSkipped )
      line += " (" + to_string( r.duration.count() ) + "ms)";

    if( isSkipped )
      line += " - " + r.llmAnswer;
    else if( r.hadError && !r.errorMessage.empty() )
      line += " - " + r.errorMessage;
    else if( r.graded && !r.correct && !r.judgementReason.empty() )
    {
      const string reason = r.judgementReason.size() > 80
        ? r.judgementReason.substr( 0, 80 ) + "..." : r.judgementReason;
      line += " - " + reason;
    }

    log( line );
  }

  log( "=== End Benchmark Results ===" );
}


void LlmBenchmarkRunner::showResultsDialog()
{
  // Build results JSON
  nlohmann::json resultsJson;
  resultsJson["benchmarkName"] = m_results.benchmarkName;
  resultsJson["modelName"] = m_results.modelName;
  resultsJson["benchmarkFile"] = m_results.benchmarkFilePath;
  resultsJson["totalQuestions"] = m_results.totalQuestions;
  resultsJson["gradedQuestions"] = m_results.gradedQuestions;
  resultsJson["correctCount"] = m_results.correctCount;
  resultsJson["errorCount"] = m_results.errorCount;
  resultsJson["totalScore"] = m_results.totalScore;
  resultsJson["maxTotalScore"] = m_results.maxTotalScore;

  nlohmann::json questionsArray = nlohmann::json::array();
  for( const BenchmarkQuestionResult &r : m_results.questionResults )
  {
    nlohmann::json qj;
    qj["problemId"] = r.problemId;
    qj["questionPart"] = r.questionPart;
    qj["prompt"] = r.prompt;
    qj["llmAnswer"] = r.llmAnswer;
    qj["expectedAnswer"] = r.expectedAnswer;
    qj["graded"] = r.graded;
    qj["correct"] = r.correct;
    qj["score"] = r.score;
    qj["maxScore"] = r.maxScore;
    qj["judgementReason"] = r.judgementReason;
    qj["durationMs"] = r.duration.count();
    qj["hadError"] = r.hadError;
    qj["errorMessage"] = r.errorMessage;
    questionsArray.push_back( qj );
  }
  resultsJson["questions"] = questionsArray;

  const string jsonStr = resultsJson.dump( 2 );

  // Show a summary dialog
  SimpleDialog *dialog = new SimpleDialog( "LLM Benchmark Results" );

  string summaryHtml = "<div style='font-family: monospace; font-size: 12px; overflow-y: auto; max-height: 75vh'>";
  summaryHtml += "<b>" + m_results.benchmarkName + "</b><br/>";
  summaryHtml += "Model: " + m_results.modelName + "<br/><br/>";

  if( m_results.gradedQuestions > 0 )
  {
    const int pct = (100 * m_results.correctCount / m_results.gradedQuestions);
    summaryHtml += "Correct: " + to_string( m_results.correctCount )
                 + "/" + to_string( m_results.gradedQuestions )
                 + " (" + to_string( pct ) + "%)<br/>";
  }

  if( m_results.maxTotalScore > 0 )
  {
    const int pct = (100 * m_results.totalScore / m_results.maxTotalScore);
    summaryHtml += "Score: " + to_string( m_results.totalScore )
                 + "/" + to_string( m_results.maxTotalScore )
                 + " (" + to_string( pct ) + "%)<br/>";
  }

  if( m_results.errorCount > 0 )
    summaryHtml += "Errors: " + to_string( m_results.errorCount ) + "<br/>";

  summaryHtml += "<br/><table style='border-collapse: collapse; width: 100%;'>";
  summaryHtml += "<tr><th style='text-align:left; padding:2px 6px;'>Status</th>"
                 "<th style='text-align:left; padding:2px 6px;'>Problem</th>"
                 "<th style='text-align:left; padding:2px 6px;'>Q</th>"
                 "<th style='text-align:left; padding:2px 6px;'>Score</th>"
                 "<th style='text-align:left; padding:2px 6px;'>Time</th></tr>";

  for( const BenchmarkQuestionResult &r : m_results.questionResults )
  {
    const bool isSkip = (r.llmAnswer.substr( 0, 9 ) == "[SKIPPED]");
    string statusColor;
    string status;
    if( isSkip )
    {
      status = "SKIP";
      statusColor = "#aaa";
    }
    else if( r.hadError )
    {
      status = "ERROR";
      statusColor = "#ff6b6b";
    }
    else if( !r.graded )
    {
      status = "---";
      statusColor = "#999";
    }
    else if( r.correct )
    {
      status = "PASS";
      statusColor = "#51cf66";
    }
    else
    {
      status = "FAIL";
      statusColor = "#ff6b6b";
    }

    summaryHtml += "<tr>"
      "<td style='padding:2px 6px; color:" + statusColor + ";'>" + status + "</td>"
      "<td style='padding:2px 6px;'>" + r.problemId + "</td>"
      "<td style='padding:2px 6px;'>" + to_string( r.questionPart ) + "</td>"
      "<td style='padding:2px 6px;'>" + (r.graded ? (to_string( r.score ) + "/" + to_string( r.maxScore )) : "-") + "</td>"
      "<td style='padding:2px 6px;'>" + to_string( r.duration.count() / 1000 ) + "s</td>"
      "</tr>";
  }

  summaryHtml += "</table></div>";

  new Wt::WText( summaryHtml, Wt::XHTMLUnsafeText, dialog->contents() );

  // Add a hidden textarea with the JSON for copy/save
  Wt::WContainerWidget *jsonContainer = new Wt::WContainerWidget( dialog->contents() );
  jsonContainer->setMargin( 8, Wt::Top );

  Wt::WPushButton *copyBtn = new Wt::WPushButton( "Copy JSON", dialog->footer() );
  copyBtn->setStyleClass( "simple-dialog-btn" );

  // Use JavaScript to copy
  const string escapedJson = Wt::WWebWidget::jsStringLiteral( jsonStr );
  copyBtn->clicked().connect(
    "function(){ navigator.clipboard.writeText(" + escapedJson + ").then(function(){"
    "alert('Results JSON copied to clipboard');}); }" );

  Wt::WPushButton *okBtn = dialog->addButton( "OK" );
  okBtn->clicked().connect( dialog, &SimpleDialog::accept );
}


#endif // USE_LLM_INTERFACE
