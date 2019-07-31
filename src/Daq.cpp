#include "InterSpec_config.h"

#include <memory>
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>
#include <utility>
#include <stdexcept>

#include <Wt/WServer>
#include <Wt/WApplication>

#include "InterSpec/Daq.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "SpecUtils/SpectrumDataStructs.h"

using namespace std;

namespace Daq
{
  /** For the example code, we'll use the passthrough example included with InterSpec
   and "add" data every few seconds, like it was comming from the DAQ.
   */
  set<int>::iterator sm_currentSampleNumIter;
  MeasurementInfo sm_exampleData;
  Measurement sm_currentSum;
  std::atomic<bool> sm_keepGoing;
  std::thread sm_exampleRunner;
  
void init_daq( int argc, char **argv )
{
  const bool loaded = sm_exampleData.load_file( "example_spectra/passthrough.n42", ParserType::kAutoParser );
  
  if( !loaded || sm_exampleData.sample_numbers().empty() )
    throw std::runtime_error( "Failed to open example file" );
  
  sm_currentSampleNumIter = std::begin( sm_exampleData.sample_numbers() );
  sm_keepGoing = true;
  
  auto worker = [](){
    while( sm_keepGoing )
    {
      std::this_thread::sleep_for( std::chrono::seconds(3) );
      if( !sm_keepGoing )
        break;
      
      ++sm_currentSampleNumIter;
      const set<int> samples( std::begin( sm_exampleData.sample_numbers() ), sm_currentSampleNumIter );
      
      if( sm_currentSampleNumIter == std::end(sm_exampleData.sample_numbers()) )
        sm_currentSampleNumIter = std::begin( sm_exampleData.sample_numbers() );
      
      vector<bool> dets_to_use( sm_exampleData.detector_names().size(), true );
      
      std::shared_ptr<Measurement> sum = sm_exampleData.sum_measurements( samples, dets_to_use );
      if( !sum )
        continue;
      
      //Just always use the first four samples as backgorund.  You can also just pass a nullptr for background.
      set<int> back_samples{ 1, 2, 3, 4 };
      std::shared_ptr<Measurement> background = sm_exampleData.sum_measurements( back_samples, dets_to_use );
      
      update_displayed_spectrum( sum, background );
    }//while( sm_keepGoing )
  };//worker
  
  sm_exampleRunner = std::thread( worker );
}//void init_daq( int argc, char **argv )

  
void update_displayed_spectrum( std::shared_ptr<const Measurement> foreground,
                                 std::shared_ptr<const Measurement> background )
{
  Wt::WServer *server = Wt::WServer::instance();
  if( !server )
  {
    cerr << "No Server!" << endl;
    return;
  }
  
  server->postAll( std::bind( [foreground,background](){
    Wt::WApplication *wtap = wApp;
    if( !wtap )
    {
      //Not sure why this happens some times.
      cerr << "No WApplication::instance() in postAll(...)" << endl;
      return;
    }
    
    InterSpecApp *app = dynamic_cast<InterSpecApp *>( wtap );
    assert( app );
    
    InterSpec *interspec = app->viewer();
    assert( interspec );
  
    interspec->update_displayed_spectrum_from_daq( foreground, background );
    
    app->triggerUpdate();
  }) );
}//update_displayed_spectrum()

  
void stop_daq()
{
  sm_keepGoing = false;
  sm_exampleRunner.join();  //May take up to 3 seconds.  Use a condition variable instead of bool to avoid this.
}//void stop_daq()
  
}//namespace Daq

