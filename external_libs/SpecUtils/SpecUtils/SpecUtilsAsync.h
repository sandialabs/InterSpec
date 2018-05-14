#ifndef SpecUtilsAsync_h
#define SpecUtilsAsync_h
/* SpecUtils: a library to parse, save, and manipulate gamma spectrum data files.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include "SpecUtils_config.h"

#include <memory>
#include <mutex>
#include <vector>
#include <thread>
#include <exception>
#include <functional>
#include <condition_variable>


#if( __APPLE__ )
  //On apple systems we can use Grand Central Dispatch to power the queue,
  //  however right now to test it out for non-apple systems, I am leaving this
  //  disabled
  #include <dispatch/dispatch.h>
  #define ThreadPool_USING_GCD 1
  //#define ThreadPool_USING_WT 1
#elif( SpecUtils_USE_WT_THREADPOOL )
  #define ThreadPool_USING_WT 1
#elif( SpecUtils_USING_NO_THREADING )
  #define ThreadPool_USING_SERIAL 1
#else
  #define ThreadPool_USING_THREADS 1
#endif


namespace SpecUtilsAsync
{
  //num_logical_cpu_cores(): if cpu has hyperthreading this is 2x physical
  //  If SpecUtils_USING_NO_THREADING is true, then returns 1
  int num_logical_cpu_cores();
  int num_physical_cpu_cores();

  
  //ThreadPool: a simple adapter class to Wt::WIOService that allows posting
  //  jobs to Wt's thread pool, then being able to call join() to ensure all
  //  jobs are completed before moving on.
  //This class also allows adapting to a non-threaded environment, as well
  //  as using Apples Grand Central Dispatch (instead of Wts thread pool).
  //The actual ThreadPool class instance is not thread safe, so dont add jobs
  //  from multiple threads, or call join from a differnt thread than you
  //  created the ThreadPool object from..
  //ToDo:
  //  For Windows could use the windows proccess threadpool - see
  //    https://devel.nuclex.org/framework/browser/Nuclex.Support.Native/trunk/Source/Threading/WindowsThreadPool.cpp?rev=1671
  //    https://devel.nuclex.org/framework/browser/Nuclex.Support.Native/trunk/Include/Nuclex/Support/Threading/WindowsThreadPool.h?rev=1671
  //    https://devel.nuclex.org/framework/browser/Nuclex.Support.Native/trunk/Include/Nuclex/Support/Threading/ThreadPool.h?rev=1671
  //    Which actually in looking around, to http://msdn.microsoft.com/en-us/library/windows/desktop/ms682478(v=vs.85).aspx
  //    looks like it will be fairly easy to interface to Windows process wide threadpool
  //    Which we could then rename ThreadPool_USING_GCD to USE_OS_THREAD_POOL
  //    (Linux will still need to be dealt with, perhaps just leave it at using
  //    do_asyncronous_work).

  
  class ThreadPool
  {
    /*
     ToDo:
       -add a constructor that takes arguments:
        std::function< void(const std::function< void(void) > &) > post_fctn
        and std::function< void(void) > join_fctn, which implement the
        internal functionality.  Then convert the existing Apple GCD, Wt 
        threadpool, and do_asyncronous_work(...) functions to use this mechanism
       -Convert do_asyncronous_work to use std::async, at least on MSVC which
        uses a threadpool (other OSs could use async as well, without much 
        penalty).
     */
  public:
    ThreadPool();
    ~ThreadPool(); //calls join()
    
    //post: post a job to the thread pool.
    //  I'm not to happy about having the fcn limited to a
    //  std::function< void(void) > (although anything with a operator()(void)
    //  operator will implicitly convert to one), but I kept running into issues
    //  trying to template it properly. I'm sure its solvable, but so many more
    //  interesting things to do than avoid this hopefully small amount of
    //  overhead for creation/copy
    //Exceptions, deriving from std::exception, thrown by 'fcn' will be caught,
    //  and re-thrown when the join() function is called.
    void post( std::function< void(void) > fcn );
    
    //join(): this function will not return untill all jobs that have been
    //  posted have completed.  You can post more jobs after calling join().
    //If an exception is thrown by any of the worker jobs, it will be re-thrown
    //  when calling this function; if multiple exceptions are thrown, then the
    //  last one last one (temporaly, as evaluated) will be re-thrown, and the
    //  others lost.
    void join();
    
  protected:
    
#if( defined(ThreadPool_USING_GCD) )
    dispatch_queue_t m_queue;
    std::mutex m_exception_mutex;
    std::exception_ptr m_exception;
#elif( defined(ThreadPool_USING_WT) )
    
    //Could switch m_nJobsLeft to a std::atomic<size_t> so this would be a
    //  lockless queue, but whatevs for now
    size_t m_nJobsLeft;
    std::mutex m_mutex;
    std::condition_variable m_condition;

    //To prevent the possibility of deadlocks when submitting to the Wt
    //  threadpool, we will only allow at most
    //  Wt::WServer::instance()->ioService().threadCount() - 3
    //  ThreadPool objects to use the Wt threadpool, any subsequent ones will
    //  use do_asyncronous_work(...) to perform the work.
    static std::mutex sm_npool_mutex;
    static int sm_npools;
    bool m_canSubmitToWtPool;
    
    std::mutex m_exception_mutex;
    std::exception_ptr m_exception;
    
    //m_nonPostedWorkers: necassarry when Wt::WServer::instance() isnt available
    //  so do_asyncronous_work(...) can be called in join(), or if
    //  SpecUtils_USE_WT_THREADPOOL isnt defined, then posted jobs are put here and ran
    //  when join() is called.
    std::vector< std::function<void(void)> > m_nonPostedWorkers;
#elif( defined(ThreadPool_USING_THREADS) )
    std::mutex m_exception_mutex;
    std::exception_ptr m_exception;
    std::vector< std::function<void(void)> > m_nonPostedWorkers;
#elif( defined(ThreadPool_USING_SERIAL) )
    std::mutex m_exception_mutex;
    std::exception_ptr m_exception;
#else
#error "Exactly on of the following must be defined: ThreadPool_USING_GCD, ThreadPool_USING_WT, ThreadPool_USING_THREADS, ThreadPool_USING_SERIAL"
#endif
    
    
  private:
#if( defined(ThreadPool_USING_WT) || defined(ThreadPool_USING_THREADS) )
    //dowork(): actually calls the user passed in function, and decrements the
    //  m_nJobsLeft counter.  To be used from the backing thread pool, and not
    //  in do_asyncronous_work(...).
    //Not const due to deccrementing m_nJobsLeft, and possibly modifying
    //  m_exception, but safe to call from multiple threads simultaneously.
    void dowork( const std::function< void(void) > &fcn );
    
    //doworkasync(): calls the user passed in function when using
    //  do_asyncronous_work(...); needed to transport any exceptions to join().
    void doworkasync( const std::function< void(void) > &fcn );
#endif  //#if( defined(ThreadPool_USING_GCD) ) / elif
  };//class ThreadPool
  
  
  //do_asyncronous_work(...) takes a vector with objects that have
  //  'void operator()(void)' calls available, typically either functionoids or
  //  or 'std::function<void()>', and then performs their task asyncronously
  //  Note exceptions are not handled
  //
  //Pass in vector of functionoids - the 'Functionoid' must copy-constructable
  //  and the container passed in will be empty upon completion
  template <class Functionoid>
  void do_asyncronous_work( std::vector< Functionoid > &jobs_to_do,
                           const bool use_only_physical_cores = false );

}//namespace SpecUtilsAsync


//implementation of of templated functions
#include <algorithm>

#if( defined(ThreadPool_USING_WT) )
#include <Wt/WServer>
#include <Wt/WIOService>
#endif


#if( defined(WIN32) )
#undef min
#undef max
#endif 


namespace SpecUtilsAsync
{
#if( !SpecUtils_USING_NO_THREADING )
  template <class Functionoid>
  void async_worker( std::vector< Functionoid > &jobs, std::mutex &gaurd )
  {
    while( 1 )
    {
      std::unique_ptr<Functionoid> work;
      
      {//begin codeblock to get the worker
        std::lock_guard<std::mutex> lock( gaurd );
        if( jobs.empty() )
          return;
        work.reset( new Functionoid( jobs.back() ) );
        jobs.pop_back();
      }//end codeblock to get the worker
      
      (*work)();
    }//while(1)
  }//async_worker(...)
#endif
  
  template <class Functionoid>
  void do_asyncronous_work( std::vector< Functionoid > &job_to_do,
                           const bool physcore )
  {
    if( job_to_do.empty() )
      return;
    
#if( SpecUtils_USING_NO_THREADING )
    for( size_t i = 0; i < job_to_do.size(); ++i )
      job_to_do[i]();
#else
  std::mutex queue_mutex;
  int nthread = (physcore ? num_physical_cpu_cores() : num_logical_cpu_cores() );
  std::vector<std::shared_ptr<std::thread> > threads( nthread );
  for( int i = 0; i < nthread; ++i )
    threads[i] = std::make_shared<std::thread>( &(async_worker<Functionoid>),
                                         std::ref(job_to_do),
                                         std::ref(queue_mutex) );
  for( int i = 0; i < nthread; ++i )
    threads[i]->join();
#endif
  }//do_asyncronous_work(...)
  
}//namespace SpecUtilsAsync





#endif  //SpecUtilsAsync_h
