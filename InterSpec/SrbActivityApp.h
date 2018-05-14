#ifndef SrbActivityApp_h
#define SrbActivityApp_h
//  Created by wcjohns on 20110321

//#include "SrbActivityFutureConfig.h"

#include <map>
#include <string>
#include <vector>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WContainerWidget>

class SrbActivityDiv;


std::string getSrbComonResourceUrl();
std::string getSrbComonResourcePath();

class SrbActivityApp : public Wt::WApplication
{
private:
  SrbActivityDiv *m_activityDiv;

  virtual void init( const Wt::WEnvironment& env ,Wt::WContainerWidget *contentDiv );
  
public:
  SrbActivityApp( const Wt::WEnvironment& env , Wt::WContainerWidget *contentDiv);
  virtual ~SrbActivityApp();
  
  bool isMobile() const;
  bool isPhone() const;
};//class SrbActivityApp


#endif
//SrbActivityApp_h
