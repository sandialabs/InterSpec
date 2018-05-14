#ifndef SRB_HEADER_FOOTER_H
#define SRB_HEADER_FOOTER_H

//#include "SrbActivityFutureConfig.h"


#include <string>

#include <Wt/WContainerWidget>

namespace Wt
{
  class WText;
}//namespace Wt

//The SrbHeader/SrbFooter classes are largely copy-pasted from anthony nm app


class SrbHeader: public Wt::WContainerWidget
{
public:
  SrbHeader( const std::string &appTitle = "", const std::string &baseurl = "",
            Wt::WContainerWidget *parentParam = 0 );
protected:
  Wt::WText *m_headerText;
  Wt::WText *m_applicationTitle;
};//class SrbHeader


class SrbFooter: public Wt::WContainerWidget
{
public:
  SrbFooter( Wt::WContainerWidget *parentParam = 0 );
protected:
  Wt::WText *m_logoText;
  Wt::WText *m_footerText;
  Wt::WText *m_footerCloseText;
};//class SrbHeader
#endif
