#ifndef FRAMPOPUPBOX_H
#define FRAMPOPUPBOX_H

#include <Wt/WPopupWidget>
#include <Wt/WString>
#include <Wt/WContainerWidget>
#include <Wt/WText>
#include <Wt/WWidget>
#include <Wt/WColor>
#include <Wt/WCssDecorationStyle>

namespace Wt 
{
class WContainerWidget;
class WText;
class WWidget;
}

class FRAMPopupBox : public Wt::WPopupWidget 
{
public:
  FRAMPopupBox();

  void showFor(Wt::WWidget *anchor, const Wt::WString& html);

private:
  Wt::WContainerWidget* m_box = nullptr;
  Wt::WText* m_text = nullptr;
};

#endif 
