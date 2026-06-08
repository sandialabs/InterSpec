#ifndef FRAMPOPUPDELEGATE_H
#define FRAMPOPUPDELEGATE_H

#include <Wt/WItemDelegate>
#include <Wt/WModelIndex>
#include <Wt/WString>
#include <Wt/WWidget>
#include <Wt/WWebWidget>

#include <nlohmann/json.hpp>
#include <boost/any.hpp>

class FRAMPopupBox;

class FRAMPopupDelegate : public Wt::WItemDelegate 
{
public:
  FRAMPopupDelegate(FRAMPopupBox& popup, int hoverColumn);

  Wt::WWidget *update(
                       Wt::WWidget *widget,
                       const Wt::WModelIndex& index,
                       Wt::WFlags<Wt::ViewItemRenderFlag> flags) override;

private:
  FRAMPopupBox& m_popup;
  int m_hoverColumn;
};


#endif