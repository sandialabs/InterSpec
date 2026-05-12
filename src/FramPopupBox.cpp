#include "InterSpec/FramPopupBox.h"

FRAMPopupBox::FRAMPopupBox()
  : Wt::WPopupWidget(new Wt::WContainerWidget())
{
  setTransient(true, 200);
  m_box = static_cast<Wt::WContainerWidget *>(implementation());
  m_box->addStyleClass("hover-popup");
  m_box->decorationStyle().setBackgroundColor(Wt::WColor("white"));
  m_box->setPadding(8);
  m_box->setAttributeValue(
    "style",
    "background-color:white;"
    "opacity:1;"
    "border:1px solid #888;"
    "border-radius:4px;"
    "box-shadow:0 2px 10px rgba(0,0,0,0.35);"
    "z-index:10000;"
  );

  m_text = new Wt::WText(m_box);
}

void FRAMPopupBox::showFor(Wt::WWidget *anchor, const Wt::WString& html)
{
  m_text->setText(html);
  setAnchorWidget(anchor, Wt::Orientation::Horizontal);
  show();
}
