#include "InterSpec/FramPopupDelegate.h"
#include "InterSpec/FRAMPopupBox.h"

namespace
{
  static const int IsotopicsJsonRole = Wt::UserRole + 100;
}

FRAMPopupDelegate::FRAMPopupDelegate(
                                      FRAMPopupBox& popup,
                                      int hoverColumn)
                                      : m_popup(popup),
                                        m_hoverColumn(hoverColumn)
{
  //nothing to do
}

Wt::WWidget* FRAMPopupDelegate::update(
                                Wt::WWidget *widget,
                                const Wt::WModelIndex& index,
                                Wt::WFlags<Wt::ViewItemRenderFlag> flags)
{
  Wt::WWidget *cell = Wt::WItemDelegate::update(widget, index, flags);

  if (index.column() == m_hoverColumn && cell) 
  {
    Wt::WInteractWidget *iw = dynamic_cast<Wt::WInteractWidget *>(cell);
    if (iw) 
    {
      iw->mouseWentOver().connect(
        [this, cell, index](const Wt::WMouseEvent&) {
          Wt::WString value;
          try 
          {
            value = boost::any_cast<Wt::WString>(index.data(Wt::DisplayRole));
          }
          catch (const boost::bad_any_cast&) 
          {
            value = "unknown";
          }

          Wt::WString isotopics_json;
          try
          {
            isotopics_json = boost::any_cast<Wt::WString>(
              index.data( IsotopicsJsonRole ) );
          }
          catch( const boost::bad_any_cast& )
          {
            isotopics_json = "";
          }

      std::string popupHtml =
        "<div style='min-width:250px'>"
        "<h3>FRAM Results</h3>";
      try
      {
        nlohmann::json results = nlohmann::json::parse(isotopics_json.toUTF8());
        const nlohmann::json* res = &results;
        if( results.is_array() )
        {
          res = nullptr;
          for( const auto &entry : results )
          {
            if( entry.is_object() && entry.contains("Isotopics") )
            {
              res = &entry;
              break;
            }
          }
        }
        if( !res || !res->is_object() || !res->contains("Isotopics") )
          throw std::runtime_error("FRAM JSON does not contain Isotopics");

        std::string parameterSet;
        if( res->contains("Header") && (*res)["Header"].is_object() )
          parameterSet = (*res)["Header"].value("parameter_set", "");
        
        popupHtml +=
          "<div><b>Parameter Set:</b> " + parameterSet + "</div><br/>";
        
        popupHtml += "<table style='border-collapse:collapse;'>";
        
        popupHtml +=
          "<tr>"
          "<th style='padding:4px 12px;text-align:left;'></th>";
      
        for( const auto &nuc : (*res)["Isotopics"] )
        {
          popupHtml +=
            "<th style='padding:4px 12px;text-align:center;'>"
            + nuc.value("nuclide", std::string("")) +
            "</th>";
        }

        popupHtml += "</tr>";
        auto addRow = [&]( const char *label, const char *jsonKey, bool percent )
        {
          popupHtml +=
            "<tr><td style='padding:4px 12px;font-weight:bold;'>"
            + std::string(label) +
            "</td>";
      
          for( const auto &nuc : (*res)["Isotopics"] )
          {
            const double value = nuc.value(jsonKey, 0.0);
            char buff[64];
            snprintf(buff, sizeof(buff), "%.4f", value);
            popupHtml +=
              "<td style='padding:4px 12px;text-align:left;'>"
              + std::string(buff)
              + (percent ? "%" : "")
              + "</td>";
          }
          popupHtml += "</tr>";
        };

        addRow("Mass %",          "mass_percent",    true);
        addRow("Sigma (stat)",    "mass_sigma_stat", false);
        addRow("RSD (stat)",      "mass_rsd_stat",   true);
        addRow("Sigma (stat+sys)","mass_sigma_sys",  false);
        addRow("RSD (stat+sys)",  "mass_rsd_sys",    true);

        popupHtml += "</table>";
      }

      catch( const std::exception &e )
      {
        popupHtml += "<div>Failed to parse FRAM JSON: " + std::string(e.what()) + "</div>";
      }

      popupHtml += "</div>";

      m_popup.showFor(cell, popupHtml);
                                                    });
      iw->mouseWentOut().connect( [this](const Wt::WMouseEvent&) {m_popup.hide();});
    }
  }

  return cell;
}


