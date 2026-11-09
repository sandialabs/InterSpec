#ifndef DrfLibraryWidgets_h
#define DrfLibraryWidgets_h
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

#include <map>
#include <set>
#include <memory>
#include <string>
#include <vector>
#include <cstdint>

#include <Wt/WSignal.h>
#include <Wt/WString.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/DrfLibrary.h"

class DrfChart;
class DetectorPeakResponse;

namespace Wt
{
  class WText;
  class WCheckBox;
  class WPushButton;
  class WInteractWidget;
}

/** The widgets of the Detector Efficiency dialog's Library (see `DrfSelect`): the list of detectors,
 and the preview of one.  They only display what they are given and say what the user did; the
 dialog decides what that means.  Text comes from the "DrfSelect" message bundle.
 */

namespace DrfLibraryText
{
  /** "Built-in", "GADRAS", ..., as on the Source badges. */
  Wt::WString originLabel( const DrfLibrary::Origin origin );

  /** The crystal, as best known: "HPGe · 85 × 30 mm", "NaI · Ø 7.6 cm", "HPGe", or empty. */
  std::string crystalLabel( const DrfLibrary::Summary &summary );

  /** "far-field only", "near-field OK", "fixed geometry" or "needs Monte-Carlo". */
  Wt::WString rangeLabel( const DrfLibrary::Range range );

  /** "44.7% rel. eff." (IEEE 325), or empty when not meaningful or the name already says it. */
  Wt::WString relEffLabel( const DrfLibrary::Summary &summary );

  /** One line describing a detector: crystal · rel. eff. · source · range (parts left out if unknown). */
  Wt::WString metaLine( const DrfLibrary::Summary &summary, const DrfLibrary::Origin origin );

  /** "Today", "Yesterday", "5 Nov", "5 Nov 2025" (in the user's time zone), or "—" for never. */
  Wt::WString lastUsedLabel( const int64_t utc );

  /** Markup for a detector's icon (its own, else the generic one), as a silhouette in the text colour. */
  std::string iconHtml( const DrfLibrary::Summary *summary, const double height_px, const double max_width_px );
}//namespace DrfLibraryText


/** The Library list: group headers, one row per detector (customized copies nested under the
 detector they were made from), and the "Show N customized copies" rows.

 The setters only record their value: call #rebuild once they are all set (rows come from
 `DrfLibrary::buildRows`, and each is a single button, so a rebuild is cheap).  #setPicked is the
 exception: it restyles the rows in place, unless picking the detector shows its customized copies.
 */
class DrfLibraryList : public Wt::WContainerWidget
{
public:
  DrfLibraryList();

  void setEntries( std::shared_ptr<const std::vector<DrfLibrary::Entry>> entries );
  void setFilter( const DrfLibrary::Filter &filter );
  void setPicked( const std::string &entryId );
  void setInUse( const std::string &entryId );

  /** The note beside the "Suggested for this spectrum" header (e.g. "Detective-EX200 · Serial 1234"). */
  void setSuggestedNote( const Wt::WString &note );

  /** Shows "Checking for changes…" in place of the rows (only while there are none). */
  void setChecking( const bool checking );

  /** Re-creates the rows for the current entries, filter and picked/in-use detectors. */
  void rebuild();

  /** Text searched in addition to each entry's own (e.g. its translated source label). */
  void setExtraSearchText( DrfLibrary::ExtraSearchText fcn );

  /** One-column rows for a phone: the name, then crystal · source · last used, and a chevron. */
  void setNarrow( const bool narrow );

  /** Scrolls the picked detector's row into view. */
  void scrollToPicked();

  /** The number of detector rows shown. */
  size_t numDetectorRows() const;

  /** The user clicked a detector; the entry id. */
  Wt::Signal<std::string> &picked();

  /** The user clicked "Import a detector" on the empty list. */
  Wt::Signal<> &importRequested();

protected:
  Wt::WString subText( const DrfLibrary::Entry &entry, const DrfLibrary::Row &row ) const;

  std::shared_ptr<const std::vector<DrfLibrary::Entry>> m_entries;
  DrfLibrary::Filter m_filter;
  DrfLibrary::ExtraSearchText m_extraText;
  std::map<std::string,bool> m_expanded;
  std::string m_pickedId, m_inUseId;
  Wt::WString m_suggestedNote;
  bool m_checking;
  bool m_narrow;
  size_t m_numDetectorRows;

  Wt::WContainerWidget *m_rows;
  Wt::WPushButton *m_pickedRow;
  /** The row of each detector shown; and the detectors with a "Show N customized copies" row, plus
   the copies listed under them - picking to or from one of these may show or hide copies. */
  std::map<std::string,Wt::WPushButton *> m_rowButtons;
  std::set<std::string> m_parentsWithCopies;

  Wt::Signal<std::string> m_picked;
  Wt::Signal<> m_importRequested;
};//class DrfLibraryList


/** The preview of one detector: what it is, its efficiency (or FWHM) against the detector in use,
 what it includes, where it came from, and what can be done with it.
 */
class DrfPreviewPanel : public Wt::WContainerWidget
{
public:
  /** Whether the previewed detector is in the Library, or a file being imported. */
  enum class Context { Library, Import };

  /** What is known about the previewed detector, beyond its summary. */
  struct Details
  {
    DrfLibrary::Origin origin = DrfLibrary::Origin::UserMade;
    bool inUse = false;
    /** Where it was found ("data/common_drfs.tsv"), or how it is kept ("Saved inside InterSpec"). */
    std::string where;
    /** Whether #where is a path (shown after "Found in" / "From"), rather than a sentence. */
    bool whereIsPath = true;
    /** For a customized copy, the name of the detector it was made from; empty if not a copy, or
     that detector is not in the Library. */
    std::string parentName;
    size_t numCopies = 0;      ///< customized copies of it in the Library
    bool deletable = false;    ///< a detector made by the user, kept only in the database
    bool canShare = true;      ///< XML / QR code (needs the full DRF, so not while it loads)
  };//struct Details

  DrfPreviewPanel( Context context );
  virtual ~DrfPreviewPanel();

  void setContext( const Context context );

  /** For a phone: the XML download becomes "Share…" (the download opens the system's share sheet). */
  void setNarrow( const bool narrow );

  /** Shows a detector from its summary; its curve follows with #setDetector once loaded. */
  void showSummary( std::shared_ptr<const DrfLibrary::Summary> summary, const Details &details );

  /** The full detector for the chart and the XML / QR actions (null while it loads). */
  void setDetector( std::shared_ptr<const DetectorPeakResponse> drf );

  /** The detector in use, drawn dashed when another is previewed (null for none). */
  void setReference( std::shared_ptr<const DetectorPeakResponse> drf );

  /** Shows the "nothing picked" placeholder. */
  void showNothing();

  /** A message in place of the preview's curve (e.g., why it could not be read). */
  void showError( const Wt::WString &message );

  std::shared_ptr<const DetectorPeakResponse> detector() const;

  Wt::Signal<> &customizeRequested();
  Wt::Signal<> &qrRequested();
  Wt::Signal<> &deleteRequested();
  Wt::Signal<> &parentRequested();
  /** "Use none" - on a phone, for the detector in use (the desktop has it on the in-use strip). */
  Wt::Signal<> &useNoneRequested();

protected:
  void updateChart();
  void updateLegendAndMessage();
  void setChartShowsFwhm( const bool fwhm );

  Context m_context;
  std::shared_ptr<const DrfLibrary::Summary> m_summary;
  Details m_details;
  std::shared_ptr<const DetectorPeakResponse> m_drf, m_reference;
  bool m_chartShowsFwhm;

  Wt::WContainerWidget *m_content, *m_empty;
  Wt::WText *m_emptyTitle, *m_emptyText;
  Wt::WText *m_tag, *m_icon, *m_name, *m_meta, *m_kidsNote;
  Wt::WContainerWidget *m_lineage;
  Wt::WPushButton *m_parentLink;
  Wt::WText *m_lineageEnd;
  Wt::WPushButton *m_effBtn, *m_fwhmBtn;
  Wt::WContainerWidget *m_chartHolder;
  DrfChart *m_chart;
  Wt::WText *m_chartMsg;
  Wt::WContainerWidget *m_legend, *m_legendSolid, *m_legendDashed;
  Wt::WText *m_legendSolidTxt, *m_legendDashedTxt;
  Wt::WContainerWidget *m_caps;
  Wt::WText *m_where;
  Wt::WContainerWidget *m_actions;
  Wt::WPushButton *m_customizeBtn, *m_qrBtn, *m_deleteBtn, *m_useNoneBtn;
  Wt::WInteractWidget *m_xmlBtn;  ///< a WAnchor on macOS/iOS (for the download), else a WPushButton

  Wt::Signal<> m_customizeRequested, m_qrRequested, m_deleteRequested, m_parentRequested, m_useNoneRequested;
};//class DrfPreviewPanel



/** Import › From a file, after a drop of several detectors' files (a folder, say): the detectors
 found, to pick one to preview and use, and whether to keep them all in the Library.
 */
class DrfFoundList : public Wt::WContainerWidget
{
public:
  /** One detector found. */
  struct Item
  {
    std::string name;
    Wt::WString format;     ///< what files it came from, e.g. "GADRAS · Detector.dat + Efficiency.csv"
    std::string crystal;    ///< see DrfLibraryText::crystalLabel; may be empty
    Wt::WString tag;        ///< "Ready", "Geometry only", ...
    bool usable = true;     ///< false styles the tag as a warning
  };

  DrfFoundList();

  /** Shows "Reading N files…" while the drop is read. */
  void setReading( const Wt::WString &message );

  /** Shows the detectors found.  `canKeep` shows the "Keep this folder in my Library" option. */
  void setFound( const Wt::WString &title, const Wt::WString &subtitle, const Wt::WString &summary,
                 const std::vector<Item> &items, const bool canKeep );

  void setPicked( const size_t index );
  bool keep() const;

  /** What "Keep this folder in my Library" will do (copy it, or list it where it is). */
  void setKeepDescription( const Wt::WString &description );

  Wt::Signal<size_t> &picked();
  Wt::Signal<> &cleared();

protected:
  Wt::WText *m_title, *m_subtitle, *m_summary, *m_hint;
  Wt::WContainerWidget *m_rows;
  Wt::WContainerWidget *m_keepRow;
  Wt::WCheckBox *m_keep;
  Wt::WText *m_keepDesc;
  std::vector<Wt::WPushButton *> m_items;

  Wt::Signal<size_t> m_picked;
  Wt::Signal<> m_cleared;
};//class DrfFoundList


/** Library › Locations: the folders and files the Library is read from - found automatically, and
 added by the user - with what was found in each, and the controls to manage them.
 */
class DrfLocationsPage : public Wt::WContainerWidget
{
public:
  DrfLocationsPage();

  /** The locations the last scan used, what it found in each (by `Location::id`), and how many
   detectors are saved inside InterSpec ("Made by you").
   */
  void setLocations( const std::vector<DrfLibrary::Location> &locations,
                     const std::vector<DrfLibrary::LocationStatus> &statuses,
                     const size_t numMadeByYou );

  Wt::Signal<> &backRequested();
  /** A location's "Show in Library" was changed: its id, and whether it is now shown. */
  Wt::Signal<std::string,bool> &shownChanged();
  Wt::Signal<std::string> &rescanRequested();
  Wt::Signal<std::string> &locateRequested();
  Wt::Signal<std::string> &removeRequested();
  Wt::Signal<> &addFolderRequested();
  Wt::Signal<> &addFileRequested();

protected:
  Wt::WContainerWidget *m_autoRows, *m_userRows;
  Wt::WText *m_yours;

  Wt::Signal<> m_backRequested;
  Wt::Signal<std::string,bool> m_shownChanged;
  Wt::Signal<std::string> m_rescanRequested, m_locateRequested, m_removeRequested;
  Wt::Signal<> m_addFolderRequested, m_addFileRequested;
};//class DrfLocationsPage

#endif //DrfLibraryWidgets_h
