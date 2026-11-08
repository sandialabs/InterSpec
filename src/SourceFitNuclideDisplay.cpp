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

#include <cmath>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <utility>
#include <iostream>
#include <optional>
#include <algorithm>
#include <functional>

#include <Wt/WAny.h>
#include <Wt/WText.h>
#include <Wt/WColor.h>
#include <Wt/WLineEdit.h>
#include <Wt/WCheckBox.h>
#include <Wt/WModelIndex.h>
#include <Wt/WPushButton.h>
#include <Wt/WContainerWidget.h>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/DecayWindow.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ShieldingSourceDisplay.h"   //for SourceFitModel
#include "InterSpec/ShieldingSourceFitCalc.h"   //for ModelSourceType
#include "InterSpec/SourceFitNuclideDisplay.h"

using namespace Wt;
using namespace std;

namespace
{
  /** RAII helper that marks the display as "refreshing" so edit handlers ignore spurious change
   signals and reconcile cannot recurse while we programmatically update widgets. */
  struct RefreshGuard
  {
    bool &m_flag;
    const bool m_prev;
    RefreshGuard( bool &flag ) : m_flag(flag), m_prev(flag){ m_flag = true; }
    ~RefreshGuard(){ m_flag = m_prev; }
  };


  /** PhysicalUnits::printToBestActivityUnits, but also for negative values (which a back-decayed
   activity can be, when a fitted ancestor's in-growth exceeds the fitted activity). */
  string print_signed_activity( const double activity, const int ndecimal, const bool curies )
  {
    return string( (activity < 0.0) ? "-" : "" )
           + PhysicalUnits::printToBestActivityUnits( fabs(activity), ndecimal, curies );
  }


  /** What to show for a model row's activity at the back-decay reference time. */
  struct BackDecayText
  {
    /** False when back-decay is not enabled (or there is no result for the row). */
    bool show = false;

    /** True if the activity was determined, in which case `value` is set. */
    bool ok = false;

    /** E.g., "2.350 mCi ± 0.04 mCi", in the same units as the rows activity. */
    string value;

    /** The reference time, and any note or error. */
    WString tooltip;
  };//struct BackDecayText


  BackDecayText back_decay_text( const SourceFitModel &model, const int row,
                                 const ShieldingSourceFitCalc::BackDecayResults *bd )
  {
    BackDecayText answer;
    if( !bd )
      return answer;

    if( !bd->error.empty() )
    {
      answer.show = true;
      answer.tooltip = WString::fromUTF8( bd->error );
      return answer;
    }//if( the whole calculation failed )

    const SandiaDecay::Nuclide * const nuc = model.nuclide( row );
    const ShieldingSourceFitCalc::BackDecaySourceResult * const res = nuc ? bd->source( nuc ) : nullptr;
    if( !res )
      return answer;

    answer.show = true;

    const string when = ShieldingSourceFitCalc::back_decay_time_str( bd->ref_time );
    const char *key = (bd->time_before_meas >= 0.0) ? "sfn-tt-act-tref-before" : "sfn-tt-act-tref-after";
    answer.tooltip = WString::tr( key )
                .arg( when.empty() ? WString::tr("sfn-tref-unknown-date") : WString::fromUTF8(when) )
                .arg( ShieldingSourceDisplay::backDecayDurationText( fabs(bd->time_before_meas) ) )
                .arg( SpecUtils::printCompact( fabs(res->half_lives), 3 ) )
                .arg( nuc->symbol );
    if( !res->note.empty() )
      answer.tooltip = WString("{1} {2}").arg( answer.tooltip ).arg( WString::fromUTF8(res->note) );

    if( res->status != ShieldingSourceFitCalc::BackDecayStatus::Ok )
      return answer;

    answer.ok = true;

    // In the same units as the card's activity (the model holds total activities, also for trace sources).
    const double value = res->activity;
    const std::optional<double> &uncert = res->activity_uncertainty;

    const bool curies = model.displayCuries();
    const string postfix = DetectorPeakResponse::det_eff_geom_type_postfix( model.detType() );
    answer.value = print_signed_activity( value, 3, curies ) + postfix;
    if( uncert.has_value() && (*uncert > 0.0) )
      answer.value += " \xC2\xB1 " + PhysicalUnits::printToBestActivityUnits( *uncert, 2, curies ) + postfix;

    return answer;
  }//BackDecayText back_decay_text(...)
}//namespace


/** One per-nuclide card.  Reads/writes everything through #SourceFitModel. */
class SourceFitNuclideDisplay::NuclideCard : public Wt::WContainerWidget
{
public:
  NuclideCard( SourceFitNuclideDisplay *disp, const SandiaDecay::Nuclide *nuc )
   : WContainerWidget(),
     m_disp( disp ),
     m_model( disp->m_model ),
     m_nuc( nuc ),
     m_nestPeaks( disp->m_nestPeaks ),
     m_open( disp->m_nestPeaks ),   //phone accordion cards start expanded (showing all peak rows)
     m_swatch( nullptr ),
     m_countLabel( nullptr ),
     m_chevron( nullptr ),
     m_summaryRow( nullptr ),
     m_bodyDiv( nullptr ),
     m_peakList( nullptr ),
     m_peakListLbl( nullptr )
  {
    addStyleClass( "NucCard" );

    // ---- header: (swatch) name, half-life, source-type tag, (count + chevron) ----
    WContainerWidget *head = addNew<WContainerWidget>();
    head->addStyleClass( "NucHead" );
    if( m_nestPeaks )
    {
      // The whole header toggles the accordion open/closed.
      head->addStyleClass( "NucHeadTap" );
      head->clicked().connect( this, &NuclideCard::toggleOpen );
      m_swatch = head->addNew<WContainerWidget>();
      m_swatch->addStyleClass( "NucSwatch" );
    }//if( m_nestPeaks )

    m_nameLabel = head->addNew<WText>();
    m_nameLabel->addStyleClass( "NucName" );
    m_halfLifeLabel = head->addNew<WText>();
    m_halfLifeLabel->addStyleClass( "NucHl" );
    m_tag = head->addNew<WText>();
    m_tag->addStyleClass( "NucTag" );
    m_tag->hide();

    if( m_nestPeaks )
    {
      WContainerWidget *grow = head->addNew<WContainerWidget>();
      grow->addStyleClass( "NucGrow" );
      m_countLabel = head->addNew<WText>();
      m_countLabel->addStyleClass( "NucCount" );
      m_chevron = head->addNew<WText>();
      m_chevron->addStyleClass( "NucChev" );
      m_chevron->setText( "›" );  //single right-pointing angle quotation mark

      // Collapsed one-line summary of which peaks are in the fit.
      m_summaryRow = addNew<WContainerWidget>();
      m_summaryRow->addStyleClass( "NucPeakSummary" );
    }//if( m_nestPeaks )

    // In nest-peaks mode the activity/age fields live in a collapsible body; otherwise they sit
    //  directly on the card (desktop behavior, unchanged).
    WContainerWidget *fields = this;
    if( m_nestPeaks )
    {
      m_bodyDiv = addNew<WContainerWidget>();
      m_bodyDiv->addStyleClass( "NucBody" );
      fields = m_bodyDiv;
    }//if( m_nestPeaks )

    // ---- activity row ----
    WContainerWidget *actRow = fields->addNew<WContainerWidget>();
    actRow->addStyleClass( "NucRow" );
    WText *actLbl = actRow->addNew<WText>( WString::tr("sfn-activity") );
    actLbl->addStyleClass( "NucRowLbl" );
    m_activityEdit = actRow->addNew<WLineEdit>();
    m_activityEdit->addStyleClass( "NucField" );
    m_activityEdit->changed().connect( this, &NuclideCard::onActivityChanged );
    m_activityEdit->enterPressed().connect( this, &NuclideCard::onActivityChanged );
    m_fitActivityCB = actRow->addNew<WCheckBox>( WString::tr("sfn-fit") );
    m_fitActivityCB->addStyleClass( "NucFitCb" );
    m_fitActivityCB->changed().connect( this, &NuclideCard::onFitActivityToggled );

    // ---- mass + act(T_ref) line ----
    WContainerWidget *massLine = fields->addNew<WContainerWidget>();
    massLine->addStyleClass( "NucMass" );
    m_massLabel = massLine->addNew<WText>();
    m_actRefLabel = massLine->addNew<WText>();
    m_actRefLabel->addStyleClass( "NucTRef" );
    m_actRefLabel->hide();

    // ---- age row ----
    m_ageRow = fields->addNew<WContainerWidget>();
    m_ageRow->addStyleClass( "NucRow" );
    WText *ageLbl = m_ageRow->addNew<WText>( WString::tr("sfn-age") );
    ageLbl->addStyleClass( "NucRowLbl" );
    m_ageEdit = m_ageRow->addNew<WLineEdit>();
    m_ageEdit->addStyleClass( "NucField" );
    m_ageEdit->changed().connect( this, &NuclideCard::onAgeChanged );
    m_ageEdit->enterPressed().connect( this, &NuclideCard::onAgeChanged );
    m_fitAgeCB = m_ageRow->addNew<WCheckBox>( WString::tr("sfn-fit") );
    m_fitAgeCB->addStyleClass( "NucFitCb" );
    m_fitAgeCB->changed().connect( this, &NuclideCard::onFitAgeToggled );

    // ---- progeny attribution row (rebuilt each refresh) ----
    m_progenyRow = fields->addNew<WContainerWidget>();
    m_progenyRow->addStyleClass( "NucProgeny" );
    m_progenyRow->hide();

    // ---- nested peak list (phone): a header + one toggleable row per peak ----
    if( m_nestPeaks )
    {
      m_peakList = addNew<WContainerWidget>();
      m_peakList->addStyleClass( "NucPeakList" );
      WContainerWidget *plHead = m_peakList->addNew<WContainerWidget>();
      plHead->addStyleClass( "NucPeakListHead" );
      m_peakListLbl = plHead->addNew<WText>();
      m_peakListLbl->addStyleClass( "NpListLbl" );
    }//if( m_nestPeaks )

    refresh();
  }//NuclideCard constructor


  const SandiaDecay::Nuclide *nuclide() const { return m_nuc; }

  bool isOpen() const { return m_open; }
  void setOpen( const bool open ){ m_open = open; applyOpenState(); }


  /** Re-reads the model row for this nuclide and updates all sub-widgets in place. */
  void refresh()
  {
    const int row = m_model->nuclideIndex( m_nuc );
    if( row < 0 )
      return;

    RefreshGuard guard( m_disp->m_refreshing );

    // Header
    m_nameLabel->setText( WString::fromUTF8(m_nuc->symbol) );
    m_halfLifeLabel->setText( WString("{1}={2}").arg( WString::tr("T1/2") )
        .arg( PhysicalUnitsLocalized::printToBestTimeUnits( m_nuc->halfLife, 2, PhysicalUnits::second ) ) );

    if( m_swatch )
    {
      const string colorCss = m_disp->peakColorForNuclide( m_nuc );
      m_swatch->setAttributeValue( "style",
                          "background:" + (colorCss.empty() ? string("#888") : colorCss) + ";" );
    }//if( m_swatch )

    const ShieldingSourceFitCalc::ModelSourceType srcType = m_model->sourceType( row );
    const bool isPoint = (srcType == ShieldingSourceFitCalc::ModelSourceType::Point);

    if( !isPoint )
    {
      const char *tagKey = (srcType == ShieldingSourceFitCalc::ModelSourceType::Trace)
                            ? "sfn-tag-trace" : "sfn-tag-shield-src";
      m_tag->setText( WString::tr(tagKey) );
      m_tag->show();
    }else
    {
      m_tag->hide();
    }

    // Activity
    const WModelIndex actIdx = m_model->index( m_nuc, SourceFitModel::kActivity );
    m_activityEdit->setText( Wt::asString( m_model->data(actIdx) ) );
    m_activityEdit->setReadOnly( !isPoint );
    m_activityEdit->toggleStyleClass( "NucFieldRo", !isPoint );

    const WModelIndex fitActIdx = m_model->index( m_nuc, SourceFitModel::kFitActivity );
    const cpp17::any fitActAny = m_model->data( fitActIdx, ItemDataRole::Checked );
    if( fitActAny.has_value() )
    {
      m_fitActivityCB->show();
      m_fitActivityCB->setChecked( cpp17::any_cast<bool>(fitActAny) );
    }else
    {
      m_fitActivityCB->hide();
    }

    // Mass + act(T_ref)
    const bool fitAct = m_model->fitActivity( row );
    m_massLabel->setText( WString("{1} {2}")
         .arg( WString::tr( fitAct ? "sfn-mass-approx" : "sfn-mass" ) )
         .arg( Wt::asString( m_model->data( m_model->index(m_nuc, SourceFitModel::kIsotopeMass) ) ) ) );

    // Activity at the back-decay reference time, when that option is enabled.
    const BackDecayText bdtxt = back_decay_text( *m_model, row, m_disp->m_backDecay.get() );
    m_actRefLabel->setHidden( !bdtxt.show );
    if( bdtxt.show )
    {
      if( bdtxt.ok )
        m_actRefLabel->setText( WString("{1}={2}").arg( WString::tr("sfn-act-tref") )
                                                  .arg( WString::fromUTF8(bdtxt.value) ) );
      else
        m_actRefLabel->setText( WString("{1}: {2}").arg( WString::tr("sfn-act-tref") )
                                                   .arg( WString::tr("sfn-act-tref-undetermined") ) );
      m_actRefLabel->toggleStyleClass( "ErrorTxt", !bdtxt.ok );
      m_actRefLabel->setToolTip( bdtxt.tooltip );
    }//if( bdtxt.show )

    const bool ageable = !PeakDef::ageFitNotAllowed( m_nuc );

    // Age row
    if( ageable )
    {
      m_ageRow->show();
      const WModelIndex ageIdx = m_model->index( m_nuc, SourceFitModel::kAge );
      m_ageEdit->setText( Wt::asString( m_model->data(ageIdx) ) );
      const bool ageEditable = m_model->flags(ageIdx).test( ItemFlag::Editable );
      m_ageEdit->setReadOnly( !ageEditable );
      m_ageEdit->toggleStyleClass( "NucFieldRo", !ageEditable );

      const WModelIndex fitAgeIdx = m_model->index( m_nuc, SourceFitModel::kFitAge );
      const cpp17::any fitAgeAny = m_model->data( fitAgeIdx, ItemDataRole::Checked );
      if( fitAgeAny.has_value() )
      {
        m_fitAgeCB->show();
        m_fitAgeCB->setChecked( cpp17::any_cast<bool>(fitAgeAny) );
      }else
      {
        m_fitAgeCB->hide();
      }
    }else
    {
      m_ageRow->hide();
    }

    updateProgeny();

    if( m_nestPeaks )
      updatePeakList();
  }//void refresh()


protected:
  /** Toggles this accordion card open/closed (phone nest-peaks mode). */
  void toggleOpen()
  {
    m_open = !m_open;
    applyOpenState();
  }

  /** Shows/hides the collapsible body + peak list and rotates the chevron to match #m_open. */
  void applyOpenState()
  {
    if( !m_nestPeaks )
      return;
    if( m_bodyDiv )    m_bodyDiv->setHidden( !m_open );
    if( m_peakList )   m_peakList->setHidden( !m_open );
    if( m_summaryRow ) m_summaryRow->setHidden( m_open );
    if( m_chevron )    m_chevron->toggleStyleClass( "NucChevOpen", m_open );
  }

  /** Gathers the peaks attributed to this nuclide (both used and unused, so the user can toggle
   them), in ascending energy order, paired with their PeakModel row index. */
  std::vector<std::pair<int,PeakModel::PeakShrdPtr>> peaksForNuclide() const
  {
    std::vector<std::pair<int,PeakModel::PeakShrdPtr>> out;
    const std::shared_ptr<PeakModel> &pm = m_disp->m_peakModel;
    if( !pm )
      return out;
    const int npeaks = static_cast<int>( pm->npeaks() );
    for( int row = 0; row < npeaks; ++row )
    {
      const PeakModel::PeakShrdPtr p = pm->peakPtr( row );
      if( p && (p->parentNuclide() == m_nuc) )
        out.push_back( std::make_pair( row, p ) );
    }
    std::sort( begin(out), end(out), []( const std::pair<int,PeakModel::PeakShrdPtr> &a,
                                         const std::pair<int,PeakModel::PeakShrdPtr> &b ){
      return a.second->mean() < b.second->mean();
    } );
    return out;
  }//peaksForNuclide()

  /** Rebuilds the nested peak list, the "sel/total" count, and the collapsed summary line. */
  void updatePeakList()
  {
    if( !m_peakList )
      return;

    const std::vector<std::pair<int,PeakModel::PeakShrdPtr>> peaks = peaksForNuclide();
    int nused = 0;
    std::vector<int> usedEnergies;
    for( const auto &rp : peaks )
    {
      if( rp.second->useForShieldingSourceFit() )
      {
        ++nused;
        usedEnergies.push_back( static_cast<int>( std::round( rp.second->mean() ) ) );
      }
    }//for( const auto &rp : peaks )

    const int total = static_cast<int>( peaks.size() );

    // Header count + summary line.
    if( m_countLabel )
    {
      m_countLabel->setText( WString("{1}/{2}").arg(nused).arg(total) );
      m_countLabel->toggleStyleClass( "NucCountHas", nused > 0 );
    }

    if( m_summaryRow )
    {
      if( nused > 0 )
      {
        string es;
        for( size_t i = 0; i < usedEnergies.size(); ++i )
          es += (i ? ", " : "") + std::to_string( usedEnergies[i] );
        m_summaryRow->clear();
        WText *t = m_summaryRow->addNew<WText>( WString::tr("sfn-peaks-in-fit").arg( WString::fromUTF8(es) ) );
        t->setTextFormat( TextFormat::XHTML );
      }else
      {
        m_summaryRow->clear();
        m_summaryRow->addNew<WText>( WString::tr("sfn-no-peaks-selected") );
      }
    }//if( m_summaryRow )

    // Peak-list label, e.g. "Peaks · 2 of 5 used".
    if( m_peakListLbl )
      m_peakListLbl->setText( WString::tr("sfn-peaks-n-of-m").arg(nused).arg(total) );

    // Rebuild rows.  We keep the header (first child) and drop the rest.
    while( m_peakList->count() > 1 )
      m_peakList->removeWidget( m_peakList->widget(1) );

    for( const auto &rp : peaks )
    {
      const int peakRow = rp.first;
      const PeakModel::PeakShrdPtr &p = rp.second;
      const bool use = p->useForShieldingSourceFit();

      WContainerWidget *row = m_peakList->addNew<WContainerWidget>();
      row->addStyleClass( "NpRow" );
      if( use )
        row->addStyleClass( "NpRowOn" );
      row->clicked().connect( this, [this,peakRow](){ onPeakToggled( peakRow ); } );

      WContainerWidget *box = row->addNew<WContainerWidget>();
      box->addStyleClass( "NpBox" );

      WContainerWidget *main = row->addNew<WContainerWidget>();
      main->addStyleClass( "NpMain" );

      char en[32] = { '\0' };
      snprintf( en, sizeof(en), "%.1f", p->mean() );
      WText *enTxt = main->addNew<WText>( WString("{1}<span class=\"NpUnit\">{2}</span>")
                                          .arg( WString::fromUTF8(en) ).arg( WString::tr("sfn-kev") ) );
      enTxt->setTextFormat( TextFormat::XHTML );
      enTxt->addStyleClass( "NpEnergy" );

      // Progeny chip: the gamma is actually emitted by a daughter of this nuclide.
      const SandiaDecay::Transition * const trans = p->nuclearTransition();
      if( trans && trans->parent && (trans->parent != m_nuc) )
      {
        WText *prog = main->addNew<WText>( WString::fromUTF8(trans->parent->symbol) );
        prog->addStyleClass( "NpProg" );
      }

      // X-ray / escape subtext.
      const SandiaDecay::RadParticle * const part = p->decayParticle();
      const char *subKey = nullptr;
      if( part && (part->type == SandiaDecay::XrayParticle) )
        subKey = "sfn-peak-xray";
      else if( (p->sourceGammaType() == PeakDef::SingleEscapeGamma)
              || (p->sourceGammaType() == PeakDef::DoubleEscapeGamma) )
        subKey = "sfn-peak-escape";
      if( subKey )
      {
        WText *sub = main->addNew<WText>( WString::tr(subKey) );
        sub->addStyleClass( "NpSub" );
      }
    }//for( const auto &rp : peaks )

    applyOpenState();
  }//void updatePeakList()

  /** Toggles whether the given peak (by PeakModel row) is used in the shielding/source fit. */
  void onPeakToggled( const int peakRow )
  {
    if( m_disp->m_refreshing )
      return;
    const std::shared_ptr<PeakModel> &pm = m_disp->m_peakModel;
    if( !pm )
      return;
    const PeakModel::PeakShrdPtr p = pm->peakPtr( peakRow );
    if( !p )
      return;
    const WModelIndex idx = pm->index( peakRow, PeakModel::kUseForShieldingSourceFit );
    pm->setData( idx, !p->useForShieldingSourceFit() );
    // The peak-model change will fan out to SourceFitModel and back to us to refresh.
  }//void onPeakToggled( const int peakRow )

  void onActivityChanged()
  {
    if( m_disp->m_refreshing )
      return;
    const bool ok = m_model->setData( m_model->index(m_nuc, SourceFitModel::kActivity),
                                      cpp17::any( m_activityEdit->text() ) );
    if( !ok )
      refresh();  // revert to the model value on a parse failure
  }

  void onFitActivityToggled()
  {
    if( m_disp->m_refreshing )
      return;
    m_model->setData( m_model->index(m_nuc, SourceFitModel::kFitActivity),
                      cpp17::any( m_fitActivityCB->isChecked() ), ItemDataRole::Checked );
  }

  void onAgeChanged()
  {
    if( m_disp->m_refreshing )
      return;
    const bool ok = m_model->setData( m_model->index(m_nuc, SourceFitModel::kAge),
                                      cpp17::any( m_ageEdit->text() ) );
    if( !ok )
      refresh();
  }

  void onFitAgeToggled()
  {
    if( m_disp->m_refreshing )
      return;
    m_model->setData( m_model->index(m_nuc, SourceFitModel::kFitAge),
                      cpp17::any( m_fitAgeCB->isChecked() ), ItemDataRole::Checked );
  }

  /** Rebuilds the "via progeny" row: for each peak attributed to this nuclide but actually
   emitted by a daughter, list the daughter and the gamma energies. */
  void updateProgeny()
  {
    m_progenyRow->clear();

    map<const SandiaDecay::Nuclide *, vector<float>> byDaughter;
    shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks
                        = m_disp->m_peakModel ? m_disp->m_peakModel->peaks() : nullptr;
    if( peaks )
    {
      for( const PeakModel::PeakShrdPtr &p : *peaks )
      {
        if( !p || !p->useForShieldingSourceFit() || (p->parentNuclide() != m_nuc) )
          continue;
        const SandiaDecay::Transition * const trans = p->nuclearTransition();
        if( !trans || !trans->parent || (trans->parent == m_nuc) )
          continue;
        const SandiaDecay::RadParticle * const part = p->decayParticle();
        const float energy = part ? part->energy : static_cast<float>( p->mean() );
        byDaughter[trans->parent].push_back( energy );
      }//for( loop over peaks )
    }//if( peaks )

    if( byDaughter.empty() )
    {
      m_progenyRow->hide();
      return;
    }

    m_progenyRow->show();
    WText *lbl = m_progenyRow->addNew<WText>( WString::tr("sfn-via-progeny") );
    lbl->addStyleClass( "NpLbl" );
    WContainerWidget *list = m_progenyRow->addNew<WContainerWidget>();
    list->addStyleClass( "NpList" );

    for( const auto &kv : byDaughter )
    {
      WContainerWidget *chip = list->addNew<WContainerWidget>();
      chip->addStyleClass( "NpChip" );
      WText *d = chip->addNew<WText>( WString::fromUTF8(kv.first->symbol) );
      d->addStyleClass( "NpDaughter" );

      string es;
      for( size_t i = 0; i < kv.second.size(); ++i )
      {
        char b[32];
        snprintf( b, sizeof(b), "%.1f", kv.second[i] );
        es += (i ? ", " : "") + string(b);
      }
      es += " keV";
      WText *e = chip->addNew<WText>( WString::fromUTF8(es) );
      e->addStyleClass( "NpEn" );
    }//for( loop over daughters )
  }//void updateProgeny()


  SourceFitNuclideDisplay *m_disp;
  SourceFitModel *m_model;
  const SandiaDecay::Nuclide *m_nuc;

  /** Phone "merged" mode: card is an accordion that nests its peak rows. */
  const bool m_nestPeaks;
  bool m_open;

  Wt::WContainerWidget *m_swatch;     //!< peak-color swatch (nestPeaks only)
  Wt::WText *m_nameLabel;
  Wt::WText *m_halfLifeLabel;
  Wt::WText *m_tag;
  Wt::WText *m_countLabel;            //!< "sel/total" peak count (nestPeaks only)
  Wt::WText *m_chevron;               //!< expand/collapse chevron (nestPeaks only)
  Wt::WContainerWidget *m_summaryRow; //!< collapsed peak summary line (nestPeaks only)
  Wt::WContainerWidget *m_bodyDiv;    //!< collapsible body wrapping the fields (nestPeaks only)
  Wt::WLineEdit *m_activityEdit;
  Wt::WCheckBox *m_fitActivityCB;
  Wt::WText *m_massLabel;
  Wt::WText *m_actRefLabel;           //!< "act(T_ref)=..." - activity at the back-decay reference time
  Wt::WContainerWidget *m_ageRow;
  Wt::WLineEdit *m_ageEdit;
  Wt::WCheckBox *m_fitAgeCB;
  Wt::WContainerWidget *m_progenyRow;
  Wt::WContainerWidget *m_peakList;   //!< nested peak rows (nestPeaks only)
  Wt::WText *m_peakListLbl;           //!< "Peaks · N of M used" (nestPeaks only)
};//class SourceFitNuclideDisplay::NuclideCard



SourceFitNuclideDisplay::SourceFitNuclideDisplay( SourceFitModel *model,
                                                  std::shared_ptr<PeakModel> peakModel,
                                                  InterSpec *interspec,
                                                  const bool nestPeaks )
  : WContainerWidget(),
    m_model( model ),
    m_peakModel( peakModel ),
    m_interspec( interspec ),
    m_nestPeaks( nestPeaks ),
    m_cardsHolder( nullptr ),
    m_backDecayFoot( nullptr ),
    m_backDecayFootTxt( nullptr ),
    m_backDecayFootSpan( nullptr ),
    m_showDecayBtn( nullptr ),
    m_refreshing( false )
{
  addStyleClass( "SourceNuclideDisplay" );
  if( m_nestPeaks )
    addStyleClass( "SourceNuclideDisplayPhone" );

  m_cardsHolder = addNew<WContainerWidget>();
  m_cardsHolder->addStyleClass( "NucCardsHolder" );

  // Only shown while back-decay is in use: the reference time, and a link that opens the decay
  //  tool from it (see updateBackDecayFooter).
  m_backDecayFoot = addNew<WContainerWidget>();
  m_backDecayFoot->addStyleClass( "NucFoot" );
  WContainerWidget * const footInfo = m_backDecayFoot->addNew<WContainerWidget>();
  footInfo->addStyleClass( "NucFootInfo" );
  m_backDecayFootTxt = footInfo->addNew<WText>();
  m_backDecayFootTxt->addStyleClass( "NucFootRef" );
  m_backDecayFootSpan = footInfo->addNew<WText>( WString(), TextFormat::Plain );
  m_backDecayFootSpan->addStyleClass( "NucFootSpan" );
  m_showDecayBtn = m_backDecayFoot->addNew<WPushButton>( WString::tr("sfn-show-decay") );
  m_showDecayBtn->addStyleClass( "SsdLink ShowDecayLink" );
  m_showDecayBtn->setToolTip( WString::tr("sfn-tt-show-decay") );
  m_showDecayBtn->clicked().connect( this, &SourceFitNuclideDisplay::showDecayFromRefTime );
  m_backDecayFoot->hide();

  // Observe the model; SourceFitModel stays the single source of truth.
  model->dataChanged().connect( this, &SourceFitNuclideDisplay::handleModelChanged );
  model->rowsInserted().connect( this, &SourceFitNuclideDisplay::handleModelChanged );
  model->rowsRemoved().connect( this, &SourceFitNuclideDisplay::handleModelChanged );
  model->layoutChanged().connect( this, &SourceFitNuclideDisplay::handleModelChanged );
  model->modelReset().connect( this, &SourceFitNuclideDisplay::handleModelChanged );

  // In nest-peaks mode the cards show every peak of a nuclide (used and unused) so the user can
  //  toggle them; toggling a currently-unused peak may not move a SourceFitModel row, so we also
  //  observe the peak model directly and refresh the nested peak lists.
  if( m_nestPeaks && m_peakModel )
  {
    m_peakModel->dataChanged().connect( this, &SourceFitNuclideDisplay::handlePeakModelChanged );
    m_peakModel->rowsInserted().connect( this, &SourceFitNuclideDisplay::handlePeakModelChanged );
    m_peakModel->rowsRemoved().connect( this, &SourceFitNuclideDisplay::handlePeakModelChanged );
    m_peakModel->layoutChanged().connect( this, &SourceFitNuclideDisplay::handlePeakModelChanged );
    m_peakModel->modelReset().connect( this, &SourceFitNuclideDisplay::handlePeakModelChanged );
  }//if( m_nestPeaks && m_peakModel )

  reconcileCards();
}//SourceFitNuclideDisplay constructor


SourceFitNuclideDisplay::~SourceFitNuclideDisplay()
{
  // m_cards point to NuclideCards owned by m_cardsHolder (a child widget), so they are torn
  //  down automatically; nothing to do here.
}


void SourceFitNuclideDisplay::handleModelChanged()
{
  reconcileCards();
}


void SourceFitNuclideDisplay::setBackDecayProvider( BackDecayProvider provider )
{
  m_backDecayProvider = std::move( provider );
  m_backDecay.reset();

  // Clearing the provider (as on destruction of its owner) does not refresh the cards.
  if( m_backDecayProvider )
    reconcileCards();
}//void setBackDecayProvider( BackDecayProvider provider )


void SourceFitNuclideDisplay::updateBackDecayFooter()
{
  // (An error - e.g., a date/time with a spectrum that has no start time - is shown by the option's
  //  status text and on each card, so there is nothing for the footer to offer.)
  const ShieldingSourceFitCalc::BackDecayResults * const bd = m_backDecay.get();
  const bool show = bd && bd->error.empty();
  m_backDecayFoot->setHidden( !show );
  toggleStyleClass( "HasNucFoot", show );
  if( !show )
    return;

  // The reference date, followed by how far it is from the measurement - or just the latter, if the
  //  spectrum has no start time.  (In date/time mode, the date shows how the entered text was read.)
  //  (Kept short, to fit on one line next to "Show Decay" at the tool's default size.)
  const double dt = bd->time_before_meas;
  const string when = ShieldingSourceFitCalc::back_decay_time_str( bd->ref_time );
  const string duration = ShieldingSourceDisplay::backDecayDurationText( fabs(dt) );
  if( when.empty() )
  {
    const WString span = WString::tr( (dt >= 0.0) ? "ssd-back-decay-before" : "ssd-back-decay-after" ).arg( duration );
    m_backDecayFootTxt->setText( WString::tr("sfn-foot-tref").arg( span ) );
    m_backDecayFootSpan->setText( WString() );
  }else
  {
    m_backDecayFootTxt->setText( WString::tr("sfn-foot-tref").arg( WString::fromUTF8(when) ) );
    m_backDecayFootSpan->setText( WString::tr( (dt >= 0.0) ? "sfn-foot-before" : "sfn-foot-after" ).arg( duration ) );
  }//if( when.empty() ) / else
  m_backDecayFootSpan->setHidden( when.empty() );

  m_showDecayBtn->setEnabled( !decayToolStart( *bd ).empty() );
}//void updateBackDecayFooter()


vector<SourceFitNuclideDisplay::DecayToolSource>
                SourceFitNuclideDisplay::decayToolStart( const ShieldingSourceFitCalc::BackDecayResults &bd )
{
  // Start at the reference time, so the decay runs through the measurement - or, for a reference
  //  time after the measurement, start at the measurement, so it runs through the reference time.
  const double dt = bd.time_before_meas;
  const bool from_ref = (dt >= 0.0);

  vector<DecayToolSource> answer;
  for( const ShieldingSourceFitCalc::BackDecaySourceResult &src : bd.sources )
  {
    DecayToolSource start;
    start.nuclide = src.nuclide;
    start.activity = from_ref ? src.activity : src.activity_now;
    start.age = from_ref ? src.age_at_ref : std::max( 0.0, src.age_at_ref + dt );

    if( start.nuclide && !start.nuclide->isStable()
       && (src.status == ShieldingSourceFitCalc::BackDecayStatus::Ok)
       && std::isfinite(start.activity) && (start.activity > 0.0) )
    {
      answer.push_back( start );
    }
  }//for( loop over sources )

  return answer;
}//decayToolStart(...)


void SourceFitNuclideDisplay::showDecayFromRefTime()
{
  // Use the current values, rather than those from the last refresh.
  const shared_ptr<const ShieldingSourceFitCalc::BackDecayResults> bd
                                              = m_backDecayProvider ? m_backDecayProvider() : nullptr;
  if( !bd || !bd->error.empty() || !m_interspec )
    return;

  const vector<DecayToolSource> sources = decayToolStart( *bd );
  if( sources.empty() )
    return;

  DecayWindow * const window = m_interspec->createDecayInfoWindow();
  if( !window )
    return;

  // The click that got us here raises the Activity/Shielding window to the front after the decay
  //  window is shown, so raise the decay window again on a delay (as LlmToolGui does).
  window->doJavaScript( "setTimeout(function(){Wt.WT.AuxWindowBringToFront('" + window->id() + "');},5);" );

  window->clearAllNuclides();

  const double span = 3.0 * fabs( bd->time_before_meas );
  const bool curies = m_model->displayCuries();
  for( const DecayToolSource &src : sources )
  {
    const SandiaDecay::Nuclide * const nuc = src.nuclide;
    window->addNuclide( nuc->atomicNumber, nuc->massNumber, nuc->isomerNumber, src.activity,
                        curies, src.age, "", ((span > 0.0) ? span : -1.0) );
  }//for( loop over sources )
}//void showDecayFromRefTime()


void SourceFitNuclideDisplay::handlePeakModelChanged()
{
  // A peak's use-in-fit (or its color/energy) changed; reconcile refreshes every card, which in
  //  nest-peaks mode rebuilds the per-card peak lists.  Guarded by m_refreshing so our own
  //  setData() during a toggle does not recurse.
  reconcileCards();
}


std::string SourceFitNuclideDisplay::peakColorForNuclide( const SandiaDecay::Nuclide *nuc ) const
{
  if( !m_peakModel || !nuc )
    return string();

  const int npeaks = static_cast<int>( m_peakModel->npeaks() );
  for( int row = 0; row < npeaks; ++row )
  {
    const PeakModel::PeakShrdPtr p = m_peakModel->peakPtr( row );
    if( p && (p->parentNuclide() == nuc) && !p->lineColor().isDefault() )
      return p->lineColor().cssText();
  }//for( loop over peaks )

  return string();
}//std::string peakColorForNuclide( const SandiaDecay::Nuclide *nuc ) const


std::vector<SourceFitNuclideDisplay::NucActivitySummary>
SourceFitNuclideDisplay::activitySummaries() const
{
  std::vector<NucActivitySummary> out;
  const int n = m_model->numNuclides();
  for( int i = 0; i < n; ++i )
  {
    const SandiaDecay::Nuclide * const nuc = m_model->nuclide( i );
    if( !nuc )
      continue;

    NucActivitySummary s;
    s.nuclide = nuc;
    s.symbol = nuc->symbol;
    s.colorCss = peakColorForNuclide( nuc );
    s.fitActivity = m_model->fitActivity( i );

    const double activity = m_model->activity( i );
    s.activity = PhysicalUnits::printToBestActivityUnits( activity, 2, m_model->displayCuries() )
                 + DetectorPeakResponse::det_eff_geom_type_postfix( m_model->detType() );

    const double actUncert = m_model->activityUncert( i );
    if( actUncert > 0.0 )
      s.activityUncert = PhysicalUnits::printToBestActivityUnits( actUncert, 2, m_model->displayCuries() );

    const BackDecayText bdtxt = back_decay_text( *m_model, i, m_backDecay.get() );
    if( bdtxt.show )
      s.backDecayActivity = bdtxt.ok ? bdtxt.value
                                     : WString::tr("sfn-act-tref-undetermined").toUTF8();

    if( !PeakDef::ageFitNotAllowed( nuc ) )
    {
      const double age = m_model->age( i );
      if( age > 0.0 )
        s.age = PhysicalUnitsLocalized::printToBestTimeUnits( age, 2, PhysicalUnits::second );
    }//if( ageable )

    out.push_back( std::move(s) );
  }//for( loop over model nuclides )

  return out;
}//std::vector<NucActivitySummary> activitySummaries() const


void SourceFitNuclideDisplay::reconcileCards()
{
  if( m_refreshing )
    return;

  RefreshGuard guard( m_refreshing );

  const int n = m_model->numNuclides();

  // 1) Ensure a card exists for every model nuclide.
  set<const SandiaDecay::Nuclide *> present;
  for( int i = 0; i < n; ++i )
  {
    const SandiaDecay::Nuclide * const nuc = m_model->nuclide( i );
    if( !nuc )
      continue;
    present.insert( nuc );
    if( m_cards.find(nuc) == m_cards.end() )
    {
      auto cardUniq = std::make_unique<NuclideCard>( this, nuc );
      m_cards[nuc] = cardUniq.get();
      m_cardsHolder->addWidget( std::move(cardUniq) );
    }
  }//for( loop over model nuclides )

  // 2) Remove cards whose nuclide is no longer in the model.
  for( auto it = m_cards.begin(); it != m_cards.end(); )
  {
    if( present.count(it->first) )
    {
      ++it;
    }else
    {
      m_cardsHolder->removeWidget( it->second );  //returned unique_ptr is destroyed here
      it = m_cards.erase( it );
    }
  }//for( loop over existing cards )

  // 3) Refresh each surviving card's contents, with the back-decay values for the model as it is now.
  m_backDecay.reset();
  try
  {
    if( m_backDecayProvider )
      m_backDecay = m_backDecayProvider();
  }catch( std::exception &e )
  {
    cerr << "SourceFitNuclideDisplay: error computing back-decay values: " << e.what() << endl;
    assert( 0 );  //compute_back_decay should not throw
  }//try / catch

  updateBackDecayFooter();

  for( int i = 0; i < n; ++i )
  {
    auto it = m_cards.find( m_model->nuclide(i) );
    if( it != m_cards.end() )
      it->second->refresh();
  }

  // 4) Reorder to match the model order (only if needed, to avoid needless re-renders).
  bool inOrder = (static_cast<int>(m_cardsHolder->count()) == n);
  for( int i = 0; inOrder && (i < n); ++i )
    inOrder = (m_cardsHolder->widget(i) == m_cards[m_model->nuclide(i)]);

  if( !inOrder )
  {
    vector<unique_ptr<WWidget>> ordered;
    for( int i = 0; i < n; ++i )
      ordered.push_back( m_cardsHolder->removeWidget( m_cards[m_model->nuclide(i)] ) );
    for( auto &w : ordered )
      m_cardsHolder->addWidget( std::move(w) );
  }//if( !inOrder )

  // Phone accordion cards default to expanded (m_open initialized to m_nestPeaks), so all peak
  //  rows are visible up front; the user can collapse individual cards and that state persists
  //  across refreshes (each card object outlives a reconcile).
}//void reconcileCards()
