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

#include <Wt/WAny.h>
#include <Wt/WText.h>
#include <Wt/WLineEdit.h>
#include <Wt/WCheckBox.h>
#include <Wt/WModelIndex.h>
#include <Wt/WContainerWidget.h>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
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
}//namespace


/** One per-nuclide card.  Reads/writes everything through #SourceFitModel. */
class SourceFitNuclideDisplay::NuclideCard : public Wt::WContainerWidget
{
public:
  NuclideCard( SourceFitNuclideDisplay *disp, const SandiaDecay::Nuclide *nuc )
   : WContainerWidget(),
     m_disp( disp ),
     m_model( disp->m_model ),
     m_nuc( nuc )
  {
    addStyleClass( "NucCard" );

    // ---- header: name, half-life, source-type tag ----
    WContainerWidget *head = addNew<WContainerWidget>();
    head->addStyleClass( "NucHead" );
    m_nameLabel = head->addNew<WText>();
    m_nameLabel->addStyleClass( "NucName" );
    m_halfLifeLabel = head->addNew<WText>();
    m_halfLifeLabel->addStyleClass( "NucHl" );
    m_tag = head->addNew<WText>();
    m_tag->addStyleClass( "NucTag" );
    m_tag->hide();

    // ---- activity row ----
    WContainerWidget *actRow = addNew<WContainerWidget>();
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

    // ---- mass + act(T0) line ----
    WContainerWidget *massLine = addNew<WContainerWidget>();
    massLine->addStyleClass( "NucMass" );
    m_massLabel = massLine->addNew<WText>();
    m_actT0Label = massLine->addNew<WText>();
    m_actT0Label->addStyleClass( "NucT0" );
    m_actT0Label->hide();

    // ---- age row ----
    m_ageRow = addNew<WContainerWidget>();
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
    m_progenyRow = addNew<WContainerWidget>();
    m_progenyRow->addStyleClass( "NucProgeny" );
    m_progenyRow->hide();

    refresh();
  }//NuclideCard constructor


  const SandiaDecay::Nuclide *nuclide() const { return m_nuc; }


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

    // Mass + act(T0)
    const bool fitAct = m_model->fitActivity( row );
    m_massLabel->setText( WString("{1} {2}")
         .arg( WString::tr( fitAct ? "sfn-mass-approx" : "sfn-mass" ) )
         .arg( Wt::asString( m_model->data( m_model->index(m_nuc, SourceFitModel::kIsotopeMass) ) ) ) );

    const bool ageable = !PeakDef::ageFitNotAllowed( m_nuc );

    if( ageable )
    {
      const double activity = m_model->activity( row );
      const double age = m_model->age( row );
      const double hl = m_nuc->halfLife;
      // Activity at T0 (production) that decays to the current activity at the given age.
      const double actT0 = (hl > 0.0) ? (activity * std::exp( age * 0.6931471805599453 / hl )) : activity;
      string t0str = PhysicalUnits::printToBestActivityUnits( actT0, 3, m_model->displayCuries() );
      t0str += DetectorPeakResponse::det_eff_geom_type_postfix( m_model->detType() );
      m_actT0Label->setText( WString("{1}={2}").arg( WString::tr("sfn-act-t0") ).arg( WString::fromUTF8(t0str) ) );
      m_actT0Label->setHidden( age <= 0.0 );
    }else
    {
      m_actT0Label->hide();
    }

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
  }//void refresh()


protected:
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

  Wt::WText *m_nameLabel;
  Wt::WText *m_halfLifeLabel;
  Wt::WText *m_tag;
  Wt::WLineEdit *m_activityEdit;
  Wt::WCheckBox *m_fitActivityCB;
  Wt::WText *m_massLabel;
  Wt::WText *m_actT0Label;
  Wt::WContainerWidget *m_ageRow;
  Wt::WLineEdit *m_ageEdit;
  Wt::WCheckBox *m_fitAgeCB;
  Wt::WContainerWidget *m_progenyRow;
};//class SourceFitNuclideDisplay::NuclideCard



SourceFitNuclideDisplay::SourceFitNuclideDisplay( SourceFitModel *model,
                                                  std::shared_ptr<PeakModel> peakModel,
                                                  InterSpec *interspec )
  : WContainerWidget(),
    m_model( model ),
    m_peakModel( peakModel ),
    m_interspec( interspec ),
    m_cardsHolder( nullptr ),
    m_refreshing( false )
{
  addStyleClass( "SourceNuclideDisplay" );

  m_cardsHolder = addNew<WContainerWidget>();
  m_cardsHolder->addStyleClass( "NucCardsHolder" );

  // Observe the model; SourceFitModel stays the single source of truth.
  model->dataChanged().connect( this, &SourceFitNuclideDisplay::handleModelChanged );
  model->rowsInserted().connect( this, &SourceFitNuclideDisplay::handleModelChanged );
  model->rowsRemoved().connect( this, &SourceFitNuclideDisplay::handleModelChanged );
  model->layoutChanged().connect( this, &SourceFitNuclideDisplay::handleModelChanged );
  model->modelReset().connect( this, &SourceFitNuclideDisplay::handleModelChanged );

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

  // 3) Refresh each surviving card's contents.
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
}//void reconcileCards()
