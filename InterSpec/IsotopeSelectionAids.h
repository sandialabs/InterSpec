#ifndef IsotopeSelectionAids_h
#define IsotopeSelectionAids_h
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

#include <vector>
#include <string>

#include <boost/any.hpp>

#include <Wt/WModelIndex>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>
#include <Wt/WAbstractItemDelegate>

#include "InterSpec/PeakDef.h"

namespace Wt
{
  class WObject;
  class WWidget;
  class WSuggestionPopup;
  class WAbstractItemModel;
}//namespace Wt



namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
  struct EnergyRatePair;
  class SandiaDecayDataBase;
}//namespace SandiaDecay


//  - Another issue is right now it is not immediately obvious to the user that
//    the isotope name field is editable.
//  - Client side highlighting of matching text does not properly handle
//    meta-stable states (e.g. <b>Co</b><b>60</b>m when user entered co60m)


namespace IsotopeSelectionAids
{
//If nuclide can obtain secular equilibrium, returns that half life, else if
//  prompt equilibrium that, and if not, the normal half life.
double halfLife( const SandiaDecay::Nuclide *nuclide );


struct NearGammaInfo
{
  /*! Abs value of distance from energy (subtracts 511/1022 for S.E./D.E. before 
      computing). 
   */
  float distance;
    
  /*! Photopeak energy (will not have 511/1022 subtracted for S.E./D.E. peaks.*/
  float gamma_energy;
    
  /*! Intensity relative to largest photopeak (parent intensity for S.E./D.E., 
      not detectors)
   */
  float relative_intensity;
    
  /*! source gamma type */
  PeakDef::SourceGammaType gamma_type;
};//struct NearGammaInfo
  
std::vector< NearGammaInfo >
  equilibriumGammasByNearestEnergy( const SandiaDecay::Nuclide *nuclide,
                                    const double energy,
                                    const double minRelativeBr,
                                    const bool includeEscapePeaks,
                                    const bool includeXrays );
}//namespace IsotopeSelectionAids


class PhotopeakDelegate : public Wt::WAbstractItemDelegate
{
public:
  enum DelegateType
  {
    NuclideDelegate, GammaEnergyDelegate
  };

  class EditWidget : public Wt::WContainerWidget
  {
  public:
    EditWidget( const Wt::WModelIndex& index,
                const Wt::WFlags<Wt::ViewItemRenderFlag> flags,
                const bool closeOnBlur,
                const DelegateType delegateType,
                PhotopeakDelegate *parent );
    virtual ~EditWidget();

    Wt::WLineEdit *edit();

    void handleBlur();
    void handleBlurWorker( boost::function<void()> worker );
    
    static void replacerJs( std::string &js );
    static void nuclideNameMatcherJs( std::string &js );
    static void gammaEnergyMatcherJs( std::string &js );

  protected:
    PhotopeakDelegate *m_parent;
    Wt::WLineEdit *m_edit;
    Wt::WSuggestionPopup *m_suggestions;
  };//class EditWidget

public:
  PhotopeakDelegate( DelegateType delegateType, bool closeOnBlur, Wt::WObject *parent = NULL );
  virtual ~PhotopeakDelegate();
  virtual Wt::WWidget *update( Wt::WWidget *widget,
                               const Wt::WModelIndex &index,
                               Wt::WFlags< Wt::ViewItemRenderFlag > flags );

protected:

  //If the editor is closed because of blurring (another field gets selcted or
  //  something), than we dont want to save the edit results if the field is
  //  blank (this can happen if the user just clicks on the field to edit it,
  //  but doesnt actually type anything).  If the field is blank and the user
  //  presses enter, than we do want to save the result
  virtual void doCloseEditor( Wt::WWidget *editor, bool save, bool isBlurr ) const;

  boost::any editState( Wt::WWidget *editor ) const;
  void setEditState( Wt::WWidget *editor, const boost::any &value ) const;
  void setModelData( const boost::any &editState,
                     Wt::WAbstractItemModel *model,
                     const Wt::WModelIndex &index ) const;

  bool m_closeOnBlur;
  DelegateType m_delegateType;
  Wt::WSuggestionPopup *m_suggestionPopup;
};//class PhotopeakDelegate



class PeakIsotopeNameFilterModel : public  Wt::WAbstractItemModel
{
public:
  PeakIsotopeNameFilterModel( const std::shared_ptr<const PeakDef> &peak,
                    Wt::WObject *parent = 0 );
  virtual ~PeakIsotopeNameFilterModel();

  virtual Wt::WModelIndex index( int row, int column,
                                 const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;

  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;

  virtual int rowCount( const Wt::WModelIndex & parent = Wt::WModelIndex() ) const;
  virtual int columnCount( const Wt::WModelIndex & parent = Wt::WModelIndex() ) const;

  virtual boost::any data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const;

  
  static Wt::WString displayText( const SandiaDecay::Nuclide *txt,
                                  double mean, PeakDef::SourceGammaType type,
                                  int role );
  static Wt::WString displayText( const SandiaDecay::Element *el,
                                  int role, double minHalfLife );
  
  //determineAndRemoveIsoLevel(...) looks for patterns such as 'Co60m', 'Co60meta',
  //  'Co 60 m', etc. to determine if the user is inputting a metts stable state.
  //  The function returns the iso level (right no just 0, 1, or 2), and
  //  removes the portion of the text indicating the meta level, from the input.
  static int determineAndRemoveIsoLevel( std::string &label );

  //All returned alpha strings will be lowercase
  //  See IsotopeNameFilterModel for a different implementation
  static void getAlphaAndNumericSubStrs( std::string label,
                                         std::vector<std::string> &alphastrs,
                                         std::vector<std::string> &numericstrs );

  //TODO: filter(...) is very poorly coded (it's a mess!) and probably quite
  //      inefficient - this should be cleaned up at some point.
  //      On my MacbookPro 2.3 GHz Intel Core i7 execution time was always
  //      less than 0.001 seconds wall time (debug compilation) on 20120607.
  //XXX - only takes into account PeakDef::CandidateNuclide::nuclide,
  //      and not PeakDef::CandidateNuclide::transition->parent
  void filter( const Wt::WString &text );
//  Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex & index ) const;

  //setPeak(...): sets m_peak, and calls filter( m_filter ) - untested as of 20130509
  void setPeak( const std::shared_ptr<const PeakDef> &peak );

protected:
  const double m_minHalfLife; //seconds
  Wt::WString m_filter;
  std::shared_ptr<const PeakDef> m_peak;
  std::vector<Wt::WString> m_displayData, m_userData;
};//class PeakIsotopeNameFilterModel


#endif  //#ifndef( IsotopeSelectionAids_h )
