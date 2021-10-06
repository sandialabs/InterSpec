#ifndef PeakModel_h
#define PeakModel_h
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

#include <deque>
#include <vector>
#include <memory>

#include <Wt/WResource>
#include <Wt/WModelIndex>
#include <Wt/Chart/WDataSeries>
#include <Wt/WAbstractItemModel>

#include <boost/any.hpp>

#include "InterSpec/PeakDef.h"

class SpecMeas;
class SpectrumDataModel;

namespace SpecUtils{ class Measurement; }

namespace SandiaDecay
{
  struct Nuclide;
  struct Transition;
  struct EnergyIntensityPair;
}//namespace SandiaDecay


class PeakModel: public Wt::WAbstractItemModel
{
  /****************************************************************************\
  | This class is used for finding, storing, and analyzing PeakDef peaks.
  | It can locate, flag, unflag, etc. It also has support for using a
  | nuclide database to figure out what each peak is likely to be signifying.
  |
  | Each 'row' is a peak whose parent index is awlsy the default index,
  | Wt::WModelIndex(), and columns represent the various peak values
  | specified by the 'Columns' enum.
  \****************************************************************************/

public:
  enum Columns
  {
    kIsotope,           //The isotope to display as the isotope responsible for this gamma line
    kUseForCalibration,
    kMean, kFwhm, kAmplitude,
    kCps,
    kPhotoPeakEnergy,
    kDifference,  //calculate the difference
    kUseForShieldingSourceFit,
    kCandidateIsotopes,
    kPeakLineColor,
    kUserLabel,
    kHasSkew, kSkewAmount,
    kType, kLowerX, kUpperX, kRoiCounts, kContinuumType,
    kNumColumns
  };//enum Columns
  
  typedef std::shared_ptr<const PeakDef> PeakShrdPtr;

public:
  PeakModel( Wt::WObject *parent = 0 );
  virtual ~PeakModel();

  //Inorder to display continuum area of peaks, you need to set the data model
  // XXX - note that with PeakDef refactorization that happened in Dec 2013,
  //       this requirement/paradyn could be eliminated probably
  void setDataModel( SpectrumDataModel *dataModel );

  //setPeakFromSpecMeas(...): when the primary spectrum is changed (either file,
  //  or the displayed sample numbers), then this function should be called to
  //  update the m_peaks pointer to be from the appropriate SpecMeas object.
  void setPeakFromSpecMeas( std::shared_ptr<SpecMeas> meas,
                            const std::set<int> &samplenums );

  //recomendUseForFit(...): determines whether or not the specified photopeak
  //  should be used for activity/shielding fit.
  static bool recomendUseForFit( const SandiaDecay::Nuclide *n,
                                 const float energy );

  //Functions add/delete/access peaks from - from these functions peaks will
  //  always be sorted by mean, in an increasing fashion
  size_t npeaks() const;
  const PeakDef &peak( const size_t peakn ) const;  //throws if invalid peak numbers
  PeakModel::PeakShrdPtr peakPtr( const size_t peakn ) const;  //throws if invalid peak numbers
  
  //nearestPeak(...): returns the nearest peak in energy.  If no peaks will
  //  return a null pointer
  PeakModel::PeakShrdPtr nearestPeak( double energy ) const;
  
  //indexOfPeak(): returns index of peak passed in; an invalid index if peak
  //  pointed to by 'peak' isnt in this model.  Note that comparison is done
  //  using the PeakDef memory location, not the PeakDef attributes.
  Wt::WModelIndex indexOfPeak( const PeakShrdPtr &peak ) const;
  
  std::shared_ptr<const std::deque< PeakModel::PeakShrdPtr > > peaks() const;
  std::vector<PeakDef> peakVec() const;

  //definePeakXRangeAndChi2(...): Inorder to save cpu (and mostly memorry access
  //  time) later on, this function will define the lower and upper energy range
  //  of peaks that have a polynomial continum defined, and have not already had
  //  the lower and upper energy values defined.  Chi2 is always computed. This
  //  function is called for all peaks in addPeaks, addPeak, and setPeaks.
  void definePeakXRange( PeakDef &peak );
  
  //addNewPeak(...): adds peak to the model, and returns index of added peak.
  //  Note that if you want the added peak to be assigned to be from the
  //  currently showing reference gamma lines, you should call
  //  InterSpec::addPeak(...)
  Wt::WModelIndex addNewPeak( const PeakDef &peak );
  
  //addPeaks(...): adds multiple peaks to the model.
  //  Does not assign currently showing reference gamma lines to peaks
  void addPeaks( const std::vector<PeakDef> &peaks );
  
  //setPeaks(...): removes all old peaks and adds peaks passed in
  void setPeaks( std::vector<PeakDef> peaks );

  //setPeaks(...): a convience function that calls the setPeaks(vector<PeakDef>)
  //  function.
  //  Note, the PeakDef objects pointed to in 'peaks' are not what are added
  //  to the model, but rather copies of them.
  void setPeaks( const std::vector<std::shared_ptr<const PeakDef> > &peaks );
  
  //removePeak( size_t ):  removes the specified `peakn`, where the peakn refers
  //  to the peak at position `peakn` when sorted by peak means (similar
  //  behavior as all other functions which refer to peaks by index).
  //  Throws if invalid peak numbers.
  void removePeak( const size_t peakn );
  
  //removePeak( WModelIndex ): removed peak at specified index.
  //  Throws if invalid peak index.
  void removePeak( Wt::WModelIndex index );
  
  //removePeak( PeakShrdPtr ): remove peak passed in
  //  Throws if the model does not contain a pointer == to the one passed in.
  void removePeak( PeakModel::PeakShrdPtr peak );
  
  /** Removes all the passed in peaks.
   
   Prefer this function call over #removePeak whenever there is more than one peak, to allow (future) undo/redo implementation.
   
   Will throw runtime exception if and of the peaks are not owned by this model.
   */
  void removePeaks( const std::vector<PeakModel::PeakShrdPtr> &peak );
  
  //removeAllPeaks(): removes all peaks
  void removeAllPeaks();
  
  
  /** Returns all peaks owned by this model that share a continuum with the peak passed in, including the peak passed in.
   
   Throws exception if peak passed in is non-null and not owned by this model.
   */
  std::vector<std::shared_ptr<const PeakDef>> peaksSharingRoi( const std::shared_ptr<const PeakDef> &peak );
  
  /** Returns all peaks owned by this model that DO NOT share a continuum with the peak passed in.
   
   Throws exception if peak passed in is non-null and not owned by this model.
   */
  std::vector<std::shared_ptr<const PeakDef>> peaksNotSharingRoi( const std::shared_ptr<const PeakDef> &peak );
  
  /** Removes the passed in 'originalPeak', then creates a new peak with values of 'newPeak'.
   
   Causes the removeRow() followed by the rowsInserted() signals to be emitted.
   
   Throws exception if peak passed in is non-null and not owned by this model.
   */
  void updatePeak( const std::shared_ptr<const PeakDef> &originalPeak, const PeakDef &newPeak );

  /** Similar to #updatePeak, but for multiple peaks.
   
   Original and new peaks do not have to be the same size.
   Will potentially cause multiple removeRow() and rowsInserted() signals to be emitted.
   
   Throws exception if and original peaks passed are not owned by this model.
   */
  void updatePeaks( const std::vector<std::shared_ptr<const PeakDef>> &originalPeaks,
                   const std::vector<PeakDef> &newPeaks );
  
  
  //setPeakFitFor(...): sets wether the specified coefficient should be fit for.
  //  The original peak is released from memory (assuming nowhere else has a
  //  shared pointer to it), and replaced with a new one, but with the
  //  m_fitFor[coef] value changed.  No changed/removed/added signals are
  //  emmitted, and nothing else is done.
  //  Note that this function is kinda against the philosphy of this class,
  //  but is used as an optimization to keep from updating the chart and every
  //  where else when the fitFor value is changed since it is only changed from
  //  the PeakEdit widget, and also only used in the PeakEdit widget.  If
  //  anything ever will be effected by changing the fitFor value, then this
  //  approach should be changed/reconsidered.
  //Throws runtime_error when an invalid index is passed in, or there is no
  //  spectrum loaded
  void setPeakFitFor( const Wt::WModelIndex index,
                      const PeakDef::CoefficientType coef,
                      const bool fitfor );
  
  //setContinuumPolynomialFitFor(...): analagous to setPeakFitFor(...), but for
  //  the continuum polynomial coefficients.
  void setContinuumPolynomialFitFor( const Wt::WModelIndex index,
                                     size_t polyCoefNum,
                                     const bool fitfor );


  //Functions for the Wt::WAbstractItemModel interface - from these functions
  //  peaks will be sorted according to m_sortColumn and m_sortOrder
  const PeakShrdPtr &peak( const Wt::WModelIndex &index ) const;
  
  virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;
  virtual boost::any data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const;
  virtual Wt::WModelIndex index( int row, int column, const
                                 Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual boost::any headerData( int section,
                                 Wt::Orientation orientation = Wt::Horizontal,
                                 int role = Wt::DisplayRole ) const;

  virtual bool removeRows( int row, int last,
                             const Wt::WModelIndex &parent = Wt::WModelIndex() );
  //removeColumns(...) should not be used and will throw std::runtime_error
  //  if called - after testing, can remove this function
  virtual bool removeColumns( int row, int count,
                              const Wt::WModelIndex &parent = Wt::WModelIndex() );

  virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const;
  virtual bool setData( const Wt::WModelIndex &index, const boost::any &value,
                        int role = Wt::EditRole );

  //setNuclideXrayReaction(): sets the nuclide, xray, or reaction specified
  //  by the 'txt' string.  If the string specifies an energy, that is used
  //  to help decide which gamma of the nuclide/xray/reaction should be
  //  assigned, otherwise the peak mean is used.
  //  'nsigma_window' is only applicable for the nuclide case (so not xray or
  //  reaction), and provides a notion of tradoff between ditance in energy, and
  //  intensity of gamma assigned; if <= 0.0, then only distance in energy will
  //  be used to decide which gamma is assigned, other wise:
  //  (0.1*nsigma_window*peak.sigma() + delta_energy)/gamma_relative_intensity;
  //  will used as; a value of 4.0 for nsigma_window yields good results for
  //  when the user clicks on a on a reasonalaby isolated peak, but can fail
  //  when identiying a small intensity photopeak near a much larger intenisty
  //  peak thats already been identified.
  enum SetGammaSource{ NoSourceChange, SourceChange, SourceAndUseChanged };
  static SetGammaSource setNuclideXrayReaction( PeakDef &peak, std::string txt,
                                                const double nsigma_window );
  

  virtual void sort( int column, Wt::SortOrder order = Wt::AscendingOrder );
  static bool compare( const PeakShrdPtr &lhs, const PeakShrdPtr &rhs,
                       Columns column, Wt::SortOrder order,
                       const std::shared_ptr<const SpecUtils::Measurement> &data );

  //isWithinRange(...) makes sure the peak mean is inside of the x-axis range
  bool isWithinRange( const PeakDef &peak ) const;
  bool isOutOfRange( const PeakDef &peak ) const;  //conveinience function

  Wt::WResource *peakCsvResource();

  //notifySpecMeasOfPeakChange(): The SpecMeas which the m_peaks actually
  //  belongs to needs to be notified when peaks are added/removed/modified, so
  //  it can mark itself as modified; this function does that notification.
  void notifySpecMeasOfPeakChange();
  
  /** Reads a CSV file generated by PeakCsvResource to peaks that are suitable
   as starting values to fit those peaks.
   
   Currently slightly forgiving that the CSV doesnt exactly match the output
   from this class, but flexible input is not a goal of this function (at this
   time - maybe in the future).
   
   Throws exception if inputs are invalid, or no peaks present in the file.
   */
  static std::vector<PeakDef> csv_to_candidate_fit_peaks(
                                                         std::shared_ptr<const SpecUtils::Measurement> meas,
                                                         std::istream &csv );
  
  
  /** Writes the peaks to a CSV file output - e.g., what the user gets when they click on the CSV download on the "Peak Manager" tab.
   */
  static void write_peak_csv( std::ostream &outstrm,
                             std::string specfilename,
                             const std::deque<std::shared_ptr<const PeakDef>> &peaks,
                             const std::shared_ptr<const SpecUtils::Measurement> &data );
  
  
protected:
  
  /** Adds a new peak to the model, returning the inserted peak and its index.
   
   Causes the rowsInserted() signal to be emitted.
   
   In the future this function will not modify undo/redo history.
   */
  std::pair<std::shared_ptr<const PeakDef>,Wt::WModelIndex> addNewPeakInternal( const PeakDef &peak );
  
  /** Removes a peak from the model.
   
   Causes the rowsRemoved() signal to be emitted.
   
   Throws exception if passed in peak is not owned by this model.
   
   In the future this function will not modify undo/redo history.
   */
  void removePeakInternal( std::shared_ptr<const PeakDef> peak );
  
  
  SpectrumDataModel *m_dataModel;

  //m_peaks and m_sortedPeaks contain the same peaks, they just differ in how
  //  they are sorted.

  //Peaks are stored sorted according to peak.mean, since some of the algorithms
  //  that work on the peaks need them to be sorted by mean
  std::shared_ptr< std::deque< PeakShrdPtr > > m_peaks;
  std::weak_ptr<SpecMeas> m_measurment;
  
  Columns m_sortColumn;
  Wt::SortOrder m_sortOrder;
  //peaks in m_sortedPeaks are sorted according to m_sortOrder and m_sortColumn
  //  for user visualization purposes
  std::deque< PeakShrdPtr > m_sortedPeaks;

  class PeakCsvResource;
  PeakCsvResource *m_csvResource;

  class PeakCsvResource : public Wt::WResource
  {
    //Simple class to allow downloading of peak info as a CSV file - compatible
    //  with other major peak analysis tool pair
  public:
    PeakCsvResource( PeakModel *parent );
    virtual ~PeakCsvResource();

  private:
    PeakModel *m_model;
    Wt::WApplication *m_app;
    
    
    virtual void handleRequest( const Wt::Http::Request &request,
                                Wt::Http::Response &response );
  };//class PeakCsvResource : public Wt::WResource

  
  friend class PeakCsvResource;
};//class PeakModel

#endif //PeakModel_h
