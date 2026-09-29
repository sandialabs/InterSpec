#ifndef PeakFileIo_h
#define PeakFileIo_h
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

#include <set>
#include <map>
#include <deque>
#include <memory>
#include <string>
#include <vector>
#include <istream>
#include <ostream>

#include "SpecUtils/SpecFile.h"

class PeakDef;
namespace RefLib{ class Library; }

/** Reading and writing peaks in N42-2012 and IAEA SPE files, the way InterSpec does, so either
 program can read the other's files.  Only built with the LIGHT_PEAK_FILE_IO option.

 - N42-2012: InterSpec keeps the peaks of each set of displayed samples, and which samples were
   displayed, in a `<DHS:InterSpec>` element (see InterSpec's `SpecMeas::appendSpecMeasStuffToXml`).
 - IAEA SPE: InterSpec appends `$PEAKLABELS:` and `$PEAK_INFO_CSV:` sections (its peak CSV) to the
   displayed spectrum (see InterSpec's `SpecMeas::write_iaea_spe`); only the CSV is read back.
 */
namespace PeakFileIo
{
  typedef std::deque<std::shared_ptr<const PeakDef>> PeakDeque;
  typedef std::map<std::set<int>,PeakDeque> SampleNumsToPeaks;

  struct FilePeaks
  {
    /** Peaks, keyed by the sample numbers they were fit on. */
    SampleNumsToPeaks peaks;

    /** For N42 files: the samples displayed when the file was saved, and as what (e.g.,
     "Foreground"); empty if not given. */
    std::set<int> displayed_samples;
    std::string display_type;

    /** Problems that did not stop the peaks being read (e.g., a peak source that couldnt be read). */
    std::vector<std::string> warnings;
  };//struct FilePeaks

  /** Where a file holds peaks InterSpec saved: an N42 `<DHS:InterSpec>` element, or an SPE
   `$PEAK_INFO_CSV:` section. */
  struct SavedPeaksLocation
  {
    enum class Type{ None, N42, Spe };
    Type type = Type::None;
    size_t offset = 0;
  };//struct SavedPeaksLocation

  /** Scans the file for saved peak information (without parsing it). */
  SavedPeaksLocation locate_saved_peaks( const std::string &path );

  /** A SpecFile that, like InterSpec's SpecMeas, keeps the sample numbers of files InterSpec wrote,
   since their saved peaks are keyed by sample number (cleanup_after_load may otherwise renumber
   e.g., gapped sample numbers).
   */
  class SampleKeepingSpecFile : public SpecUtils::SpecFile
  {
  public:
    explicit SampleKeepingSpecFile( const bool keep ) : m_keep_sample_numbers( keep ) {}
    void cleanup_after_load( const unsigned int flags = StandardCleanup ) override;

  private:
    bool m_keep_sample_numbers;
  };//class SampleKeepingSpecFile

  /** Reads the peaks saved at `loc` in the file at `path` (already parsed as `spec`).  Peak sources
   in SPE files are matched to the reference-line library `lib`.  Throws if the peak information is
   invalid.
   */
  FilePeaks read_file_peaks( const std::string &path, const SavedPeaksLocation &loc,
                             const SpecUtils::SpecFile &spec, const RefLib::Library &lib );

  /** Writes `spec` as N42-2012, adding InterSpec's `<DHS:InterSpec>` element with `peaks`, and the
   samples and detectors displayed.
   */
  void write_n42( std::ostream &out, const SpecUtils::SpecFile &spec, const SampleNumsToPeaks &peaks,
                  const std::set<int> &displayed_samples,
                  const std::vector<std::string> &displayed_detectors );

  /** Returns the IAEA SPE file `spe` (as SpecUtils writes it) with InterSpec's peak sections added
   before its final `$ENDRECORD:`; `data` is the spectrum written in the file.
   */
  std::string add_spe_peaks( const std::string &spe, const PeakDeque &peaks,
                             const std::shared_ptr<const SpecUtils::Measurement> &data );

  /** Port of InterSpec's `PeakModel::csv_to_candidate_fit_peaks`: reads a peak CSV (as InterSpec,
   and `PeakCsv::write`, write them); sources are matched to the reference-line library `lib`.
   Throws on error, or if no peaks are found.
   */
  std::vector<PeakDef> read_peak_csv( std::istream &csv,
                                      const std::shared_ptr<const SpecUtils::Measurement> &meas,
                                      const RefLib::Library &lib );
}//namespace PeakFileIo

#endif //PeakFileIo_h
