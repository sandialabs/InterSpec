#ifndef DrfImport_h
#define DrfImport_h
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

#include <string>
#include <vector>
#include <memory>
#include <cstdint>

#include "InterSpec/DetectorEffG2kPar.h"

struct AngleOutxContents;
class DetectorPeakResponse;
class DetectorEfficiencyUncert;

/** Identifies and imports detector-efficiency (DRF) files, independent of any GUI.

 Used by `DrfImportWidget`, which is shown both on the "Import" tab of the DRF Select dialog and
 when a DRF file is dropped onto the app.  The flow is:
   1. `parseFile` - identify one file's kind and parse it (cheap; works on in-memory bytes, since
      uploaded spool files are transient).
   2. `makeSource` - combine a file with its companion, if any (Efficiency.csv + Detector.dat, or
      a .par grid + DETECTOR.txt), and work out which interpretations are allowed.  Slow only for
      the .par pair (see `isSlow`).
   3. `build` - apply the user's interpretation/options; always returns a NEW DRF, and never
      modifies the `Source` or any DRF it returned before.
 */
namespace DrfImport
{
  /** Largest file accepted; far above any real DRF file. */
  const size_t sm_maxFileBytes = 32*1024*1024;

  enum class FileKind
  {
    DrfXml,              ///< InterSpec `<DetectorPeakResponse>` XML
    MakeDrfCsv,          ///< CSV exported by the "Make Detector Response" tool
    MultiDrfCsv,         ///< one-DRF-per-line relative-efficiency CSV/TSV
    GammaQuantCsv,       ///< GammaQuant column-format efficiency CSV
    IsocsEcc,            ///< ISOCS .ecc
    Angle,               ///< ANGLE XML (.outx, .detx, ...)
    EfficiencyCsv,       ///< energy/efficiency CSV (gamEff, Run_effoutput, ...)
    GadrasEfficiencyCsv, ///< GADRAS Efficiency.csv
    GadrasDetectorDat,   ///< GADRAS Detector.dat
    ParGrid,             ///< binary .par full-energy-peak efficiency grid
    ParDetectorTxt       ///< DETECTOR.txt geometry records for .par grids
  };//enum class FileKind

  /** Kinds whose header matches, most likely first, from the start of a file.  Only nominates;
   `parseFile` is the authority.
   */
  std::vector<FileKind> candidateKinds( const uint8_t *header, size_t headerLen, size_t fileSize );

  /** Whether a kind is a complete DRF (or list of them) that needs no interpretation. */
  bool isCompleteDrfKind( const FileKind kind );

  /** Whether the two kinds are the two halves of one detector definition. */
  bool areCompanions( const FileKind a, const FileKind b );

  /** Reads a whole file (UTF-8 path; handles Windows).  Throws if it can't, or it is over
   `sm_maxFileBytes`.
   */
  std::shared_ptr<const std::string> readFile( const std::string &path );


  /** One identified and parsed file. */
  struct ParsedFile
  {
    std::string displayName;
    FileKind kind = FileKind::DrfXml;
    std::shared_ptr<const std::string> data;

    /** The DRF(s) the file holds: complete DRFs, the efficiency curve to interpret, or (for a
     Detector.dat) a geometry-only seed.  Empty for .par/DETECTOR.txt.
     */
    std::vector<std::shared_ptr<const DetectorPeakResponse>> drfs;

    /** ISOCS .ecc source quantities (PhysicalUnits; zero if absent) and per-energy uncertainties. */
    double sourceArea = 0.0, sourceMass = 0.0;
    std::vector<float> uncertEnergies, baselineFrac, convergenceFrac;

    std::shared_ptr<const AngleOutxContents> angle;
    std::shared_ptr<const DetEffG2kPar::ParFile> par;
    std::vector<DetEffG2kPar::DetectorDef> detectorDefs;

    /** Warnings and parse notes from the file (plain text, English as the parsers give them). */
    std::vector<std::string> notes;

    /** The "#credit:" lines of a list of DRFs, which are HTML. */
    std::vector<std::string> credits;
  };//struct ParsedFile

  /** Identifies and parses a file.  With `anyTextAsCsv`, a file matching no header pattern is
   still tried as an energy/efficiency CSV, and as a list of DRFs (the DRF Select "Import" tab is
   only ever given DRF files; an app-wide drop is not).  Throws std::exception if it is not a DRF
   file.
   */
  std::shared_ptr<const ParsedFile> parseFile( const std::string &displayName,
                                               std::shared_ptr<const std::string> data,
                                               const bool anyTextAsCsv );


  enum class Interpretation
  {
    AsIs,               ///< use the file's DRF unchanged
    FarFieldIntrinsic,  ///< efficiencies are intrinsic; needs a diameter
    FarFieldAbsolute,   ///< efficiencies are absolute at a distance; needs diameter + distance
    FixedTotal,         ///< fixed geometry, total activity
    FixedPerCm2,
    FixedPerM2,
    FixedPerGram,
    GenericDetector     ///< ANGLE detector model + reference curve, as a geometry-modeled DRF
  };//enum class Interpretation

  enum class Status
  {
    Ready,
    NeedsCompanion,         ///< a .par without its DETECTOR.txt, or vice-versa
    NeedsCharacterization,  ///< geometry only (Detector.dat, ANGLE .detx); see `Source::base`
    NeedsDiameter,
    NeedsDistance,
    Error                   ///< see `Result::error`
  };//enum class Status

  /** A file, with its companion if any, ready to be built into a DRF. */
  struct Source
  {
    Status status = Status::Error;

    /** The DRF to interpret or copy; for `NeedsCharacterization`, the geometry-only seed. */
    std::shared_ptr<const DetectorPeakResponse> base;

    /** For `Interpretation::GenericDetector`: the ANGLE reference curve with the detector
     geometry and a curve-transfer response attached.
     */
    std::shared_ptr<const DetectorPeakResponse> generic;

    /** Allowed interpretations; empty means `base` is used as-is. */
    std::vector<Interpretation> interpretations;
    Interpretation defaultInterpretation = Interpretation::AsIs;

    double sourceArea = 0.0, sourceMass = 0.0;
    std::vector<float> uncertEnergies, baselineFrac, convergenceFrac;

    /** Choices when a file holds several DRFs, or a DETECTOR.txt several records none of which
     name the .par; `record` is the one used.
     */
    std::vector<std::string> recordNames;
    size_t record = 0;

    /** Suggested DRF name; `nameIsFileStem` when it came from the file name, not the file. */
    std::string name;
    bool nameIsFileStem = false;

    std::vector<std::string> notes;    ///< plain text
    std::vector<std::string> credits;  ///< HTML
    std::string error;
  };//struct Source

  /** Whether `makeSource` for these files takes long enough (seconds) to belong off the session
   thread - i.e., building from a .par grid.
   */
  bool isSlow( const ParsedFile *a, const ParsedFile *b );

  /** Combines `a` with its companion `b` (either may be the companion; `b` may be null).
   `record` selects among `Source::recordNames`.  Never throws; failures give `Status::Error`.
   */
  Source makeSource( std::shared_ptr<const ParsedFile> a,
                     std::shared_ptr<const ParsedFile> b,
                     const size_t record );

  struct Options
  {
    Interpretation interpretation = Interpretation::AsIs;
    double diameter = 0.0;   ///< PhysicalUnits
    double distance = 0.0;   ///< PhysicalUnits
    double setback = -1.0;   ///< PhysicalUnits; negative keeps the file's
    bool setUncert = false;  ///< replace the efficiency uncertainty with `uncert` (null clears)
    std::shared_ptr<const DetectorEfficiencyUncert> uncert;
    std::string name;        ///< empty keeps the DRF's name
  };//struct Options

  struct Result
  {
    Status status = Status::Error;
    std::shared_ptr<DetectorPeakResponse> drf;  ///< a new object; set only when `Ready`
    std::string error;
  };//struct Result

  Result build( const Source &src, const Options &options );
}//namespace DrfImport

#endif //DrfImport_h
