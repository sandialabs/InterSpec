#ifndef GadrasDetectorDat_h
#define GadrasDetectorDat_h
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

#include <map>
#include <array>
#include <string>
#include <cstdint>
#include <iosfwd>
#include <vector>

/** \brief In-memory representation of a GADRAS `Detector.dat` file.

 GADRAS describes a gamma detector's response with a `Detector.dat` file that comes in two
 on-disk variants:
   - a legacy fixed-column TEXT format: a flat list of numbered parameter lines,
     `<index> <value> <fit_flag>   <label text>`, for indices 1..84 (older files stop at 64
     or 80).  A single line `Version <x.y.z>` is wedged in between params 64 and 65 of the
     newer text files.
   - an XML format: `<gamma_detector>` with human-readable nested sub-blocks (`<dimensions>`,
     `<peak_shape>`, `<lower_level_discrimination>`, attenuators, shields, ...).  Every numeric
     leaf is a `<value>` element with an optional `<varying>` (the fit flag).

 This struct reads either variant into a single representation and can write the XML variant
 back out.  The 84 GADRAS parameters are the single source of truth: they are stored verbatim
 in #m_params (both the value column and the integer/fit-flag column), so a parsed file can be
 losslessly round-tripped and only the fields a caller owns need be overwritten.  Named
 accessors read #m_params by index and return values in GADRAS native units (cm, g/cm2, keV,
 %, dimensionless).

 The parameter meanings, defaults, and the material table below are transcribed from
 `scratch/GADRAS_DRF_parameter_documentation.md`; section references (e.g. "doc §5") point back
 into that document.

 \note The legacy text format carries the crystal material only as an integer INDEX in the
       integer (fit-flag) column of param 59, which indexes into the default GADRAS material
       array (the #materialTable below).  The text file has no material name.  The XML format
       instead stores a material NAME (`<material>NaI</material>`).  Hence #materialName resolves
       an index through #materialByIndex for text files, and uses the stored name for XML files.
 */
struct GadrasDetectorDat
{
  /** Which on-disk variant this instance was parsed from. */
  enum class SourceFormat
  {
    LegacyText,   //!< Fixed-column numbered-parameter text file (doc §1.1).
    Xml,          //!< `<gamma_detector>` XML file (doc §1.2).
    Unknown       //!< Default-constructed, not yet parsed.
  };

  /** One GADRAS parameter line.

   In the text format @p value is the second (value) column and @p int_col is the third
   (fit-flag / integer) column.  Several parameters overload the integer column to carry an
   actual integer datum (e.g. the material index at param 59, min/max energy at 62/63, channel
   count at 64) - see doc §2.  In the XML format @p int_col holds the `<varying>` flag (0/1)
   for most leaves, or the unpacked integer datum for the overloaded parameters. */
  struct Param
  {
    float value = 0.0f;   //!< Value column (native units).
    int   int_col = 0;    //!< Integer / fit-flag column.
  };

  /** Number of GADRAS parameters (doc: `NUMBER_PARAMETERS = 84`). */
  static constexpr int sm_num_params = 84;

  SourceFormat m_format = SourceFormat::Unknown;

  /** The 84 parameters, 1-based (index 0 is unused so param N lives at m_params[N]). */
  std::array<Param, sm_num_params + 1> m_params{};

  /** Highest parameter index actually present in the source file (for round-trip: older files
   stop at 64 or 80).  0 until a file is parsed. */
  int m_highest_param = 0;

  /** The GADRAS response-model version, e.g. "19.2.3".  From the text `Version <x>` line, or the
   XML `<response_version>` element.  Empty implies an old (<=64 param) legacy file. */
  std::string m_response_version;

  /** XML `<file_version>` (the XML schema version, e.g. "1.2.0").  Empty for text files. */
  std::string m_file_version;

  /** Crystal material name from the XML `<material>` element (e.g. "NaI").  Empty for text files
   (which carry only the material index at param 59) - use #materialName to resolve either. */
  std::string m_material_name;

  /** Verbatim label text for each parameter line as read from a text file (1-based).  Empty
   entries fall back to #canonicalLabel when writing text.  Lets #toText reproduce a text file
   byte-for-byte even when a particular install used non-standard labels. */
  std::array<std::string, sm_num_params + 1> m_labels{};

  /** String / non-numeric XML leaves preserved for round-trip, keyed by their slash-delimited
   node path (e.g. "photon_scatter/environment", "timing/compute_pileup",
   "anti_coincidence_shield/material").  Populated only when parsing XML. */
  std::map<std::string,std::string> m_xml_extras;


  //=============================== Raw parameter access ===============================

  /** Value column of parameter @p idx (1..84); returns 0 if out of range. */
  float value( int idx ) const;

  /** Integer / fit-flag column of parameter @p idx (1..84); returns 0 if out of range. */
  int intCol( int idx ) const;

  /** Set the value column of parameter @p idx (1..84); no-op if out of range.  Also grows
   #m_highest_param so a subsequently-written file includes the parameter. */
  void setValue( int idx, float v );

  /** Canonical GADRAS label for parameter @p idx (1..84), used when writing a text file if no
   verbatim label was captured.  Returns "" if out of range. */
  static const char *canonicalLabel( int idx );


  //=============================== Named accessors ===================================
  // All return GADRAS native units.  The trailing comment gives the parameter index and, where
  // useful, the meaning/units adapted from GADRAS_DRF_parameter_documentation.md.

  float length()        const { return value(10); }  //!< param 10: crystal box depth along axis, cm (doc §5).
  float width()         const { return value(11); }  //!< param 11: crystal box transverse width, cm (doc §5).
  float heightToWidth() const { return value(12); }  //!< param 12: box height = ratio*width (doc §5).
  float height()        const { return width() * heightToWidth(); } //!< derived transverse height, cm.
  float shapeFactor()   const { return value(13); }  //!< param 13: crystal-profile/chord model; big=box, small=rounded (doc §6).
  float solidAnglePercent() const { return value(9); } //!< param 9: fractional solid angle (%) at ref geometry (doc §7).
  float distanceCm()    const { return value(17); }  //!< param 17: reference source distance, cm.
  float setbackCm()     const { return value(40); }  //!< param 40: crystal recess behind reference point, cm (doc §7.3).
  float deadLayerMm()   const { return value(51); }  //!< param 51: crystal dead-layer thickness, mm (doc §8.3/§10.2).
  float effScalar()     const { return value(18); }  //!< param 18: constant efficiency grounding scalar (doc §7.3).

  float resOffset()   const { return value(6); }  //!< param 6: FWHM model, resolution @ E=0 (keV); sign selects low-E model (doc §11.1).
  float resFWHM661()  const { return value(7); }  //!< param 7: FWHM model, % FWHM @ 661 keV (doc §11.1).
  float resPower()    const { return value(8); }  //!< param 8: FWHM model, resolution power (doc §11.1).

  float lldKeV()       const { return value(35); } //!< param 35: lower-level discriminator energy, keV (doc §11.2).
  float lldSharpness() const { return value(72); } //!< param 72: LLD roll-off sharpness (doc §11.2).

  float lowerEnergyKeV() const { return static_cast<float>( intCol(62) ); } //!< param 62 int col: min valid energy, keV (doc §7.3).
  float upperEnergyKeV() const { return static_cast<float>( intCol(63) ); } //!< param 63 int col: max valid energy, keV (doc §7.3).

  // Peak skew / tailing (store-and-forward for the GADRAS peak-shape code; doc §11.5).
  float lowSkew()        const { return value(47); } //!< param 47: low-E-side skew amplitude.
  float highSkew()       const { return value(29); } //!< param 29: high-E-side skew amplitude.
  float lowSkewPower()   const { return value(57); } //!< param 57: low-E skew energy power.
  float highSkewPower()  const { return value(73); } //!< param 73: high-E skew energy power.
  float lowSkewExtent()  const { return value(71); } //!< param 71: low-E skew extent.
  float highSkewExtent() const { return value(77); } //!< param 77: high-E skew extent.

  /** InterSpec's round-face-equivalent crystal diameter, in cm (native, NOT PhysicalUnits):
   `2*sqrt(width*width*heightToWidth/pi)` - the diameter of a circle whose area equals the box
   transverse face (width x height).  Callers multiply by PhysicalUnits::cm. */
  float equivalentCircularDiameterCm() const;


  //=============================== Shielding (doc §8) ================================

  /** An attenuating housing layer: a single effective atomic number, areal density (g/cm2), and
   porosity (%).  See doc §8.2. */
  struct Attenuator
  {
    float atomicNumber = 0.0f;
    float arealDensity = 0.0f;   //!< g/cm2
    float porosityPercent = 0.0f;
  };

  Attenuator innerAttenuator() const;  //!< Front window/endcap: params 14 (Z), 15 (AD), 16 (porosity).
  Attenuator outerAttenuator() const;  //!< Outer housing:      params 19 (Z), 20 (AD), 32 (porosity).

  /** Directional side shield: one effective Z/AD plus per-face coverage percentages (doc §8.4). */
  struct SideShield
  {
    float atomicNumber = 0.0f;   //!< param 66
    float arealDensity = 0.0f;   //!< param 67, g/cm2
    float covPlusX = 0.0f;       //!< param 65, %
    float covMinusX = 0.0f;      //!< param 81, % (defaults to +X in old files)
    float covPlusY = 0.0f;       //!< param 82, %
    float covMinusY = 0.0f;      //!< param 83, %
  };
  SideShield sideShield() const;

  /** Back shield: one effective Z/AD plus a single coverage percentage (doc §8.4). */
  struct BackShield
  {
    float atomicNumber = 0.0f;   //!< param 69
    float arealDensity = 0.0f;   //!< param 70, g/cm2
    float coverage = 0.0f;       //!< param 68, %
  };
  BackShield backShield() const;


  //=============================== Material (doc §3) ================================

  /** GADRAS detector-material index (1..36) into the default material array (#materialTable).
   For text files this is the integer column of param 59 (0 = none/unknown "perfect" detector).
   For XML files it is resolved from the stored #m_material_name.  Returns 0 if unknown. */
  int materialIndex() const;

  /** Resolved crystal-material name: the XML `<material>` name if present, otherwise the name of
   #materialIndex from #materialTable.  Empty if neither is available. */
  std::string materialName() const;

  /** One row of the default GADRAS material table (`DetectorData.gadras`, doc §3.1.2).

   \note The shipped table is user-editable and can differ between GADRAS installs, so prefer
         resolving a material by NAME when one is available; use the numeric index only with the
         matching install. */
  struct MaterialInfo
  {
    int          index;          //!< 1-based table index.
    const char  *name;           //!< e.g. "NaI".
    const char  *formula;        //!< GADRAS formula syntax, e.g. "Cd0.96Zn0.04Te".
    float        density;        //!< g/cm3.
    float        resolution661;  //!< nominal % FWHM @ 661 keV.
    bool         specialCased;   //!< true if GADRAS has a name-keyed physics branch (doc §3.1.3).
  };

  /** The 36-entry default GADRAS material table (doc §3.1.2). */
  static const std::array<MaterialInfo, 36> &materialTable();

  /** Table row for a 1-based index, or nullptr if out of range / 0. */
  static const MaterialInfo *materialByIndex( int idx );

  /** Table row for a material name (case-insensitive), or nullptr if not found. */
  static const MaterialInfo *materialByName( const std::string &name );

  //=========================== Shape inference (doc §5.1) ===========================

  /** Recovered physical crystal shape class. */
  enum class Shape
  {
    Rectangular,       //!< Bar/log/panel (HeightToWidth != 1).
    Cylinder,          //!< Right circular cylinder (scintillators, square cross-section).
    CoaxialCylinder,   //!< Coaxial HPGe (bore geometry not recoverable from the file).
    Box,               //!< Cube/box solid-state block (CZT, CdTe, TlBr, Si).
    Unknown
  };

  /** A recovered shape plus its three characteristic dimensions (cm).  Meaning of dimA/B/C:
   - Rectangular: depth(length) x width x height(=heightToWidth*width)
   - Box:         length x width x width
   - Cylinder / CoaxialCylinder: dimA = radius (=width/2), dimB = height (=length), dimC unused (0). */
  struct InferredShape
  {
    Shape shape = Shape::Unknown;
    float dimA = 0.0f, dimB = 0.0f, dimC = 0.0f;
  };

  /** Guess the original physical crystal shape from the stored box dimensions, material, and
   shape factor, following the decision recipe in doc §5.1.  @p materialOverride, if non-empty,
   is used instead of #materialName for the material-family branch.  Recovered sizes are nominal
   (GADRAS often stores fitted effective dimensions). */
  InferredShape inferShape( const std::string &materialOverride = std::string() ) const;

  /** The normalized crystal chord-length distribution (doc §6.2), used by GADRAS's analytic
   interaction-probability integral.  `lengths[i]` (i = 1..count) is the fraction of incident
   photons whose chord through the crystal box is `length*i/32` cm.  Provided as a utility; not
   used during parsing. */
  struct ChordDistribution
  {
    std::array<double,128> lengths{};
    int count = 0;
  };

  /** Reproduces the GADRAS `PathLengthComputer` chord-length ray trace (doc §6.2).  All lengths
   in cm; @p distance = param17 + param40(setback). */
  static ChordDistribution computeChordDistribution( double length, double width,
                                                     double heightToWidthRatio,
                                                     double shapeFactor, double distance );


  //=============================== Read / write ====================================

  /** Parse a Detector.dat from @p input, auto-detecting the text vs XML variant by sniffing the
   first line for "xml".  Throws std::runtime_error on failure.  The stream must be seekable
   (the first line is read to detect the variant, then the stream is rewound).

   */
  static GadrasDetectorDat fromStream( std::istream &input );

  /** Convenience: open @p path and parse it via #fromStream. */
  static GadrasDetectorDat fromFile( const std::string &path );

  /** Whether @p input plausibly IS a `Detector.dat`, for classifying an uploaded
   or dropped file of unknown type.

   Deliberately stricter than "#fromStream did not throw": the text parser needs
   only one well-formed `<index> <value>` line and skips everything else, so it
   accepts a great many CSVs and logs.  This requires either the XML variant's
   `<gamma_detector>` root, or a text file with a substantial number of parameter
   lines AND the two values every real detector has - a positive %FWHM at 661 keV
   (param 7) and a positive crystal width (param 11).

   Never throws; leaves @p input's read position where it found it.
   */
  static bool isCandidateDetectorDat( std::istream &input );

  /** Write this detector as a GADRAS XML `Detector.dat` (the real GADRAS schema).  Numeric
   leaves are emitted as `<value>` + `<varying>` under their proper sub-blocks. */
  void toXml( std::ostream &out ) const;

  /** Write this detector as a legacy fixed-column text `Detector.dat`.  Uses verbatim labels
   (#m_labels) where available, else #canonicalLabel; inserts the `Version` line after param 64
   when #m_response_version is set. */
  void toText( std::ostream &out ) const;

private:
  /** Parse the legacy fixed-column text variant from @p input into this instance. */
  void parseText( std::istream &input );

  /** Parse the XML variant from @p input into this instance. */
  void parseXml( std::istream &input );
};//struct GadrasDetectorDat

#endif //GadrasDetectorDat_h
