#ifndef AngleOutxImport_h
#define AngleOutxImport_h
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
#include <istream>
#include <utility>

class DetectorPeakResponse;

/** Reading of ANGLE 5 XML files (".outx", ".savx", ".detx", ".recx", ...).

 The format is documented at https://www.angle.me/file-formats.html; every file
 type shares one `<angle generator="ANGLE" version=… build=… units="mm|cm|in">`
 root, and the same sub-element grammar (`<detector>`, `<container>`,
 `<geometry>`, `<source>`, `<energies>`, `<referenceEfficiencyCurve>`,
 `<results>`, ...) appears in whichever files carry that data.  The structures
 below mirror that grammar, so a single parser handles all of the file types.

 Conventions:
   - all lengths are in InterSpec `PhysicalUnits` units (the file's mm/cm/inch
     value multiplied by the root element's unit); the two fields whose names
     end in "Cm" are the exception and are plain centimeters;
   - materials are kept exactly as the file gives them (a bare name for one of
     ANGLE's predefined materials, or a full inline definition), so the parser
     needs no `MaterialDB` - name resolution happens in `CeeLoUtils`;
   - counts are all left generic: nothing is assumed about how many layers,
     coatings, absorbing layers, regions or points a given file contains.
 */

/** The eight predefined ANGLE detector types.

 Many `<detector>` attributes are only defined for a subset of these, and the
 mapping onto a CeeLo geometry differs per type, so the raw string is turned
 into this enum once at parse time (`AngleOutxContents::detTypeEnum`). */
enum class AngleDetectorType
{
  ClosedEndCoaxHPGe,   //!< "Closed-end coaxial HPGe"   (spec type 1)
  TrueCoaxHPGe,        //!< "True coaxial HPGe"         (spec type 2)
  ClosedEndCoaxGeLi,   //!< "Closed-end coaxial Ge(Li)" (spec type 3)
  OpenEndCoaxGeLi,     //!< "Open-end coaxial Ge(Li)"   (spec type 4)
  PlanarLEPD,          //!< "Planar LEPD"               (spec type 5)
  Well,                //!< "Well"                      (spec type 6)
  NaI,                 //!< "NaI"                       (spec type 7)
  NaIWell,             //!< "NaI Well"                  (spec type 8)
  Unknown              //!< absent, or a string this version does not know
};//enum class AngleDetectorType

/** The exact ANGLE spelling of `type`, or "Unknown". */
const char *to_string( const AngleDetectorType type );

/** Maps an ANGLE `<detector type="...">` string onto #AngleDetectorType.
 Matching is case-insensitive and ignores surrounding whitespace; anything
 unrecognized gives `Unknown`. */
AngleDetectorType angleDetectorTypeFromString( const std::string &type );


/** One element of an ANGLE `<elements>` material: `massFraction` is a fraction
 in [0,1] (the file's mass-percent divided by 100). */
struct AngleElementFraction
{
  std::string symbol;
  double massFraction = 0.0;
};//struct AngleElementFraction

/** One element of an ANGLE `<compound>`: a stoichiometric atom count. */
struct AngleCompoundElement
{
  std::string symbol;
  int atoms = 0;
};//struct AngleCompoundElement

/** One ANGLE `<compound>`.  `massFraction` is a fraction in [0,1]; it is 1.0
 for the single-compound form (which has no `massFraction` attribute). */
struct AngleCompound
{
  double massFraction = 1.0;
  std::string chemicalFormula;
  std::vector<AngleCompoundElement> elements;
};//struct AngleCompound

/** An ANGLE `<material>`, in any of the contexts it appears in.

 It is either a bare reference to one of ANGLE's 19 predefined materials
 (`<material name="Water"/>` - name only, resolved downstream), or a full
 inline definition: a density plus exactly one of `<elements>`, `<compound>` or
 `<compounds>`.  Inline definitions are what let a layer's attenuation be
 modeled without depending on the name being in InterSpec's material database. */
struct AngleMaterial
{
  std::string name;
  /** Density in g/cm3; zero for a bare predefined reference. */
  double density_g_cm3 = 0.0;
  /** From `<elements>`; fractions in [0,1], normalized to sum to 1. */
  std::vector<AngleElementFraction> elements;
  /** From `<compound>` (one entry) or `<compounds>` (N entries, normalized). */
  std::vector<AngleCompound> compounds;

  /** True when the file's mass percentages did not add to 100 and had to be
   rescaled - the spec requires them to, so this means a hand edit or a corrupt
   file, and the composition may not be what was intended. */
  bool fractionsRescaled = false;

  /** True when nothing at all was given (no name and no composition). */
  bool empty() const;
  /** True when only a name was given, so the composition must come from the
   predefined-material table / InterSpec's `MaterialDB`. */
  bool isBareReference() const;
};//struct AngleMaterial


/** Where a concentric detector layer sits in the physical stack.  Kept so the
 CeeLo mapping can order and interpret the layers rather than relying on the
 order they happened to be parsed in. */
enum class AngleLayerRole
{
  Contact,                //!< `<contact>`, when it is not folded into the dead layer
  ReflectingLayer,        //!< `<reflectingLayer>`  (NaI / NaI Well)
  HousingInner,           //!< `<housing><topLower>` or `<sideInner>` (one each)
  HousingOuter,           //!< `<housing><topUpper>` or `<sideOuter>` (one each)
  VacuumInner,            //!< `<vacuum>` inside an antimicrophonic shield
  AntimicrophonicShield,  //!< `<antimicrophonicShield>`
  VacuumOuter,            //!< `<vacuum>` (the whole gap when there is no shield)
  EndCap,                 //!< `<endCap>`
  EndCapWindow,           //!< `<endCap><window>`
  EndCapCoating,          //!< `<endCap><coatings><coating>`
  OpticalCoupling,        //!< `<opticalCoupling>`      (behind the crystal)
  PmtWall,                //!< `<photomultiplierTube><wall>`   (behind)
  PmtWindow,              //!< `<photomultiplierTube><window>` (behind)
  Other
};//enum class AngleLayerRole

const char *to_string( const AngleLayerRole role );


/** The `<well>` cavity of a well-type detector (spec types 6 and 8). */
struct AngleWell
{
  bool present = false;
  double depth = 0.0;
  double radius = 0.0;
};//struct AngleWell

/** `<contact>` - the electrical contact on an HPGe crystal (spec types 1, 2).
 A thin implanted contact of the crystal material reads as dead crystal; a
 thick one of a different material is an attenuating layer. */
struct AngleContact
{
  AngleMaterial material;
  double front = 0.0;   //!< `topThickness` (type 1 only)
  double side = 0.0;    //!< `sideThickness`
  /** True when this contact should be counted as dead crystal rather than as
   an attenuating layer: it is unspecified, made of the crystal material, or
   thin enough (< 10 um) to be an implanted contact. */
  bool readsAsDeadCrystal = true;
};//struct AngleContact

/** `<contactPin>` - the pin down the crystal's bore (spec types 1-5). */
struct AngleContactPin
{
  AngleMaterial material;
  double radius = 0.0;
};//struct AngleContactPin


/** An ANGLE `<container>` - the cup the sample sits in.  This is between the
 source and the detector, NOT around the crystal. */
struct AngleContainer
{
  bool present = false;
  /** "Cylindrical" or "Marinelli". */
  std::string type;
  std::string name;
  std::string description;

  // <shape> - which attributes are used depends on `type`
  double innerRadius = 0.0;
  double sideThickness = 0.0;
  double bottomThickness = 0.0;   //!< Cylindrical (file may spell it "botomThickness")
  double footHeight = 0.0;        //!< Cylindrical
  double cavityRadius = 0.0;      //!< Marinelli
  double cavityDepth = 0.0;       //!< Marinelli
  double sideInnerThickness = 0.0;  //!< Marinelli
  double bottomUpperThickness = 0.0; //!< Marinelli
  double bottomLowerThickness = 0.0; //!< Marinelli

  AngleMaterial material;

  /** Up to two `<coatings><coating>` layers. */
  struct Coating
  {
    AngleMaterial material;
    double side = 0.0;
    double bottom = 0.0;        //!< Cylindrical
    double bottomUpper = 0.0;   //!< Marinelli
    double bottomLower = 0.0;   //!< Marinelli
  };//struct Coating

  std::vector<Coating> coatings;
};//struct AngleContainer

/** An ANGLE `<geometry>` - the holder the container rests on, plus any extra
 absorbing layers.  Also between the source and the detector. */
struct AngleHolderGeometry
{
  bool present = false;
  /** "General", "Marinelli" or "Well". */
  std::string type;
  std::string name;
  std::string description;

  double holderOuterRadius = 0.0;
  /** On-axis height of the holder, i.e. the standoff from the detector face to
   the bottom of whatever rests on it. */
  double holderHeight = 0.0;

  double capThickness = 0.0;
  AngleMaterial capMaterial;
  double wallThickness = 0.0;
  AngleMaterial wallMaterial;

  /** Up to five `<absorbingLayers><absorbingLayer>`. */
  struct AbsorbingLayer
  {
    AngleMaterial material;
    double top = 0.0;
    double bottom = 0.0;   //!< well-detector geometry only
    double side = 0.0;
  };//struct AbsorbingLayer

  std::vector<AbsorbingLayer> absorbingLayers;
};//struct AngleHolderGeometry

/** An ANGLE `<source>`.  A zero `height` means a point or disk source, in which
 case ANGLE does not require a material. */
struct AngleSource
{
  bool present = false;
  double radius = 0.0;
  double height = 0.0;
  AngleMaterial material;
};//struct AngleSource


/** One `<results><result>` row.  `efficiency`/`efficiencyPrecision` are only
 written when the calculation had a reference efficiency curve to transfer
 from; `solidAngle` is always present. */
struct AngleResultRow
{
  double energyKeV = 0.0;
  double solidAngle = 0.0;
  double solidAnglePrecision = 0.0;
  double efficiency = 0.0;
  double efficiencyPrecision = 0.0;
  bool hasEfficiency = false;
};//struct AngleResultRow

/** One `<cascadeSummingCorrections><nuclide><correction>` row: ANGLE's true
 coincidence summing correction at a single energy. */
struct AngleCascadeCorrection
{
  double energyKeV = 0.0;
  double value = 1.0;
  double branchingRatio = 0.0;
  double correctedBranchingRatio = 0.0;
};//struct AngleCascadeCorrection

struct AngleCascadeNuclide
{
  std::string nuclide;
  std::vector<AngleCascadeCorrection> corrections;
};//struct AngleCascadeNuclide


/** Plain-data result of parsing an ANGLE XML file.

 ANGLE files carry far more than the geometry-transferred efficiency curve
 (`<results>`) InterSpec historically imported: a full physical detector model
 (`<detector>`), the sample container and holder it was measured in
 (`<container>`, `<geometry>`, `<source>`), and a measured reference efficiency
 curve (`<referenceEfficiencyCurve>`).  This struct holds the raw contents so
 callers can either (Mode B) use `fixedGeomDrf` as today, or (Mode A) map the
 geometry onto a `ceelo::GeometryDescriptor` and anchor an efficiency transfer
 / MC on the reference curve.
 */
struct AngleOutxContents
{
  /** The fixed-geometry DRF built from `<results>` - the historical import
   product (Mode B).  Null when the file carries no per-energy efficiencies
   (e.g. a ".detx" detector definition, or an ".outx" computed without a
   reference curve, where `<result>` has solid angles only). */
  std::shared_ptr<DetectorPeakResponse> fixedGeomDrf;

  // --- provenance, from the `<angle>` root element ---
  std::string generator;   //!< always "ANGLE"
  std::string version;     //!< minimum ANGLE version needed to open the file
  std::string build;       //!< ANGLE version that wrote the file
  std::string unitsStr;    //!< "mm", "cm" or "in" as written in the file
  /** The `PhysicalUnits` length factor `unitsStr` resolved to; every length in
   this struct has already been multiplied by it. */
  double lengthUnitPU = 0.0;

  // --- detector geometry (PhysicalUnits length units) ---
  std::string detName;
  std::string detDescription;
  std::string detType;                 //!< raw `type` attribute
  AngleDetectorType detTypeEnum = AngleDetectorType::Unknown;

  double crystalRadius = 0.0;
  double crystalLength = 0.0;
  /** ANGLE `<crystal bulletizingRadius>`: the outer front edge is rounded off
   with a quarter-torus fillet of this radius, as HPGe crystals usually are.
   Zero means a sharp 90-degree edge.  Maps to
   `ceelo::GeometryDescriptor::bullet_radius_cm`. */
  double bulletizingRadius = 0.0;

  bool hasCore = false;
  double coreRadius = 0.0;
  double coreDepth = 0.0;
  /** ANGLE `<core rounded="yes">`: the core was drilled with a round-tipped
   bit, so its closed end is a hemisphere of the core radius rather than a flat
   bottom.  The depth is measured to the apex either way.  Maps to
   `ceelo::BoreHoleConfig::rounded_tip`. */
  bool coreRounded = false;

  /** The `<well>` cavity of a well-type detector. */
  AngleWell well;

  // --- inactive (dead) germanium, `<inactiveGe>` ---
  /** `topThickness` / `sideThickness` (spec types 1-5).  For a well detector
   (types 6, 8) these are the OUTWARD-facing pair - `topUpperThickness` and
   `sideOuterThickness` - since only those are on the path from a source in
   front of the detector. */
  double deadLayerFront = 0.0;
  double deadLayerSide = 0.0;
  /** `bottomThickness`, planar LEPD only.  Maps to `DeadLayerConfig::back`. */
  double deadLayerBack = 0.0;
  /** The raw well-detector spellings.  A well crystal is annular, so these are
   FOUR DISTINCT SURFACES, not a stack: `topUpper`/`sideOuter` face outward and
   become `deadLayerFront`/`deadLayerSide`, while `topLower`/`sideInner` line
   the well itself and are recorded but not modeled. */
  double deadLayerTopUpper = 0.0;
  double deadLayerTopLower = 0.0;
  double deadLayerSideInner = 0.0;
  double deadLayerSideOuter = 0.0;

  AngleContact contact;
  AngleContactPin contactPin;

  /** Vacuum gap in front of / beside the crystal - a spacing, not an
   attenuator (attenuation ~ 0).  These are the TOTALS over every `<vacuum>`
   attribute in that direction; the gap also appears in `layers` (with no
   material) so that it holds its place in the stack. */
  double vacuumFront = 0.0;
  double vacuumSide = 0.0;
  /** The inner/outer split of the gap, non-zero only when the file has an
   `<antimicrophonicShield>` (ANGLE only writes these attributes then), so the
   shield can be placed at the right radius / depth.  For a well detector the
   gap against the well floor (`topLower`/`bottomLower`) is a different surface
   and is not included here. */
  double vacuumFrontInner = 0.0;
  double vacuumFrontOuter = 0.0;
  double vacuumSideInner = 0.0;
  double vacuumSideOuter = 0.0;

  /** `<endCap><window>` geometry that is not carried by a layer: the window's
   own radius, the radius of the end-cap hole it covers, and the thickness of
   the annulus that clamps it.  The holder is a ring, not a disc, so it is not
   on the on-axis path and is reported rather than modeled. */
  double endCapWindowRadius = 0.0;
  double endCapWindowHoleRadius = 0.0;
  double endCapWindowHolderThickness = 0.0;

  /** One concentric endcap/window/housing layer.  `front`/`side` are the
   thicknesses (PhysicalUnits units). */
  struct Layer
  {
    AngleLayerRole role = AngleLayerRole::Other;
    AngleMaterial material;
    double front = 0.0;
    double side = 0.0;
    /** True for pieces that sit BEHIND the crystal (the optical coupling and
     photomultiplier tube of a NaI detector).  CeeLo models no back layers, so
     these are carried for fidelity but dropped - with a warning - when a
     `ceelo::GeometryDescriptor` is built. */
    bool behindCrystal = false;
  };//struct Layer

  /** The concentric stack around the crystal - reflecting layer, housing,
   vacuum gap, antimicrophonic shield, endcap, window and coatings - ordered
   innermost -> outermost, followed by any pieces behind the crystal.  A layer
   whose `role` is `VacuumInner`/`VacuumOuter` carries no material: it is the
   cryostat gap, present only so everything outside it sits at the right depth
   and radius. */
  std::vector<Layer> layers;

  // --- reference efficiency curve ---
  /** Measured reference efficiency points: (energy in keV, absolute
   efficiency), ascending in energy. */
  std::vector<std::pair<float,float>> referencePoints;

  /** Reference on-axis source distance, in cm (derived from the reference
   curve's `<holder height>`); approximate, meant to be user-editable. */
  double referenceDistanceCm = 0.0;

  std::string referenceCurveName;
  std::string referenceCurveDescription;
  /** `<referenceEfficiencyCurve><detector name>` - the detector the reference
   measurement was made on (a name reference only, no geometry). */
  std::string referenceDetectorName;

  /** The container / holder / source of the REFERENCE measurement.  Together
   with the sample ones below these give the source-side attenuation on both
   sides of an efficiency transfer (it cancels when they are the same, which is
   the common case). */
  AngleContainer referenceContainer;
  AngleHolderGeometry referenceGeometry;
  AngleSource referenceSource;

  // --- sample geometry (the top-level `<container>`/`<geometry>`/`<source>`) ---
  /** These describe the SAMPLE the file's `<results>` efficiency is for (not
   the reference curve).  They are not used to build the detector geometry -
   they sit between the source and the detector, not around the crystal. */
  AngleContainer sampleContainer;
  AngleHolderGeometry sampleGeometry;
  AngleSource sampleSource;

  /** On-axis standoff from the endcap front (i.e. the outer detector face) to
   the sample's NEAR (bottom) face - where the sample begins, not its center.
   This is the empty gap between detector and sample, taken from the sample
   holder's `<height>` (the sample rests on top of the holder).  For an extended
   source, add half of the source height to reach the source center; for a
   point/thin source the two coincide.  In cm. */
  double sampleDistanceCm = 0.0;

  /** On-axis standoff from the endcap front to the SOURCE's near face, in cm:
   `sampleDistanceCm` plus everything the container puts between the holder and
   the source (its foot, its bottom wall, and any coatings on it).  Equal to
   `sampleDistanceCm` when there is no container. */
  double sampleSourceDistanceCm = 0.0;

  // --- calculation inputs / outputs ---
  /** The `<energies>` the calculation was run at, in keV, and the name of that
   energy set. */
  std::vector<float> energies;
  std::string energiesName;
  /** `<precision>`: the requested calculation precision, in percent. */
  double precision = 0.0;
  /** `<elapsedTime>`: ANGLE's calculation time, in seconds. */
  double elapsedTimeSec = 0.0;

  /** Every `<results><result>` row as written, including the solid angles that
   are present even when no efficiencies are. */
  std::vector<AngleResultRow> results;

  /** `<cascadeSummingCorrections>`: ANGLE's true-coincidence-summing factors,
   per nuclide.  Not used to build a DRF; kept as a comparison dataset for
   InterSpec's own cascade-summing calculations. */
  std::vector<AngleCascadeNuclide> cascadeCorrections;

  // --- fitted region equations ---
  /** One fitted log-log efficiency region.  `equation` is the raw text of the
   XML comment ANGLE writes after each `<region>` (e.g.
   "log10 e = c0 + c1 * log10 Eg + c2 * (log10 Eg)^2 + ...").  Extracted only
   for a self-consistency dev-check; not interpreted otherwise. */
  struct FitRegion
  {
    double startKeV = 0.0;
    double endKeV = 0.0;
    int polynomOrder = 0;
    std::string equation;
  };//struct FitRegion

  std::vector<FitRegion> fitRegions;

  /** Whether a usable detector geometry / reference curve was parsed; gate
   Mode A availability on these. */
  bool hasGeometry = false;
  bool hasReference = false;

  /** Whether the parsed detector can be represented as a
   `ceelo::GeometryDescriptor` at all.  False for the detector types CeeLo has
   no shape for (a well cavity, or a true-coaxial through-bore); the
   fixed-geometry (Mode B) import still works for those. */
  bool modeASupported = false;
  /** When `modeASupported` is false, a short user-facing explanation of why. */
  std::string modeAObstruction;

  /** Non-fatal issues noticed while parsing (an unusual unit, an attribute the
   declared detector type does not define, a feature we cannot model, ...).
   Meant to be shown to the user, not logged and forgotten. */
  std::vector<std::string> parseNotes;
};//struct AngleOutxContents


namespace AngleOutx
{
  /** Parses any ANGLE 5 XML file into an #AngleOutxContents.

   Throws only when the input is not usable ANGLE XML at all (too large, not
   XML, no `<angle>` root, an unrecognized `units` value, a malformed number
   where one is required).  A well-formed file that simply lacks a section -
   a ".detx" with no `<results>`, an ".outx" computed without a reference curve
   - parses successfully with the corresponding fields left empty; check
   `hasGeometry` / `hasReference` / `results` / `fixedGeomDrf` for what is
   actually present.  `fixedGeomDrf` is NOT filled in by this function; see
   `DetectorPeakResponse::parseAngleOutxFileFull`.
   */
  AngleOutxContents parse( std::istream &input,
                           const size_t max_bytes = 10u * 1024u * 1024u );
}//namespace AngleOutx

#endif //AngleOutxImport_h
