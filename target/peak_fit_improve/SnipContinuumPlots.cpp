/* InterSpec: dev tool - export foreground + production SNIP continuum for visual assessment.

 For each configured (detector x live-time) directory of injected test spectra, load each PCF,
 take the first gamma measurement as the injected foreground, run the SAME SNIP continuum
 estimator used in the app (`estimateContinuum` - both the legacy fixed-window form and the
 FWHM-window variants under evaluation), downsample to a manageable display resolution, and
 emit one compact JSON file per config.

 The JSON is consumed by a small standalone HTML viewer (snip_plots_viewer) so the SNIP continuum
 can be eyeballed against the data across the whole test set.

 Not part of the app build; added as its own executable in this directory's CMakeLists.
 */

#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>

#include <Wt/WColor>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"
#include "SpecUtils/D3SpectrumExport.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/FitPeaksForNuclides.h"

using namespace std;

namespace
{
  // One (detector x live-time) config to export.
  struct Config
  {
    string label;      // short id used for the output file name
    string title;      // human title shown in the viewer
    string dir;        // directory of PCFs to walk
    float (*fwhm_fcn)( const float energy );  // nominal FWHM(E) for this detector class
    string d3_variant; // which continuum variant to draw on this config's D3 page
  };

  // One FWHM-window SNIP variant to export alongside the legacy fixed-window continuum.
  struct SnipVariant
  {
    string key;      // series key in the JSON ("sn" object)
    double num_fwhm; // clip half-window, in local FWHM
    int order;       // clipping filter order (2/4/6)
    int smooth_hw;   // presmooth half-width (0=off, 1=3ch, 2=5ch, ...) before clipping
    bool lls;        // clip in log-log-sqrt space (classic SNIP transform)
    bool restrict;   // restrict clipping to the valid spectroscopic-extent energy range
  };

  // Round to keep the JSON small; counts at display resolution are large so integer rounding is fine.
  string num( const double v )
  {
    char buf[32];
    snprintf( buf, sizeof(buf), "%.0f", v );
    return buf;
  }

  string energy_str( const double v )
  {
    char buf[32];
    snprintf( buf, sizeof(buf), "%.2f", v );
    return buf;
  }

  string json_escape( const string &s )
  {
    string out;
    for( const char c : s )
    {
      if( c == '"' || c == '\\' )
        out += '\\';
      if( static_cast<unsigned char>(c) >= 0x20 )
        out += c;
    }
    return out;
  }

  // One spectrum's contents for the D3 (SpectrumChartD3.js) page.
  struct D3ChartInput
  {
    string name;
    shared_ptr<const SpecUtils::Measurement> foreground;
    shared_ptr<const SpecUtils::Measurement> continuum;  // SNIP continuum, shown as a 2nd line
    double extent_lo = 0.0;
    double extent_hi = 0.0;
  };

  // Draw a vertical marker at `energy` on the chart via the reference-line mechanism (same
  // approach write_spectroscopic_extent_html uses).
  bool add_extent_ref_line( std::map<std::string,std::string> &lines, const string &key,
                            const double energy )
  {
    if( !(energy > 0.0) )
      return false;

    RefLineInput input;
    input.m_input_txt = std::to_string( energy ) + " keV";
    input.m_color = Wt::WColor( 255, 165, 0 );  // orange
    input.m_showGammas = true;
    input.m_showXrays = false;
    input.m_showAlphas = false;
    input.m_showBetas = false;
    input.m_promptLinesOnly = false;
    input.m_lower_br_cutt_off = 0.0;

    const std::shared_ptr<ReferenceLineInfo> ref = ReferenceLineInfo::generateRefLineInfo( input );
    if( !ref || (ref->m_validity != ReferenceLineInfo::InputValidity::Valid) )
      return false;

    string json;
    ref->toJson( json );
    lines[key] = json;
    return true;
  }

  // Write one InterSpec-style HTML page (SpectrumChartD3.js) for a config: one zoomable chart per
  // spectrum, each showing the foreground plus the SNIP continuum as a second spectrum, with
  // vertical orange lines at the lower/upper spectroscopic extent.  Mirrors
  // FitPeaksForNuclideDev::write_spectroscopic_extent_html.
  void write_snip_d3_html( const string &title, const string &variant_desc,
                           const vector<D3ChartInput> &charts, const string &outpath )
  {
    ofstream html_output( outpath );
    if( !html_output.is_open() )
    {
      cerr << "Error: could not open '" << outpath << "' for writing." << endl;
      return;
    }

    try
    {
      D3SpectrumExport::write_html_page_header( html_output, title, "InterSpec_resources" );
    }catch( std::exception &e )
    {
      cerr << "\n\nError: " << e.what() << "\nSymlink InterSpec_resources into the CWD." << endl;
      return;
    }

    html_output << "<body>\n"
                << "<style>fieldset{width: 90vw; margin: 20px auto 0;}"
                << " legend{font-size: 13px;}</style>\n"
                << "<h2>SNIP continuum &mdash; " << title << "</h2>\n"
                << "<p style=\"width:90vw;margin:0 auto;color:#555;font-size:13px\">"
                << "Foreground (black) + SNIP continuum <b>" << variant_desc << "</b>"
                << " (blue); orange vertical lines = lower/upper spectroscopic extent."
                << " Drag to zoom, double-click to reset.</p>\n";

    size_t chart_counter = 0;
    for( const D3ChartInput &ci : charts )
    {
      if( !ci.foreground || (ci.foreground->num_gamma_channels() < 16) )
        continue;

      std::map<std::string,std::string> reference_lines_json;
      const bool have_extent = (ci.extent_lo > 0.0) && (ci.extent_hi > ci.extent_lo);
      if( have_extent )
      {
        add_extent_ref_line( reference_lines_json, "Extent Min", ci.extent_lo );
        add_extent_ref_line( reference_lines_json, "Extent Max", ci.extent_hi );
      }

      const float xMin = static_cast<float>( ci.foreground->gamma_energy_min() );
      const float xMax = static_cast<float>( ci.foreground->gamma_energy_max() );

      const D3SpectrumExport::D3SpectrumChartOptions options(
        ci.name, "Energy (keV)", "Counts/Channel",
        /*dataTitle=*/"", /*useLogYAxis=*/true,
        /*showVerticalGridLines=*/false, /*showHorizontalGridLines=*/false,
        /*legendEnabled=*/true, /*compactXAxis=*/true,
        /*showPeakUserLabels=*/false, /*showPeakEnergyLabels=*/false,
        /*showPeakNuclideLabels=*/false, /*showPeakNuclideEnergyLabels=*/false,
        /*showEscapePeakMarker=*/false, /*showComptonPeakMarker=*/false,
        /*showComptonEdgeMarker=*/false, /*showSumPeakMarker=*/false,
        /*backgroundSubtract=*/false,
        xMin, xMax, reference_lines_json
      );

      D3SpectrumExport::D3SpectrumOptions fg_opts;
      fg_opts.line_color = "black";
      fg_opts.title = "data";
      fg_opts.display_scale_factor = 1.0;
      fg_opts.spectrum_type = SpecUtils::SpectrumType::Foreground;

      D3SpectrumExport::D3SpectrumOptions cont_opts;
      cont_opts.line_color = "#2a78d6";
      cont_opts.title = "SNIP continuum";
      cont_opts.display_scale_factor = 1.0;
      cont_opts.spectrum_type = SpecUtils::SpectrumType::Background;

      const string div_id = "chart_" + std::to_string( chart_counter );
      chart_counter++;

      string legend_text = ci.name;
      if( have_extent )
        legend_text += " [extent: " + std::to_string( static_cast<int>( ci.extent_lo + 0.5 ) )
                     + " - " + std::to_string( static_cast<int>( ci.extent_hi + 0.5 ) ) + " keV]";

      html_output << "<fieldset>\n<legend>" << legend_text << "</legend>\n"
                  << "<div id=\"" << div_id << "\" class=\"chart\" oncontextmenu=\"return false;\"></div>\n"
                  << "<script>\n";

      D3SpectrumExport::write_js_for_chart( html_output, div_id, options.m_dataTitle,
                                            options.m_xAxisTitle, options.m_yAxisTitle );

      std::vector< std::pair<const SpecUtils::Measurement *, D3SpectrumExport::D3SpectrumOptions> > measurements;
      measurements.emplace_back( ci.foreground.get(), fg_opts );
      if( ci.continuum )
        measurements.emplace_back( ci.continuum.get(), cont_opts );

      D3SpectrumExport::write_and_set_data_for_chart( html_output, div_id, measurements );

      html_output << "\nconst resizeChart" << (chart_counter - 1) << " = function(){\n"
                  << "  let height = window.innerHeight;\n"
                  << "  let width = document.documentElement.clientWidth;\n"
                  << "  let el = spec_chart_" << div_id << ".chart;\n"
                  << "  el.style.width = 0.8*width + \"px\";\n"
                  << "  el.style.height = Math.min(500,Math.max(250, Math.min(0.4*width,height-175))) + \"px\";\n"
                  << "  el.style.marginLeft = 0.05*width + \"px\";\n"
                  << "  el.style.marginRight = 0.05*width + \"px\";\n"
                  << "  spec_chart_" << div_id << ".handleResize();\n"
                  << "};\n"
                  << "window.addEventListener('resize', resizeChart" << (chart_counter - 1) << ");\n";

      D3SpectrumExport::write_set_options_for_chart( html_output, div_id, options );
      html_output << "spec_chart_" << div_id << ".setReferenceLines( reference_lines_" << div_id << " );\n";
      html_output << "resizeChart" << (chart_counter - 1) << "();\n</script>\n</fieldset>\n";
    }//for( each chart )

    html_output << "</body>\n</html>\n";
    html_output.close();

    cout << "  wrote " << outpath << " (" << chart_counter << " charts)" << endl;
  }//write_snip_d3_html
}//namespace


int main( int argc, char **argv )
{
  // Optional flags anywhere on the command line: --no-d3 (skip the heavy SpectrumChartD3 pages),
  // --d3-variant=<key> (which continuum variant to draw on the D3 pages).
  bool write_d3 = true;
  string d3_variant_override;  // empty = use each config's per-class d3_variant
  vector<string> positional;
  for( int i = 1; i < argc; ++i )
  {
    const string a = argv[i];
    if( a == "--no-d3" )
      write_d3 = false;
    else if( a.rfind("--d3-variant=",0) == 0 )
      d3_variant_override = a.substr( 13 );
    else
      positional.push_back( a );
  }

  const string base = (positional.size() > 0)
    ? positional[0]
    : "/Users/wcjohns/coding/InterSpec_peak_fit_improve/peak_fit_accuracy_inject_compact";
  const string outdir = (positional.size() > 1) ? positional[1] : "snip_plots";

  SpecUtils::create_directory( outdir );

  // D3 continuum per class: HPGe -> k2.0slx (LLS, wider), NotHPGe -> k1.5s7x (7ch smooth, no LLS).
  const string hpge_d3 = "k2.0slx";
  const string lowmed_d3 = "k1.5s7x";
  const vector<Config> configs = {
    { "detx_30s",  "Detective-X (HPGe) - Livermore - 30 s",   base + "/Detective-X/Livermore/30_seconds",  &PeakFitUtils::hpge_fwhm_fcn, hpge_d3 },
    { "detx_300s", "Detective-X (HPGe) - Livermore - 300 s",  base + "/Detective-X/Livermore/300_seconds", &PeakFitUtils::hpge_fwhm_fcn, hpge_d3 },
    { "nai_30s",   "identiFINDER R500 (NaI) - Livermore - 30 s",  base + "/IdentiFINDER-R500-NaI/Livermore/30_seconds",  &PeakFitUtils::nai_fwhm_fcn, lowmed_d3 },
    { "nai_300s",  "identiFINDER R500 (NaI) - Livermore - 300 s", base + "/IdentiFINDER-R500-NaI/Livermore/300_seconds", &PeakFitUtils::nai_fwhm_fcn, lowmed_d3 },
    { "labr_30s",  "identiFINDER LaBr3 - Livermore - 30 s",   base + "/IdentiFINDER-LaBr3/Livermore/30_seconds",  &PeakFitUtils::labr_fwhm_fcn, lowmed_d3 },
    { "labr_300s", "identiFINDER LaBr3 - Livermore - 300 s",  base + "/IdentiFINDER-LaBr3/Livermore/300_seconds", &PeakFitUtils::labr_fwhm_fcn, lowmed_d3 },
    { "czt_30s",   "H3D M400 (CZT) - Livermore - 30 s",   base + "/CZT_H3D_M400_ORNL_25cm/Livermore/30_seconds",  &PeakFitUtils::czt_fwhm_fcn, lowmed_d3 },
    { "czt_300s",  "H3D M400 (CZT) - Livermore - 300 s",  base + "/CZT_H3D_M400_ORNL_25cm/Livermore/300_seconds", &PeakFitUtils::czt_fwhm_fcn, lowmed_d3 },
  };

  // FWHM-window variants (all order-2, the right filter for FWHM-scaled windows).  The candidate
  // per-class favorites from the first round were k1.0 (low/med) and k1.5+smooth (HPGe); the "x"
  // variants additionally restrict clipping to the valid spectroscopic extent, to test whether
  // that fixes the sub-extent detector-turn-on problem and lets ONE variant serve all classes.
  // Per-class candidates (LLS helps HPGe tall peaks but wrecks NaI Compton-edge tracking; wider
  // smoothing helps sparse-high-E CZT; NaI/LaBr want the narrow no-LLS estimator).
  const vector<SnipVariant> variants = {
    { "k1.5sx",   1.5, 2, 1, false, true },  // NaI/LaBr pick (3ch smooth, no LLS)
    { "k1.5slx",  1.5, 2, 1, true,  true },  // HPGe pick (+ LLS)
    { "k2.0slx",  2.0, 2, 1, true,  true },  // HPGe alt (wider + LLS)
    { "k1.5s5x",  1.5, 2, 2, false, true },  // CZT pick? (5ch smooth, no LLS)
    { "k1.5s7x",  1.5, 2, 3, false, true },  // CZT alt (7ch smooth, no LLS)
    { "k2.0s5x",  2.0, 2, 2, false, true },  // CZT alt (wider + 5ch smooth, no LLS)
  };

  const size_t target_display_bins = 2048;

  for( const Config &cfg : configs )
  {
    vector<string> files = SpecUtils::recursive_ls( cfg.dir, ".pcf" );
    std::sort( begin(files), end(files) );

    if( files.empty() )
    {
      cerr << "No PCFs under '" << cfg.dir << "' - skipping " << cfg.label << endl;
      continue;
    }

    // Unique energy axes (spectra in a config almost always share one calibration).
    vector<vector<double>> axes;

    // Per-spectrum records, built as JSON fragments.
    vector<string> spectra_json;

    // Per-spectrum inputs for the zoomable SpectrumChartD3.js page (foreground + the unified
    // candidate continuum k1.5sx + extent).
    vector<D3ChartInput> d3_charts;

    size_t n_ok = 0;
    for( const string &filename : files )
    {
      SpecUtils::SpecFile sf;
      if( !sf.load_file( filename, SpecUtils::ParserType::Auto ) )
      {
        cerr << "  failed to load " << filename << endl;
        continue;
      }

      // First gamma measurement is the injected foreground (matches the fitter harness).
      shared_ptr<const SpecUtils::Measurement> fg;
      for( size_t i = 0; i < sf.num_measurements(); ++i )
      {
        const shared_ptr<const SpecUtils::Measurement> m = sf.measurement( i );
        if( m && m->gamma_counts() && !m->gamma_counts()->empty()
            && m->energy_calibration() && m->energy_calibration()->valid() )
        {
          fg = m;
          break;
        }
      }

      if( !fg )
        continue;

      const size_t nchan = fg->num_gamma_channels();
      if( nchan < 16 )
        continue;

      // Production SNIP continuum (same call the app uses).
      shared_ptr<SpecUtils::Measurement> snip;
      try
      {
        snip = estimateContinuum( fg );
      }catch( const std::exception &e )
      {
        cerr << "  estimateContinuum failed for " << filename << ": " << e.what() << endl;
        continue;
      }

      if( !snip || (snip->num_gamma_channels() != nchan) )
        continue;

      // FWHM-window variants (the tuned estimator under evaluation).
      const std::function<double(double)> fwhm_at
                    = [&cfg]( const double energy ) -> double {
        return cfg.fwhm_fcn( static_cast<float>(energy) );
      };

      // Valid spectroscopic extent (keV), for the "restrict"-flagged variants.
      const std::pair<double,double> extent = FitPeaksForNuclides::find_valid_energy_range( fg );

      vector<shared_ptr<SpecUtils::Measurement>> variant_snips( variants.size() );
      bool variant_failed = false;
      for( size_t v = 0; v < variants.size(); ++v )
      {
        const double lo = variants[v].restrict ? extent.first : 0.0;
        const double hi = variants[v].restrict ? extent.second : 0.0;
        try
        {
          variant_snips[v] = estimateContinuum( fg, fwhm_at, variants[v].num_fwhm,
                                    variants[v].order, variants[v].smooth_hw, variants[v].lls,
                                    lo, hi );
        }catch( const std::exception &e )
        {
          cerr << "  estimateContinuum(" << variants[v].key << ") failed for "
               << filename << ": " << e.what() << endl;
          variant_failed = true;
          break;
        }

        if( !variant_snips[v] || (variant_snips[v]->num_gamma_channels() != nchan) )
        {
          variant_failed = true;
          break;
        }
      }//for( each variant )

      if( variant_failed )
        continue;

      const vector<float> &fg_counts = *fg->gamma_counts();
      const vector<float> &snip_counts = *snip->gamma_counts();

      const size_t factor = std::max<size_t>( 1, static_cast<size_t>(
          std::llround( static_cast<double>(nchan) / static_cast<double>(target_display_bins) ) ) );
      const size_t ndisp = (nchan + factor - 1) / factor;

      vector<double> disp_e( ndisp ), disp_fg( ndisp, 0.0 ), disp_snip( ndisp, 0.0 );
      vector<vector<double>> disp_var( variants.size(), vector<double>(ndisp, 0.0) );
      for( size_t d = 0; d < ndisp; ++d )
      {
        const size_t c0 = d * factor;
        const size_t c1 = std::min( nchan, c0 + factor );  // exclusive
        double fsum = 0.0, ssum = 0.0;
        for( size_t c = c0; c < c1; ++c )
        {
          fsum += fg_counts[c];
          ssum += snip_counts[c];
        }
        disp_fg[d] = fsum;
        disp_snip[d] = ssum;
        disp_e[d] = 0.5 * ( fg->gamma_channel_lower(c0) + fg->gamma_channel_upper(c1 - 1) );

        for( size_t v = 0; v < variants.size(); ++v )
        {
          const vector<float> &vc = *variant_snips[v]->gamma_counts();
          double vsum = 0.0;
          for( size_t c = c0; c < c1; ++c )
            vsum += vc[c];
          disp_var[v][d] = vsum;
        }
      }

      // Find or add a matching energy axis.
      int axis_index = -1;
      for( size_t a = 0; a < axes.size(); ++a )
      {
        if( axes[a].size() == disp_e.size()
            && std::fabs( axes[a].front() - disp_e.front() ) < 1.0e-3
            && std::fabs( axes[a].back() - disp_e.back() ) < 1.0e-3 )
        {
          axis_index = static_cast<int>( a );
          break;
        }
      }
      if( axis_index < 0 )
      {
        axis_index = static_cast<int>( axes.size() );
        axes.push_back( disp_e );
      }

      string name = SpecUtils::filename( filename );
      const size_t dot = name.rfind( ".pcf" );
      if( dot != string::npos )
        name = name.substr( 0, dot );

      string rec = "{\"n\":\"" + json_escape(name) + "\",\"lt\":"
                   + num( fg->live_time() ) + ",\"ax\":" + std::to_string(axis_index)
                   + ",\"fg\":[";
      for( size_t d = 0; d < ndisp; ++d )
        rec += (d ? "," : "") + num( disp_fg[d] );
      rec += "],\"bg\":[";
      for( size_t d = 0; d < ndisp; ++d )
        rec += (d ? "," : "") + num( disp_snip[d] );
      rec += "],\"sn\":{";
      for( size_t v = 0; v < variants.size(); ++v )
      {
        rec += (v ? ",\"" : "\"") + variants[v].key + "\":[";
        for( size_t d = 0; d < ndisp; ++d )
          rec += (d ? "," : "") + num( disp_var[v][d] );
        rec += "]";
      }
      // Spectroscopic extent (keV): [lower, upper] vertical markers for the viewer.
      rec += "},\"ex\":[" + energy_str( extent.first ) + "," + energy_str( extent.second ) + "]}";

      spectra_json.push_back( std::move(rec) );

      // Retain the class-appropriate continuum variant for the zoomable D3 page.
      const string d3_key = d3_variant_override.empty() ? cfg.d3_variant : d3_variant_override;
      shared_ptr<const SpecUtils::Measurement> unified_cont;
      for( size_t v = 0; v < variants.size(); ++v )
      {
        if( variants[v].key == d3_key )
          unified_cont = variant_snips[v];
      }
      d3_charts.push_back( D3ChartInput{ name, fg, unified_cont, extent.first, extent.second } );

      ++n_ok;
    }//for( each file )

    // Zoomable SpectrumChartD3.js page for this config (foreground + SNIP continuum + extent).
    if( write_d3 )
    {
      const string d3_key = d3_variant_override.empty() ? cfg.d3_variant : d3_variant_override;
      write_snip_d3_html( cfg.title, d3_key, d3_charts,
                          outdir + "/snip_d3_" + cfg.label + ".html" );
    }

    const string outpath = outdir + "/" + cfg.label + ".json";
    ofstream out( outpath.c_str() );
    if( !out )
    {
      cerr << "Could not open " << outpath << " for writing" << endl;
      continue;
    }

    out << "{\"title\":\"" << json_escape(cfg.title) << "\",\"variants\":[";
    for( size_t v = 0; v < variants.size(); ++v )
      out << (v ? ",\"" : "\"") << variants[v].key << "\"";
    out << "],\"axes\":[";
    for( size_t a = 0; a < axes.size(); ++a )
    {
      out << (a ? "," : "") << "[";
      for( size_t i = 0; i < axes[a].size(); ++i )
        out << (i ? "," : "") << energy_str( axes[a][i] );
      out << "]";
    }
    out << "],\"spectra\":[";
    for( size_t s = 0; s < spectra_json.size(); ++s )
      out << (s ? ",\n" : "\n") << spectra_json[s];
    out << "\n]}";
    out.close();

    cout << cfg.label << ": " << n_ok << " spectra -> " << outpath << endl;
  }//for( each config )

  return 0;
}
