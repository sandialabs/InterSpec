#ifndef BatchGuiInputFile_h
#define BatchGuiInputFile_h
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
#include <memory>

#include <Wt/WContainerWidget.h>


// Forward declarations
class SpecMeas;
class NativeFloatSpinBox;
class D3SpectrumDisplayDiv;
class DetectorPeakResponse;
class FileDragUploadResource;

namespace Wt
{
  class WText;
  class WLabel;
  class WCheckBox;
  class WLineEdit;
}

namespace SpecUtils
{
  class Measurement;
}

class BatchGuiInputSpectrumFile : public Wt::WContainerWidget
{
public:
  enum class ShowPreviewOption
  {
    Show,
    DontShow
  };

protected:
  const std::string m_filename;
  const std::string m_display_name;
  const bool m_should_cleanup;
  const ShowPreviewOption m_show_preview;

  Wt::WContainerWidget *m_preview_container;
  D3SpectrumDisplayDiv *m_spectrum;
  
  bool m_preview_created;
  std::shared_ptr<SpecMeas> m_spec_meas;
  
  bool m_is_peaks_csv;
  
  Wt::Signal<bool> m_preview_created_signal;
  Wt::Signal<BatchGuiInputSpectrumFile *> m_remove_self_request_signal;

  /** Creates the widgets common to both constructors. */
  void init_gui( const std::string &display_name, const std::string &path_to_file );

public:
  BatchGuiInputSpectrumFile( const std::string display_name,
                  const std::string path_to_file,
                  const bool should_cleanup,
                  const ShowPreviewOption preview );

  /** Constructor for a spectrum file that is already parsed and resident in memory - e.g., a file
   the user has open in InterSpec.

   `path_to_file` is the on-disk path, if one is known; it is used only for the tool-tip and
   `path_to_file()`, and may be empty.  The file is never deleted by this class.

   Note that, like the disk-loading constructor, the preview is not created synchronously - it is
   created from a delayed callback, since creating the `D3SpectrumDisplayDiv` in the same event
   loop that creates this widget gives JS errors inside a `SimpleDialog` on the macOS app.
   */
  BatchGuiInputSpectrumFile( const std::string display_name,
                  const std::string path_to_file,
                  std::shared_ptr<SpecMeas> spec_meas,
                  const ShowPreviewOption preview );

  void set_spectrum( std::shared_ptr<SpecMeas> spec_meas, std::shared_ptr<int> status_ptr );
  
  Wt::Signal<bool> &preview_created_signal();
  
  ~BatchGuiInputSpectrumFile();
  
  std::shared_ptr<SpecMeas> spec_meas() const;
  
  bool is_peaks_csv() const;
  
  void requestRemoveSelf();
  
  const std::string &display_name() const;
  
  const std::string &path_to_file() const;
  
  Wt::Signal<BatchGuiInputSpectrumFile *> &remove_self_request();
};// class BatchGuiInputSpectrumFile


class BatchGuiInputFile : public Wt::WContainerWidget
{
protected:
  const std::string m_filename;
  const std::string m_display_name;
  const bool m_should_cleanup;
  
  Wt::WText *m_txt;
  
  Wt::Signal<BatchGuiInputFile *> m_remove_self_request_signal;
  
public:
  BatchGuiInputFile( const std::string display_name,
                      const std::string path_to_file,
                      const bool should_cleanup );
  
  ~BatchGuiInputFile();
  
  void requestRemoveSelf();
  
  const std::string &display_name() const;
  
  const std::string &path_to_file() const;
  
  Wt::Signal<BatchGuiInputFile *> &remove_self_request();
};// class BatchGuiInputFile

#endif // BatchGuiInputFile_h 
