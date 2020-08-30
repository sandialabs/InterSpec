#ifndef EnergyCalMultiFile_h
#define EnergyCalMultiFile_h
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
#include <memory>
#include <utility>

#include <Wt/WDialog>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>

// Forward Declarations
class PeakDef;
class AuxWindow;
class EnergyCalTool;
class SpectraFileModel;
class EnergyCalMultiFileModel;


/** Tool to fit peaks from multiple files to come up with the best energy calibration.
 
 Currently, for the case of multiple detectors in the system, this tool just uses the display energy
 calibration for each set of peaks, and assumes it should be the same across all files.  This isnt
 totally correct, but good enough for the majority of use cases problably.
 */

/*
class EnergyCalMultiFile : public Wt::WContainerWidget
{
public:
  EnergyCalMultiFile( EnergyCalTool *cal, AuxWindow *parent );
  virtual ~EnergyCalMultiFile();
  
  void doFit();
  
  void handleFinish( Wt::WDialog::DialogCode result );
protected:
  void updateCoefDisplay();
  
  EnergyCalTool *m_calibrator;
  EnergyCalMultiFileModel *m_model;
  std::vector<Wt::WCheckBox *> m_fitFor;
  std::vector<Wt::WLineEdit *> m_coefvals;
  Wt::WPushButton *m_use;
  Wt::WPushButton *m_cancel;
  Wt::WPushButton *m_fit;
  Wt::WTextArea *m_fitSumary;
  
  
  std::vector<double> m_calVal;
  std::vector<double> m_calUncert;
};//class EnergyCalMultiFile



class EnergyCalMultiFileModel : public  Wt::WAbstractItemModel
{
public:
  EnergyCalMultiFileModel( EnergyCalTool *calibrator, Wt::WObject *parent = 0 );
  virtual ~EnergyCalMultiFileModel();
  
  virtual Wt::WModelIndex index( int row, int column,
                                 const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;
  virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  virtual boost::any data( const Wt::WModelIndex &index,
                           int role = Wt::DisplayRole ) const;
  virtual bool setData( const Wt::WModelIndex &index,
                        const boost::any &value, int role = Wt::EditRole );
  virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const;
  virtual boost::any headerData( int section,
                                 Wt::Orientation orientation = Wt::Horizontal,
                                 int role = Wt::DisplayRole) const;
  void refreshData();
  
protected:
  std::vector<std::vector< std::pair<bool,std::shared_ptr<const PeakDef> > > > m_peaks;
  EnergyCalTool *m_calibrator;
  SpectraFileModel *m_fileModel;
  
  friend class EnergyCalMultiFile;
};//class EnergyCalMultiFileModel
*/


#endif //EnergyCalMultiFile_h
