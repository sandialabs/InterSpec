#!/usr/bin/python2.7

# SpecUtils: a library to parse, save, and manipulate gamma spectrum data files.
# 
# Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
# Government retains certain rights in this software.
# For questions contact William Johnson via email at wcjohns@sandia.gov, or
# alternative emails of interspec@sandia.gov, or srb@sandia.gov.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 

import SpecUtils

filename = "Cal.pcf"
info = SpecUtils.MeasurementInfo()

try:
    info.loadFile( filename, SpecUtils.ParserType.kAutoParser )
except RuntimeError as e:
    print "Failed to decode file: {0}".format( e )
    exit( 1 )

meass = info.measurements()
nmeas = info.numMeasurements()

if nmeas < 1:
    print filename, "didnt contain any spectroscopic data"
    exit( 1 )

if len(meass) != nmeas :
    print "Fatal Error: len(meass) != nmeas"
    exit( 1 )

print "There were", nmeas, "measurements (or records) in the file"

meas = meass[0]
#meas = info.measurement(0)

counts = meas.gammaCounts()
startime = meas.startTime()

numchannel = len(counts)
print "For measurment started at", startime, ":"
print numchannel, "channels, with a few mid channels of the first measurement having counts:"

print "\tChannel\tCounts"
for i in range(numchannel/2,min(numchannel/2+10,numchannel)):
    print "\t", i, "\t", counts[i]
print "With live time:", meas.liveTime(), "seconds, and total counts:", meas.gammaCountSum()


nenergy = 511
channel = meas.findGammaChannel(nenergy)
content = meas.gammaChannelContent( channel );
lenergy = meas.gammaChannelLower( channel );
uenergy = meas.gammaChannelUpper( channel );
print nenergy, "keV corresponds to channel ", channel, "which has", content, "counts, and energy range (", lenergy, ",", uenergy, ")"


lenergy = 400
uenergy = 800
gammaint = meas.gammaIntegral(lenergy,uenergy)
print "Between", lenergy, "and", uenergy, "keV, the sum of gamma counts is", gammaint

lchannel = 20
uchannel = 30
gammasum = meas.gammaChannelsSum(lchannel,uchannel)
print "Channels", lchannel, "through", uchannel, "summed give", gammasum, "gamma sums"

print "DetectorNumbers:", info.detectorNumbers()
print "SampleNumbers:", info.sampleNumbers()


summedmeas = info.sumMeasurements( [1,2], info.detectorNumbers() )
print "Summed measurement has liveTime=", summedmeas.liveTime()



savetoname = "Cal_pyconverted.chn"
f = open( savetoname, 'w' );

sampleNums = info.sampleNumbers()
detToUse = [0]
try:
    info.writeIntegerChn( f, sampleNums, detToUse )
except RuntimeError as e:
    print "Error writing Integer CHN file: {0}.".format( e )
    exit( 1 )
    
f.close()
print "Wrote", savetoname

savetoname = "Cal_pyconverted.pcf"
f = open( savetoname, 'w' );

try:
    info.writePcf( f )
except RuntimeError as e:
    print "Error writing PCF file: {0}.".format( e )
    exit( 1 )
    
f.close()
print "Wrote", savetoname


savetoname = "Cal_pyconverted.n42"
f = open( savetoname, 'w' );

try:
    info.write2011N42Xml( f )
except RuntimeError as e:
    print "Error writing 2011 N42 file: {0}.".format( e )
    exit( 1 )
    
f.close()
print "Wrote", savetoname

#still having trouble reading from python source when seeking is done by the reader
f = open( "Cal_pyconverted.pcf", 'r' );
rereadinfo = SpecUtils.MeasurementInfo()
try:
    rereadinfo.setInfoFromPcfFile( f );
except RuntimeError as e:
    print "Failed to decode the converted PCF file: {0}.".format( e )
    exit( 1 )

print "Read converted PCF file, that has", rereadinfo.numMeasurements(), " measurements"


