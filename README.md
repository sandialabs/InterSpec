# InterSpec
InterSpec is a web-based application to assist in analyzing spectral nuclear radiation data, using a peak-based methodology.
Common uses include identifying nuclides present, determining source activity, shielding amounts, source age, or other nuclear reactions present.
InterSpec also provides a number of other tools useful for analyzing radiation data including spectral file format converting,
dose rate calculations, interactive nuclide decay and reference information, gamma cross section calculations, and more.
InterSpec can open data files from most common spectral radiation detectors (e.g., most NaI, HPGe, LaBr, 
CLYC, CsI, etc. based systems) and assist in their analysis.
InterSpec can also be built as a stand-alone (e.g., not Internet connection necessary) native app (using either OS provided WebView, or [Electron](https://github.com/electron/electron) for display) for Windows, Linux, macOS, iOS, or Android.


## Getting Started

If you just care about using a pre-compiled executable, head over to 
[releases](https://github.com/sandialabs/interspec/releases) and download the latest.
A brief getting started guide can be found in [tutorials/brief_analysis_intro](https://github.com/sandialabs/InterSpec/tree/master/tutorials/brief_analysis_intro/); a more comprehsive tutorial is planned for the future.

If you would like to build the app, some instructions are bellow.

If you are developer of radiation analysis software, there is, in particular some code that 
may be useful for you:
1. [sandia_decay](https://github.com/sandialabs/interspec/tree/master/external_libs/sandia_decay): This small library reads in [sandia.decay.xray.xml](https://github.com/sandialabs/interspec/tree/master/data/sandia.decay.xray.xml) and allows you to retrieve the half-lives, xrays, gammas, alpha, betas, etc given off by any certain nuclide, or you can easily decay nuclides and find activities of the decay chain nuclides, the particles that will be given off at any time, and more.  An emphasis is placed on easy to use, but computationally efficient code.
2. [SpecUtils](https://github.com/sandialabs/interspec/tree/master/external_libs/SpecUtils): If you have to read in spectral files, including [N42](https://www.nist.gov/programs-projects/ansiieee-n4242-standard), [PCF](http://prod.sandia.gov/techlib/access-control.cgi/2017/179107.pdf), SPC, DAT, CHN, etc., this library may do what you want.  It decodes many file format variants, and saves to about a dozen different formats.  This library can be called from c++, Python, or Java, or if you want a very simple gui to do this, or a command line application, check out [cambio](https://hekili.ca.sandia.gov/cambio/) (which is powered by this library)


### Prerequisites

To compile, you need a c++11 compiler, and:
* [Wt](https://www.webtoolkit.eu/wt) versions 3.3.1 through 3.3.4 should work, but 3.3.4 is recommended.
* [boost](https://www.boost.org/) versions 1.48 through 1.65.1 will probably work, but development is done using 1.65.1.
* [cmake](https://cmake.org/) 

To make building a little easier, other required libraries, including 
[Cuba](http://www.feynarts.de/cuba), 
[Minuit2](https://github.com/root-project/root/tree/master/math/minuit2), 
[muparserx](https://github.com/beltoforion/muparserx), 
[rapidxml](http://rapidxml.sourceforge.net/), 
and [websocketpp](https://github.com/zaphoyd/websocketpp), are included in this git repository already in the  
[external_libs](https://github.com/sandialabs/interspec/tree/master/external_libs) directory, while 
the [js](https://github.com/sandialabs/interspec/tree/master/js) directory contains some ECMAScript libraries, including
[jQuery](https://jquery.org),
[qTip2](http://qtip2.com),
[SILK ICONS](http://www.famfamfam.com/lab/icons/silk), and
[D3](https://d3js.org).


### Installing

After installing boost, Wt, and cmake using either your package manager, or building from source, clone the 
InterSpec repository, and from the terminal, run cmake:

```
git clone https://github.com/sandialabs/interspec/
cd interspec
mkdir build
cd build
cmake ..
```

If you are on Windows, or prefer a GUI, running the CMake for building InterSpec should be like
most other CMake project; you will probably have to fill in paths to boost and Wt manually.


And then make:
```
make -j8
```

If all goes well, InterSpec default to building an executable that when ran from the command line, starts a local web server, you can access then access from your browser.  To run InterSpec, use a command like:
```
./bin/InterSpec.exe --docroot . --http-address 0.0.0.0 --http-port 8080 -c ./data/config/wt_config_localweb.xml
```
and then point your browser to [http://localhost:8080](http://localhost:8080).


If you would like build as a native-ish app, see the *BUILD_AS_OSX_APP*, or *BUILD_AS_ELECTRON_APP* CMake options. 
For building as a WebApp behind something like nginx or apache, see *BUILD_FOR_WEB_DEPLOYMENT*.
There are also a number of CMake options available to set to control which features get included in InterSpec.


Building for iOS and Android are both possible (see the [target](https://github.com/sandialabs/interspec/tree/master/target) directory), and these instructions will be updated in the future.

## Authors
Ethan Chan, **William Johnson**, Christian Morte, with
extensive additional support provided by Noel Nachtigal and Edward Walsh.

## License
This project is licensed under the LGPL v2.1 License - see the [LICENSE.md](LICENSE.md) file for details

## Copyright
Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.

 ## Disclaimer
```
DISCLAIMER OF LIABILITY NOTICE:
The United States Government shall not be liable or responsible for any maintenance,
updating or for correction of any errors in the SOFTWARE or subsequent approved version
releases.


THE INTERSPEC (SOFTWARE) AND ANY OF ITS SUBSEQUENT VERSION
RELEASES, SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF
ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT
NOT LIMITED TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO
SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY
WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY
THAT THE DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE
SOFTWARE. IN NO EVENT SHALL THE UNITED STATES GOVERNMENT OR ITS
CONTRACTORS OR SUBCONTRACTORS BE LIABLE FOR ANY DAMAGES,
INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR
CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY
WAY CONNECTED WITH THE SOFTWARE OR ANY OTHER PROVIDED
DOCUMENTATION, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT,
TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY
PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
SUSTAINED FROM, OR AROSE OUT OF THE RESULT OF, OR USE OF, THE
SOFTWARE OR ANY PROVIDED DOCUMENTATION. THE UNITED STATES
GOVERNMENT DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING
THIRD PARTY SOFTWARE, IF PRESENT IN THE SOFTWARE, AND DISTRIBUTES
IT "AS IS."
```

 ## Acknowledgement
This InterSpec Software was developed with funds from the Science and Technology Directorate of the Department of Homeland Security.