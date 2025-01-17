@echo off
setlocal enableDelayedExpansion

rem This scripts builds the Boost, Wt, zlib, Eigen, Ceres-solver, and wxWidgets dependencies for InterSpec on Windows
rem This script is meant to be used from within the "x86 Native Tools Command Prompt for VS 2022" terminal
rem Only the Release build of libraries will be built by this script (the x64 script builds both)9
rem Usage: dep_build_msvc2022_x86.bat "Path to build at" "Path to install to"

rem This is my first cmd.exe .bat script, so I'm sure it could be done much better

IF "%2"=="" GOTO invalidarg
IF not "%3"=="" GOTO invalidarg

pushd .

rem We will assume this script is in the patch directory
set PATCH_DIR=%~dp0

set ORIG_DIR=%CD%

echo Checking of tgz/zip file hashes is temporarily disabled

if not exist "%1" (
  MKDIR "%1" && (
    echo "Created build dir %1"
  ) || (
    echo "Failed to make dir %1"
    popd
    exit /b -2
  )
) else (
    echo "Build directory %1 already exists."
)

pushd .
cd %1
set BUILD_DIR=%CD%
popd


if not exist "%2" (
  MKDIR "%2" && (
    echo "Created install dir %2"
  ) || (
    echo "Failed to make dir %2"
    popd
    exit /b -3
  )
) else (
    echo "Install directory %2 already exists."
)

pushd .
cd %2
set MY_PREFIX=%CD%
popd

echo "Will build in %BUILD_DIR%"
echo "Will install to %MY_PREFIX%"
echo "Will patch from %PATCH_DIR%"

cd %BUILD_DIR%

set errorlevel=

rem Build/install boost
set BOOST_TAR=boost_1_78_0.tar.gz
set BOOST_DIR=boost_1_78_0
set BOOST_BUILT_FILE=built_%BOOST_DIR%
set BOOST_REQUIRED_SHA256=f22143b5528e081123c3c5ed437e92f648fe69748e95fa6e2bd41484e2986cc3

if not exist %BOOST_BUILT_FILE% (

    if not exist %BOOST_TAR% (
        curl -L https://sourceforge.net/projects/boost/files/boost/1.78.0/boost_1_78_0.zip/download --output %BOOST_TAR% && (
            echo Downloaded Boost
        ) || (
            echo Error downloading boost
            GOTO :cmderr
        )
    ) else (
        echo %BOOST_TAR% already downloaded
    )
    
    set "BOOST_SHA256="
    for /f %%A in ('certutil -hashfile "%BOOST_TAR%" SHA256 ^| find /i /v ":" ') do set "BOOST_SHA256=%%A"
    
    if not "!BOOST_SHA256!"=="%BOOST_REQUIRED_SHA256%" (
        echo Invalid hash of boost.  Expected "%BOOST_REQUIRED_SHA256%" and got "!BOOST_SHA256!"
rem     GOTO :cmderr
    )

    if not exist %BOOST_DIR% (
        tar -xzvf %BOOST_TAR% && (
            echo Un-tarred %BOOST_TAR%
        ) || (
            echo "Failed to unzip boost"
            GOTO :cmderr
        )
    ) else (
        echo "Boost was already unzipped"
    )

    cd "%BOOST_DIR%"

    echo "Running boost bootstrap.bat"
    (
        set foundErr=1
        cmd /c .\bootstrap.bat

        if errorlevel 0 if not errorlevel 1 set "foundErr="
        if defined foundErr (
            echo "Failed to run bootstrap.bat"
            GOTO :cmderr
        )
    )


    echo "Building boost release"
    .\b2.exe runtime-link=static link=static threading=multi variant=release address-model=32 architecture=x86 --prefix="%MY_PREFIX%" --build-dir=win_build_release -j8 install && (
        rmdir /s /q win_build_release
        echo Done building boost release
    ) || (
        echo "Failed to run build boost release"
        GOTO :cmderr
    )

    echo "Built boost!"

    cd %BUILD_DIR%

    rmdir /s /q "%BOOST_DIR%" && (
        echo Removed "%BOOST_DIR%" directory
    ) || (
        echo Failed to remove "%BOOST_DIR%" directory
    )

    echo "Built boost" > %BOOST_BUILT_FILE%
) else (
    echo "Boost was already built (%BOOST_BUILT_FILE% existed)"
)


rem Build/install Zlib
set ZLIB_TAR="zlib-1.2.13.tar.gz"
set ZLIB_DIR="zlib-1.2.13"
set ZLIB_BUILT_FILE=built_%ZLIB_DIR%
set ZLIB_REQUIRED_SHA256=b3a24de97a8fdbc835b9833169501030b8977031bcb54b3b3ac13740f846ab30

if not exist %ZLIB_BUILT_FILE% (
    rem if not exist %ZLIB_TAR% (
        rem curl -L https://zlib.net/fossils/zlib-1.2.13.tar.gz --output %ZLIB_TAR%
    rem    curl -L https://zlib.net/zlib-1.2.13.tar.gz --output %ZLIB_TAR%
    rem ) else (
    rem    echo %ZLIB_TAR% already downloaded
    rem )

    rem set "ZLIB_SHA256="
    rem for /f %%A in ('certutil -hashfile "%ZLIB_TAR%" SHA256 ^| find /i /v ":" ') do set "ZLIB_SHA256=%%A"
    
    rem if not "!ZLIB_SHA256!"=="%ZLIB_REQUIRED_SHA256%" (
    rem    echo Invalid hash of zlib.  Expected "%ZLIB_REQUIRED_SHA256%" and got "!ZLIB_SHA256!"
rem     GOTO :cmderr
    rem )

    rem tar -xzvf %ZLIB_TAR% && (
    rem    echo Untarred zlib
    rem ) || (
    rem    echo "Failed to untar zlib"
    rem    GOTO :cmderr
    rem )

    git clone https://github.com/madler/zlib.git %ZLIB_DIR% && (
        echo Cloned into zlib
    ) || (
        echo "Failed to clone into zlib"
        GOTO :cmderr
    )

    cd %ZLIB_DIR%

    rem checkout zlib 1.2.13
    git checkout 04f42ceca40f73e2978b50e93806c2a18c1281fc && (
        echo Checked out zlib 1.2.13
    ) || (
        echo "Failed to checkout wanted zlib commit"
        GOTO :cmderr
    )

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX="%MY_PREFIX%" -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>" -A Win32 .. && (
        echo Configured zlib Release
    ) || (
        echo "Failed to cmake configure zlib"
        GOTO :cmderr
    )


    cmake --build . --config Release --target install && (
        echo Built/installed zlib Release
    ) || (
        echo "Failed to build zlib Release"
        GOTO :cmderr
    )

    rem Get rid of zlib dll's, otherwise CMake used to have InterSpec link
    rem against those, instead of zlibstatic.lib and zlibstaticd.lib, although
    rem I think I fixed this - but didnt test
    del "%MY_PREFIX%\bin\zlib.dll"
    del "%MY_PREFIX%\bin\zlibd.dll"
    del "%MY_PREFIX%\lib\zlib.lib"
    del "%MY_PREFIX%\lib\zlibd.lib"

    echo "Built zlib!"
    cd %BUILD_DIR%

    rmdir /s /q "%ZLIB_DIR%" && (
        echo Removed "%ZLIB_DIR%" directory
    ) || (
        echo Failed to remove "%ZLIB_DIR%" directory
    )


    echo "Built zlib" > %ZLIB_BUILT_FILE%
) else (
    echo "Zlib was already built (%ZLIB_BUILT_FILE% existed)"
)
    
    
rem Build/install Wt 
set WT_DIR=wt-3.7.1
set WT_BUILT_FILE=built_%WT_DIR%

rem git clone https://github.com/emweb/wt.git %WT_DIR%
rem cd %WT_DIR%
rem git checkout b84925215d2b45879cf20c0cb340c4e7960d0c53
rem 

if not exist %WT_BUILT_FILE% (
    if not exist %WT_DIR% (
        git clone https://github.com/emweb/wt.git %WT_DIR% && (
            echo Cloned into Wt
        ) || (
            echo "Failed to clone into Wt"
            GOTO :cmderr
        )

        cd %WT_DIR%

        rem checkout Wt 3.7.1
        git checkout b84925215d2b45879cf20c0cb340c4e7960d0c53 && (
            echo Checked out Wt 3.7.1
        ) || (
            echo "Failed to checkout wanted commit"
            GOTO :cmderr
        )

        echo About to patch Wt
        git apply --reject --ignore-space-change --ignore-whitespace "%PATCH_DIR%\wt\3.7.1\NormalBuild\wt_3.7.1_git.patch" && (
            echo Patched Wt
        ) || (
            echo "Failed to patch Wt"
            GOTO :cmderr
        )
    ) else (
        echo %WT_DIR% already cloned and patched
        cd %WT_DIR%
    )
    

    mkdir build_msvc2022
    cd build_msvc2022

    cmake -A Win32 -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="%MY_PREFIX%" -DBoost_INCLUDE_DIR="%MY_PREFIX%/include" -DBOOST_PREFIX="%MY_PREFIX%" -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="%MY_PREFIX%" -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DINSTALL_FINDWT_CMAKE_FILE=ON -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++11" -DINSTALL_FINDWT_CMAKE_FILE=OFF -DCONFIGURATION=data/config/wt_config_electron.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR="%MY_PREFIX%/etc/wt" -DBoost_USE_STATIC_RUNTIME=ON .. && (
        echo configured Wt
    ) || (
        echo "Failed to cmake configure Wt"
        GOTO :cmderr
    )

    echo Will build Wt
    cmake --build . --config Release --target install && (
        echo Built and installed Wt Release
    ) || (
        echo "Failed to build Release Wt"
        GOTO :cmderr
    )

    echo "Built Wt!"
    cd %BUILD_DIR%

    rmdir /s /q "%WT_DIR%" && (
        echo Removed "%WT_DIR%" directory
    ) || (
        echo Failed to remove "%WT_DIR%" directory
    )

    echo "Built Wt" > %WT_BUILT_FILE%
) else (
    echo Wt was already built (%WT_BUILT_FILE% existed)
)

rem Build/install Eigen
set EIGEN_DIR="eigen-3.x"
set EIGEN_BUILT_FILE=built_%EIGEN_DIR%
set EIGEN_GIT_HASH=2e76277bd049f7bec36b0f908c69734a42c5234f

if not exist %EIGEN_BUILT_FILE% (
    git clone --recursive https://gitlab.com/libeigen/eigen.git --branch master --single-branch --depth 1 %EIGEN_DIR% && (
        echo Cloned into Eigen
    ) || (
        echo "Failed to clone into Eigen"
        GOTO :cmderr
    )

    cd %EIGEN_DIR%

    git fetch --depth 1 origin %EIGEN_GIT_HASH% && (
        echo Fetched Eigen
    ) || (
        echo "Failed to fetch wanted Eigen commit"
        GOTO :cmderr
    )

    git checkout %EIGEN_GIT_HASH% && (
        echo Checked out wanted Eigen version
    ) || (
        echo "Failed to checkout wanted Eigen commit"
        GOTO :cmderr
    )

    git submodule update --init --recursive && (
        echo Updated Eigen submodules
    ) || (
        echo "Failed to update Eigen submodules"
        GOTO :cmderr
    )
    
    cmake -A Win32 -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DCMAKE_BUILD_TYPE=Release -DEIGEN_MPL2_ONLY=1 -DEIGEN_BUILD_SHARED_LIBS=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_TESTING=OFF .. && (
        echo Configured Eigen
    ) || (
        echo Failed to configure Eigen
        GOTO :cmderr
    )

    cmake --build . --config Release --target install && (
        echo Installed Eigen
    ) || (
        echo Failed to install Eigen
        GOTO :cmderr
    )

    echo "Built Eigen!"
    cd %BUILD_DIR%

    rmdir /s /q "%EIGEN_DIR%" && (
        echo Removed "%EIGEN_DIR%" directory
    ) || (
        echo Failed to remove "%EIGEN_DIR%" directory
    )

    echo "Built Eigen" > %EIGEN_BUILT_FILE%
) else (
    echo Eigen was already built (%EIGEN_BUILT_FILE% existed)
)

rem Build/install Ceres Solver
set CERES_DIR="ceres-solver-2.2.0"
set CERES_BUILT_FILE=built_%CERES_DIR%
set CERES_GIT_HASH=85331393dc0dff09f6fb9903ab0c4bfa3e134b01

if not exist %CERES_BUILT_FILE% (

    git clone --recursive https://github.com/ceres-solver/ceres-solver.git --branch master --single-branch --depth 1 %CERES_DIR% && (
        echo Cloned into ceres-solver
    ) || (
        echo "Failed to clone into ceres-solver"
        GOTO :cmderr
    )

    cd %CERES_DIR%

    rem checkout version 2.2.0, Oct 12, 2023.  We need to fetch, then checkout, then bring submodules to listed versions
    git fetch --depth 1 origin %CERES_GIT_HASH% && (
        echo Fetched ceres 2.2.0
    ) || (
        echo "Failed to fetch wanted ceres commit"
        GOTO :cmderr
    )

    git checkout %CERES_GIT_HASH% && (
        echo Checked out ceres 2.2.0
    ) || (
        echo "Failed to checkout wanted ceres commit"
        GOTO :cmderr
    )

    git submodule update --init --recursive && (
        echo Updated Ceres submodules
    ) || (
        echo "Failed to update Ceres submodules"
        GOTO :cmderr
    )

    mkdir build_msvc_rel
    cd build_msvc_rel

    cmake -A Win32 -DCMAKE_PREFIX_PATH=%MY_PREFIX% -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DMINIGLOG=ON -DGFLAGS=OFF -DACCELERATESPARSE=OFF -DUSE_CUDA=OFF -DEXPORT_BUILD_DIR=ON -DBUILD_TESTING=ON -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded -DMSVC_USE_STATIC_CRT=ON .. && (
        echo Configured Release Ceres
    ) || (
        echo Failed to configure Release ceres-solver
        GOTO :cmderr
    )

    cmake --build . --config Release --target install && (
        echo Built and installed Release Ceres
    ) || (
        echo Failed to build Release ceres-solver
        GOTO :cmderr
    )

    echo "Built Ceres!"
    cd %BUILD_DIR%

    rmdir /s /q "%CERES_DIR%" && (
        echo Removed "%CERES_DIR%" directory
    ) || (
        echo Failed to remove "%CERES_DIR%" directory
    )

    echo "Built Ceres" > %CERES_BUILT_FILE%
) else (
    echo Ceres-Solver was already built (%CERES_BUILT_FILE% existed)
)

rem Build/install wxWidgets
set WX_TAR="wxWidgets-3.2.1.zip"
set WX_DIR="wxWidgets-3.2.1"
set WX_BUILT_FILE=built_%WX_DIR%
set WX_REQUIRED_SHA256=cc8868c3c8ec4eddaf659a8b81589a3d83126d5afde012350f61031a607a56d8

if not exist %WX_BUILT_FILE% (
    mkdir %WX_DIR%
    cd %WX_DIR%
    curl -L https://github.com/wxWidgets/wxWidgets/releases/download/v3.2.1/%WX_TAR% --output %WX_TAR% && (
        echo Downloaded wxWidgets
    ) || (
        echo Failed to download wxWidgets
        GOTO :cmderr
    )

    set "WX_SHA256="
    for /f %%A in ('certutil -hashfile "%WX_TAR%" SHA256 ^| find /i /v ":" ') do set "WX_SHA256=%%A"
    
    echo Temporarily not checking hash
    rem if not "!WX_SHA256!"=="%WX_REQUIRED_SHA256%" (
    rem     echo Invalid hash of ceres.  Expected "%WX_REQUIRED_SHA256%" and got "!WX_SHA256!"
    rem     GOTO :cmderr
    rem )

    tar -xzvf %WX_TAR% && (
        echo Unzipped wxWidgets
    ) || (
        echo Failed to unzip wxWidgets
        GOTO :cmderr
    )

    mkdir build_msvc
    cd build_msvc

    rem TODO: we could/should turn off a lot of wxWidgets components
    cmake -A Win32 -DCMAKE_PREFIX_PATH=%MY_PREFIX% -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DwxUSE_WEBVIEW_EDGE=ON -DwxUSE_WEBVIEW_EDGE_STATIC=ON -DwxBUILD_USE_STATIC_RUNTIME=ON -DwxBUILD_SHARED=OFF .. && (
        echo Configured wxWidgets
    ) || (
        echo Failed to configure wxWidgets
        GOTO :cmderr
    )

    cmake --build . --config Release --target install && (
        echo Installed Debug wxWidgets
    ) || (
        echo Failed to install Release wxWidgets
        GOTO :cmderr
    )

    rem wxWidgets doesnt also install WebView2 stuff, so we'll do it manually
    copy packages\Microsoft.Web.WebView2.1.0.705.50\build\native\x86\WebView2LoaderStatic.lib "%MY_PREFIX%\lib" && (
        echo Successfully copied WebView2LoaderStatic.lib to lib dir
    ) || (
        echo Failed to copy WebView2LoaderStatic.lib to lib dir
        GOTO :cmderr
    )

    rem I dont think we actually need these next three things, so I wont bother doing error checking
    copy packages\Microsoft.Web.WebView2.1.0.705.50\build\native\x86\WebView2Guid.lib "%MY_PREFIX%\lib"
    copy packages\Microsoft.Web.WebView2.1.0.705.50\build\native\include\WebView2.h "%MY_PREFIX%\include"
    copy packages\Microsoft.Web.WebView2.1.0.705.50\build\native\include\WebView2EnvironmentOptions.h "%MY_PREFIX%\include"

    echo "Built wxWidgets!"
    cd %BUILD_DIR%

    rmdir /s /q "%WX_DIR%" && (
        echo Removed "%WX_DIR%" directory
    ) || (
        echo Failed to remove "%WX_DIR%" directory
    )

    
    echo "Built wxWidgets" > %WX_BUILT_FILE%
) else (
    echo wxWidgets was already built (%WX_BUILT_FILE% existed)
)

popd
echo "Completed Succesfully"
exit /b 0
goto :EOF

:invalidarg
echo "Usage: build_dep_msvc2022.bat <temp path> <install path>"
exit /b 1
goto :EOF

:cmderr
echo "Error building %1, code %2"
popd
exit /b 2
goto :EOF
