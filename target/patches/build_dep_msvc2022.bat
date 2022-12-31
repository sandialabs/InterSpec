@echo off

rem This scripts builds the Boost, Wt, zlib, Eigen, Ceres-solver, and wxWidgets dependencies for InterSpec on Windows
rem This script is meant to be used from within the "x64 Native Tools Command Prompt for VS 2022" terminal
rem Usage: dep_build_msvc2022.bat "Path to build at" "Path to install to"

rem This is my first cmd.exe .bat script, so I'm sure it could be done much better

IF "%2"=="" GOTO invalidarg
IF not "%3"=="" GOTO invalidarg

pushd .

rem We will assume this script is in the patch directory
set PATCH_DIR=%~dp0

set ORIG_DIR=%CD%


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
set BOOST_REQUIRED_MD5=e193e5089060ed6ce5145c8eb05e67e3

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

    echo Im not sure boost hash checking is working
    set "BOOST_MD5=" & for /F "skip=1 delims=" %%H in ('
       2^> nul CertUtil -hashfile "%BOOST_TAR%" MD5
    ') do if not defined BOOST_MD5 set "BOOST_MD5=%%H"

    if not "%BOOST_MD5%"=="%BOOST_REQUIRED_MD5%" (
        echo Invalid hash of boost.  Expected "%BOOST_REQUIRED_MD5%" and got "%BOOST_MD5%"
        GOTO :cmderr
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


    echo "Building boost"
    .\b2.exe runtime-link=static link=static threading=multi variant=release address-model=64 architecture=x86 --prefix="%MY_PREFIX%" --build-dir=win_build -j8 install && (
        echo Done building boost
    ) || (
        echo "Failed to run build boost"
        GOTO :cmderr
    )

    echo "Built boost!"

    cd %BUILD_DIR%
    echo "Built boost" > %BOOST_BUILT_FILE%
) else (
    echo "Boost was already built (%BOOST_BUILT_FILE% existed)"
)


rem Build/install Zlib
set ZLIB_TAR="zlib-1.2.13.tar.gz"
set ZLIB_DIR="zlib-1.2.13"
set ZLIB_BUILT_FILE=built_%ZLIB_DIR%
set ZLIB_REQUIRED_MD5=9b8aa094c4e5765dabf4da391f00d15c

if not exist %ZLIB_BUILT_FILE% (
    if not exist %ZLIB_TAR% (
        curl -L https://zlib.net/fossils/%ZLIB_TAR% --output %ZLIB_TAR%
    ) else (
        echo %ZLIB_TAR% already downloaded
    )

    echo Im not sure zlib hash checking is working
    set "ZLIB_MD5="
    for /f "skip=1 delims=" %%a in ('certutil -hashfile %ZLIB_TAR% MD5') do if not defined ZLIB_MD5 set "ZLIB_MD5=%%a"

    if not "%ZLIB_MD5%"=="%ZLIB_REQUIRED_MD5%" (
        echo Invalid hash of ZLIB.  Expected "%ZLIB_REQUIRED_MD5%" and got "%ZLIB_MD5%"
        GOTO :cmderr
    )

    tar -xzvf %ZLIB_TAR%
    if %errorlevel% neq 0 (
        echo "Failed to untar zlib"
        GOTO :cmderr
    )

    cd %ZLIB_DIR%
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX="%MY_PREFIX%" -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>" ..
    if %errorlevel% neq 0 (
        echo "Failed to cmake configure zlib"
        GOTO :cmderr
    )

    cmake --build . --config Debug --target install
    if %errorlevel% neq 0 (
        echo "Failed to build zlib debug"
        GOTO :cmderr
    )


    cmake --build . --config Release --target install
    if %errorlevel% neq 0 (
        echo "Failed to build zlib Release"
        GOTO :cmderr
    )

    echo "Built zlib!"
    cd %BUILD_DIR%
    echo "Built zlib" > %ZLIB_BUILT_FILE%
) else (
    echo "Zlib was already built (%ZLIB_BUILT_FILE% existed)"
)
    
    
rem Build/install Wt 
set WT_DIR=wt-3.7.1
set WT_BUILT_FILE=built_%WT_DIR%
set WT_REQUIRED_MD5=0f6635a213cc310fc4ec69ca48f2f46d

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

        set WT_PATCH_FILE="%PATCH_DIR%\wt\3.7.1\NormalBuild\wt_3.7.1_git.patch"
        git apply --reject --ignore-space-change --ignore-whitespace %WT_PATCH_FILE% && (
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

    cmake -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="%MY_PREFIX%" -DBoost_INCLUDE_DIR="%MY_PREFIX%/include" -DBOOST_PREFIX="%MY_PREFIX%" -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="%MY_PREFIX%" -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DINSTALL_FINDWT_CMAKE_FILE=ON -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++11" -DINSTALL_FINDWT_CMAKE_FILE=OFF -DCONFIGURATION=data/config/wt_config_electron.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR="%MY_PREFIX%/etc/wt" -DBoost_USE_STATIC_RUNTIME=ON .. && (
        echo configured Wt
    ) || (
        echo "Failed to cmake configure Wt"
        GOTO :cmderr
    )

    cmake --build . --config Release --target install && (
        echo Built and installed Wt Release
    ) || (
        echo "Failed to build Release Wt"
        GOTO :cmderr
    )

    cmake --build . --config Debug --target install && (
        echo Built and installed Wt Debug
    ) || (
        echo "Failed to build Debug Wt"
        GOTO :cmderr
    )

    echo "Built Wt!"
    cd %BUILD_DIR%
    echo "Built Wt" > %WT_BUILT_FILE%
) else (
    echo Wt was already built (%WT_BUILT_FILE% existed)
)

rem Build/install Eigen
set EIGEN_TAR="eigen-3.4.0.tar.gz"
set EIGEN_DIR="eigen-3.4.0"
set EIGEN_BUILT_FILE=built_%EIGEN_DIR%
set EIGEN_REQUIRED_MD5=9b8aa094c4e5765dabf4da391f00d15c

if not exist %EIGEN_BUILT_FILE% (
    curl -L https://gitlab.com/libeigen/eigen/-/archive/3.4.0/%EIGEN_TAR% --output %EIGEN_TAR% && (
        echo Downloaded Eigen
    ) || (
        echo Failed to download Eigen
        GOTO :cmderr
    )

    echo Im not sure Eigen hash checking is working
    set "EIGEN_MD5="
    for /f "skip=1 delims=" %%a in ('certutil -hashfile %EIGEN_MD5% MD5') do if not defined EIGEN_MD5 set "EIGEN_MD5=%%a"

    if not "%EIGEN_MD5%"=="%EIGEN_REQUIRED_MD5%" (
        echo Invalid hash of Eigen.  Expected "%EIGEN_REQUIRED_MD5%" and got "%EIGEN_MD5%"
        GOTO :cmderr
    )

    tar -xzvf %EIGEN_TAR% && (
        echo "Untarred Eigen"
    ) || (
        echo Failed to untar Eigen
        GOTO :cmderr
    )

    cd %EIGEN_DIR%
    mkdir build
    cd build
    
    cmake -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DCMAKE_BUILD_TYPE=Release -DEIGEN_MPL2_ONLY=1 -DEIGEN_BUILD_SHARED_LIBS=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_TESTING=OFF .. && (
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
    echo "Built Eigen" > %EIGEN_BUILT_FILE%
) else (
    echo Eigen was already built (%EIGEN_BUILT_FILE% existed)
)

rem Build/install Ceres Solver
set CERES_TAR="ceres-solver-2.1.0.tar.gz"
set CERES_DIR="ceres-solver-2.1.0"
set CERES_BUILT_FILE=built_%CERES_DIR%
set CERES_REQUIRED_MD5=9b8aa094c4e5765dabf4da391f00d15c

if not exist %CERES_BUILT_FILE% (

    curl -L http://ceres-solver.org/%CERES_TAR% --output %CERES_TAR% && (
        echo Downloaded Eigen
    ) || (
        echo Failed to download ceres-solver
        GOTO :cmderr
    )

    echo Need to check Ceres hash

    tar -xzvf %CERES_TAR% && (
        echo Untarred Ceres
    ) || (
        echo Failed to untar ceres-solver
        GOTO :cmderr
    )

    cd %CERES_DIR%
    mkdir build_msvc
    cd build_msvc

    cmake -DCMAKE_PREFIX_PATH=%MY_PREFIX% -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DMINIGLOG=ON -DGFLAGS=OFF -DCXSPARSE=OFF -DACCELERATESPARSE=OFF -DCUDA=OFF -DEXPORT_BUILD_DIR=ON -DBUILD_TESTING=ON -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDebug -DMSVC_USE_STATIC_CRT=ON .. && (
        echo Configured Debug Ceres
    ) || (
        echo Failed to configure Debug ceres-solver
        GOTO :cmderr
    )

    cmake --build . --config Debug --target install && (
        echo Built and installed Debug Cerres
    ) || (
        echo Failed to build Debug ceres-solver
        GOTO :cmderr
    )

    cmake -DCMAKE_PREFIX_PATH=%MY_PREFIX% -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DMINIGLOG=ON -DGFLAGS=OFF -DCXSPARSE=OFF -DACCELERATESPARSE=OFF -DCUDA=OFF -DEXPORT_BUILD_DIR=ON -DBUILD_TESTING=ON -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded -DMSVC_USE_STATIC_CRT=ON .. && (
        echo Configured Release Ceres
    ) || (
        echo Failed to configure Release ceres-solver
        GOTO :cmderr
    )

    cmake --build . --config Release --target install && (
        echo Built and installed Release Cerres
    ) || (
        echo Failed to build Release ceres-solver
        GOTO :cmderr
    )

    echo "Built Ceres!"
    cd %BUILD_DIR%
    echo "Built Ceres" > %CERES_BUILT_FILE%
) else (
    echo Ceres-Solver was already built (%CERES_BUILT_FILE% existed)
)

rem Build/install wxWidgets
set WX_TAR="wxWidgets-3.2.1.zip"
set WX_DIR="wxWidgets-3.2.1"
set WX_BUILT_FILE=built_%WX_DIR%
set WX_REQUIRED_MD5=9b8aa094c4e5765dabf4da391f00d15c

if not exist %WX_BUILT_FILE% (
    mkdir %WX_DIR%
    cd %WX_DIR%
    curl -L https://github.com/wxWidgets/wxWidgets/releases/download/v3.2.1/%WX_TAR% --output %WX_TAR% && (
        echo Downloaded wxWidgets
    ) || (
        echo Failed to download wxWidgets
        GOTO :cmderr
    )

    rem TODO: check hash
    echo Need to check wxWidgets hash

    tar -xzvf %WX_TAR% && (
        echo Unzipped wxWidgets
    ) || (
        echo Failed to unzip wxWidgets
        GOTO :cmderr
    )

    mkdir build_msvc
    cd build_msvc

    rem # TODO: we could/should turn off a lot of wxWidgets components
    cmake -DCMAKE_PREFIX_PATH=%MY_PREFIX% -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DwxUSE_WEBVIEW_EDGE=ON -DwxUSE_WEBVIEW_EDGE_STATIC=ON -DwxBUILD_USE_STATIC_RUNTIME=ON -DwxBUILD_SHARED=OFF .. && (
        echo Configured wxWidgets
    ) || (
        echo Failed to configure wxWidgets
        GOTO :cmderr
    )

    cmake --build . --config Debug --target install && (
        echo Installed Debug wxWidgets
    ) || (
        echo Failed to install Debug wxWidgets
        GOTO :cmderr
    )

    cmake --build . --config Release --target install && (
        echo Installed Debug wxWidgets
    ) || (
        echo Failed to install Release wxWidgets
        GOTO :cmderr
    )

    echo "Built wxWidgets!"
    cd %BUILD_DIR%
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
