rem @echo off

rem This is my first cmd.exe .bat script, so I'm sure it could be done much better

IF "%2"=="" GOTO invalidarg
IF not "%3"=="" GOTO invalidarg

pushd .

rem We will assume this script is in the patch directory
set PATCH_DIR=%~dp0

set ORIG_DIR=%CD%


if not exist "%1" (
  MKDIR "%1"
  if %errorlevel% neq 0 (
    echo "Failed to make dir %1"
    popd
    exit /b %errorlevel%
  ) else (
    echo "Created build dir %1"
  )
) else (
    echo "Build directory %1 already exists."
)

pushd .
cd %1
set BUILD_DIR=%CD%
popd



if not exist "%2" (
  MKDIR "%2"
  if %errorlevel% neq 0 (
    echo "Failed to make dir %2"
    popd
    exit /b %errorlevel%
  ) else (
    echo "Created install dir %2"
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

set errorlevel=0

set BOOST_TAR=boost_1_78_0.tar.gz
set BOOST_DIR=boost_1_78_0
set BOOST_BUILT_FILE=built_%BOOST_DIR%
set BOOST_REQUIRED_MD5=e193e5089060ed6ce5145c8eb05e67e3

if not exist %BOOST_BUILT_FILE% (

    if not exist %BOOST_TAR% (
        curl -L https://sourceforge.net/projects/boost/files/boost/1.78.0/boost_1_78_0.zip/download --output %BOOST_TAR%
    ) else (
        echo %BOOST_TAR% already downloaded
    )

    set "BOOST_MD5=" & for /F "skip=1 delims=" %%H in ('
       2^> nul CertUtil -hashfile "%BOOST_TAR%" MD5
    ') do if not defined BOOST_MD5 set "BOOST_MD5=%%H"

    if not "%BOOST_MD5%"=="%BOOST_REQUIRED_MD5%" (
        echo Invalid hash of boost.  Expected "%BOOST_REQUIRED_MD5%" and got "%BOOST_MD5%"
        GOTO :cmderr
    )

    if not exist %BOOST_DIR% (
        tar -xzvf boost_1_78_0.tar.gz
    ) else (
        echo "Boost was already unzipped"
    )

    if %errorlevel% neq 0 (
        echo "Failed to unzip boost"
        GOTO :cmderr
    )

    cd boost_1_78_0

    echo "Running boost bootstrap.bat"
    cmd /c .\bootstrap.bat

    if %errorlevel% neq 0 (
        echo "Failed to run bootstrap.bat"
        GOTO :cmderr
    )


    echo "Building boost"
    .\b2.exe runtime-link=static link=static threading=multi variant=release address-model=64 architecture=x86 --prefix="%MY_PREFIX%" --build-dir=win_build -j8 install


    if %errorlevel% neq 0 (
        echo "Failed to run build boost"
        GOTO :cmderr
    )

    echo "Built boost!"

    cd %BUILD_DIR%
    echo "Built boost" > %BOOST_BUILT_FILE%
) else (
    echo "Boost was already built (%BOOST_BUILT_FILE% existed)"
)

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
    
    
rem Build Wt 
set WT_TAR=wt-3.7.1.tar.gz
set WT_DIR=wt-3.7.1
set WT_BUILT_FILE=built_%WT_DIR%
set WT_REQUIRED_MD5=0f6635a213cc310fc4ec69ca48f2f46d

if not exist %WT_BUILT_FILE% (
    if not exist %WT_TAR% (
        curl -L https://github.com/emweb/wt/archive/refs/tags/3.7.1.tar.gz --output "%WT_TAR%"
    ) else (
        echo %WT_TAR% already downloaded
    )

    if %errorlevel% neq 0 (
        echo "Failed to download Wt"
        GOTO :cmderr
    )

    set "WT_MD5="
    for /f "skip=1 delims=" %%a in ('certutil -hashfile "%WT_TAR%" MD5') do if not defined WT_MD5 set "WT_MD5=%%a"

    if not "%WT_MD5%"=="%WT_REQUIRED_MD5%" (
        echo Invalid hash of Wt.  Expected "%WT_REQUIRED_MD5%" and got "%WT_MD5%"
        GOTO :cmderr
    )

    tar -xzvf %WT_TAR%
    if %errorlevel% neq 0 (
        echo "Failed to untar Wt"
        GOTO :cmderr
    )

    cd %WT_DIR%

    set WT_PATCH_FILE="%PATCH_DIR%\wt\3.7.1\NormalBuild\wt_3.7.1_git.patch"
    git apply --ignore-space-change --ignore-whitespace %WT_PATCH_FILE%
    if %errorlevel% neq 0 (
        echo "Failed to patch Wt"
        GOTO :cmderr
    )

    mkdir build_msvc2022
    cd build_msvc2022

    cmake -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="%MY_PREFIX%" -DBoost_INCLUDE_DIR="%MY_PREFIX%/include" -DBOOST_PREFIX="%MY_PREFIX%" -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="%MY_PREFIX%" -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DINSTALL_FINDWT_CMAKE_FILE=ON -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++11" -DCONFIGURATION=data/config/wt_config_electron.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR="%MY_PREFIX%/etc/wt" -DBoost_USE_STATIC_RUNTIME=ON ..
    
    if %errorlevel% neq 0 (
        echo "Failed to cmake configure Wt"
        GOTO :cmderr
    )

    cmake --build . --config Release --target install
    if %errorlevel% neq 0 (
        echo "Failed to build Release Wt"
        GOTO :cmderr
    )

    cmake --build . --config Debug --target install
    if %errorlevel% neq 0 (
        echo "Failed to build Debug Wt"
        GOTO :cmderr
    )

    echo "Built Wt!"
    cd %BUILD_DIR%
    echo "Built Wt" > %WT_BUILT_FILE%
) else (
    echo "Wt was already built (%WT_BUILT_FILE% existed)"
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
