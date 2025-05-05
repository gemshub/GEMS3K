if "%APPVEYOR_BUILD_WORKER_IMAGE%"=="Visual Studio 2017" (
    call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars64.bat"
)
if "%APPVEYOR_BUILD_WORKER_IMAGE%"=="Visual Studio 2019" (
    call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
)

echo "Configuring..."
cmake -G"Visual Studio 16 2019" -A x64 -S . -DCMAKE_INSTALL_PREFIX:PATH="%CONDA_PREFIX%\Library" -DBUILD_SHARED_LIBS=OFF -B build
echo "Building..."
cmake --build build --config %CONFIGURATION% --target install
if errorlevel 1 exit 1
