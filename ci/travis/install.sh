if [ ! -f $HOME/miniconda/bin/conda ]; then
    echo "Downloading and installing miniconda"
    if [ $TRAVIS_OS_NAME = "linux" ]; then OS=Linux-x86_64; else OS=MacOSX-x86_64; fi
    wget -O miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-$OS.sh
    rm -rf $HOME/miniconda
    bash miniconda.sh -b -p $HOME/miniconda
fi
if [ ! -f $HOME/miniconda/bin/conda ]; then
    echo ERROR: conda was not installed.
    exit 1
fi
bash $HOME/miniconda/etc/profile.d/conda.sh
export PATH=$HOME/miniconda/bin/:$PATH
conda config --set always_yes yes --set changeps1 no
conda config --add channels conda-forge
conda install conda-devenv
conda update -q conda
conda info -a
conda devenv
source activate GEMS3K
mkdir build
cd build
# Configure step
cmake -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_LIBDIR=lib \
    ..
ninja install
if [ $? -eq 1 ]
then
echo "The install failed" >&2
exit 1
fi
conda list
