sudo apt-get install -qq -y g++ gfortran csh
sudo apt-get install -qq -y g++-multilib gcc-multilib openbabel
wget http://repo.continuum.io/miniconda/Miniconda-3.0.5-Linux-x86_64.sh
bash Miniconda-3.0.5-Linux-x86_64.sh -b
PIP_ARGS="-U"

export PATH=$HOME/miniconda/bin:$PATH
echo "PYTHON IS "
echo ${python}
conda update --yes conda
conda config --add channels http://conda.binstar.org/omnia
conda config --add channels https://conda.binstar.org/ric
conda create --yes --no-default-packages -n ${python} python=${python} --file tools/ci/requirements-conda.txt 
source activate $python
