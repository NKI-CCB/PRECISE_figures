# Create environment
conda create -n precise_figures python=3.6
conda activate precise_figures

# Activate Jupyter notebook
pip install ipykernel
python -m ipykernel install --user --name precise_figures --display-name "Python (precise_figures)"
conda install -c r rpy2
conda install tzlocal


# For PRECISE
pip install numpy
pip install pandas
pip install scikit-learn
pip install joblib

# For notebooks
pip install matplotlib
pip install seaborn
conda install xarray
conda install netcdf4

# Install PRECISE
cd ../precise/
python setup.py install

# Install edgeR
python install_edgeR.py