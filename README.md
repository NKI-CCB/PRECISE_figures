# PRECISE_figures
Jupyter notebooks and code to reproduce figures.

## TO DO before launching the notebooks

- Install all the dependencies
sh setup.sh

- Install GSEA package (Java implementation) [Subramanian et al 2005] available at the following address:
<br/>
http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/software/gsea-3.0.jar
This has to be put in the root folder.

- Clone ComBat version using Brent Pedersen implementation donwloadable as this location:
<br/>
git clone https://github.com/brentp/combat.py
<br/>
mv combat.py combat

- Install statannot using Marc Weber implementation:
conda activate precise_figures
pip install git+https://github.com/webermarcolivier/statannot
conda deactivate
