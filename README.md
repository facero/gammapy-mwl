Caveat: this is a work in progress. More details and examples will be added here

Gammapy-mwl 
=======
A Python tool for the analysis of Gamma-ray and MWL data in a unified framework, using physically-motivated spectral models.

## What this package can do and can not do 

TBC



## Installation and Set-up

These instructions assume that you have previously installed a version of `conda` or `mamba` on your machine.
To set-up the work environment with conda::

## Create the environment 
```
  conda create -n gammapy-mwl-0.1
  conda activate gammapy-mwl-0.1
  conda install -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge ciao sherpa
  conda install -c conda-forge gammapy
```

## Clone the repository
```
git clone https://github.com/yourusername/gammapy-mwl.git
cd gammapy-mwl
```
## Install the package (editable mode for development)
```
pip install -e .
```

Additionally, if you wish to combine the absorption models provided in Sherpa with the physical models provided by Naima, you have to `install Naima <https://naima.readthedocs.io/en/latest/installation.html>`_.
TBD: show how this is done!



Citing
+++++++++++++++++++++++++++++++++++++++++++++

A software description is provided in the following publication: TBD


If you use gammapyXray for work/research presented in a publication (whether directly, or as a dependency to another package), we ask that you please cite it using the following links

???

We encourage you to also include citations to the paper in the main text
wherever appropriate, using the recommended BibTeX entry:


Licence
+++++++
This folder is licensed under a 3-clause BSD style license - see the
`LICENSE.rst <https://github.com/gammapy/gammapy/blob/master/LICENSE.rst>`_ file.

.. image:: https://anaconda.org/conda-forge/gammapy/badges/license.svg
    :target: TBD
    :alt: Licence









--------OLDER NOTES FOR LATER--------------

# gammapy-mwl
A repository for tools to convert various MWL data to gammapy

## Repositories with MWL tools using gammapy

###  A MWL workflow from optical to GeV 
- from @mireianievas : https://github.com/mireianievas/gammapy_mwl_workflow

### X-ray OGIP data manipulation 

- Latest fork of gammapyXray : https://github.com/mireianievas/gammapy-ogip-spectra

## Test data
- HEACIT curated list of X-ray test data : https://github.com/HEACIT/curated-test-data

