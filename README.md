# UNH-RVAT Reynolds number dependence experiment

This repository contains the processing and plotting code, as well as the 
derived dataset from the UNH-RVAT Reynolds number dependence experiment
performed in Spring 2014.

**Note:** This repository is still in active development and is not ready
for use.

Documents
---------

  * [Uncertainty calculations](http://nbviewer.ipython.org/github/UNH-CORE/RVAT-Re-dep/blob/master/Documents/IPython%20notebooks/uncertainty.ipynb)

Download/usage
--------------

To download, use the git clone URL to the right. For example, in a terminal

    git clone https://github.com/UNH-CORE/RVAT-Re-dep.git

See `run.py` for functions to process and visualize the data.

To contribute back to the main repository, use GitHub's fork/pull mechanism.

### Dependencies

  * NumPy
  * SciPy
  * matplotlib
  * pandas
  * [PXL](https://github.com/petebachant/PXL)

Dependencies can be installed automatically by running

    pip install -r requirements.txt

## How to cite
Please cite 

```bibtex
@Misc{Bachant2014_Re-dep-data,
  Title                    = {UNH-RVAT Reynolds number dependence experiment: Reduced dataset and processing code, rev1},
  Author                   = {Peter Bachant and Martin Wosnik},
  HowPublished             = {fig\textbf{share}. http://dx.doi.org/10.6084/m9.figshare.1286960,
  Month                    = {January},
  Year                     = {2015},
  Doi                      = {10.6084/m9.figshare.1286960},
  Url                      = {http://dx.doi.org/10.6084/m9.figshare.1286960}
}
```

Publications
------------
These data were used in the following publications:

```bibtex
@INPROCEEDINGS{Bachant2014
  author = {Bachant, P. and Wosnik, M.},
  title = {Reynolds Number Dependence of Cross-Flow Turbine Performance and Near-Wake Characteristics},
  booktitle = {Proceedings of the 2nd Marine Energy Technology Symposium METS2014},
  year = {2014},
}
```

Other resources
---------------

Turbine CAD (STEP) files are available at http://figshare.com/articles/UNH_RVAT_CAD_models/1062009

License
-------
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
<img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by/4.0/88x31.png" />
</a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
Creative Commons Attribution 4.0 International License</a>.
