# UNH-RVAT Reynolds number dependence experiment

This repository contains the processing and plotting code, as well as the
derived dataset from the UNH-RVAT Reynolds number dependence experiment
performed in Spring 2014.

**Note:** This repository is still in active development and is not ready
for use.

Getting started
---------------

Clone this repository with

    git clone https://github.com/UNH-CORE/RVAT-Re-dep.git

To run the processing/plotting code we recommend the
[Anaconda Python distribution](https://store.continuum.io/cshop/anaconda/)
(Python 3.4) since it includes most dependencies. The remaining
can be installed by executing

    pip install -r requirements.txt

After installing all dependencies, execute `python plot.py` to generate
figures from the experiment.

Contributing
------------

Pull requests welcome!

Documents and other resources
-----------------------------

  * [Uncertainty calculations](http://nbviewer.ipython.org/github/UNH-CORE/RVAT-Re-dep/blob/master/Documents/IPython%20notebooks/uncertainty.ipynb)
  * [UNH-RVAT CAD models](http://figshare.com/articles/UNH_RVAT_CAD_models/1062009)

## How to cite
Please cite

```bibtex
@Misc{Bachant2014-Re-dep-data,
  Title                    = {UNH-RVAT Reynolds number dependence experiment: Reduced dataset and processing code},
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
@INPROCEEDINGS{Bachant2014-METS
  Author                   = {Bachant, P. and Wosnik, M.},
  Title                    = {Reynolds Number Dependence of Cross-Flow Turbine Performance and Near-Wake Characteristics},
  Booktitle                = {Proceedings of the 2nd Marine Energy Technology Symposium METS2014},
  Month                    = {April},
  Year                     = {2014},
  Address                  = {Seattle, WA}
  Url                      = {http://hdl.handle.net/10919/49210}
}
```

License
-------
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
<img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by/4.0/88x31.png" />
</a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
Creative Commons Attribution 4.0 International License</a>.
