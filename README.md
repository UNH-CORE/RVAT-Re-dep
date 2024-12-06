# UNH-RVAT Reynolds number dependence experiment

This repository contains the processing and plotting code, as well as the
derived dataset from the UNH-RVAT Reynolds number dependence experiment
performed in Spring 2014.

## Reproducing these results

You will need
[Calkit](https://github.com/calkit/calkit)
and its dependencies installed in order to run the pipeline that generates
all of the figures.
You will also need to have a
[token set in your config](https://github.com/calkit/calkit?tab=readme-ov-file#cloud-integration)
to be able to pull the raw data (necessary for one figure.)

Clone this repository with:

```sh
calkit clone https://github.com/UNH-CORE/RVAT-Re-dep.git
```

After that, call `calkit run`.

## Reusing these materials

If you'd like to reuse some of the data loading and plotting functionality
in other Python code outside this project,
add it as a Git submodule in your own repo with:

```sh
git submodule add https://github.com/UNH-CORE/RVAT-Re-dep rvat-re-dep
```

then install the Python package in local mode with:

```sh
pip install -e ./rvat-re-dep
```

You'll then be able to run code like:

```python
import pyrvatrd

perf_curve = pyrvatrd.load_perf_curve(tow_speed=1.0)
perf_curve.plotcp()
```

See `notebook.ipynb` for more examples.
Also see
[this repo](https://github.com/petebachant/reuse-rvat-re-dep)
for an example of how this data can be reused in your own project.

## Contributing

To contribute to this repository, please create a fork, add changes on a
descriptively-named branch, and submit a pull request.

## Documents and other resources

- [Uncertainty calculations](http://nbviewer.ipython.org/github/UNH-CORE/RVAT-Re-dep/blob/master/Documents/IPython%20notebooks/uncertainty.ipynb)
- [UNH-RVAT CAD models](http://figshare.com/articles/UNH_RVAT_CAD_models/1062009)

## How to cite

Please cite

Bachant, P. and Wosnik, M. (2016) UNH-RVAT Reynolds number dependence experiment: Reduced dataset and processing code. figshare.
DOI: [10.6084/m9.figshare.1286960.v5](https://dx.doi.org/10.6084/m9.figshare.1286960.v5)

## Publications

These data were used in the following publications:

Bachant, P. and Wosnik, M. (2016) [Effects of Reynolds Number on the Energy Conversion and Near-Wake Dynamics of a High Solidity
Vertical-Axis Cross-Flow Turbine](http://doi.org/10.3390/en9020073). _Energies_, 9

Bachant, P. and Wosnik, M. (2014) [Reynolds Number Dependence of Cross-Flow Turbine Performance and Near-Wake
Characteristics](http://hdl.handle.net/10919/49210). _In Proceedings of the 2nd Marine Energy Technology Symposium METS2014_, April
15--18, Seattle, WA

## License

Code licensed under the MIT license. See `LICENSE` for details.
All other materials licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
Creative Commons Attribution 4.0 International License</a>.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
<img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by/4.0/88x31.png" />
</a>
