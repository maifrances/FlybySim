# FlybySim 
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Overview
Python script to simulate the field of view of the Thermal Infrared Imager (TIRI) onboard the Modular Infrared Molecules and Ices Sensor (MIRMIS) during the comet encounter.


## Setup
### Requirements
This code was developed and tested in [Python 3.10.9](https://docs.python.org/release/3.10.9/). The comet encounter is simulated using [PyQt v5](https://www.riverbankcomputing.com/static/Docs/PyQt5/) bindings as the graphical interface backend engine. The mesh handling (download and display) sections are heavily inspired by [Part 12 of Thomas Albin's space science tutorials](https://github.com/ThomasAlbin/SpaceScienceTutorial)

The Python packages required for running all functions are:

* [Argparse](https://docs.python.org/3/library/argparse.html)
* [Astropy](http://www.astropy.org/)
* [imageio](https://imageio.readthedocs.io/en/stable/)
* [Matplotlib](https://matplotlib.org/)
* [Numpy](http://www.numpy.org/)
* [Pandas](https://pandas.pydata.org/)
* [Pathlib](https://docs.python.org/3/library/pathlib.html)
* [Seaborn](https://seaborn.pydata.org/)
* [Shapely](https://shapely.readthedocs.io/en/stable/)
* [Urllib](https://docs.python.org/3/library/urllib.html)
* [Visvis](https://github.com/almarklein/visvis)

### Usage
```bash
$ FlybySim --help
usage: FlybySim [-h] [-rest {tiri,sc}] [-ff] [-V] [-Vc] [-ca] [-o] [-readfreq] [-fps] [--nosave] [--noann] [--saveA] [--saveraw] {path to mesh}

positional arguments:
  {path to mesh}
  path to mesh     path to comet shape model .OBJ file stored on local disk.

options:
  -h, --help       show this help message and exit.
  -rest {tiri,sc}  Rest frame of simulation view (default: tiri).
  -ff      		   Fast forward non-datacube simulation frames by this amount (default: 20).
  -V               Encounter velocity vector in km/s (default: 70).
  -Vc              Comet velocity in km/s (default: 15).
  -ca              Closest approach distance in km (default: 1000).
  -o       		   Path to local directory to save output. If not provided,defaults to current working directory (default:current working directory).
  -readfreq        Readout frequency of the TIRI detector in Hz (default:0.05).
  -fps             Number of frames per second in output gif (default: 20).
  --nosave         Don't save model gif and table output (default: True).
  --noann          Don't annotate simulation frames with TIRI filters, etc (default: True).
  --saveA          Don't save annotated simulation frames as .npy arrays (default: False).
  --saveraw        Save raw simulation frames to local disk (default: False).
```


## Documentation
Documentation is included in the [docstrings](https://www.python.org/dev/peps/pep-0257/) (strings enclosed by triple quotes `"""..."""`) within the source code.
Or for those with access, in the TIRI Flyby Model Description (COMET-OXF-MIR-TN-130).


## Outputs
Simulation data is output in two formats:
1. .csv - an ascii table with raw data e.g:
2. .gif - An animation of the annotated simulation frames with a user defined frame rate:

| <img src="/output_examples/flyby-SC-0.05.gif" width="300"> | <img src="/output_examples/flyby-TIRI-0.05.gif" width="300"> |
| -- | -- |
| *Spacecraft rest frame* | *TIRI rest frame* |








