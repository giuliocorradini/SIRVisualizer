# SIRVisualizer

Spread model simulator and Italy COVID-19 data plotting facility.

## What is SIRVisualizer?

SIRVisualizer is a Python 3.x software that provides common facilities to visualize,
interpolate, plot, generate, simulate virus spread data; and visualize Italian COVID-19
data directly from the official repository.

## What it isn't

This program uses a very simple approach, thus must not be used if you want reliable
data or predictions.

I DON'T TAKE ANY RESPONSIBILITY for the correctness of ANY assumption, result, inference,
prediction, interpolation this program may produce.

## Usage

`data_fetcher.py` and `sir.py` provides the basic facilities to download data from the Italian repository
and to interpolate them.
You can run these files directly from the command line, or you can use the classes `SIRInterpolation` and `SIRPredict`
in your code.

`COVID19.ipynb` uses facilities from all the modules to simulate the spread of COVID-19 in Italy.
(This is done by doing massive semplification, any result mustn't be taken for accurate nor truthful).


`cellular_automata.py` aims to simulate the spread using a CA-based model, results aren't accurate yet.

`Cellular Automata Virus Spread.ipynb` uses the CA module to simulate the further spread of the COVID-19
in Italy. (This is done by doing massive semplification, any result mustn't be taken for accurate nor truthful).

## Disclaimer

This program was made for academic research sake, it's not a tool to predict the COVID-19 spread in a truthful way.
It's responsibility for anything you do with it.

### License

The program is licensed under the MIT license.