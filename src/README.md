# PPGraw

*PPGraw* is a Python tool developed for the analytical quality review of supposedly raw, unfiltered photoplethysmography signals. This directory contains an excerpt from the raw, unfiltered reference data [ref_sample.p](ref_sample.p), the Python tool itself [PPGraw.py](PPGraw.py), a Jupyter notebook [ref_show.ipynb](ref_show.ipynb) to visualize the reference signal, and a Jupyter notebook [ref_review.ipynb](ref_review.ipynb) that runs a quality review on the reference data and details the 7 multi-varied decision metrics.

* **ref_sample.p** is a Python Pickle file containing a 3-minute excerpt of the raw, unfiltered photoplethysmography reference data
* **PPGraw.py** contains the Python class of the analytical tool *PPGraw* for the automated quality review of raw photoplethysmography signals
* **ref_show.ipynb** is a Jupyter notebook that visualizes the provided reference data
* **ref_review.ipynb** is a Jupyter notebook that runs the quality review on the reference data and details the 7 multi-varied decision metrics

### Download
This GitHub repository provides the developed analytical tool *PPGraw*.
The raw photoplethysmography reference signals can be downloaded via the following link:
https://ubicomp.eti.uni-siegen.de/home/datasets/data20/index.html.en

### Citation
"[The Quest for Raw Signals: A Quality Review of Publicly Available Photoplethysmography Datasets](https://ubicomp.eti.uni-siegen.de/home/datasets/data20/index.html.en)", <a href="https://ubicomp.eti.uni-siegen.de/home/team/fwolling.html.en" target="_blank">Florian Wolling</a> and <a href="https://ubicomp.eti.uni-siegen.de/home/team/kristof.html.en" target="_blank">Kristof Van Laerhoven</a>. In *DATA'20: Proceedings of the 3rd Workshop on Data Acquisition To Analysis, DATA 2020, Virtual Event, Japan, November 2020*, ACM, 2020. <a href="https://doi.org/10.1145/3419016.3431485" target="_blank">https://doi.org/10.1145/3419016.3431485</a>

### Disclaimer
You may use the source code of the developed analytical tool *PPGraw* for scientific, non-commercial purposes, provided that you give credit to the owners when publishing any work based on it. We would also be very interested to hear back from you if you use our tool or metrics in any way and are happy to answer any questions or address any remarks related to it.
