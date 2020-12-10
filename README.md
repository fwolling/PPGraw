# PPGraw
## The Quest for Raw Signals:<br>A Quality Review of Publicly Available Photoplethysmography Datasets

In this GitHub repository, we present an analytical tool for the quality review of raw photoplethysmography (PPG) signals, based on 7 multi-varied decision metrics. It has been applied in the review of 10 publicly available photoplethysmography datasets, referred below in [Citation](#citation). Although all [evaluated datasets](#evaluated-datasets) were advertised to contain raw signals, the characteristics of the PPG data look quite diverse. Our developed tool enables to automatically analyze the suitability and applicability of datasets and helps to identify preprocessed and filtered signals with a limited evidence. The [raw reference data](#reference-data), recorded with the MAX86140EVSYS# evaluation system, as well as the implemented [Python tool](/python/), based on the presented 7 decision metrics, are available for [download](#download), to support the reproducibility and the review of new datasets.

### Download
This GitHub repository provides the developed analytical tool *PPGraw*.
The raw photoplethysmography reference signals can be downloaded via the following link:
https://ubicomp.eti.uni-siegen.de/home/datasets/data20/index.html.en

### Citation
"[The Quest for Raw Signals: A Quality Review of Publicly Available Photoplethysmography Datasets](https://ubicomp.eti.uni-siegen.de/home/datasets/data20/index.html.en)", <a href="https://ubicomp.eti.uni-siegen.de/home/team/fwolling.html.en" target="_blank">Florian Wolling</a> and <a href="https://ubicomp.eti.uni-siegen.de/home/team/kristof.html.en" target="_blank">Kristof Van Laerhoven</a>. In *DATA'20: Proceedings of the 3rd Workshop on Data Acquisition To Analysis, DATA 2020, Virtual Event, Japan, November 2020*, ACM, 2020. <a href="https://doi.org/10.1145/3419016.3431485" target="_blank">https://doi.org/10.1145/3419016.3431485</a>

### Disclaimer
You may use the source code of the developed analytical tool *PPGraw* for scientific, non-commercial purposes, provided that you give credit to the owners when publishing any work based on it. We would also be very interested to hear back from you if you use our tool or metrics in any way and are happy to answer any questions or address any remarks related to it.

### Presentation Video
<a href="https://www.youtube.com/watch?v=RshKMVtH7P0" target="_blank"><img src="https://raw.githubusercontent.com/fwolling/PPGraw/main/fig/youtube.png" alt="DATA'20 - The Quest for Raw Signals - A Quality Review of Photoplethysmography Datasets" width="600" style="float: center;" /></a>

<br>

### Reference Data
**REF:** MAX86140EVSYS# [<a href="https://ubicomp.eti.uni-siegen.de/home/datasets/data20/index.html.en" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/00_max86140.png" alt="Reference Data: MAX86140EVSYS#" style="float: left; margin-right: 10px;" />

### Evaluated Datasets

The following figures show excerpts from the reviewed datasets: Short close-up of few pulses on the left, a 30-second window in the middle, and its respective frequency spectrum (FFT) on the right. Note that the PPG-BP dataset (S05) contains only snippets of 2.1 s length. Frequency bands: very low frequency (VLF, < 0.167 Hz, red), low frequency (LF, 0.167 to 0.667 Hz, orange), and intermediate frequency (IF, 0.5 to 3.0 Hz , green) while the high frequency (HF, > 3.0 Hz) noise and harmonics are clipped.

**S01:** MAXREFDES100# [<a href="#ref_s01">Biagetti et al.</a>, <a href="https://www.sciencedirect.com/science/article/pii/S2352340919314003" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/01_maxrefdes100.png" alt="S01: MAXREFDES100#" style="float: left; margin-right: 10px;" />
**S02:** PPG-DaLiA [<a href="#ref_s02">Reiss et al.</a>, <a href="https://ubicomp.eti.uni-siegen.de/home/datasets/sensors19/index.html.en" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/02_dalia.png" alt="S02: PPG-DaLiA" style="float: left; margin-right: 10px;" />
**S03:** WESAD [<a href="#ref_s03">Schmidt et al.</a>, <a href="https://ubicomp.eti.uni-siegen.de/home/datasets/icmi18/index.html.en" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/03_wesad.png" alt="S03: WESAD" style="float: left; margin-right: 10px;" />
**S04:** BloodLossSVM [<a href="#ref_s04">Reljin et al.</a>, <a href="https://figshare.com/articles/NR_bloodlosssvm_zip/5594644" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/04_nrblsvm.png" alt="S04: BloodLossSVM" style="float: left; margin-right: 10px;" />
**S05:** PPG-BP [<a href="#ref_s05">Liang et al.</a>, <a href="https://figshare.com/articles/PPG-BP_Database_zip/5459299" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/05_ppgbp.png" alt="S05: PPG-BP" style="float: left; margin-right: 10px;" />
**S06:** BIDMC [<a href="#ref_s06">Pimentel et al.</a>, <a href="https://physionet.org/content/bidmc/1.0.0/" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/06_bidmc.png" alt="S06: BIDMC" style="float: left; margin-right: 10px;" />
**S07:** Wrist PPG During Exercise [<a href="#ref_s07">Jarchi et al.</a>, <a href="https://physionet.org/content/wrist/1.0.0/" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/07_wpde.png" alt="S07: Wrist PPG During Exercise" style="float: left; margin-right: 10px;" />
**S08:** Cuff-Less Blood Pressure Estimation [<a href="#ref_s08">Kachuee et al.</a>, <a href="https://archive.ics.uci.edu/ml/datasets/Cuff-Less+Blood+Pressure+Estimation" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/08_clbpe.png" alt="S08: Cuff-Less Blood Pressure Estimation" style="float: left; margin-right: 10px;" />
**S09:** IEEE SPC 2015 (TROIKA) [<a href="#ref_s09">Zhang et al.</a>, <a href="https://sites.google.com/site/researchbyzhang/ieeespcup2015" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/09_spc2015.png" alt="S09: IEEE SPC 2015 (TROIKA)" style="float: left; margin-right: 10px;" />
**S10:** IEEE SPC 2013 [<a href="#ref_s10">Karlen et al.</a>, <a href="http://www.capnobase.org/index.php?id=857" target="_blank">data</a>]
<img src="https://github.com/fwolling/PPGraw/blob/main/fig/10_spc2013.png" alt="S10: IEEE SPC 2013" style="float: left; margin-right: 10px;" />

<br>

### References

<a id="ref_s01">**S01:**</a> Giorgio Biagetti, Paolo Crippa, Laura Falaschetti, Leonardo Saraceni, AndreaTiranti, and Claudio Turchetti. 2020. "Dataset from PPG wireless sensor for activity monitoring". *Data in brief 29 (2020), 105044*. https://doi.org/10.1016/j.dib.2019.105044

<a id="ref_s02">**S02:**</a> Attila Reiss, Ina Indlekofer, Philip Schmidt, and Kristof Van Laerhoven. 2019. "DeepPPG: Large-Scale Heart Rate Estimation with Convolutional Neural Networks". *Sensors (Basel, Switzerland) 19, 14 (2019)*. https://doi.org/10.3390/s19143079

<a id="ref_s03">**S03:**</a> Philip Schmidt, Attila Reiss, Robert Duerichen, Claus Marberger, and Kristof Van Laerhoven. 2018. "Introducing WESAD, a Multimodal Dataset for Wearable Stressand Affect Detection". https://doi.org/10.1145/3242969.3242985

<a id="ref_s04">**S04:**</a> Natasa Reljin, Gary Zimmer, Yelena Malyuta, Kirk Shelley, Yitzhak Mendel-son, David  J. Blehar, Chad E. Darling, and Ki H. Chon. 2018. "Using support vector machines on photoplethysmographic signals to discriminate between hypovolemia and euvolemia". *PloS one 13, 3 (2018), e0195087*. https://doi.org/10.1371/journal.pone.0195087

<a id="ref_s05">**S05:**</a> Yongbo Liang, Zhencheng Chen, Guiyong Liu, and Mohamed Elgendi. 2018. "A new, short-recorded photoplethysmogram dataset for blood pressure monitoring in China". *Scientific Data 5, 1 (2018), 180020*. https://doi.org/10.1038/sdata.2018.20

<a id="ref_s06">**S06:**</a> Marco A. F. Pimentel, Alistair E. W. Johnson, Peter H. Charlton, Drew Birrenkott, Peter J. Watkinson, Lionel Tarassenko, and David A. Clifton. 2017. "Toward a Robust Estimation of Respiratory Rate From Pulse Oximeters". *IEEE transactions on bio-medical engineering 64, 8 (2017), 1914–1923*.  https://doi.org/10.1109/TBME.2016.2613124

<a id="ref_s07">**S07:**</a> Delaram Jarchi and Alexander Casson. 2017. "Description of a Database Containing Wrist PPG Signals Recorded during Physical Exercise with Both Accelerometer and Gyroscope Measures of Motion". *Data 2, 1 (2017), 1*. https://doi.org/10.3390/data2010001

<a id="ref_s08">**S08:**</a> Mohamad Kachuee, Mohammad Mahdi Kiani, Hoda Mohammadzade, and Mahdi Shabany. 2015. "Cuff-less high-accuracy calibration-free blood pressure estimation using pulse transit time". *(2015), 1006–1009*. https://doi.org/10.1109/ISCAS.2015.7168806

<a id="ref_s09">**S09:**</a> Zhilin Zhang, Zhouyue Pi, and Benyuan Liu. 2015. "TROIKA: A General Framework for Heart Rate Monitoring Using Wrist-Type Photoplethysmographic Signals During Intensive Physical Exercise". *IEEE transactions on bio-medical engineering 62, 2 (2015), 522–531*. https://doi.org/10.1109/TBME.2014.2359372

<a id="ref_s10">**S10:**</a> Walter Karlen, Srinivas Raman, J. Mark Ansermino, and Guy A. Dumont. 2013. "Multiparameter Respiratory Rate Estimation from the Photoplethysmogram". *IEEE transactions on bio-medical engineering 60, 7 (2013), 1946–1953*. https://doi.org/10.1109/TBME.2013.2246160
