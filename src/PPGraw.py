# Florian Wolling
# Ubiquitous Computing Lab
# https://ubicomp.eti.uni-siegen.de
# University of Siegen, Germany
#
# "The Quest for Raw Signals: A Quality Review of Publicly Available Photoplethysmography Datasets"
# by Florian Wolling and Kristof Van Laerhoven. In DATA'20: Proceedings of the 3rd Workshop on Data Acquisition To
# Analysis, DATA 2020, Virtual Event, Japan, November 2020, ACM, 2020. https://doi.org/10.1145/3419016.3431485
#
# Version 1.0 (November 2020)


import pickle
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
from tabulate import tabulate

__author__ = "Florian Wolling, https://github.com/fwolling/"
__version__ = "1.0"


class PPGraw:
    """PPGraw - analytical tools for the quality review of raw photoplethysmography datasets

    PPGraw provides analytical tools for the automatic quality review of photoplethysmography datasets, containing
    preferably raw, unfiltered, and non-preprocessed signals.

    The tools have been developed in context of the following publication for the DATA'20 workshop:
    "The Quest for Raw Signals: A Quality Review of Publicly Available Photoplethysmography Datasets"
    by Florian Wolling and Kristof Van Laerhoven. In DATA'20: Proceedings of the 3rd Workshop on Data Acquisition To
    Analysis, DATA 2020, Virtual Event, Japan, November 2020, ACM, 2020. https://doi.org/10.1145/3419016.3431485

    Information: https://ubicomp.eti.uni-siegen.de/home/datasets/data20/index.html.en
    GitHub repository: https://github.com/fwolling/PPGraw

    Version 1.0 (November 2020)

    Attributes
    ----------
    time : None or list of floats, optional, default = None
        The vector that contains the timestamps.
    signal : list of floats, required
        The vector that contains the signal samples.
    fs : float
        The desired sampling frequency in Hz.

    Methods
    -------
    review_timebase
        ...
    review_amplitude
        ...
    review_granularity
        ...
    review_normalization
        ...
    review_flipped
        ...
    review_frequency
        ...
    review_artifacts
        ...
    review
        ...

    Notes
    -----


    References
    ----------
    [1] "The Quest for Raw Signals: A Quality Review of Publicly Available Photoplethysmography Datasets"
    by Florian Wolling and Kristof Van Laerhoven. In DATA'20: Proceedings of the 3rd Workshop on Data Acquisition To
    Analysis, DATA 2020, Virtual Event, Japan, November 2020, ACM, 2020. https://doi.org/10.1145/3419016.3431485
    """

    # Parameters

    print_table_format = None

    # Error margin in % for detection of zero-centered signals.
    ampl_zc_err_margin = 1.5

    # Parameters for the flip detection.
    flip_bp_order = 3  # Order of the band-pass filter.
    flip_bp_fc_low = 0.4  # Lower corner frequency.
    flip_bp_fc_high = 15.0  # Upper corner frequency.

    # Frequency Analysis
    frq_win_len = 60.0  # Length in s of Hanning window for PSD calculation.

    def __init__(self, time=None, signal=None, fs=None,
                 table_format="pretty",
                 debug=False):
        
        # Store signal data in instance variables if provided.
        if signal is not None:
            self.signal = signal

            # Store time stamps if provided.
            if time is not None:
                self.time = time
            # Store sampling rate if provided.
            if fs is not None:
                self.fs = fs

            # Check whether either time base or sampling rate are provided.
            if time is None and fs is None:
                raise ValueError("Either time base 'time' or sampling rate 'fs' must be passed to the instance.")

            # Forward arguments to the instance's parameters.
            self.print_table_format = table_format

    def lowpass(self, data, fs, order, fc,
                debug=False):
        nyq = 0.5 * fs  # Calculate the Nyquist frequency.
        cut = fc / nyq  # Calculate the cutoff frequency (-3 dB).
        lp_b, lp_a = sig.butter(order, cut, btype="lowpass")  # Design and apply the low-pass filter.
        lp_data = list(sig.filtfilt(lp_b, lp_a, data))  # Apply forward-backward filter with linear phase.
        return lp_data

    def highpass(self, data, fs, order, fc,
                 debug=False):
        nyq = 0.5 * fs  # Calculate the Nyquist frequency.
        cut = fc / nyq  # Calculate the cutoff frequency (-3 dB).
        hp_b, hp_a = sig.butter(order, cut, btype="highpass")  # Design and apply the high-pass filter.
        hp_data = list(sig.filtfilt(hp_b, hp_a, data))  # Apply forward-backward filter with linear phase.
        return hp_data

    def bandpass(self, data, fs, order, fc_low, fc_hig,
                 debug=False):
        nyq = 0.5 * fs  # Calculate the Nyquist frequency.
        cut_low = fc_low / nyq  # Calculate the lower cutoff frequency (-3 dB).
        cut_hig = fc_hig / nyq  # Calculate the upper cutoff frequency (-3 dB).
        bp_b, bp_a = sig.butter(order, (cut_low, cut_hig), btype="bandpass")  # Design and apply the band-pass filter.
        bp_data = list(sig.filtfilt(bp_b, bp_a, data))  # Apply forward-backward filter with linear phase.
        return bp_data

    def psd(self, data, fs, n_seg, n_fft, win="hanning", mode="densitiy",
            debug=False):
        """ Welch's Method for Power Spectrum Estimation

        :param data:
        :param fs:
        :param n_seg:
        :param n_fft:
        :param win:
        :param mode: Power spectral density "density" where Pxx has units of V**2/Hz or power spectrum "spectrum" where Pxx has units of V**2.
        :return:
        """

        # default: noverlap=n_seg//2
        bin, psd = sig.welch(data, fs, window=win, nperseg=n_seg, nfft=n_fft, return_onesided=True, scaling=mode)
        return bin, psd

    def detrend(self, data,
                debug=False):
        return data - np.mean(data)  # Zero-Centered

    def norm01(self, data,
                debug=False):
        return (data - np.min(data)) / (np.max(data) - np.min(data))  # [0,1]-Normalization

    def norm11(self, data,
                debug=False):
        # Calculate error margin based on the given percentage.
        err_margin = (np.max(data) - np.min(data)) * (self.ampl_zc_err_margin/100.0)
        # [0,1]-Normalization
        if -err_margin <= np.mean(data) <= err_margin:
            if np.max(data) > abs(np.min(data)):
                return data / np.max(data)
            else:
                return data / np.min(data)
        else:
            return 2 * (data - np.min(data)) / (np.max(data) - np.min(data)) - 1.0

    def flip(self, data,
             debug=False):
        return [-d for d in data]  # Flipping

    def review_timebase(self, time=None, signal=None, fs=None,
                       factor=0,  # 0: seconds, 3: milliseconds, 6: microseconds
                       debug=False):
        """Analyse the provided time base regarding real sampling rate and jitter.

        :param time:
        :param signal:
        :param fs:
        :param factor:
        :param debug:

        :return: a dict containing:
            "time_base": {False, True}, if time base provided;
            "time_off": float of time in seconds offset of the time base;
            "time_len": length of the recording;
            "fs": float, desired sampling rate in Hz;
            "fs_real": float, real average sampling rate in Hz;
            "fs_std": float, jitter of the real sampling rate in seconds s;
        """

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        # Check whether either time base or sampling rate are provided.
        if time is None and fs is None:
            raise ValueError("Either time base 'time' or sampling rate 'fs' must be provided.")

        # Check whether a time base is provided.
        time_off = 0.0
        if time is None:
            # Time base not provided.
            time_base = False
            # Generate the time base according to the sampling rate.
            time = np.arange(0.0, 1.0/fs * len(signal), 1.0/fs)
        else:
            # Time base provided.
            time_base = True
            # Adjust time base according to unit factor.
            unit_fact = pow(10, factor)
            time = [t/unit_fact for t in time]
            # Align time base to zero.
            if time[0] != 0.0:
                time_off = time[0]
                time = [t - time[0] for t in time]

        # Determine length of the dataset in seconds.
        time_len = time[-1] - time[0]

        # Calculate the steps between the particular timestamps.
        t_dif = [time[i] - time[i-1] for i in range(1, len(time))]
        # Calculate mean and standard deviation of the sampling rate.
        fs_mean = 1.0/np.mean(t_dif)
        fs_std = np.std(t_dif) if np.std(t_dif) > 0.0 else None

        # Adapt time base from milli or even micro seconds.
        if abs(fs_mean) < 1.0:
            # Reversely call the analysis with different unit factor.
            factor += 3
            res = self.review_timebase(time, signal, fs, factor=factor, debug=False)
            time_off = res["time_off"]
            time_len = res["time_len"]
            fs_mean = res["fs_real"]
            fs_std = res["fs_std"]

        # Store measures of real sampling rate.
        fs_real = fs_mean if time_base is not None else None
        fs_std = fs_std if time_base is not None and fs_std is not None else None

        # Print the debug information.
        if debug:
            print(" TIME BASE")
            print(tabulate([
                ["time base provided", str("no" if time_base is None else "yes")],
                ["time base", str("s" if factor == 0 else "ms" if factor == 3 else "$\mu$s" if factor == 9 else "/")],
                ["length", str(round(time_len, 1))+" s"],
                ["offset", str(time_off)+" s"],
                ["desired sampling rate", str(fs) + " Hz"],
                ["real sampling rate", ((str(fs_real) + " Hz") if time_base is not None else str(None))],
                ["jitter", (str(fs_std) + str(" s") if time_base is not None and fs_std is not None else str(None))],
            ], tablefmt="pretty", colalign=("left", "right")))

        # Return a dictionary containing the collected results.
        return {
            "time_base":                time_base,
            "time_len":                 time_len,
            "time_off":                 time_off,
            "fs":                       fs,
            "fs_real":                  fs_real,
            "fs_std":                   fs_std,
        }

    def review_amplitude(self, time=None, signal=None, fs=None,
                        debug=False):

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        # Calculate minimum, maximum, span, median, and mean of the signal.
        ampl_min = np.min(signal)
        ampl_max = np.max(signal)
        ampl_span = np.max(signal) - np.min(signal)
        ampl_med = np.median(signal)
        ampl_mean = np.mean(signal)

        # Print the debug information.
        if debug:
            print(" AMPLITUDE")
            print(tabulate([
                ["minimum", np.min(signal)],
                ["maximum", np.max(signal)],
                ["span", (np.max(signal) - np.min(signal))],
                ["median", np.median(signal)],
                ["mean", np.mean(signal)],
            ], tablefmt="pretty", colalign=("left", "right")))

        # Return a dictionary containing the collected results.
        return {
            "ampl_min":                 ampl_min,
            "ampl_max":                 ampl_max,
            "ampl_span":                ampl_span,
            "ampl_med":                 ampl_med,
            "ampl_mean":                ampl_mean,
        }

    def review_granularity(self, time=None, signal=None, fs=None,
                          debug=False):

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        # Determine granularity respectively amplitude resolution.
        signal_sort = np.unique(signal)  # Sort the values and remove duplicates.
        num_val = len(signal_sort)  # Determine number of amplitude values.
        signal_step = [signal_sort[i]-signal_sort[i-1] for i in
                       range(1, len(signal_sort))]  # Determine step sizes.
        granularity = np.min(signal_step)  # Determine granularity.

        # Print the debug information.
        if debug:
            print(" GRANULARITY")
            print(tabulate([
                ["number of values", num_val],
                ["granularity", (int(1) if granularity == 1.0 else granularity)],
            ], tablefmt="pretty", colalign=("left", "right")))

        # Return a dictionary containing the collected results.
        return {
            "ampl_num":                 num_val,
            "ampl_gran":                granularity,
        }

    def review_normalization(self, time=None, signal=None, fs=None,
                            debug=False):

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        # Calculate error margin based on the given percentage.
        err_margin = (np.max(signal) - np.min(signal)) * (self.ampl_zc_err_margin/100.0)

        # Zero-Centered
        ampl_zc = True if -err_margin <= np.mean(signal) <= err_margin else False

        # [0,1]-Normalized
        ampl_norm01 = True if np.min(signal) >= 0.0 and np.max(signal) <= 1.0 else "quasi" \
            if np.min(signal) >= 0.0 - err_margin and np.max(signal) <= 1.0 + err_margin else False

        # [-1,1]-Normalized
        ampl_norm11 = True if (-1.0 <= np.min(signal) < 0.0) and np.max(signal) <= 1.0 else "quasi" \
            if (-1.0 - err_margin <= np.min(signal) < 0.0) and np.max(signal) <= 1.0 + err_margin else False

        # Print the debug information.
        if debug:
            print(" NORMALIZATION")
            print(tabulate([
                ["zero-centered", str("yes" if ampl_zc else "no")],
                ["[0,1]-normalized", str("yes" if ampl_norm01 is True else
                                         "quasi (\u03B5=" + str(round(err_margin, 3)) + ")"
                                         if ampl_norm01 == "quasi" else "no")],
                ["[-1,1]-normalized", str("yes" if ampl_norm11 is True else
                                          "quasi (\u03B5=" + str(round(err_margin, 3)) + ")"
                                          if ampl_norm11 == "quasi" else "no")],
            ], tablefmt="pretty", colalign=("left", "right")))

        # Return a dictionary containing the collected results.
        return {
            "ampl_zc":                  ampl_zc,
            "ampl_norm01":              ampl_norm01,
            "ampl_norm11":              ampl_norm11,
        }

    def review_flip(self, time=None, signal=None, fs=None,
                    plot=False,
                    debug=False):

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        # Apply band-pass filter before the analysis.
        signal_bp = signal - np.mean(signal)
        signal_bp = self.bandpass(signal_bp, fs=fs, order=self.flip_bp_order,
                                  fc_low=self.flip_bp_fc_low, fc_hig=self.flip_bp_fc_high)

        # Check if the signal has been flipped.
        granularity = self.review_granularity(time, signal, fs, debug=False)["ampl_gran"]

        # Determine the number of bins for the histogram.
        bins_n = int(round((np.max(signal_bp) - np.min(signal_bp)) / granularity))
        # Limit the number of bins to 500.
        bins_n = (501 if bins_n > 500 else bins_n)
        bins_range = (int(np.min(signal_bp)), int(np.max(signal_bp)))
        hist_ys, hist_bins = np.histogram(signal_bp, bins=bins_n, range=bins_range)
        hist_max = (hist_bins[int(np.argmax(hist_ys))] + hist_bins[int(np.argmax(hist_ys)) + 1]) / 2

        # Plot the histogram.
        if plot:
            plt.figure()
            plt.title("Amplitude Flipping")
            plt.hist(signal_bp, bins=bins_n, range=bins_range)
            plt.xlabel("amplitude")
            plt.ylabel("count")

        # Determine the amplitudes' center of mass.
        ampl_flip = True if hist_max > 0 else False

        # Print the debug information.
        if debug:
            print(" AMPLITUDE FLIPPING")
            print(tabulate([
                ["center of mass", str(hist_max)],
                ["flipped", str("yes" if ampl_flip else "no")],
            ], tablefmt="pretty", colalign=("left", "right")))

        # Return a dictionary containing the collected results.
        return {
            "ampl_flip":                ampl_flip,
            "ampl_flip_hist_max":       hist_max,
        }

    def review_frequency(self, time=None, signal=None, fs=None,
                         plot=False,
                         debug=False):

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        nperseg = int(fs * self.frq_win_len)

        # fft_bin, fft_mag = fourier.fft(signal, fs=fs, win="hanning", scale="mag")

        psd_bin, psd_mag = self.psd(signal, fs=fs, n_seg=nperseg, n_fft=nperseg,
                                    win="hanning", mode="density")  # Calculate Power Spectral Density in V**2/Hz
        # psd_bin, psd_mag = fourier.psd(signal, fs=fs, n_seg=nperseg, n_fft=nperseg, win="hanning", mode="spectrum")  # Calculate Power Spectrum in V**2

        # Determine the indices associated with the frequency ranges.
        range_vlf = (0, max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= 0.166667]))
        range_lf = (max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= 0.166667]),
                    max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= 0.666667]))
        range_if = (max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= 0.5]),
                    max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= 3.0]))

        # Maximum Magnitude
        frq_vlf_max = max(psd_mag[range_vlf[0]:range_vlf[1]])
        frq_lf_max = max(psd_mag[range_lf[0]:range_lf[1]])
        frq_if_max = max(psd_mag[range_if[0]:range_if[1]])

        # Mean Magnitude
        frq_vlf_mean = np.mean(psd_mag[range_vlf[0]:range_vlf[1]])
        frq_lf_mean = np.mean(psd_mag[range_lf[0]:range_lf[1]])
        frq_if_mean = np.mean(psd_mag[range_if[0]:range_if[1]])

        # Median Magnitude
        frq_vlf_med = np.median(psd_mag[range_vlf[0]:range_vlf[1]])
        frq_lf_med = np.median(psd_mag[range_lf[0]:range_lf[1]])
        frq_if_med = np.median(psd_mag[range_if[0]:range_if[1]])

        # Determine the ratios of the features.
        frq_ratio_vlf = frq_vlf_max / frq_if_max  # max(VLF) / max(IF)
        frq_ratio_lf = frq_lf_max / frq_if_mean  # max(LF) / mean(IF)

        if plot:
            fft_dc_drop = 1
            plt.figure()
            plt.title("Frequency Spectrum")
            # plt.plot(fft_bin[fft_dc_drop:], fft_mag[fft_dc_drop:])
            # Calculate number of window intervals in entire signal period with 50% window overlap.
            # psd_n_fft = 2 * int(np.floor(len(signal) / nperseg)) - 1
            plt.plot(psd_bin[fft_dc_drop:], psd_mag[fft_dc_drop:])  # /(psd_n_fft*ds_fs)
            plt.xlim(0.0, 3.3)
            plt.xlabel("frequency")
            plt.ylabel("magnitude")
            plt.axvspan(0.0, 0.166667, color="red", alpha=0.15)
            plt.axvspan(0.166667, 0.666667, color="yellow", alpha=0.15)
            plt.axvspan(0.5, 3.0, color="green", alpha=0.15)
            plt.gca().text(0.98, 0.98, str("VLF/IF: "+"{0:.3f}".format(frq_ratio_vlf)), fontsize=10,
                           ha="right", va="top", transform=plt.gca().transAxes)
            plt.gca().text(0.98, 0.91, str("LF/IF: "+"{0:.3f}".format(frq_ratio_lf)), fontsize=10,
                           ha="right", va="top", transform=plt.gca().transAxes)

        # Print the debug information.
        if debug:
            print(" FREQUENCY DOMAIN")
            print(tabulate([
                ["VLF max", frq_vlf_max],
                ["LF max", frq_lf_max],
                ["IF max", frq_if_max],
                ["VLF mean", frq_vlf_mean],
                ["LF mean", frq_lf_mean],
                ["IF mean", frq_if_mean],
                ["VLF med", frq_vlf_med],
                ["LF med", frq_lf_med],
                ["IF med", frq_if_med],
                ["VLF / IF ratio", frq_ratio_vlf],
                ["LF / IF ratio", frq_ratio_lf],
            ], tablefmt="pretty", colalign=("left", "right")))

        # Return a dictionary containing the collected results.
        return {
            "frq_vlf_max":              frq_vlf_max,
            "frq_lf_max":               frq_lf_max,
            "frq_if_max":               frq_if_max,
            "frq_vlf_mean":             frq_vlf_mean,
            "frq_lf_mean":              frq_lf_mean,
            "frq_if_mean":              frq_if_mean,
            "frq_vlf_med":              frq_vlf_med,
            "frq_lf_med":               frq_lf_med,
            "frq_if_med":               frq_if_med,
            "frq_ratio_vlf":            frq_ratio_vlf,
            "frq_ratio_lf":             frq_ratio_lf,
        }

    def review_artifacts(self, time=None, signal=None, fs=None,
                        debug=False):

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        # Print the debug information.
        if debug:
            print("\nARTIFACTS")

        # Print the debug information.
        if debug:
            print("Not yet implemented.")

        # Return a dictionary containing the collected results.
        return {}

    def review(self, time=None, signal=None, fs=None,
               plot=False,
               debug=True):

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        results = {
            "time_base":                None,
            "fs":                       None,
            "fs_real":                  None,
            "fs_std":                   None,
            "ampl_min":                 None,
            "ampl_max":                 None,
            "ampl_span":                None,
            "ampl_med":                 None,
            "ampl_mean":                None,
            "ampl_num":                 None,
            "ampl_gran":                None,
            "ampl_zc":                  None,
            "ampl_norm-01":             None,
            "ampl_norm-11":             None,
            "ampl_flip":                None,
            "frq_vlf_max":              None,
            "frq_vlf_med":              None,
            "frq_vlf_mean":             None,
            "frq_lf_max":               None,
            "frq_lf_med":               None,
            "frq_lf_mean":              None,
            "frq_if_max":               None,
            "frq_if_med":               None,
            "frq_if_mean":              None,
            "frq_ratio_vlf":            None,
            "frq_ratio_lf":             None,
            "artf_clip":                None,
        }

        import warnings
        warnings.filterwarnings("ignore")

        # Time Base
        res = self.review_timebase(debug=debug)
        results = {key: res.get(key, results[key]) for key in results}  # Update results.

        # Amplitude Characteristics
        res = self.review_amplitude(debug=debug)
        results = {key: res.get(key, results[key]) for key in results}  # Update results.

        # Granularity
        res = self.review_granularity(debug=debug)
        results = {key: res.get(key, results[key]) for key in results}  # Update results.

        # Normalization
        res = self.review_normalization(debug=debug)
        results = {key: res.get(key, results[key]) for key in results}  # Update results.

        # Flipped Amplitude
        res = self.review_flip(plot=plot, debug=debug)
        results = {key: res.get(key, results[key]) for key in results}  # Update results.

        # Frequency Domain
        res = self.review_frequency(plot=plot, debug=debug)
        results = {key: res.get(key, results[key]) for key in results}  # Update results.

        # Artifacts
        res = self.review_artifacts(debug=debug)
        results = {key: res.get(key, results[key]) for key in results}  # Update results.

        # Return a dictionary containing the collected results.
        return results
