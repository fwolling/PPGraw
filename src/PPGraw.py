# Florian Wolling
# Ubiquitous Computing Lab
# https://ubicomp.eti.uni-siegen.de
# University of Siegen, Germany
#
# "The Quest for Raw Signals: A Quality Review of Publicly Available Photoplethysmography Datasets"
# by Florian Wolling and Kristof Van Laerhoven. In DATA'20: Proceedings of the 3rd Workshop on Data Acquisition To
# Analysis, DATA 2020, Virtual Event, Japan, November 2020, ACM, 2020. https://doi.org/10.1145/3419016.3431485
#
# Version 1.1 (February 2021)


import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
from tabulate import tabulate

__author__ = "Florian Wolling, https://github.com/fwolling/"
__version__ = "1.1"


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

    Version 1.1 (February 2021)

    Attributes
    ----------
    time : None or list of floats, optional, default=None
        The vector that contains the timestamps.
    signal : list of floats, required
        The vector that contains the signal samples.
    fs : float
        The desired sampling frequency in Hz.

    Methods
    -------
    review_timebase
        The regular sampling of devices is often assumed to be constant at a desired rate, but often the sampling period
        is slightly deviating, resulting in a jittering frequency. The tool review_timebase considers the metrics mean
        and standard deviation to review the provided timebase, whether the timestamps are absent, original and raw, or
        artificial and subsequently added. The determined mean is also compared to the desired reference fs.
    review_amplitude
        The amplitude of raw PPG signals naturally shows a large DC offset. Hence, the span of the amplitude, minimum
        and maximum respectively, are utilized by the review_amplitude tool to reveal whether the signal has been
        preprocessed, shifted towards zero, or has even been filtered. Furthermore, the mean and median of the signal
        are derived to provide insights into the characteristics of the PPG signal.
    review_granularity
        Ideally, the granularity of a sampled signal is identical to the amplitude resolution. The metric helps to
        unveil applied preprocessing as the values obtained from an ADC are integers by nature, preprocessed ones are,
        however, usually represented by floats. The review_granularity tool determines the granularity from the sorted
        list of unique values without duplicates by seeking the minimum Euclidean distance.
    review_normalization
        While the scope of a raw PPG signal is only a minor fraction of the overall signal's extent, it is commonly
        cropped at its minimum to reduce the memory demands. Furthermore, it is common to scale and normalize the
        scope to a range of [0,1] or [-1,1]. The review_normalization tool analyzes the signal's minimum and maximum.
    review_flipped
        Most traditional pulse oximetry sensors monitor the PPG signal proportional to the course of the arterial
        blood pressure (ABP) and hence flipped the signal to enable this analogy. The raw signals of both PPG modes,
        however, originally show an inversely proportional course. The review_flipped tool identifies the pulse
        direction by 1) determining the pulses' center of mass and 2) comparing the steepness of the down and up slopes.
    review_frequency
        As all physiological signals, the raw PPG signal is a non-stationary one and dominated by baseline wandering.
        Most approaches apply high-pass filters to remove the low-frequency components and to limit the pulsatile signal
        in a constant boundary envelope. The review_frequency tool analyzes the frequency bands associated with activity
        in the autonomic nervous system (VLF 0.0 to 0.167 Hz) and respiration (LF 0.167 to 0.667 Hz), and compares their
        maxima with the maximum and mean of the frequency band typical for the natural heart rate (IF 0.5 to 3.0 Hz).
    review_artifacts
        Clipping artifacts are common for aged datasets where the signal has been normalized while the caps of the
        lowest or highest peaks are cut due to inaccuracies. The review_artifacts tool detects these flat tops by means
        of multiple successive samples that stay at the constant boundary values for a longer period.
    review
        The review tool sequentially runs all previously provided tools.

    Future Work
    -----------
    The detection of typical motion artifacts is not yet implemented, would however be an interesting quality metric.

    References
    ----------
    [1] "The Quest for Raw Signals: A Quality Review of Publicly Available Photoplethysmography Datasets"
    by Florian Wolling and Kristof Van Laerhoven. In DATA'20: Proceedings of the 3rd Workshop on Data Acquisition To
    Analysis, DATA 2020, Virtual Event, Japan, November 2020, ACM, 2020. https://doi.org/10.1145/3419016.3431485
    """

    # DEFAULT PARAMETERS

    # Default format of the printed result tables.
    # Other format options: https://pypi.org/project/tabulate/
    PRINT_TABLE_FORMAT = "pretty"

    # Error margin in % for the detection of zero-centered signals.
    AMPL_ZC_ERR_MARGIN = 1.5

    # Parameters required for the flip detection.
    AMPL_FLIP_BP_ORDER = 3  # The order of the band-pass filter.
    AMPL_FLIP_BP_FC_LOW = 0.4  # The lower corner frequency.
    AMPL_FLIP_BP_FC_HIGH = 15.0  # The upper corner frequency.

    # Mode of the frequency spectral analysis:
    # * "density", calculate Power Spectral Density in V**2/Hz.
    # * "spectrum", calculate Power Spectrum in V**2.
    FRQ_MODE = "density"
    # Window function applied to the signal.
    FRQ_WIN = "hanning"
    # The length in s of the Hanning window applied in the PSD calculation.
    FRQ_WIN_LEN = 60.0
    FRQ_VLF = [0.0, 0.166667]  # Very low frequency band 0.0 to 0.166667 Hz.
    FRQ_LF = [0.166667, 0.666667]  # Low frequency band 0.166667 to 0.666667 Hz.
    FRQ_IF = [0.5, 3.0]  # Intermediate frequency band 0.5 to 3.0 Hz.
    # HF High frequency band > 3.0 Hz.

    # INITIALIZATION

    def __init__(self, time=None, signal=None, fs=None,
                 table_format=PRINT_TABLE_FORMAT,
                 debug=False):
        """Initialization of the PPGraw review tools.

        :param time: None or list of ints or floats, optional, default=None; The vector that contains the timestamps.
        :param signal: list of ints or floats, optional, default=None; The vector that contains the signal samples.
        :param fs: float, optional, default=None; The desired sampling frequency in Hz.
        :param table_format: String, default=PRINT_TABLE_FORMAT; The default format of the printed result tables.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.
        """
        
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
            self.PRINT_TABLE_FORMAT = table_format

    # PREPROCESSING METHODS

    @staticmethod
    def lowpass(signal, fs, order, fc,
                debug=False):
        """Butterworth forward-backward low-pass filter.

        :param signal: list of ints or floats; The vector containing the signal samples.
        :param fs: float; The sampling frequency in Hz.
        :param order: int; The order of the filter.
        :param fc: int or float; The cutoff frequency of the filter.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: list of floats; The filtered signal.
        """

        nyq = 0.5 * fs  # Calculate the Nyquist frequency.
        cut = fc / nyq  # Calculate the cutoff frequency (-3 dB).
        lp_b, lp_a = sig.butter(order, cut, btype="lowpass")  # Design and apply the low-pass filter.
        lp_data = list(sig.filtfilt(lp_b, lp_a, signal))  # Apply forward-backward filter with linear phase.
        return lp_data

    @staticmethod
    def highpass(signal, fs, order, fc,
                 debug=False):
        """Butterworth forward-backward high-pass filter.

        :param signal: list of ints or floats; The vector containing the signal samples.
        :param fs: float; The sampling frequency in Hz.
        :param order: int; The order of the filter.
        :param fc: int or float; The cutoff frequency of the filter.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: list of floats; The filtered signal.
        """
        nyq = 0.5 * fs  # Calculate the Nyquist frequency.
        cut = fc / nyq  # Calculate the cutoff frequency (-3 dB).
        hp_b, hp_a = sig.butter(order, cut, btype="highpass")  # Design and apply the high-pass filter.
        hp_data = list(sig.filtfilt(hp_b, hp_a, signal))  # Apply forward-backward filter with linear phase.
        return hp_data

    @staticmethod
    def bandpass(signal, fs, order, fc_low, fc_hig,
                 debug=False):
        """Butterworth forward-backward band-pass filter.

        :param signal: list of ints or floats; The vector containing the signal samples.
        :param fs: float; The sampling frequency in Hz.
        :param order: int; The order of the filter.
        :param fc_low: int or float; The lower cutoff frequency of the filter.
        :param fc_hig: int or float; The upper cutoff frequency of the filter.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: list of floats; The filtered signal.
        """
        nyq = 0.5 * fs  # Calculate the Nyquist frequency.
        cut_low = fc_low / nyq  # Calculate the lower cutoff frequency (-3 dB).
        cut_hig = fc_hig / nyq  # Calculate the upper cutoff frequency (-3 dB).
        bp_b, bp_a = sig.butter(order, (cut_low, cut_hig), btype="bandpass")  # Design and apply the band-pass filter.
        bp_data = list(sig.filtfilt(bp_b, bp_a, signal))  # Apply forward-backward filter with linear phase.
        return bp_data

    @staticmethod
    def detrend(signal,
                debug=False):
        """Tool to detrend the signal, hence to remove the DC offset.

        :param signal: list of ints or floats; The vector containing the signal samples.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: list of floats; The detrended signal.
        """

        return signal - np.mean(signal)  # Zero-Centered

    @staticmethod
    def norm01(signal,
               debug=False):
        """Tool to normalize the signal in the range of [0,1].

        :param signal: list of ints or floats; The vector containing the signal samples.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: list of floats; The normalized signal in the range of [0,1].
        """

        return (signal - np.min(signal)) / (np.max(signal) - np.min(signal))  # [0,1]-Normalization

    @staticmethod
    def norm11(signal,
               zc_err_margin=AMPL_ZC_ERR_MARGIN,
               debug=False):
        """Tool to normalize the signal in the range of [-1,1].

        :param signal: list of ints or floats; The vector containing the signal samples.
        :param zc_err_margin: float; The error margin in % for the detection of zero-centered signals.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: list of floats; The normalized signal in the range of [-1,1].
        """
        # Calculate error margin based on the given percentage.
        err = (np.max(signal) - np.min(signal)) * (zc_err_margin/100.0)
        # [0,1]-Normalization
        if -err <= np.mean(signal) <= err:  # Check whether signal is zero-centered.
            if np.max(signal) > abs(np.min(signal)):  # Check whether maximum or minimum is larger.
                return signal / np.max(signal)  # Adjust zero-centered signal according to maximum.
            else:
                return signal / np.min(signal)  # Adjust zero-centered signal according to minimum.
        else:  # Normalize non-zero-centered signal.
            return 2 * (signal - np.min(signal)) / (np.max(signal) - np.min(signal)) - 1.0

    @staticmethod
    def flip(signal,
             debug=False):
        """Tool to flip the signal.

        :param signal: list of ints or floats; The vector containing the signal samples.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: list of floats; The flipped signal.
        """
        return [-v for v in signal]  # Flipping

    # ANALYTICAL METHODS

    @staticmethod
    def psd(signal, fs,
            n_seg, n_fft, win="hanning", mode="densitiy",
            debug=False):
        """Tool to apply Welch's method for the estimation of the power spectral density or power spectrum.

        :param signal: list of ints or floats; The vector containing the signal samples.
        :param fs: float; The sampling frequency in Hz.
        :param n_seg: int; The length of the window segments.
        :param n_fft: int; The number of samples taken for the FFT of each window segment, usually n_fft = n_seg.
        :param win: String, default="hanning";
        :param mode: String, default="density"; The power spectral density "density" where Pxx has units of V**2/Hz or
            the power spectrum "spectrum" where Pxx has units of V**2.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return bin: ndarray; The array of sample frequencies.
        :return psd: ndarray; The power spectral density or power spectrum of the signal.
        """

        # default: noverlap=n_seg//2
        bin, psd = sig.welch(signal, fs, window=win, nperseg=n_seg, nfft=n_fft, return_onesided=True, scaling=mode)
        return bin, psd

    # REVIEW TOOL METHODS

    def review_timebase(self, time=None, signal=None, fs=None,
                        factor=0,  # 0: seconds, 3: milliseconds, 6: microseconds
                        debug=False):
        """Tool to review the timebase regarding real sampling rate and jitter.

        The regular sampling of devices is often assumed to be constant at a desired rate, but often the sampling period
        is slightly deviating, resulting in a jittering frequency. The tool review_timebase considers the metrics mean
        and standard deviation to review the provided timebase, whether the timestamps are absent, original and raw, or
        artificial and subsequently added. The determined mean is also compared to the desired reference fs.

        :param time: None or list of floats, optional, default=None; The vector containing the timestamps.
        :param signal: None or list of floats, optional, default=None; The vector containing the signal samples.
        :param fs: None or float, optional, default=None; The desired sampling frequency in Hz.
        :param factor: int, required, default=0; The factor of the timebase: 0: seconds, 3: milliseconds, 6: microseconds.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: a dict containing:
            "time_base": {False, True}, if time base provided.
            "time_off": float of time in seconds offset of the time base.
            "time_len": length of the recording.
            "fs": float, desired sampling rate in Hz.
            "fs_real": float, real average sampling rate in Hz.
            "fs_std": float, jitter of the real sampling rate in seconds s.
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
            ], tablefmt=self.PRINT_TABLE_FORMAT, colalign=("left", "right")))

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
        """Tool to review the amplitude regarding a natural DC offset.

        The amplitude of raw PPG signals naturally shows a large DC offset. Hence, the span of the amplitude, minimum
        and maximum respectively, are utilized by the review_amplitude tool to reveal whether the signal has been
        preprocessed, shifted towards zero, or has even been filtered. Furthermore, the mean and median of the signal
        are derived to provide insights into the characteristics of the PPG signal.

        :param time: None or list of floats, optional, default=None; The vector containing the timestamps.
        :param signal: None or list of floats, optional, default=None; The vector containing the signal samples.
        :param fs: None or float, optional, default=None; The desired sampling frequency in Hz.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: a dict containing:
            "ampl_min": float; The minimum amplitude value.
            "ampl_max": float; The maximum amplitude value.
            "ampl_span": float; The span of the amplitude.
            "ampl_med": float; The median of the amplitude values.
            "ampl_mean": float; The mean of the amplitude values.
        """

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
            ], tablefmt=self.PRINT_TABLE_FORMAT, colalign=("left", "right")))

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
        """Tool to review the amplitude's granularity respectively resolution.

        Ideally, the granularity of a sampled signal is identical to the amplitude resolution. The metric helps to
        unveil applied preprocessing as the values obtained from an ADC are integers by nature, preprocessed ones are,
        however, usually represented by floats. The review_granularity tool determines the granularity from the sorted
        list of unique values without duplicates by seeking the minimum Euclidean distance.

        :param time: None or list of floats, optional, default=None; The vector containing the timestamps.
        :param signal: None or list of floats, optional, default=None; The vector containing the signal samples.
        :param fs: None or float, optional, default=None; The desired sampling frequency in Hz.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: a dict containing:
            "ampl_num": int; The number of unique values.
            "ampl_gran": int or float; The granularity measure.
        """

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        # Determine granularity respectively amplitude resolution.
        signal_sort = np.unique(signal)  # Sort the values and remove duplicates.
        num_val = len(signal_sort)  # Determine number of unique amplitude values.
        signal_step = [signal_sort[i]-signal_sort[i-1] for i in
                       range(1, len(signal_sort))]  # Determine step sizes.
        granularity = np.min(signal_step)  # Determine granularity.

        # Print the debug information.
        if debug:
            print(" GRANULARITY")
            print(tabulate([
                ["number of values", num_val],
                ["granularity", (int(1) if granularity == 1.0 else granularity)],
            ], tablefmt=self.PRINT_TABLE_FORMAT, colalign=("left", "right")))

        # Return a dictionary containing the collected results.
        return {
            "ampl_num":                 num_val,
            "ampl_gran":                granularity,
        }

    def review_normalization(self, time=None, signal=None, fs=None,
                             debug=False):
        """Tool to review the amplitude regarding cropping or normalization.

        While the scope of a raw PPG signal is only a minor fraction of the overall signal's extent, it is commonly
        cropped at its minimum to reduce the memory demands. Furthermore, it is common to scale and normalize the
        scope to a range of [0,1] or [-1,1]. The review_normalization tool analyzes the signal's minimum and maximum.

        :param time: None or list of floats, optional, default=None; The vector containing the timestamps.
        :param signal: None or list of floats, optional, default=None; The vector containing the signal samples.
        :param fs: None or float, optional, default=None; The desired sampling frequency in Hz.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: a dict containing:
            "ampl_zc": bool; Zero-centered or not.
            "ampl_norm01": bool; [0,1]-normalized or not.
            "ampl_norm11": bool; [-1,1]-normalized or not.
        """

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        # Calculate error margin based on the given percentage.
        err_margin = (np.max(signal) - np.min(signal)) * (self.AMPL_ZC_ERR_MARGIN/100.0)

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
            ], tablefmt=self.PRINT_TABLE_FORMAT, colalign=("left", "right")))

        # Return a dictionary containing the collected results.
        return {
            "ampl_zc":                  ampl_zc,
            "ampl_norm01":              ampl_norm01,
            "ampl_norm11":              ampl_norm11,
        }

    def review_flip(self, time=None, signal=None, fs=None,
                    plot=False,
                    debug=False):
        """Tool to review the amplitude regarding pulse direction and flipping.

        Most traditional pulse oximetry sensors monitor the PPG signal proportional to the course of the arterial
        blood pressure (ABP) and hence flipped the signal to enable this analogy. The raw signals of both PPG modes,
        however, originally show an inversely proportional course. The review_flipped tool identifies the pulse
        direction by 1) determining the pulses' center of mass and 2) comparing the steepness of the down and up slopes.

        :param time: None or list of floats, optional, default=None; The vector containing the timestamps.
        :param signal: None or list of floats, optional, default=None; The vector containing the signal samples.
        :param fs: None or float, optional, default=None; The desired sampling frequency in Hz.
        :param plot: bool, default=False; Flag to enable the plot of the analysis' illustrations.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: a dict containing:
            "ampl_flip": bool; Pulse amplitude direction flipped or not.
            "ampl_flip_hist_max": int; The position of the center of mass in the histogram.
        """

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        # Apply band-pass filter before the analysis.
        signal_bp = signal - np.mean(signal)
        signal_bp = self.bandpass(signal_bp, fs=fs, order=self.AMPL_FLIP_BP_ORDER,
                                  fc_low=self.AMPL_FLIP_BP_FC_LOW, fc_hig=self.AMPL_FLIP_BP_FC_HIGH)

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
            ], tablefmt=self.PRINT_TABLE_FORMAT, colalign=("left", "right")))

        # Return a dictionary containing the collected results.
        return {
            "ampl_flip":                ampl_flip,
            "ampl_flip_hist_max":       hist_max,
        }

    def review_frequency(self, time=None, signal=None, fs=None,
                         plot=False,
                         debug=False):
        """Tool to review the signal's frequency spectrum regarding low-frequency components.

        As all physiological signals, the raw PPG signal is a non-stationary one and dominated by baseline wandering.
        Most approaches apply high-pass filters to remove the low-frequency components and to limit the pulsatile signal
        in a constant boundary envelope. The review_frequency tool analyzes the frequency bands associated with activity
        in the autonomic nervous system (VLF 0.0 to 0.167 Hz) and respiration (LF 0.167 to 0.667 Hz), and compares their
        maxima with the maximum and mean of the frequency band typical for the natural heart rate (IF 0.5 to 3.0 Hz).

        :param time: None or list of floats, optional, default=None; The vector containing the timestamps.
        :param signal: None or list of floats, optional, default=None; The vector containing the signal samples.
        :param fs: None or float, optional, default=None; The desired sampling frequency in Hz.
        :param plot: bool, default=False; Flag to enable the plot of the analysis' illustrations.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: a dict containing:
            "frq_vlf_max": int or float; The maximum value in the VLF frequency band.
            "frq_lf_max": int or float; The maximum value in the LF frequency band.
            "frq_if_max": int or float; The maximum value in the IF frequency band.
            "frq_vlf_mean": float; The mean of the values in the VLF frequency band.
            "frq_lf_mean": float; The mean of the values in the LF frequency band.
            "frq_if_mean": float; The mean of the values in the IF frequency band.
            "frq_vlf_med": int or float; The median of the values in the VLF frequency band.
            "frq_lf_med": int or float; The median of the values in the LF frequency band.
            "frq_if_med": int or float; The median of the values in the IF frequency band.
            "frq_ratio_vlf": float; The ratio of max(VLF) / max(IF).
            "frq_ratio_lf": float; The ratio of max(LF) / mean(IF).
        """

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        nperseg = int(fs*self.FRQ_WIN_LEN)

        # fft_bin, fft_mag = fourier.fft(signal, fs=fs, win="hanning", scale="mag")

        psd_bin, psd_mag = self.psd(signal, fs=fs, n_seg=nperseg, n_fft=nperseg, win=self.FRQ_WIN, mode=self.FRQ_MODE)  # Calculate Power Spectral Density in V**2/Hz
        # psd_bin, psd_mag = fourier.psd(signal, fs=fs, n_seg=nperseg, n_fft=nperseg, win="hanning", mode="spectrum")  # Calculate Power Spectrum in V**2

        # Determine the indices associated with the frequency ranges.
        range_vlf = (max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= self.FRQ_VLF[0]]),
                     max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= self.FRQ_VLF[1]]))
        range_lf = (max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= self.FRQ_LF[0]]),
                    max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= self.FRQ_LF[1]]))
        range_if = (max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= self.FRQ_IF[0]]),
                    max([i_bin for i_bin, k_bin in enumerate(psd_bin) if k_bin <= self.FRQ_IF[1]]))

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
            plt.axvspan(self.FRQ_VLF[0], self.FRQ_VLF[1], color="red", alpha=0.15)
            plt.axvspan(self.FRQ_LF[0], self.FRQ_LF[1], color="yellow", alpha=0.15)
            plt.axvspan(self.FRQ_IF[0], self.FRQ_IF[1], color="green", alpha=0.15)
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
            ], tablefmt=self.PRINT_TABLE_FORMAT, colalign=("left", "right")))

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
                         win=30.0,
                         debug=False):
        """Tool to review the signal amplitude regarding slipping artifacts.

        review_artifacts
        Clipping artifacts are common for aged datasets where the signal has been normalized while the caps of the
        lowest or highest peaks are cut due to inaccuracies. The review_artifacts tool detects these flat tops by means
        of multiple successive samples that stay at the constant boundary values for a longer period.

        :param time: None or list of floats, optional, default=None; The vector containing the timestamps.
        :param signal: None or list of floats, optional, default=None; The vector containing the signal samples.
        :param fs: None or float, optional, default=None; The desired sampling frequency in Hz.
        :param win: int or float, optional, default=30.0; The length of the window in samples (int) or seconds (float).
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: a dict containing:
            "atf_clip_n": int; The total number of clipping artifacts counted.
            "atf_clip_win": int or float; The average number of clipping artifacts per window of length win.
        """

        # Fetch data from instance if no arguments are passed.
        if time is None and signal is None and fs is None:
            time = self.time
            signal = self.signal
            fs = self.fs

        atf_n = 0  # TODO Not yet implemented.
        atf_win = 0.0  # TODO Not yet implemented.

        # Print the debug information.
        if debug:
            print("\nARTIFACTS")

        # Print the debug information.
        if debug:
            pass

        print("Not yet implemented.")  # TODO Not yet implemented.

        # Return a dictionary containing the collected results.
        return {
            "atf_clip_n":               atf_n,
            "atf_clip_win":             atf_win,
        }

    def review(self, time=None, signal=None, fs=None,
               plot=False,
               debug=True):
        """Tool to run all provided review tools sequentially.

        The review tool sequentially runs all previously provided tools.

        :param time: None or list of floats, optional, default=None; The vector containing the timestamps.
        :param signal: None or list of floats, optional, default=None; The vector containing the signal samples.
        :param fs: None or float, optional, default=None; The desired sampling frequency in Hz.
        :param plot: bool, default=False; Flag to enable the plot of the analysis' illustrations.
        :param debug: bool, default=False; Flag to enable the debug mode that prints additional information.

        :return: a dict containing:
            "time_base": {False, True}, if time base provided.
            "fs": float, desired sampling rate in Hz.
            "fs_real": float, real average sampling rate in Hz.
            "fs_std": float, jitter of the real sampling rate in seconds s.
            ---
            "ampl_min": float; The minimum amplitude value.
            "ampl_max": float; The maximum amplitude value.
            "ampl_span": float; The span of the amplitude.
            "ampl_med": float; The median of the amplitude values.
            "ampl_mean": float; The mean of the amplitude values.
            ---
            "ampl_num": int; The number of unique values.
            "ampl_gran": int or float; The granularity measure.
            ---
            "ampl_zc": bool; Zero-centered or not.
            "ampl_norm01": bool; [0,1]-normalized or not.
            "ampl_norm11": bool; [-1,1]-normalized or not.
            ---
            "ampl_flip": bool; The pulse amplitude direction flipped or not.
            ---
            "frq_vlf_max": int or float; The maximum value in the VLF frequency band.
            "frq_vlf_med": int or float; The median of the values in the VLF frequency band.
            "frq_vlf_mean": float; The mean of the values in the VLF frequency band.
            "frq_lf_max": int or float; The maximum value in the LF frequency band.
            "frq_lf_med": int or float; The median of the values in the LF frequency band.
            "frq_lf_mean": float; The mean of the values in the LF frequency band.
            "frq_if_max": int or float; The maximum value in the IF frequency band.
            "frq_if_med": int or float; The median of the values in the IF frequency band.
            "frq_if_mean": float; The mean of the values in the IF frequency band.
            "frq_ratio_vlf": float; The ratio of max(VLF) / max(IF).
            "frq_ratio_lf": float; The ratio of max(LF) / mean(IF).
            ---
            "atf_clip_n": int; The total number of clipping artifacts counted.
            "atf_clip_win": int or float; The average number of clipping artifacts per window of length win.
        """

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
            "ampl_norm01":              None,
            "ampl_norm11":              None,
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
            "atf_clip":                 None,
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
