# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 22:07:17 2025

@author: jianzou
"""

import numpy as np
from scipy.signal import freqz
import matplotlib.pyplot as plt

def iirbandc(V0, fc, fs, Q):
    # Band cut
    K = np.pi * fc / fs
    V0 = 10 ** (V0 / 20)

    b0 = (1 + K / Q + K * K) / (1 + V0 * K / Q + K * K)
    b1 = 2 * (K * K - 1) / (1 + V0 * K / Q + K * K)
    b2 = (1 - K / Q + K * K) / (1 + V0 * K / Q + K * K)

    a0 = 1
    a1 = b1
    a2 = (1 - V0 * K / Q + K * K) / (1 + V0 * K / Q + K * K)

    b = [b0, b1, b2]
    a = [a0, a1, a2]

    return b, a

def iirbandb(V0, fc, fs, Q):
    # Band boost
    K = np.pi * fc / fs
    V0 = 10 ** (V0 / 20)

    b0 = (1 + V0 * K / Q + K * K) / (1 + K / Q + K * K)
    b1 = 2 * (K * K - 1) / (1 + K / Q + K * K)
    b2 = (1 - V0 * K / Q + K * K) / (1 + K / Q + K * K)

    a0 = 1
    a1 = b1
    a2 = (1 - K / Q + K * K) / (1 + K / Q + K * K)

    b = [b0, b1, b2]
    a = [a0, a1, a2]

    return b, a


def iirbb(V0, fc, fs):
    # Base Boost
    K = 2 * fc / fs
    V0 = 10 ** (V0 / 20)

    b0 = (1 - np.sqrt(2 * V0) * K + V0 * K * K) / (1 + np.sqrt(2) * K + K * K)
    b1 = 2 * (V0 * K * K - 1) / (1 + np.sqrt(2) * K + K * K)
    b2 = (1 + np.sqrt(2 * V0) * K + V0 * K * K) / (1 + np.sqrt(2) * K + K * K)

    a0 = 1
    a1 = 2 * (K * K - 1) / (1 + np.sqrt(2) * K + K * K)
    a2 = (1 - np.sqrt(2) * K + K * K) / (1 + np.sqrt(2) * K + K * K)

    b = [b0, b1, b2]
    a = [a0, a1, a2]

    return b, a

def iirbc(V0, fc, fs):
    # Base Cut
    K = 2 * fc / fs
    V0 = 10 ** (V0 / 20)

    b0 = (1 + np.sqrt(2) * K + K * K) / (1 + np.sqrt(2 * V0) * K + V0 * K * K)
    b1 = 2 * (K * K - 1) / (1 + np.sqrt(2 * V0) * K + V0 * K * K)
    b2 = (1 - np.sqrt(2) * K + K * K) / (1 + np.sqrt(2 * V0) * K + V0 * K * K)

    a0 = 1
    a1 = 2 * (V0 * K * K - 1) / (1 + np.sqrt(2 * V0) * K + V0 * K * K)
    a2 = (1 - np.sqrt(2 * V0) * K + V0 * K * K) / (1 + np.sqrt(2 * V0) * K + V0 * K * K)

    b = [b0, b1, b2]
    a = [a0, a1, a2]

    return b, a

def iirtb(V0, fc, fs):
    # Treble Boost
    K = 2 * fc / fs
    V0 = 10 ** (V0 / 20)

    b0 = (V0 + np.sqrt(2 * V0) * K + K * K) / (1 + np.sqrt(2) * K + K * K)
    b1 = 2 * (K * K - V0) / (1 + np.sqrt(2) * K + K * K)
    b2 = (V0 - np.sqrt(2 * V0) * K + K * K) / (1 + np.sqrt(2) * K + K * K)

    a0 = 1
    a1 = 2 * (K * K - 1) / (1 + np.sqrt(2) * K + K * K)
    a2 = (1 - np.sqrt(2) * K + K * K) / (1 + np.sqrt(2) * K + K * K)

    b = [b0, b1, b2]
    a = [a0, a1, a2]

    return b, a

def iirtc(V0, fc, fs):
    # Treble Cut
    K = 2 * fc / fs
    V0 = 10 ** (V0 / 20)

    b0 = (1 - np.sqrt(2) * K + K * K) / (V0 + np.sqrt(2 * V0) * K + V0 * K * K)
    b1 = 2 * (K * K - 1) / (V0 + np.sqrt(2 * V0) * K + V0 * K * K)
    b2 = (1 + np.sqrt(2) * K + K * K) / (V0 + np.sqrt(2 * V0) * K + V0 * K * K)

    a0 = 1
    a1 = 2 * (K * K - V0) / (V0 + np.sqrt(2 * V0) * K + V0 * K * K)
    a2 = (V0 - np.sqrt(2 * V0) * K + K * K) / (V0 + np.sqrt(2 * V0) * K + V0 * K * K)

    b = [b0, b1, b2]
    a = [a0, a1, a2]

    return b, a
##
### Example usage
##V0 = -7  # Example value
##fc = 5000  # Example value
##fs = 48000  # Example value
##Q = 3 # Example value
##
###b, a = iirbandc(V0, fc, fs, Q)
###b, a = iirbandb(V0, fc, fs, Q)
###b, a = iirbb(V0, fc, fs)
###b, a = iirbc(V0, fc, fs)
##b, a = iirtc(V0, fc, fs)
##
##
### Plotting the frequency response
##w, h = freqz(b, a, worN=8000, fs=fs)
##plt.plot(w, 20 * np.log10(abs(h)))
##plt.xscale('log')
##plt.xlim([100, fs / 2])
##plt.title(f'V0={V0} fc={fc} Q={Q}')
##plt.xlabel('Frequency (Hz)')
##plt.ylabel('Magnitude (dB)')
##plt.grid()
##plt.show()
