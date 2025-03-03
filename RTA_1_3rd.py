# -*- coding: utf-8 -*-
"""
Get the 1/3 RTA spectrum

Arguments:

    -- data is 1D array
    -- fs is sample rate of the input data
return data:
    
    -- centerFrequency_Hz   1/3 center frequency of 1/3 oct
    -- S_F                  RMS value of the spectrum at 1/3 frequency band             

Created on Tue Apr 21 21:20:33 2020

@author: jianzo
"""

from scipy import signal
from scipy.signal import lfilter
import matplotlib.pyplot as plt
import numpy as np  
from scipy.io import wavfile
from A_weighting import A_weighting
import time

def RTA_1_3rd(data,fs):

    nyquistRate=fs/2.0
    centerFrequency_Hz = np.array([50,63,80,100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000,
    12500, 16000])
    
    factor = np.power(2.0,1.0/6.0)
    lowerCutoffFrequency_Hz=centerFrequency_Hz/factor;
    upperCutoffFrequency_Hz=centerFrequency_Hz*factor;
    
    S_F=[];
    
    for lower,upper in zip(lowerCutoffFrequency_Hz, upperCutoffFrequency_Hz):
        # Design filter
        sos = signal.butter( N=8, Wn=np.array([ lower, 
        upper])/nyquistRate, btype='bandpass', analog=False, 
        output='sos');
    
        # Compute frequency response of the filter.
        w, h = signal.sosfreqz(sos,worN=4096,whole=False,fs=fs)
        
    
        #plt.semilogx(w, 20 * np.log10(abs(h)), 'b')
        filteredSpeech = signal.sosfiltfilt(sos, data)
        S_F.append(np.sqrt(np.mean(filteredSpeech**2)))
    
    # show the 1/3 RTA spectrum    
    #plt.semilogx(centerFrequency_Hz, 20 * np.log10(S_F), 'b')
    #ax=plt.gca()
    #ax.yaxis.grid(True,'major',linestyle='--')
    #ax.xaxis.grid(True,'major',linestyle='--')
    #plt.minorticks_on()
    #plt.grid(True,'minor','both')
    #plt.xticks([20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000])
    #plt.xlim([20, 20000])
   # ticks=np.array([-100:0:10])
    #plt.yticks(np.linspace(-100,0,11))
    #plt.ylim([-100,0])
    #plt.title('1/3 RTA spectrum')
    #plt.show()
    
    return centerFrequency_Hz, S_F
    
    

if __name__ == "__main__":
    fs, data = wavfile.read('C:\\D\\01-Work\\00-standard\\06-TestSignal\\IEEE_269-2010_Male_mono_48_kHz.wav')
    data=data/32768
    startT1=time.time()
    RTA_1_3rd(data,fs)
    
    elapsed_time = time.time() - startT1
    print(' RTA_1_3rd Analysis Time is %.2f second' %(elapsed_time))
   
   
    
   # plt.show()



