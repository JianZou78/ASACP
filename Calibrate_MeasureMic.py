# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 22:01:16 2025

@author: jianzou
"""

import sounddevice as sd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import signal
from datetime import datetime

MIN_LOG_OFFSET = 1e-21

import tkinter as tk

from tkinter import filedialog

def get_input(prompt, default_value):
    user_input = input(f"{prompt} (default: {default_value}): ")
    return user_input if user_input else default_value

def select_sounddevice():
    

    deviceSelect = get_input('Please input DUT name keyword, for example: Fireface\n','Fireface')
    
    apiselect = int(get_input('Please input hostapi index,0=> MME; 1==> Windows DirectSound ; 2==>ASIO ; 3==>Windows WASAPI ; 4==> Windows WDM-KS\n','2'))#3 # 0=> MME; 1==> Windows DirectSound ; 2==>ASIO ; 3==>Windows WASAPI ; 4==> Windows WDM-KS
    
    devName = sd.query_devices()
    NofDev  = len(devName)
    #print(devName)
    ##print(NofDev)
    #print(deviceSelect)
    
    #print(apiselect)
    print('-------------Selected device contains %s, API = %d-----------------\n'%(deviceSelect,apiselect))
    
    
    InChannel= -1;
    OutChannel= -1;
    for LL in range(0,NofDev):
        #print(devName[LL]["name"])
        
        if (devName[LL]["hostapi"] == apiselect and (deviceSelect.lower() in devName[LL]["name"].lower()) and devName[LL]["max_input_channels"] > 0):
            print('Index %d, Api = %d,%s'%(LL,devName[LL]["hostapi"],devName[LL]["name"]))
            InChannel = LL
            print("Input Channel index = %d\n"%(LL))
            break
                  
        
    for LL in range(0,NofDev):        
                
        if( devName[LL]["hostapi"] == apiselect and (deviceSelect.lower() in devName[LL]["name"].lower()) and devName[LL]["max_output_channels"] > 0):
            print('Index %d, Api = %d,%s'%(LL,devName[LL]["hostapi"],devName[LL]["name"]))
            print("Output Channel index = %d\n"%(LL))
            OutChannel = LL
            break
       
    if InChannel== -1 :
        print("Selected Input device not found!")
    if OutChannel==-1 :
        print("Selected Output device not found!")
        
    sd.default.device = [InChannel,OutChannel]
            
    defaultDev = sd.default.device  
    print('-------------Listing the selected device support API -----------------')
         
    for LL in range(0,NofDev):
        if devName[defaultDev[0]]["name"] in devName[LL]["name"]:
            print('Input device: Index %d, Api = %d,%s\n'%(LL,devName[LL]["hostapi"],devName[LL]["name"]))
        elif devName[defaultDev[1]]["name"] in devName[LL]["name"]:
            print('Output device: Index %d, Api = %d,%s\n'%(LL,devName[LL]["hostapi"],devName[LL]["name"]))
            
    return sd.default.device,devName[defaultDev[0]]["name"], devName[defaultDev[0]]["max_input_channels"],devName[defaultDev[1]]["name"]


#======================================================================Start Calibration=============================
print ("Calibration start, insert the measurement mic to the standard sound source, e.g. B&K4231, select/create a txt file to save test output....\n" )

CalSetFile = filedialog.askopenfilename(filetypes=[("Txt files", "*.txt")])
CalFrequency = int(get_input('Calibration Frequency [Hz]\n','1000'))
CalLevel= float(get_input('Calibration Level [dBSPL]\n','94'))
    

print("================Select Audio Frontend  \n====================")
asio_devCH,asio_InName,asio_maxCH,asio_OutName= select_sounddevice();

# Parameters
duration = 10  # seconds
sample_rate = sd.query_devices(sd.default.device, 'input')["default_samplerate"] # Hz

nFFT = 4096;
blocksize_ms = nFFT*1000/sample_rate #
chunk_size = int(blocksize_ms * sample_rate/1000)  # samples per chunk

MIN_LOG = 1e-24

# Initialize plot
fig, ax = plt.subplots()
x = np.arange(0, chunk_size)
line, = ax.plot(x, np.ones(chunk_size))
ax.set_ylim(-1, 1)
ax.set_xlim(0, chunk_size)

figF, axF = plt.subplots()
x =  np.random.rand(chunk_size)
##faxis = np.linspace(0, sample_rate,nFFT)
##yf = np.fft.fft(x, nFFT ) /nFFT
#yfdB = 20*np.log10(abs(yf *2 + MIN_LOG))  # add mirror frequency energy back

faxis,yf=signal.welch(x,sample_rate,'hann', nFFT, nFFT*0.75,nFFT,detrend=False,scaling='density',return_onesided = True)
yfdB = 20*np.log10(yf*sample_rate)
ymaxV = np.max(yfdB)
xmaxHz = faxis[np.argmax(yfdB)]
lineF, = axF.plot(faxis, yfdB)

#ax.set_xlim(0, chunk_size)

def audio_callback(indata, frames, time, status):
    if status:
        print(status)
    global audio_data
    audio_data = indata[:, 0]

def update(frame):
    global audio_data
    line.set_ydata(audio_data)
    rmsT =20*np.log10( np.sqrt(np.mean(audio_data**2))+MIN_LOG)
    ax.set_title(f'RMS={rmsT:.1f},blocksize: {blocksize_ms:.1f} ms ')

    return line,

def update_F(frame):
    global audio_data,max_value,xmaxHz
   # line.set_ydata(audio_data)
    #yf = np.fft.fft(audio_data, nFFT )/nFFT
    #yfdB = 20*np.log10(abs(yf *2+ MIN_LOG))
    faxis,yf=signal.welch(audio_data,sample_rate,'hann', nFFT, nFFT*0.75,nFFT,detrend=False,scaling='density',return_onesided = True)
    yfdB = 20*np.log10(yf*sample_rate)
    lineF.set_ydata(yfdB)
    #max_value = np.max(yfdB)
    max_value = 20*np.log10(np.sqrt( sample_rate/nFFT* np.sum(yf)))
    xmaxHz = faxis[np.argmax(yfdB)]
    axF.set_title(f'nFFT={nFFT:.0f},Max: {max_value:.1f} dB @ {xmaxHz:.0f} Hz,{CalLevel:.0f}dBSPL')
    axF.yaxis.grid(True)
    axF.xaxis.grid(True)
    axF.set_xscale('log')
    axF.minorticks_on()
    axF.grid(which='both',axis='x')
    axF.set_xlabel('Frequency:[Hz]')
    axF.set_ylabel('Amplitude:[dB]')
    axF.set_xlim([100,24000])
    return lineF,

# Start audio stream
stream = sd.InputStream(callback=audio_callback, channels=1, samplerate=sample_rate, blocksize=chunk_size)
with stream:
    ani = FuncAnimation(fig, update, interval=100)
    aniF = FuncAnimation(figF, update_F, interval=100)
    
    
    plt.show()
    aniF.save(os.path.split(CalSetFile)[0]+"\Calibration.gif",writer='pillow')
    
with open(CalSetFile, "a", encoding="utf-8") as file:
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S\n")
    file.write(f'nFFT={nFFT:.0f},Max: {max_value:.1f} dB @ {xmaxHz:.0f} Hz,{CalLevel:.0f}dBSPL\n\n')




