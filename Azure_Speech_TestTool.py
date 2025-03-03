import sounddevice as sd
import soundfile as sf
import os
import scipy
from scipy import signal
from scipy.signal import lfilter
from datetime import datetime
from RTA_1_3rd import RTA_1_3rd
from IIRParametric import iirbandc
 
import matplotlib.pyplot as plt
import numpy as np
from A_weighting import A_weighting
import yulewalker as yw
import time

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

def LoadCalibration(fPath=r"\\sox-bgn\_shared_folder\Tools_to_share\SpeechPlatformToolchain_V10.1-WithTeamsMeetingRoomsDevices-Release\SpeechPlatform-OnLabPC\CalibrationWAV"):
    fnameEQ = r"EQCalibration-8ch.wav"
    fnameLvl = r"LevelCalibration-8ch.wav"
    d_EQ,fs = sf.read(os.path.join(fPath,fnameEQ))
    d_Lvl,fs = sf.read(os.path.join(fPath,fnameLvl))
    return d_EQ, d_Lvl, fs

# Create a Tkinter root window (it will be hidden)

print (" select the calibration wavfile folder ....\n" )
root = tk.Tk()

root.withdraw()
 
# Open a folder dialog and get the selected folder path

Calibrationfolder_path = filedialog.askdirectory()
 
# Print the selected folder path

print(f"Selected folder: {Calibrationfolder_path}")

root = tk.Tk()

root.withdraw()

# Open a file dialog to select a .txt file

print (" select a txt file to save test output....\n" )

CalSetFile = filedialog.askopenfilename(filetypes=[("Txt files", "*.txt")])


d_EQ, d_Lvl, fs = LoadCalibration(fPath= Calibrationfolder_path)

asio_devCH,asio_InName,asio_maxCH,asio_OutName= select_sounddevice();
#sd.default.device = [14,14] # the number may change later


#Feb 23,2025, 94dBSPL 1kHz, calibration on fireface 802 default input channel 1
#np.sqrt(np.mean(myrecord**2))
#
MicSdB=20*np.log10(0.0011411147801350208) +3 # ~3dB difference vs ACQUA



#asio_out = sd.AsioSettings(channel_selectors=[0,1,2,3]);sd.play(d_EQ[0:10*fs,0:3]*0.1, fs, extra_settings=asio_out)

Lvl_0 = np.zeros(8) # to store intial level from each BGN speaker
RecT = int(get_input('Please input time duration for calibration file playback and record in seconds, max=30s\n','10')) # time to record for analysis, max = 30s for calibration
Scale_48dB = np.ones(8)
FineTuneG = 0;
#Scale_48dB=np.power(10,(48 - FineTuneG -Lvl_0)/20)  #scaling factor to get 48dBSPL(A), fine tune

Lvl_48dB = np.zeros(8);
##print('-------------Playing 8 channel Level Calibration wavfile on the BGN speaker -----------------')
##Lvl_target_Lvl_1spk = float(get_input("Each BGN speaker target dBSPL(A) Level with level calibration wav\n","57"))
##for CH in range(0,8):
##    asio_out = sd.AsioSettings(channel_selectors=[CH])
##    asio_in = sd.AsioSettings(channel_selectors=[0])
##    sd.default.extra_settings = asio_in, asio_out
##    #sd.play(d_EQ[0:3*fs,CH]*0.1, fs, extra_settings=asio_out)
##    
##    myrecord = sd.playrec(d_Lvl[0:RecT*fs,CH], fs, channels=1)
##    
##    time.sleep(RecT +1 )
##    myrecord = myrecord.flatten()
##    #add A weighting
##    b, a = A_weighting(fs)
##    y = lfilter(b, a, myrecord)
##    L=4096;
##    
##    f,Pxx_den=signal.welch(y,fs,'hann', L, L*0.75,L,detrend=False,scaling='density',return_onesided = True)
##    Pxx_den = Pxx_den + MIN_LOG_OFFSET
##    
##    frange=[20, 20000]
##    fN1=np.argmin(np.abs(f-frange[0]))
##    fN2=np.argmin(np.abs(f-frange[1]))
##    ResultLevel=20*np.log10(np.sqrt( fs/L* np.sum(Pxx_den[fN1:fN2:1]))) - MicSdB +94
##    print('BGN speaker %d Level is %.1f dBSPL(A) with 0dB Gain \n'%(CH +1,ResultLevel))
##    Lvl_0[CH] = ResultLevel
##
##with open(CalSetFile, "a", encoding="utf-8") as file:
##    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
##    file.write(f"Time:{current_time}\n")
##    file.write(f"Initial Gain=0dB, BGN speaker Level: \n{Lvl_0}\n")
##
##
###print('-------------Playing 8 channel Level calibration wavfile on the BGN speaker to get %.1f dBSPL(A)-----------------'%(Lvl_target_Lvl_1spk))
##
##
##
##with open(CalSetFile, "a", encoding="utf-8") as file:
##    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
##    file.write(f"Time:{current_time}\n")
##    file.write(f"Recording Time {RecT} seconds...\n")
##    file.write(f"Final Scaling Factor to get target level:\n {Scale_48dB}\n")
##    file.write(f"Final Check level, BGN speaker Level:\n {Lvl_48dB}\n")
##
########################################################## EQ 
Lvl_EQ0 = np.zeros(8) # to store intial level from each BGN speaker
Lvl_EQ_eq = np.zeros(8) # to store intial level from each BGN speaker
d_EQ_eq1 = np.zeros(d_EQ.shape) # to store time signal of EQ'ed

N_tap = int(get_input("IIR filter order\n","16")) # filter tap to design yulewalk filter

Lvl_target_EQ_1spk = float(get_input("Each BGN speaker target dBSPL(A) Level\n","48"))

b_eq = np.ones((N_tap +1, 8)) # yulewalk filter numerator 分子
a_eq = np.ones((N_tap +1, 8)) # yulewalk filter Denominator 分子
FR_EQ = np.zeros((26,8))   # 1/3oct 50 Hz to 16kHz, total 26 frequencies, 8 Channels

Q = np.power(2,4/3) - np.power(2,-4/3) # the last IIR filter Quality factor

b_eq_add_iir = np.ones((N_tap + 3, 8))
a_eq_add_iir = np.ones((N_tap + 3, 8))

FR_EQ = np.zeros((26,8))   # 1/3oct 50 Hz to 16kHz, total 26 frequencies, 8 Channels, to store final EQ'ed FR

print('-------------Playing 8 channel Equlization wavfile on the BGN speaker -----------------')
for CH in range(0,8):
    asio_out = sd.AsioSettings(channel_selectors=[CH])
    asio_in = sd.AsioSettings(channel_selectors=[0])
    sd.default.extra_settings = asio_in, asio_out
   
    myrecord = sd.playrec(d_EQ[0:RecT *fs,CH], fs, channels=1)
    
    time.sleep(RecT +1)
    myrecord = myrecord.flatten()
    #add A weighting
    b, a = A_weighting(fs)
    y = lfilter(b, a, myrecord)
    L=4096;
    
    f,Pxx_den=signal.welch(y,fs,'hann', L, L*0.75,L,detrend=False,scaling='density',return_onesided = True)
    Pxx_den = Pxx_den + MIN_LOG_OFFSET
    
    frange=[20, 20000]
    fN1=np.argmin(np.abs(f-frange[0]))
    fN2=np.argmin(np.abs(f-frange[1]))
    ResultLevel=20*np.log10(np.sqrt( fs/L* np.sum(Pxx_den[fN1:fN2:1]))) - MicSdB +94
    print('BGN speaker  %d playing EQCalibration file, Level is %.1f dBSPL(A), Gain 0dB\n'%(CH +1,ResultLevel))
    Lvl_EQ0[CH] = ResultLevel

    count = 0
    while True:
        Scale_48dB[CH]=np.power(10,(Lvl_target_EQ_1spk - FineTuneG -Lvl_EQ0[CH])/20)  #scaling factor to get 48dBSPL(A), fine tune
        myrecord = sd.playrec(d_EQ[0:RecT*fs,CH] * Scale_48dB[CH] , fs, channels=1)
        time.sleep(RecT +1 )
        myrecord = myrecord.flatten()
        #add A weighting
        b, a = A_weighting(fs)
        y = lfilter(b, a, myrecord)
        L=4096;
        f,Pxx_den=signal.welch(y,fs,'hann', L, L*0.75,L,detrend=False,scaling='density',return_onesided = True)
        Pxx_den = Pxx_den + MIN_LOG_OFFSET
    
        frange=[20, 20000]
        fN1=np.argmin(np.abs(f-frange[0]))
        fN2=np.argmin(np.abs(f-frange[1]))
        ResultLevel=20*np.log10(np.sqrt( fs/L* np.sum(Pxx_den[fN1:fN2:1]))) - MicSdB +94
        print('BGN speaker %d Level is %.1f dBSPL(A), Try %d , Gain = %.1fdB ...\n'%(CH +1,ResultLevel, count +1, 20*np.log10(Scale_48dB[CH])))
        Lvl_48dB[CH]= ResultLevel
        FineTuneG = FineTuneG + Lvl_48dB[CH] - Lvl_target_EQ_1spk
        

        CF,F_data = RTA_1_3rd(myrecord[0:RecT*fs],fs);

        fig,ax0 = plt.subplots(1,figsize=(16,9))
        ax0.semilogx(CF, 20*np.log10(np.array(F_data)) + 94 - MicSdB,'r',label='BGN Speaker '+str(CH+1)+' level='+f"{ResultLevel:.1f}"+" dBSPL(A)")
        ax0.yaxis.grid(True,'major',linestyle='--')
        ax0.xaxis.grid(True,'minor')
        plt.legend(loc="upper left")
        plt.xlim([100, 20000])
        plt.title('BGN speaker RTA 1/3oct noise spectrum')
        plt.savefig(os.path.split(CalSetFile)[0]+"\BGN_FR_"+str(CH+1)+".png")#"\\sox-bgn\_shared_folder\Tools_to_share\SpeechPlatformToolchain_V10.1-WithTeamsMeetingRoomsDevices-Release\SpeechPlatform-OnLabPC\CalibrationWAV\BGN_FR_"+str(CH+1)+".png")
    
        count = count + 1

        if abs( Lvl_48dB[CH]- Lvl_target_EQ_1spk)<=1:
            break
        if count >=2:
            break

 ### get target filter for yulewalk
    frange=[160, 8000] # in Hz
    fN1=np.argmin(np.abs(CF-frange[0]))
    fN2=np.argmin(np.abs(CF-frange[1]))

    m_dB = 20*np.log10(np.array(F_data)) + 94 - MicSdB
    m_dB1 = np.mean(m_dB[fN1:fN2])-m_dB

    m1 = np.concatenate((np.concatenate (([0,0], m_dB1[fN1:fN2+1])),[0]))
    f1 = np.concatenate((np.concatenate (([0], CF[fN1-1:fN2+1])),[fs/2]))

    a1, b1 = yw.yulewalk(N_tap, f1/fs*2, np.power(10,m1/20))
    w, h = scipy.signal.freqz(b1, a1)
    b_eq[:,CH] = b1
    a_eq[:,CH] = a1
##    print("BGN Speaker %d EQ coefficient Numerator b:\n"%(CH+1))
##    print(a1)
##    print("\n")
##    print("BGN Speaker %d EQ coefficient Denominator a:\n"%(CH+1))
##    print(b1)
##    print("\n")
 
# show target and yulewalk estimated
    #fig, ax = plt.subplots(2, figsize=(7, 10));ax[0].set_ylabel('Magnitude');ax[1].set_ylabel('Phase');ax[1].set_xlabel('Hz');ax[0].plot(w/np.pi*48000/2, 20*np.log10(np.abs(h)));ax[0].plot(f1, m1);ax[0].legend(['yw-estimated', 'desired']);ax[0].set_xscale('log');plt.show()

### get yulewalk filtered d_EQ signal
    d_EQ_eq1[:,CH] = lfilter(b1, a1, d_EQ[:,CH])
    
with open(CalSetFile, "a", encoding="utf-8") as file:
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    file.write(f"Time:{current_time}\n")
    file.write(f"IIR filter tap :\t {N_tap}\n")
    file.write(f"EQ coefficient Demonimator a :\n {a_eq}\n")
    file.write(f"EQ coefficient Numerator b :\n {b_eq}\n") 


for CH in range(0,8):
    asio_out = sd.AsioSettings(channel_selectors=[CH])
    asio_in = sd.AsioSettings(channel_selectors=[0])
    sd.default.extra_settings = asio_in, asio_out
    #sd.play(d_EQ[0:3*fs,CH]*0.1, fs, extra_settings=asio_out)
    
    myrecord_EQ = sd.playrec(d_EQ_eq1[0:RecT*fs,CH] * Scale_48dB[CH], fs, channels=1)
    
    time.sleep(RecT +1 )
    myrecord_EQ = myrecord_EQ.flatten()
    #add A weighting
    b, a = A_weighting(fs)
    y = lfilter(b, a, myrecord_EQ)
    L=4096;
    
    f,Pxx_den=signal.welch(y,fs,'hann', L, L*0.75,L,detrend=False,scaling='density',return_onesided = True)
    Pxx_den = Pxx_den + MIN_LOG_OFFSET
    
    frange=[20, 20000]
    fN1=np.argmin(np.abs(f-frange[0]))
    fN2=np.argmin(np.abs(f-frange[1]))
    ResultLevel=20*np.log10(np.sqrt( fs/L* np.sum(Pxx_den[fN1:fN2:1]))) - MicSdB +94
    print('BGN speaker  %d playing EQCalibration file with EQ, Level is %.1f dBSPL(A)\n'%(CH+1,ResultLevel))
    Lvl_EQ_eq[CH] = ResultLevel

    #f1,Pxx_den1=signal.welch(y,fs,'flattop', L1,L1*0.75,L1, detrend=False,scaling='spectrum',return_onesided = True)
    # calculate BGN speaker A weighting level

    CF,F_data = RTA_1_3rd(myrecord_EQ[0:RecT*fs],fs);

    fig,ax = plt.subplots(1,figsize=(16,9))
    FR_EQ[:,CH] = 20*np.log10(np.array(F_data)) + 94 - MicSdB
    ax.semilogx(CF, FR_EQ[:,CH],'r',label='BGN Speaker '+str(CH+1)+' level='+f"{ResultLevel:.1f}"+" dBSPL(A)")

    frange=[160, 8000] # in Hz
    fN1=np.argmin(np.abs(CF-frange[0]))
    fN2=np.argmin(np.abs(CF-frange[1]))

    #check the max failure frequency point, if fail
    fr1 = FR_EQ[:,CH]- np.min(FR_EQ[fN1:fN2+1,CH])
    
    if fr1[np.argmax(fr1[fN1:fN2+1])+fN1] > 6:
        print("FR_EQ is fail, adding 2nd order parametric IIR band cut filter... \n")

        b_bandc, a_bandc = iirbandc( fr1[np.argmax(fr1[fN1:fN2+1])+fN1] -3, CF[np.argmax(fr1[fN1:fN2+1])+fN1], fs, Q)
        b_eq_add_iir[:,CH] = np.convolve(b_eq[:,CH],b_bandc)
        a_eq_add_iir[:,CH] = np.convolve(a_eq[:,CH],a_bandc)

        #recalculate time signal, and play analyze again
        d_EQ_eq1[:,CH] = lfilter(b_eq_add_iir[:,CH], a_eq_add_iir[:,CH], d_EQ[:,CH])
  

        myrecord_EQ = sd.playrec(d_EQ_eq1[0:RecT*fs,CH] * Scale_48dB[CH], fs, channels=1)
    
        time.sleep(RecT +1 )
        myrecord_EQ = myrecord_EQ.flatten()
        #add A weighting
        b, a = A_weighting(fs)
        y = lfilter(b, a, myrecord_EQ)
        L=4096;
        
        f,Pxx_den=signal.welch(y,fs,'hann', L, L*0.75,L,detrend=False,scaling='density',return_onesided = True)
        Pxx_den = Pxx_den + MIN_LOG_OFFSET
        
        frange=[20, 20000]
        fN1=np.argmin(np.abs(f-frange[0]))
        fN2=np.argmin(np.abs(f-frange[1]))
        ResultLevel=20*np.log10(np.sqrt( fs/L* np.sum(Pxx_den[fN1:fN2:1]))) - MicSdB +94
        print('BGN speaker  %d playing EQCalibration file with EQ and parametric IIR, Level is %.1f dBSPL(A)\n'%(CH+1,ResultLevel))
        Lvl_EQ_eq[CH] = ResultLevel

        #f1,Pxx_den1=signal.welch(y,fs,'flattop', L1,L1*0.75,L1, detrend=False,scaling='spectrum',return_onesided = True)
        # calculate BGN speaker A weighting level

        CF,F_data = RTA_1_3rd(myrecord_EQ[0:RecT*fs],fs);
        frange=[160, 8000] # in Hz
        fN1=np.argmin(np.abs(CF-frange[0]))
        fN2=np.argmin(np.abs(CF-frange[1]))

        fig,ax = plt.subplots(1,figsize=(16,9))
        FR_EQ[:,CH] = 20*np.log10(np.array(F_data)) + 94 - MicSdB
        ax.semilogx(CF, FR_EQ[:,CH],'r',label='BGN Speaker '+str(CH+1)+' level='+f"{ResultLevel:.1f}"+" dBSPL(A)")
        ax.semilogx([160,8000], np.min(FR_EQ[:,CH][fN1:fN2+1])*np.array([1,1]),linestyle = '--',color='k',label='upper')
        ax.semilogx([160,8000], np.min(FR_EQ[:,CH][fN1:fN2+1])*np.array([1,1]) + 6,linestyle = '--',color='k',label='lower')
        
        #plt.semilogx(f1,10*np.log10(Pxx_den1)+20*np.log10(L1/L),'g', label='Receiving Idle Noise SFI');
        #ax=plt.gca()
        ax.yaxis.grid(True,'major',linestyle='--')
        ax.xaxis.grid(True,'minor')
        plt.legend(loc="upper left")
        ax.set_xlim((100, 10000))
        #plt.ylim([-120,0])
        plt.title('BGN speaker RTA 1/3oct noise spectrum with EQ')
        plt.savefig(os.path.split(CalSetFile)[0]+ "\BGN_FR_EQ_iirP"+str(CH+1)+".png")
        plt.close()
        #check the max failure frequency point, if fail
        fr1 = FR_EQ[:,CH]- np.min(FR_EQ[fN1:fN2+1,CH])
        if fr1[np.argmax(fr1[fN1:fN2+1])+fN1] > 6:
            print(f"FR_EQ with parametric IIR is still fail {(fr1[np.argmax(fr1[fN1:fN2+1])+fN1] -6):.1f}dB... \n")
        else:
            print("FR_EQ with parametric IIR is Pass.. \n")
 
    
    else:
        print("FR_EQ is pass...\n")
        b_eq_add_iir[:,CH] = np.convolve(b_eq[:,CH],np.array([1,0,0]))
        a_eq_add_iir[:,CH] = np.convolve(a_eq[:,CH],np.array([1,0,0]))
        ax.semilogx([160,8000], np.min(FR_EQ[:,CH][fN1:fN2+1])*np.array([1,1]),linestyle = '--',color='k',label='upper')
        ax.semilogx([160,8000], np.min(FR_EQ[:,CH][fN1:fN2+1])*np.array([1,1]) + 6,linestyle = '--',color='k',label='lower')
        
        #plt.semilogx(f1,10*np.log10(Pxx_den1)+20*np.log10(L1/L),'g', label='Receiving Idle Noise SFI');
        #ax=plt.gca()
        ax.yaxis.grid(True,'major',linestyle='--')
        ax.xaxis.grid(True,'minor')
        plt.legend(loc="upper left")
        ax.set_xlim((100, 10000))
        #plt.ylim([-120,0])
        plt.title('BGN speaker RTA 1/3oct noise spectrum with EQ')
        plt.savefig(os.path.split(CalSetFile)[0]+"\BGN_FR_EQ"+str(CH+1)+".png")
        plt.close()

  
with open(CalSetFile, "a", encoding="utf-8") as file:
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    file.write(f"Time:{current_time}\n")
    file.write(f"IIR filter tap :\t {N_tap}\n")
    file.write(f"EQ coefficient Demonimator a :\n {a_eq_add_iir}\n")
    file.write(f"EQ coefficient Numerator b :\n {b_eq_add_iir}\n") 

  
#--------------play d_EQ_eq1 from 8 BGN speakers and check level

asio_out = sd.AsioSettings(channel_selectors=[0,1,2,3,4,5,6,7])
asio_in = sd.AsioSettings(channel_selectors=[0])
sd.default.extra_settings = asio_in, asio_out
#sd.play(d_EQ[0:3*fs,CH]*0.1, fs, extra_settings=asio_out)

myrecord_EQ = sd.playrec(d_EQ_eq1[0:RecT*fs,:] * Scale_48dB, fs, channels=1)

time.sleep(RecT +1 )
myrecord_EQ = myrecord_EQ.flatten()
#add A weighting
b, a = A_weighting(fs)
y = lfilter(b, a, myrecord_EQ)
L=4096;

f,Pxx_den=signal.welch(y,fs,'hann', L, L*0.75,L,detrend=False,scaling='density',return_onesided = True)
Pxx_den = Pxx_den + MIN_LOG_OFFSET

frange=[20, 20000]
fN1=np.argmin(np.abs(f-frange[0]))
fN2=np.argmin(np.abs(f-frange[1]))
ResultLevel_8CH_EQ=20*np.log10(np.sqrt( fs/L* np.sum(Pxx_den[fN1:fN2:1]))) - MicSdB +94
print('BGN speaker  1-8 playing EQCalibration file with EQ, Level is %.1f dBSPL(A)\n'%(ResultLevel_8CH_EQ))

if abs(ResultLevel_8CH_EQ - 57) <=1:
    print("EQ calibration done, 8 BGN Spk level within +/-1 dB target!")

else:
    Scale_48dB = Scale_48dB * np.log10(10, (57-ResultLevel_8CH_EQ)/20)
    myrecord_EQ = sd.playrec(d_EQ_eq1[0:RecT*fs,:] * Scale_48dB, fs, channels=1)

    time.sleep(RecT +1 )
    myrecord_EQ = myrecord_EQ.flatten()
    #add A weighting
    b, a = A_weighting(fs)
    y = lfilter(b, a, myrecord_EQ)
    L=4096;

    f,Pxx_den=signal.welch(y,fs,'hann', L, L*0.75,L,detrend=False,scaling='density',return_onesided = True)
    Pxx_den = Pxx_den + MIN_LOG_OFFSET

    frange=[20, 20000]
    fN1=np.argmin(np.abs(f-frange[0]))
    fN2=np.argmin(np.abs(f-frange[1]))
    ResultLevel_8CH_EQ=20*np.log10(np.sqrt( fs/L* np.sum(Pxx_den[fN1:fN2:1]))) - MicSdB +94
    print('BGN speaker  1-8 playing EQCalibration file with EQ, Level is %.1f dBSPL(A)\n'%(ResultLevel_8CH_EQ))

###### 
        
              
# ----------------- ask if start Lab PC test wav playback and select wavfile
# Ask  user for input

user_input = input("Do you want to select the test wavfile to start test ? (Y/N): ")

if user_input.lower() == 'y':

    # Create a Tkinter root window (it will be hidden)

    root = tk.Tk()

    root.withdraw()

    # Open a file dialog to select a .wav file

    file_path = filedialog.askopenfilename(filetypes=[("WAV files", "*.wav")])

    if file_path:

        print(f"Selected file: {file_path}")
        d_TestWav,fs = sf.read(file_path)
        d_EQ_TestWav = np.zeros(d_TestWav.shape) ## to store time signal of EQ'ed Test wavfile

        for CH in range (0, 8):
            d_EQ_TestWav[:,CH] = lfilter(b_eq_add_iir[:,CH], a_eq_add_iir[:,CH], d_TestWav[:,CH])

        
        asio_out = sd.AsioSettings(channel_selectors=[0,1,2,3,4,5,6,7])
        asio_in = sd.AsioSettings(channel_selectors=[0])
        sd.default.extra_settings = asio_in, asio_out
        print("Testing, waiting for %.1f minutes...\n"%(d_EQ_TestWav.shape[0]/fs/60))

        sd.play(d_TestWav* Scale_48dB, fs)

        

        time.sleep(d_EQ_TestWav.shape[0]/fs)

        print("test signal playback finish...\n")
        

        
        

    else:

        print("No file selected.")

else:

    print("Quitting the program.")

   

    
