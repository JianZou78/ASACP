# ASACP
Automatically calibrate the background noise level and equalize the frequency response for Microsoft Azure Speech Spec. This is achieved through a single execution of a Python script. After the calibration process, the tool prompts the user to select the necessary .wav signal to complete the tests.

# Usuage

- python Calibrate_MeasureMic.py
  --
  use calibration sound source to calibrate mic sensitivity and store calibration result to txt file

- python Azure_Speech_TestTool.py
  --
  Follow tool UI to auto calibrate and equalize 8x BGN (BackGround Noise) Speakers, then follow instruction to select or skip real test signal playback from 8 BGN speakers
  
  

