# -*- coding: utf-8 -*-
'''
Created on Tue Jul 14 12:28:39 2015
@author: Yaoyao
'''
#Measurement Settings  
pointCount = 2000
timeRange = 1e-6
position = 16.15e-6
VRange = 1
avgCount = 16
feature = 'Wall'
savePath = '..\\'
CH_list = ['CH0','CH7']
channels = [1,2]
VOffset = [-0.9, -1.4]
VRange = [1, 0.1]

#Data Process Settings
dataProcess = True
FFT_points = 2**14
gain = 22
attn = 20
drive = 10
Vin = 2 * drive * 10 ** (-attn/20)

import visa
import sys
import numpy as np
import time
import winsound
import scipy
import scipy.signal
import numpy.fft
import matplotlib.pyplot as plt
import os

###############################################################
#Measurement Functions
###############################################################  
def scope_init():
    rm = visa.ResourceManager()
    res_list = rm.list_resources()
    if res_list == None:
        exit()
    scope = rm.get_instrument(res_list[0])
    IDN = scope.ask('*IDN?')
    model = IDN.split(',')[1]
    print model
    if model != 'MSO-X 3052A':
        print 'Wrong instrument!'
        sys.exit(0)
    return scope
  
def setup(scope):
    scope.write('*RST')
    scope.write('*CLS')
    scope.write(':AUT')
    scope.write(':CHAN1:PROB 1')
    scope.write(':CHAN1:RANG 8')
    scope.write(':CHAN1:SCAL 0.04')
    scope.write(':CHAN1:INP ONEM')
    scope.write(':MEAS:VPP CHAN1')
    
    scope.write(':CHAN2:PROB 1')
    scope.write(':CHAN2:RANG 8')
    scope.write(':CHAN2:SCAL 0.04')
    scope.write(':CHAN2:INP ONEM')
    scope.write(':MEAS:VPP CHAN2')
    
    scope.write(':TIM:RANG 5e-6')
    scope.write(':ACQ:MODE RTIM')
    scope.write(':ACQ:TYPE AVER')
    scope.write(':ACQ:COUN 4')
    scope.write(':TIM:REF CENT')
    scope.write(':TIM:POS 15e-6')
    scope.write(':TRIG:SOUR EXT')
    scope.write(':TRIG:MODE EDGE')
    scope.write(':TRIG:EDGE:LEV 1.5') 
    scope.write(':WAV:POIN 2000')
    time.sleep(1)


def setVoltage(scope, Range = 0.3, offset = 0, source = 1, autoRange = False):
    if autoRange:
        scope.write(':CHAN%s:RANG 5' %(source))
        Vmax = float(scope.ask(':MEAS:VMAX? CHAN1'))
        Vmin = float(scope.ask(':MEAS:VMIN? CHAN1'))
        Range = round((Vmax-Vmin)*1.5,1)
    scope.write(':CHAN%s:RANG %s' %(source, Range))
    scope.write(':CHAN%s:OFFS %s' %(source, offset))

def setTime(scope, timeRange = 5e-6, position = 0):
    scope.write(':TIM:RANG %s' %(timeRange))
    scope.write(':TIM:POS %s' %(position))
    
def get_waveform(scope,  fmt = 'ASCII', source = 1, points = 2000, avg = 0):
    scope.write(':WAV:POIN:MODE NORM')
    scope.write(':WAV:FORM %s' %(fmt))
    scope.write(':WAV:SOUR CHAN%s' %(source))
    scope.write(':WAV:POIN %s' %(points))
    if avg > 1 and avg <= 65532:
        scope.write(':ACQ:TYPE AVER')
        scope.write(':ACQ:COUN %s' %(avg))
       
    scope.write(':WAV:DATA?')
    time.sleep(1)
    rawData = scope.read()[10:]
    temp = rawData.split(',')
    data = np.array(temp,float)
    points = len(data)
    return points, data

def get_timeData(timeRange, timeCenter, time_points):
#    timeRange = float(scope.ask(':TIM:RANG?'))
#    timeCenter = float(scope.ask(':TIM:POS?'))
    timeStep = timeRange/(time_points-1)
    time_data = np.linspace(timeCenter-timeRange/2, timeCenter+timeRange/2, time_points)
    return time_data

def saveData(timeData, waveform, filePath):
    data = np.transpose(np.vstack([timeData, waveform]))
    np.savetxt(filePath,data)
    
###############################################################
#Data Process Functions
###############################################################
def loadTraceFile(fname):
#Select a new file
    data = None
    #try:
    if True:
        data = np.loadtxt(fname)    
    traceTime = 0.
    traceVal = 0.
    
    try:
       if data is not None:
           traceTime = data[:,0]
           traceVal = data[:,1]
    except:
        print 'Data Process Error:  ' + str(fname)
        pass
    
    return traceTime, traceVal
    
def calcNormFactor(V):
    K = np.max(V) 
    return K
    
def calcTau(V0,V1,dt):
    xCorrel = scipy.signal.correlate(V0,V1)
    tau = dt * (np.argmax(xCorrel)-len(V0)-1)
    return tau

def calcSens(Vpp, drive, gain, attn):
    Sens = 20*np.log10(Vpp/(2*drive))+attn*2-gain
    return Sens
    
def calcBW(V,dt,cutoff=0):
    x_data = np.fft.fftfreq(FFT_points,dt)/1.0E+06
    y_data = 20.*np.log10(numpy.fft.fft(V, n=FFT_points))
    y_data = y_data - 60*(x_data<5)
    y_data -= np.max(y_data)      
    fl = min(x_data[y_data>-6.0])
    fh = max(x_data[y_data>-6.0])
    fp = x_data[y_data==0]
    bw = 100.*2.*(fh/fl-1.)/(fh/fl+1.)
    fc = (fh+fl)/2
    return x_data, y_data, fl, fh, fp, fc, bw
    

#Init scope   
scope = scope_init()
setup(scope)

#setVoltage(scope, VRange)
time.sleep(1)


#Start measurement
print 'Start pulst-echo test on %s' %(feature)
raw_input('Align the transducers. Press Enter to continue...')
setTime(scope, timeRange, 16e-6)
setVoltage(scope, Range = VRange[0], offset = VOffset[0], source = channels[0], autoRange = False)
setVoltage(scope, Range = VRange[1], offset = VOffset[1], source = channels[1], autoRange = False)
filePath = savePath + feature + '\\'
scope.write(':TRIG:EDGE:LEV 3')
if not os.path.exists(filePath):
    os.mkdir(filePath)
for i in range(2):
    CH = CH_list[i]
    raw_input('Press Enter to collect %s...' %(CH))
    setTime(scope, timeRange, position)
    setVoltage(scope, Range = VRange[i], offset = VOffset[i], source = channels[i], autoRange = False)
    timeCount, waveform = get_waveform(scope,  source = channels[i], points = pointCount, avg = avgCount)
    offset = np.average(waveform[:timeCount/4])
    waveform -= offset
    timeData = get_timeData(timeRange, position, timeCount)
    print len(waveform), timeCount, len(timeData)
    fileName = feature+'_%s.txt' %(CH)
    saveData(timeData, waveform, filePath + fileName)
    print '%s collected' %(CH)
    winsound.Beep(4500,300)
scope.write (':ACQ:COUN 4')
scope.write (':CHAN1:RANG 0.3')
#Data Process
if dataProcess == False:
    exit()
    
fig1 = plt.figure(figsize=(12, 12), dpi=60)
f = open(filePath+'Summary_%s.txt'%(feature),'w')
f.write('Transducer: %s\tStimulus: +/- %sV\tGain: %sdB\tAttenuation: %sdB x2\n' %(feature, drive, gain,  attn)) 
f.write('CH\tVpp(V)\tIL (dB)\tFl(MHz)\tFh(MHz)\tFp(MHz)\tFc(MHz)\tBandwidth(%)\n')   
print feature 
for CH in CH_list:
    #Time domain
    if CH == 'CH0':
        plt.subplot(311)
        fname = filePath + feature + '_%s.txt' %(CH)
        traceT, traceV = loadTraceFile(fname)    
        Vmax = np.max(np.abs(traceV))    
        Vpp = np.max(traceV)-np.min(traceV)
        sens = calcSens(Vpp, drive, gain, attn)
        dt = traceT[1] - traceT[0]
        feq, pw, fl, fh, fp, fc, bw = calcBW(traceV,dt)
        f.write('%s\t%.3f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'%(CH, Vpp, sens, fl, fh, fp, fc,bw))
        print '%.1f\t%.2f'%(sens, bw)
        plt.plot(traceT*1.0E+06,traceV,label='%s' %(CH))
        plt.title(feature+'\nPulse-Echo\nTime (us)')
        plt.xlim(position*1e6-0.5, position*1e6+0.5)
#    plt.ylim(-1.2,1.2)
        plt.legend()   
        plt.ylabel('Amplitude')
    else:    
        plt.subplot(312)
        fname = filePath + feature + '_%s.txt' %(CH)
        traceT, traceV = loadTraceFile(fname)    
        Vmax = np.max(np.abs(traceV))    
        Vpp = np.max(traceV)-np.min(traceV)
        sens = calcSens(Vpp, drive, gain, attn)
        dt = traceT[1] - traceT[0]
        feq, pw, fl, fh, fp, fc, bw = calcBW(traceV,dt)
        f.write('%s\t%.3f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'%(CH, Vpp, sens, fl, fh, fp, fc,bw))
        print '%.1f\t%.2f'%(sens, bw)
        plt.plot(traceT*1.0E+06,traceV,label='%s' %(CH), color = 'g')
        plt.xlim(position*1e6-0.5, position*1e6+0.5)
    #    plt.ylim(-1.2,1.2)
        plt.legend()   
        plt.ylabel('Amplitude')
    #Frequency domain
    plt.subplot(313)
    plt.xlim(10,40)
    plt.ylim(-25,5)
    plt.plot(feq,pw,label='%s' %(CH))
    plt.title('Power Spectrum')
    plt.ylabel('Power (dB)')
    plt.xlabel('Frequency (MHz)')
    plt.legend()

f.close()
plt.savefig(filePath + feature +'_PulseEcho_CH0_Ch7.png')
#plt.show()
winsound.Beep(5000,500)

