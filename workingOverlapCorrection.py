"""
Overlap correction functions - probability and correction

These functions work in conjunction to correct for the fact that there is an overlap between detection of events and reading the data from a Timepix pixel. Until the frame has been read out, the pixel cannot detect further events, and so we calculate the probability of a pixel being occupied and use this to produce a corrected count number.
"""

from astropy.io import fits
import numpy as np
import copy
import Tkinter as Tk
from tkFileDialog import askdirectory
import glob
import os
import time



"""
this function allows the user to select the directory containing their data.
parameters: None
Returns: working directory
"""
#@profile
def openDirectory():
    
    return askdirectory()

path = openDirectory()
    
"""
This function finds all .fits files in directory chosen by openDirectory()
parameters: None
Returns: List of filenames contained within the directory
"""
#@profile
def listFits(path):
    
    return glob.glob(os.path.join(path, '*[0-9][0-9][0-9][0-9][0-9].fits'))

"""
open two consecutive files at a time as it's too difficult to open them all
parameters: file indices
Returns: opened .fits files data and headers
"""
#@profile
def openFits(a,b):
    
    files = listFits(path)
    hdua = fits.getheader(files[a])
    hdub = fits.getheader(files[b])
    fileaData = fits.getdata(files[a])
    filebData = fits.getdata(files[b])
    return fileaData, filebData, hdua, hdub

"""
reads shutter count/time from text files in chosen directory
parameters: path to directory
Returns: shutter counts, shutter times and TOF spectra
"""
#@profile
def readShutter(path):
    
    countFile = glob.glob(os.path.join(path, "*ShutterCount.txt"))
    timeFile = glob.glob(os.path.join(path, "*ShutterTimes.txt"))
    spectraFile = glob.glob(os.path.join(path, "*Spectra.txt"))
    readCount = open(str(countFile[0]))
    readTime = open(str(timeFile[0]))
    readSpectra = open(str(spectraFile[0]))
    countData = []
    for line in readCount:
        counts = line.split()
        if counts[1] == '0':
            break
        countData.append(counts)
    timeData = []
    for line in readTime:
        time = line.split()
        if time[2] == '0':
            break
        timeData.append(time)
    spectraData = []
    for line in readSpectra:
        spectra = line.split()
        spectraData.append(spectra)
        
        
    return countData, timeData, spectraData

"""
this function predetermines the data that belongs to each shutter, by comparing TOF data with shutter times. The ranges computed are then used to compute the necessary corrections
Paramters: None
Returns: shutter indices
"""
#@profile
def preBinData():
    
    shutterData = readShutter(path)
    shutterTimes = []
    x = 0
    for i in shutterData[0]:
        time = float(shutterData[1][x][1]) + float(shutterData[1][x][2])
        shutterTimes.append(time)
        x += 1
    x = 0
    for item in shutterTimes:
        if x == 0:
            shutterTimes[x] += 0
        else:
            shutterTimes[x] += shutterTimes[x-1]
        x += 1
    number_of_shutters = len(shutterData[0])
    shutterIndices = []
    for i in range(0,number_of_shutters):
        shutterIndices.insert(0,[])
    x = 0
    for i in range(0,number_of_shutters):
        for line in shutterData[2][x:]:
            if float(line[0]) < shutterTimes[i]:
                shutterIndices[i].append(x)
                x += 1
            else:
                break
    return shutterIndices

t0 = float(time.time())

"""
uses prebinned data to map each slice to a specific shutter and apply the correction algorithms for the MCP detector from page 3 of "optimisation of timepix count rate capabilities for the applications with a periodic input signal" - AS Tremsin et al, with the minor modification that the summation runs from i = 0 to i = j-1

This function will check to see if the overlap corrected data already exists, and if it does not, it will write the corrected .fits files to a new directory inside the directory containing the uncorrected data.

parameters: None
returns: None
"""
def overlapCorrection():
    
    a = 0
    b = 1
    binnedData = preBinData()
    shutterCounts = readShutter(path)
    shutterIndices = binnedData
    shutterValues = shutterCounts[0]
    if not os.path.exists(path+"/overlapCorrected"):
        os.mkdir(path+"/overlapCorrected")
    
    x = 0
    for subIndices in shutterIndices:
        runningTot = np.zeros((512,512),dtype=np.float16)
        while b <= subIndices[-1]:

            fitsFiles = openFits(a,b)
            shutters = float(shutterValues[x][1])
            prob = runningTot / shutters
            runningTot += fitsFiles[0]
            correctedPix = np.round(fitsFiles[0] / (1 - prob)).astype(np.int16)
            hdu = fits.PrimaryHDU()
            hdu.data = correctedPix
            hdu.header = fitsFiles[2]
            hdu.header['N_COUNTS'] = sum(sum(correctedPix))
            hdu.writeto(path+'/overlapCorrected/corrected'+str(a)+'.fits')
            
            prob = runningTot / shutters
            runningTot += fitsFiles[1]
            correctedPix = np.round(fitsFiles[1] / (1 - prob)).astype(np.int16)
            hdu = fits.PrimaryHDU()
            hdu.data = correctedPix
            hdu.header = fitsFiles[3]
            hdu.header['N_COUNTS'] = sum(sum(correctedPix))
            hdu.writeto(path+'/overlapCorrected/corrected'+str(b)+'.fits')
            a += 2
            b += 2
            print b
        print shutters
        x += 1

overlapCorrection()

t1 = float(time.time())
timeElapsed = (t1 - t0) / 60.0
print "%f minutes" % timeElapsed
