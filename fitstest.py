from astropy.io import fits
import numpy as np
import copy
import Tkinter as Tk
from tkFileDialog import askdirectory
import glob
import os
import time

def openDirectory():
    return askdirectory()

pathToOpen = openDirectory()
pathToSample = openDirectory()

def listFits(path):
    return glob.glob(os.path.join(path, '*[0-9][0-9][0-9][0-9][0-9].fits'))

def openFits(a):

    filesOpen = listFits(pathToOpen)
    filesSample = listFits(pathToSample)

    hduOpena = fits.getheader(filesOpen[a])
    fileOpenaData = fits.getdata(filesOpen[a])

    hduSamplea = fits.getheader(filesSample[a])
    fileSampleaData = fits.getdata(filesSample[a])

    return hduOpena, fileOpenaData, hduSamplea, fileSampleaData

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

def preBinData():
    
    shutterData = readShutter(pathToSample)
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

def overlapCorrection():
    
    a = 0
    binnedData = preBinData()
    shutterCounts = readShutter(pathToSample)
    shutterIndices = binnedData
    shutterValues = shutterCounts[0]
    if not os.path.exists(pathToOpen+"/1overlapCorrected"):
        os.mkdir(pathToOpen+"/1overlapCorrected")
    if not os.path.exists(pathToSample+"/1overlapCorrected"):
        os.mkdir(pathToSample+"/1overlapCorrected")
    
    x = 0
    for subIndices in shutterIndices:
        runningTotOpen = np.zeros((512,512),dtype=np.float32)
        runningTotSample = np.zeros((512,512),dtype=np.float32)
        while a <= subIndices[-1]:

            fitsFiles = openFits(a)
            shutters = float(shutterValues[x][1])
            prob = runningTotOpen / shutters
            runningTotOpen += fitsFiles[1]
            correctedPix = (fitsFiles[1] / (1 - prob)).astype(np.int16)
            hdu = fits.PrimaryHDU()
            hdu.data = correctedPix
            hdu.header = fitsFiles[0]
            hdu.header['N_COUNTS'] = sum(sum(correctedPix))
            hdu.writeto(pathToOpen+'/1overlapCorrected/corrected'+str(a)+'.fits')
            
            prob = runningTotSample / shutters
            runningTotSample += fitsFiles[3]
            correctedPix = (fitsFiles[3] / (1 - prob)).astype(np.int16)
            hdu = fits.PrimaryHDU()
            hdu.data = correctedPix
            hdu.header = fitsFiles[2]
            hdu.header['N_COUNTS'] = sum(sum(correctedPix))
            hdu.writeto(pathToSample+'/1overlapCorrected/corrected'+str(a)+'.fits')
            a += 1
            print a
        print shutters
        x += 1

overlapCorrection()
