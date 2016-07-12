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

"""
this function allows the user to select the directory containing their data.
parameters: None
Returns: working directory
"""
def openDirectory():
    name = askdirectory()
    return name

path = openDirectory()
    
"""
This function finds all .fits files in directory chosen by openDirectory()
parameters: None
Returns: List of filenames contained within the directory
"""
def listFits(path):
    #path = openDirectory()
    files = glob.glob(os.path.join(path, '*[0-9][0-9][0-9][0-9][0-9].fits'))
    return files

"""
open two consecutive files at a time as it's too difficult to open them all
parameters: file indices
Returns: opened .fits files
"""
def openFits(a,b):
    files = listFits(path)
    filea = fits.open(files[a])
    fileahdu = filea[0]
    fileaData = fileahdu.data
    fileb = fits.open(files[b])
    filebhdu = fileb[0]
    filebData = filebhdu.data
    return fileaData, filebData

"""
reads shutter count/time from text files in chosen directory
parameters: path to directory
Returns: shutter counts and shutter times
"""
def readShutter(path):
    countFile = glob.glob(os.path.join(path, "*ShutterCount.txt"))
    timeFile = glob.glob(os.path.join(path, "*ShutterTimes.txt"))
    readCount = open(str(countFile[0]))
    readTime = open(str(timeFile[0]))
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
    return countData, timeData
    
"""
parameters: shutters acquired (float), sum (float)
returns: probability of occupied pixel (float)
"""
def overlapProbability(shutters, sum):
    probability = sum / shutters
    return probability

"""
paramters: probability (float), current pixel value (float)
returns: pixel intensity corrected for overlap (float)
"""
def correctedCount(probability, value):
    correctedPixel = value / (1 - probability)
    return correctedPixel

"""
this function should handle all the previous functions in order to apply the overlap correction to the data.
"""
def overlapCorrection():
    a = 0
    b = 1
    runningTot = np.zeros((512,512),dtype=int)
    
    while b < len(listFits(path)):
        fitsFiles = openFits(a,b)
        runningTot = runningTot + fitsFiles[0]
        i = 0
        for row in fitsFiles[1]:
            j = 0
            for item in row:
                runningTot[i,j] = runningTot[i,j] + item
                j = j + 1
            i = i + 1
        print b

        a = a + 1
        b = b + 1
    print runningTot
                
overlapCorrection()
