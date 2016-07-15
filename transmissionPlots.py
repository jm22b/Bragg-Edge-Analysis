import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits
import Tkinter as tk
from tkFileDialog import askdirectory
import glob
import os

def openSample():
    return askdirectory()


def listSampleFiles(pathToSamples):
    return glob.glob(os.path.join(pathToSamples, '*[0-9][0-9][0-9][0-9][0-9].fits'))

def openOpenBeam():
    pathToOpenBeam = askdirectory()
    return pathToOpenBeam

def listOpenBeamFiles(pathToOpenBeam):
    return glob.glob(os.path.join(pathToOpenBeam, '*[0-9][0-9][0-9][0-9][0-9].fits'))

def compute():
    pathToSample = openSample()
    pathToOpenBeam = openOpenBeam()

    sampleFiles = listSampleFiles(pathToSample)
    openBeamFiles = listOpenBeamFiles(pathToOpenBeam)

    zipped = zip(sampleFiles,openBeamFiles)

    transmissionData = []

    for x,y in zipped:
        hdux = fits.getheader(x)
        N_COUNTSsample = hdux['N_COUNTS']
        TOF = hdux['TOF']
        hduy = fits.getheader(y)
        N_COUNTSopen = hduy['N_COUNTS']
        transmitted = float(N_COUNTSsample)/float(N_COUNTSopen)
        transmissionData.insert(-1,(transmitted,TOF))

    return transmissionData

def plot():
    data = compute()
    y_list = [x[0] for x in data]
    x_list = [x[1] for x in data]
    x_list = sorted(x_list)
    
    plt.plot(x_list, y_list,'-')
    plt.show()
    plt.close()


plot()
