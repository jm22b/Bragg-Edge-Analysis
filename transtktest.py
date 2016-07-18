import Tkinter as tk
from astropy.io import fits
import numpy as np
import copy
from tkFileDialog import askdirectory
import glob
import os
import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot as plt

class tkintertest():

    def __init__(self,master):
        self.frame = tk.Frame(master)
        self.frame.pack()

        msg = tk.Label(root, text = "tkinter test")
        msg.pack()

        self.widgets()

    def widgets(self):
        root.option_add("*tearoff", "FALSE")
        self.menubar = tk.Menu(root)
        self.menubar.add_command(label = "Quit", command = root.quit())
        root.config(menu = self.menubar)

        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        self.filemenu.add_separator()
        self.filemenu.add_command(label = "Load Open Beam", command = self.openOpenDirectory)
        self.filemenu.add_separator()
        self.filemenu.add_command(label = "Load Sample", command = self.openSampleDirectory)
        
        self.filemenu.add_separator()
        self.filemenu.add_command(label = "Exit", command = root.quit)
        self.menubar.add_cascade(label = "File", menu = self.filemenu)

        
        self.label1 = tk.Label(self.frame, text="Overlap Correction (open beam)")
        self.label1.grid(row=1)

        self.buttonOverlapCorrectionOpen = tk.Button(
            self.frame, text="Apply",width=10, command=self.overlapOpenCorrection
            )
        self.buttonOverlapCorrectionOpen.grid(row=1,column=1)

        self.label2 = tk.Label(self.frame, text="Overlap Correction (sample)")
        self.label2.grid(row=2)
        self.buttonOverlapCorrectionSample = tk.Button(
            self.frame, text="Apply", width=10, command=self.overlapSampleCorrection)
        self.buttonOverlapCorrectionSample.grid(row=2,column=1)
        
        self.label3 = tk.Label(self.frame, text="Transmission Plot")
        self.label3.grid(row=3)
        self.buttonTransmissionPlot = tk.Button(
            self.frame, text="Plot", width=10, command=self.plot)
        self.buttonTransmissionPlot.grid(row=3,column=1)


    def openOpenDirectory(self):

        global openPath
        openPath = askdirectory()
        return openPath

    def openSampleDirectory(self):

        global samplePath
        samplePath = askdirectory()
        return samplePath

    def listFits(self, path):
    
        return glob.glob(os.path.join(path, '*[0-9][0-9][0-9][0-9][0-9].fits'))

    def openOpenFits(self,a,b):
    
        files = self.listFits(openPath)
        hdua = fits.getheader(files[a])
        hdub = fits.getheader(files[b])
        fileaData = fits.getdata(files[a])
        filebData = fits.getdata(files[b])
        return fileaData, filebData, hdua, hdub

    def openSampleFits(self,a,b):

        files = self.listFits(samplePath)
        hdua = fits.getheader(files[a])
        hdub = fits.getheader(files[b])
        fileaData = fits.getdata(files[a])
        filebData = fits.getdata(files[b])
        return fileaData, filebData, hdua, hdub

    def readOpenShutter(self,path):
    
        countFile = glob.glob(os.path.join(openPath, "*ShutterCount.txt"))
        timeFile = glob.glob(os.path.join(openPath, "*ShutterTimes.txt"))
        spectraFile = glob.glob(os.path.join(openPath, "*Spectra.txt"))
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

    def readSampleShutter(self,path):
    
        countFile = glob.glob(os.path.join(samplePath, "*ShutterCount.txt"))
        timeFile = glob.glob(os.path.join(samplePath, "*ShutterTimes.txt"))
        spectraFile = glob.glob(os.path.join(samplePath, "*Spectra.txt"))
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

    def preBinOpenData(self):
    
        shutterData = self.readOpenShutter(openPath)
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

    def preBinSampleData(self):
    
        shutterData = self.readSampleShutter(samplePath)
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

    def overlapOpenCorrection(self):
    
        a = 0
        b = 1
        binnedData = self.preBinOpenData()
        shutterCounts = self.readOpenShutter(openPath)
        shutterIndices = binnedData
        shutterValues = shutterCounts[0]
        if not os.path.exists(openPath+"/overlapCorrected"):
            os.mkdir(openPath+"/overlapCorrected")
        f = open(openPath+'/overlapCorrected/TOFData','wb+')
        x = 0
        for subIndices in shutterIndices:
            runningTot = np.zeros((512,512),dtype=np.float32)
            while b <= subIndices[-1]:

                fitsFiles = self.openOpenFits(a,b)
                shutters = float(shutterValues[x][1])
                
                prob = runningTot / shutters
                runningTot += fitsFiles[0]
                correctedPix = np.round(fitsFiles[0] / (1 - prob)).astype(np.int16)
                hdu = fits.PrimaryHDU()
                hdu.data = correctedPix
                hdu.header = fitsFiles[2]
                hdu.header['N_COUNTS'] = sum(sum(correctedPix))
                counts = hdu.header['N_COUNTS']
                TOF = hdu.header['TOF']
                line = '%.16f %d \n' % (TOF, counts)
                f.write(line)
                hdu.writeto(openPath+'/overlapCorrected/corrected'+str(a)+'.fits')
            
                prob = runningTot / shutters
                runningTot += fitsFiles[1]
                correctedPix = np.round(fitsFiles[1] / (1 - prob)).astype(np.int16)
                hdu = fits.PrimaryHDU()
                hdu.data = correctedPix
                hdu.header = fitsFiles[3]
                hdu.header['N_COUNTS'] = sum(sum(correctedPix))
                counts = hdu.header['N_COUNTS']
                TOF = hdu.header['TOF']
                line = '%.16f %d \n' % (TOF, counts)
                f.write(line)
                hdu.writeto(openPath+'/overlapCorrected/corrected'+str(b)+'.fits')
                a += 2
                b += 2
                print b
            print shutters
            x += 1
        f.close()

    def overlapSampleCorrection(self):
    
        a = 0
        b = 1
        binnedData = self.preBinSampleData()
        shutterCounts = self.readSampleShutter(samplePath)
        shutterIndices = binnedData
        shutterValues = shutterCounts[0]
        if not os.path.exists(samplePath+"/overlapCorrected"):
            os.mkdir(samplePath+"/overlapCorrected")
        f = open(samplePath+'/overlapCorrected/TOFData.txt','ab+')
        x = 0
        for subIndices in shutterIndices:
            runningTot = np.zeros((512,512),dtype=np.float32)
            while b <= subIndices[-1]:

                fitsFiles = self.openSampleFits(a,b)
                shutters = float(shutterValues[x][1])
                
                prob = runningTot / shutters
                runningTot += fitsFiles[0]
                correctedPix = np.round(fitsFiles[0] / (1 - prob)).astype(np.int16)
                hdu = fits.PrimaryHDU()
                hdu.data = correctedPix
                hdu.header = fitsFiles[2]
                hdu.header['N_COUNTS'] = sum(sum(correctedPix))
                counts = hdu.header['N_COUNTS']
                TOF = hdu.header['TOF']
                line = '%.16f %d \n' % (TOF, counts)
                f.write(line)
                hdu.writeto(samplePath+'/overlapCorrected/corrected'+str(a)+'.fits')
            
                prob = runningTot / shutters
                runningTot += fitsFiles[1]
                correctedPix = np.round(fitsFiles[1] / (1 - prob)).astype(np.int16)
                hdu = fits.PrimaryHDU()
                hdu.data = correctedPix
                hdu.header = fitsFiles[3]
                hdu.header['N_COUNTS'] = sum(sum(correctedPix))
                counts = hdu.header['N_COUNTS']
                TOF = hdu.header['TOF']
                line = '%.16f %d \n' % (TOF, counts)
                f.write(line)
                hdu.writeto(samplePath+'/overlapCorrected/corrected'+str(b)+'.fits')
                a += 2
                b += 2
                print b
            print shutters
            x += 1
        f.close()
    
    def getData(self):
        
        f = open('/Users/jacobmaresca/python/python2/BraggEdgeAnalysis/IMATdatajuly2016/OpenBeam/overlapCorrected/TOFData.txt', 'ab+')
        TOF = []
        openData = []
        for line in f:
            a = line.split()
            openData.append(a)
            TOF.append(a[0])
        g = open('/Users/jacobmaresca/python/python2/BraggEdgeAnalysis/IMATdatajuly2016/Weld1/Data/overlapCorrected/TOFData.txt', 'ab+')
        sampleData = []
        for line in g:
            b = line.split()
            sampleData.append(b)

        zipped = zip(openData, sampleData)
        transmissionData = []
        for x,y in zipped:
            transmitted = float(y[1])/float(x[1])
            TOF = float(x[0])
            transmissionData.insert(-1,(transmitted,TOF))

        f.close()
        g.close()

        return transmissionData
    
    
    def plot(self):

        data = self.getData()
        y_list = [x[0] for x in data]
        x_list = [x[1] for x in data]
    
        plt.plot(x_list, y_list,'x')
        plt.xlabel('Time of Flight (micro seconds)')
        plt.ylabel('Transmission')
        plt.title('Neutron Transmission')
        plt.show()
        plt.close()
    
    
    
    def hello(self):

        self.overlapCorrection()

if __name__ == '__main__':

    root = tk.Tk()
    root.geometry("500x500")
    root.title("Tkinter test")
    app = tkintertest(root)
    root.mainloop()
    root.destroy()
