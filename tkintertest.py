import Tkinter as tk
from astropy.io import fits
import numpy as np
import copy
from tkFileDialog import askdirectory
import glob
import os

class tkintertest():

    def __init__(self,master):
        self.frame = tk.Frame(master, bg = "#f2f2f2")
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
        self.filemenu.add_command(label = "Open", command = self.openDirectory)
        self.filemenu.add_separator()
        self.filemenu.add_command(label = "Test", command = self.hello)
        
        self.filemenu.add_separator()
        self.filemenu.add_command(label = "Exit", command = root.quit)
        self.menubar.add_cascade(label = "File", menu = self.filemenu)

        
        self.label1 = tk.Label(self.frame, text="Overlap Correction")
        self.label1.grid(row=5)

        self.button_overlapCorrection = tk.Button(
            self.frame, text="Apply",width=10, command=self.hello
            )
        self.button_overlapCorrection.grid(row=5,column=1)


    def openDirectory(self):

        global path
        path = askdirectory()
        return path

    def listFits(self, path):
    
        return glob.glob(os.path.join(path, '*[0-9][0-9][0-9][0-9][0-9].fits'))

    def openFits(self,a,b):
    
        files = self.listFits(path)
        hdua = fits.getheader(files[a])
        hdub = fits.getheader(files[b])
        fileaData = fits.getdata(files[a])
        filebData = fits.getdata(files[b])
        return fileaData, filebData, hdua, hdub

    def readShutter(self,path):
    
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

    def preBinData(self):
    
        shutterData = self.readShutter(path)
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

    def overlapCorrection(self):
    
        a = 0
        b = 1
        binnedData = self.preBinData()
        shutterCounts = self.readShutter(path)
        shutterIndices = binnedData
        shutterValues = shutterCounts[0]
        if not os.path.exists(path+"/overlapCorrected"):
            os.mkdir(path+"/overlapCorrected")
    
        x = 0
        for subIndices in shutterIndices:
            runningTot = np.zeros((512,512),dtype=np.float32)
            while b <= subIndices[-1]:

                fitsFiles = self.openFits(a,b)
                shutters = float(shutterValues[x][1])
                prob = runningTot / shutters
                runningTot += fitsFiles[0]
                correctedPix = (fitsFiles[0] / (1 - prob)).astype(np.int16)
                hdu = fits.PrimaryHDU()
                hdu.data = correctedPix
                hdu.header = fitsFiles[2]
                hdu.header['N_COUNTS'] = sum(sum(correctedPix))
                hdu.writeto(path+'/overlapCorrected/corrected'+str(a)+'.fits')
            
                prob = runningTot / shutters
                runningTot += fitsFiles[1]
                correctedPix = (fitsFiles[1] / (1 - prob)).astype(np.int16)
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

    
    def hello(self):

        self.overlapCorrection()

if __name__ == '__main__':

    root = tk.Tk()
    root.geometry("300x300")
    root.title("Tkinter test")
    app = tkintertest(root)
    root.mainloop()
    root.destroy()
