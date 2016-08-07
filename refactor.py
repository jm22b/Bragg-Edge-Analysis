import Tkinter as tk
import glob
import os
import numpy as np
from tkFileDialog import askdirectory
from astropy.io import fits


class BraggEdgeAnalysisGUI:

    def __init__(self, root_):

        self.root = root_
        self.frame = tk.Frame(self.root)
        self.frame.pack()

        self.directory = GetDirectories()
        self.test = Test(self.directory)
        self.correction = OverlapCorrectionAndScaling(self.directory)

        self.menubar = tk.Menu(self.root)
        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        self.actionmenu = tk.Menu(self.menubar, tearoff=0)

        self.openButton = tk.Button(self.frame, text="open", command=self.directory.getOpenPath)
        self.sampleButton = tk.Button(self.frame, text="sample", command=self.directory.getSamplePath)

        self.testButton = tk.Button(self.frame, text="test", command= self.test.do)

        self.widgets()

    def widgets(self):

        root.option_add("*tearoff", "FALSE")

        self.filemenu.add_command(
            label="Load Open Beam", command=self.directory.getOpenPath)
        self.filemenu.add_separator()
        self.filemenu.add_command(
            label="Load Sample Beam", command=self.directory.getSamplePath)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=root.destroy)
        self.menubar.add_cascade(label="File", menu=self.filemenu)

        self.actionmenu.add_command(label="Correct & Scale Data", command=self.correction.doBoth)
        self.menubar.add_cascade(label="Actions", menu=self.actionmenu)

        root.config(menu=self.menubar)
        
        self.openButton.pack()
        self.sampleButton.pack()
        self.testButton.pack()


class FitsData:

    def __init__(self, names=None, headers=None, arrays=None):

        if names == None:
            names = []
        if headers == None:
            headers = []
        if arrays == None:
            arrays = []

        self.names = names
        self.headers = headers
        self.arrays = arrays
        

class DirectoryHandler:

    def __init__(self):
        self.openPath = None
        self.samplePath = None

    def openOpenDirectory(self):
        self.openPath = askdirectory()
        return self.openPath

    def openSampleDirectory(self):
        self.samplePath = askdirectory()
        return self.samplePath

class GetDirectories:

    def __init__(self):

        self.directory = DirectoryHandler()
        self.openFits = FitsData()
        self.sampleFits = FitsData()
        
        self.openPath = None
        self.samplePath = None

    def getOpenPath(self):
        
        self.directory.openOpenDirectory()
        self.openPath = self.directory.openPath
        if os.path.exists(os.path.join(self.openPath, "scaledOpenBeam")):
            path = os.path.join(self.openPath, "scaledOpenBeam")
            self.loadData(path, self.openFits)
        else:
            self.loadData(self.openPath, self.openFits)
        
    def getSamplePath(self):

        self.directory.openSampleDirectory()
        self.samplePath = self.directory.samplePath
        if os.path.exists(os.path.join(self.samplePath, "overlapCorrected")):
            path = os.path.join(self.samplePath, "overlapCorrected")
            self.loadData(path, self.sampleFits)
        else:
            self.loadData(self.samplePath, self.sampleFits)

    def loadData(self, path, container):
        f = glob.glob(os.path.join(path, "*[0-9][0-9][0-9][0-9][0-9].fits"))
        for fitsFile in f:
            hdulist = fits.open(fitsFile, memmap=False)
            name = hdulist.filename().split("\\")[-1]
            header = hdulist[0].header
            data = hdulist[0].data
            hdulist.close()
            container.names.append(name)
            container.headers.append(header)
            container.arrays.append(data)
            print name

            
class OverlapCorrectionAndScaling:

    def __init__(self, directory):

        self.directory = directory
        self.openPath = self.directory.openPath
        self.samplePath = self.directory.samplePath
        
        self.sampleFits = self.directory.sampleFits
        self.openFits = self.directory.openFits

        #self.openData = openData
        #self.sampleData = sampleData

    def readShutter(self, path):
        # finds the ShutterCount file in openbeam folder
        countFile = glob.glob(os.path.join(path, "*ShutterCount.txt"))
        print countFile
        # finds the ShutterTimes file in openbeam folder
        timeFile = glob.glob(os.path.join(path, "*ShutterTimes.txt"))
        # finds the spectra file in openbeam folder
        spectraFile = glob.glob(os.path.join(path, "*Spectra.txt"))
        # opens the above files
        readCount = open(str(countFile[0]))
        readTime = open(str(timeFile[0]))
        readSpectra = open(str(spectraFile[0]))
        countData = []
        for line in readCount:
            counts = line.split()  # splits the line into two parts on the empty space
            if counts[1] == '0':  # if the second part of the line is zero, stop
                break
            countData.append(counts)  # insert shutter counts into countData
        timeData = []
        for line in readTime:  # same as above
            time = line.split()
            if time[2] == '0':
                break
            timeData.append(time)
        spectraData = []
        for line in readSpectra:  # same as above but no need to check for zeros
            spectra = line.split()
            spectraData.append(spectra)

        return countData, timeData, spectraData

    def preBinData(self, path):

        shutterData = self.readShutter(path)  # calls the read shutter function
        shutterTimes = []
        x = 0
        for i in shutterData[0]:  # computes the shutter intervals
            time = float(shutterData[1][x][1]) + float(shutterData[1][x][2])
            shutterTimes.append(time)
            x += 1
        x = 0
        for item in shutterTimes:
            if x == 0:
                shutterTimes[x] += 0
            else:
                shutterTimes[x] += shutterTimes[x - 1]
            x += 1
        number_of_shutters = len(shutterData[0])
        shutterIndices = []
        for i in range(0, number_of_shutters):  # constructs a list of lists containing how to map an image slice
            # to a shutter value. e.g if '700' is contained within the second
            # list, then the second shutter value will be used for correction
            shutterIndices.insert(0, [])
        x = 0
        for i in range(0, number_of_shutters):
            for line in shutterData[2][x:]:
                if float(line[0]) < shutterTimes[i]:
                    shutterIndices[i].append(x)
                    x += 1
                else:
                    break
        return shutterIndices

    def overlapCorrection(self, path):

        if path == self.directory.samplePath:

            shutterIndices = self.preBinData(path)
            shutterValues = self.readShutter(path)[0]

            if os.path.exists(path + "/overlapCorrected"):

                return ctypes.windll.user32.MessageBoxA(0, "Corrected files already exist", "Error", 1)

            os.mkdir(path + "/overlapCorrected")
            f = open(path + "/overlapCorrected/TOFData.csv", "wb")
            zipped = zip(self.sampleFits.arrays, self.sampleFits.headers, self.sampleFits.names)
            s = 0
            for subIndex in shutterIndices:
                subList = zipped[subIndex[0]:subIndex[-1]+1]
                runningTot = np.zeros((512, 512))

                for data, header, name in subList:

                    shutter = float(shutterValues[s][1])
                    prob = np.divide(runningTot, shutter)
                    runningTot += data
                    correction = np.round(np.divide(data, (1 - prob))).astype(np.int16)

                    # hdu = fits.PrimaryHDU()
                    # hdu.data = correction
                    # hdu.header = header
                    # counts = sum(sum(correction))
                    # hdu.header["N_COUNTS"] = counts
                    # TOF = hdu.header["TOF"]

                    # line = "%.16f, %d\n" % (TOF, counts)
                    # f.writelines(line)

                    #hdu.writeto(path + "/overlapCorrected/corrected"+name)
                    print name
                print s
                s += 1
            f.close()

        else:
            shutterIndices = self.preBinData(path)
            shutterValuesOpen = self.readShutter(path)[0]
            shutterValuesSample = self.readShutter(path)[0]

            zipShutters = zip(shutterValuesOpen, shutterValuesSample)
            ratio = []
            for svo, svs in zipShutters:
                ratio.append(float(svs[1]) / float(svo[1]))

            if os.path.exists(path + "/scaledOpenBeam"):
                return ctypes.windll.user32.MessageBoxA(0, "Scaled files already exist", "Error", 1)

            os.mkdir(path + "/scaledOpenBeam")
            f = open(path + "/scaledOpenBeam/TOFData.csv", "wb")
            zipped = zip(self.openFits.arrays, self.openFits.headers, self.openFits.names)
            s = 0

            for subIndex in shutterIndices:
                sublist = zipped[subIndex[0]:subIndex[-1]+1]
                runningTot = np.zeros((512, 512))
                scaleFactor = ratio[s]

                for data, header, name in sublist:
                    shutter = float(shutterValuesOpen[s][1])
                    prob = np.divide(runningTot, shutter)
                    runningTot += data
                    correction = np.round(np.divide(data, (1 - prob))).astype(np.int16)
                    scaled = correction * scaleFactor

                    # hdu = fits.PrimaryHDU()
                    # hdu.data = scaled
                    # hdu.header = header
                    # counts = sum(sum(scaled))
                    # hdu.header["N_COUNTS"] = counts
                    # TOF = hdu.header["TOF"]

                    # line = "%.16f, %d\n" % (TOF, counts)
                    # f.writelines(line)

                    # hdu.writeto(path + "/scaledOpenBeam/scaled" + name)
                    print name
                print s
                s += 1
            f.close()

    def doBoth(self):
        self.overlapCorrection(self.directory.openPath)
        self.overlapCorrection(self.directory.samplePath)

    def do(self):
        pass

class Test:

    def __init__(self, dir):
        self.a = dir

    def do(self):
        print self.a.openFits.arrays

    

if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("500x600")
    root.title("Bragg Edge Analysis")
    app = BraggEdgeAnalysisGUI(root)
    root.mainloop()
    root.destroy()
