import Tkinter as tk
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.widgets import Slider, RectangleSelector
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
        self.transplot = tk.Menu(self.menubar)

        self.flightpath = tk.Entry(self.frame, width=30)

        self.showDataButton = tk.Button(self.frame, text="Show Sample", width=10, command=lambda: ShowData(self.root, self.directory).plot())

        self.testButton = tk.Button(self.frame, text="test", command=lambda: Test(self.directory).do())

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
        self.actionmenu.add_separator()
        self.actionmenu.add_cascade(label="Transmission", menu=self.transplot)
        self.transplot.add_command(label="Plot (TOF)", command=lambda: TransPlot(self.directory, self.flightpath).plotTransTOF())
        self.transplot.add_separator()
        self.transplot.add_command(label="Plot (Wavelength)", command=lambda: TransPlot(self.directory, self.flightpath).plotTransWavelength())
        self.menubar.add_cascade(label="Actions", menu=self.actionmenu)

        root.config(menu=self.menubar)

        self.flightpath.insert(0, "Default flight path: 56m")
        self.flightpath.pack()
        self.showDataButton.pack()
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
        
        #self.openPath = None
        #self.samplePath = None

    def getOpenPath(self):
        
        self.directory.openOpenDirectory()
        #self.openPath = self.directory.openPath
        if os.path.exists(os.path.join(self.directory.openPath, "scaledOpenBeam")):
            path = os.path.join(self.directory.openPath, "scaledOpenBeam")
            self.loadData(path, self.openFits)
        else:
            self.loadData(self.directory.openPath, self.openFits)
        
    def getSamplePath(self):

        self.directory.openSampleDirectory()
        #self.samplePath = self.directory.samplePath
        if os.path.exists(os.path.join(self.directory.samplePath, "overlapCorrected")):
            path = os.path.join(self.directory.samplePath, "overlapCorrected")
            self.loadData(path, self.sampleFits)
        else:
            self.loadData(self.directory.samplePath, self.sampleFits)

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
        #self.openPath = self.directory.openPath
        #self.samplePath = self.directory.samplePath

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
            #zipped = zip(self.sampleFits.arrays, self.sampleFits.headers, self.sampleFits.names)
            s = 0
            for subIndex in shutterIndices:
                #subList = zipped[subIndex[0]:subIndex[-1]+1]
                runningTot = np.zeros((512, 512))

                for i in range(subIndex[0], subIndex[0]+len(subIndex)):
                    shutter = float(shutterValues[s][1])
                    prob = np.divide(runningTot, shutter)
                    runningTot += self.directory.sampleFits.arrays[i]
                    self.directory.sampleFits.arrays[i] = np.round(np.divide(self.directory.sampleFits.arrays[i], (1 - prob))).astype(np.int16)
                    #hdu = fits.PrimaryHDU()
                    #hdu.data = self.directory.sampleFits.arrays[i]
                    #counts = sum(sum(self.directory.sampleFits.arrays[i]))
                    #hdu.header = self.directory.sampleFits.headers[i]
                    #hdu.header["N_COUNTS"] = counts
                    #TOF = hdu.header["TOF"]
                    
                    #line = "%.16f, %d\n" % (TOF, counts)
                    #f.writelines(line)

                    #hdu.writeto(path + "/overlapCorrected/corrected"+self.directory.sampleFits.names[i])
                    
            

                #for data, header, name in subList:

                    #shutter = float(shutterValues[s][1])
                    #prob = np.divide(runningTot, shutter)
                    #runningTot += data
                    #correction = np.round(np.divide(data, (1 - prob))).astype(np.int16)

                    # hdu = fits.PrimaryHDU()
                    # hdu.data = correction
                    # hdu.header = header
                    # counts = sum(sum(correction))
                    # hdu.header["N_COUNTS"] = counts
                    # TOF = hdu.header["TOF"]

                    # line = "%.16f, %d\n" % (TOF, counts)
                    # f.writelines(line)

                    #hdu.writeto(path + "/overlapCorrected/corrected"+name)
                    #print name
                print s
                s += 1
            f.close()
            print self.directory.sampleFits.arrays[100]

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
            #zipped = zip(self.openFits.arrays, self.openFits.headers, self.openFits.names)
            s = 0

            for subIndex in shutterIndices:
                #sublist = zipped[subIndex[0]:subIndex[-1]+1]
                runningTot = np.zeros((512, 512))
                scaleFactor = ratio[s]

                for i in range(subIndex[0], subIndex[0]+len(subIndex)):
                    shutter = float(shutterValuesOpen[s][1])
                    prob = np.divide(runningTot, shutter)
                    runningTot += self.directory.openFits.arrays[i]
                    self.directory.openFits.arrays[i] = np.round(np.divide(self.directory.openFits.arrays[i], (1 - prob))).astype(np.int16) * scaleFactor
                    
                    #hdu = fits.PrimaryHDU()
                    #hdu.data = self.directory.openFits.arrays[i]
                    #counts = sum(sum(self.directory.openFits.arrays[i]))
                    #hdu.header = self.directory.openFits.headers[i]
                    #hdu.header["N_COUNTS"] = counts
                    #TOF = hdu.header["TOF"]
                    
                    #line = "%.16f, %d\n" % (TOF, counts)
                    #f.writelines(line)

                    #hdu.writeto(path + "/scaledOpenBeam/scaled"+self.directory.sampleFits.names[i])
                   

                #for data, header, name in sublist:
                    #shutter = float(shutterValuesOpen[s][1])
                    #prob = np.divide(runningTot, shutter)
                    #runningTot += data
                    #correction = np.round(np.divide(data, (1 - prob))).astype(np.int16)
                    #scaled = correction * scaleFactor

                    # hdu = fits.PrimaryHDU()
                    # hdu.data = scaled
                    # hdu.header = header
                    # counts = sum(sum(scaled))
                    # hdu.header["N_COUNTS"] = counts
                    # TOF = hdu.header["TOF"]

                    # line = "%.16f, %d\n" % (TOF, counts)
                    # f.writelines(line)

                    # hdu.writeto(path + "/scaledOpenBeam/scaled" + name)
                    #print name
                print s
                s += 1
            f.close()
            print self.directory.openFits.arrays[100] 

    def doBoth(self):
        self.overlapCorrection(self.directory.openPath)
        self.overlapCorrection(self.directory.samplePath)

    
class ShowData:

    def __init__(self, root, directory):
        
        self.root = root
        self.directory = directory
        self.fig = Figure(figsize=(6, 6))
        self.ax = self.fig.add_subplot(111)
        self.plotted = False
        self.l = None
        self.canvas = None
        plt.show()

    def onSelect(self, eclick, erelease):
        print "Start position: (%f, %f)" % (eclick.xdata, eclick.ydata)
        print "End position: (%f, %f)" % (erelease.xdata, erelease.ydata)
        global a
        a = eclick.xdata
        global b
        b = erelease.xdata
        global c
        c = eclick.ydata
        global d
        d = erelease.ydata
        return a, b, c, d

    def plot(self, **kwargs):
        self.slider = tk.Scale(
            self.root, from_=0, to=len(self.directory.sampleFits.arrays)-1, resolution=1, orient=tk.HORIZONTAL, command=self.update
            )
        self.slider.pack()
        self.plotted = True
        self.s = 0
        im = self.histeq(self.directory.sampleFits.arrays[self.s])[0]
        self.l = self.ax.imshow(im, cmap=plt.cm.gray, **kwargs)
        self.canvas = FigureCanvasTkAgg(self.fig, self.root)
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()
        self.myrectsel = MyRectangleSelector(self.ax, self.onSelect, drawtype="box", rectprops=dict(
            facecolor="red", edgecolor="black", alpha=0.2, fill=True))

    def update(self, val):
        ind = int(self.slider.get())
        if self.plotted:
            im = self.histeq(self.directory.sampleFits.arrays[ind])[0]
            self.l.set_data(im)
            self.canvas.draw()

    def histeq(self, im, nbr_bins=256):
        # get image histogram
        imhist, bins = np.histogram(im.flatten(), nbr_bins, normed=True)
        cdf = imhist.cumsum()  # cumulative distribution function
        cdf = 255 * cdf / cdf[-1]  # normalize

        # use linear interpolation of cdf to find new pixel values
        im2 = np.interp(im.flatten(), bins[:-1], cdf)

        return im2.reshape(im.shape), cdf

    
class MyRectangleSelector(RectangleSelector):

    def release(self, event):
        super(MyRectangleSelector, self).release(event)
        self.to_draw.set_visible(True)
        self.canvas.draw()


class TransPlot:

    def __init__(self, directory, val):

        self.directory = directory
        try:
            self.val = val.get()
            if self.val == "Default flight path: 56m":
                self.L = float(self.val.split(':')[1].strip('m'))

            else:
                self.L = float(self.val.strip('m'))

        except NameError:
            return ctypes.windll.user32.MessageBoxA(0, "You must select an ROI first", "Error", 1)

    def produceTransData(self):

        scaledIntensities = []
        for scaled in self.directory.openFits.arrays:
            scaledIntensities.append(sum(sum(scaled[a:b,c:d])))

        sampleIntensities = []
        for sample in self.directory.sampleFits.arrays:
            sampleIntensities.append(sum(sum(sample[a:b,c:d])))

        transmitted = []
        zipped = zip(sampleIntensities, scaledIntensities)
        for sample, scaled in zipped:
            transmitted.append(float(sample)/float(scaled))

        TOF = []
        for header in self.directory.sampleFits.headers:
            TOF.append(header["TOF"])

        return TOF, transmitted

    def convertToWavelength(self, data):

        wavelength = []
        h = 6.6E-34
        m = 1.67E-27
        A = 10**10
        for point in data:
            wavelength.append(((h*float(point))/(self.L*m))*A)
        return wavelength

    def plotTransTOF(self):

        xyData = self.produceTransData()
        global timeOF
        timeOF = xyData[0]

        ymin = min(xyData[1]) - 0.05
        ymax = max(xyData[1]) + 0.05

        self.fig = plt.figure(1)
        self.ax = self.fig.add_subplot(111)
        self.ax.autoscale(enable=True, axis="both", tight=True)
        self.myrectsel = MyRectangleSelector(
            self.ax, self.onSelect, drawtype='box',rectprops=dict(
                facecolor='red', edgecolor='black', alpha=0.5, fill=True))
        plt.plot(xyData[0], xyData[1])
        plt.ylim(ymin, ymax)
        plt.xlabel("Time of Flight (s)")
        plt.ylabel("Neutron Transmission")
        plt.show()
        plt.close()

        return timeOF

    def plotTransWavelength(self):

        xyData = self.produceTransData()
        global transW
        transW = xyData[1]
        wavelength = self.convertToWavelength(xyData[0])

        ymin = min(xyData[1]) - 0.05
        ymax = max(xyData[1]) + 0.05

        self.fig = plt.figure(1)
        self.ax = self.fig.add_subplot(111)
        self.ax.autoscale(enable=True, axis="both", tight=True)
        self.myrectsel = MyRectangleSelector(
            self.ax, self.onSelect, drawtype='box', rectprops=dict(
                facecolor='red', edgecolor='black', alpha=0.5, fill=True))

        plt.plot(wavelength, xyData[1])
        plt.xlabel(u"Wavelength (\u00C5)")
        plt.ylabel("Neutron Transmission")
        plt.ylim(ymin, ymax)
        plt.show()
        plt.close()
        data = None
        return wavelength, transW

    def onSelect(self, eclick, erelease):
        print "Start position: (%f, %f)" % (eclick.xdata, eclick.ydata)
        print "End position: (%f, %f)" % (erelease.xdata, erelease.ydata)
        global atp
        atp = eclick.xdata
        global btp
        btp = erelease.xdata
        global ctp
        ctp = eclick.ydata
        global dtp
        dtp = erelease.ydata
        return atp, btp, ctp, dtp    
        
class Test:

    def __init__(self, dir):
        self.a = dir

    def do(self):
        print self.a.openFits.arrays[100]

    

if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("500x600")
    root.title("Bragg Edge Analysis")
    app = BraggEdgeAnalysisGUI(root)
    root.mainloop()
    root.destroy()
