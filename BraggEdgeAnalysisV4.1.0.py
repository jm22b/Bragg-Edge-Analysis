import Tkinter as tk
import glob
import os
import numpy as np
import ctypes
import scipy.special
import scipy.signal
import warnings
import matplotlib
import csv
# this backend must be used
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib import path
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg  # Note: add toolbar
from matplotlib.figure import Figure
from matplotlib.widgets import Slider, RectangleSelector, LassoSelector
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.signal import convolve2d
from tkFileDialog import askdirectory, asksaveasfilename
from astropy.io import fits
#from skimage.filters.rank import mean

# Main page of the GUI
class BraggEdgeAnalysisGUI:
    def __init__(self, root_):
        self.root = root_
        self.frame = tk.Frame(self.root)
        self.frame.pack() #pack is used to manage the position of widgets

        self.directory = GetDirectories()  # instantiate the class here to be passed to other classes, avoids multiples
        self.correction = OverlapCorrectionAndScaling(self.directory) # instance of overlapcorrection takes the aforemention instance of self.directory so that it has access

        self.menubar = tk.Menu(self.root)  # creates top level menus to be populated in widgets()
        self.filemenu = tk.Menu(self.menubar, tearoff=0) # add a file menu
        self.actionmenu = tk.Menu(self.menubar, tearoff=0) # add an action menu
        self.transplot = tk.Menu(self.menubar, tearoff=0) # add a sub menu for plotting functions
        self.bits = tk.Menu(self.menubar, tearoff=0) # sub menu for choosing 16/32 bit options
        self.results = tk.Menu(self.menubar, tearoff=0)
        # button for showing the sample images
        self.showDataButton = tk.Button(
            self.frame, text="Show Sample", width=10, command=lambda: ShowData(self.root, self.directory).plot())
        # text field that is used for specifying the flight path of the instrument
        self.flightpath = tk.Entry(self.frame, width=30)
        # calls the widgets function, populating the GUI with it's objects
        self.widgets()

    def widgets(self):
        
        root.option_add("*tearoff", "FALSE")
        # populates top level menus and maps the commands to them
        
        self.filemenu.add_command(
            label="Load Open Beam", command=self.directory.getOpenPath)
        
        self.filemenu.add_separator()
        self.filemenu.add_command(
            label="Load Sample Beam", command=self.directory.getSamplePath)
        
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=root.destroy)
        self.menubar.add_cascade(label="File", menu=self.filemenu)

        self.actionmenu.add_cascade(label="Correct & Scale Data", menu=self.bits)
        self.bits.add_command(
            label="16 Bit Integer Data", command=lambda: self.correction.doBoth(np.int16))
        
        self.bits.add_separator()
        self.bits.add_command(
            label="32 Bit Float Data", command=lambda: self.correction.doBoth(np.float32))
        
        self.actionmenu.add_separator()
        self.actionmenu.add_cascade(label="Plotting", menu=self.transplot)
        self.transplot.add_command(
            label="Transmission Plots", command=lambda: TransPlot(
                self.directory, self.flightpath).combinedTransPlot())
        
        self.transplot.add_separator()
        self.transplot.add_command(label="Z-Axis Profile", command=lambda: TransPlot(self.directory, self.flightpath).ZAxisProfile())

        self.actionmenu.add_separator()
        self.actionmenu.add_command(label="Fit Bragg Edge", command=lambda: EdgeFitting().subplotCall())
        
        self.actionmenu.add_separator()
        self.actionmenu.add_command(label="2D Strain Mapping", command=lambda: StrainMapping(self.directory).do())
        
        self.actionmenu.add_separator()
        self.actionmenu.add_command(label="Principal Component Analysis", command=lambda: PrincipalComponentAnalysis(self.directory))

        self.menubar.add_cascade(label="Actions", menu=self.actionmenu)

        self.results.add_command(label="Results", command=lambda: ResultsTable().populateTable())
        self.menubar.add_cascade(label="Results", menu=self.results)

        root.config(menu=self.menubar)
        # supplies a default value to the flight path
        self.flightpath.insert(0, "Default flight path: 56m")
        self.flightpath.pack()
        self.showDataButton.pack()
        # self.contrastButton.pack()


class FitsData:
    """This class acts as a data model for the open beam and sample data.
    creating an instance of it will produce a blank template to be filled with the
    relevant data by the loading functions
    """

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
    # holds the methods for opening file dialogs and finding path variables

    def __init__(self):
        self.openPath = None
        self.samplePath = None

    def openOpenDirectory(self):
        self.openPath = askdirectory() # tk file dialog
        return self.openPath

    def openSampleDirectory(self):
        self.samplePath = askdirectory()
        return self.samplePath


class GetDirectories:
    def __init__(self):
        self.directory = DirectoryHandler()
        self.openFits = FitsData() # instances of the FitsData() class are blank templates to be filled
        self.sampleFits = FitsData()

        self.openPath = None
        self.samplePath = None

    def getOpenPath(self):
        """if scaled data exists, load it. Otherwise load original data"""

        self.directory.openOpenDirectory()
        self.openPath = self.directory.openPath
        self.loadData(self.directory.openPath, self.openFits) # load data handles the logic of what to load from where

    def getSamplePath(self):
        """if overlap corrected data exists, load it. otherwise load the original data"""

        self.directory.openSampleDirectory()
        self.samplePath = self.directory.samplePath
        self.loadData(self.directory.samplePath, self.sampleFits) # the same load data function is called

    def loadData(self, path, container):
        """handles the loading of the data, filling up the 'container' (FitsData() instance)
        behaves identically for both open beam and sample data"""

        f = glob.glob(os.path.join(path, "*[0-9][0-9][0-9][0-9][0-9].fits"))  # relies on 5-digit numbering of files
        for fitsFile in f:
            hdulist = fits.open(fitsFile, memmap=False)  # memmap tries to keep things open after closing them
            name = hdulist.filename().split("\\")[-1] # does soom foo to extract filename
            header = hdulist[0].header # accesses header object
            data = hdulist[0].data # accesses data object
            hdulist.close() # close the file
            container.names.append(name)  # populate container with name, header, data
            container.headers.append(header)
            container.arrays.append(data)
            print name  # debugging


class OverlapCorrectionAndScaling:
    """deals with overlap correction and scaling of the selected data. Needs refactoring and the reasoning about
    what these functions will do needs pushing up towards the GUI."""

    def __init__(self, directory):

        self.directory = directory

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

        """
        This function is responsible for determining the indices of the files that
        belong to each shutter
        """

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

    def overlapCorrection(self, path, bits):

        """most hopeful suspect for refactoring and pushing of logic upwards."""

        shutterIndices = self.preBinData(path)
        shutterValues = self.readShutter(path)[0]

        if bits == np.int16:
            os.mkdir(os.path.join(path, "16-bit-overlapCorrected"))
            self.f = open(os.path.join(path, "16-bit-overlapCorrected", "TOFData.csv"), "wb")
            # zipped = zip(self.sampleFits.arrays, self.sampleFits.headers, self.sampleFits.names)
            s = 0
            for subIndex in shutterIndices:
                # subList = zipped[subIndex[0]:subIndex[-1]+1]
                runningTot = np.zeros((512, 512))

                for i in range(subIndex[0], subIndex[0] + len(subIndex)):
                    shutter = float(shutterValues[s][1])
                    prob = np.divide(runningTot, shutter)
                    runningTot += self.directory.sampleFits.arrays[i]
                    self.directory.sampleFits.arrays[i] = np.round(
                        np.divide(self.directory.sampleFits.arrays[i], (1 - prob))).astype(np.int16)
                    self.writeToFolder(
                        self.directory.sampleFits.arrays[i], self.directory.sampleFits.headers[i],
                        self.directory.sampleFits.names[i], path, "16-bit-overlapCorrected", "16-bit-corrected")
                    print i
                print s
                s += 1
        else:
            os.mkdir(os.path.join(path, "32-bit-overlapCorrected"))
            self.f = open(os.path.join(path, "32-bit-overlapCorrected", "TOFData.csv"), "wb")
            s = 0
            for subIndex in shutterIndices:
                runningTot = np.zeros((512, 512))
                for i in range(subIndex[0], subIndex[0] + len(subIndex)):
                    shutter = float(shutterValues[s][1])
                    prob = np.divide(runningTot, shutter)
                    runningTot += self.directory.sampleFits.arrays[i]
                    self.directory.sampleFits.arrays[i] = (
                        np.divide(self.directory.sampleFits.arrays[i], (1 - prob))).astype(bits)
                    self.writeToFolder(
                        self.directory.sampleFits.arrays[i], self.directory.sampleFits.headers[i],
                        self.directory.sampleFits.names[i], path, "32-bit-overlapCorrected", "32-bit-corrected")
                    print i
                print s
                s += 1

        self.f.close()

    def overlapCorrectionScaling(self, path, bits):

        shutterIndices = self.preBinData(path)
        shutterValuesOpen = self.readShutter(path)[0]
        shutterValuesSample = self.readShutter(self.directory.samplePath)[0]

        zipShutters = zip(shutterValuesOpen, shutterValuesSample)
        ratio = []
        for svo, svs in zipShutters:
            ratio.append(float(svs[1]) / float(svo[1]))

        if bits == np.int16:
            os.mkdir(os.path.join(path, "16-bit-scaledOpenBeam"))
            self.f = open(os.path.join(path, "16-bit-scaledOpenBeam", "TOFData.csv"), "wb")
            # fmod = str(bits).split('.')[-1][0:-2]
            s = 0
            for subIndex in shutterIndices:
                # sublist = zipped[subIndex[0]:subIndex[-1]+1]
                runningTot = np.zeros((512, 512))
                scaleFactor = ratio[s]

                for i in range(subIndex[0], subIndex[0] + len(subIndex)):
                    shutter = float(shutterValuesOpen[s][1])
                    prob = np.divide(runningTot, shutter)
                    runningTot += self.directory.openFits.arrays[i]
                    self.directory.openFits.arrays[i] = np.round(
                        (np.divide(self.directory.openFits.arrays[i], (1 - prob))) * scaleFactor).astype(bits)
                    self.writeToFolder(
                        self.directory.openFits.arrays[i], self.directory.openFits.headers[i],
                        self.directory.openFits.names[i], path, "16-bit-scaledOpenBeam", "16-bit-scaled")
                    print i
                print s
                s += 1

        else:
            os.mkdir(os.path.join(path, "32-bit-scaledOpenBeam"))
            self.f = open(os.path.join(path, "32-bit-scaledOpenBeam", "TOFData.csv"), "wb")
            # fmod = str(bits).split('.')[-1][0:-2]
            # zipped = zip(self.openFits.arrays, self.openFits.headers, self.openFits.names)
            s = 0

            for subIndex in shutterIndices:
                # sublist = zipped[subIndex[0]:subIndex[-1]+1]
                runningTot = np.zeros((512, 512))
                scaleFactor = ratio[s]

                for i in range(subIndex[0], subIndex[0] + len(subIndex)):
                    shutter = float(shutterValuesOpen[s][1])
                    prob = np.divide(runningTot, shutter)
                    runningTot += self.directory.openFits.arrays[i]
                    self.directory.openFits.arrays[i] = (
                        (np.divide(self.directory.openFits.arrays[i], (1 - prob))) * scaleFactor).astype(bits)
                    self.writeToFolder(
                        self.directory.openFits.arrays[i], self.directory.openFits.headers[i],
                        self.directory.openFits.names[i], path, "32-bit-scaledOpenBeam", "32-bit-scaled")
                    print i
                print s
                s += 1
        self.f.close()

    def writeToFolder(self, array, header, name, path, mod1, mod2):

        """
        generic function for saving the data once it has been corrected/scaled
        """

        hdu = fits.PrimaryHDU() # create fits file
        hdu.data = array # populate file with data
        counts = sum(sum(array)) # get new number of counts
        hdu.header = header # copy over old headers
        hdu.header["N_COUNTS"] = counts # update with new number of counts
        TOF = hdu.header["TOF"] # extract TOF
        line = "%.16f, %d\n" % (TOF, counts) # line to be written to .csv file
        self.f.writelines(line) # write line to.csv
        hdu.writeto(os.path.join(path, mod1, mod2 + name)) # saves the file based on the args of the function

    def doBoth(self, bits):

        """
        This function calls the overlap correction functions and catches the exception
        caused if the data already exists
        """
        
        try:
            fmod = str(bits).split('.')[-1][0:-2][-2:] + "-bit-" # foo to extract a string representing the number of bits being used (used to identify the saved data)
            if not os.path.exists(os.path.join(self.directory.samplePath, fmod + "overlapCorrected")):
                self.overlapCorrection(self.directory.samplePath, bits)
            else:
                ctypes.windll.user32.MessageBoxA(0, "Corrected files already exist.", "Error", 1)
            if not os.path.exists(os.path.join(self.directory.openPath, fmod + "scaledOpenBeam")):
                self.overlapCorrectionScaling(self.directory.openPath, bits)
            else:
                ctypes.windll.user32.MessageBoxA(0, "Scaled and corrected files already exist.", "Error", 1)
        except TypeError:
            return ctypes.windll.user32.MessageBoxA(0, "You need to select some data", "Error", 1)


class ShowData:

    """
    This class handles the data visualisation of the sample
    """
    
    def __init__(self, root, directory):
        self.root = root # want to show in the 'root' page
        self.directory = directory

        self.fig = Figure(figsize=(7, 7)) 
        self.ax = self.fig.add_subplot(111)
        self.plotted = False
        self.l = None
        self.canvas = None
        plt.show() # shows the figure upon initialisation

        self.slider = tk.Scale(
            self.root, from_=0, to=(len(self.directory.sampleFits.arrays) - 1), resolution=1, orient=tk.HORIZONTAL,
            command=self.update) # adds a slider to move through the image stack
        self.slider.pack()

        self.vmax = tk.Entry(self.root, width=10) # Text entry takes the vmax argument
        self.vmax.insert(0, "100") # 100 is a reasonable default value for most imgs
        self.vmax.pack()
        # button for accessing the histogram equalisation function
        self.button = tk.Button(text="Histogram Equalisation", command=self.contrast)
        self.button.pack()

    def onSelect(self, eclick, erelease):

        """
        This function handles the selection of an ROI, returning the values of the
        rectangles corners
        """
        print "Start position: (%f, %f)" % (eclick.xdata, eclick.ydata)
        print "End position: (%f, %f)" % (erelease.xdata, erelease.ydata)
        global a
        a = int(eclick.xdata)
        global b
        b = int(erelease.xdata)
        global c
        c = int(eclick.ydata)
        global d
        d = int(erelease.ydata)
        return a, b, c, d

    def plot(self):

        """
        shows first image of stack
        """
        self.plotted = True
        self.s = 0
        self.canvas = FigureCanvasTkAgg(self.fig, self.root)
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()
        self.myrectsel = MyRectangleSelector(self.ax, self.onSelect, drawtype="box", rectprops=dict(
            facecolor="red", edgecolor="black", alpha=0.2, fill=True)) #sub class to force the rectangle to persist

    def update(self, val):

        """
        This function deals with updating the canvas when the slider moves
        """
        global sliderInd
        sliderInd = int(self.slider.get()) # position of slider for indexing stack
        if self.plotted:
            im = self.directory.sampleFits.arrays[sliderInd]
            self.l = self.ax.imshow(im, cmap=plt.cm.gray, interpolation="nearest", vmin=0, vmax=int(self.vmax.get()))
            # self.l.set_data(im)
            # print self.directory.sampleFits.arrays[ind]
            self.canvas.draw()
        return sliderInd

    def histeq(self, im, nbr_bins=256):
        # get image histogram
        imhist, bins = np.histogram(im.flatten(), nbr_bins, normed=True)
        cdf = imhist.cumsum()  # cumulative distribution function
        cdf = 255 * cdf / cdf[-1]  # normalize

        # use linear interpolation of cdf to find new pixel values
        im2 = np.interp(im.flatten(), bins[:-1], cdf)

        return im2.reshape(im.shape), cdf

    def contrast(self):

        """
        applies the histogram equalisation to the current slice
        """
        ind = int(self.slider.get())
        im = self.histeq(self.directory.sampleFits.arrays[ind])[0]
        self.l.set_data(im)
        self.canvas.draw()


class MyRectangleSelector(RectangleSelector):
    """
    This class overrides the default behavior of RectangleSelector in order to make it
    remain visible after releasing the mouse
    """
    def release(self, event):
        super(MyRectangleSelector, self).release(event)
        self.to_draw.set_visible(True)
        self.canvas.draw()
        

class TransPlot:
    def __init__(self, directory, val):

        self.directory = directory
        self.val = val.get()
        if self.val == "Default flight path: 56m":
            self.L = float(self.val.split(':')[1].strip('m'))

        else:
            self.L = float(self.val.strip('m'))

        self.curr_pos = 0
        self.currT_pos = 0

    def produceTransData(self):

        scaledIntensities = []
        for scaled in self.directory.openFits.arrays:
            scaledIntensities.append(sum(sum(scaled[c:d, a:b])))

        sampleIntensities = []
        for sample in self.directory.sampleFits.arrays:
            sampleIntensities.append(sum(sum(sample[c:d, a:b])))

        transmitted = []
        zipped = zip(sampleIntensities, scaledIntensities)
        for sample, scaled in zipped:
            transmitted.append(float(sample) / float(scaled))

        TOF = []
        for header in self.directory.sampleFits.headers:
            TOF.append(header["TOF"])

        return TOF, transmitted

    def convertToWavelength(self, data):
        convertedwavelength = []
        h = 6.6E-34
        m = 1.67E-27
        A = 10 ** 10
        for point in data:
            convertedwavelength.append(((h * float(point)) / (self.L * m)) * A)
        return convertedwavelength
    
    def combinedTransPlot(self):
        XYData = self.produceTransData()
        global Transmitted
        Transmitted = XYData[1]
        global TimeOfFlight
        TimeOfFlight = XYData[0]
        global wavelength
        wavelength = self.convertToWavelength(TimeOfFlight)
        
        self.TransPlots = [(TimeOfFlight, Transmitted),(wavelength, Transmitted)]
        self.fig2 = plt.figure(2)
        
        self.ax2 = self.fig2.add_subplot(111)
        self.fig2.canvas.mpl_connect("Rectangle Select", MyRectangleSelector(
            self.ax2, self.onSelect, drawtype='box', rectprops=dict(
                facecolor='red', edgecolor='black', alpha=0.5, fill=True)))
        self.fig2.canvas.mpl_connect("key_press_event", self.key_event_Transmission)
        self.ax2.ymin = np.min(Transmitted) - 0.05
        self.ax2.ymax = np.max(Transmitted) + 0.05
        self.ax2.plot(TimeOfFlight, Transmitted, 'x', ms=3)
        self.xTlabels = ["TOF (s)", u"Wavelength (\u00C5)"]
        self.ax2.set_xlabel(self.xTlabels[0])
        self.ax2.set_ylabel("Neutron Transmission")
        plt.show()
        return Transmitted, TimeOfFlight, wavelength
        
    def key_event_Transmission(self, e):
        
        if e.key == "right":
            self.currT_pos += 1
        elif e.key == "left":
            self.currT_pos -= 1
        else:
            return
        self.currT_pos %= len(self.TransPlots)
        self.ax2.cla()
        self.ax2.plot(self.TransPlots[self.currT_pos][0], self.TransPlots[self.currT_pos][1], 'x', ms=3)
        self.ax2.set_xlabel(self.xTlabels[self.currT_pos])
        self.ax2.set_ylabel("Neutron Transmission")
        self.myrectsel = MyRectangleSelector(
            self.ax2, self.onSelect, drawtype='box', rectprops=dict(
                facecolor='red', edgecolor='black', alpha=0.5, fill=True))
        self.fig2.canvas.draw()

    def ZAxisProfile(self):

        TOF = []
        for header in self.directory.sampleFits.headers:
            TOF.append(header['TOF'])
        print len(TOF)

        wavelengthZ = self.convertToWavelength(TOF)
        print len(wavelengthZ)

        avg = []
        for data in self.directory.sampleFits.arrays:
            avg.append(np.mean(data[c:d, a:b]))

        self.plots = [(TOF, avg), (wavelengthZ, avg)]
        self.fig1 = plt.figure(1)
        self.fig1.canvas.mpl_connect("key_press_event", self.key_event_ZAxis)
        self.ax = self.fig1.add_subplot(111)
        self.ax.plot(TOF, avg)

        self.xlabels = ["TOF (s)", u"Wavelength (\u00C5)"]
        self.ax.set_xlabel(self.xlabels[0])
        self.ax.set_ylabel("Average Number of Counts")
        plt.show()

    def key_event_ZAxis(self, e):

        if e.key == "right":
            self.curr_pos += 1
        elif e.key == "left":
            self.curr_pos -= 1
        else:
            return
        self.curr_pos %= len(self.plots)
        self.ax.cla()
        self.ax.plot(self.plots[self.curr_pos][0], self.plots[self.curr_pos][1])
        self.ax.set_xlabel(self.xlabels[self.curr_pos])
        self.ax.set_ylabel("Average Number of Counts")
        self.fig1.canvas.draw()

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


class ResultsTable:

    def __init__(self):

        self.frame = tk.Toplevel()
        self.textFrame = tk.Frame(self.frame, width = 400, height = 600)

        self.textwidget = tk.Text(self.textFrame, borderwidth=3, relief="sunken")
        
        self.scrollbar = tk.Scrollbar(self.textFrame, command = self.textwidget.yview)
        
        self.menubar = tk.Menu(self.frame)
        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        
        self.TOFlabel = tk.Label(self.frame, text=u"TOF (s)      wavelength (\u00C5) Transmisson")
        #self.wavelengthlabel = tk.Label(self.frame, text=u"wavelength (\u00C5)")
        #self.transmissionlabel = tk.Label(self.frame, text="Transmisson")
        
        self.widgets()

    def widgets(self):
        
        self.TOFlabel.pack()
        #self.wavelengthlabel.pack(side="left")
        #self.transmissionlabel.pack(side="left")
        
        self.textFrame.pack(side='bottom', fill="both", expand=True)
        self.textFrame.grid_propagate(False)
        self.textFrame.grid_rowconfigure(0, weight=1)
        self.textFrame.grid_columnconfigure(0, weight=1)
        
        self.textwidget.config(font=("consolas", 12), undo=True, wrap='word')
        self.textwidget.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        
        self.scrollbar.grid(row=0, column=1, sticky="nsew")
        self.textwidget['yscrollcommand'] = self.scrollbar.set
        
        self.filemenu.add_command(label="Save As", command=self.save)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.frame.config(menu=self.menubar)
        
    def save(self):
        
        name = asksaveasfilename()
        print name
        if name == None:
            return
        contents = self.textwidget.get("1.0", "end-1c")
        f = open(name, "wb")
        f.writelines(contents.replace(" ", ","))
        f.close()
        
    def populateTable(self):

        zipped = zip(TimeOfFlight, wavelength, Transmitted)
        #results = []
        for x,y,z in zipped:
            xyzstr = "%f \t%f \t%f\n" % (x,y,z)
            self.textwidget.insert(tk.END, xyzstr)


class EdgeFitting:
    def __init__(self):

        # self.xvalsW = wavelength
        # self.trans = transW
        # self.xvalsT = timeOF
        self.subx = []
        self.suby = []

        self.frame = tk.Toplevel()
        self.fig = Figure(figsize=(5, 5))
        self.ax = self.fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0)

        # self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame)
        # self.toolbar.update()

        self.plotButton = tk.Button(self.frame, text="Fit Curve", command=self.fitCurve)

        self.lambda0var = tk.Entry(self.frame, width=10)
        self.lambda0label = tk.Label(self.frame, text=u"      Lambda Edge \u03BB \u2080")
        self.sigmavar = tk.Entry(self.frame, width=10)
        self.sigmalabel = tk.Label(self.frame, text=u" Bragg Edge Width \u03C3")
        self.tauvar = tk.Entry(self.frame, width=10)
        self.taulabel = tk.Label(self.frame, text=u"   Edge Asymmetry \u03C4")

        self.widgets()

        self.lambda0 = self.lambda0var.get()
        self.sigma = self.sigmavar.get()
        self.tau = self.tauvar.get()

    def widgets(self):

        self.plotButton.grid(row=1)


        self.lambda0label.grid(sticky="W")
        self.lambda0var.grid(row=2)
        #self.lambda0var.insert(0, "1")
        self.sigmalabel.grid(sticky="W")
        self.sigmavar.grid(row=3)
        self.sigmavar.insert(0, "1")
        self.taulabel.grid(sticky="W")
        self.tauvar.grid(row=4)
        self.tauvar.insert(0, "0.01")

    def noiseFiltering(self, y):
        """
        low frequency filter to remove noise for parameter estimation
        """
        b, a = scipy.signal.butter(3, 0.05)
        zi = scipy.signal.lfilter_zi(b, a)
        z, _ = scipy.signal.lfilter(b, a, y, zi = zi*y[0])
        z2, _ = scipy.signal.lfilter(b, a, y, zi = zi*z[0])
        ynew = scipy.signal.filtfilt(b,a,y)
        
        return ynew
    
    def parameterEstimation(self, x, y):
        
        self.ynew = self.noiseFiltering(y)
        """
        Use first derivative to find edge position (maximal gradient of filtered data indicates position)
        """
        dy = np.diff(self.ynew)
        dx = np.diff(x)
        dydx = dy/dx
        global edgeIndex
        edgeIndex = np.argmax(dydx)
        lambda_0 = x[edgeIndex]
        return lambda_0, edgeIndex

    def subPlot(self, XData, YData):
        
        zipped = zip(XData, YData)
        #pos = 0
        global posList 
        posList = []
        
        for xval, yval in zipped:
            if xval >= atp and xval <=btp:

                self.subx.append(xval)
                self.suby.append(yval)
                posList.append(zipped.index((xval,yval)))

        self.ax.plot(self.subx, self.suby, 'x', ms=3)
        self.arrx = np.array(self.subx)
        self.arry = np.array(self.suby)
        lambda_0, self.edgeIndex = self.parameterEstimation(self.arrx, self.arry)
        self.lambda0var.insert(0, lambda_0)
        return posList
    
    def subplotCall(self):
        
        if atp < 1:
            self.subPlot(TimeOfFlight, Transmitted)
        else:
            self.subPlot(wavelength, Transmitted)
            
    def rightFunc(self, x, m, c):
        """
        function describing bragg edge curve for x >> lambda_0
        """
        return np.exp(-(m*x + c))
    
    def leftFunc(self, x, b, a):
        """
        function describing bragg edge curve for x << lambda_0
        """
        return np.exp(-(self.m_right*x + self.a_right))*np.exp(-(a + b*x))
    
    def centralFunc(self, x, lambda_0, sigma, tau):
        """
        function describing bragg edge curve in the vicinity of lambda_0
        """
        x_sig = -(x - lambda_0) / (np.sqrt(2)*sigma)
        x_tau = -(x - lambda_0) / tau
        sig_tau = sigma / tau
        b = 0.5*(
            scipy.special.erfc(x_sig) - np.exp(x_tau + (0.5*sig_tau**2))*scipy.special.erfc(x_sig + sig_tau))
                  
        return np.exp(
            -(self.a_right + self.m_right*x))*(
            np.exp(-(self.a_left + self.m_left*x)) + (1-np.exp(-(self.a_left + self.m_left*x)))*b)

    def fitCurve(self):
        global fit_params #m_r, a_r, m_l, a_l
        fit_params = []
        warnings.simplefilter("error", OptimizeWarning)

        try:

            self.ax.cla()
            self.ax.plot(self.subx, self.suby, 'x')
            global initial_guess
            initial_guess = [float(self.lambda0var.get()), float(self.sigmavar.get()), float(self.tauvar.get())]

            self.m_right = (self.ynew[self.edgeIndex+20:][-1] - self.ynew[self.edgeIndex+20:][0])/(self.subx[self.edgeIndex+20:][-1] - self.subx[self.edgeIndex+20:][0])
            self.a_right = np.median(self.suby[self.edgeIndex:])
            popt_r, pcov = curve_fit(self.rightFunc, self.subx[self.edgeIndex+30:], self.suby[self.edgeIndex+30:], p0=[self.m_right, self.a_right])
            self.m_right = popt_r[0]
            self.a_right = popt_r[1]
            fit_params.append(self.m_right)
            fit_params.append(self.a_right)
            
            self.m_left = (self.ynew[0:self.edgeIndex-20][-1] - self.ynew[0:self.edgeIndex-20][0])/(self.subx[0:self.edgeIndex-20][-1] - self.subx[0:self.edgeIndex-20][0])
            self.a_left = np.median(self.suby[0:self.edgeIndex])
            popt_l, pcov = curve_fit(self.leftFunc, self.subx[:self.edgeIndex-20], self.suby[:self.edgeIndex-20], p0=[self.m_left, self.a_left])
            self.m_left = popt_l[0]
            self.a_left = popt_l[1]
            fit_params.append(self.m_left)
            fit_params.append(self.a_left)
            
            popt_c, pcov = curve_fit(self.centralFunc, self.arrx, self.arry, p0=initial_guess)
            
            print " a_hkl: %f \n b_hkl: %f \n a_0: %f \n b_0: %f \n lambda_0: %f \n sigma: %f \n tau: %f \n" % (
                self.a_right, self.m_right, self.a_left, self.m_left, popt_c[0], popt_c[1], popt_c[2])
            
            self.ax.plot(self.arrx, self.centralFunc(self.arrx, popt_c[0], popt_c[1], popt_c[2]))
            self.canvas.show()
            self.clearText()
            self.lambda0var.insert(0, popt_c[0])
            self.sigmavar.insert(0, popt_c[1])
            self.tauvar.insert(0, popt_c[2])
            initial_guess = [float(self.lambda0var.get()), float(self.sigmavar.get()), float(self.tauvar.get())]
            return initial_guess, fit_params
         
        except (RuntimeError, OptimizeWarning):
            self.ax.cla()
            self.ax.plot(self.subx, self.suby, 'x')
            x = np.linspace(self.subx[0], self.subx[-1], 100)
            self.ax.plot(
                x, self.centralFunc(
                    x, initial_guess[0], initial_guess[1], initial_guess[2]))
            self.canvas.show()
            return ctypes.windll.user32.MessageBoxA(0, "Please refine your parameters", "Error", 1)
        
    

    def clearText(self):
        fields = [self.lambda0var, self.sigmavar, self.tauvar]
        for field in fields:
            field.delete(0, "end")
        
        
class StrainMapping:
        
    def __init__(self, directory):
        
        self.directory = directory
        self.sampleArray = self.directory.sampleFits.arrays
        self.openArray = self.directory.openFits.arrays
        self.im = self.sampleArray[sliderInd]
        
        self.frame = tk.Toplevel()
        self.fig = Figure(figsize=(7, 7))
        self.ax = self.fig.add_subplot(111)
        
        self.strainButton = tk.Button(
            self.frame, text="Strain Map", command=self.strainMap)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0)
        self.canvas.mpl_connect('key_press_event', self.onKey)
        self.strainButton.grid(row=1)
        
        self.mask = np.zeros((512,512))
        self.pix = np.arange(512)
        self.XX, self.YY = np.meshgrid(self.pix, self.pix)
        self.pix = np.vstack((self.XX.flatten(), self.YY.flatten())).T
        self.lasso = LassoSelector(self.ax, self.onselect)
        
        
        
    def do(self):
        self.canvas.mpl_connect('key_press_event', self.onKey)
        self.ax.imshow(self.im, cmap = plt.cm.gray)
        

    def strainMap(self):

        beg = posList[0]
        end = posList[-1]
        zipped = zip(self.sampleArray[beg:end], self.openArray[beg:end])
        transmitted = np.zeros((len(zipped),1,512*512)).astype(np.float32)
        
        l = 0
        kernel = np.ones((50,50))
        kernel = kernel / kernel.sum()
        
        
        for sample, empty in zipped:
            sample = sample * self.mask.reshape(512,512)
            empty = empty * self.mask.reshape(512,512)
            
            transmitted[l] = np.nan_to_num(convolve2d(sample / empty, kernel, mode='same')).flatten() #convolve2d(sample, kernel, mode='same').flatten() / convolve2d(empty, kernel, mode='same').flatten()
            l += 1
            print l
        
        
        lambdas = []
        for c in range(512*512):
            if transmitted[:,:,c].all() == False:
                lambdas.append(0)
                print 'empty'

            else:
                try:
                    fit_params = []
                    self.m_right =(np.dstack(transmitted[:,:,c][edgeIndex+20:][-1])[0][0] - np.dstack(transmitted[:,:,c][edgeIndex+20:][0])[0][0]) / (np.array(wavelength)[beg:end][edgeIndex+20:][-1] - np.array(wavelength)[beg:end][edgeIndex+20:][0])
                    self.a_right = np.median(np.dstack(transmitted[:,:,c][edgeIndex:])[0][0])
                    popt_r, pcov = curve_fit(
                        self.rightFunc, np.array(wavelength)[beg:end][edgeIndex+30:], np.dstack(transmitted[:,:,c][edgeIndex+30:])[0][0], p0=[self.m_right, self.a_right])
                    
                    self.m_right = popt_r[0]
                    self.a_right = popt_r[1]
                    fit_params.append(self.m_right)
                    fit_params.append(self.a_right)
                    
                    self.m_left = (np.dstack(transmitted[:,:,c][0:edgeIndex-20][-1])[0][0] - np.dstack(transmitted[:,:,c][0:edgeIndex-20][0])[0][0]) / (np.array(wavelength)[beg:end][0:edgeIndex-20][-1] - np.array(wavelength)[beg:end][0:edgeIndex-20][0])
                    self.a_left = np.median(transmitted[:,:,c][0:edgeIndex])
                    popt_l, pcov = curve_fit(self.leftFunc, np.array(wavelength)[beg:end][:edgeIndex-20], transmitted[:,:,c][:edgeIndex-20], p0=[self.m_left, self.a_left])
            
                    self.m_left = popt_l[0]
                    self.a_left = popt_l[1]
                    fit_params.append(self.m_left)
                    fit_params.append(self.a_left)
            
                    popt, pcov = curve_fit(self.centralFunc, np.array(wavelength)[beg:end], np.dstack(transmitted[:,:,c])[0][0], p0=initial_guess)
                
                    lambdas.append((popt[2] - initial_guess[2])/initial_guess[2])
                    print 'full'
                except (OptimizeWarning, RuntimeError):
                    #lambdas.append((initial_guess[2] - popt[2])/initial_guess[2])
                    a = 1
                    lambdas.append(a)
                    print 'Exception'
                
                
        strainMap = np.array(lambdas).reshape(512,512)*self.mask.reshape(512,512) 
        strainMap = np.ma.masked_where(strainMap == 0, strainMap)
        minVal = strainMap.min()
        maxVal = strainMap.max()
        cmap = plt.cm.bwr
        cmap.set_bad(color='black')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cax = ax.imshow(strainMap, interpolation='None', cmap=cmap)
        cbar = fig.colorbar(cax, ticks=[minVal, 0, maxVal])
        cbar.ax.set_yticklabels([minVal, '0', maxVal])
        plt.show()
        plt.close()
        
     
    def rightFunc(self, x, m, const):
        """
        function describing bragg edge curve for x >> lambda_0
        """
        return np.exp(-(m*x + const))
    
    def leftFunc(self, x, b, a):
        """
        function describing bragg edge curve for x << lambda_0
        """
        return np.exp(-(self.m_right*x + self.a_right))*np.exp(-(a + b*x))
        
    def centralFunc(self, x, lambda_0, sigma, tau):
        """
        function describing bragg edge curve in the vicinity of lambda_0
        """
        x_sig = -(x - lambda_0) / (np.sqrt(2)*sigma)
        x_tau = -(x - lambda_0) / tau
        sig_tau = sigma / tau
        b = 0.5*(
            scipy.special.erfc(x_sig) - np.exp(x_tau + (0.5*sig_tau**2))*scipy.special.erfc(x_sig + sig_tau))
                  
        return np.exp(
            -(fit_params[1] + fit_params[0]*x))*(
            np.exp(-(fit_params[3] + fit_params[2]*x)) + (1-np.exp(-(fit_params[3] + fit_params[2]*x)))*b)
        
    def updateArray(self, im, indices,mask):
        lin = np.arange(self.im.size)
        self.mask = self.mask.flatten()
        self.mask[lin[indices]] = 1
        newArray = im.flatten()
        #newArray[lin[indices]] = 1
        newArray = newArray*self.mask
        self.ax.imshow(newArray.reshape(self.im.shape), cmap=plt.cm.gray)
        self.canvas.draw()
        return newArray.reshape(self.im.shape)

    def onselect(self, verts):
        p = path.Path(verts)
        ind = p.contains_points(self.pix, radius=1)
        array = self.updateArray(self.im, ind, self.mask)
        self.canvas.draw_idle()
        
    def onKey(self, event):
        print "key pressed"
        if event.key == 'r':
            print "resetting mask"
            self.ax.imshow(self.im, cmap=plt.cm.gray)
            self.canvas.draw()
            
class PrincipalComponentAnalysis:
    
    def __init__(self, directory):
        self.sampleArrays = directory.sampleFits.arrays
        self.frame = tk.Toplevel()
        self.slicesEntry = tk.Entry(self.frame, width = 30)
        self.PCAbutton = tk.Button(
            self.frame, text="Perform PCA", width=10, command=self.controller)
        
        self.fig = Figure(figsize=(7,7))
        self.ax = self.fig.add_subplot(111)
        self.plotted = False
        self.l = None
        self.canvas = None
        plt.show()
        
        
        
        self.widgets()
        
    def widgets(self):
        self.slicesEntry.insert(0, "number of slices to combine")
        self.slicesEntry.pack()
        self.PCAbutton.pack()
    
    def plot(self):
        
        self.plotted = True
        self.s = 0
        self.canvas = FigureCanvasTkAgg(self.fig, self.frame)
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()
    
    def update(self, val):
        global PCAslider
        PCAslider = int(self.slider.get())
        if self.plotted:
            im = self.PCAimages[PCAslider]
            self.l = self.ax.imshow(im, cmap=plt.cm.gray, interpolation=None)
            self.l.set_data(im)
            # print self.directory.sampleFits.arrays[ind]
            self.canvas.draw()
        return PCAslider
    
    def pca(self, X):
        
        numDat, dims = X.shape
        meanX = X.mean(axis = 0)
        for i in range(numDat):
            X[i] = X[i] - meanX
        

        print "compacting"
        M = np.dot(X, X.T)
        e, EV = np.linalg.eigh(M)
        tmp = np.dot(X.T, EV).T
        V = tmp[::-1]
        S = np.sqrt(e)[::-1]
            

       
        return V, S, meanX
    
    def creator(self, a, b):
        immatrix = np.array([self.sampleArrays[i].flatten() for i in range(a, b)], 'f')
        V, S, immean = self.pca(immatrix)
        immean = immean.reshape(self.m, self.n)
        mode = V[0].reshape(self.m, self.n)
        
        return mode
        #self.fig1 = plt.figure(1)
        #self.ax1 = self.fig1.add_subplot(111)
        #self.ax1.imshow(immean, cmap=plt.cm.gray)
        
        #self.fig2 = plt.figure(2)
        #self.ax2 = self.fig2.add_subplot(111)
        #self.ax2.imshow(mode, cmap=plt.cm.gray)
        #plt.show()
        
    def controller(self):
        self.im = self.sampleArrays[0]
        self.m, self.n = self.im.shape[0:2]
        self.imnbr = len(self.sampleArrays)
        
        a, b = 0, 1
        slicesToCombine = int(self.slicesEntry.get())
        newLength = self.imnbr / slicesToCombine
        self.PCAimages = []
        slices = []
        
        self.slider = tk.Scale(
            self.frame, from_=0, to=newLength-1, resolution=1, orient=tk.HORIZONTAL, command=self.update)
        self.slider.pack()
         
        for i in range(newLength):
            slices.append(i*slicesToCombine)
            
        slices.append(self.imnbr)
        
        while b < newLength:
            self.PCAimages.append(self.creator(slices[a], slices[b]))
            a += 1
            b +=1
        
        self.plot()
        

if __name__ == "__main__":
    root = tk.Tk()
    # root.geometry("600x700")
    root.title("Bragg Edge Analysis")
    app = BraggEdgeAnalysisGUI(root)
    root.mainloop()
    root.destroy()
