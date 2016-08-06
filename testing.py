import Tkinter as tk
import glob
import os
import numpy as np
import ctypes
import matplotlib
import scipy.special
import warnings
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.widgets import Slider, RectangleSelector
from tkFileDialog import askdirectory
from astropy.io import fits
from scipy.optimize import curve_fit, OptimizeWarning

class BraggEdgeAnalysisGUI:

    def __init__(self, root_):
        self.root = root_
        self.frame = tk.Frame(self.root)
        self.frame.pack()

        self.getFiles = FileHandler()
        # self.openFiles = OpenFiles()
        self.overlap = OverlapCorrectionAndScaling()
        # self.plot = TransPlot()

        self.menubar = tk.Menu(root)
        self.filemenu = tk.Menu(self.menubar, tearoff=0)
        self.actionmenu = tk.Menu(self.menubar, tearoff=0)
        self.transplot = tk.Menu(self.menubar)

        # self.testButton = tk.Button(self.frame, text="Test", width=10, command=lambda: EdgeFitting().subPlot())

        self.showButton = tk.Button(
            self.frame, text="Show Data", width=10, command=lambda: ShowData(self.root).plot())

        self.flightpath = tk.Entry(self.frame, width=30)

        # self.test = Test(self.flightval)

        self.widgets()

    def widgets(self):
        root.option_add("*tearoff", "FALSE")
        self.filemenu.add_command(
            label="Load Open Beam", command=self.getFiles.getOpenPath)

        self.filemenu.add_separator()

        self.filemenu.add_command(
            label="Load Sample", command=self.getFiles.getSamplePath)

        self.filemenu.add_separator()

        self.filemenu.add_command(label="Exit", command=root.destroy)

        self.menubar.add_cascade(label="File", menu=self.filemenu)

        self.actionmenu.add_command(
            label="Correct and Scale data", command=lambda: self.overlap.overlapCorrection(
                openPath) & self.overlap.overlapCorrection(samplePath))
        self.actionmenu.add_separator()
        self.actionmenu.add_cascade(label="Transmission", menu=self.transplot)
        self.transplot.add_command(label="Plot (TOF)", command=lambda: TransPlot(self.flightpath).plotTrans(True))
        self.transplot.add_separator()
        self.transplot.add_command(label="Plot (wavelength)", command=lambda: TransPlot(self.flightpath).plotTrans())
        self.actionmenu.add_separator()
        self.actionmenu.add_command(label="Fit Bragg Edge", command=lambda: EdgeFitting().subPlot())
        self.menubar.add_cascade(label="Actions", menu=self.actionmenu)
        root.config(menu=self.menubar)

        self.flightpath.insert(0, "Default flight path: 56m")
        self.flightpath.pack()
        # self.testButton.pack()
        self.showButton.pack()


class DirectoryHandler:

    def __init__(self):
        self.openPath = None
        self.samplePath = None

    def openOpenDirectory(self):
        global openPath
        openPath = askdirectory()
        return openPath

    def openSampleDirectory(self):
        global samplePath
        samplePath = askdirectory()
        return samplePath

    def pathToOpen(self):
        # print self.openPath
        return self.openPath

    def pathToSample(self):
        # print self.samplePath
        return self.samplePath


class FileHandler:

    # noinspection PyDefaultArgument,PyDefaultArgument,PyDefaultArgument,PyDefaultArgument,PyDefaultArgument,PyDefaultArgument,PyDefaultArgument,PyDefaultArgument,PyDefaultArgument,PyDefaultArgument
    def __init__(
            self, openFiles=[], sampleFiles=[], sampleArray=[], openArray=[], scaledArray=[], openHeaders=[], sampleHeaders=[], scaledHeaders=[], openNames=[], sampleNames=[]):
        self.getFiles = DirectoryHandler()
        self.openFiles = openFiles
        self.sampleFiles = sampleFiles
        self.scaledArray = scaledArray
        self.scaledHeaders = scaledHeaders
        self.openArray = openArray
        self.openHeaders = openHeaders
        self.sampleHeaders = sampleHeaders
        self.sampleArray = sampleArray
        self.sampleNames = sampleNames
        self.openNames = openNames
        self.pathToOpen = self.getFiles.openPath
        self.pathToSample = self.getFiles.samplePath

    def getOpenPath(self):
        self.getFiles.openOpenDirectory()
        # self.pathToOpen = self.getFiles.openPath
        self.listFits(openPath)

    def getSamplePath(self):
        self.getFiles.openSampleDirectory()
        # self.pathToSample = self.getFiles.samplePath
        self.listFits(samplePath)

    def listFits(self, path):

        if path == openPath:

            if os.path.exists(openPath + "/scaledOpenBeam"):
                openf = glob.glob(os.path.join(path+"/scaledOpenBeam", "*[0-9][0-9][0-9][0-9][0-9].fits"))
                for fitsfile in openf:

                    data = fits.getdata(fitsfile)
                    header = fits.getheader(fitsfile)
                    self.scaledArray.append(data)
                    self.scaledHeaders.append(header)
                    print fitsfile

            else:
                openf = glob.glob(os.path.join(path, "*[0-9][0-9][0-9][0-9][0-9].fits"))
                for fitsfile in openf:
                    hdulist = fits.open(fitsfile)
                    name = hdulist.filename().split("\\")[-1]
                    data = fits.getdata(fitsfile)
                    header = fits.getheader(fitsfile)
                    hdulist.close()
                    self.openArray.append(data)
                    self.openHeaders.append(header)
                    self.openNames.append(name)
                    print name

        else:

            if os.path.exists(samplePath + "/overlapCorrected"):
                samplef = glob.glob(os.path.join(samplePath+"/overlapCorrected", "*[0-9][0-9][0-9][0-9][0-9].fits"))
                for sf in samplef:
                    hdulist = fits.open(sf)
                    name = hdulist.filename().split("\\")[-1]
                    data = fits.getdata(sf)
                    header = fits.getheader(sf)
                    hdulist.close()
                    self.sampleArray.append(data)
                    self.sampleHeaders.append(header)
                    self.sampleNames.append(name)
                    print name

            else:
                samplef = glob.glob(os.path.join(samplePath, "*[0-9][0-9][0-9][0-9][0-9].fits"))
                for sf in samplef:
                    hdulist = fits.open(sf)
                    name = hdulist.filename().split("\\")[-1]
                    data = fits.getdata(sf)
                    header = fits.getheader(sf)
                    self.sampleArray.append(data)
                    self.sampleHeaders.append(header)
                    self.sampleNames.append(name)
                    print name

    def getDir(self):
        print self.scaledArray[0]
        print self.scaledHeaders[0]
        print self.openArray[0]
        print self.openHeaders[0]
        print self.sampleArray[0]
        print self.sampleHeaders[0]


class OverlapCorrectionAndScaling:

    def __init__(self):
        self.directories = DirectoryHandler()
        self.files = FileHandler()
        # self.pathToOpen = openPath
        # self.pathTosample = samplePath
        self.openArray = self.files.openArray
        self.openHeaders = self.files.openHeaders
        self.sampleArray = self.files.sampleArray
        self.sampleHeaders = self.files.sampleHeaders
        self.openNames = self.files.openNames
        self.sampleNames = self.files.sampleNames
        # self.ctuple = tuple(self.sampleArray)

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

        if path == samplePath:

            shutterIndices = self.preBinData(samplePath)
            shutterValues = self.readShutter(samplePath)[0]

            if os.path.exists(path + "/overlapCorrected"):

                return ctypes.windll.user32.MessageBoxA(0, "Corrected files already exist", "Error", 1)

            os.mkdir(path + "/overlapCorrected")
            f = open(path + "/overlapCorrected/TOFData.csv", "wb")
            zipped = zip(self.sampleArray, self.sampleHeaders, self.sampleNames)
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
            shutterValuesSample = self.readShutter(samplePath)[0]

            zipShutters = zip(shutterValuesOpen, shutterValuesSample)
            ratio = []
            for svo, svs in zipShutters:
                ratio.append(float(svs[1]) / float(svo[1]))

            if os.path.exists(path + "/scaledOpenBeam"):
                return ctypes.windll.user32.MessageBoxA(0, "Scaled files already exist", "Error", 1)

            os.mkdir(path + "/scaledOpenBeam")
            f = open(path + "/scaledOpenBeam/TOFData.csv", "wb")
            zipped = zip(self.openArray, self.openHeaders, self.openNames)
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

    def do(self):

        print self.sampleNames




class ShowData:

    def __init__(self, root):

        # global a
        # a = None
        # global b
        # b = None
        # global c
        # c = None
        # global d
        # d = None

        self.transdata = FileHandler()
        self.scaledArray = self.transdata.scaledArray

        self.files = FileHandler()
        self.sampleArray = self.files.sampleArray
        """
        if self.sampleArray == []:
            self.sa = np.random.rand(512, 512)
            self.sb = np.random.rand(512, 512)
            self.sc = np.random.rand(512, 512)
            self.ctuple = (self.sa, self.sb, self.sc)
            self.cube = np.dstack(self.ctuple)
        else:

            self.cube = np.empty((len(self.sampleArray), 512, 512))
            for n in range(len(self.sampleArray)):
                self.cube[n] = self.sampleArray[n]
                print n
        """
            # self.ctuple = tuple(self.sampleArray)
            # self.cube = np.dstack(self.ctuple)
        self.root = root
        self.fig = Figure(figsize=(6, 6))
        self.ax = self.fig.add_subplot(111)
        self.plotted = False
        self.l = None
        self.canvas = None
        plt.show()

    def grayifyCmap(self, cmap):
        cmap = plt.cm.get_cmap(cmap)
        colors = cmap(np.arange(cmap.N))
        RGB_weight = [0.299, 0.587, 0.114]
        luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
        colors[:, :3] = luminance[:, np.newaxis]

        return cmap.from_list(cmap.name + "_grayscale", colors, cmap.N)

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
            self.root, from_=0, to=len(self.sampleArray)-1, resolution=1, orient=tk.HORIZONTAL, command=self.update
            )
        self.slider.pack()
        self.plotted = True
        self.s = 0
        im = self.histeq(self.sampleArray[self.s])[0]
        self.l = self.ax.imshow(im, cmap=plt.cm.gray, **kwargs)
        self.canvas = FigureCanvasTkAgg(self.fig, self.root)
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()
        self.myrectsel = MyRectangleSelector(self.ax, self.onSelect, drawtype="box", rectprops=dict(
            facecolor="red", edgecolor="black", alpha=0.2, fill=True))

    def update(self, val):
        ind = int(self.slider.get())
        if self.plotted:
            im = self.histeq(self.sampleArray[ind])[0]
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


class TransPlot:

    def __init__(self, val):

        try:

            self.val = val.get()
            if self.val == "Default flight path: 56m":
                self.L = float(self.val.split(':')[1].strip('m'))

            else:
                self.L = float(self.val.strip('m'))

            self.data = FileHandler()
            self.scaledArray = self.data.scaledArray
            self.sampleArray = self.data.sampleArray

            self.a = a
            self.b = b
            self.c = c
            self.d = d

        except NameError:
            return ctypes.windll.user32.MessageBoxA(0, "You must select an ROI first", "Error", 1)

    def produceTransData(self):

        scaledIntensities = []
        for scaled in self.scaledArray:
            sliced = scaled[self.a:self.b, self.c:self.d]
            summedScaled = sum(sum(sliced))
            scaledIntensities.append(summedScaled)

        sampleIntensities = []
        for sample in self.sampleArray:
            sliced = sample[self.a:self.b, self.c:self.d]
            summedSample = sum(sum(sliced))
            sampleIntensities.append(summedSample)

        transmitted = []
        zipped = zip(sampleIntensities, scaledIntensities)
        for sample, scaled in zipped:
            transmitted.append(float(sample)/float(scaled))

        TOF = []

        TOFdata = open(openPath + '/scaledOpenBeam/TOFData.csv', 'rb+')  # grabs the TOF data
        for line in TOFdata:
            TOF.append(line.split(',')[0])

        return TOF, transmitted

    def plotTrans(self, *args):

        if len(args) > 0:
            data = self.produceTransData()

            global timeOF
            timeOF = data[0]

            ymin = min(data[1]) - 0.05
            ymax = max(data[1]) + 0.05

            self.fig = plt.figure(1)
            self.ax = self.fig.add_subplot(111)
            # self.ax.set_ylim([ymin, ymax])
            self.ax.autoscale(enable=True, axis="both", tight=True)
            self.myrectsel = MyRectangleSelector(self.ax, self.onSelect, drawtype='box', rectprops=dict(facecolor='red', edgecolor='black', alpha=0.5, fill=True))

            plt.plot(data[0], data[1])
            plt.ylim(ymin, ymax)
            plt.xlabel("Time of Flight (s)")
            plt.ylabel("Neutron Transmission")
            plt.show()
            plt.close()
            data = None

            return timeOF

        else:
            data = self.produceTransData()
            global transW
            transW = data[1]
            print self.L
            global wavelength
            wavelength = []

            h = 6.6E-34
            m = 1.67E-27
            for point in data[0]:
                wavelength.append(((h*float(point))/(self.L*m))*(10**10))

            ymin = min(data[1]) - 0.05
            ymax = max(data[1]) + 0.05

            self.fig = plt.figure(1)
            self.ax = self.fig.add_subplot(111)
            # self.ax.set_ylim([ymin, ymax])
            self.ax.autoscale(enable=True, axis="both", tight=True)
            self.myrectsel = MyRectangleSelector(self.ax, self.onSelect, drawtype='box',
                                                 rectprops=dict(facecolor='red', edgecolor='black', alpha=0.5,
                                                                fill=True))

            plt.plot(wavelength, data[1])
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


class MyRectangleSelector(RectangleSelector):

    def release(self, event):
        super(MyRectangleSelector, self).release(event)
        self.to_draw.set_visible(True)
        self.canvas.draw()


class EdgeFitting:

    def __init__(self):

        self.xvalsW = wavelength
        self.trans = transW
        self.subx = []
        self.suby = []

        self.frame = tk.Toplevel()
        self.fig = Figure(figsize=(5, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0)

        #self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame)
        #self.toolbar.update()

        self.plotButton = tk.Button(self.frame, text="Fit Curve", command=self.fitCurve)

        self.coeff1 = tk.Entry(self.frame, width=10)
        self.coeff1label = tk.Label(self.frame, text=u"      Edge Pedestal C \u2081")
        self.coeff2 = tk.Entry(self.frame, width=10)
        self.coeff2label = tk.Label(self.frame, text=u"        Edge Height C \u2082")
        self.lambda0var = tk.Entry(self.frame, width=10)
        self.lambda0label = tk.Label(self.frame, text=u"      Lambda Edge \u03BB \u2080")
        self.sigmavar = tk.Entry(self.frame, width=10)
        self.sigmalabel = tk.Label(self.frame, text=u" Bragg Edge Width \u03C3")
        self.tauvar = tk.Entry(self.frame, width=10)
        self.taulabel = tk.Label(self.frame, text=u"   Edge Asymmetry \u03C4")

        self.widgets()

        self.c1 = self.coeff1.get()
        self.c2 = self.coeff2.get()
        self.lambda0 = self.lambda0var.get()
        self.sigma = self.sigmavar.get()
        self.tau = self.tauvar.get()

    def widgets(self):

        self.plotButton.grid(row=1)

        self.coeff1label.grid(sticky="W")
        self.coeff1.grid(row=2)
        self.coeff1.insert(0, "1")
        self.coeff2label.grid(sticky="W")
        self.coeff2.grid(row=3)
        self.coeff2.insert(0, "1")
        self.lambda0label.grid(sticky="W")
        self.lambda0var.grid(row=4)
        self.lambda0var.insert(0, "1")
        self.sigmalabel.grid(sticky="W")
        self.sigmavar.grid(row=5)
        self.sigmavar.insert(0, "1")
        self.taulabel.grid(sticky="W")
        self.tauvar.grid(row=6)
        self.tauvar.insert(0, "1")

    def func(self, x, c_1, c_2, lambda0, sigma, tau):
        return c_1 * (scipy.special.erfc((lambda0 - x) / (np.sqrt(2) * sigma)) - np.exp(
            ((lambda0 - x) / tau) + (sigma ** 2 / (2 * tau ** 2))) * scipy.special.erfc(
            ((lambda0 - x) / (np.sqrt(2) * sigma)) + sigma / (np.sqrt(2) * tau))) + c_2

    def subPlot(self):

        if self.xvalsW != []:

            zipped = zip(self.xvalsW, self.trans)
            for xval, yval in zipped:
                if xval >= atp and xval <= btp:
                    self.subx.append(xval)
                    self.suby.append(yval)

            self.ax.plot(self.subx, self.suby, 'x')

    def fitCurve(self):

        warnings.simplefilter("error", OptimizeWarning)

        try:

            self.ax.cla()
            self.ax.plot(self.subx, self.suby, 'x')

            initial_guess = [float(self.coeff1.get()), float(self.coeff2.get()), float(self.lambda0var.get()), float(self.sigmavar.get()), float(self.tauvar.get())]

            popt, pcov = curve_fit(self.func, self.subx, self.suby, p0=initial_guess)
            self.ax.plot(self.subx, self.func(self.subx, popt[0], popt[1], popt[2], popt[3], popt[4]))
            self.canvas.show()
            self.clearText()
            self.coeff1.insert(0, popt[0])
            self.coeff2.insert(0, popt[1])
            self.lambda0var.insert(0, popt[2])
            self.sigmavar.insert(0, popt[3])
            self.tauvar.insert(0, popt[4])
        except (RuntimeError, OptimizeWarning):
            self.ax.cla()
            self.ax.plot(self.subx, self.suby, 'x')
            x = np.linspace(self.subx[0], self.subx[-1], 100)
            self.ax.plot(x, self.func(x, initial_guess[0], initial_guess[1], initial_guess[2], initial_guess[3], initial_guess[4]))
            self.canvas.show()
            return ctypes.windll.user32.MessageBoxA(0, "Please refine your parameters", "Error", 1)

    def clearText(self):
        fields = [self.coeff1, self.coeff2, self.lambda0var, self.sigmavar, self.tauvar]
        for field in fields:
            field.delete(0, "end")


class Test:

    def __init__(self, val):

        self.val = val.get()

    def do(self):
        print self.val

if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("500x600")
    root.title("Bragg Edge Analysis")
    app = BraggEdgeAnalysisGUI(root)
    root.mainloop()
    root.destroy()
