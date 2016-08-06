import Tkinter as tk
import glob
import os
from tkFileDialog import askdirectory
from astropy.io import fits

class BraggEdgeAnalysisGUI:

    def __init__(self, root_):

        self.root = root_
        self.frame = tk.Frame(self.root)
        self.frame.pack()

        self.directory = GetDirectories()
        self.test = Test(self.directory)

        self.openButton = tk.Button(self.frame, text="open", command=self.directory.getOpenPath)

        self.testButton = tk.Button(self.frame, text="test", command= self.test.do

        self.widgets()

    def widgets(self):
        
        self.openButton.pack()


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

    def getOpenPath(self):
        
        self.directory.openOpenDirectory()
        openPath = self.directory.openPath
        if os.path.exists(os.path.join(openPath, "scaledOpenBeam")):
            path = os.path.join(openPath, "scaledOpenBeam")
            self.loadData(path, self.openFits)
        else:
            self.loadData(openPath, self.openFits)
            
        print self.openFits.headers
        
        

    def getSamplePath(self):

        self.directory.openSampleDirectory()
        samplePath = self.directory.samplePath
        if os.path.exists(os.path.join(samplePath, "overlapCorrected")):
            path = os.path.join(samplePath, "overlapCorrected")
            self.loadData(path, self.sampleFits)
        else:
            self.loadData(samplePath, self.sampleFits)

        print self.directory.samplePath

    def loadData(self, path, container):
        f = glob.glob(os.path.join(path, "*[0-9][0-9][0-9][0-9][0-9].fits"))
        for fitsFile in f:
            hdulist = fits.open(fitsFile)
            name = hdulist.filename().split("\\")[-1]
            header = hdulist[0].header
            data = hdulist[0].data
            hdulist.close()
            container.names.append(name)
            container.headers.append(header)
            container.arrays.append(data)
            print name

class Test:

    def __init__(self, dir):
        self.a = dir

    def do(self):
        print self.a

    

if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("500x600")
    root.title("Bragg Edge Analysis")
    app = BraggEdgeAnalysisGUI(root)
    root.mainloop()
    root.destroy()
