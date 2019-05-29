""" Save user settings and analyse images

Refer to comments in function. 

"""

from __future__ import with_statement, division
from ij import IJ, ImagePlus, WindowManager, Prefs
from ij.plugin import ImageCalculator, ChannelSplitter, Duplicator
from ij.process import ImageConverter, ImageProcessor
import os
import csv
from ij import IJ, ImagePlus
from ij.gui import DialogListener, GenericDialog, ImageWindow, \
                   MessageDialog, YesNoCancelDialog
from java.awt.event import ActionListener, ActionEvent, \
     ItemListener, ItemEvent, WindowAdapter
from ij.plugin.frame.Editor import actionPerformed
from java.awt import Button, Checkbox, TextField, \
                     BorderLayout, Font, Color, \
                     Dimension, Cursor
from java.util import Random, Vector
from java.lang import String
from javax.swing import JFrame, JButton, JOptionPane, \
                        BorderFactory, JFrame, JLabel, JPanel, \
                        JTextArea, JSlider, JCheckBox, JComboBox, \
                        JTextField, JScrollBar
from loci.common.DebugTools import enableIJLogging
from ij.plugin.filter import RankFilters
from javax.swing.AbstractButton import setBorderPainted
from ij.io import FileSaver
import time
from javax.swing import ImageIcon
import sys
from ij.macro import Interpreter as IJ1
import mpicbg.ij.clahe.Flat
from java.lang import Runtime, System
from inra.ijpb.morphology.attrfilt import BoxDiagonalOpeningQueue
from inra.ijpb.morphology import Morphology
w = WindowManager
OS = System.getProperty("os.name")


def closeimage():
        """ Close active image
        """
        
        images = w.getImageCount()
        if images != 0:
                IJ.run("Close")


def closeallimages():
        """ Close all active windows
        """
        
        images = w.getImageCount()
        if images != 0:
                IJ.run("Close All")


def newUser(cwd, greyscaleMinVal, g, r, backgroundsubRolling, cellbodycb, user, SN, mCLAHE, backgroundsubNeurite, setpixels,
            statcb, threshChoice, despeckle, contrast, Min, Max, radius, Min2, Max2, 
            threshChoice2, mCLAHE2, Sbgcbstate):
            
            """ Save user settings
            Creates a comma separated values (csv) file called "username".csv which contains
            all of the image analysis settings defines by the user. The file is saved in the
            MyelinJ folder within ImageJ.
            """
           
            closeimage()
            # get all user defined settings and make a CSV file, where its name is the user name
            imagesettings1 = []
            imagesettings1.append(Min)
            imagesettings1.append(Max)
            imagesettings1.append(threshChoice)
            imagesettings1.append(despeckle)
            imagesettings1.append(g)
            imagesettings1.append(r)
            imagesettings1.append(backgroundsubRolling)
            imagesettings2 = []
            if cellbodycb is True:
                imagesettings2.append(radius)
            else:
                imagesettings2.append("0")
            imagesettings2.append(mCLAHE)
            imagesettings2.append(backgroundsubNeurite)
            imagesettings2.append(setpixels)
            imagesettings2.append(greyscaleMinVal)
            if SN is False:
                imagesettings2.append(contrast)
            else:
                imagesettings2.append("0")
            imagesettings2.append(cellbodycb)
            
            if SN is True:
                imagesettings3 = []
                imagesettings3.append(Sbgcbstate)
                imagesettings3.append(mCLAHE2)
                imagesettings3.append(threshChoice2)
                imagesettings3.append(Min2)
                imagesettings3.append(Max2)
                imagesettings3.append("0")
                imagesettings3.append("0")
            
            totalsettings = []
            totalsettings.append(imagesettings1)
            totalsettings.append(imagesettings2)
            if SN is True:
                totalsettings.append(imagesettings3)
            
            root = cwd
            filename = user
            fullpath = os.path.join(root, filename)
            f = open(fullpath, 'wb')
            writer = csv.writer(f)
            with open(user, "wb"):
                    writer = csv.writer(f)
                    writer.writerows(totalsettings)
            f.close()


class Finished(JFrame):

                def __init__(self):
                    super(Finished, self).__init__()
                    self.initUI()

                def initUI(self):
                    """ Finished Dialog box
                    Simple dialog box that says "Finished", to
                    bw displayed when all image analysis has
                    finished. When OK button is pressed all
                    ImageJ windows are closed.
                    """
                    
                    panel = JPanel()
                    self.getContentPane().add(panel)
                    panel.setBackground(Color.WHITE)
                    panel.setLayout(None)
                    self.setTitle("Analysis has finished")
                    self.setSize(300, 150)
                    OKbutton = JButton("OK", actionPerformed=self.onOK)
                    OKbutton.setBackground(Color.BLACK)
                    OKbutton.setBounds(80, 50, 100, 30)
                    panel.add(OKbutton)
                    Title = JTextArea("Analysis has finised!! :-)")
                    Title.setBounds(15, 10, 250, 20)
                    panel.add(Title)
                    self.setLocationRelativeTo(None)
                    self.setLocation(int(IJ.getScreenSize().width * 0.01),
                                        int(IJ.getScreenSize().height *3 /10))
                    self.setVisible(True)

                def onOK(self, cb2):
                        self.dispose()
                        closeallimages()


def analyse(cwd, user, imagefolder, stats, experiments, multi, Rloc2, subfoldernames, names, statsfolderPath, cwdR):
        """ Main image analysis
        Gets user image analysis settings from the .csv file.
        If multiple experiments have been selected by the user
        (multi) each subfolder will be looped through. A nested
        loop will then interate through each .tif image and
        analyse. A .csv file will be produced for each folder
        analysed with the name of each image and its % neurite
        density and % myelination. A summary csv file will also
        be produced with the average % neurite density and %
        myelination for each subfolder. If statistical analysis
        has been selected (stats) then MyelinJ's Rscript will be
        run via the command line. If multple experiments is not
        selected then all of the images within the selected
        folder will be analysed together and no summary .csv will
        be produced.
        Independ of the analysis settings defined, a processed
        myelin channel image and a processed neurite channel
        image will be saved. The images can be any number of
        subdirectories (folders within folders).
        Parameters
        ----------
        cwd : string
            Path for current working directory (location of
            MyelinJ folder in Fiji).
        user: string
            User name
        imagefolder: string
            Path to .tiff image folder(s) defined by user.
        stats: boolean
            Perform statistical analysing using R?
        experiments: 2D list of strings
            list of all the subfolders (experiments) that are in each
            experimental condition.
        multi: boolean
            Analyse multiple experiments?
        Rloc2: string
            file path to Rscript location
        subfoldernames: string
            name of each subfolder which denoates each individual
            experiment, if multple experiments are being analysed.
        names: array
            array of textfields for each experimental condition defined by
            user. User will enter the name of each experimental condition.
        statsfolderPath: string
            file path to the create statsfolder.
        cwdR: string
            file path to MyelinJstats.R
        """
        # read settings from the user name CSV
        bg = False
        readsettings = []
        imagenames = []
        neuritedensity = []
        myelinoverlay = []
        myelinaverage2 = []
        neuriteaverage2 = []
        root = cwd
        filename = user
        fullpath = os.path.join(root, filename)
        f = open(fullpath, 'rb')
        readCSV = csv.reader(f)
        for row in readCSV:
                readsettings.append(row[0])
                readsettings.append(row[1])
                readsettings.append(row[2])
                readsettings.append(row[3])
                readsettings.append(row[4])
                readsettings.append(row[5])
                readsettings.append(row[6])
        f.close()
        i = 0
        
        for i in range(len(subfoldernames)):
            # if multiple experimental conditions has been selected each folder is treated as a
            # separate experiment and looped through separately otherwise all folders will be
            # treated as one experiment this only works for sub directories within the main folder.
            # Further folders will be ignored (each image can be in its own folder for example)
            if multi is True:
                 # if multiple experiments are being analysed the file path is changed to the
                 # current subfolder
                 settings2 = os.path.join(imagefolder, subfoldernames[i])
                 if "Windows" in OS:
                     settings2 = settings2 + "\\"
                 elif "Mac" in OS:
                    settings2 = settings2 + "/"
            else:
                 settings2 = imagefolder
             # loop through all .tiff files in location
            for root, dirs, files in os.walk(settings2):
              for name in files:
                if name.endswith((".tif")):
                    imagenames.append(os.path.join(name))
                    # open .tiff image, split channels and
                    # convert to 8bit grey scale.
                    imp = IJ.openImage(os.path.join(root, name))
                    g = int(readsettings[4])
                    r = int(readsettings[5])
                    imp = ChannelSplitter.split(imp)
                    green = imp[g]
                    red = imp[r]
                    conv = ImageConverter(red)
                    conv.convertToGray8()
                    conv = ImageConverter(green)
                    conv.convertToGray8()
                    
                    # thresholding to select cell bodies
                    green2 = green.duplicate()
                    if (readsettings[0] != "0") or (readsettings[1] != "0"):
                        bg = True    
                        IJ.setAutoThreshold(green2, readsettings[2])
                        IJ.setRawThreshold(green2, int(readsettings[0]), int(readsettings[1]), None)
                        Prefs.blackBackground = True
                        IJ.run(green2, "Convert to Mask", "")
                        IJ.run(green2, "Invert LUT", "")
                        if readsettings[7] != "0":
                               IJ.run(green2, "Make Binary", "")
                               IJ.run(green2, "Remove Outliers...", "radius="+readsettings[7]+" threshold=50 which=Bright")
                    
                    # CLAHE and background subtraction
                    if readsettings[8] == "True":
                        mpicbg.ij.clahe.Flat.getFastInstance().run(green, 127, 256, 3, None, False)
                    if readsettings[9] == "True":
                        calc = ImageCalculator()
                        green = calc.run("Subtract create", green, red)
                        IJ.run(green, "Subtract...", "value="+readsettings[10])
                    elif readsettings[6] == "True":
                        IJ.run(green, "Subtract Background...", "rolling=50")
                    
                    # run frangi vesselness
                    IJ.run(green, "Frangi Vesselness (imglib, experimental)", "number=1 minimum=0.454009 maximum=0.454009")
                    green = IJ.getImage()
                    
                    # convert frangi vesselness image to 8bit grey scale
                    conv = ImageConverter(green)
                    conv.convertToGray8()
                    IJ.run(green, "Convert to Mask", "")
                    
                    # remove cell bodies
                    if bg is "True":
                        green = ImageCalculator().run("Subtract create", green, green2)
                    
                    # run grey scale morphology filter from MorpholibJ
                    if readsettings[11] != "0":
                        green = green.getProcessor()
                        algo = BoxDiagonalOpeningQueue()
                        algo.setConnectivity(4)
                        result = algo.process(green, int(readsettings[11]))
                        green = ImagePlus("result", result)
                    IJ.run(green, "Invert LUT", "")
                    
                    
                    if len(readsettings) > 14:
                        # sparse neurite image analysis
                        if readsettings[15] == "True":                        
                            IJ.run(red, "Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)")
                        if readsettings[14] == "True":
                            IJ.run(red, "Subtract Background...", "rolling=50")
                        IJ.setAutoThreshold(red, readsettings[16])
                        IJ.setRawThreshold(red, int(readsettings[17]), int(readsettings[18]), None)
                        IJ.run(red, "Convert to Mask", "")
                    else:
                        # dense neurite image analysis
                        IJ.run(red, "Normalize Local Contrast", "block_radius_x=40 block_radius_y=40 standard_deviations="+readsettings[12]+" center stretch")
                        IJ.run(red, "Auto Threshold", "method=Default white")
                        IJ.run(red, "Invert LUT", "")
                    if readsettings[3] == "True":
                            IJ.run(red, "Despeckle", "")
                            
                    IJ.saveAs(red, "Jpeg", settings2+name+"neurites")
                    # get number of neurite pixels
                    
                    # get number of neurite pixels
                    statsneurite = red.getProcessor()
                    statsneurite = statsneurite.getHistogram()
                    neuritedensity.append(statsneurite[255])
                    IJ.saveAs(green, "Jpeg", settings2+name+"myelinFinal")
                    
                    # get number of myelin pixels
                    statsmyelin = green.getProcessor()
                    statsmyelin = statsmyelin.getHistogram()
                    myelinoverlay.append(statsmyelin[255])
                    closeallimages()
                    
                    # get pixel total of image
                    whitepixels = (statsneurite.histogram[0])
                    blackpixels = (statsneurite.histogram[255])
                    
            totalpixels = whitepixels + blackpixels
            totalpixels = [totalpixels]*len(neuritedensity)
            
            # for each image calculate % myelination as number of myelin pixels
            # divided by the number of neurite pixels * 100
            myelinoverlay = [x1/x2*100 for (x1, x2) in zip(myelinoverlay, neuritedensity)]
            myelinaverage = sum(myelinoverlay)/len(myelinoverlay)
            myelinaverage2.append(myelinaverage)
            
            # for each image calculate % neurite density as neurite pixels divided
            # by the total number of pixels in the image * 100.
            neuritedensity = [x1/x2*100 for (x1, x2) in zip(neuritedensity, totalpixels)]
            neuriteaverage = sum(neuritedensity)/len(neuritedensity)
            neuriteaverage2.append(neuriteaverage)
            name = "Image names"
            green = "% myelination"
            red = "% neurite density"
            imagenames = [name]+imagenames
            neuritedensity = [red]+neuritedensity
            myelinoverlay = [green]+myelinoverlay
            result = []
            result.append(imagenames)
            result.append(neuritedensity)
            result.append(myelinoverlay)
            
            root = settings2
            filename = "Results.csv"
            fullpath = os.path.join(root, filename)
            f = open(fullpath, 'wb')
            writer = csv.writer(f)
            for d in range(len(result)):
                row = [result[d]]
                writer.writerows(row)
            f.close()
            
            # must be reset to 0 for each iteration.
            y = 0
            r = 0
            
            # if statistical analysis is being performed the results .csv file
            # is also saved to a subfolder within the statistical analysis folder
            # which denotes the experimental condition the results belong to.
            if stats is True:
                # nested for loop to identify correct experimental condition
                # for the current subfolder being analysed.
                for y in range(0, len(experiments)):
                    for r in range(0, len(experiments[0])):
                        if experiments[y][r] == subfoldernames[i]:
                            if "Windows" in OS:
                                root = imagefolder+"\\statistical analysis\\"+names[y].getText()
                            elif "Mac" in OS:
                                root = imagefolder+"/statistical analysis/"+names[y].getText()
                            filename = subfoldernames[i]+".csv"
                            fullpath = os.path.join(root, filename)
                            f = open(fullpath, 'wb')
                            writer = csv.writer(f)
                            for e in range(len(result)):
                                row = [result[e]]
                                writer.writerows(row)
                            f.close()
            
                            break
            cwd2 = os.getcwd()
            for files in os.listdir(cwd2):
                    if files.endswith(".csv"):
                        os.remove(os.path.join(cwd2, files))
            imagenames = []
            myelinoverlay = []
            neuritedensity = []
            
        # create .csv summary sheet with average % neurite density
        # and average % myelination for each subfolder (experiment).
        if multi is True:
            name = "Folder name"
            imagenames = [name]+subfoldernames
            neuritedensity = [red]+neuriteaverage2
            myelinoverlay = [green]+myelinaverage2
            result = []
            result.append(imagenames)
            result.append(neuritedensity)
            result.append(myelinoverlay)
            if "Windows" in OS:
                root = imagefolder+"\\"
            elif "Mac" in OS:
                root = imagefolder+"/"
            filename = "Result-Summary.csv"
            fullpath = os.path.join(root, filename)
            f = open(fullpath, 'wb')
            writer = csv.writer(f)
            for p in range(len(result)):
                row = [result[p]]
                writer.writerows(row)
            f.close()
            imagenames = []
            myelinoverlay = []
            neuritedensity = []
            
        # Run Rscript for statistical analysis via the command line
        if stats is True:
            cmd = Rloc2+" "+cwdR+" "+statsfolderPath
            Runtime.getRuntime().exec(cmd)
        Finished()
