
"""MyelinJ graphical user interface
This script provides a GUI to set and adjust settings for the analysis of
micrographs of myelinating cultures. The GUI accepts images in .tif format
where the neurite and myelin micrographs have been merged together.
The GUI saves all of the user defined settings in a comma separated values
(.csv) file using the imported module MyelinJanalysis.
Some basic functions are also imported from the module basicfunctions
This script contains the following functions:
    * analysed - saves settings in .csv file and runs image analysis
    * getNext - this allows the user to select the next image in a large
    file of images, so analysis settings can be checked on mutiple images
    * getimage - this opens an image from the file selected by the user.
    Can be used in conjunction with getNext to open the next image in a file.
This scipt contains the folowing Dialog boxes:
    * Dialog1 - if a username has already been made then this can be
    selected and the settings will be used for image analysis. Otherwise a
    new user can be made.
    * DialogStats - the user can define the number of experimental
    conditions and which experimental condition each folder of images
    belongs to. This is required because statistical analysis requires
    at least three repeats for each experimental condition.
    * Dialog2 - enter a username
    * Dialog3 - define which channel in merged .tif corresponds to neurites
    and which corresponds to myelin.
    * Dialog4 - define analysis settings for myelin micrographs
    * Dialog5 - define analysis settings for neurite micropgraphs
    * Dialog6 - alternative settings for the analysis of neurite micrographs
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
SN = False
w = WindowManager
OS = System.getProperty("os.name")
cwd = os.getcwd()

# set folder containing the macro as the current working directory i.e
# location for username CSV (which  should be in Fiji's plugin file) based
# on the operating system bein used
if "Windows" in OS:
    cwdR = cwd+"\\plugins\\MyelinJ\\MyelinJstats.R"
    cwdR = '"'+cwdR+'"'
    cwd = cwd+"\\plugins\\MyelinJ\\"
elif "Mac" in OS:
    cwdR = cwd + "/plugins/MyelinJ/MyelinJstats.R"
    cwd = cwd+"/plugins/MyelinJ/"
sys.path.append(cwd)

import MyelinJanalysis
import basicfunctions
import config

# select file containing images for analysis
imagefolder = IJ.getDirectory("Choose a Directory")
basicfunctions.closeallimages()  # close any images already open
username = []
usernamepath = []
g = 0
r = 0

def neuritesubtract():
        """ Subtract neurites from myelin channel
        duplicates the current myelin channel image, subtracts neurites
        (red) using the ImageCalculator and then displays image.
        """
        green3 = IJ.getImage()
        green3 = green3.duplicate()
        green3 = ImageCalculator().run("Subtract create", green3, red)
        ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                        int(IJ.getScreenSize().height * 1/14))
        basicfunctions.closeimagebg()
        green3.show()


def analysed():
        """
        Runs image analysis using settings defined by the user. If a new
        user name has been created a .csv file will first be created.
        Analyse will then read settings from the user name .csv file and
        perform the analysis. The imported module MyelinJ analysis is used
        for this function.
        """
        if config.newusercb is True:
            MyelinJanalysis.newUser(cwd, config.greyscaleMinVal, g, r, config.backgroundsubRolling,
                                     config.cellbodycb, config.user, config.SN, config.mCLAHE,
                                     config.backgroundsubNeurite, config.setpixels, config.stats,
                                     config.threshChoice, config.despeckle, config.contrast,
                                     config.Min, config.Max, config.radius, config.Min2,
                                     config.Max2, config.threshChoice2, config.mCLAHE2,
                                     config.Sbgcbstate)
                                     
        MyelinJanalysis.analyse(cwd, config.user, imagefolder, config.stats,
                                 config.experiments, config.multi, config.RscriptPath,
                                 config.subfoldernames, config.names, config.statsfolderPath,
                                 cwdR)

      
def getNext():
    """
        Gets the next image in the selected folder. Goes back to the first
        image if all images have been opened.
        Parameters
        ----------
        imagecount: int
            Total number of images in selected folder. From Dialog1.onOK.
        returns (global in config)
        ----------
        imageposition : int
            Position of the current image in the selected folder.
        """
    config.imageposition = config.imageposition + 1
    if config.imageposition > config.imagecount - 1:
        config.imageposition = 0


def getimage():
        """
        Opens image in folder. Position of the image opened in the folder
        is defined by the function getNext. Alternativly, if the user has
        chosen to use an image they have opened themselves, then this image
        will be selected as the current image.
        Parameters
        ----------
        userimage2 : checkbox
            if True then the user has opened an image themselves to test
            the settings on rather than using an image from the folder they
            selected.
        userimagename : imp
            reference to image opened by the user. This image is the current
            active image.
        g: int
            position of the myelin micrograph after the image has been
            split into separate channels.
        r: int
            position of the myelin micrograph after the image has been
            split into separate channels.
        """
        global red, green
        if config.userimage2 is True:
            imp = userimagename
        else:
            imp = IJ.openImage(config.listAllImages[config.imageposition])
        imp = ChannelSplitter.split(imp)
        green = imp[g]
        conv = ImageConverter(green)
        conv.convertToGray8()
        red = imp[r]
        conv = ImageConverter(red)
        conv.convertToGray8()
        IJ.run(red, "Grays", "")
        green.setTitle("original")
        ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                            int(IJ.getScreenSize().height * 1/14))


def applybackground():
         """ background subtraction
          Performs any background settings defined by the user the image
          currently being displayed. Changes title of the image.
          Parameters
          ----------
          mCLAHE : bool
           Perform CLAHE?
          backgroundsubRolling: bool
              perform rolling ball background subtraction?
          backgroundsubNeurite: bool
              perform neurite subtraction
          setpixels: string
              value for pixel subtraction?
         """
         if config.mCLAHE is True:
                basicfunctions.CLAHE()
         if config.backgroundsubRolling is True:
                basicfunctions.rollingsubtract()
         elif config.backgroundsubNeurite is True:
                neuritesubtract()
         if config.setpixels != "0":
                IJ.run("Subtract...", "value="+config.setpixels)
         green = IJ.getImage()
         green.setTitle("original + background subtraction")


def removecellbodies():
            if config.cellbodycb is True:
                getimage()
                green.show()
                if (config.Min != 0) or (config.Max != 255):
                    IJ.setAutoThreshold(green, config.threshChoice+" dark")
                    IJ.setRawThreshold(green, config.Min, config.Max, None)
                    Prefs.blackBackground = True
                    IJ.run(green, "Convert to Mask", "")
                if config.radius != "0" or "":
                        IJ.run(green, "Remove Outliers...",
                                "radius="+config.radius+" threshold=50 which=Bright")
                green.setTitle("cell body selection")


def frangifilter(self):
            """ Frangi vesselness
            Performs frangi vesselness for myelin channel image.If cell
            bodies have been selected this is assumed to be the
            activeimage. A new image with any background subtraction is
            performed, any cellbodies selected are removed using image
            calculator and frangi vesselness is performed.
            Parameters
            ----------
            cellbodycb : bool
                  remove cell bodies?
            attrival2: string
                  value for performing grey scale morphology filtering.
            """
            global green
            self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
            # if cell bodies removal is selected then a new a new myelin channel
            # image is opened and any background subtraction selected by user
            # is performed.
            if config.cellbodycb is True:
                bg = IJ.getImage()
                bg = bg.duplicate()
                getimage()
                green.show()
                applybackground()
            ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                            int(IJ.getScreenSize().height * 1/14))
            self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
            IJ.run("Frangi Vesselness (imglib, experimental)",
                    "number=1 minimum=0.454009 maximum=0.454009")
            g = IJ.getImage()
            g.setTitle("original + background subtraction + vesselness")
            g = g.duplicate()
            conv = ImageConverter(g)
            conv.convertToGray8()
            IJ.run(g, "Convert to Mask", "")
            self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
            if config.cellbodycb is True:
                green = ImageCalculator().run("Subtract create", g, bg)
                ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                                int(IJ.getScreenSize().height * 1/14))
                green.show()
                green.setTitle("original + background subtraction + vesselness + cell body subtraction")
                self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
            ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                            int(IJ.getScreenSize().height * 1/14))
            if config.greyscaleMinVal != "0" or "":
                IJ.run(green, "Gray Scale Attribute Filtering", "operation=Opening attribute=[Box Diagonal] minimum="+config.greyscale.MinVal+" connectivity=4")
                green = IJ.getImage()
                green.setTitle("original + background subtraction + vesselness + cell body subtraction")
            self.setCursor(Cursor.getDefaultCursor())


basicfunctions.closeallimages()  # close any images already open
username = []
usernamepath = []
# find all of the CSV files for user names (if any) - gets the name of the
# file (user name) and the location of the file
for root, dirs, files in os.walk(cwd):
    for name in files:
        if name.endswith((".csv")):
            username.append(name)
            usernamepath.append(os.path.join(root, name))


class Dialog1(JFrame):
    def __init__(self):
        super(Dialog1, self).__init__()
        self.initUI()
        """
        Select previously made user settings and run image analysis or make
        a new user name. Select whether multiple experimental conditions or
        statistical analysis will be performed (this is not saved in the user
        name settings and so has to be defined everytime an analysis is
        performed - this is to allow more flexbility so less user names are
        required for the same analysis settings).
        Attributes
        ----------
         selectuser: combobox
            Drop down box listing all previous user names as .csv files.
        newuser : checkbox
           check if making a new user name.
        multicb: checkbox
            check if multiple experimental conditions are to be analysed at
            the same time.
        statcb: checkbox
            check if statistical analysis will be performed
        Okbutton: button
            if a user name has been selected image analysis will be started.
            Otherwise more Dialogs will open to imput required analysis
            settings.
        Cancelbutton:
            Dialog box will close and the macro will stop.
        Raises
        ------
        error - if a folder is not selected
        """
    def initUI(self):
        panel = JPanel()
        self.getContentPane().add(panel)
        panel.setBackground(Color.WHITE)
        panel.setLayout(None)
        self.setTitle("Choose user name")
        self.setSize(300, 200)

        self.selectuser = JComboBox(username)
        self.selectuser.setBounds(50, 10, 180, 20)
        panel.add(self.selectuser)

        newusercb = JCheckBox("New user?", True,
                                actionPerformed=self.onNewuser)
        newusercb.setBounds(20, 50, 100, 20)
        newusercb.setSelected(False)
        panel.add(newusercb)

        self. multicb = JCheckBox("Multiple experimental conditions?", True)
        self.multicb.setBounds(20, 70, 250, 20)
        self.multicb.setSelected(False)
        panel.add(self.multicb)

        self.statcb = JCheckBox("Perform statistics?", True)
        self.statcb.setBounds(20, 90, 170, 20)
        self.statcb.setSelected(False)
        panel.add(self.statcb)

        OKbutton = JButton("OK", actionPerformed=self.onOK)
        OKbutton.setBackground(Color.BLACK)
        OKbutton.setBounds(20, 120, 100, 30)
        panel.add(OKbutton)

        Cancelbutton = JButton("Cancel", actionPerformed=self.onCancel)
        Cancelbutton.setBackground(Color.BLACK)
        Cancelbutton.setBounds(150, 120, 100, 30)
        panel.add(Cancelbutton)

        self.setLocationRelativeTo(None)
        self.setLocation(int(IJ.getScreenSize().width * 0.01),
                              int(IJ.getScreenSize().height * 3/10))
        self.setVisible(True)
        # error if user has not selected a folder
        global imagefolder
        if imagefolder is None:
            IJ.showMessage("Error: please selected a folder first")
            imagefolder = IJ.getDirectory("Choose a Directory")
 
    def onNewuser(self, cb2):
        newusercb = cb2.getSource()
        config.newusercb = newusercb.isSelected()

    def onCancel(self, b):
        self.dispose()

    def onOK(self, cb2):
        """Runs username analysis or starts input of new user
        Starts image analysis if a previously defined user name has been
        selected. Otherwise starts the input of a new user name and
        settings. Multicb allows analysis of multiple experimental
        conditiosn simultaneously (where each conditionis analysed
        individually). statcb = True will start. Dialogstat to define
        settings for statistical analysis with R.
        Parameters
        ----------
        newuser: bool
            create a new user profile?
        returns (global in config.py)
        --------------------------
        multi: bool
            multiple experiments to be analysed?
        stats: bool
            statistical analysis to be performed by R?
        subfoldernames: list
            list of all subdirectories within a folder
        test: list
            list of all .tif images in selected folder
        imagecount: int
            total number of .tif images within selected folder
        user: string
            predefined user name selected.
        Raises
        ------
        A folder must be selected.
        The selected folder must contain at least one .tif image.
        For statistical analysis multi (for mulitple
        experimental conditions) must also be True.
        """
        global imagefolder
        config.multi = self.multicb.isSelected()
        config.stats = self.statcb.isSelected()
        if (config.stats is True) and (config.multi is False):
            IJ.showMessage("Error: multiple experimental conditions are required for statistical analysis")
        else:
            # get a list of all the images in the user selected folder if .tif
            config.subfoldernames = list()
        # error if user has not selected a folder
        if imagefolder is None:
            IJ.showMessage("Error: Please select a folder first")
            imagefolder = IJ.getDirectory("Choose a Directory")
        else:
          for root, dirs, files in os.walk(imagefolder):
                config.subfoldernames.append(dirs)
                for name in files:
                        if name.endswith((".tif")):
                            config.listAllImages.append(os.path.join(root, name))
        # error if no .tif files in selected folder
        if config.listAllImages is []:
                IJ.showMessage("Error: no .tif files found in folder")
                imagefolder = IJ.getDirectory("Choose a Directory")
        else:
            # select subdirectory names within folder. This is for analysing
            # multiple experimental conditions where each subdirectory denotes
            # a separate experimental condition
            config.subfoldernames = config.subfoldernames[0]
            if config.multi is False:
                # sets subfoldernames to 1 if multiple experiment analysis is
                # not required (means can iterate over once in a for loop).
                config.subfoldernames = [1]
            config.imagecount = len(config.listAllImages)
            if config.newusercb is True:
                self.dispose()
                config.newuser = True
                if config.stats is False:
                    Dialog2()
                else:
                        DialogStats()
            else:
                     config.user = self.selectuser.getSelectedItem()
                     self.dispose()
                     if config.stats is False:
                         analysed()
                     else:
                            DialogStats()


class DialogStats(JFrame):
    def __init__(self):
        super(DialogStats, self).__init__()
        self.initUI()
        """ Defines which experimental condition each subfolder belongs to.
        For statistical analysis each experimental condition requires a
        minimum od three experimental repeats. This Dialog box allows the
        user to define the names of each experimental repeat (which will
        appear in the grpahs produced by R) and select which folder belongs
        to which experimental condition.
        Attributes
        ----------
        Rloc: textfield
            path to Rscript
        noConditions : textfield
           Number of experimental conditions
        Methods
        -------
        A dynamic dialog box is produced depending on the number of
        subfolders (experiments) and the number of experimental conditions
        defined by user. A matrix is dispayed where subfolder names are row
        names, column names are experimental conditions (defined by the
        user) and checkboxes allow the user to define which experiments
        belong to which experimental conditions. onOK creates a 2D list that
        defines which subfolders belong to which exerimental conditions -
        based on the checkboxes selected by the user.
        Raises
        ------
        Must be at least six folders i.e two experimetal conditions with
        three repeats each.
        """
    def initUI(self):
        global panel
        panel = JPanel()
        self.getContentPane().add(panel)
        panel.setBackground(Color.WHITE)
        panel.setLayout(None)
        self.setTitle("Define experimental conditions")
        self.setSize(config.width1, config.height1)
        Min2 = JTextField("Rscript location:")
        Min2.setBounds(20, 10, 85, 20)
        Min2.setEditable(False)
        panel.add(Min2)
        if "Windows" in OS:
            self.Rloc = JTextField("C:\\Program Files\\R\\R-3.5.1\\bin\\Rscript")
        elif "Mac" in OS:
            self.Rloc = JTextField("/Library/Frameworks/R.framework/Versions/3.3/Resources/Rscript")
        self.Rloc.setBounds(110, 10, 350, 20)
        self.Rloc.setEditable(True)
        panel.add(self.Rloc)
        Min1 = JTextField("No. experimental conditions (press enter)")
        Min1.setBounds(20, 40, 250, 20)
        Min1.setEditable(False)
        panel.add(Min1)
        self.noConditions = TextField("", actionPerformed=self.onClick)
        self.noConditions.setBounds(290, 41, 40, 20)
        panel.add(self.noConditions)
        self.setLocationRelativeTo(None)
        self.setLocation(int(IJ.getScreenSize().width * 0.01),
                           int(IJ.getScreenSize().height * 3/10))
        self.setVisible(True)

    def onClick(self, cb2):
        """ Array of texfields and checkboxes for each folder (experimental
         condition)
        The user enters the number fo experimental conditions and then
        presses OK on the keyboard.
        This displays a matrix of textfields and checkboxes for each
        folder. Row names are the name of each subfolder (experiment). The
        headings of the columns are empty texfields corresponding to the
        number of experimental conditions defined by the user.The user can
        enter the names of each experimental condition (which will appear
        on any graphs produced by R later). The matrix is made up of an
        array of checkboxes which allows the user to define which
        experimental condition each folder belongs to.
        Parameters
        ----------
        repeats: int
            number of experimental repeats defind by the user.
        subfoldernames : list
           list of all subfolders within main folder (each folder
           corresponds to an experimental condition)
        heigh1: int
            height on the dialog box.
        width1:int
            width of the dialog box.
        text: array
            array of textfields for the name of each subfolder
        accross: int
            position of text (array). Increments by 30 for each textfield
            required (repeats).
        check: checkbox
            array of checkboxes Where the number will be: number of
            experimental conditions * number of subfolders.
        down: int
            position of checkbox array going down. Increments for each
            subfolder
        accross2: int
            position of checkbox array going acrross. Increments for each
            experimental condition.
        names: array
            array of textfields for each experimental condition defined by
            user. User will enter the name of each experimental condition.
        OKbutton: button
            starts image analysis if a predefined user name has already
            been selected. Otherwise continues to create a new user name.
        resetbutton: button
            sets the dialogbox back to the start - before the number of
            experimental conditions was provided.
           Raises
        ------
        Must be at least six folders i.e two experimetal conditions with
        three repeats each. Must be at least two experimental conditions.
        """
        global  panel, check, repeats
        accross = 85
        text = [0] * len(config.subfoldernames)
        repeats = int(self.noConditions.getText())

        # need to be at least six subfolders (experiments and
        # two experimental conditions defined.
        if len(config.subfoldernames) < 6:
            IJ.showMessage("Error: must be at least six experiments (i.e 3 repeats of two conditions)")
        elif repeats < 2:
            IJ.showMessage("Error: must be at least two different experimental conditions")
        elif (len(config.subfoldernames) >= 6) and (repeats >= 2):
            repeats2 = len(config.subfoldernames) - 3
            repeats2 = repeats2 * 44
            config.height1 = config.height1 + repeats2
            self.setSize(config.width1, config.height1)

            # creates an array of texfields for the name of
            # each subfolder (experiment) and adds then to
            # the dialog box as a list.
            for i in range(0, len(config.subfoldernames)):
                text[i] = JTextField(config.subfoldernames[i])
                text[i].setBounds(20, accross, 100, 20)
                panel.add(text[i])
                accross = accross + 30

            if repeats > 3:
                repeats2 = repeats - 3
                repeats2 = repeats2 * 110
                config.width1 = config.width1 + repeats2
                self.setSize(config.width1, config.height1)

            check = [[0 for x in range(len(config.subfoldernames))] for y in range(repeats)]
            down = 85
            accross2 = 25
            config.names = [0] * repeats

            # adds an array of textfields for the number of experimental
            # conditions defined by the user as headings.
            for z in range(0, repeats):
                config.names[z] = JTextField()
                accross2 = accross2 + 105
                config.names[z].setBounds(accross2, 63, 100, 20)
                panel.add(config.names[z])
                down = 85

                # adds an array of checkboxes to make a matrix
                for x in range(0, len(config.subfoldernames)):
                    check[z][x] = JCheckBox()
                    check[z][x].setBounds(accross2, down, 40, 20)
                    panel.add(check[z][x])
                    down = down + 30
            panel.repaint()

            # Place buttons at bottom of the dialog box.
            OKbutton = JButton("OK", actionPerformed=self.onOK)
            OKbutton.setBackground(Color.BLACK)
            OKbutton.setBounds(80, down, 100, 30)
            panel.add(OKbutton)
            Resetbutton = JButton("Reset", actionPerformed=self.onReset)
            Resetbutton.setBackground(Color.BLACK)
            Resetbutton.setBounds(180, down, 100, 30)
            panel.add(Resetbutton)

    def onReset(self, cb2):
        # closes and reopens dialog box.
        self.dispose()
        DialogStats()

    def onOK(self, cb2):
        """
        Continuation of onClick. retrieves the names of each experimental
        condition and which subfolder (experiment) belongs to which
        experimental condition. creates a folder called statistical
        analysis and subfolders for each experimental condition. Retrives
        the location of Rscript, as defined by the user. This is required
        for running R via the command line later on.
        Parameters
        ----------
        subfoldernames : list
           list of all subfolders within main folder (each folder
           corresponds to an experimental condition).
        repeats2: int
            used to calculate the change in width depending on the number
            of experimental conditions defined.
        text: array
            array of textfields for the name of each subfolder
        same: int
            count to determine whether any subfolders (experiments) have
            been assigned more than one experimental condition.
        check: checkbox
            array of checkboxes Where the number will be: number of
            experimental conditions * number of subfolders.
        names: array
            array of textfields for each experimental condition defined by
            user. User will enter the name of each experimental condition.
        Returns (global in config.py)
        -------
        experiments: 2D list of strings
             list of all the subfolders (experiments) that are in each
             experimental condition.
        allimageNames: list of strings
            name of each experimental condition.
        allUniqueimageNames: list of strings
            same as allimageNames but is a list of unique names
            any duplicates have been removed.
        statsfolderPath: string
            file path to the create statsfolder.
        RscriptPath: string
            file path for R script (defined by user).
           Raises
        ------
        If any subfolders (experiments) have been assigned more than one
        experimental condition. If any of the experimental conditions have
        been assigened the same name.
        """

        config.experiments = [[0 for x in range(0)] for y in range(repeats)]
        same = 0

        # creates a 2D list (experiments) for each experimental conditions
        # with the names of each subfolder. as defined by the user (check).
        for i in range(0, repeats):
            for x in range(0, len(config.subfoldernames)):
                if check[i][x].isSelected() is True:
                    config.experiments[i].append(config.subfoldernames[x])
                    same = same + 1

        # if same is more than the number of subfolders then some of
        # the subfolders must have been assigned to multiple
        # experimental conditions.
        if same > len(config.subfoldernames):
            IJ.showMessage("Error: at least one experiment has been assigned two conditions")
        elif same != len(config.subfoldernames):
            IJ.showMessage("Error: all experiments must be assigned a condition")
        else:
            allimageNames = []

            # get the name of each experimental condition from names (allimageNames)
            for s in range(0, len(config.names)):
                allimageNames.append(config.names[s].getText())
            # remove any duplicate names
            allUniqueimageNames = set(allimageNames)

            # if allimageNames and allUniqueimageNames are not the same then some of the experimental
            # conditions must have the same name.
            if len(allUniqueimageNames) != len(allimageNames):
                IJ.showMessage("Error: all experimental conditions must have a unique name")
            # save file path to statistical analysis folder

            if "Windows" in OS:
                statsfolder = imagefolder+"\\statistical analysis"
                config.statsfolderpPath = imagefolder+"statistical analysis"
                config.statsfolderPath = statsfolderPath.replace("\\", "/")
            elif "Mac" in OS:
                statsfolder = imagefolder+"/statistical analysis"
                config.statsfolderPath = statsfolder

            # create folder called statistical analysis (if it does not
            # already exist).
            if not os.path.exists(statsfolder):
                    os.makedirs(statsfolder)

            # create subfolders for each experimental conditions,
            # within the statistical analysis folder.
            for f in range(0, len(config.names)):
                if config.names[f].getText() == "":
                    IJ.showMessage("Error: each experimental condition must be given a name")
                    break
                if "Windows" in OS:
                    if not os.path.exists(statsfolder+"\\"+config.names[f].getText()):
                        os.makedirs(statsfolder+"\\"+config.names[f].getText())
                elif "Mac" in OS:
                    if not os.path.exists(statsfolder+"/"+config.names[f].getText()):
                        os.makedirs(statsfolder+"/"+config.names[f].getText())

            self.dispose()
            config.RscriptPath = self.Rloc.getText()
            if "Windows" in OS:
                config.RscriptPath = '"'+RscriptPath +'"'
            if config.newuser is True:
                Dialog2()
            else:
                analysed()


class Dialog2(JFrame):
    def __init__(self):
        super(Dialog2, self).__init__()
        self.initUI()
        """ Enter new user name
        Attributes
        ----------
         Title: textfield
            to enter new user name
        OKbutton : button
           continues to next dialog box for entering user name settings.
           Raises
        ------
        If the user name already exists
        """

    def initUI(self):
        panel = JPanel()
        self.getContentPane().add(panel)
        panel.setBackground(Color.WHITE)
        panel.setLayout(None)
        self.setTitle("Enter user name")
        self.setSize(280, 150)
        
        self.Title = JTextArea("")
        self.Title.setBounds(15, 10, 250, 20)
        border = BorderFactory.createEtchedBorder()
        self.Title.setBorder(border)
        panel.add(self.Title)
        
        OKbutton = JButton("Enter", actionPerformed=self.onEnter)
        OKbutton.setBackground(Color.BLACK)
        OKbutton.setBounds(80, 50, 100, 30)
        panel.add(OKbutton)
        
        self.setLocationRelativeTo(None)
        self.setLocation(int(IJ.getScreenSize().width * 0.01),
                          int(IJ.getScreenSize().height * 3/10))
        self.setVisible(True)

    def onEnter(self, e):
        """ Enter new user name
        Parameters
        ----------
         Title: textfield
             to enter user name
        user : string
           new user name
        username: list of string
            list of all predefined user names
           Raises
        ------
        If the user name entered already exists. User is given the choice
        of either entering a new user name or overwriting the existing username.
        """
        config.user = self.Title.getText()
        config.user = config.user+".csv"
        
        # if the new user name entered already exists the user is asked
        # whether to overwrite the existing user name or choose a different
        # name
        if config.user in username:
            frame = JFrame()
            usernameExists = YesNoCancelDialog(frame, "Multi usernames defined", "Multiple usernames defined. Would you like to overwrite the settings?")
            if usernameExists.yesPressed() is True:
                Dialog3()
                self.dispose()               
        else:
             self.dispose()
             Dialog3()
# user selects which channels represent green and red in the images


class Dialog3(JFrame):
    def __init__(self):
        super(Dialog3, self).__init__()
        self.initUI()
        """ select channels for myelin and neurites
           Defines channels for myelin and neurites. Opens the first image in
           selected folder and displays the myelin channel only. For use by
           Dialod4.
           Raises
        ------
        Selected channels for neurites and myelin cannot be the same.
        """
    def initUI(self):
        panel = JPanel()
        self.getContentPane().add(panel)
        panel.setBackground(Color.WHITE)
        panel.setLayout(None)
        self.setTitle("Image settings")
        self.setSize(300, 200)

        Title2 = JTextArea("Neurite:")
        Title2.setBounds(15, 5, 70, 20)
        Title2.setEditable(False)
        background = Color.GRAY
        Title2.setBackground(background)
        panel.add(Title2)

        self.neurite = JComboBox(config.neurites)
        self.neurite.setBounds(15, 30, 250, 20)
        panel.add(self.neurite)

        Title3 = JTextArea("Myelin:")
        Title3.setBounds(15, 60, 80, 20)
        Title3.setEditable(False)
        background = Color.GRAY
        Title3.setBackground(background)
        panel.add(Title3)

        self.myelin = JComboBox(config.neurites)
        self.myelin.setBounds(15, 90, 250, 20)
        panel.add(self.myelin)

        OKbutton = JButton("Enter", actionPerformed=self.onEnter)
        OKbutton.setBackground(Color.BLACK)
        OKbutton.setBounds(80, 120, 100, 30)
        panel.add(OKbutton)

        self.setLocationRelativeTo(None)
        self.setLocation(int(IJ.getScreenSize().width * 0.01),
                             int(IJ.getScreenSize().height * 3/10))
        self.setVisible(True)

    def onEnter(self, e):
            global r, g
            r = self.neurite.getSelectedItem()
            r = int(config.neurites.index(r))
            g = self.myelin.getSelectedItem()
            g = int(config.neurites.index(g))
            if r == g:
                IJ.showMessage("Error: colours for myelin and neurites needs to be different")
            else:
                    self.dispose()
                    getimage()
                    green.show()
                    Dialog4()


class Dialog4(JFrame):
    def __init__(self):
        super(Dialog4, self).__init__()
        self.initUI()
    """ Graphical interface for defining myelin channel    image analysis
    settings
    Interactive interface where the user can change image analysis settings
    and visualise how these settings effect the image. The user can press
    next image at any time, which will open the next image in the selected
    folder and apply any analysis settings that have already been defined.
    Alternatively, the user can open a specific image and check "use
    selected image". The settings can be "Reset" back to defaul at any time.
    
    bug
    ---
    When image processing takes place a wait cursor is displayed. There is
    a lag between when the analysis or function has finished and when the
    image(s) is displayed. This means the wait cursor stops before the
    analysis has finished. This is especially pronounced when several
    images are being displayed.
    """
    def initUI(self):
        panel = JPanel()
        self.getContentPane().add(panel)
        panel.setBackground(Color.WHITE)
        panel.setLayout(None)

        Title3 = JTextArea("Background")
        Title3.setBounds(30, 0, 320, 20)
        Title3.setEditable(False)
        background = Color.GRAY
        Title3.setBackground(background);
        panel.add(Title3)

        backgroundlabel = JLabel()
        backgroundlabel.setBounds(30, 0, 320, 160)
        border = BorderFactory.createEtchedBorder()
        backgroundlabel.setBorder(border)
        panel.add(backgroundlabel)

        self.bCLAHE = JCheckBox("CLAHE?", False, actionPerformed=self.onCLAHE)
        self.bCLAHE.setFocusable(False)
        self.bCLAHE.setBounds(60, 20, 250, 20)
        panel.add(self.bCLAHE)

        self.rollingballcb = JCheckBox("Subtract background (rolling ball)?", 
                                       False, actionPerformed=self.onSubtractrollingball)
        self.rollingballcb.setFocusable(False)
        self.rollingballcb.setBounds(60, 40, 250, 20)
        panel.add(self.rollingballcb)

        self.neurireSubcb = JCheckBox("Subtract background (using neurites)?", 
                                      False, actionPerformed=self.onSubtractneurites)
        self.neurireSubcb.setFocusable(False)
        self.neurireSubcb.setBounds(60, 60, 250, 20)
        panel.add(self.neurireSubcb)

        Title4 = JTextArea("Remove back ground pixels")
        Title4.setBounds(30, 80, 320, 20)
        Title4.setEditable(False)
        background = Color.GRAY
        Title4.setBackground(background)
        panel.add(Title4)

        subpixels = JTextField("Subtract pixels using 'Math' - press enter")
        subpixels.setBounds(90, 110, 210, 20)
        border = BorderFactory.createLineBorder(Color.black)
        subpixels.setBorder(border)
        background = Color.GRAY
        subpixels.setBackground(background)
        panel.add(subpixels)

        self.setpixels2 = JTextField("0", actionPerformed=self.onSubtractpixels)
        self.setpixels2.setBounds(90, 130, 210, 20)
        border = BorderFactory.createLineBorder(Color.black)
        self.setpixels2.setBorder(border)
        panel.add(self.setpixels2)

        self.cellbodycb = JCheckBox("Remove cell bodies?", 
                                    False, actionPerformed=self.onRemovecellbodies)
        self.cellbodycb.setFocusable(False)
        self.cellbodycb.setBounds(120, 160, 170, 20)
        panel.add(self.cellbodycb)

        Title = JTextArea("Threshold method::")
        Title.setBounds(30, 180, 320, 20)
        Title.setEditable(False)
        background = Color.GRAY
        Title.setBackground(background)
        panel.add(Title)

        threshlabel = JLabel()
        threshlabel.setBounds(30, 180, 321, 140)
        border = BorderFactory.createEtchedBorder()
        threshlabel.setBorder(border)
        panel.add(threshlabel)

        autoThresh = ("Default", "Huang", "Intermodes", "IsoData", "Li", "MaxEntropy", 
                     "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", 
                     "RenyiEntropy", "Shanbhag", "Triangle", "Yen")
        self.autoT = JComboBox(autoThresh)
        self.autoT.setBounds(125, 205, 120, 20)
        self.autoT.setEnabled(False)
        panel.add(self.autoT)

        Min1 = JTextArea("Min:")
        Min1.setBounds(90, 220, 40, 20)
        Min1.setEditable(False)
        panel.add(Min1)

        self.tMin = JTextArea("200")
        self.tMin.setBounds(90, 240, 40, 20)
        border = BorderFactory.createLineBorder(Color.black)
        self.tMin.setBorder(border)
        self.tMin.setEditable(False)
        panel.add(self.tMin)

        Max1 = JTextArea("Max:")
        Max1.setBounds(240, 220, 40, 20)
        Max1.setEditable(False)
        panel.add(Max1)

        self.tMax = JTextArea("255")
        self.tMax.setBounds(240, 240, 40, 20)
        border = BorderFactory.createLineBorder(Color.black)
        self.tMax.setBorder(border)
        self.tMax.setEditable(False)
        panel.add(self.tMax)

        setThresholdbutton = JButton("Apply", actionPerformed=self.onSetThreshold)
        setThresholdbutton.setBackground(Color.BLACK)
        setThresholdbutton.setBounds(140, 280, 80, 30)
        panel.add(setThresholdbutton)

        threshlabel2 = JLabel()
        threshlabel2.setBounds(30, 330, 321, 80)
        border = BorderFactory.createEtchedBorder()
        threshlabel2.setBorder(border)
        panel.add(threshlabel2)

        Title3 = JTextField("Remove cell bodies")
        Title3.setBounds(120, 330, 150, 20)
        Title3.setEditable(False)
        background = Color.GRAY
        Title3.setBackground(background)
        panel.add(Title3)

        setCells2 = JTextField("Cell body radius - press enter")
        setCells2.setBounds(100, 355, 190, 20)
        border = BorderFactory.createLineBorder(Color.black)
        setCells2.setBorder(border)
        background = Color.GRAY
        setCells2.setBackground(background)
        setCells2.setEditable(False)
        panel.add(setCells2)

        self.setOutliers = JTextField("0", actionPerformed=self.onOutlier)
        self.setOutliers.setBounds(100, 375, 190, 20)
        border = BorderFactory.createLineBorder(Color.black)
        self.setOutliers.setBorder(border)
        self.setOutliers.setEditable(False)
        panel.add(self.setOutliers)

        threshlabel3 = JLabel()
        threshlabel3.setBounds(30, 440, 321, 65)
        border = BorderFactory.createEtchedBorder()
        threshlabel3.setBorder(border)
        panel.add(threshlabel3)

        self.frangicb = JCheckBox("Apply frangi analysis?", 
                                  False, actionPerformed=self.onFrangi)
        self.frangicb.setFocusable(False)
        self.frangicb.setBounds(120, 415, 170, 20)
        panel.add(self.frangicb)

        setCells3 = JTextField("Grey scale filter - press enter")
        setCells3.setBounds(97, 450, 200, 20)
        border = BorderFactory.createLineBorder(Color.black)
        setCells3.setBorder(border)
        background = Color.GRAY
        setCells3.setBackground(background)
        setCells3.setEditable(False)
        panel.add(setCells3)

        self.greyscale = JTextField("0", actionPerformed=self.onGreyscale)
        self.greyscale.setBounds(97, 470, 200, 20)
        border = BorderFactory.createLineBorder(Color.black)
        self.greyscale.setBorder(border)
        self.greyscale.setEditable(True)
        panel.add(self.greyscale)

        userimage = JCheckBox("Use selected image", False, actionPerformed=self.onUserimage)
        userimage.setFocusable(False)
        userimage.setBounds(180, 530, 170, 20)
        panel.add(userimage)

        Nextimagebutton = JButton("Next Image", actionPerformed=self.onNextimage)
        Nextimagebutton.setBackground(Color.BLACK)
        Nextimagebutton.setBounds(240, 560, 100, 30)
        panel.add(Nextimagebutton)

        Cancelbutton = JButton("Cancel", actionPerformed=self.onCancel)
        Cancelbutton.setBackground(Color.BLACK)
        Cancelbutton.setBounds(140, 560, 100, 30)
        panel.add(Cancelbutton)

        Nextdialogbutton = JButton("Next", actionPerformed=self.onNextdialog)
        Nextdialogbutton.setBounds(40, 560, 100, 30)
        Nextdialogbutton.setBackground(Color.BLACK)
        panel.add(Nextdialogbutton)

        Resetbutton = JButton("Reset", actionPerformed=self.onReset)
        Resetbutton.setBounds(70, 530, 100, 30)
        Resetbutton.setBackground(Color.BLACK)
        panel.add(Resetbutton)

        self.setTitle("Settings for myelin")
        self.setSize(395, 640)
        self.setLocationRelativeTo(None)
        self.setLocation(int(IJ.getScreenSize().width * 0.01),
                             int(IJ.getScreenSize().height * 3/10))
        self.setVisible(True)
        global c, d, f, m, myelinpanel
        c = 0
        d = 0
        f = False
        m = False
        myelinpanel = self

    def onCLAHE(self, cb1):
        """
        Applies CLAHE to image if checked. If unchecked closes image and
        reverts back to the original unaltered image (if other background
        settings have not been applied). For the below background
        subtraction methods (onSubtractrollingball() and
        onSubtractneurites()) CLAHE will always be applied to the image
        first. Consequently, if one of the bacground methods is applied
        before checking CLAHE, the image will be closed and a new image
        will be opened where CLAHE will be applied to the image and then
        the background subtracted. If CLAHE is unchecked and background
        subtraction has already been set, the image will be closed and
        the same image with only background subtraction will be displayed.
        Parameters
        ----------
        backgroundsubRolling : bool
           state of rolling ball background subtraction checkbox
           (rollingballcb)
        backgroundsubNeurite : bool
           state of neurite background subtraction checkbox (neurireSubcb)
        Returns (global in config.py)
        -------
        mCLAHE: bool
             state of CLAHE checkbox (bCLAHE)
        """
        self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
        mCLAHE = cb1.getSource()
        config.mCLAHE = mCLAHE.isSelected()
        if config.mCLAHE is True:
            # remove unwanted images so that the number of images open
            # does not pile up.
            basicfunctions.toOriginal()
            basicfunctions.green2()
            basicfunctions.CLAHE()
            if config.backgroundsubRolling is True:
                basicfunctions.rollingsubtract()
            elif config.backgroundsubNeurite is True:
                neuritesubtract()
        else:
            basicfunctions.closeimage()
            if (config.backgroundsubRolling is True) or (config.backgroundsubNeurite is True):
                basicfunctions.green2()
            if config.backgroundsubRolling is True:
                basicfunctions.rollingsubtract()
            elif config.backgroundsubNeurite is True:
                neuritesubtract()
        self.setCursor(Cursor.getDefaultCursor())

    def onSubtractrollingball(self, cb1):
        """ Rolling ball background subtraction
        Refer to onCLAHE for full description of background subtraction.
        Parameters
        ----------
        backgroundsubNeurite : bool
           state of neurite background subtraction checkbox (neurireSubcb)
        mCLAHE: bool
             state of CLAHE checkbox (bCLAHE)
        Returns (global in config.py)
        -------
        backgroundsubRolling: bool
             state of rolling ball background subtraction checkbox
             (rollingballcb)
        Raises
        ------
        Only one background subtraction method can be be selected at any
        one time (rolling ball or neurite background subtraction). The
        user must first deselect one before the other can be selected.
        """
        if config.backgroundsubNeurite is False:
            self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
            bgcbstate = cb1.getSource()
            config.backgroundsubRolling = bgcbstate.isSelected()
            if config.backgroundsubRolling is True:
                basicfunctions.ifOriginal()
                basicfunctions.rollingsubtract()
            elif config.backgroundsubRolling is False:
                 basicfunctions.closeimage()
                 if config.mCLAHE is True:
                     basicfunctions.green2()
                     basicfunctions.CLAHE()
            self.setCursor(Cursor.getDefaultCursor())
        else:
            IJ.showMessage("Error: only one background subtraction can be selected")
            self.bgcb.setSelected(False)

    def onSubtractneurites(self, cb1):
        """ neurite background subtraction
        Refer to onCLAHE for full description of background subtraction.
        Parameters
        ----------
         backgroundsubRolling : bool
           state of neurite background subtraction checkbox (neurireSubcb)
        mCLAHE: bool
             state of CLAHE checkbox (bCLAHE)
        Returns (global in config.py)
        -------
        backgroundsubNeurite: bool
             state of rolling ball background subtraction checkbox
             (rollingballcb)
        Raises
        ------
        same as onSubtractrollingball()
        """
        if config.backgroundsubRolling is False:
            self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
            bgcbstate2 = cb1.getSource()
            config.backgroundsubNeurite = bgcbstate2.isSelected()
            if config.backgroundsubNeurite is True:
                basicfunctions.ifOriginal()
                neuritesubtract()
            elif config.backgroundsubRolling is False:
                 basicfunctions.closeimage()
                 if config.mCLAHE is True:
                     basicfunctions.green2()
                     basicfunctions.CLAHE()
            self.setCursor(Cursor.getDefaultCursor())
        else:
            IJ.showMessage("Error: only one background subtraction can be selected")
            self.bgcb.setSelected(False)

    def onSubtractpixels(self, e):
         """ Math: pixel subtraction
         if performed for the first time the current image will be
         duplicated, displayed and pixel subtraction performed. If this is
         repeated the current image will  first be closed (stops number of
         displayed images mounting up).
         Parameters
         ----------
         m : bool
           marker for whether analysis has already been performed
         Returns (global in config.py)
         -------
         setpixles: string
             pixel subtraction value (has to be a string to run).
         """
         global m
         config.setpixels = self.setpixels2.getText()
         if m is True:
            basicfunctions.closeimage()
         basicfunctions.subpixels(config.setpixels)
         m = True

    def onRemovecellbodies(self, cb3):
        """ cell body selection
         closes any images from defining background subtraction to only the
         unprocessed original myelin channel image is displayed. Disables
         checkboxes for background subtraction enables checkboxes etc for
         selecting cell bodies. If deselected the opposite will be
         performed. Any previous settings for background subtraction will
         be performed.
        """
        config.cellbodycb = self.cellbodycb.isSelected()
        if self.cellbodycb.isSelected() is True:
            self.setOutliers.setEditable(True)
            self.setpixels2.setEditable(False)
            self.rollingballcb.setEnabled(False)
            self.neurireSubcb.setEnabled(False)
            self.bCLAHE.setEnabled(False)
            self.tMin.setEditable(True)
            self.tMax.setEditable(True)
            self.autoT.setEnabled(True)
            basicfunctions.closeimagebg()
            if config.setpixels != "0":
                basicfunctions.closeimage()
            getimage()
            green.show()
            
        else:
            self.setOutliers.setEditable(False)
            self.rollingballcb.setEnabled(True)
            self.neurireSubcb.setEnabled(True)
            self.bCLAHE.setEnabled(True)
            self.setpixels2.setEditable(True)
            self.tMin.setEditable(False)
            self.tMax.setEditable(False)
            basicfunctions.closeimage()
            getimage()
            green.show()
            applybackground()
            self.autoT.setEnabled(False)

    def onSetThreshold(self, e):
        """ Applies threshold
         Applies thresholding method selected by the user and crops the
         histogram according to user defined settings. First time this
         button is pressed the threshold will not be performed.This is
         because the user is meant to perform the thresholding/histogram
         cropping themselves using ImageJ - try and find the correct
         settings. The second the button is pressed thresholding will be
         performed - so the user can still change settings.
         Parameters
         ----------
         cellbodycb : bool
           checkbox to be selected before cell bodies are selected.
         Returns (global in config.py)
         -------
         c: integer
             counter for how many times button has been pressed by the user.
         Min: string
             value for cropping histogram
         Max: strong
             value for cropping histogram
         threshChoise: strong
             thresholding method selected by user
         Raises
         ------
         cellbodycb must be selected.
        """
        
        global d,c
        self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
        if config.cellbodycb is True:
            config.Min = int(self.tMin.getText())
            config.Max = int(self.tMax.getText())
            config.threshChoice = self.autoT.getSelectedItem()
            # perform threshold is button has been pressed more
            # than once.
            if c != 0:
                # if button has been prssed more than once close
                # current image (stops lots of images piling up).
                if c > 0:
                    basicfunctions.closeimage()
                getimage()
                ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                                int(IJ.getScreenSize().height * 1/14))
                green.show()
                IJ.setAutoThreshold(green, config.threshChoice+" dark")
                IJ.setRawThreshold(green, config.Min, config.Max, None)
                Prefs.blackBackground = True
                IJ.run(green, "Convert to Mask", "")
                basicfunctions.bgTitle2()
                threshResult = green.duplicate()
                d = 0
        else:
            IJ.showMessage("Error: first select checkbox for: Remove cell bodies?")
        c = c + 1
        self.setCursor(Cursor.getDefaultCursor())

    def onOutlier(self, e):
            """ Remove outliers
            Runs remove outliers to select cellbodies using user defined
            value.
             Parameters
             ----------
             cellbodycb : bool
                  checkbox to be selected before cell bodies are selected.
             Returns (global in config.py)
             -------
             d: integer
                 counter for how many times analysis has been performed.
             radius: string
                 value for removing outliers.
             Raises
             ------
             cellbodycb must be selected.
            """
            global d
            if config.cellbodycb is True:
                    sender = e.getSource()
                    value = sender.getText()
                    config.radius = str(value)
                    self.setOutliers.setText(config.radius)
                    # if remove outliers has been performed previously
                    # current image is closed (to stop images piling up).
                    if d > 0:
                        basicfunctions.closeimage()
                    threshResult = IJ.getImage()
                    threshResult = threshResult.duplicate()
                    ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                                    int(IJ.getScreenSize().height * 1/14))
                    IJ.run(threshResult, "Remove Outliers...", "radius="+config.radius+" threshold=50 which=Bright")
                    threshResult.show()
                    basicfunctions.bgTitle2()
                    d = d + 1
            else:
                    IJ.showMessage("Error: Please check the remove cell bodies checkbox first")

    def onFrangi(self, e):
        """ Frangi vesselness
            Performs frangi vesselness using the function frangifilter().
             Returns (global in config.py)
             -------
             d: integer
                 counter for how many times analysis has been performed.
             radius: string
                 value for removing outliers.
        """
        config.frangicb = self.frangicb.isSelected()
        if config.frangicb is True:
             frangifilter(self)
        else:
            basicfunctions.closeimage()

    def onGreyscale(self, e):
        """ grey scale morphology filter
         Returns (global in config.py)
         -------
         attrival2: string
             grey scale filter value ("Minimum").
         Raises
         ------
         Frangi vesselness must be performed first.
        """
        if config.frangicb is True:
            self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
            sender = e.getSource()
            attrival2 = sender.getText()
            config.greyscaleMinVal = str(attrival2)
            self.greyscale.setText(config.greyscaleMinVal)
            ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                                int(IJ.getScreenSize().height * 1/14))
            
            # if title of image does not contain vesselness then close
            basicfunctions.ifVesselness()
            IJ.run("Gray Scale Attribute Filtering", "operation=Opening attribute=[Box Diagonal] minimum="+config.greyscaleMinVal+" connectivity=4")
            g = IJ.getImage()
            g.setTitle("grey scale")
            self.setCursor(Cursor.getDefaultCursor())
        else:
            IJ.showMessage("Error: Please apply frangi analysis first")

    def onUserimage(self, e):
        """
        Allow user to open an image and select as current image so that it
        can used to optimise analysis settings. This must be a .tif
        composite of myelin and neurires channels. The image will be split
        and only the myelin channel displayed.
        """
        global userimagename
        sender = e.getSource()
        config.userimage2 = True
        userimagename = IJ.getImage()
        if sender.isSelected() is True:
            getimage()
            green.show()
        else:
            config.userimage2 = False

    def onCancel(self, e):
        basicfunctions.closeallimages()
        self.dispose()

    def onNextdialog(self, e):
        config.userimage2 = False
        self.dispose()
        basicfunctions.closeallimages()
        Dialog5()

    def onReset(self, e):
        """
        Reset all settings back to default.
        Open unprocessed myelin channel image.
           """
        global c, d, m
        basicfunctions.closeallimages()
        getimage()
        green.show()
        self.setOutliers.setEditable(False)
        self.cellbodycb.setSelected(False)
        self.rollingballcb.setEnabled(True)
        self.rollingballcb.setSelected(False)
        self.neurireSubcb.setEnabled(True)
        self.neurireSubcb.setSelected(False)
        self.bCLAHE.setEnabled(True)
        self.bCLAHE.setSelected(False)
        self.setpixels2.setEditable(True)
        self.tMin.setEditable(False)
        self.tMax.setEditable(False)
        self.frangicb.setSelected(False)
        self.autoT.setEnabled(True)
        self.setpixels2.setText("0")
        self.tMin.setText("200")
        self.tMax.setText("255")
        self.setOutliers.setText("0")
        self.greyscale.setText("0")
        config.mCLAHE = False
        config.backgroundsubRolling = False
        config.backgroundsubNeurite = False
        config.setpixels = "0"
        config.cellbodycb = False
        config.Min = "0"
        config.Max = "0"
        config.radius = "0"
        config.frangicb = False
        config.greyscaleMinVal = "0"
        c = 0
        d = 0
        m = False

    def onNextimage(self, e):
            """ Next image
            Dislays the next image in folder. Any settings that have been
            defined will be performed.
            """
            self.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR))
            getNext()
            getimage()
            green2 = green.duplicate()
            green.show()
            if (config.cellbodycb is False) and (config.frangicb is False):
                ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                                int(IJ.getScreenSize().height * 1/14))
                green2.show()
                applybackground()
            if (config.cellbodycb is True) and (config.frangicb is False):
                removecellbodies()
            elif config.frangicb is True:
                if config.cellbodycb is False:
                    frangifilter(self)
                if config.cellbodycb is True:
                        removecellbodies()
                        frangifilter(self)
            self.setCursor(Cursor.getDefaultCursor())


class Dialog5(JFrame):
    def __init__(self):
        super(Dialog5, self).__init__()
        self.initUI()
        
    def initUI(self):
        """ Graphical interface for defining neurite channel analysis
        settings
         Provides choices for normalise local contrast (NLC) and despeckle.
         Also has a sparse neurites option which opens another dialog box.
        
        """
        
        panel = JPanel()
        self.getContentPane().add(panel)
        panel.setBackground(Color.WHITE)
        # show neuite channel and run NLC if sparse neurite (SN) settings
        # have not been defined (Dialog6).
        if config.SN is False:
              getimage()
              red.show()
              IJ.run(red, "Normalize Local Contrast", "block_radius_x=40 block_radius_y=40 standard_deviations="+config.contrast+" center stretch")
              IJ.run(red, "Auto Threshold", "method=default white")
              IJ.run(red, "Invert LUT", "")

        Title = JTextArea("Normalise local contrast:")
        Title.setBounds(30, 0, 320, 20)
        Title.setEditable(False)
        background = Color.GRAY
        Title.setBackground(background)
        panel.add(Title)

        border1 = JLabel()
        border1.setBounds(30, 20, 321, 150)
        border = BorderFactory.createEtchedBorder()
        border1.setBorder(border)
        panel.add(border1)

        sd = JTextArea("Standard deviation")
        sd.setBounds(125, 45, 160, 20)
        sd.setEditable(False)
        panel.add(sd)

        tMin = JTextField("0.5", actionPerformed=self.onEnter)
        tMin.setBounds(165, 70, 40, 20)
        border = BorderFactory.createLineBorder(Color.black)
        tMin.setBorder(border)
        if config.SN is True:
           tMin.setEditable(False)
        else:
           tMin.setEditable(True)
        panel.add(tMin)

        self.despecklecb = JCheckBox("Despeckle?", True, actionPerformed=self.onDespeckle)
        self.despecklecb.setFocusable(False)
        self.despecklecb.setBounds(140, 110, 110, 20)
        self.despecklecb.setSelected(False)        
        if config.SN is True:
            self.despecklecb.setEnabled(False)
            
        panel.add(self.despecklecb)
        startbutton = JButton("Sparse neurites", actionPerformed=self.onSN)
        startbutton.setBackground(Color.BLACK)
        startbutton.setBounds(130, 135, 120, 30)
        panel.add(startbutton)
        
        Nextimagebutton = JButton("Next Image", actionPerformed=self.onNextimage)
        Nextimagebutton.setBackground(Color.BLACK)
        Nextimagebutton.setBounds(240, 65, 100, 30)
        panel.add(Nextimagebutton)
        panel.setLayout(None)
        
        self.setTitle("Settings for dense neurites")
        self.setLocationRelativeTo(None)
        self.setLocation(int(IJ.getScreenSize().width * 0.01),
                             int(IJ.getScreenSize().height * 3/10))
        self.setVisible(True)
        startbutton = JButton("Start analysis", actionPerformed=self.onStartanalysis)
        startbutton.setBackground(Color.BLACK)
        startbutton.setBounds(130, 180, 120, 30)
        panel.add(startbutton)
        
        cancelbutton = JButton("Cancel", actionPerformed=self.onCancel)
        cancelbutton.setBackground(Color.BLACK)
        cancelbutton.setBounds(10, 180, 100, 30)
        panel.add(cancelbutton)

        userimage = JCheckBox("Use selected image", False, actionPerformed=self.onUserimage)
        userimage.setFocusable(False)
        userimage.setBounds(260, 185, 170, 20)
        panel.add(userimage)
        
        self.setSize(410, 250)

    def onEnter(self, e):
        """ Normalise local contrast
        Parameters
        ----------
        contrast: string
            value for normalise local contrast
        """
        sender = e.getSource()
        config.contrast = sender.getText()
        basicfunctions.closeimage()
        getimage()
        red.show()
        IJ.run(red, "Normalize Local Contrast", "block_radius_x=40 block_radius_y=40 standard_deviations="+config.contrast+" center stretch")
        IJ.run(red, "Auto Threshold", "method=default white")
        IJ.run(red, "Invert LUT", "")

    def onDespeckle(self, cb1):
        despecklecbstate = cb1.getSource()
        config.despeckle = despecklecbstate.isSelected()
        if config.despeckle is True:
                IJ.run("Despeckle", "")
        else:
                basicfunctions.closeimage()
                getimage()
                IJ.run(red, "Normalize Local Contrast", "block_radius_x=40 block_radius_y=40 standard_deviations="+config.contrast+" center stretch")
                IJ.run(red, "Auto Threshold", "method=default white")
                IJ.run(red, "Invert LUT", "")
                red.show()

    def onSN(self, setb):
        """ Open dialog for sparse neurite settings.
        """

        config.userimage2 = False
        self.dispose()
        basicfunctions.closeallimages()
        getimage()
        red.show()
        Dialog6()

    def onCancel(self, setb):
         self.dispose()
         basicfunctions.closeimage()

    def onNextimage(self, setb):
         """ Get next image in selected folder
          Open the next image in user selected image folder, display the
          red channel and run NLC.
          Parameters
          ----------
          constrast: string
              "standard deviations" for NLC.
        """
         getNext()
         getimage()
         red.show()
         IJ.run(red, "Normalize Local Contrast", "block_radius_x=40 block_radius_y=40 standard_deviations="+config.contrast+" center stretch")
         IJ.run(red, "Invert LUT", "")
         if config.despeckle is True:
             IJ.run(red, "Despeckle", "")

    def onStartanalysis(self, b):
        self.dispose()
        start_time = time.time()
        analysed()

    def onUserimage(self, e):
        """
        Allow user to open an image and select as current image so that it
        can used to optimise analysis settings. This must be a .tif
        composite of myelin and neurires channels. The image will be split
        and only the myelin channel displayed.
        """
        global userimagename
        sender = e.getSource()
        config.userimage2 = True
        userimagename = IJ.getImage()
        if sender.isSelected() is True:
            getimage()
            red.show()
        else:
            config.userimage2 = False

class Dialog6(JFrame):

    def __init__(self):
            super(Dialog6, self).__init__()
            self.initUI()

    def initUI(self):
            """ Graphical interface for defining sparse neurite analysis

            GUI to define background subtraction, enhancement using CLAHE
            and all ImageJ threshold methods. 
  
            """

            panel = JPanel()
            self.getContentPane().add(panel)
            panel.setBackground(Color.WHITE)
            self.setTitle("Settings sparse neurites")
            panel = JPanel()
            self.getContentPane().add(panel)
            
            panel.setBackground(Color.WHITE)
            panel.setLayout(None)
            
            Title3 = JTextArea("Background")
            Title3.setBounds(30, 0, 320, 20)
            Title3.setEditable(False)
            background = Color.GRAY
            Title3.setBackground(background)
            panel.add(Title3)
            
            backgroundlabel = JLabel()
            backgroundlabel.setBounds(30, 0, 320, 60)
            border = BorderFactory.createEtchedBorder()
            backgroundlabel.setBorder(border)
            panel.add(backgroundlabel)

            cellbodycb = JCheckBox("CLAHE?", False, actionPerformed=self.onCLAHE)
            cellbodycb.setFocusable(False)
            cellbodycb.setBounds(110, 20, 80, 20)
            panel.add(cellbodycb)
            
            bgcb = JCheckBox("Subtract background?", False, actionPerformed=self.onbg)
            bgcb.setFocusable(False)
            bgcb.setBounds(110, 40, 180, 20)
            panel.add(bgcb)
            
            Title = JTextArea("Threshold method:")
            Title.setBounds(30, 60, 320, 20)
            Title.setEditable(False)
            background = Color.GRAY
            Title.setBackground(background)
            panel.add(Title)
            
            test = JLabel()
            test.setBounds(30, 60, 321, 170)
            border = BorderFactory.createEtchedBorder()
            test.setBorder(border)
            panel.add(test)

            autoThresh = ("Default", "Huang", "Intermodes", "IsoData", "Li", "MaxEntropy", 
                          "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", 
                          "RenyiEntropy", "Shanbhag", "Triangle", "Yen")
            self.autoT2 = JComboBox(autoThresh)
            self.autoT2.setBounds(125, 90, 120, 20)
            panel.add(self.autoT2)
            
            panel.repaint()
            Min = JTextArea("Min:")
            Min.setBounds(110, 130, 40, 20)
            Min.setEditable(False)
            panel.add(Min)
            
            self.tMin = JTextArea("0")
            self.tMin.setBounds(110, 150, 40, 20)
            border = BorderFactory.createLineBorder(Color.black)
            self.tMin.setBorder(border)
            self.tMin.setEditable(True)
            panel.add(self.tMin)
            
            Max = JTextArea("Max:")
            Max.setBounds(210, 130, 40, 20)
            Max.setEditable(False)
            panel.add(Max)
            
            self.tMax = JTextArea("255")
            self.tMax.setBounds(210, 150, 40, 20)
            border = BorderFactory.createLineBorder(Color.black)
            self.tMax.setBorder(border)
            self.tMax.setEditable(True)
            panel.add(self.tMax)
            
            setbutton = JButton("Apply", actionPerformed=self.onSet)
            setbutton.setBackground(Color.BLACK)
            setbutton.setBounds(140, 190, 80, 30)
            panel.add(setbutton)

            userimage = JCheckBox("Use selected image", False, actionPerformed=self.onUserimage)
            userimage.setFocusable(False)
            userimage.setBounds(100, 240, 170, 20)
            panel.add(userimage)

            cancelbutton = JButton("Cancel", actionPerformed=self.onCancel)
            cancelbutton.setBackground(Color.BLACK)
            cancelbutton.setBounds(30, 270, 80, 30)
            panel.add(cancelbutton)
            
            cancelbutton = JButton("Next image", actionPerformed=self.onNext)
            cancelbutton.setBackground(Color.BLACK)
            cancelbutton.setBounds(250, 270, 100, 30)
            panel.add(cancelbutton)
            
            cancelbutton = JButton("Set for analysis", actionPerformed=self.onSetanalysis)
            cancelbutton.setBackground(Color.BLACK)
            cancelbutton.setBounds(115, 270, 130, 30)
            panel.add(cancelbutton)
            
            global c
            c = 0
            self.setSize(400, 330)
            self.setLocationRelativeTo(None)
            self.setLocation(int(IJ.getScreenSize().width * 0.01),
                                 int(IJ.getScreenSize().height * 2/10))
            self.setVisible(True)

    def onSetanalysis(self, cb1):
            if (config.Min2 == 0) and (config.Max2 == 0):
               IJ.showMessage("Error: please select threshold settings")            
            elif (config.Min2 == 0) and (config.Max2 == 255):
               IJ.showMessage("Error: please select threshold settings")
            else:
               config.userimage2 = False
               config.SN = True
               self.setVisible(False)
               Dialog5()

    def onCLAHE(self, cb1):
            """ CLAHE

            Run CLAHE if checked. If unchecked close image and run 
            background subtraction if selected.
            
            Parameters
            ----------
            cellbodycb : bool
                  checkbox to be selected before cell bodies are selected.
            Returns (global in config.py)
            -------
            mCLAHE2: bool
                 Perform CLAHE?
                     
            """
            
            CLAHEcbstate = cb1.getSource()
            config.mCLAHE2 = CLAHEcbstate.isSelected()
            if config.mCLAHE2 is True:
                IJ.run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)")
            else:
                basicfunctions.closeimage()
                getimage()
                red.show()
                if config.Sbgcbstate is True:
                    IJ.run("Subtract Background...", "rolling=50")

    def onbg(self, cb1):
            """ rolling ball background subtraction. 

           
                     
            """
            
            bgcbstate2 = cb1.getSource()
            config.Sbgcbstate = bgcbstate2.isSelected()
            if config.Sbgcbstate is True:
                if w.getImageCount() == 0:
                    getimage()
                    red.show()
                    IJ.run("Subtract Background...", "rolling=50")
                else:
                    IJ.run("Subtract Background...", "rolling=50")
            else:
                basicfunctions.closeimage()
                getimage()
                red.show()
                if config.mCLAHE2 is True:
                    IJ.run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)")

    def onCancel(self, e):
        """ Cancel and close dialog box. 
                   
        """

        config.userimage2 = False
        basicfunctions.closeallimages()
        self.dispose()
        config.SN = False
        Dialog5()

    def onSet(self, e):
            current = IJ.getImage()
            global c
            config.Min2 = int(self.tMin.getText())
            config.Max2 = int(self.tMax.getText())
            config.threshChoice2 = self.autoT2.getSelectedItem()
            if c != 0:
                basicfunctions.closeimage()
                getimage()
                red.show()
                if config.mCLAHE2 is True:
                    IJ.run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)")
                if config.Sbgcbstate is True:
                    IJ.run("Subtract Background...", "rolling=50")
                IJ.setAutoThreshold(red, config.threshChoice2)
                IJ.run(red, "Invert LUT", "")
                IJ.setRawThreshold(red, config.Min2, config.Max2, None)
                IJ.run(red, "Convert to Mask", "")
            c = c + 1

    def onNext(self, cb1):
        basicfunctions.closeimage()
        getNext()
        getimage()
        red.show()
        if config.mCLAHE2 is True:
                IJ.run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate)")
        if config.Sbgcbstate is True:
            IJ.run("Subtract Background...", "rolling=50")
        if (config.Min2 != "0") or (config.Max2 != "0"):
                IJ.setAutoThreshold(red, config.threshChoice2)
                IJ.setRawThreshold(red, config.Min2, config.Max2, None)
                IJ.run("Convert to Mask", "")

    def onUserimage(self, e):
        """
        Allow user to open an image and select as current image so that it
        can used to optimise analysis settings. This must be a .tif
        composite of myelin and neurires channels. The image will be split
        and only the myelin channel displayed.
        """
        global userimagename
        sender = e.getSource()
        config.userimage2 = True
        userimagename = IJ.getImage()
        if sender.isSelected() is True:
            getimage()
            red.show()
        else:
            config.userimage2 = False
Dialog1()
