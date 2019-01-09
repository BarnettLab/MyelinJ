
from ij import IJ, WindowManager
from ij.gui import ImageWindow
w = WindowManager

def closeimage():
        """
        Close active window
        """

        images = w.getImageCount()

        if images != 0:
                IJ.run("Close")


def closeallimages():
        """
        Close all active windows
        """

        images = w.getImageCount()
        if images != 0:
                IJ.run("Close All")


def closeimagebg():
        """
        If more than one image is open, close.
        """
        images = w.getImageCount()
        if images > 1:
                IJ.run("Close")


def green2():
    """
      Displays a duplicate of the active image.
    """
    green2 = IJ.getImage()
    green2 = green2.duplicate()
    ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                    int(IJ.getScreenSize().height * 1/14))
    green2.show()


def CLAHE():
        IJ.run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None* fast_(less_accurate) process_as_composite")
        bgTitle()


def rollingsubtract():
        IJ.run("Subtract Background...", "rolling=50")
        bgTitle()


def subpixels(pixel):
     """ Math: subtraction
         duplicates and displayes current image.
        Performs puxel subtraction using Math.
     """
     green = IJ.getImage()
     green = green.duplicate()
     IJ.run(green, "Subtract...", "value="+pixel)
     ImageWindow.setNextLocation(int(IJ.getScreenSize().width * 1/3),
                                     int(IJ.getScreenSize().height * 1/14))
     green.show()
     bgTitle()


def toOriginal():
    """ Go back to unprocessed image
    The original unprocessed image displayed is called
    "original". If the active image is not called
    "original" then it is closed. This stops the
    number of images being displayed building up.
        """
    green = IJ.getImage()
    title = green.getTitle()
    if title not in "original":
        closeimage()


def ifOriginal():
    green = IJ.getImage()
    title = green.getTitle()
    if title in "original":
        green2()


def ifVesselness():
    """
        If image title does not contain "vesselness" then close image.
    """
    green = IJ.getImage()
    title = green.getTitle()
    if "vesselness" not in title:
        closeimage()


def bgTitle():
    """
        sets current image title
    """
    g = IJ.getImage()
    g.setTitle("original + background subtraction")


def bgTitle2():
    """
        sets current image title
    """
    g = IJ.getImage()
    g.setTitle("original + cell body selection")
