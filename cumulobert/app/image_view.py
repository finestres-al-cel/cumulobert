"""Define ImageView widget as an extension of pg.PlotWidget"""
from astropy.visualization import ImageNormalize, ZScaleInterval, SinhStretch
import numpy as np

from PyQt6.QtCore import Qt, QPointF, QRectF
from PyQt6.QtGui import QTransform, QFont
import pyqtgraph as pg


class ImageView(pg.PlotWidget):
    """ Manage image plotting

    Methods
    -------
    (see pg.PlotWidget)
    __init__

    Attributes
    ----------
    (see pg.PlotWidget)
    """
    def __init__(self, image):
        """Initialize instance

        Arguments
        ---------
        image: Image
        The Image to be shown
        """
        # initialize plotting
        super().__init__()
        self.show()

        # tranform to RA/Dec
        # TODO: fix the transformation from pixels to radius in cluster_image.py
        # then change this to True
        self.pix_to_radec = False

        # keep image
        self.image = image
        self.imageData = image.data
        self.imageShape = image.data.shape
        self.imageRa = image.ra
        self.imageDec = image.dec

        self.aperRadiusList = None
        self.annularRadiusList = None

        # plot image
        self.imageItem = None
        self.colorBar = None
        self.colorBarValues = None
        if self.pix_to_radec:
            self.invertXFlag = True
            self.invertYFlag = True
        else:
            self.invertXFlag = False
            self.invertYFlag = False
        self.updatePlot()

    def setInvertXFlag(self, status):
        """Set the invert flag for X axis

        Arguments
        ---------
        status: bool
        Set the invertXFlag status.
        """
        self.invertXFlag = status

    def setInvertYFlag(self, status):
        """Set the invert flag for Y axis

        Arguments
        ---------
        status: bool
        Set the invertYFlag status.
        """
        self.invertYFlag = status

    def resetPlot(self):
        """Reset plot"""
        # reset labels
        self.setLabel(axis='left', text='')
        self.setLabel(axis='bottom', text='')
        # Remove existing items if they exists
        if self.colorBar is not None:
            self.getPlotItem().layout.removeItem(self.colorBar)
            self.colorBar = None
        if self.imageItem is not None:
            self.getPlotItem().removeItem(self.imageItem)
            self.imageItem = None
        # reset plot
        self.clear()

    def setPlot(self):
        """Load plot settings"""

        # add label
        self.showAxes(True)
        if self.pix_to_radec:
            self.setLabel(axis='left', text='Dec')
            self.setLabel(axis='bottom', text='RA')
        else:
            self.setLabel(axis='left', text='pix')
            self.setLabel(axis='bottom', text='pix')

        # Set To Larger Font
        leftAxis = self.getAxis('left')
        bottomAxis = self.getAxis('bottom')
        font = QFont("Helvetica", 18)
        leftAxis.label.setFont(font)
        bottomAxis.label.setFont(font)
        leftAxis.setTickFont(font)
        bottomAxis.setTickFont(font)

        # invert axis
        if self.invertXFlag:
            self.invertX(True)
        else:
            self.invertX(False)
        if self.invertYFlag:
            self.invertY(True)
        else:
            self.invertY(False)

    def updatePlot(self):
        """Update plot"""
        # keep colorbar values before reset
        if self.colorBar is not None:
            self.colorBarValues = self.colorBar.values

        # reset plot
        self.resetPlot()

        # plot image
        self.imageItem = pg.ImageItem(self.imageData.transpose())
        self.addItem(self.imageItem)

        # set axis range
        if self.pix_to_radec:
            # QRectF(aleft: float, atop: float, awidth: float, aheight: float)
            rect = QRectF(
                self.imageRa[0],
                self.imageDec[0],
                self.imageRa[-1] - self.imageRa[0],
                self.imageDec[-1] - self.imageDec[0]
            )
            self.imageItem.setRect(rect)

        colorMap = pg.colormap.get("CET-L2")  # choose perceptually uniform, diverging color map

        # generate an adjustabled color bar
        if self.colorBarValues is None:
            self.colorBar = pg.ColorBarItem(
                # we multiply the minimum and maximum levels by 1.08 and 0.017
                # these are arbitrary values that can be changed in the display
                # but have been found to have a good balance in some images
                values=(np.min(self.imageData)*1.055, np.max(self.imageData)*0.017),
                colorMap=colorMap)
        else:
            self.colorBar = pg.ColorBarItem(
                # we use previous limits
                values=self.colorBarValues,
                colorMap=colorMap)
        # link color bar and color map to correlogram, and show it in plotItem:
        self.colorBar.setImageItem(self.imageItem, insert_in=self.getPlotItem())

        # plot aperture radius
        if self.aperRadiusList is not None:
            apertureRadiusItem = pg.ScatterPlotItem(
                size=1, movable=False, pen=pg.mkPen("r"))
            if self.pix_to_radec:
                apertureRadiusItem.addPoints([
                    {"pos": (item.ra.value, item.dec.value), 'data': 1}
                    for item in self.aperRadiusList
                ])
            else:
                apertureRadiusItem.addPoints([
                    {"pos": (x, y), 'data': 1}
                    for x, y in self.aperRadiusList
                ])
            self.addItem(apertureRadiusItem)

        # plot annular radius
        if self.annularRadiusList is not None:
            annularRadiusItem = pg.ScatterPlotItem(
                size=1, movable=False, pen=pg.mkPen("g"))
            if self.pix_to_radec:
                annularRadiusItem.addPoints([
                    {"pos": (item.ra.value, item.dec.value), 'data': 1}
                    for item in self.annularRadiusList
                ])
            else:
                annularRadiusItem.addPoints([
                    {"pos": (x, y), 'data': 1}
                    for x, y in self.annularRadiusList
                ])
            self.addItem(annularRadiusItem)

        # load plot settings
        self.setPlot()
