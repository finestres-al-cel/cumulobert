"""cumulobert main window"""
import os

from PyQt6.QtCore import QSize, Qt, pyqtSlot
from PyQt6.QtWidgets import (
    QCheckBox,
    QFileDialog,
    QLabel,
    QMainWindow,
    QPushButton,
    QStatusBar,
    QTableView,
    QToolBar,
    QMessageBox,
)

from cumulobert.app.environment import (
    HEIGHT, ICON_SIZE, WIDTH
)
from cumulobert.app.error_dialog import ErrorDialog
from cumulobert.app.find_star_settings_dialog import FindStarSettingsDialog
from cumulobert.app.image_view import ImageView
from cumulobert.app.load_actions import (
    loadFileMenuActions,
    reduceActions,
)
from cumulobert.app.query_catalogue_dialog import QueryCatalogueDialog
from cumulobert.app.success_dialog import SuccessDialog
from cumulobert.app.stars_table_view import StarsTableView
from cumulobert.app.utils import getFileType
from cumulobert.cluster_image import ClusterImage
from cumulobert.errors import ClusterImageError
from cumulobert.external_catalogue import ExternalCatalogue

class MainWindow(QMainWindow):
    """Main Window

    Methods
    -------
    (see QMainWindow)
    __init__
    _createToolBar
    _createMenuBar
    _createStatusBar
    _loadActions

    Attributes
    ----------
    (see QMainWindow)

    centralWidget: QtWidget
    Central widget

    cluster_image: ClusterImage
    Opened image

    ext_cat: ExternalCatalogue
    External catalogue with sources in our field

    fileActions: list of QAction
    List of file action items. They are plotted in the menu and also in the toolbar

    imageView: ImageView
    Image window to show the cluster image

    reduceActions: list of QAction
    List of reducction action items. They are plotted in the menu and also in the toolbar

    star_finder_options_set: bool
    Flag specifying if the star finder options are set
    """
    def __init__(self):
        """Initialize class instance """
        super().__init__()

        self.setWindowTitle("Cumulobert")
        self.resize(WIDTH, HEIGHT)

        self.centralWidget = QLabel("Welcome to Cumulobert")
        self.centralWidget.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setCentralWidget(self.centralWidget)

        self.reduceActions = reduceActions(self)
        self.fileActions = loadFileMenuActions(self)

        self._createToolBar()
        self._createStatusBar()
        self._createMenuBar()

        self.cluster_image = None
        self.ext_cat = None
        self.imageView = None
        self.star_finder_options_set = False

    def _createToolBar(self):
        """Create tool bars"""
        fileToolBar = QToolBar("File toolbar")
        fileToolBar.setIconSize(QSize(ICON_SIZE, ICON_SIZE))
        for menuAction in self.fileActions:
            fileToolBar.addAction(menuAction)
            fileToolBar.addSeparator()
        self.addToolBar(fileToolBar)

        reduceToolBar = QToolBar("Reduce toolbar")
        reduceToolBar.setIconSize(QSize(ICON_SIZE, ICON_SIZE))
        for menuAction in self.reduceActions:
            reduceToolBar.addAction(menuAction)
            reduceToolBar.addSeparator()
        self.addToolBar(reduceToolBar)

    def _createMenuBar(self):
        """Create menu bars"""
        menu = self.menuBar()

        fileMenu = menu.addMenu("&File")
        for menuAction in self.fileActions:
            fileMenu.addAction(menuAction)
            fileMenu.addSeparator()

        reduceMenu = menu.addMenu("&Reduce")
        for menuAction in self.reduceActions:
            reduceMenu.addAction(menuAction)
            reduceMenu.addSeparator()

    def _createStatusBar(self):
        """Create status bar"""
        self.setStatusBar(QStatusBar(self))

    @pyqtSlot()
    def extractStars(self):
        """Open dialog to query Vizier catalogue"""
        if self.cluster_image.star_cat is None:
            error = "Missing stars"
            errorDialog = ErrorDialog(
                "An error occurred when extracting stellar information:\n" + str(error))
            errorDialog.exec()
            return

        # sort table
        star_cat = self.cluster_image.star_cat
        # sort by RA/Dec
        if self.imageView.pix_to_radec:
            star_cat.sort(['ra', 'dec'])
        # sort by x/y
        else:
            star_cat.sort(['xcentroid', 'ycentroid'])


        windowWidth = self.frameGeometry().width()
        windowHeight = self.frameGeometry().height()

        # plot table
        self.starsTableView = StarsTableView(
            self.cluster_image.star_cat)
        self.starsTableView.resize(windowWidth, windowHeight)
        self.starsTableView.show()

    @pyqtSlot()
    def findStars(self):
        """Find stars in the image"""
        if self.ext_cat is None:
            error = "Missing external catalogue"
            errorDialog = ErrorDialog(
                "An error occurred when finding stars:\n" + str(error))
            errorDialog.exec()
            return

        try:
            # Add query to change these two values accordingly to the
            # cluster/target you have observed
            findStarSettingDialog = FindStarSettingsDialog()

            if findStarSettingDialog.exec():
                try:
                    threshold = float(findStarSettingDialog.thresholdQuestion.text())
                except ValueError as error:
                    errorDialog = ErrorDialog("Selected threshold must be a float")
                    errorDialog.exec()
                    return
                try:
                    fwhm = float(findStarSettingDialog.fwhmQuestion.text())
                except ValueError as error:
                    errorDialog = ErrorDialog("Selected fwhm must be a float")
                    errorDialog.exec()
                    return
                try:
                    minsep_fwhm = float(findStarSettingDialog.minsepFwhmQuestion.text())
                except ValueError as error:
                    errorDialog = ErrorDialog("Selected minsep_fwhm must be a float")
                    errorDialog.exec()
                    return

                self.cluster_image.find_stars(
                    self.ext_cat,
                    threshold=threshold,
                    fwhm=fwhm,
                    minsep_fwhm=minsep_fwhm)

                aper_radius_list = self.cluster_image.aper_radius_arrays(
                    self.imageView.pix_to_radec)
                self.imageView.aperRadiusList = aper_radius_list
                annular_radius_list = self.cluster_image.annular_radius_arrays(
                    self.imageView.pix_to_radec)
                self.imageView.annularRadiusList = annular_radius_list
                self.imageView.updatePlot()

                successDialog = SuccessDialog(
                    f"Detected {len(self.cluster_image.star_cat)} stars")
                successDialog.exec()

        except Exception as error:
            errorDialog = ErrorDialog(
                "An error occurred when finding stars:\n" + str(error))
            errorDialog.exec()
            return

    @pyqtSlot()
    def openFile(self):
        """Open dialog to select and open file"""
        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Open File",
            "${HOME}",
            "Fits Files (*fit *.fits *.fits.gz);; Data (*dat);; All files (*)",
        )

        # figure out whether to open an Image or a Spectrum
        file_type = getFileType(filename)

        # Unkonw extension, report message in the status bar
        if file_type is None:
            message = "Unrecognized file extension"
            self.statusBar().showMessage(message)

        # open Image
        elif file_type == "Image":
            try:
                # load image
                self.cluster_image = ClusterImage(filename)

                # plot image
                self.imageView = ImageView(self.cluster_image)
                self.setCentralWidget(self.imageView)

                # enable stellar extraction options
                for action in self.reduceActions:
                    action.setEnabled(True)

                # Add a status bar to show mouse position
                self.status_bar = QStatusBar()
                self.setStatusBar(self.status_bar)

                # Connect mouse move event to a custom slot
                self.imageView.scene().sigMouseMoved.connect(self.on_mouse_move)

            except ClusterImageError as error:
                errorDialog = ErrorDialog(
                    "An error occurred when opening an image:\n" + str(error))
                errorDialog.exec()

        # open Table
        elif file_type == "Table":
            error = "Not implemented"
            errorDialog = ErrorDialog(
                "An error occurred when opening an Table:\n" + str(error))
            errorDialog.exec()

    def on_mouse_move(self, pos):
        # Get the mouse position in the view coordinates
        mouse_point = self.imageView.plotItem.vb.mapSceneToView(pos)
        x = mouse_point.x()
        y = mouse_point.y()

        # Update the status bar with the mouse position
        # TODO: fix this if RA/Dec is show instead of X/Y
        self.status_bar.showMessage(f"X: {x:.2f}, Y: {y:.2f}")

    @pyqtSlot()
    def queryCatalogue(self):
        """Open dialog to query Vizier catalogue"""

        try:
            # Add query to change these two values accordingly to the
            # cluster/target you have observed
            addRaDecDialog = QueryCatalogueDialog()

            if addRaDecDialog.exec():
                try:
                    ra = float(addRaDecDialog.raQuestion.text())
                except ValueError as error:
                    errorDialog = ErrorDialog("Selected RA must be a float")
                    errorDialog.exec()
                    return
                try:
                    dec = float(addRaDecDialog.decQuestion.text())
                except ValueError as error:
                    errorDialog = ErrorDialog("Selected Dec must be a float")
                    errorDialog.exec()
                    return
                try:
                    field_size = float(addRaDecDialog.sizeQuestion.text())
                    # convert from minutes to degrees
                    field_size /= 60.0
                except ValueError as error:
                    errorDialog = ErrorDialog("Selected field size must be a float")
                    errorDialog.exec()
                    return

                self.ext_cat = ExternalCatalogue(ra, dec, field_size)

                successDialog = SuccessDialog("Query success")
                successDialog.exec()

        except Exception as error:
            errorDialog = ErrorDialog(
                "An error occurred when querying Vizier catalogue:\n" + str(error))
            errorDialog.exec()
            return
