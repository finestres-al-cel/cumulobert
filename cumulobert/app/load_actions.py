""" Functions to load Actions"""
from PyQt6.QtGui import QAction, QIcon

from cumulobert.app.environment import BUTTONS_PATH

def loadFileMenuActions(window):
    """Load file menu actions

    Arguments
    ---------
    window: MainWindow
    Window where the actions will act

    Return
    ------
    menuAction: list of QAction
    List of actions in the file menu
    """
    menuActions = []

    load_spectrum_option = QAction(
        QIcon(f"{BUTTONS_PATH}/load_image.png"),
        "&Load Image",
        window)
    load_spectrum_option.setStatusTip("Load Image")
    load_spectrum_option.triggered.connect(window.openFile)
    menuActions.append(load_spectrum_option)

    return menuActions

def reduceActions(window):
    """Load reduce menu actions

    Arguments
    ---------
    window: MainWindow
    Window where the actions will act

    Return
    ------
    menuAction: list of QAction
    List of actions in the spectral extraction menu
    """
    menuActions = []

    query_catalogue_option = QAction(
        QIcon(f"{BUTTONS_PATH}/query_catalogue.jpg"),
        "&Query Catalogue",
        window)
    query_catalogue_option.setStatusTip("Query Catalogue")
    query_catalogue_option.triggered.connect(window.queryCatalogue)
    query_catalogue_option.setCheckable(True)
    query_catalogue_option.setEnabled(False)
    menuActions.append(query_catalogue_option)

    find_stars_option = QAction(
        QIcon(f"{BUTTONS_PATH}/find_stars.png"),
        "&Find Stars",
        window)
    find_stars_option.setStatusTip("Find Stars")
    find_stars_option.triggered.connect(window.findStars)
    find_stars_option.setEnabled(False)
    menuActions.append(find_stars_option)

    extract_stars_option = QAction(
        QIcon(f"{BUTTONS_PATH}/extract_stars.png"),
        "&Extract Stars",
        window)
    extract_stars_option.setStatusTip("Extract Stars")
    extract_stars_option.triggered.connect(window.extractStars)
    extract_stars_option.setEnabled(False)
    menuActions.append(extract_stars_option)

    invert_axis_option = QAction(
        QIcon(f"{BUTTONS_PATH}/invert_axis.png"),
        "&Invert Axis",
        window)
    invert_axis_option.setStatusTip("Invert Axis")
    invert_axis_option.triggered.connect(window.invertAxis)
    invert_axis_option.setEnabled(False)
    menuActions.append(invert_axis_option)



    return menuActions
