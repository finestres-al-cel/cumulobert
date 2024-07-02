"""Define ImageView widget as an extension of pg.PlotWidget"""
from astropy.visualization import ImageNormalize, ZScaleInterval, SinhStretch
import numpy as np

from PyQt6.QtWidgets import (
    QTableWidget, QTableWidgetItem, QWidget, QVBoxLayout, QFileDialog,
    QMessageBox, QPushButton, QHBoxLayout
)
import pyqtgraph as pg

from cumulobert.app.success_dialog import SuccessDialog

SHOW_COLUMNS = [
    "xcentroid", "ycentroid", "ra", "dec", "filter",
    "inst. mag", "TIC", "GAIA", "Bmag", "Vmag", "is reference", "calib mag"]
MIN_REF_STARS = 3

class StarsTableView(QWidget):
    """ Manage table plotting

    Methods
    -------
    (see QtCore.QAbstractTableModel)
    __init__

    Attributes
    ----------
    (see QtCore.QAbstractTableModel)

    _data: astropy.table.Table
    The data table
    """
    def __init__(self, data):
        """Initialize instance

        Arguments
        ---------
        data: astropy.table.Table
        The data to be shown
        """
        # initialize plotting
        super().__init__()

        # load data
        self.data = data

        # window settings
        self.setWindowTitle("Stellar properties")
        # Create a layout
        layout = QVBoxLayout(self)
        # Create the table widget
        self.table_widget = QTableWidget()
        layout.addWidget(self.table_widget)
        # Create a horizontal layout for the buttons
        button_layout = QHBoxLayout()
        # Create a button to compute the magnitudes
        self.find_mags_button = QPushButton("Compute magnitudes")
        self.find_mags_button.clicked.connect(self.find_mags)
        button_layout.addWidget(self.find_mags_button)
        # Create a button to save the table
        self.save_button = QPushButton("Save Table")
        self.save_button.clicked.connect(self.save_table)
        button_layout.addWidget(self.save_button)

        # Add the button layout to the main layout
        layout.addLayout(button_layout)

        # Set up the table widget
        self.setup_table_widget()

    def find_mags(self):
        """Fit the instrumental magnitudes using the refernce magnitudes"""
        pos = np.where(self.data["is reference"])
        if len(pos[0]) < MIN_REF_STARS:
            errorDialog = ErrorDialog(
                f"Error: I need at least {MIN_REF_STARS} reference stars")
            errorDialog.exec()
            return
        if np.unique(self.data["filter"]).size > 1:
            errorDialog = ErrorDialog(
                "Error: I see there are many filters in this table."
                "I expected only one to be present")
            errorDialog.exec()
            return

        mag_colname = self.data["filter"][0].split(" ")[1]
        mag_colname = f"{mag_colname}mag"
        if mag_colname not in self.data.colnames:
            errorDialog = ErrorDialog(
                "Error: I could not find the magnitude for this filter")
            errorDialog.exec()
            return

        ref_mags = self.data[mag_colname][pos]
        inst_mags = self.data["inst. mag"][pos]

        # The full formula should be something like this
        # m_B = A_1 + B_1*B_inst + C_1*(B_inst-V_inst) + D_1*B_inst^2 + E_1*(B_inst-V_inst)^2
        # m_V = A_2 + B2"*V_inst + C_2*(B_inst-V_inst) + D_2*V_inst^2 + E_2*(B_inst-V_inst)^2
        # However, we have just one filter here so we do a simpler fit
        # m = A + B m_inst
        mag_fit = np.poly1d(np.polyfit(inst_mags, ref_mags, 1))

        self.data["calib mag"] = mag_fit(self.data["inst. mag"])

        successDialog = SuccessDialog("Magnitudes computed!")
        successDialog.exec()

        self.setup_table_widget()


        return

    def setup_table_widget(self):
        """Setup table widget.

        We fill this widget with the table data and connect it with the
        update_table function when changes are made
        """
        # select columns to show
        self.show_columns = [col for col in SHOW_COLUMNS if col in self.data.colnames]

        n_rows = len(self.data)
        n_cols = len(self.show_columns)
        self.table_widget.setRowCount(n_rows)
        self.table_widget.setColumnCount(n_cols)
        self.table_widget.setHorizontalHeaderLabels(self.show_columns)


        # Fill the table with data from the astropy table
        for row in range(len(self.data)):
            for col, colname in enumerate(self.show_columns):
                item = QTableWidgetItem(str(self.data[colname][row]))
                self.table_widget.setItem(row, col, item)

         # Connect cell change signal to update Astropy table
        self.table_widget.itemChanged.connect(self.update_table)

    def save_table(self):
        """Save the current table"""
        file_name, _ = QFileDialog.getSaveFileName(
            self, "Save Table", "", "CSV Files (*.csv);;All Files (*)",
        )
        if file_name:
            self.data.write(file_name, format='csv', overwrite=True)
            QMessageBox.information(self, "Save Table", f"Table saved to {file_name}")

    def update_table(self, item):
        """Update table"""
        row = item.row()
        col = item.column()
        colname = self.show_columns[col]
        self.data[colname][row] = item.text()
