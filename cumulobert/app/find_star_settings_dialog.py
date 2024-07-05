""" Dialog to add calibration points"""
from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QGridLayout, QLabel, QLineEdit
)

class FindStarSettingsDialog(QDialog):
    """ Class to define the settings for the stellar finder

    Methods
    -------
    (see QDialog)
    __init__

    Arguments
    ---------
    (see QDialog)

    buttonBox: QDialogButtonBox
    Accept/cancel button

    fwhmQuestion: QLineEdit
    Field to modify the FWHM for the Gaussian Kernel used for the detection

    minsepFwhmQuestion: QLineEdit
    Field to modify the minimum separation between detected objects (in units
    of the FWHM)

    thresholdQuestion: QLineEdit
    Field to modify the detection threshold
    """
    def __init__(self):
        """Initialize instance"""
        super().__init__()

        self.setWindowTitle("Find stars")

        QButtons = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel

        self.buttonBox = QDialogButtonBox(QButtons)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        self.thresholdQuestion = QLineEdit()
        self.thresholdQuestion.setMaxLength(10)

        self.fwhmQuestion = QLineEdit()
        self.fwhmQuestion.setMaxLength(10)

        self.minsepFwhmQuestion = QLineEdit()
        self.minsepFwhmQuestion.setMaxLength(10)

        layout = QGridLayout()
        layout.addWidget(QLabel("Threshold"), 0, 0)
        layout.addWidget(self.thresholdQuestion, 0, 1)
        layout.addWidget(QLabel("FWHM Gaussian Kernel"), 1, 0)
        layout.addWidget(self.fwhmQuestion, 1, 1)
        layout.addWidget(QLabel("Min. obj. sep. (in units of FWHM)"), 2, 0)
        layout.addWidget(self.minsepFwhmQuestion, 2, 1)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)
