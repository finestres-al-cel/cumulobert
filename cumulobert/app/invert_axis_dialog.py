""" Dialog to add calibration points"""
from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QGridLayout, QLabel, QLineEdit
)

class InvertAxisDialog(QDialog):
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

        self.setWindowTitle("Invert Axis")

        QButtons = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel

        self.buttonBox = QDialogButtonBox(QButtons)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        self.whichAxisQuestion = QLineEdit()
        self.whichAxisQuestion.setMaxLength(4)

        layout = QGridLayout()
        layout.addWidget(QLabel("Axis to invert (x, y, both)"), 0, 0)
        layout.addWidget(self.whichAxisQuestion, 0, 1)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)
