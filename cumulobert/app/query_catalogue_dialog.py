""" Dialog to add calibration points"""
from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QGridLayout, QLabel, QLineEdit
)

class QueryCatalogueDialog(QDialog):
    """ Class to define the Vizier query dialog to select RA/Dec

    Methods
    -------
    (see QDialog)
    __init__

    Arguments
    ---------
    (see QDialog)

    buttonBox: QDialogButtonBox
    Accept/cancel button

    raQuestion: QLineEdit
    Field to modify the right ascension

    decQuestion: QLineEdit
    Field to modify the declination

    sizeQuestion: QLineEdit
    Field to modify the field size
    """
    def __init__(self):
        """Initialize instance"""
        super().__init__()

        self.setWindowTitle("Query Vizier")

        QButtons = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel

        self.buttonBox = QDialogButtonBox(QButtons)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        self.raQuestion = QLineEdit()
        self.raQuestion.setMaxLength(10)

        self.decQuestion = QLineEdit()
        self.decQuestion.setMaxLength(10)

        self.sizeQuestion = QLineEdit()
        self.sizeQuestion.setMaxLength(10)

        layout = QGridLayout()
        layout.addWidget(QLabel("Center RA (degrees)"), 0, 0)
        layout.addWidget(self.raQuestion, 0, 1)
        layout.addWidget(QLabel("Center Dec (degrees)"), 1, 0)
        layout.addWidget(self.decQuestion, 1, 1)
        layout.addWidget(QLabel("Field size (arcminutes)"), 2, 0)
        layout.addWidget(self.sizeQuestion, 2, 1)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)
