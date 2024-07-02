"""Error Dialog """
from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QGridLayout, QLabel, QLineEdit, QPushButton, QVBoxLayout
)

class SuccessDialog(QDialog):
    """ Class to define the dialog to report success

    Methods
    -------
    (see QDialog)
    __init__

    Arguments
    ---------
    (see QDialog)
    """
    def __init__(self, message):
        """Initialize instance

        Arguments
        ---------
        message: str
        Error message
        """
        super().__init__()

        self.setWindowTitle("Success!")

        QButtons = QDialogButtonBox.StandardButton.Ok

        self.buttonBox = QDialogButtonBox(QButtons)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(QLabel(message))
        layout.addWidget(self.buttonBox)

        self.setLayout(layout)
