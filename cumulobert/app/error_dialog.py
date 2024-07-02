"""Error Dialog """
from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QGridLayout, QLabel, QLineEdit, QPushButton, QVBoxLayout
)

class ErrorDialog(QDialog):
    """ Class to define the dialog to report errors

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

        self.setWindowTitle("Error")

        QButtons = QDialogButtonBox.StandardButton.Ok

        self.buttonBox = QDialogButtonBox(QButtons)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(QLabel(message))
        layout.addWidget(self.buttonBox)

        self.setLayout(layout)
