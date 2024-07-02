"""Utils file contianing useful functions"""

from cumulobert.cluster_image import ACCEPTED_FORMATS as ACCEPTED_FORMATS_IMAGE
#from cumulobert.table import ACCEPTED_FORMATS as ACCEPTED_FORMATS_TABLE

def getFileType(filename):
    """Figure out if the file is an Image or a Spectrum

    Make the decision based on the filename extension

    Attributes
    ----------
    filename: str
    Name of the file

    Return
    ------
    file_type: str or None
    'Image' or 'Table' if the extension is in ACCEPTED_FORMATS_IMAGE or
    ACCEPTED_FORMATS_TABLE respectively. None otherwise
    """
    file_type = None
    for format in ACCEPTED_FORMATS_IMAGE:
        if filename.endswith(format):
            file_type = "Image"
    #for format in ACCEPTED_FORMATS_TABLE:
    #    if filename.endswith(format):
    #        file_type = "Table"

    return file_type
