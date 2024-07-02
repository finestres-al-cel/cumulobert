""" External Catalogue """
from astroquery.vizier import Vizier
from astropy import units
from astropy.coordinates import Angle, SkyCoord

class ExternalCatalogue:
    """ External catalogue with sources in the field

    Methods
    -------
    __init__

    Attributes
    ----------

    """
    def __init__(self, ra, dec, field_size):
        """Construct cross-match catalogue

        Arguments
        ---------
        ra: float
        Central field Right Ascention (in degrees)

        dec: float
        Central field declination (in degrees)
        """
        # Here we set the columns we want to keep from the catalogue.
        # This way we do not download all the info which, in some cases could
        # be quite memory demanding.
        # You can add/remove any columns you think might be of use or not.
        vizier = Vizier(
            columns=["RAJ2000", "DEJ2000", "TIC", "GAIA", "Tmag", "Bmag", "Vmag"],
            catalog="IV/39/tic82",
            row_limit=-1)
        self.cat = vizier.query_region(
            SkyCoord(ra=ra, dec=dec, unit=(units.deg, units.deg), frame='icrs'),
            width=Angle(field_size, "deg"),
            height=Angle(field_size, "deg"),
            catalog="IV/39/tic82",
            column_filters={'Tmag': '<13.0'})

        # These are the coordinates of the stars from the external catalogue
        # we have found in the previous query.
        # They will be used to crossmatch our detections in each frame with
        # the external catalogue.
        self.sky_coords = SkyCoord(self.cat[0]['RAJ2000'], self.cat[0]['DEJ2000'])
        
