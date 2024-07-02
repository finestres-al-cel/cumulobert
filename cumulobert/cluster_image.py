""" Basic Cluster Image """
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.table import Table, hstack, join, vstack, unique
from astropy.wcs import WCS
from photutils.detection import find_peaks
from photutils.detection import IRAFStarFinder
from photutils.aperture import (
    aperture_photometry, ApertureStats, CircularAnnulus, CircularAperture)
from scipy.ndimage import rotate
import numpy as np

from cumulobert.errors import ClusterImageError

ACCEPTED_FORMATS = [".fit", ".fits", ".fits.gz"]

NUM_POINTS_CIRCLE = 1000
ANGLES = np.linspace(0, 2*np.pi, 1000)

class ClusterImage:
    """ Cluster

    Methods
    -------
    __init__
    load_image
    star_finder

    Attributes
    ----------
    filename: str
    Name of the file containing the cluster image

    data: array of float
    The 2D cluster image data

    dec: array of float
    The declinations corresponding to the y-pixels in data

    header: astropy.io.fits.header.Header
    The image header

    image_extension: str
    Extension of the loaded file

    ra: array of float
    The right ascension corresponding to the x-pixels in data

    wcs: astropy.wcs.WCS
    This loads the astrometric information stored in the header.
    Allows to convert from pixel to sky coordinates
    """
    def __init__(self, filename):
        """Initialize instance

        Arguments
        ---------
        filename: str
        Filename to open

        Raise
        -----
        ClusterImageError if filename is not a string
        ClusterImageError if filename does not have the correct format
        """
        # check filename type
        if not isinstance(filename, str):
            raise ClusterImageError(
                f"Image: Argument 'filename' has incorrect type. Expected string "
                f"found {type(filename)}. {filename}")

        # check filename extension
        self.image_extension = None
        for format_check in ACCEPTED_FORMATS:
            if filename.endswith(format_check):
                self.image_extension = format_check
        if self.image_extension is None:
            raise ClusterImageError(
                "Image: 'filename' has incorrect extension. Valid"
                "extensions are " + ", ".join(ACCEPTED_FORMATS)
                )

        self.filename = filename

        # load image
        self.data = None
        self.dec = None
        self.header = None
        self.ra = None
        self.wcs = None
        self.load_image()

        self.star_cat = None

    def load_image(self):
        """ Load the cluster image"""
        try:
            hdu = fits.open(self.filename)

            # read data
            self.header = hdu[0].header
            self.data = hdu[0].data

            # This loads the astrometric information stored in the header.
            # Allows to convert from pixel to sky coordinates
            self.wcs = WCS(self.header)

            # pixel-to-sky transformation
            x_pix = np.arange(self.data.shape[0]) # x pixels
            y_pix = np.arange(self.data.shape[1]) # y pixels

            y_aux = np.array(y_pix[0]*x_pix.size)
            aux = self.wcs.pixel_to_world(x_pix, y_aux)
            self.ra = [item.ra.value for item in aux]
            x_aux = np.array(y_pix[0]*x_pix.size)
            aux = self.wcs.pixel_to_world(x_aux, y_pix)
            self.dec = [item.dec.value for item in aux]

        except IOError as error:
            raise ClusterImageError("ClusterImage:", str(error)) from error

        except KeyError as error:
            raise ClusterImageError("ClusterImage:", str(error)) from error

        finally:
            hdu.close()

    def find_stars(self, ext_cat, threshold=50, fwhm=8, minsep_fwhm=2.5):
        """Find stars in the image

        Arguments
        ---------
        ext_cat: ExternalCatalogue
        External catalogue with information about the sources in the field

        threshold: float
        The absolute image value above which to select sources.
        (see https://photutils.readthedocs.io/en/stable/api/photutils.detection.IRAFStarFinder.html
        for details)

        fwhm: float
        The full-width half-maximum (FWHM) of the 2D circular Gaussian kernel
        in units of pixels.
        (see https://photutils.readthedocs.io/en/stable/api/photutils.detection.IRAFStarFinder.html
        for details)

        minsep_fwhm: float
        The separation (in units of fwhm) for detected objects. The minimum
        separation is calculated as int((fwhm * minsep_fwhm) + 0.5) and is
        clipped to a minimum value of 2. Note that large values may result in
        long run times.
        (see https://photutils.readthedocs.io/en/stable/api/photutils.detection.IRAFStarFinder.html
        for details)
        """
        # In this part of the code we use Photutils
        # (https://photutils.readthedocs.io/en/stable/) to detect all possible
        # sources within our image.
        # IRAFStarFinder searches images for local density maxima that have a
        # peak amplitude greater than threshold above the local background and
        # have a PSF full-width at half-maximum similar to the input fwhm
        # (https://photutils.readthedocs.io/en/stable/api/photutils.detection.IRAFStarFinder.html)
        # NOTE: you can change the values of threshold and fwhm to change the
        # number of detected objects. You can keep the other parameters set to
        # their default values.
        find = IRAFStarFinder(
            threshold,
            fwhm,
            minsep_fwhm=minsep_fwhm,
            roundlo=0.01,
            roundhi=1,
            sharphi=1,
            exclude_border=True)

        sources = find(self.data)
        # We add filter, exptime and JD info to the sources dataframe
        sources['filter'] = self.header['FILTER']
        sources['exptime'] = self.header['EXPTIME']
        sources['dateobs'] = self.header['DATE-OBS']

        peak_median = np.median(sources['peak'])
        starcat = []
        for source in sources:
            # Here we set the positions in the image (i.e. the centroid of
            # the detected objects) where we will do the aperture photometry
            star_pos = np.transpose((source['xcentroid'], source['ycentroid']))
            aper_radius = source["fwhm"]*2 + 0.09*source['peak']/peak_median
            circular_aper = CircularAperture(star_pos, r=aper_radius)
            circular_annulus_inner = aper_radius + 2
            circular_annulus_outer = aper_radius + 4
            circular_annulus = CircularAnnulus(
                star_pos, r_in=circular_annulus_inner, r_out=circular_annulus_outer)
            star_data = aperture_photometry(self.data, circular_aper)
            star_data["aper_radius"] = aper_radius
            star_data["circular_annulus_inner"] = circular_annulus_inner
            star_data["circular_annulus_outer"] = circular_annulus_outer

            # Here the local background inside the annulus is computed
            aperstats = ApertureStats(self.data, circular_annulus)
            bkg_mean = aperstats.mean
            aperture_area = circular_aper.area_overlap(self.data)
            total_bkg = bkg_mean * aperture_area

            # Here the final background-substracted flux, its error and the
            # instrumental magnitude are computed
            star_data['aperture_sum_bkgsub'] = star_data['aperture_sum'] - total_bkg
            star_data['fluxerr'] = np.sqrt(star_data['aperture_sum_bkgsub']) # We assume the flux error as sqrt(flux)
            star_data['inst. mag'] = -2.5* np.log10(star_data['aperture_sum_bkgsub'] / source['exptime'])
            star_data['filter'] = source['filter']
            star_data['exptime'] = source['exptime']
            star_data['dateobs'] = source['dateobs']

            # Here we add the corresponding (RA, Dec) coordinates for each
            # target taking into account their (X,Y) coordinates and the WCS
            # info from the header
            sky = self.wcs.pixel_to_world(source['xcentroid'], source['ycentroid'])
            star_data['xcentroid'] = source['xcentroid']
            star_data['ycentroid'] = source['ycentroid']
            star_data['ra'] = sky.ra
            star_data['dec'] = sky.dec

            starcat.append(star_data)

        # Here we crossmatch the list of detected sources with the external
        # catalogue (i.e. TESS Input Catalogue, Gaia, ...)
        starcat = vstack(starcat, join_type='exact')

        frame_cat = SkyCoord(starcat['ra'], starcat['dec'])
        idx_ext_cat, d2d_ext_cat, _ = frame_cat.match_to_catalog_sky(
            ext_cat.sky_coords)

        cross_frame_cat = hstack([starcat, ext_cat.cat[0][idx_ext_cat]])
        cross_frame_cat['distance'] = d2d_ext_cat

        # To avoid spurious detections (i.e. hot pixels, wrong crossmatches,
        # etc) we only keep those objects whose angular distance to the external
        # catalogue ones are less than 0.001 deg
        mask = cross_frame_cat['distance'].value < 0.001
        self.star_cat = cross_frame_cat[mask]

        # add flag to set reference stars
        self.star_cat["ref mag"] = np.nan

    def annular_radius_arrays(self, pix_to_radec=False):
        """ Find the points building the annular radius used for the background
        subtraction in each of the identified stars

        Arguments
        ---------
        pix_to_radec: bool - Default: False
        If True, return the points in RA/Dec coordinates. Otherwise, return
        pixel positions

        Return
        ------
        annular_radius_list: list of (float, float) or list of SkyCoord
        The found points. Format for each item in (x,y) if pix_to_radec is False
        and SkyCoord otherwise
        """
        annular_radius_list = []
        for item in self.star_cat:

            x_center = item["xcentroid"]
            y_center = item["ycentroid"]

            # inner circle
            inner_radius = item["circular_annulus_inner"]
            x_circle_inner = inner_radius*np.cos(ANGLES) + x_center
            y_circle_inner = inner_radius*np.sin(ANGLES) + y_center

            # outer circle
            outer_radius = item["circular_annulus_outer"]
            x_circle_outer = outer_radius*np.cos(ANGLES) + x_center
            y_circle_outer = outer_radius*np.sin(ANGLES) + y_center

            if pix_to_radec:
                annular_radius_list += [
                    item for item in self.wcs.pixel_to_world(
                        x_circle_inner, y_circle_inner)
                ]
                annular_radius_list += [
                    item for item in self.wcs.pixel_to_world(
                        x_circle_outer, y_circle_outer)
                ]
            else:
                annular_radius_list += [
                    (x, y) for (x, y) in zip(x_circle_inner, y_circle_inner)
                ]
                annular_radius_list += [
                    (x, y) for (x, y) in zip(x_circle_outer, y_circle_outer)
                ]

            """
            # inner circle
            inner_radius = item["circular_annulus_inner"]
            x_circle_inner = np.linspace(
                x_center - inner_radius, x_center + inner_radius, NUM_POINTS_CIRCLE)
            y_inner_aux = np.sqrt(inner_radius*inner_radius - (x_circle_inner - x_center)**2)
            y_circle_inner_top = y_center + y_inner_aux
            y_circle_inner_bottom = y_center - y_inner_aux

            # outer circle
            outer_radius = item["circular_annulus_outer"]
            x_circle_outer = np.linspace(
                x_center - outer_radius, x_center + outer_radius, NUM_POINTS_CIRCLE)
            y_outer_aux = np.sqrt(outer_radius*outer_radius - (x_circle_outer - x_center)**2)
            y_circle_outer_top = y_center + y_outer_aux
            y_circle_outer_bottom = y_center - y_outer_aux

            # add to list
            if pix_to_radec:
                annular_radius_list += [
                    item for item in self.wcs.pixel_to_world(
                        x_circle_inner, y_circle_inner_top)
                ]
                annular_radius_list += [
                    item for item in self.wcs.pixel_to_world(
                        x_circle_inner, y_circle_inner_bottom)
                ]
                annular_radius_list += [
                    item for item in self.wcs.pixel_to_world(
                        x_circle_outer, y_circle_outer_top)
                ]
                annular_radius_list += [
                    item for item in self.wcs.pixel_to_world(
                        x_circle_outer, y_circle_outer_bottom)
                ]
            else:
                annular_radius_list += [
                    (x, y) for (x, y) in zip(
                        x_circle_inner, y_circle_inner_top)
                ]
                annular_radius_list += [
                    (x, y) for (x, y) in zip(
                        x_circle_inner, y_circle_inner_bottom)
                ]
                annular_radius_list += [
                    (x, y) for (x, y) in zip(
                        x_circle_outer, y_circle_outer_top)
                ]
                annular_radius_list += [
                    (x, y) for (x, y) in zip(
                        x_circle_outer, y_circle_outer_bottom)
                ]
            """

        return annular_radius_list

    def aper_radius_arrays(self, pix_to_radec=False):
        """ Find the points building the aperture radius used for the flux
        estimation for each of the identified stars

        Arguments
        ---------
        pix_to_radec: bool - Default: False
        If True, return the points in RA/Dec coordinates. Otherwise, return
        pixel positions

        Return
        ------
        aper_radius_list: list of (float, float) or list of SkyCoord
        The found points. Format for each item in (x,y) if pix_to_radec is False
        and SkyCoord otherwise
        """
        aper_radius_list = []
        for item in self.star_cat:

            x_center = item["xcentroid"]
            y_center = item["ycentroid"]
            radius = item["aper_radius"]

            x_circle = radius*np.cos(ANGLES) + x_center
            y_circle = radius*np.sin(ANGLES) + y_center
            if pix_to_radec:
                aper_radius_list += [
                    item for item in self.wcs.pixel_to_world(x_circle, y_circle)
                ]
            else:
                aper_radius_list += [
                    (x, y) for (x, y) in zip(x_circle, y_circle)
                ]


        return aper_radius_list
