import sys, glob, os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt


from astropy.io import fits
from astropy.wcs import WCS
from astroquery.vizier import Vizier
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table, hstack, join, vstack, unique
from astropy.visualization import ZScaleInterval, SinhStretch, ImageNormalize

from photutils.detection import find_peaks
from photutils.detection import IRAFStarFinder
from photutils.aperture import aperture_photometry, ApertureStats, CircularAnnulus, CircularAperture

target_files = glob.glob('new*.fits')
print('Number of available files: ', len(target_files), '\n')

ra, dec = (276.771, 6.399) #Change these two values accordingly to the cluster/target you have observed

print("Querying catalogue around target position: ", ra, dec, "...\n")

"""
Here we set the columns we want to keep from the catalogue. This way we do not download all the info which, in some cases could be quite memory demanding.
You can add/remove any columns you think might be of use or not. See NOTE at the end of STEP 3.
"""
v = Vizier(columns=["RAJ2000", "DEJ2000", "TIC", "GAIA", "Tmag", "Bmag", "Vmag"], catalog="IV/39/tic82", row_limit=-1)
result = v.query_region(SkyCoord(ra=ra, dec=dec,\
                                 unit=(u.deg, u.deg), frame='icrs'), width=Angle(21.4/60, "deg"), height=Angle(21.4/60, "deg"), catalog="IV/39/tic82", column_filters={'Tmag': '<13.0'})

"""
These are the coordinates of the stars from the external catalogue we have found in the previous query.
They will be used to crossmatch our detections in each frame with the external catalogue.
"""

vizier_cat = SkyCoord(result[0]['RAJ2000'], result[0]['DEJ2000'])

allcats = []
for i in target_files:
    print('===================================================================================')
    print('Running aperture photometry for file: ', i)

    """
    Here we open the files and load their corresponding data, header and WCS (astrometry) info (check https://docs.astropy.org/en/stable/io/fits/)
    """
    hdulist = fits.open(i, ignore_missing_end = True)
    data = hdulist[0].data #This is the actual pixel array of the image
    hdr = hdulist[0].header #This loads the header info
    wcs = WCS(hdr) #This loads the astrometric information stored in the header. Allows to convert from pixel to sky coordinates
    hdulist.close()

    """
    In this part of the code we use Photutils (https://photutils.readthedocs.io/en/stable/) to detect all possible sources within our images.
    IRAFStarFinder searches images for local density maxima that have a peak amplitude greater than threshold above the local background and
    have a PSF full-width at half-maximum similar to the input fwhm (https://photutils.readthedocs.io/en/stable/api/photutils.detection.IRAFStarFinder.html#photutils.detection.IRAFStarFinder).
    NOTE: you can change the values of threshold and fwhm to change the number of detected objects. You can keep the other parameters set to their default values.
    """
    find = IRAFStarFinder(50, 8, minsep_fwhm=2, roundlo=0.01, roundhi=1, sharphi=1, exclude_border=True)
    sources = find(data)
    sources['filter'] = hdr['FILTER']    # We add filter, exptime and JD info to the sources dataframe
    sources['exptime'] = hdr['EXPTIME']  #
    sources['dateobs'] = hdr['DATE-OBS']       #

    """
    This part is OPTIONAL but can be used to check if the objects we have detected are correct or not using SAOImageDS9. The next step creates a "region" file
    that can be overploted together with your frame. Open the file in DS9 and then "region --> open" and select the .reg file created here.
    """
    outfile = open('./' + os.path.basename(i.replace('.fits','.reg')), 'w') #Change the path accordingly to your local path where you want to store the files
    outfile.write('global color=green dashlist=8 3 width=2 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    outfile.write("physical\n")
    for i in range(0, len(sources)):
        outfile.write("box(%s,%s,%s,%s,0)\n" %(sources["xcentroid"][i], sources["ycentroid"][i], 10, 10))
    outfile.close()

    """
    Here is where the aperture photometry is performed using the procedure explained in https://photutils.readthedocs.io/en/stable/aperture.html. In the example below
    we are following the Local Background Subtraction method using a circular annulus around each of our targets.

    NOTE 1: we iteratively set the aperture radius of each target. The radius of the circular aperture for each star is defined in "aper_radius". We use the fwhm and
    peak data to set a slightly different aperture radius for each target depending on their size and magnitude. You can play with this "aper_radius" as well as the
    "CircularAnnulus" radii to see how the photometry for each target changes.
    NOTE 2: you might also want to change the code to use a Global Background Subtraction method and compare results.
    NOTE 3: the code also plots the apertures and annulus around each target so you can check if they are correctly set.
    """

    fig = plt.figure(figsize=(15,15))                           # You can comment this part of the code if
    ax = plt.gca()                                              # you do not want to plot the apertures
    norm = ImageNormalize(data, interval=ZScaleInterval(),      #
                          stretch=SinhStretch())                #
    plt.imshow(data, cmap='gray', norm=norm, origin='lower')    #

    starcat = []
    plot = True
    for i in sources:
        """
        Here we set the positions in the image (i.e. the centroid of the detected objects) where we will do the aperture photometry
        """
        star_pos = np.transpose((i['xcentroid'], i['ycentroid']))
        aper_radius = i["fwhm"]*2 + 0.09*i['peak']/np.median(sources['peak'])
        star_aper = CircularAperture(star_pos, r=aper_radius)
        star_annulus = CircularAnnulus(star_pos, r_in=aper_radius+2, r_out=aper_radius+4)
        star_data = aperture_photometry(data, star_aper)
        star_data["aper_radius"] = aper_radius

        """
        Here the local background inside the annulus is computed
        """
        aperstats = ApertureStats(data, star_annulus)
        bkg_mean = aperstats.mean
        aperture_area = star_aper.area_overlap(data)
        total_bkg = bkg_mean * aperture_area

        """
        Here the final background-substracted flux, its error and the instrumental magnitude are computed
        """
        star_data['aperture_sum_bkgsub'] = star_data['aperture_sum'] - total_bkg
        star_data['fluxerr'] = np.sqrt(star_data['aperture_sum_bkgsub']) # We assume the flux error as sqrt(flux)
        star_data['inst. mag'] = -2.5*np.log10(star_data['aperture_sum_bkgsub']/i['exptime'])
        star_data['filter'] = i['filter']
        star_data['exptime'] = i['exptime']
        star_data['dateobs'] = i['dateobs']

        """
        Here we add the corresponding (RA, Dec) coordinates for each target taking into account their (X,Y) coordinates and the WCS info from the header
        """
        sky = wcs.pixel_to_world(i['xcentroid'], i['ycentroid'])
        star_data['ra'] = sky.ra
        star_data['dec'] = sky.dec

        starcat.append(star_data)

        star_aper.plot(color='red', lw=1.5, alpha=0.5)
        star_annulus.plot(color='green', lw=1.5, alpha=0.5)

    #ax.set_xlim(1900, 2100) # You can change this limits to check different parts of the image
    #ax.set_ylim(1900, 2300) # You can comment this part of the code if
    plt.show()              # you do not want to plot the apertures
    plt.close()             #

    """
    Here we crossmatch the list of detected sources with the external catalogue (i.e. TESS Input Catalogue, Gaia, ...)
    """
    starcat = vstack(starcat, join_type='exact')

    frame_cat = SkyCoord(starcat['ra'], starcat['dec'])
    idx_vizier, d2d_vizier, d3d_vizier = frame_cat.match_to_catalog_sky(vizier_cat)

    cross_frame_cat = hstack([starcat, result[0][idx_vizier]])
    cross_frame_cat['distance'] = d2d_vizier
    #cross_frame_cat.keep_columns(["TIC", "GAIA", "dateobs", "xcenter", "ycenter" ,"ra", "dec", "RAJ2000", "DEJ2000", "aperture_sum_bkgsub", "fluxerr", "inst. mag", "exptime", "filter", "Tmag", "Bmag", "Vmag", "distance"])
    """
    NOTE: in the previous line, "keep_columns" limits the total number of final columns that our data table will have. You can comment this line to keep all the info from the external catalogue
    or you can add any extra columns you see fit. For example, for TIC, keeping the proper motions "pmRA" & "pmDEC" or the parallax "plx" can help you discriminate objects which are really part
    of the cluster or not. WARNING: if you use TIC, the proper motions, and parallax come from Gaia DR2. So it is recommended to crossmatch again the catalogue using the Gaia ids
    to obtain more up-to-date Gaia DR3 data. WARNING 2: You will have to add any extra value in the columns parameter in here too (see STEP 2): v = Vizier(columns=["RAJ2000", "DEJ2000", "TIC", "GAIA", "Tmag", "Bmag", "Vmag"], catalog="IV/39/tic82", row_limit=-1)
    """

    """
    To avoid spurious detections (i.e. hot pixels, wrong crossmatches, etc) we only keep those objects whose angular distance to the external catalogue ones are less than 0.001 deg
    """
    mask = cross_frame_cat['distance'].value < 0.001 *1e10
    filt_cat = cross_frame_cat[mask] #This contains all the "good" detections in the frame that have been correctly crossmatched with the external catalogue
                                     #with all their associated info. You can print it to check its contents.
    allcats.append(filt_cat)

    filt_cat.write("results_dani.csv")
