Usage instructions:
1. Go to [astrometry.net](https://nova.astrometry.net/) 
2. Upload your image. Then download the `new_image.fits` file and rename it appropriately.
3. Keep the center coordinates (in degrees) and the field size (in arcmin)
4. Redo 2 and 3 for the other filter
5. Open a new terminal
6. Go to cumulobert `bin` folder:
    `cd <path to your cumulobert repo>/bin`
7. Load your conda environment (see README.md):
    `conda activate my_cumulobert_env`
8. Open the cumulobert app to extract the photometry of your cluster
    `python cumulobert_app.py`
    1. Load the file
    2. Query the catalogue (you will need to input the center coordinates and field size you saved in step 3)
    3. Find stars
    4. Extract stars
    5. Save the catalogue
9. Redo 8 for the other filter
10. Open jupyter lab
    `jupyter-lab`
11. Open the CumulObert.ipynb file and finish its instructions to finish the analysis