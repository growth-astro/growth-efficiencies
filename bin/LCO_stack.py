#!/usr/bin/env python

import os, sys, optparse, shutil, glob, copy
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from scipy.ndimage import median_filter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import aplpy

import ztfsub.utils

from skimage.feature import register_translation
import image_registration

# run LCO_stack.py -i /Users/tahumada/growth-too-papers/data/GW200115/ZTF20aafqvyc/r/elp*fits.fz -o /Users/tahumada/growth-too-papers/data/GW200115/ZTF20aafqvyc/r/stack_r.fits -r 56.9926 -d 38.4422 --doAstrometryNet --doOverwrite
def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()

    # parser.add_option("-i","--inputfiles",default="/Users/mcoughlin/Downloads/lcogtdata-20200118-18/g/*fits.fz",type="string") # didn't work, it only loads 1 file
    parser.add_option("-i","--inputfolder",default="/Users/mcoughlin/Downloads/lcogtdata-20200118-18/g/",type="string") #added to get all files
    parser.add_option("-t","--telescope",default="1m",type="string") # added to get all files , use 'm' to get all data from all telescopes
    parser.add_option("-o","--outputfile",default="/Users/mcoughlin/Downloads/lcogtdata-20200118-18/g/stack_g_91.fits") # need the '91' to run photometry later!
    

    parser.add_option("-r","--ra",default=62.9711,type=float)
    parser.add_option("-d","--declination",default=43.7496,type=float)

    parser.add_option("--defaultsDir",default="../defaults")

    parser.add_option("--doAstrometryNet",  action="store_true", default=False)
    parser.add_option("--doOverwrite",  action="store_true", default=False)
    parser.add_option("--doMedian",  action="store_true", default=False)
    parser.add_option("--doSextractor",  action="store_true", default=False)

    opts, args = parser.parse_args()

    return opts

opts = parse_commandline()
# inputfiles = opts.inputfiles
inputfolder = opts.inputfolder

print(inputfolder+'*'+opts.telescope+'*fits*')
outputfile = opts.outputfile
# fitsfiles = sorted(glob.glob(inputfiles))
fitsfiles = sorted(glob.glob(inputfolder+'*'+opts.telescope+'*fits*'))
defaultsDir = opts.defaultsDir

if opts.doOverwrite:
    rm_command = "rm -rf %s"%outputfile
    os.system(rm_command)

if not os.path.isfile(outputfile):
    hdulist2 = []
    cnt = 1
    median_size = 40
    for jj in range(len(fitsfiles)):
        print(jj,fitsfiles[jj])
        fitsfile = fitsfiles[jj]
        hdulist = fits.open(fitsfile)
        data_tmp = hdulist[1].data

        if cnt == 1:
            hdulist_hold = copy.copy(hdulist[1])
            xshape, yshape = hdulist_hold.data.shape
            if opts.doMedian:
                data = np.empty([xshape,yshape,0])
            else:
                data = np.zeros([xshape,yshape])
                reference = hdulist[1].data
        shift, error, diffphase = register_translation(reference, data_tmp, upsample_factor=1)
        shifted = image_registration.fft_tools.shiftnd(data_tmp, (shift[0], shift[1]))
        if opts.doMedian:
            data = np.append(data,np.expand_dims(shifted,axis=2),axis=2)
        else:
            data = data + shifted
        cnt = cnt + 1            
        
        if opts.doMedian:
            hdulist_hold.data = np.median(data,axis=2)
        else:
            hdulist_hold.data = data
        #hdulist_hold.data = np.mean(data,axis=2)
        hdulist2 = fits.HDUList(hdus=hdulist_hold)
        hdulist2.writeto(outputfile,output_verify='fix',overwrite=True)
else:
    print('Using existing %s... run --doOverwrite if unwanted.' % outputfile)

if opts.doAstrometryNet:
    if opts.doSextractor:
        hdulist = fits.open(outputfile)
        xsize,ysize = hdulist[1].data.shape

        catfile = outputfile.replace(".fits",".cat")
        ztfsub.utils.sextractor(outputfile, defaultsDir,
                                doSubtractBackground = False,
                                catfile = catfile)

        cat = np.loadtxt(catfile)
        psfs = cat[:,20]
        psfthresh = 3.0/0.26
        idx = np.where(cat[:,20] <= psfthresh)[0]
        cat = cat[idx,:]

        magthresh = np.min(cat[:,4])+8.0
        idx = np.where(cat[:,4] <= magthresh)[0]
        cat = cat[idx,:]
        np.savetxt(catfile,cat,fmt='%.5f')

        index_xyls = catfile.replace(".cat",".cat.fits")
        tbl = Table([cat[:,0],cat[:,1],cat[:,4]], names=('XIMAGE','YIMAGE','MAG'))
        tbl.write(index_xyls, format='fits',overwrite=True)

        ztfsub.utils.astrometrynet(outputfile,pixel_scale=0.26,ra=opts.ra,dec=opts.declination,radius=5.0,depth=100,index_xyls=index_xyls)
        xs, ys = cat[:,0],cat[:,1]

    else:
        ztfsub.utils.astrometrynet(outputfile,pixel_scale=0.26,ra=opts.ra,dec=opts.declination,radius=5.0,depth=100,ext=1)    

        hdulist = fits.open(outputfile.replace(".fits",".axy"))
        data_out = hdulist[1].data
        xs, ys = [], []
        for data_tmp in data_out:
            xs.append(data_tmp[0])
            ys.append(data_tmp[1])
        xs, ys = np.array(xs), np.array(ys)

    plotName = outputfile.replace(".fits",".png")
    fig = plt.figure(figsize=(12,10))
    f1 = aplpy.FITSFigure(outputfile,figure=fig)
    #f1.set_tick_labels_font(size='xx-large')
    #f1.set_axis_labels_font(size='xx-large')
    f1.show_grayscale(invert=False)
    f1.show_circles(xs,ys,5,coords_frame='pixel',zorder=99,edgecolor='red')
    #f1.axis_labels.set_xtext('Right Ascension')
    #f1.axis_labels.set_ytext('Declination')

    fig.canvas.draw()
    plt.savefig(plotName)
    plt.close()   
 
