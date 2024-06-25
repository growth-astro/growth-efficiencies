
import matplotlib
#matplotlib.rc('text', usetex=True)
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16, 'legend.fontsize': 10})
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import matplotlib.pyplot as plt
from matplotlib import rcParams

import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.table import Table
from astropy.io import fits, ascii
from astropy.time import Time

#Below are all the needed functions for smoothing the spectrum.
def autorange(start, end, increment):
    while start <=end:
        yield start
        start +=increment

def endpt(data, pointA, pointB):
    mx=np.size(data)
    endpt_smooth=np.ones(mx)
    
    for x in range(pointA, pointB):
        if x in range(0,2):
            endpt_smooth[x]=data[x]
        elif x in range(mx-1,mx+1):
            endpt_smooth[x]=data[x]
        else:
            A=data[x-2:x+2]
            average5=np.mean(A)
            endpt_smooth[x]=average5
    return endpt_smooth

def smooth2(data, start, end, increment):
    first_smooth=np.ones(np.size(data))
    endpoint=endpt(data,0,start)
    ep1=endpoint[0:start]
    first_smooth[0:start]=ep1
    
    for x in autorange(start, end, increment):
        #data[np.isnan(data)]=1
        A=data[x-1:x+1]
        average3=np.mean(A)
        first_smooth[x]=average3
        
    endpoint2=endpt(data,end,np.size(data))
    ep2=endpoint2[end:np.size(data)]
    first_smooth[end:np.size(data)]=ep2
    return first_smooth

def repeat_smooth(data, start, end, increment, repeat):
    for y in range(0,repeat+1):
        if y==0:
            data1=data
        else:
            data1=spectrum_smooth
            start=start+2
            end=end-2
        spectrum_smooth=smooth2(data1, start, end, increment)
    return spectrum_smooth

def normalize(data, repeat_smooth, xmin, xmax):
    normalized_spectrum=data[xmin:xmax]/repeat_smooth[xmin:xmax]
    return normalized_spectrum

def mask_nan(data):
    mask=np.ones(np.size(data))
    nanvals=np.isnan(data)
    for x in range(0,np.size(data)):
        if nanvals[x]==True: 
            mask[x]=-9999
        else:
            mask[x]=data[x]
    return mask

path = '../data/GW190425/spectra/'
files = np.array(os.listdir(path))
candidates = ['AT2019ebq', 'ZTF19aarykkb', 'ZTF19aasckkq', 'ZTF19aarzaod']
classifications = ['Ib/c', 'II', 'IIb' ,'II']
PS1_transient = 'AT2019ebq'
file_cands = [f[0:12] for f in files]
event_date = Time('2019-04-25')

# smoothing the spectra for each of the candidates
filenames = []
redshift = [0.037, 0.024, 0.0528, 0.028]
spectra = {}
for i in range(len(files)):
    data = Table(ascii.read(path+files[i]))
    candidate, date, telescope, v = files[i].replace(".ascii","").split('_')
    key = "%s_%s" % (telescope, v)
    j = [j for j in range(len(candidates)) if candidates[j] == candidate][0]
    wavelength, flux = np.array(data['col1']), np.array(data['col2'])
    wavelength = wavelength[~np.isnan(flux)]
    flux = flux[~np.isnan(flux)]
    # smoothed_flux = repeat_smooth(flux, 10, len(flux)-1, 1, 0)
    if not candidate in spectra:
        spectra[candidate] = {}
    if not key in spectra[candidate]:
        spectra[candidate][key] = {}
    spectra[candidate][key]['date'] = date[0:4]+'-'+date[4:6]+'-'+date[6:8] 
    spectra[candidate][key]['telescope'] = telescope 
    spectra[candidate][key]['wavelength'] = wavelength
    spectra[candidate][key]['flux'] = flux
    spectra[candidate][key]['classification'] = classifications[j]
    spectra[candidate][key]['redshift'] = redshift[j]

candidates = ['AT2019ebq', 'ZTF19aarykkb', 'ZTF19aasckkq', 'ZTF19aarzaod']
color = ['indigo', 'b', 'k', 'goldenrod']
# color = ['k', 'darkgray', 'gray', 'k']
# linestyle = ['-', ':', '--', '-']
offset = -1.6
ii = 0    
plt.figure(1)
plt.rc('font', size=15)
fig, ax = plt.subplots(2, 1, sharey='row', gridspec_kw={'height_ratios': [3,1]})
# He I Lines:
ax[0].axvline(3790, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(3740, 5.45 + offset, 'He I', fontsize=10, rotation='vertical')
ax[0].axvline(5700, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(5650, 5.45 + offset,'He I', fontsize=10, rotation='vertical')
ax[0].axvline(5930, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(5880, 5.45 + offset,'He I', fontsize=10, rotation='vertical')
ax[0].axvline(6700, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(6650, 5.45 + offset,'He I', fontsize=10, rotation='vertical')
ax[0].axvline(7065, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(7015, 5.45 + offset,'He I', fontsize=10, rotation='vertical')


# H I lines
ax[0].axvline(3970.1, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(3920.1, 5.35 + offset, 'H$\epsilon$', fontsize=10, rotation='vertical')
ax[0].axvline(4101.7, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(4060, 5.35 + offset,'H$\delta$', fontsize=10, rotation='vertical')
ax[0].axvline(4340.5, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(4290.5, 5.35 + offset, 'H$\gamma$', fontsize=10, rotation='vertical')
ax[0].axvline(4861.3, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(4811, 5.35 + offset,r'H$\beta$', fontsize=10, rotation='vertical')
ax[0].axvline(6562.8, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(6510, 5.4 + offset,r'H$\alpha$', fontsize=10, rotation='vertical')

# Telluric
ax[0].axvline(7430, color='lightgray', linestyle='--', linewidth=0.5)
ax[0].text(7380, 5.65 + offset,'Telluric', fontsize=10, rotation='vertical')


# marking important features on AT2019ebq
ax[1].axvline(10330, color='lightgray', linestyle='--', linewidth=0.5)
ax[1].axvline(10897, color='lightgray', linestyle='--', linewidth=0.5)
ax[1].axvline(14079, color='lightgray', linestyle='--', linewidth=0.5)
ax[1].text(10330, 0.32,'He I', fontsize=10, rotation='vertical')
ax[1].text(10897, 0.32, 'O I', fontsize=10, rotation='vertical')
ax[1].text(14079, 0.32,'Mg I', fontsize=10, rotation='vertical')

for i, candidate in enumerate(candidates):
    candidate_dict = spectra[candidate]
    if candidate == 'ZTF19aarykkb': key = list(candidate_dict.keys())[3]
    else: key = list(candidate_dict.keys())[-1]
    classification = candidate_dict[key]['classification']
    if not "wavelength" in candidate_dict[key]: continue
    wavelength, flux =  candidate_dict[key]['wavelength'] / (1 + candidate_dict[key]['redshift']), \
        candidate_dict[key]['flux']
    tobs = Time(candidate_dict[key]['date']) - event_date
    wavelength = wavelength[~np.isnan(flux)]
    flux = flux[~np.isnan(flux)]
    label = "%s %s" % (candidate, candidate_dict[key]['telescope'])
    if candidate in ['AT2019ebq']: 
        ax[1].plot(wavelength, (flux / np.amax(flux)), 'k', label=label, linewidth=0.75)
    else:
        ax[0].plot(wavelength, ii + (flux / np.amax(flux)), 'k', label=label, linewidth=0.75)
    if candidate in ['ZTF19aasckkq']: ii = ii + 0.5
    else: ii = ii + 1.0 

ax[0].text(6750, 1.75, 'ZTF19aarykkb DCT', fontsize=10) 
ax[0].text(6750, 2.5, 'ZTF19aasckkq P200', fontsize=10) 
ax[0].text(6750, 3.4, 'ZTF19aarzaod SALT', fontsize=10) 
ax[1].text(18000, 0.20, 'AT2019ebq/PS19qp', fontsize=10)
ax[1].text(20500, 0.12, 'Keck2', fontsize=10)
fig.text(0.5, 0.0, 'Rest Wavelength [$\AA$]', ha='center')
fig.text(0.04, 0.5, '$F_\lambda$ [erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$] + constant', va='center', rotation='vertical')

ax[0].set_xlim(3000, 8500)
ax[1].set_xlim(9500, 22000)
# ax[0].set_ylim(0.0, 0.4)
ax[1].set_ylim(0.0, 0.4)
ax[0].yaxis.set_major_locator(plt.NullLocator())
ax[1].yaxis.set_major_locator(plt.NullLocator())
ax[0].xaxis.set_minor_locator(plt.MultipleLocator(250))
ax[1].xaxis.set_minor_locator(plt.MultipleLocator(500))
# ax[0].grid(color='k')
# ax[1].grid(color='k')
# ax[0].legend(loc='upper left')
# ax[1].legend(loc='upper right')
plt.savefig('../output/GW190425/spectra.pdf', bbox_inches='tight')
plt.close()

