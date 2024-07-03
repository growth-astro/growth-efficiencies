
import os
import glob
import optparse

import pickle
import numpy as np
from astropy.time import Time

import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 24}
matplotlib.rc('xtick', labelsize=22)
matplotlib.rc('ytick', labelsize=22)
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()

    parser.add_option("-p","--python",default="python")

    parser.add_option("-s", "--skymap", help="GW skymap.", default='/home/michael.coughlin/ZTF/growth-too-papers/data/GW190426/LALInference1.fits.gz')
    parser.add_option("-i", "--inputDir", help="input directory",default="/home/michael.coughlin/ZTF/growth-too-papers/input/")
    parser.add_option("-o", "--outputDir", help="output directory",default="/home/michael.coughlin/ZTF/growth-too-papers/output/serendipitous/GW190426/")
    parser.add_option("--possisDir", help="possis directory",default="/home/michael.coughlin/ZTF/growth-too-papers/data/possis/")
    parser.add_option("--secrets", help="secrets filename.", default='/home/michael.coughlin/ZTF/growth-too-papers/bin/secrets.csv')

    parser.add_option("--observations", help="output file",default="/home/michael.coughlin/ZTF/growth-too-papers/data/serendipitous/2019-10-05-2019-10-11_normal.txt")

    parser.add_option("-m", "--modelType", help="(Bulla, Tophat, afterglow)",default="Tophat")

    parser.add_option("--thetadist", help="(sine, uniform)",default="uniform")
    parser.add_option("--thetai", help="Initial theta", default=0, type=float)
    parser.add_option("--thetaf", help="Final theta", default=90, type=float)
    parser.add_option("--ntransient", help="Final theta", default=1000, type=float)

    opts, args = parser.parse_args()

    return opts

# Parse command line
opts = parse_commandline()
baseoutputDir = opts.outputDir
skymap = opts.skymap

if opts.modelType == "Bulladynwind":
    baseoutputDir = os.path.join(baseoutputDir,"%d_%d"%(opts.thetai, opts.thetaf))

data = np.empty((0,3), float)
mags, dmags = [], []
if opts.modelType == "Bulladynwind":
    mejwinds, mejdyns = [], []
ntransients_1_19 = []
ntransients_1_17 = []
ntransients_1_18 = []
ntransients_2 = []
sim = np.empty(0, float)
if opts.modelType == "Bulladynwind":
    outputDirs = glob.glob(os.path.join(baseoutputDir, '*_*_*_*_*'))
elif opts.modelType in ["Tophat","ntransients"]:
    outputDirs = glob.glob(os.path.join(baseoutputDir, '*_*_*_*_*'))
else:
    outputDirs = glob.glob(os.path.join(baseoutputDir, '*_*_*_*'))

for outputDir in outputDirs:
    outsplit = list(filter(None,outputDir.split("/")[-1].split("_")))
    if opts.modelType == "Bulladynwind":
        jd, mejwind, mejdyn, phiangle, theta, phi = np.array(outsplit, dtype=float)
    elif opts.modelType in ["Tophat","ntransients"]:
        jd, mag, dmag, phi, theta = np.array(outsplit, dtype=float)
    else:
        jd, mag, phi, theta = np.array(outsplit, dtype=float)
    filename = os.path.join(outputDir, 'prob.dat')
    if os.path.isfile(filename):
        data_out = np.loadtxt(filename)
        if data_out.size == 0: continue
        data = np.append(data, np.atleast_2d(data_out), axis=0)
    filename = os.path.join(outputDir, 'transients_1_19.dat')
    if os.path.isfile(filename):
        lines = [line.rstrip('\n') for line in open(filename)]
        line = lines[0].replace(" ","").replace("set()","").replace("'","").replace("{","").replace("}","")
        transients = list(filter(None,line.split(",")))
        ntransients_1_19.append(len(transients))
    filename = os.path.join(outputDir, 'transients_1_17.dat')
    if os.path.isfile(filename):
        lines = [line.rstrip('\n') for line in open(filename)]
        line = lines[0].replace(" ","").replace("set()","").replace("'","").replace("{","").replace("}","")
        transients = list(filter(None,line.split(",")))
        ntransients_1_17.append(len(transients))
    filename = os.path.join(outputDir, 'transients_1_18.dat')
    if os.path.isfile(filename):
        lines = [line.rstrip('\n') for line in open(filename)]
        line = lines[0].replace(" ","").replace("set()","").replace("'","").replace("{","").replace("}","")
        transients = list(filter(None,line.split(",")))
        ntransients_1_18.append(len(transients))
    filename = os.path.join(outputDir, 'transients_2.dat')
    if os.path.isfile(filename):
        lines = [line.rstrip('\n') for line in open(filename)]
        line = lines[0].replace(" ","").replace("set()","").replace("'","").replace("{","").replace("}","")
        transients = list(filter(None,line.split(",")))
        ntransients_2.append(len(transients))
    filename = os.path.join(outputDir, 'sim.pkl')
    if os.path.isfile(filename):
        f1 = open(outputDir+'/sim.pkl', 'rb')
        lcs = pickle.load(f1)
        f1.close()               
        data_out = len(lcs.lcs)/opts.ntransient
        sim = np.append(sim, data_out)
        if opts.modelType == "Bulladynwind":
            mejwinds.append(mejwind)
            mejdyns.append(mejdyn)
        else:
            mags.append(mag)
        if opts.modelType == "Tophat":
            dmags.append(dmag)

outputDir = os.path.join(baseoutputDir, 'combined')
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)
 
color2 = 'coral'
color1 = 'cornflowerblue'
color3 = 'palegreen'
color4 = 'darkmagenta'

if opts.modelType == "Bulladynwind":
    mejwinds = np.array(mejwinds)
    mejdyns = np.array(mejdyns)
else:
    mags = np.array(mags)
    if opts.modelType == "afterglow":
        mags = 10**mags
    mags_unique = np.unique(mags)

if opts.modelType == "Tophat":
    models = [['SN Ia', -0.200, -0.186, -19.00, -16.61],
              ['SN IIP', -0.157, -0.166, -16.84, -17.03],
              ['SN IIb', -0.169, -0.261, -17.44, -17.32],
              #['SN IIn', 0.105, 0.115, -17.58, -15.32],
              ['SN Ib', -0.142, -0.146, -17.4, -17.25],
              ['SN Ic', -0.120, -0.143, -17.48, -14.92],
              ['GW170817', 0.5, 0.5, -16.0, -16.0]]

    dmags_unique = np.unique(dmags)
    X, Y = np.meshgrid(mags_unique,dmags_unique)
    Z = np.zeros(X.shape)
    for ii, mag in enumerate(mags_unique):
        for jj, dmag in enumerate(dmags_unique):
            if dmag < 0:
                mag_zero = mag - 3*dmag
                diff = np.abs(mags - mag_zero)
                mindiff = np.amin(diff)
                idx = np.where((dmag == dmags) & (diff == mindiff))[0]
                #idx = np.where((mag_zero == mags) & (dmag == dmags))[0]
                if len(idx) == 0: continue
                if mag_zero > -10.0: continue
            else:
                idx = np.where((mag == mags) & (dmag == dmags))[0]
                if len(idx) == 0: continue

            Z[jj, ii] = np.median(sim[idx])

    diffs = [[0.44,-0.15], [0.4,0.25],
             [0.47,0.32], #[0.5,-0.25],
             [0.4,-0.2], [0.45, 0.3],
             [0.45, 0.12]]

    Z[Z==0.0] = 1e-3
    plotName = os.path.join(outputDir,'tophat.pdf')
    fig = plt.figure(figsize=(8,6))
    c = plt.pcolor(X,Y,Z,vmin=1e-2,vmax=np.max(Z),cmap='coolwarm')
    cbar = fig.colorbar(c)
    plt.xlabel('Absolute Magnitude')
    plt.ylabel('Change in magnitude per day')
    for model, diff in zip(models,diffs):
        #plt.text(model[3], model[1], model[0])
        plt.plot(model[3], model[1], marker='o', markerfacecolor='none',
                    markeredgecolor='white', markersize=30)

        diffx, diffy = diff
        x, y = model[3], model[1]

        x1text = x*diffx*2.0
        y1text = y*diffy*15.0

        x1arrow = x*diffx*2.2
        y1arrow = y*diffy*12.0

        plt.text(x1text, y1text, model[0])
        plt.annotate('', xy=(x,y), xytext=(x1arrow,y1arrow),
                     arrowprops=dict(facecolor='black', arrowstyle='->'),
                    )

    cbar.set_label('Efficiency')
    plt.ylim([-1,1])
    plt.xlim([-20,-12])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.savefig(plotName,bbox_inches='tight')
    plt.close()


elif opts.modelType == "Bulladynwind":
    mejwinds_unique = np.unique(mejwinds)
    mejwinds_unique = np.append(mejwinds_unique,1)
    mejdyns_unique = np.unique(mejdyns)
    mejdyns_unique = np.append(mejdyns_unique,1)
    X, Y = np.meshgrid(mejwinds_unique,mejdyns_unique)
    Z = np.zeros(X.shape)
    for ii, mejwind in enumerate(mejwinds_unique):
        for jj, mejdyn in enumerate(mejdyns_unique):
            idx = np.where((mejwind == mejwinds) & (mejdyn == mejdyns))[0]
            if len(idx) == 0: continue
            Z[jj,ii] = np.median(sim[idx])
            print(mag, dmag, sim[idx])

    Z[Z==0.0] = 1e-3
    plotName = os.path.join(outputDir,'Bulla.pdf')
    fig = plt.figure(figsize=(8,6))
    c = plt.pcolor(X,Y,Z,norm=LogNorm(),vmin=1e-3,vmax=np.max(Z),cmap='coolwarm')
    cbar = fig.colorbar(c)
    plt.xlabel('Wind ejecta [solar mass]')
    plt.ylabel('Dynamical ejecta [solar mass]')
    cbar.set_label('Efficiency')
    plt.xlim([0.01, 0.13])
    plt.ylim([0.001, 0.02])
    plt.savefig(plotName,bbox_inches='tight')
    plt.close()

plotName = os.path.join(outputDir,'sim.pdf')
fig = plt.figure(figsize=(19,12))
ax = fig.add_subplot(1, 1, 1)
labels = []
for mag in mags_unique:
    idx = np.where(mag == mags)[0]
    magsi = sim[idx]
    mag_10, mag_50, mag_90 = np.percentile(magsi,10), np.percentile(magsi,50), np.percentile(magsi,90)
    #asymmetric_error = np.atleast_2d([mag_50-mag_10, mag_90-mag_50]).T
    #plt.errorbar(mag, mag_50, yerr=asymmetric_error, c='r',zorder=1,fmt='o')

    if opts.modelType == "Bulla":
        widths=0.01
    elif opts.modelType == "Tophat":
        widths=0.05
    elif opts.modelType == "afterglow":
        widths=mag*0.1

    parts = plt.violinplot(magsi.T,[mag],widths=widths)
    for partname in ('cbars','cmins','cmaxes'):
        vp = parts[partname]
        vp.set_edgecolor(color1)
        vp.set_linewidth(1)
    for pc in parts['bodies']:
        pc.set_facecolor(color1)
        pc.set_edgecolor(color1)

    plt.plot([mag-widths,mag+widths],[mag_10,mag_10],'--',color=color1)
    plt.plot([mag-widths,mag+widths],[mag_90,mag_90],'--',color=color1)

    labels.append(mag)

#ax.set_xticks(np.arange(len(mags_unique)))
#ax.set_xticklabels(labels)

if opts.modelType == "afterglow":
    plt.xlabel(r"Isotropic Energy [erg/s]",fontsize=24)
elif opts.modelType == "Tophat":
    plt.xlabel(r"Absolute Magnitude",fontsize=24)
elif opts.modelType == "Bulla":
    plt.xlabel(r"Ejecta mass",fontsize=24)
plt.ylabel('Efficiency',fontsize=24)
#plt.xlim([0,1])
if opts.modelType == "afterglow":
    plt.xlim([1e48,1e50])
    ax.set_xscale('log')
elif opts.modelType == "Tophat":
    plt.xlim([-19.0,-15.0])
    ax.invert_xaxis()
elif opts.modelType == "Bulla":
    ax.set_yscale('log')
#plt.legend()
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.savefig(plotName)
plt.close()

bins = np.arange(-0.05,1.05,0.05)
hist_1, bin_edges_1 = np.histogram(data[:,0],bins)
hist_2, bin_edges_2 = np.histogram(data[:,1],bins)
hist_3, bin_edges_3 = np.histogram(data[:,2],bins)
bins = (bin_edges_1[1:] + bin_edges_1[:-1])/2.0

histcum_1 = np.cumsum(hist_1)/float(np.sum(hist_1))
histcum_2 = np.cumsum(hist_2)/float(np.sum(hist_2))
histcum_3 = np.cumsum(hist_3)/float(np.sum(hist_3))

bins_ntran = np.arange(-0.5,100.5,1)
hist_ntran_1_19, bin_edges_ntran_1_19 = np.histogram(ntransients_1_19, bins_ntran)
hist_ntran_1_17, bin_edges_ntran_1_17 = np.histogram(ntransients_1_17, bins_ntran)
hist_ntran_1_18, bin_edges_ntran_1_18 = np.histogram(ntransients_1_18, bins_ntran)
hist_ntran_2, bin_edges_ntran_2 = np.histogram(ntransients_2, bins_ntran)
bins_ntran = (bin_edges_ntran_2[1:] + bin_edges_ntran_2[:-1])/2.0

histcum_ntran_1_19 = np.cumsum(hist_ntran_1_19)/float(np.sum(hist_ntran_1_19))
histcum_ntran_1_17 = np.cumsum(hist_ntran_1_17)/float(np.sum(hist_ntran_1_17))
histcum_ntran_1_18 = np.cumsum(hist_ntran_1_18)/float(np.sum(hist_ntran_1_18))
histcum_ntran_2 = np.cumsum(hist_ntran_2)/float(np.sum(hist_ntran_2))

plotName = os.path.join(outputDir,'coverage_hist.pdf')
plt.figure(figsize=(19,12))
plt.step(bins, hist_1, 'k-', linewidth=2, label='1', where='mid')
plt.step(bins, hist_2, 'b--',linewidth=2, label='2', where='mid')
plt.step(bins, hist_3, 'g.-',linewidth=2, label='3+', where='mid')
plt.xlabel(r"Skymap Coverage",fontsize=24)
plt.ylabel('Probability Density Function',fontsize=24)
plt.xlim([0,1])
plt.legend()
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.savefig(plotName)
plt.close()

plotName = os.path.join(outputDir,'coverage_cumhist.pdf')
plt.figure(figsize=(19,12))
plt.step(bins, histcum_1, 'k-', linewidth=2, label='1', where='mid')
plt.step(bins, histcum_2, 'b--',linewidth=2, label='2', where='mid')
plt.step(bins, histcum_3, 'g.-',linewidth=2, label='3+', where='mid')
plt.xlabel(r"Skymap Coverage",fontsize=24)
plt.ylabel('Cumulative Density Function',fontsize=24)
plt.xlim([0,1])
plt.ylim([0,1])
plt.legend()
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.savefig(plotName)
plt.close()

plotName = os.path.join(outputDir,'ntransients_hist.pdf')
plt.figure(figsize=(19,12))
plt.step(bins_ntran, hist_ntran_1_19, 'k-', label='1+ mag<19', linewidth=2, where='mid')
plt.step(bins_ntran, hist_ntran_1_17, 'k-', label='1+ mag<17', linewidth=2, where='mid')
plt.step(bins_ntran, hist_ntran_1_18, 'k-', label='1+ mag<18', linewidth=2, where='mid')
plt.step(bins_ntran, hist_ntran_2, 'b--', label='2+', linewidth=2, where='mid')
plt.xlabel(r"Number of Transients",fontsize=24)
plt.ylabel('Probability Density Function',fontsize=24)
plt.legend()
plt.xlim([0,10])
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.savefig(plotName)
plt.close()

plotName = os.path.join(outputDir,'ntransients_cumhist.pdf')
fig = plt.figure(figsize=(19,12))
ax = fig.add_subplot(1, 1, 1)
plt.step(bins_ntran, histcum_ntran_1_19, 'k-', linewidth=2, label='1+ mag<19', where='mid')
plt.step(bins_ntran, histcum_ntran_1_17, 'r:', linewidth=2, label='1+ mag<17', where='mid')
plt.step(bins_ntran, histcum_ntran_1_18, 'g.-', linewidth=2, label='1+ mag<18', where='mid')
plt.step(bins_ntran, histcum_ntran_2, 'b--', linewidth=2, label='2+', where='mid')
plt.xlabel(r"Number of Transients",fontsize=24)
plt.ylabel('Cumulative Density Function',fontsize=24)
#plt.xlim([0,10])
plt.ylim([0,1])
ax.set_xscale('log')
plt.legend()
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.savefig(plotName)
plt.close()

sim = sim * 100 # percentages
bins = np.arange(0.000,1.000,0.001)
hist_1, bin_edges_1 = np.histogram(sim[:],bins)
hist_2, bin_edges_2 = np.histogram(sim[:],bins)
bins = (bin_edges_1[1:] + bin_edges_1[:-1])/2.0

histcum_1 = np.cumsum(hist_1)/float(np.sum(hist_1))
histcum_2 = np.cumsum(hist_2)/float(np.sum(hist_2))

plotName = os.path.join(outputDir,'sim_hist.pdf')
plt.figure(figsize=(19,12))
plt.step(bins, hist_1, 'k-', linewidth=2, label='1', where='mid')
plt.step(bins, hist_2, 'b--',linewidth=2, label='2', where='mid')
plt.xlabel(r"Percentage of KNe",fontsize=24)
plt.ylabel('Probability Density Function',fontsize=24)
#plt.xlim([0,1])
plt.legend()
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.savefig(plotName)
plt.close()

plotName = os.path.join(outputDir,'sim_cumhist.pdf')
plt.figure(figsize=(19,12))
plt.step(bins, histcum_1, 'k-', linewidth=2, label='1', where='mid')
plt.step(bins, histcum_2, 'b--',linewidth=2, label='2', where='mid')
plt.xlabel(r"Percentage of KNe",fontsize=24)
plt.ylabel('Cumulative Density Function',fontsize=24)
#plt.xlim([0,1])
#plt.ylim([0,1])
plt.legend()
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.savefig(plotName)
plt.close()
