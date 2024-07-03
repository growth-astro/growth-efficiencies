"""Run Tophat model for a number of absolute magnitudes and plot the results."""

import numpy as np
import subprocess
import argparse
import sys
import os
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--skymap", help="GW skymap.", default='../data/GW190425/LALInference.fits.gz')
    parser.add_argument("-o", "--outputDir", help="output file",default="../output/GW190425/simsurvey/tophat/")
    parser.add_argument("-n", "--ntransient", help="int number", default=10000, type=int)
    parser.add_argument("--minmag", help="minimum absolute magnitude", default=-10.0, type=float)
    parser.add_argument("--maxmag", help="maximum absolute magnitude", default=-20.0, type=float)
    parser.add_argument("--doPlotEfficiency",  action="store_true", default=False)
    parser.add_argument("--outputFile", help="file to store terminal output", default=sys.stdout, type=str)
    parser.add_argument("--nbins", help="number of luminosity bins", default=50, type=int)
    opts = parser.parse_args()
    return opts

def simsurvey_outputtxt(txtfile, mag):
    knsim = {}
    f = open(txtfile,'r')
    i = 0
    mag[i] = round(mag[i],1)
    knsim[mag[i]] = {}
    for line in f.readlines():
        if '1 exposure' in line: probexp1 = float(line.strip(' ').split(':')[1])
        elif '2 exposure' in line: probexp2 = float(line.strip(' ').split(':')[1])
        elif '3+ exposure' in line: probexp3 = float(line.strip(' ').split(':')[1])
        elif 'Passing kNe' in line: knsim[mag[i]]['detected'] = int(line.strip(' ').split(':')[1])
        elif 'Number of created kNe falling in the covered area' in line: 
            knsim[mag[i]]['observed'] = int(line.strip(' ').split(':')[1])
            knsim[mag[i]]['det_efficiency'] = float(knsim[mag[i]]['detected'] / opts.ntransient)
            knsim[mag[i]]['obs_efficiency'] = float(knsim[mag[i]]['detected'] / knsim[mag[i]]['observed'])
            i += 1
            if i < len(mag):
                mag[i] = round(mag[i],1)
                knsim[mag[i]] = {}
        else: pass
    f.close()
    return knsim, [probexp1, probexp2, probexp3] 


opts = parse_commandline()
absmag = np.linspace(opts.maxmag,opts.minmag, opts.nbins)
outputname = os.path.join(opts.outputDir, opts.outputFile)
output = open(outputname, 'w+')

for mag in absmag:
    process = subprocess.run(['python', 'kilonova_sim.py', '--skymap', str(opts.skymap), '--observation_type', 'obsfile' \
        '--outputDir', str(opts.outputDir), '--ntransient', str(opts.ntransient), '--mag', str(mag), '--doSimSurvey', \
        '--pickleFile', 'm'+str(np.abs(mag))+'_sim.pkl', '1>&2'], stdout=subprocess.PIPE)
    rows = str(process.stdout).split('\\n')
    for row in rows: output.write(row+'\n')

output.close()

if opts.doPlotEfficiency:
    figname = os.path.join(str(opts.outputDir), 'efficiency_tophat.pdf')
    knsim, probs = simsurvey_outputtxt(outputname, absmag)
    plt.figure()
    for mag in knsim.keys():
        plt.scatter(mag, knsim[mag]['det_efficiency'], color='m', label=None)
        plt.scatter(mag, knsim[mag]['obs_efficiency'], color='goldenrod', label=None)
    mline = Line2D([],[],color='m', marker='o',label='overall efficiency')
    gline = Line2D([],[],color='goldenrod', marker='o', label='efficiency within observed area')
    plt.xlim(opts.minmag+0.5,opts.maxmag-0.5)
    plt.ylim(0.0, 1.0)
    plt.xlabel('absolute magnitude')
    plt.ylabel('Fraction of KNe detected')
    plt.legend(handles=[mline, gline])
    plt.grid()
    plt.savefig(figname)
