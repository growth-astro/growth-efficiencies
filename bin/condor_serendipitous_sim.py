
import os, sys
import glob
import optparse
from astropy.time import Time

import tables
import pandas as pd
import numpy as np
import h5py

import ztfperiodic.utils

try:
    from penquins import Kowalski
except:
    print("penquins not installed... need to use matchfiles.")

np.random.seed(0)

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()

    parser.add_option("-p","--python",default="python")

    parser.add_option("-i", "--inputDir", help="input directory",default="/home/michael.coughlin/ZTF/growth-too-papers/input/")
    parser.add_option("-o", "--outputDir", help="output directory",default="/home/michael.coughlin/ZTF/growth-too-papers/output/serendipitous/GW190426/")
    parser.add_option("--possisDir", help="possis directory",default="/home/michael.coughlin/ZTF/growth-too-papers/data/possis/")
    parser.add_option("--secrets", help="secrets filename.", default='/home/michael.coughlin/ZTF/growth-too-papers/bin/secrets.csv')

    parser.add_option("--observations", help="output file",default="/home/michael.coughlin/ZTF/growth-too-papers/data/serendipitous/2019-04-25-2019-10-03_subdata.txt")

    parser.add_option("-n", "--Nsamples", help="number of samples",default=50,type=int)

    parser.add_option("-m", "--modelType", help="(Bulla, Tophat, afterglow)",default="Tophat")

    parser.add_option("-t", "--starttime", help="time",default="2019-05-01T00:00:00")
    parser.add_option("-e", "--endtime", help="time",default="2019-07-01T00:00:00")

    parser.add_option("--doSerendipitousSim",  action="store_true", default=False)
    parser.add_option("--thetadist", help="(sine, uniform)",default="uniform")
    parser.add_option("--thetai", help="Initial theta", default=0, type=float)
    parser.add_option("--thetaf", help="Final theta", default=90, type=float)

    parser.add_option("--observation_type", help="observation type.", default='obsfile')

    opts, args = parser.parse_args()

    return opts

# Parse command line
opts = parse_commandline()
baseoutputDir = opts.outputDir
Nsamples = opts.Nsamples

if opts.modelType == "Bulladynwind":
    baseoutputDir = os.path.join(baseoutputDir,"%d_%d"%(opts.thetai, opts.thetaf))

condorDir = os.path.join(baseoutputDir,'condor')
if not os.path.isdir(condorDir):
    os.makedirs(condorDir)

logDir = os.path.join(condorDir,'logs')
if not os.path.isdir(logDir):
    os.makedirs(logDir)

job_number = 0 
job_number_query = 0

dir_path = os.path.dirname(os.path.realpath(__file__))

if Nsamples == 1:
    jd_min = Time(opts.starttime, format='isot', scale='utc').jd
    jd_max = Time(opts.endtime, format='isot', scale='utc').jd
    ttstart, ttend = Time(jd_min, format='jd').isot, Time(jd_max, format='jd').isot
else:
    jd_min = Time('2019-10-05T00:00:00', format='isot', scale='utc').jd
    jd_max = Time('2019-10-08T00:00:00', format='isot', scale='utc').jd
phi_min, phi_max = -90.0, 90.0
theta_min, theta_max = 0.0, 360.0
ntransient = 1000
NPoints = 1
mag_min, mag_max = -20, -12
mags = np.linspace(mag_min, mag_max, 31)
dmags = np.arange(-0.05,1.10,0.05)

E0_min, E0_max = 1e45, 1e53
E0s = np.logspace(np.log10(E0_min), np.log10(E0_max), 51)
#phis = [0, 15, 30, 45, 60, 75, 90]
phis = [15, 30]
mejwinds = [0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13]
mejdyns = [0.001, 0.005, 0.01, 0.02]

thetai, thetaf = opts.thetai, opts.thetaf

if opts.modelType == "afterglow":
    Nmag = len(E0s)
elif opts.modelType == "Tophat":
    Nmag = len(mags)*len(dmags)
elif opts.modelType == "Bulladynwind":
    Nmag  = len(mejwinds)*len(mejdyns)*len(phis)

if opts.doSerendipitousSim:
    condordag = os.path.join(condorDir,'condor.dag')
    fid = open(condordag,'w') 
    condorsh = os.path.join(condorDir,'condor.sh')
    fid1 = open(condorsh,'w') 

for ii in range(Nsamples):
    for jj in range(Nmag):
        if opts.modelType == "afterglow":
            E0 = E0s[jj] 
            outputDir = os.path.join(baseoutputDir, '%.3f_%.3f_%.2f' % (jd_min, jd_max, np.log10(E0)))
        elif opts.modelType == "Tophat":
            mm, dind = np.divmod(jj, len(mags))
            mag = mags[dind]
            dmag = dmags[mm]
            outputDir = os.path.join(baseoutputDir, '%.3f_%.3f_%.2f_%.2f' % (jd_min, jd_max,
                                                                         mag,
                                                                         dmag
                                                                         ))    
        elif opts.modelType == "Bulladynwind":
            mm, rem = np.divmod(jj, len(mejdyns)*len(phis))
            pp, tind = np.divmod(rem, len(mejdyns))
            mejwind = mejwinds[mm]
            phiangle = phis[pp]
            mejdyn = mejdyns[tind]
            outputDir = os.path.join(baseoutputDir, '%.3f_%.3f_%.3f_%.3f_%d' % (
                                                    jd_min, jd_max, mejwind, mejdyn,
                                                    phiangle))    

        if opts.doSerendipitousSim:
            if opts.modelType == "afterglow":
                fid1.write('%s %s/serendipitous_sim.py --doPlots --doFilter --doSimSurvey --doRotate -o %s --start_time %s --end_time %s --observation_type %s --observations %s --inputDir %s --possisDir %s --ntransient %d --E0 %.5e --NPoints %d --modelType afterglow\n' % (opts.python, dir_path, outputDir, ttstart, ttend, opts.observation_type, opts.observations, opts.inputDir, opts.possisDir,ntransient,E0,NPoints))
            elif opts.modelType == "Tophat":
                fid1.write('%s %s/serendipitous_sim.py --doPlots --doFilter --doSimSurvey --doRotate -o %s --start_time %s --end_time %s --observation_type %s --observations %s --inputDir %s --possisDir %s --ntransient %d --mag %.5f --dmag %.5f --NPoints %d --modelType Tophat\n' % (opts.python, dir_path, outputDir, ttstart, ttend, opts.observation_type, opts.observations, opts.inputDir, opts.possisDir,ntransient,mag,dmag,NPoints)) 
            elif opts.modelType == "Bulladynwind":
                fid1.write('%s %s/serendipitous_sim.py --doPlots --doFilter --doSimSurvey --doRotate -o %s -s %s --start_time %s -end_time %s --observation_type %s --observations %s --inputDir %s --possisDir %s --ntransient %d --mwin %.3f --opening_angle %d --mdyn %.3f --NPoints %d --modelType Bulladynwind --thetai %.3f --thetaf %.3f\n' % (opts.python, dir_path, outputDir, ttstart, ttend, opts.observation_type, opts.observations, opts.inputDir, opts.possisDir,ntransient,mejwind,phiangle,mejdyn,NPoints,thetai,thetaf))

            fid.write('JOB %d condor.sub\n'%(job_number))
            fid.write('RETRY %d 3\n'%(job_number))
            if opts.modelType == "afterglow":
                fid.write('VARS %d jobNumber="%d" outputDir="%s" start_time="%s" end_time="%s" observations="%s" inputDir="%s" possisDir="%s" E0="%.5e"\n'%(job_number,job_number,outputDir,ttstart,ttend,opts.observations,opts.inputDir,opts.possisDir,E0))
            elif opts.modelType == "Tophat":
                fid.write('VARS %d jobNumber="%d" outputDir="%s" start_time="%s" end_time="%s" observations="%s" inputDir="%s" possisDir="%s" mag="%.5f" dmag="%.5f"\n'%(job_number,job_number,outputDir,ttstart,ttend,opts.observations,opts.inputDir,opts.possisDir,mag,dmag))
            elif opts.modelType == "Bulladynwind":
                fid.write('VARS %d jobNumber="%d" outputDir="%s" start_time="%s" end_time="%s" observations="%s" inputDir="%s" possisDir="%s" mwin="%.3f" opening_angle="%d" mdyn="%.3f" thetai="%.3f" thetaf="%.3f"\n'%(job_number,job_number,outputDir,ttstart,ttend,opts.observations,opts.inputDir,opts.possisDir,mejwind,phiangle,mejdyn,thetai,thetaf))
            fid.write('\n\n')
            job_number = job_number + 1
       
if opts.doSerendipitousSim:
    fid1.close()
    fid.close()

if opts.doSerendipitousSim:    
    fid = open(os.path.join(condorDir,'condor.sub'),'w')
    fid.write('executable = %s/serendipitous_sim.py\n'%dir_path)
    fid.write('output = logs/out.$(jobNumber)\n');
    fid.write('error = logs/err.$(jobNumber)\n');

    if opts.modelType == "afterglow":
        fid.write('arguments = --doFilter --doSimSurvey -o $(outputDir) --start_time $(start_time) --end_time $(end_time) --observation_type %s --observations $(observations) --inputDir $(inputDir) --possisDir $(possisDir) --E0 $(E0) --modelType afterglow\n' % opts.observation_type)
    elif opts.modelType == "Tophat":
        fid.write('arguments = --doFilter --doSimSurvey -o $(outputDir) --start_time $(start_time) --end_time $(end_time) --observation_type %s --observations $(observations) --inputDir $(inputDir) --possisDir $(possisDir) --mag $(mag) --dmag $(dmag) --modelType Tophat --NPoints 1\n' % opts.observation_type)
    elif opts.modelType == "Bulladynwind":
        fid.write('arguments = --doFilter --doSimSurvey -o $(outputDir) --start_time $(start_time) --end_time $(end_time) --observation_type %s --observations $(observations) --inputDir $(inputDir) --possisDir $(possisDir) --mwin $(mwin) --mdyn $(mdyn) --opening_angle $(opening_angle) --modelType Bulladynwind --thetai $(thetai) --thetaf $(thetaf)\n' % opts.observation_type)
    fid.write('requirements = OpSys == "LINUX"\n');
    fid.write('request_memory = 16384\n');
    fid.write('request_cpus = 1\n');
    fid.write('accounting_group = ligo.dev.o2.burst.allsky.stamp\n');
    fid.write('notification = never\n');
    fid.write('getenv = true\n');
    fid.write('log = /local/michael.coughlin/simsurvey.log\n')
    fid.write('+MaxHours = 24\n');
    fid.write('universe = vanilla\n');
    fid.write('queue 1\n');
    fid.close()
    
    
