
import os, sys
import glob
import optparse
from astropy.time import Time

#import tables
import pandas as pd
import numpy as np
import h5py

#try:
#    from penquins import Kowalski
#except:
#    print("penquins not installed... need to use matchfiles.")

np.random.seed(0)

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
    parser.add_option("--sfdDir", help="SFD directory", default="/home/michael.coughlin/ZTF/growth-too-papers/input/sfd98/")
    parser.add_option("--observations", help="output file",default="/home/michael.coughlin/ZTF/growth-too-papers/data/serendipitous/2019-10-05-2019-10-11_normal.txt")
    parser.add_option("--logDir", help="directory to log file",default="/local/michael.coughlin/")
    parser.add_option("-n", "--Nsamples", help="number of samples",default=50,type=int)

    parser.add_option("-m", "--modelType", help="(Bulladynwind, BulladynwindNSBH, Tophat, afterglow)",default="Tophat")

    parser.add_option("-t", "--starttime", help="time",default="2020-01-05T16:24:26.057208")

    parser.add_option("--doSerendipitousSim",  action="store_true", default=False)
    parser.add_option("--doQuerySkymap",  action="store_true", default=False)

    parser.add_option("--thetadist", help="(sine, uniform)",default="uniform")
    parser.add_option("--thetai", help="Initial theta", default=0, type=float)
    parser.add_option("--thetaf", help="Final theta", default=90, type=float)

    opts, args = parser.parse_args()

    return opts

# Parse command line
opts = parse_commandline()
baseoutputDir = opts.outputDir
Nsamples = opts.Nsamples
skymap = opts.skymap

if opts.modelType == "Bulladynwind" or opts.modelType == "BulladynwindNSBH":
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
    jd_max = jd_min
else:
    jd_min = Time('2019-10-05T00:00:00', format='isot', scale='utc').jd
    jd_max = Time('2019-10-08T00:00:00', format='isot', scale='utc').jd
phi_min, phi_max = -90.0, 90.0
theta_min, theta_max = 0.0, 360.0
ntransient = 10000
NPoints = 1
mag_min, mag_max = -19, -12
mags = np.linspace(mag_min, mag_max, 50)
dmags = np.arange(-0.5,1.6,0.1)

E0_min, E0_max = 1e45, 1e53
E0s = np.logspace(np.log10(E0_min), np.log10(E0_max), 51)
#phis = [0, 15, 30, 45, 60, 75, 90]
phis = [15, 30]
mejwinds = [0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13]
mejdyns = [0.001, 0.005, 0.01, 0.02]

phis_nsbh = [30]
mejwinds_nsbh = np.arange(0.01, 0.1, 0.01)
mejdyns_nsbh = np.arange(0.01, 0.1, 0.01)

thetai, thetaf = opts.thetai, opts.thetaf

if opts.modelType == "afterglow":
    Nmag = len(E0s)
elif opts.modelType == "Tophat":
    Nmag = len(mags)*len(dmags)
elif opts.modelType == "Bulladynwind":
    Nmag  = len(mejwinds)*len(mejdyns)*len(phis)
elif opts.modelType == "BulladynwindNSBH":
    Nmag  = len(mejwinds_nsbh)*len(mejdyns_nsbh)*len(phis_nsbh)

if opts.doSerendipitousSim:
    condordag = os.path.join(condorDir,'condor.dag')
    fid = open(condordag,'w') 
    condorsh = os.path.join(condorDir,'condor.sh')
    fid1 = open(condorsh,'w') 

if opts.doQuerySkymap:
    condordag = os.path.join(condorDir,'condor_query.dag')
    fid2 = open(condordag,'w')
    condorsh = os.path.join(condorDir,'condor_query.sh')
    fid3 = open(condorsh,'w')

for ii in range(Nsamples):
    if Nsamples == 1:
        jd = jd_min
        tt = Time(jd, format='jd', scale='utc').isot
        phi = 0.0
        theta = 0.0
    else:
        jd = np.random.uniform(jd_min, jd_max)
        tt = Time(jd, format='jd', scale='utc').isot
        phi = np.random.uniform(phi_min, phi_max)
        theta = np.random.uniform(theta_min, theta_max)

    for jj in range(Nmag):
        if opts.modelType == "afterglow":
            E0 = E0s[jj] 
            outputDir = os.path.join(baseoutputDir, '%.3f_%.2f_%d_%d' % (jd, np.log10(E0),
                                                                         theta, phi))
        elif opts.modelType == "Tophat":
            mm, dind = np.divmod(jj, len(mags))
            mag = mags[dind]
            dmag = dmags[mm]
            outputDir = os.path.join(baseoutputDir, '%.3f_%.2f_%.2f_%d_%d' % (jd,
                                                                         mag,
                                                                         dmag,
                                                                         theta, phi))    
        elif opts.modelType == "Bulladynwind":
            mm, rem = np.divmod(jj, len(mejdyns)*len(phis))
            pp, tind = np.divmod(rem, len(mejdyns))
            mejwind = mejwinds[mm]
            phiangle = phis[pp]
            mejdyn = mejdyns[tind]
            outputDir = os.path.join(baseoutputDir, '%.3f_%.3f_%.3f_%d_%d_%d' % (
                                                    jd, mejwind, mejdyn,
                                                    phiangle, theta, phi))    
        elif opts.modelType == "BulladynwindNSBH":
            mm, rem = np.divmod(jj, len(mejdyns_nsbh)*len(phis_nsbh))
            pp, tind = np.divmod(rem, len(mejdyns_nsbh))
            mejwind = mejwinds_nsbh[mm]
            phiangle = phis_nsbh[pp]
            mejdyn = mejdyns_nsbh[tind]
            outputDir = os.path.join(baseoutputDir, '%.3f_%.3f_%.3f_%d_%d_%d' % (
                                                    jd, mejwind, mejdyn,
                                                    phiangle, theta, phi))


        if opts.doSerendipitousSim:
            if opts.modelType == "afterglow":
                fid1.write('%s %s/kilonova_sim.py --doPlots --doFilter --doSimSurvey --doRotate -o %s -s %s --start_time %s --theta %.5f --phi %.5f --observation_type obsfile --observations %s --sfdDir %s --inputDir %s --possisDir %s --ntransient %d --E0 %.5e --NPoints %d --modelType afterglow\n' % (opts.python, dir_path, outputDir, skymap, tt, theta, phi, opts.observations, opts.sfdDir, opts.inputDir, opts.possisDir,ntransient,E0,NPoints))
            elif opts.modelType == "Tophat":
                fid1.write('%s %s/kilonova_sim.py --doPlots --doFilter --doSimSurvey --doRotate -o %s -s %s --start_time %s --theta %.5f --phi %.5f --observation_type obsfile --observations %s --sfdDir %s --inputDir %s --possisDir %s --ntransient %d --mag %.5f --dmag %.5f --NPoints %d --modelType Tophat\n' % (opts.python, dir_path, outputDir, skymap, tt, theta, phi, opts.observations, opts.sfdDir, opts.inputDir, opts.possisDir,ntransient,mag,dmag,NPoints)) 
            elif opts.modelType == "Bulladynwind":
                fid1.write('%s %s/kilonova_sim.py --doPlots --doFilter --doSimSurvey --doRotate -o %s -s %s --start_time %s --theta %.5f --phi %.5f --observation_type obsfile --observations %s --sfdDir %s --inputDir %s --possisDir %s --ntransient %d --mwin %.3f --opening_angle %d --mdyn %.3f --NPoints %d --modelType Bulladynwind --thetai %.3f --thetaf %.3f\n' % (opts.python, dir_path, outputDir, skymap, tt, theta, phi, opts.observations, opts.sfdDir, opts.inputDir, opts.possisDir,ntransient,mejwind,phiangle,mejdyn,NPoints,thetai,thetaf))
            elif opts.modelType == "BulladynwindNSBH":
                fid1.write('%s %s/kilonova_sim.py --doPlots --doFilter --doSimSurvey --doRotate -o %s -s %s --start_time %s --theta %.5f --phi %.5f --observation_type obsfile --observations %s --sfdDir %s --inputDir %s --possisDir %s --ntransient %d --mwin %.3f --opening_angle %d --mdyn %.3f --NPoints %d --modelType BulladynwindNSBH --thetai %.3f --thetaf %.3f\n' % (opts.python, dir_path, outputDir, skymap, tt, theta, phi, opts.observations, opts.sfdDir, opts.inputDir, opts.possisDir,ntransient,mejwind,phiangle,mejdyn,NPoints,thetai,thetaf))

            fid.write('JOB %d condor.sub\n'%(job_number))
            fid.write('RETRY %d 3\n'%(job_number))
            if opts.modelType == "afterglow":
                fid.write('VARS %d jobNumber="%d" outputDir="%s" skymap="%s" start_time="%s" theta="%.5f" phi="%.5f" observations="%s" sfdDir="%s" inputDir="%s" possisDir="%s" E0="%.5e"\n'%(job_number,job_number,outputDir,skymap,tt,theta,phi,opts.observations,opts.sfdDir,opts.inputDir,opts.possisDir,E0))
            elif opts.modelType == "Tophat":
                fid.write('VARS %d jobNumber="%d" outputDir="%s" skymap="%s" start_time="%s" theta="%.5f" phi="%.5f" observations="%s" sfdDir="%s" inputDir="%s" possisDir="%s" mag="%.5f" dmag="%.5f"\n'%(job_number,job_number,outputDir,skymap,tt,theta,phi,opts.observations,opts.sfdDir,opts.inputDir,opts.possisDir,mag,dmag))
            elif opts.modelType == "Bulladynwind":
                fid.write('VARS %d jobNumber="%d" outputDir="%s" skymap="%s" start_time="%s" theta="%.5f" phi="%.5f" observations="%s" sfdDir="%s" inputDir="%s" possisDir="%s" mwin="%.3f" opening_angle="%d" mdyn="%.3f" thetai="%.3f" thetaf="%.3f"\n'%(job_number,job_number,outputDir,skymap,tt,theta,phi,opts.observations,opts.sfdDir,opts.inputDir,opts.possisDir,mejwind,phiangle,mejdyn,thetai,thetaf))
            elif opts.modelType == "BulladynwindNSBH":
                fid.write('VARS %d jobNumber="%d" outputDir="%s" skymap="%s" start_time="%s" theta="%.5f" phi="%.5f" observations="%s" sfdDir="%s" inputDir="%s" possisDir="%s" mwin="%.3f" opening_angle="%d" mdyn="%.3f" thetai="%.3f" thetaf="%.3f"\n'%(job_number,job_number,outputDir,skymap,tt,theta,phi,opts.observations,opts.sfdDir,opts.inputDir,opts.possisDir,mejwind,phiangle,mejdyn,thetai,thetaf))
            fid.write('\n\n')
            job_number = job_number + 1
       
    if opts.doQuerySkymap: 
        min_mags = [19, 17, 18]
        for min_mag in min_mags:
            ndethist, min_days = 1, 0
            filename = os.path.join(outputDir, 'transients_1_%d.dat' % min_mag) 
            fid3.write('%s %s/query_skymap.py --skymap %s --level 90 --after-trigger True --jd-trigger %.5f --min-days 0 --within-days 3 --theta %.5f --phi %.5f --ndethist 1 --out %s --min-mag %.1f --secrets %s\n' % (opts.python, dir_path, skymap, jd, theta, phi, filename, min_mag, opts.secrets))
    
            fid2.write('JOB %d condor_query.sub\n'%(job_number_query))
            fid2.write('RETRY %d 3\n'%(job_number_query))
            fid2.write('VARS %d jobNumber="%d" out="%s" skymap="%s" jd="%.5f" theta="%.5f" phi="%.5f" min_mag="%.1f" ndethist="%d" min_days="%.2f"\n'%(job_number_query,job_number_query,filename,skymap,jd,theta,phi,min_mag,ndethist,min_days))
            fid2.write('\n\n')
            job_number_query = job_number_query + 1
    
        ndethist, min_days = 2, 0.02
        filename = os.path.join(outputDir, 'transients_2.dat')
        fid3.write('%s %s/query_skymap.py --skymap %s --level 90 --after-trigger True --jd-trigger %.5f --min-days 0.02 --within-days 3 --theta %.5f --phi %.5f --ndethist 2 --out %s --secrets %s\n' % (opts.python, dir_path, skymap, jd, theta, phi, filename,opts.secrets))
    
        fid2.write('JOB %d condor_query.sub\n'%(job_number_query))
        fid2.write('RETRY %d 3\n'%(job_number_query))
        fid2.write('VARS %d jobNumber="%d" out="%s" skymap="%s" jd="%.5f" theta="%.5f" phi="%.5f" min_mag="%.1f" ndethist="%d" min_days="%.2f"\n'%(job_number_query,job_number_query,filename,skymap,jd,theta,phi,30,ndethist,min_days))
        fid2.write('\n\n')
        job_number_query = job_number_query + 1

if opts.doSerendipitousSim:
    fid1.close()
    fid.close()

if opts.doQuerySkymap:
    fid2.close()
    fid3.close()

if opts.doSerendipitousSim:    
    fid = open(os.path.join(condorDir,'condor.sub'),'w')
    fid.write('executable = %s/kilonova_sim.py\n'%dir_path)
    fid.write('output = logs/out.$(jobNumber)\n');
    fid.write('error = logs/err.$(jobNumber)\n');

    if opts.modelType == "afterglow":
        fid.write('arguments = --doFilter --doSimSurvey --doRotate -o $(outputDir) -s $(skymap) --start_time $(start_time) --theta $(theta) --phi $(phi) --observation_type obsfile --observations $(observations) --sfdDir $(sfdDir) --inputDir $(inputDir) --possisDir $(possisDir) --E0 $(E0) --modelType afterglow\n')
    elif opts.modelType == "Tophat":
        fid.write('arguments = --doFilter --doSimSurvey --doRotate -o $(outputDir) -s $(skymap) --start_time $(start_time) --theta $(theta) --phi $(phi) --observation_type obsfile --observations $(observations) --sfdDir $(sfdDir) --inputDir $(inputDir) --possisDir $(possisDir) --mag $(mag) --dmag $(dmag) --modelType Tophat --NPoints 1\n')
    elif opts.modelType == "Bulladynwind":
        fid.write('arguments = --doFilter --doSimSurvey --doRotate -o $(outputDir) -s $(skymap) --start_time $(start_time) --theta $(theta) --phi $(phi) --observation_type obsfile --observations $(observations) --sfdDir $(sfdDir) --inputDir $(inputDir) --possisDir $(possisDir) --mwin $(mwin) --mdyn $(mdyn) --opening_angle $(opening_angle) --modelType Bulladynwind --thetai $(thetai) --thetaf $(thetaf)\n')
    elif opts.modelType == "BulladynwindNSBH":
        fid.write('arguments = --doFilter --doSimSurvey --doRotate -o $(outputDir) -s $(skymap) --start_time $(start_time) --theta $(theta) --phi $(phi) --observation_type obsfile --observations $(observations) --sfdDir $(sfdDir) --inputDir $(inputDir) --possisDir $(possisDir) --mwin $(mwin) --mdyn $(mdyn) --opening_angle $(opening_angle) --modelType BulladynwindNSBH --thetai $(thetai) --thetaf $(thetaf)\n')
    fid.write('requirements = OpSys == "LINUX"\n');
    fid.write('request_memory = 16384\n');
    fid.write('request_cpus = 1\n');
    fid.write('accounting_group = ligo.dev.o2.burst.allsky.stamp\n');
    fid.write('notification = never\n');
    fid.write('getenv = true\n');
    fid.write('log =' + opts.logDir + 'simsurvey.log\n')
    fid.write('+MaxHours = 24\n');
    fid.write('universe = vanilla\n');
    fid.write('queue 1\n');
    fid.close()
    
if opts.doQuerySkymap:
    fid = open(os.path.join(condorDir,'condor_query.sub'),'w')
    fid.write('executable = %s/query_skymap.py\n'%dir_path)
    fid.write('output = logs/out.$(jobNumber)\n');
    fid.write('error = logs/err.$(jobNumber)\n');
    fid.write('arguments = --skymap $(skymap) --level 90 --after-trigger True --jd-trigger $(jd) --min-days $(min_days) --within-days 3 --theta $(theta) --phi $(phi) --ndethist $(ndethist) --out $(out) --secrets %s --min-mag $(min_mag)\n' % opts.secrets)
    fid.write('requirements = OpSys == "LINUX"\n');
    fid.write('request_memory = 16384\n');
    fid.write('request_cpus = 1\n');
    fid.write('accounting_group = ligo.dev.o2.burst.allsky.stamp\n');
    fid.write('notification = never\n');
    fid.write('getenv = true\n');
    fid.write('log = /local/michael.coughlin/query_skymap.log\n')
    fid.write('+MaxHours = 24\n');
    fid.write('universe = vanilla\n');
    fid.write('queue 1\n');
    fid.close()
    
