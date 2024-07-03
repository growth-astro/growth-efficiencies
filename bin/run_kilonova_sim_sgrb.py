import os
import numpy as np

import datetime

from astropy import units as u
from astropy.time import Time, TimeDelta
from scipy import integrate

t0 = datetime.datetime.now()

path_repo = '/Users/anasaguescarracedo/Dropbox/PhD/github_repositories/'
path_grb = '/Users/anasaguescarracedo/Dropbox/PhD/SGRB/'

use_rate = True
use_doPlot = False

os.chdir(path_repo+'growth-too-papers/bin/')

events = [
{'n':'GRB180523B',
'o': path_grb+'GRB180523B/',
'start_time': Time(2458262.2823, format='jd', scale='utc').isot,
'skymap': path_grb+'GRB180523B/glg_healpix_all_bn180523782.fit',
'f': path_grb+'ztf_obsfile_status1_v2.csv'
},
{'n':'GRB180626C',
'o': path_grb+'GRB180626C/',
'start_time': Time(2458295.8916, format='jd', scale='utc').isot,
'skymap': path_grb+'GRB180626C/glg_healpix_all_bn180626392.fit',
'f': path_grb+'ztf_obsfile_status1_v2.csv'
},
{'n':'GRB180715B',
'o': path_grb+'GRB180715B/',
'start_time': Time(2458315.2412, format='jd', scale='utc').isot,
'skymap': path_grb+'GRB180715B/glg_healpix_all_bn180715741.fit',
'f': path_grb+'ztf_obsfile_status1_v2.csv'
},
{'n':'GRB180728B',
'o': path_grb+'GRB180728B/',
'start_time': Time(2458328.3819, format='jd', scale='utc').isot,
'skymap': path_grb+'GRB180728B/glg_healpix_all_bn180728882.fit',
'f': path_grb+'ztf_obsfile_status1_v2.csv'
},
{'n':'GRB180913A',
'o': path_grb+'GRB180913A/',
'start_time': Time(2458375.2834, format='jd', scale='utc').isot,
'skymap': path_grb+'GRB180913A/glg_healpix_all_bn180913783.fit',
'f': path_grb+'ztf_obsfile_status1_v2.csv'
},
{'n':'GRB181126B',
'o': path_grb+'GRB181126B/',
'start_time': Time(2458448.6617, format='jd', scale='utc').isot,
'skymap': path_grb+'GRB181126B/glg_healpix_all_bn181126162.fit',
'f': path_grb+'ztf_obsfile_status1_v2.csv'
}]


LF = np.genfromtxt('LF_sgrb.txt',delimiter=',',unpack=True)
n,edges = np.histogram(LF[0],bins=20)

masks = [np.array((LF[0]>edges[i]) * (LF[0]<edges[i+1])) for i in range(len(edges)-1) ]
prob = np.array([integrate.simps(LF[1][mask] ,LF[0][mask]) for mask in masks])

N_transient = 1e5
ntran =  np.round(N_transient*prob,0).astype(int)

E0s = 10**((edges[:-1]+edges[1:])/2)
# E0s = ['1e46', '1e47', '1e48', '1e49', '1e50', '1e51', '1e52' , '1e53']

for event in events:
    print(event['n'])
    outputdir = event['o']
    start_time = event['start_time']
    skymap = event['skymap']
    observinglog = event['f']
    dmax = 500
    NPoints = 1

    for i,E0 in enumerate(E0s):
        outputfile = '{:.2e}'.format(E0)+'.pkl'
        ntransient = ntran[i]
        rate = float(1e-3)

        system_command = 'python3 kilonova_sim.py --doSimSurvey --NPoints='+str(NPoints)+' -a '+path_grb+' --dmax='+str(dmax)+' -n '+str(ntransient)+' --thetadist=uniform --thetai=0 --thetaf=20 -o '+outputdir+' --start_time '+start_time+' --dt=2.0 --skymap '+skymap+' -f '+observinglog+' -m afterglow --E0='+str(E0)
        if use_doPlot:
            system_command += ' --doPlots '

        if use_rate:
            system_command += ' --useRate -r '+str(rate)
            outputfile = '{:.2e}'.format(E0)+'_rate'+'{:.0e}'.format(rate)+'.pkl'

        system_command += ' --pickleFile='+outputfile

        os.system(system_command)
