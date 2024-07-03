import os
import numpy as np

import datetime

from astropy import units as u
from astropy.time import Time, TimeDelta
from scipy import integrate

def eiso(ep): # eq 2 https://arxiv.org/abs/1208.0429
    z = lambda e,k:10**(52.42+ (k*-0.15)) * (e/774.5)**(1.58+ k*0.28)
    return [z(ep,-1,),z(ep,1),z(ep,0)]



# workspaces
growth_bin = './'
output = '../data/SGRB/new3/'
maps = '../data/SGRB/maps/'
ztf_tiles = '../data/'

path_repo = '../'
path_grb = '../data/SGRB/'


os.chdir(growth_bin)

events = [
{'n':'GRB180523B',
'o': output+'GRB180523B/',
'start_time': Time(2458262.2823, format='jd', scale='utc').isot,
'skymap': maps+'glg_healpix_all_bn180523782.fit',
'f': ztf_tiles+'ztf_obsfile_status1_v2.csv',
 'ep':1430,
 'ep_err':687
},
{'n':'GRB180626C',
'o': output+'GRB180626C/',
'start_time': Time(2458295.8916, format='jd', scale='utc').isot,
'skymap':  maps+'glg_healpix_all_bn180626392.fit',
'f': ztf_tiles+'ztf_obsfile_status1_v2.csv',
 'ep':446,
 'ep_err':98
},
{'n':'GRB180715B',
'o': output+'GRB180715B/',
'start_time': Time(2458315.2412, format='jd', scale='utc').isot,
'skymap':  maps+'glg_healpix_all_bn180715741.fit',
'f': ztf_tiles+'ztf_obsfile_status1_v2.csv',
 'ep':559,
 'ep_err':112
},
{'n':'GRB180728B',
'o': output+'GRB180728B/',
'start_time': Time(2458328.3819, format='jd', scale='utc').isot,
'skymap':  maps+'glg_healpix_all_bn180728882.fit',
'f': ztf_tiles+'ztf_obsfile_status1_v2.csv',
 'ep':504,
 'ep_err':61
},
{'n':'GRB180913A',
'o': output+'GRB180913A/',
'start_time': Time(2458375.2834, format='jd', scale='utc').isot,
'skymap':  maps+'glg_healpix_all_bn180913783.fit',
'f': ztf_tiles+'ztf_obsfile_status1_v2.csv',
 'ep':444,
 'ep_err':175
},
{'n':'GRB181126B',
'o': output+'GRB181126B/',
'start_time': Time(2458448.6617, format='jd', scale='utc').isot,
'skymap':  maps+'glg_healpix_all_bn181126162.fit',
'f': ztf_tiles+'ztf_obsfile_status1_v2.csv',
 'ep':1049,
 'ep_err':241
}]

LF = np.genfromtxt('LF_sgrb.txt',delimiter=',',unpack=True)
n,edges = np.histogram(LF[0],bins=20)

masks = [np.array((LF[0]>edges[i]) * (LF[0]<edges[i+1])) for i in range(len(edges)-1) ]
prob = np.array([integrate.simps(LF[1][mask] ,LF[0][mask]) for mask in masks])
    
N_transient = 1e5
ntran =  np.round(N_transient*prob,0).astype(int)

E0s = 10**((edges[:-1]+edges[1:])/2)
E0s = 10**np.linspace(52,54.2,5)
E_ISO = True

distances =  [(0.1,500),
#               (500,1000),
#               (1000,2500),
#               (2500,5000),
#              (5000,7500),
#              (7500,10000),
#              (10000,12500),
#              (12500,15000),
             ]
sys_commands = []

import time
start = time.time()

for dmin,dmax in distances:
    for event in events:
        outputdir = event['o']
        start_time = event['start_time']
        skymap = event['skymap']
        observinglog = event['f']
        NPoints = 1
        if event['ep'] != 0 and E_ISO:
#             print(event['n'])
            eiso_max,eiso_min = np.log10(np.amax([eiso(event['ep']+event['ep_err']),eiso(event['ep']-event['ep_err'])])),\
                        np.log10(np.amin([eiso(event['ep']+event['ep_err']),eiso(event['ep']-event['ep_err'])]))
#             E0s = 10**np.linspace(eiso_min,eiso_max,10)
            E0s = 10**np.array([eiso_min,np.log10(eiso(event['ep'])[-1]),eiso_max])
#             print(np.log10(eiso(event['ep'])),eiso_max,eiso_min)
            print(event['n'],'$',np.round(np.log10(eiso(event['ep'])[-1]),2),\
                              '\ _{-',np.round(np.log10(eiso(event['ep'])[-1])-eiso_min,2),'}',\
                               '^{+',np.round(-np.log10(eiso(event['ep'])[-1])+eiso_max,2),'}$')
        elif E_ISO:
            E0s = 10**np.linspace(49,54,10)

        for (i,E0) in enumerate(E0s):
    #         E0 = E0s[7]    
            outputfile = str(dmax)+'_{:.2e}'.format(E0)+'.pkl'
        #         ntransient = ntran[i]
            ntransient = 7000

            system_command = 'python3 kilonova_sim.py --doSimSurvey --NPoints='+str(NPoints)+\
            ' -a '+path_grb+' --dmax='+str(dmax)+' --dmin='+str(dmin)+' -n '+str(ntransient)+\
            ' --thetadist=uniform --thetai=0 --thetaf=20 -o '+outputdir+' --start_time '+start_time+\
            ' --dt=2.0 --skymap '+skymap+' -f '+observinglog+' -m afterglow --E0='+str(E0)+' --pickleFile='+\
            outputfile
#             print(system_command)
            sys_commands.append(system_command)
            os.system(system_command)


end = time.time()
t1 = end-start
print( '-----'*20)
print( 'Time', t1)
print( '-----'*20)


