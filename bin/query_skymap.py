#!/usr/bin/env python

'''
Query Kowalski searching for transients
given a set of constraints.

'''
import os
import numpy as np
import json
import pdb

from astropy.time import Time
from astropy.io import ascii
#from astropy.io import fits
from ligo.skymap.io import fits
from astropy.table import Table
import requests

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
import healpy as hp

def rotate_map(hmap, rot_theta, rot_phi):
    """
    Take hmap (a healpix map array) and return another healpix map array 
    which is ordered such that it has been rotated in (theta, phi) by the 
    amounts given.
    """
    nside = hp.npix2nside(len(hmap))

    # Get theta, phi for non-rotated map
    t,p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside))) #theta, phi

    # Define a rotator
    r = hp.Rotator(deg=False, rot=[rot_phi,rot_theta])

    # Get theta, phi under rotated co-ordinates
    trot, prot = r(t,p)

    # Interpolate map onto these co-ordinates
    rot_map = hp.get_interp_val(hmap, trot, prot)

    return rot_map


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'Yes', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'No', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def check_args(args):
    ''' Check that the jd of the trigger is provided
    if only detections after a certain date are required '''
    if args.after_trigger == True:
        if args.jd_trigger == -1:
            print('>>> Expected --jd-trigger')
            print('>>> Exiting...')
            exit()

    '''Check required duration'''
    if args.min_days >= args.max_days:
        print('>>> --min-days must be smaller than --max-days')
        print('>>> Exiting...')
        exit()
    '''Check that RA and Dec arrays have the same length  '''

    return


def print_query_params(args, ra_center, dec_center):
    '''Print a summary of the query parameters'''

    print("#-----")
    print("Cone search parameters:")
    print(f"A list of {len(ra_center)} coordinate pairs will be explored")
    print(f"Search radius {args.radius} arcmin")
    if args.after_trigger or args.jd_trigger > 0:
        print(f"Only sources detected for the first time after {Time(args.jd_trigger, format='jd').iso} will be considered")
    print(f"Minimum time between the first and last alert {args.min_days} days")
    print(f"Maximum time between the first and last alert {args.max_days} days")    
    print(f"Query divided in {args.slices} slices")
    print("#-----")
    print(" ")

    return


def get_programidx(program_name, username, password):
    ''' Given a marshal science program name, it returns its programidx'''

    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_programs.cgi', auth=(username, password))
    programs=json.loads(r.text)
    program_dict={p['name']:p['programidx'] for i,p in enumerate(programs)}

    try:
        return program_dict[program_name]
    except KeyError:
        print(f'The user {username} does not have access to the program {program_name}')
        return None


def get_candidates_growth_marshal(program_name, username, password):
    ''' Query the GROWTH db for the science programs '''

    programidx=get_programidx(program_name, username, password)
    if programidx==None:
        return None
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi', \
        auth=(username, password), data={'programidx':str(programidx)})
    sources=json.loads(r.text)
    sources_out=[]
    for s in sources:
            coords=SkyCoord(ra=s['ra']*u.deg, dec=s['dec']*u.deg, frame='icrs')
            sources_out.append({"name":s['name'], "ra":coords.ra, "dec":coords.dec, \
	        "classification":s['classification'], "redshift":s['redshift'], "creation_date":s['creationdate']})

    return sources_out


def query_kowalski_clu(username, password, clu):
    '''Query kowalski to get a table of CLU galaxies. '''

    k = Kowalski(username=username, password=password, verbose=False)
    q = {"query_type": "general_search", 
        "query": "db['CLU_20180513'].find({},{'distmpc': 1})" 
        }
    r = k.query(query=q)

    return r


def check_clu_transients(sources_kowalski, clu_sources):
    '''Check if the selected sources are present in the 
    CLU science program.  If so, print out the relevant information.'''

    sources_in_clu = []
    sources_not_in_clu = []
    list_clu_sources = list(s['name'] for s in clu_sources)

    for source in sources_kowalski:
        print("-------")
        if source in list_clu_sources:
            clu_source = clu_sources[np.where(np.array(list_clu_sources) == source)[0][0]]
            try:
                for k in clu_source.keys():
                    print(f"{k}: {clu_source[k]}")
                sources_in_clu.append(source)
            except:
                pdb.set_trace()
        else:
            print(f"{source} was not saved in CLU")
            sources_not_in_clu.append(source)
        print("-------")
    print("Summary:")
    print(f"Sources saved in CLU: {sources_in_clu}")
    print(f"Sources not saved in CLU: {sources_not_in_clu}")

    return


def read_skymap(skymap_filename, theta=0, phi=0):
    '''Read the healpix skymap'''

    skymap, metadata = fits.read_sky_map(skymap_filename, nest=True)
    nside = hp.get_nside(skymap)
    npix = hp.nside2npix(nside)
    skymap = skymap[hp.ring2nest(nside, np.arange(npix))]
    if (not theta==0) or (not phi == 0):
        skymap = rotate_map(skymap, np.deg2rad(theta), np.deg2rad(phi))

    return skymap


def tesselation_spiral(FOV, scale=0.80):
    FOV = np.pi*FOV*FOV*scale

    area_of_sphere = 4*np.pi*(180/np.pi)**2
    n = int(np.ceil(area_of_sphere/FOV))
    print("Using %d points to tile the sphere..."%n)

    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(n)
    z = np.linspace(1 - 1.0 / n, 1.0 / n - 1, n)
    radius = np.sqrt(1 - z * z)

    points = np.zeros((n, 3))
    points[:,0] = radius * np.cos(theta)
    points[:,1] = radius * np.sin(theta)
    points[:,2] = z

    ra, dec = hp.pixelfunc.vec2ang(points, lonlat=True)

    return ra, dec


def do_getfields(healpix, FOV=60/3600.0, ra=None, dec=None, radius=None, level=None):
    from ligo.skymap import postprocess
    import matplotlib

    ras, decs = tesselation_spiral(FOV, scale=0.80)
    if (not ra is None) and (not dec is None):
        dist = angular_distance(ras, decs, ra, dec) 
        idx = np.where(dist <= radius)[0]
        ras, decs = ras[idx], decs[idx]
    elif (not level is None):
        cls = 100 * postprocess.find_greedy_credible_levels(healpix)
        paths = postprocess.contour(cls, [level], degrees=True, simplify=True)
        paths = paths[0]

        pts = np.vstack((ras, decs)).T
        idx = np.zeros((len(ras)))
        for path in paths:
            polygon = matplotlib.path.Path(path)
            check = polygon.contains_points(pts)
            check = list(map(int, check))
            idx = np.maximum(idx, check)
        idx = np.where(idx == 1)[0] 
        ras, decs = ras[idx], decs[idx]

    return ras, decs


def query_kowalski(username, password, ra_center, dec_center, radius,
                   jd_trigger, min_days, max_days, slices, ndethist_min,
                   within_days, after_trigger=True, min_mag=np.inf):
    '''Query kowalski and apply the selection criteria'''

    from penquins import Kowalski

    k = Kowalski(username=username, password=password, verbose=False)
    #Initialize a set for the results
    set_objectId_all = set([])
    slices = slices + 1
    for slice_lim,i in zip(np.linspace(0,len(ra_center),slices)[:-1], np.arange(len(np.linspace(0,len(ra_center),slices)[:-1]))):
        #if slice_lim < 2013:
        #    continue
        try:
            ra_center_slice = ra_center[int(slice_lim):int(np.linspace(0,len(ra_center),slices)[:-1][i+1])]
            dec_center_slice = dec_center[int(slice_lim):int(np.linspace(0,len(dec_center),slices)[:-1][i+1])]
        except IndexError:
            ra_center_slice = ra_center[int(slice_lim):]
            dec_center_slice = dec_center[int(slice_lim):]
        coords_arr = []
        for ra, dec in zip(ra_center_slice, dec_center_slice):
            try:
                #Remove points too far south for ZTF.  Say, keep only Dec>-40 deg to be conservative
                if dec < -40.:
                    continue
                coords=SkyCoord(ra=float(ra)*u.deg, dec=float(dec)*u.deg)
                coords_arr.append((coords.ra.deg,coords.dec.deg))
            except ValueError:
                print("Problems with the galaxy coordinates?")
                continue
        #Correct the minimum number of detections
        ndethist_min_corrected = int(ndethist_min - 1)
        #Correct the jd_trigger if the user specifies to query also before the trigger
        if after_trigger == False:
            jd_trigger = 0
        try: 
            print(f"slice: {int(slice_lim)}:{int(np.linspace(0,len(ra_center),slices)[:-1][i+1])}" )
        except:
            print(f"slice: {int(slice_lim)}:{int(len(ra_center))}" )
        q = {"query_type": "cone_search",
        "object_coordinates": {
             "radec": f"{coords_arr}", 
             "cone_search_radius": f"{radius}",
             "cone_search_unit": "arcmin"
         },
         "catalogs": {
             "ZTF_alerts": {
                 "filter": {
		     "candidate.jd": {'$gt': jd_trigger},
		     "candidate.drb": {'$gt': 0.5},
		     "candidate.ndethist": {'$gt': ndethist_min_corrected},
		     "candidate.jdstarthist": {'$gt': jd_trigger}
		     },
                 "projection": {
                     "objectId": 1,
                     "candidate.rcid": 1,
                     "candidate.ra": 1,
                     "candidate.dec": 1,
                     "candidate.jd": 1,
                     "candidate.ndethist": 1,
                     "candidate.jdstarthist": 1,
                     "candidate.jdendhist": 1,
                     "candidate.jdendhist": 1,
                     "candidate.magpsf": 1,
                     "candidate.sigmapsf": 1,
                     "candidate.fid": 1,
                     "candidate.programid": 1,
                     "candidate.isdiffpos": 1,
                     "candidate.ndethist": 1,
                     "candidate.ssdistnr": 1,
                     "candidate.rb": 1,
                     "candidate.drb": 1,
                     "candidate.distpsnr1": 1,   
                     "candidate.sgscore1": 1,
                     "candidate.srmag1": 1,
                     "candidate.distpsnr2": 1,   
                     "candidate.sgscore2": 1,
                     "candidate.srmag2": 1,
                     "candidate.distpsnr3": 1,   
                     "candidate.sgscore3": 1,
                     "candidate.srmag3": 1
                 }
             }
         },
        "kwargs": {"hint": "gw01"}
         }

        #Perform the query
        r = k.query(query=q)
        print('Search completed for this slice.')

#        #Dump the results in a json file
#        with open(f'results_clu25Mpc_1week_{i+1}.json', 'w') as j:
#            json.dump(r, j)

        #Identify 'candid' for all relevant candidates
        objectId_list = []
        with_neg_sub = []
        old = []
        out_of_time_window = []
        stellar_list = []
        try:
            keys_list=list(r['result_data']['ZTF_alerts'].keys())
        except:
            print("Error in the keys list?? Check 'r' ")
            #Try one more time 
            print("Trying to query again the same slice")
            try:
                r = k.query(query=q)
                keys_list=list(r['result_data']['ZTF_alerts'].keys())
            except:
                print("The query failed again, skipping slice..")
                continue

        for key in keys_list:
            all_info=r['result_data']['ZTF_alerts'][key]
            
            for info in all_info:
                if info['objectId'] == 'ZTF19abyfbii':
                    pdb.set_trace()
                if info['objectId'] in old:
                    continue
                if info['objectId'] in stellar_list:
                    continue
                try:
                    if info['candidate']['drb'] < 0.5:
                        continue
                except:
                    do = 'do nothing.'
                if np.abs(info['candidate']['ssdistnr']) < 10:
                    continue
                if info['candidate']['isdiffpos'] in ['f',0]:
                    with_neg_sub.append(info['objectId'])
                if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) < min_days:
                    continue
                if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) > max_days:
                    old.append(info['objectId'])
                if (info['candidate']['jdstarthist'] - jd_trigger) > within_days:
                    old.append(info['objectId'])
                if after_trigger == True:
                    if (info['candidate']['jdendhist'] - jd_trigger) > max_days:
                        out_of_time_window.append(info['objectId'])
                else:
                    if (info['candidate']['jdendhist'] - info['candidate']['jdstarthist']) > max_days:
                        out_of_time_window.append(info['objectId'])
                try:
                    if (np.abs(info['candidate']['distpsnr1']) < 3. and info['candidate']['sgscore1'] > 0.0):
                        stellar_list.append(info['objectId'])
                except:
                    do = 'do nothing.'
                try:
                    if (np.abs(info['candidate']['distpsnr1']) < 15. and info['candidate']['srmag1'] < 15. and info['candidate']['sgscore1'] >= 0.5):
                        continue
                except:
                    do = 'do nothing.'
                try:
                    if (np.abs(info['candidate']['distpsnr2']) < 15. and info['candidate']['srmag2'] < 15. and info['candidate']['sgscore2'] >= 0.5):
                        continue
                except:
                    do = 'do nothing.'
                try:
                    if (np.abs(info['candidate']['distpsnr3']) < 15. and info['candidate']['srmag3'] < 15. and info['candidate']['sgscore3'] >= 0.5):
                        continue
                except:
                    do = 'do nothing.'
                if np.abs(info['candidate']['magpsf']) > min_mag:
                    print('removing this object...')
                    continue

                objectId_list.append(info['objectId'])

        set_objectId = set(objectId_list)

        #Remove those objects with negative subtraction
        for n in set(with_neg_sub):
            try:
                set_objectId.remove(n)
            except:
                do = 'do nothing'

        #Remove stellar objects
        for n in set(stellar_list):
            try:
                set_objectId.remove(n)
            except:
                do = 'do nothing'

        #Remove those objects considered old
        for n in set(old):
            try:
                set_objectId.remove(n)
            except:
                do = 'do nothing'

        #Remove those objects whole alerts go bejond jd_trigger+max_days
        for n in set(out_of_time_window):
            try:
                set_objectId.remove(n)
            except:
                do = 'do nothing'
        print(set_objectId)

        set_objectId_all = set_objectId_all | set_objectId
        print("Cumulative:", set_objectId_all)

        '''
        print('----stats-----')
        print('Number of sources with negative sub: ', len(set(with_neg_sub)))
        print('Number of sources with only pos subtraction: ', len(set_objectId))
        print(f"Number of sources older than {max_days} days: {len(set(old))}, specifically {set(old)}")
        '''

    return set_objectId_all


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Query kowalski.')
    parser.add_argument('--skymap', dest='skymap_filename', type=str, required=False, \
    help='Skymap filename', default=None)
    parser.add_argument('--secrets', dest='secrets', type=str, required=False, \
    help='secrets filename', default="secrets.csv")
    parser.add_argument('--level', dest='level', type=float, required=False, \
    help='Enclosed probability', default = 90)
    parser.add_argument('--fov', dest='fov', type=float, required=False, \
    help='Field of view of each cone (radius, in arcmin)', default = 60)
    parser.add_argument('--ra-center', dest='ra_center', nargs='+', required=False, \
    help='Right ascension of the center (array, in degrees)')
    parser.add_argument('--dec-center', dest='dec_center', nargs='+', required=False, \
    help='Declination of the center (array, in degrees)')
    parser.add_argument('--radius', dest='radius', type=float, required=False, \
    help='Search radius (min), by default radius = fov', default = None)
    parser.add_argument('--after-trigger', dest='after_trigger', type=str2bool, required=False, \
    help='Query only alerts whose first detection occurred after a certain date. \
    If this boolean value is True, then --jd-trigger must be given.', default=True)  
    parser.add_argument('--jd-trigger', dest='jd_trigger', type=float, required=False, \
    help='Julian Day of the trigger', default = -1.)            
    parser.add_argument('--min-days', dest='min_days', type=float, required=False, \
    help='Minimum time (days) between the first and last alert', default = 0.)
    parser.add_argument('--max-days', dest='max_days', type=float, required=False, \
    help='Maximum time (days) between the first and last alert', default = 10000.)
    parser.add_argument('--within-days', dest='within_days', type=float, required=False, \
    help='Maximum time (days) between the jd-trigger and the first alert', default = 1000.)
    parser.add_argument('--min-mag', dest='min_mag', type=float, required=False, \
    help='Minimum mag at peak', default = np.inf)
    parser.add_argument('--ndethist', dest='ndethist_min', type=int, required=False, \
    help='Minimum number of detections', default=2)
    parser.add_argument('--slices', dest='slices', type=int, required=False, \
    help='Number (integer) of slices in which the query will be devided', default = 10)
    parser.add_argument('--out', dest='out', type=str, required=False, \
    help='Output filename', default = 'results.txt')

    parser.add_argument('--theta', dest='theta', type=float, required=False, \
    help='Theta rotation of skymap', default = 0.0)
    parser.add_argument('--phi', dest='phi', type=float, required=False, \
    help='Phi rotation of skymap', default = 0.0)
    
    args = parser.parse_args()

    #Check that the input args make sense
    check_args(args)

    #Read the skymap and create the tessellation
    if args.skymap_filename != None:
        healpix = read_skymap(args.skymap_filename, theta=args.theta, phi=args.phi)
        ra_center, dec_center = do_getfields(healpix, FOV = args.fov/60., level = args.level)
    else:
        ra_center, dec_center = args.ra_center, args.dec_center 

    if args.radius is None:
        args.radius = args.fov
    #Pre-select coordinates above Dec = -30
    ra_center = ra_center[np.where(dec_center > -30)]
    dec_center = dec_center[np.where(dec_center > -30)]

    #Print a summary of the query input
    print_query_params(args, ra_center, dec_center)

    #Read the secrets
    secrets = ascii.read(args.secrets, format = 'csv')

    username = secrets['kowalski_user'][0]
    password = secrets['kowalski_pwd'][0]

    #Query kowalski
    sources_kowalski = query_kowalski(username, password, ra_center, dec_center, args.radius, args.jd_trigger, args.min_days, args.max_days, args.slices, args.ndethist_min, args.within_days, after_trigger=args.after_trigger, min_mag=args.min_mag)

    outdir = "/".join(args.out.split("/")[:-1])
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    #Print results to an output text file
    with open(args.out, 'w') as f:
        f.write(f"{sources_kowalski} \n")

    #Check the CLU science program on the Marshal
    username_marshal = secrets['marshal_user'][0]
    password_marshal= secrets['marshal_pwd'][0]
    
    program_name='Census of the Local Universe'
    #clu_sources = get_candidates_growth_marshal(program_name, username_marshal, password_marshal)    

    #For each transient check if it is present in the CLU science program
    #check_clu_transients(sources_kowalski, clu_sources)

    print("Done.")
