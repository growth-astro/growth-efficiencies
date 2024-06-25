import numpy as np
import glob
from read_avro import read_avro

# to download alerts:
#  https://ztf.uw.edu/alerts/partnership/
# ask for user/pass -- http://www.oir.caltech.edu/twiki_ptf/pub/ZTF/ZTFCollaborationMeetingAugust2018/ 


def compiledFunction(current_observation,t_gw = 2458598.846134): # can be changed, just copy the filter from the GROWTH marshal page
    
    filteron = False
    annotations={}
    calccount=10000
    bright = False
    nopointunderneath = True
    mover = True
    real = False
    tdiff = (-99.0)
    magdiff = (-99)
    riserate = (-99)
    decayrate = (-99)
    hostgr = (-99)
    hostri = (-99)
    positivesubtraction = False
    brightstar = False
    highlum = False
    tooflag = 0
    jdendhist = 0
    jdstarthist = 0
    new = True
    idname = current_observation['objectId']
    prevcandidates = current_observation['prv_candidates']
    m_now = current_observation['candidate']['magpsf']
    m_app = current_observation['candidate']['magap']
    t_now = current_observation['candidate']['jd']
    fid_now = current_observation['candidate']['fid']
    sgscore = current_observation['candidate']['sgscore1']
    sgscore2 = current_observation['candidate']['sgscore2']
    sgscore3 = current_observation['candidate']['sgscore3']
    srmag = current_observation['candidate']['srmag1']
    srmag2 = current_observation['candidate']['srmag2']
    srmag3 = current_observation['candidate']['srmag3']
    sgmag = current_observation['candidate']['sgmag1']
    simag = current_observation['candidate']['simag1']
    rbscore = current_observation['candidate']['rb']
    magnr = current_observation['candidate']['magnr']
    distnr = current_observation['candidate']['distnr']
    distpsnr1 = current_observation['candidate']['distpsnr1']
    distpsnr2 = current_observation['candidate']['distpsnr2']
    distpsnr3 = current_observation['candidate']['distpsnr3']
    scorr = current_observation['candidate']['scorr']
    fwhm = current_observation['candidate']['fwhm']
    elong = current_observation['candidate']['elong']
    nbad = current_observation['candidate']['nbad']
    chipsf = current_observation['candidate']['chipsf']
    tooflag = current_observation['candidate']['tooflag']
    jdstarthist = current_observation['candidate']['jdstarthist']
    jdendhist = current_observation['candidate']['jdendhist']
    psfminap = m_now - m_app
    bright = m_now < 99.0
    if (current_observation['candidate']['isdiffpos'] and (current_observation['candidate']['isdiffpos'] == 't' or current_observation['candidate']['isdiffpos'] == '1')):
        positivesubtraction = True
        calccount -= 2
    if (rbscore and rbscore > 0.25):
        real = True
        calccount -= 2
    if (sgscore and distpsnr1 and sgscore > 0.6 and distpsnr1 < 2):
        nopointunderneath = False
        calccount -= 2
    if ((distpsnr1 and srmag and distpsnr1 < 20 and srmag < 15.0 and srmag > 0 and sgscore > 0.49) or (distpsnr2 and srmag2 and distpsnr2 < 20 and srmag2 < 15.0 and srmag2 > 0 and sgscore2 > 0.49) or (distpsnr3 and srmag3 and distpsnr3 < 20 and srmag3 < 15.0 and srmag3 > 0 and sgscore3 > 0.49)):
        brightstar = True
        calccount -= 2
    if prevcandidates != None:
        for candidate in prevcandidates:
            calccount -= 2
            if (candidate['jd'] and candidate['magpsf'] and candidate['fid'] and candidate['isdiffpos'] and candidate['isdiffpos'] == 't'):
                dt = t_now - candidate['jd']
                if (dt > 0.01 and candidate['magpsf'] < 99):
                    mover = False
                    calccount -= 2
                calccount -= 3
            if calccount < 0:
                break

    histjd = jdendhist - jdstarthist
    if (histjd > 3):
        new = False
        calccount -= 2

    filteron = bright and nopointunderneath and ((not mover)) and real and positivesubtraction and tooflag and ((not brightstar)) and new# and (jdstarthist > 2458600.14022)
    return [filteron,bright , nopointunderneath , ((not mover)) , real , positivesubtraction , tooflag , ((not brightstar)) , new, jdstarthist > t_gw]

def apply_filter(alert_path,t_gw = 2458598.846134, test=False): # change t_gw according to evet (jd time)
    
    alert_folder = glob.glob(alert_path+'/*') # getting alerts
    print ('# of alerts', len(alert_folder))
    
    n=1.
    if test:
        n=0.05
    
    # analyzing alerts
    pass_filter=[]
    for i in range(int(len(alert_folder)*n)):
        if i%10000==0:
            if i==0:
                print('Starting ... \nAlert #\n')
            else:
                print(i)
        current_observation = read_avro(alert_folder[i])
        pass_filter.append(compiledFunction(current_observation,t_gw=t_gw))
    pass_filter = np.asarray(pass_filter)
    
    # printing info
    # [0] filteron, [1] bright , [2] nopointunderneath , [3] not mover , [4] real , [5] positivesubtraction , [6] tooflag , [7] not brightstar , [8] new , [9] jd > t_gw
    print ('Finish: ',str(i),' alerts analized')
    print ('too',np.sum(pass_filter.T[6]))
    print ('positive',np.sum(pass_filter.T[6]*pass_filter.T[5]))
    print ('real',np.sum(pass_filter.T[6]*pass_filter.T[5]*pass_filter.T[4]))
    print ('not stellar',np.sum(pass_filter.T[6]*pass_filter.T[5]*pass_filter.T[4]*pass_filter.T[2]))
    print ('far from star',np.sum(pass_filter.T[6]*pass_filter.T[5]*pass_filter.T[4]*pass_filter.T[2]*pass_filter.T[7]))
    print ('not mover',np.sum(pass_filter.T[6]*pass_filter.T[5]*pass_filter.T[4]*pass_filter.T[2]*pass_filter.T[7]*pass_filter.T[3]))
    print ('new',np.sum(pass_filter.T[6]*pass_filter.T[5]*pass_filter.T[4]*pass_filter.T[2]*pass_filter.T[7]*pass_filter.T[3]*pass_filter.T[8]))
    print ('jd',np.sum(pass_filter.T[6]*pass_filter.T[5]*pass_filter.T[4]*pass_filter.T[2]*pass_filter.T[7]*pass_filter.T[3]*pass_filter.T[8]*pass_filter.T[9]))

    return alert_folder,pass_filter


import healpy as hp
from astropy.table import Table, Column

def in_out(skymap,ra_obj,dec_obj):
    # Read skymap, calculate top pixels
    top_fraction = 0.90 # limit skymap top 90% region
    npix = len(skymap)
    nside = hp.npix2nside(npix)

    # Convert to astropy Table, easier for manipulation
    indices = np.arange(len(skymap))
    tm = Table(data=(indices, skymap), names=('id', 'prob'))
    tm.sort('prob')
    cs = np.cumsum(tm['prob'])
    cs.name='cumsum'
    tm.add_column(cs)

    top_pix = (tm['cumsum'] > 1 - top_fraction)
    tp = Column(data=top_pix, name="top")
    tm.add_column(tp)

    # Cast as a set for easier comparison below
    top_subset = set(tm['id'][tm['top']])

    inside = False
    pix = hp.ang2pix(nside,ra,dec, lonlat=True)
    if type(pix)==np.int64:
        if pix in top_subset:
            inside = True 
    else:
        inside =[]
        for i in range(len(pix)):
            if pix[i] in top_subset:
                inside.append(True)
            else: 
                inside.append(False)
    return inside

t_gw = 2458598.846134
alert_path = '/Users/tahumada/Downloads/ztf_partnership_20190425/' # path to the uncompressed folder
alerts,passed = apply_filter(alert_path,t_gw = t_gw, test=True)

#checkin inside outside -- change map

cond = passed.T[6]*passed.T[5]*passed.T[4]*passed.T[2]*passed.T[7]*passed.T[3]*passed.T[8]*passed.T[9]
ra,dec,inside = [],[],[]
map_path = '/Users/tahumada/Downloads/LALInference1.fits.gz,0'
skymap = hp.read_map(map_path)

for i,file in enumerate(np.asarray(alerts)[cond.astype('bool')]):
    current_observation = read_avro(file)
    ra.append(current_observation['candidate']['ra'])
    dec.append(current_observation['candidate']['dec'])

inside=(in_out(skymap,ra,dec))
print('inside ',sum(inside))

