from penquins import Kowalski
import io
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np 
import astropy
import astropy.stats


k = Kowalski(username='sanand', password='PrettyStr0ngPa$$w0rd')

# k = Kowalski(host='localhost', port=8000, protocol='http', username='admin', password='admin')

connection_ok = k.check_connection()
print(f'Connection OK: {connection_ok}')
collection_ZTF_alerts = 'ZTF_alerts'
collection_ZTF_alerts_aux = 'ZTF_alerts_aux'

def plot_cutouts(name,s_min,s_max):
    q = {"query_type": "general_search", "query": "db['ZTF_alerts'].find({'objectId': {'$eq': '"+name+"'}})" }
    r = k.query(query=q,timeout=30)
    # getting cutouts from the arvo packages kowalski
    
    candidate = r['result_data']['query_result'][0]
    # type(candidate['cutoutScience']['stampData']) is bytes
    
    sdir = './temp/'
    f1=open(sdir+candidate['objectId']+'-sci.fits','wb').write(io.BytesIO(candidate['cutoutScience']['stampData']).getvalue()) 
    f2=open(sdir+candidate['objectId']+'-ref.fits','wb').write(io.BytesIO(candidate['cutoutTemplate']['stampData']).getvalue())
    f3=open(sdir+candidate['objectId']+'-dif.fits','wb').write(io.BytesIO(candidate['cutoutDifference']['stampData']).getvalue())
    
    
    cm = plt.cm.cubehelix
    stamp_ext = 'fits'
    
    fig = plt.figure(figsize=(7, 9))

    ax1 = plt.subplot(121)
    ax1.margins(0.05)           # Default margin is 0.05, value 0 means fit, Values >0.0 zoom out
    itype='-sci.'
    img_data = fits.getdata(sdir+candidate['objectId']+itype+stamp_ext)
    img_data[np.isnan(img_data)]=np.nanmedian(img_data)

    ax1.set_xticks([])
    ax1.set_yticks([])

    imgd_scaled = np.log10(img_data)
    vmax = np.sort(imgd_scaled[28:36,28:36].flatten())[-3]
    npixel = (len(img_data)+1)**2
    imgd_flat = img_data.flatten()
    imgd_scaled[imgd_scaled<0]=np.nanmedian(img_data)
    v_onesig = np.log10(np.nanmedian(img_data) - astropy.stats.mad_std(imgd_flat[np.isnan(imgd_flat)==False])*1.4826)
    vmin= max(v_onesig, np.nanmin(imgd_scaled))
    ax1.imshow(imgd_scaled, cmap=cm, vmax=s_max[0]*vmax, vmin=vmin*s_min[0])

    ax2 = plt.subplot(122)
    ax2.margins(2, 2)       
    itype='-ref.'
    img_data = fits.getdata(sdir+candidate['objectId']+itype+stamp_ext)
    img_data[np.isnan(img_data)]=np.nanmedian(img_data)
    imgd_scaled = np.log10(img_data)
    vmin= max(v_onesig, np.nanmin(imgd_scaled))
    ax2.imshow(imgd_scaled, cmap=cm, vmax=s_max[1]*vmax, vmin=vmin*s_min[1])

    ax2.set_xticks([])
    ax2.set_yticks([])
    
    # marks for the candidate
    ax2.axhline(y=31,xmin=0.5-0.15,xmax=0.5-0.05, c='white',linewidth=1.75,alpha=0.85)
    ax2.axhline(y=31,xmin=0.5+0.05,xmax=0.5+0.15, c='white',linewidth=1.75,alpha=0.85)
    ax2.axvline(x=31,ymin=0.5-0.15,ymax=0.5-0.05, c='white',linewidth=1.75,alpha=0.85)
    ax2.axvline(x=31,ymin=0.5+0.05,ymax=0.5+0.15, c='white',linewidth=1.75,alpha=0.85)

    ax1.axhline(y=31,xmin=0.5-0.15,xmax=0.5-0.05, c='white',linewidth=1.75,alpha=0.85)
    ax1.axhline(y=31,xmin=0.5+0.05,xmax=0.5+0.15, c='white',linewidth=1.75,alpha=0.85)
    ax1.axvline(x=31,ymin=0.5-0.15,ymax=0.5-0.05, c='white',linewidth=1.75,alpha=0.85)
    ax1.axvline(x=31,ymin=0.5+0.05,ymax=0.5+0.15, c='white',linewidth=1.75,alpha=0.85)


    ax1.set_title(candidate['objectId'],fontsize=14,loc= 'left',fontweight='bold')
    ax2.set_title('REF',fontsize=14,loc= 'left',fontweight='bold')
    plt.subplots_adjust(wspace=.1, hspace=0)
    plt.savefig('./cuts/'+candidate['objectId']+'_cut.pdf',bbox_inches='tight')

def get_cutouts_path(name):
    q = {"query_type": "general_search", "query": "db['ZTF_alerts'].find({'objectId': {'$eq': '"+name+"'}})" }
    r = k.query(query=q,timeout=30)
    # getting cutouts from the arvo packages kowalski
    
    candidate = r['result_data']['query_result'][0]
    # type(candidate['cutoutScience']['stampData']) is bytes
    
    sdir = './temp/'
    f1=open(sdir+candidate['objectId']+'-sci.fits','wb').write(io.BytesIO(candidate['cutoutScience']['stampData']).getvalue()) 
    f2=open(sdir+candidate['objectId']+'-ref.fits','wb').write(io.BytesIO(candidate['cutoutTemplate']['stampData']).getvalue())
    f3=open(sdir+candidate['objectId']+'-dif.fits','wb').write(io.BytesIO(candidate['cutoutDifference']['stampData']).getvalue())
    
    
    cm = plt.cm.cubehelix
    stamp_ext = 'fits'
    
    return sdir+candidate['objectId']+'-sci.'+stamp_ext,sdir+candidate['objectId']+'-ref.'+stamp_ext
    
#     fig = plt.figure(figsize=(7, 9))

#     ax1 = plt.subplot(121)
#     ax1.margins(0.05)           # Default margin is 0.05, value 0 means fit, Values >0.0 zoom out
#     itype='-sci.'
#     img_data = fits.getdata(sdir+candidate['objectId']+itype+stamp_ext)
#     img_data[np.isnan(img_data)]=np.nanmedian(img_data)

#     ax1.set_xticks([])
#     ax1.set_yticks([])

#     imgd_scaled = np.log10(img_data)
#     vmax = np.sort(imgd_scaled[28:36,28:36].flatten())[-3]
#     npixel = (len(img_data)+1)**2
#     imgd_flat = img_data.flatten()
#     imgd_scaled[imgd_scaled<0]=np.nanmedian(img_data)
#     v_onesig = np.log10(np.nanmedian(img_data) - astropy.stats.mad_std(imgd_flat[np.isnan(imgd_flat)==False])*1.4826)
#     vmin= max(v_onesig, np.nanmin(imgd_scaled))
#     ax1.imshow(imgd_scaled, cmap=cm, vmax=s_max[0]*vmax, vmin=vmin*s_min[0])

#     ax2 = plt.subplot(122)
#     ax2.margins(2, 2)       
#     itype='-ref.'
#     img_data = fits.getdata(sdir+candidate['objectId']+itype+stamp_ext)
#     img_data[np.isnan(img_data)]=np.nanmedian(img_data)
#     imgd_scaled = np.log10(img_data)
#     vmin= max(v_onesig, np.nanmin(imgd_scaled))
#     ax2.imshow(imgd_scaled, cmap=cm, vmax=s_max[1]*vmax, vmin=vmin*s_min[1])

#     ax2.set_xticks([])
#     ax2.set_yticks([])
    
#     # marks for the candidate
#     ax2.axhline(y=31,xmin=0.5-0.15,xmax=0.5-0.05, c='white',linewidth=1.75,alpha=0.85)
#     ax2.axhline(y=31,xmin=0.5+0.05,xmax=0.5+0.15, c='white',linewidth=1.75,alpha=0.85)
#     ax2.axvline(x=31,ymin=0.5-0.15,ymax=0.5-0.05, c='white',linewidth=1.75,alpha=0.85)
#     ax2.axvline(x=31,ymin=0.5+0.05,ymax=0.5+0.15, c='white',linewidth=1.75,alpha=0.85)

#     ax1.axhline(y=31,xmin=0.5-0.15,xmax=0.5-0.05, c='white',linewidth=1.75,alpha=0.85)
#     ax1.axhline(y=31,xmin=0.5+0.05,xmax=0.5+0.15, c='white',linewidth=1.75,alpha=0.85)
#     ax1.axvline(x=31,ymin=0.5-0.15,ymax=0.5-0.05, c='white',linewidth=1.75,alpha=0.85)
#     ax1.axvline(x=31,ymin=0.5+0.05,ymax=0.5+0.15, c='white',linewidth=1.75,alpha=0.85)


#     ax1.set_title(candidate['objectId'],fontsize=14,loc= 'left',fontweight='bold')
#     ax2.set_title('REF',fontsize=14,loc= 'left',fontweight='bold')
#     plt.subplots_adjust(wspace=.1, hspace=0)
#     plt.savefig('./cuts/'+candidate['objectId']+'_cut.pdf',bbox_inches='tight')

# plot_cutouts('ZTF20abwysqy',[1.001,1.001],[1.001,1.001],k)
