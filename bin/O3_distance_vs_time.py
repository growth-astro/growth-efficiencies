
from astropy import units as u, constants as c
from datetime import datetime
import numpy as np
from astropy.table import Table
from astropy.time import Time, TimeDelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.pyplot import rc,rcParams
from matplotlib import pylab

# some font set up.  Others may want different values
def computer_modern(ht_mult=1.0,wid_mult=1.0,font_mult=1.0):
    rc("font",**{"family":"serif","serif":["Computer Modern Roman"]})
    fig_width = 6.0*wid_mult
    fig_height = fig_width*0.8*ht_mult
    fig_size = [fig_width*wid_mult,fig_height*ht_mult]
    params = {"backend": "pdf",
              "font.size"      : 14*fig_width/8.5*font_mult,
              "axes.labelsize" : 14*fig_width/8.5*font_mult,
#              "text.fontsize"  : 12*fig_width/8.5*font_mult,
              "legend.fontsize": 12*fig_width/8.5*font_mult,
              "xtick.labelsize": 14*fig_width/8.5*font_mult,
              "ytick.labelsize": 14*fig_width/8.5*font_mult,
              "text.usetex"    : True,
              "figure.figsize" : fig_size,
              "axes.unicode_minus": True}
    rcParams.update(params)

fontsize=16
# pylab.rc('font',**{'family':'serif','serif':['Times New Roman']})
pylab.rc('xtick',labelsize=fontsize)
pylab.rc('ytick',labelsize=fontsize)
pylab.rc('text',usetex=True)

computer_modern(font_mult=1.4)

def drawPieMarker(xs, ys, ratios, sizes, colors, alpha, ax):
    assert sum(ratios) <= 1, 'sum of ratios needs to be < 1'

    markers = []
    previous = 0
    # calculate the points of the pie pieces
    for color, ratio in zip(colors, ratios):
        if ratio>1e-2:
            this = 2 * np.pi * ratio + previous
            x  = [0] + np.cos(np.linspace(previous, this, 50)).tolist() + [0]
            y  = [0] + np.sin(np.linspace(previous, this, 50)).tolist() + [0]
            xy = np.column_stack([x, y])
            previous = this
            markers.append({'marker':xy, 's':np.abs(xy).max()**2*np.array(sizes), 'facecolor':color, 'alpha': alpha, 'edgecolor': 'none'})

    # scatter each of the pie pieces to create pies
    for marker in markers:
        ax.scatter(xs, ys, **marker)

def get_size(area):
    # this is somewhat arbitrary
    # may need to be changed 
    return 40e3/area


def get_alpha(FAR,logFARmin=-15, logFARmax=-7, alphamin=0.25):
    # this is somewhat arbitrary
    # may need to be changed 
    y=1-(np.log10(FAR.value)-logFARmin)/(logFARmax-logFARmin)
    if isinstance(y, np.ndarray):
        y[y>1]=1
        y[y<alphamin]=alphamin
    else:
        y=min(y,1)
        y=max(y,alphamin)
    return y

# whether or not we also plot GW170817 for comparison
gw170817=True

# https://en.wikipedia.org/wiki/List_of_gravitational_wave_observations#Observation_candidates_from_O3/2019

filename = '../data/O3/gw_O3.dat'
data=Table.read(filename,format='ascii.csv',comment='#')
t=Time(data['Detectiontime (UTC)'])
name=data['GW event']
area=data['Locationarea[n 13](deg2)']*u.deg**2
D_L=np.array([float(x.split()[0]) for x in data['Luminositydistance(Mpc)[n 14]']])*u.Mpc
FAR=data['False Alarm Rate (Hz)']*u.Hz
BNS=data['BNS']
BHNS=data['NSBH']
BBH=data['BBH']
MassGap=data['Mass Gap']
Terr=data['Terr']

probs=np.c_[data['BNS'].data,data['NSBH'].data,data['BBH'].data,data['Mass Gap'].data,data['Terr'].data]
probs=(probs.T/probs.sum(axis=1)).T
# colors for BNS, BHNS, BBH, Mass Gap, Terrestrial
colors=['b','r','g','c','k']

size=get_size(area.value)

plt.figure()
fig, ax = plt.subplots(figsize=(8, 5))
for i in range(len(data)):
    drawPieMarker(t[i].plot_date, D_L[i].value, probs[i], size[i], colors, get_alpha(FAR[i]), plt.gca())

    if (np.max(probs[i,0:2]) > np.max(probs[i,2:4])) or (i==42):

        if i == 24:
            scale = -0.2
        elif i == 25:
            scale = 0.7
        elif i == 35:
            scale = 0.3
        elif i == 36:
            scale = -0.1
        elif i == 42:
            scale = 0.5
        elif i == 46:
            scale = 0.1
        elif np.mod(i,3) == 0:
            scale = 0.5
        elif np.mod(i,3) == 1:
            scale = 0.3
        else:
            scale = 0.1

        if i == 46:
            x, y = (t[i] + TimeDelta(30*u.day) + TimeDelta(1.3*np.log(size[i])*u.day)).plot_date, (1.0+scale)*D_L[i].value
        elif i == 42:
            x, y = (t[i] + TimeDelta(50*u.day) + TimeDelta(1.3*np.log(size[i])*u.day)).plot_date, (1.0+scale)*D_L[i].value
        else:
            x, y = (t[i] + TimeDelta(40*u.day) + TimeDelta(1.3*np.log(size[i])*u.day)).plot_date, (1.0+scale)*D_L[i].value

        print(i, name[i])

        plt.text(x,y,name[i],
             horizontalalignment='center',
             verticalalignment='top',fontsize=12)

        if i == 46:
            plt.annotate('', xy=(t[i].plot_date, D_L[i].value), xytext=(x*0.999969,y*0.95), color='k',
                         arrowprops=dict(facecolor='black', arrowstyle='->'),
                         )
        elif i == 42:
            plt.annotate('', xy=(t[i].plot_date, D_L[i].value), xytext=(x*0.999973,y*0.95), color='k',
                         arrowprops=dict(facecolor='black', arrowstyle='->'),
                         )
        else:
            plt.annotate('', xy=(t[i].plot_date, D_L[i].value), xytext=(x*0.99996,y*0.95), color='k',
                         arrowprops=dict(facecolor='black', arrowstyle='->'),
                         )

#plt.gcf().autofmt_xdate()
#plt.gca().set_xlim([Time('2019-04-01').datetime,Time('2020-05-01').datetime])
plt.gca().set_yscale('log')


for c,label in zip(colors,['BNS','NSBH','BBH','Mass Gap','Terrestrial']):
    plt.scatter(Time('2019-12-01').plot_date, 1,c=c, edgecolors='none',label=label)    
if gw170817:
    plt.gca().set_ylim([10,6e3])
else:
    plt.gca().set_ylim([100,6e3])

# these are the area -> size legend
areas=[10e3,5e3,1e3,5e2,1e2]
y=20
for i,a in zip(range(len(areas)),areas):
    # the date here may need to be shifted as more events are found
    plt.scatter(Time('2019-04-15').plot_date,y,c='r',s=get_size(a),
                alpha=0.5,edgecolors='none',)
    # and same here
    plt.text(Time('2019-04-25').plot_date,y,'%d deg$^2$' % a,
             verticalalignment='center',fontsize=10)
    y*=1.3

if gw170817:
    plt.scatter(Time('2019-08-17').plot_date,40,c='b',s=40e3/40,alpha=1,edgecolors='none')
    plt.text(Time('2019-08-17').plot_date,20,'GW170817+2yr',
             horizontalalignment='center',
             verticalalignment='top',fontsize=12)

ax.get_xaxis().set_major_locator(mdates.MonthLocator(interval=2))
ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%b %Y"))
fig.autofmt_xdate()
    
plt.ylabel('Luminosity Distance (Mpc)',fontsize=16)
plt.legend(loc=4)

plotName = '../output/O3_distance_vs_time.pdf'

plt.savefig(plotName, bbox_inches='tight')
plt.close()

names = []
dates = []
idxs = []
for i in range(len(data)):
    if (np.max(probs[i,0:2]) > np.max(probs[i,2:4])) or (i==42):
        names.append(name[i])
        dates.append(t[i].isot.split(".")[0])
        idxs.append(i)
dates = [datetime.strptime(ii, "%Y-%m-%dT%H:%M:%S") for ii in dates]

plt.figure()
bns_levels = np.array([5, 3, 1, 5, 1, 5, 1, 5, 1, 5, 3, 1, 5, 3, 1])
nsbh_levels = np.array([-5, -3, -5, -4, -3, -1, -5, -3, -1, -5, -3, -1, -5, -3, -1])

fig, ax = plt.subplots(figsize=(8, 5))

# Create the base line
start = min(dates)
stop = max(dates)
ax.plot((start, stop), (0, 0), 'k', alpha=.5)

bns_count, nsbh_count = 0, 0

color1 = 'cornflowerblue'
color2 = 'coral'
color3 = 'darkgreen'
color4 = 'pink'
color5 = 'cyan'

names2 = ["GW190425", "GW190814", "S200105ae", "S200115j"]
colors2 = [color1, color2, color3, color4, color5]

# Iterate through releases annotating each one
for ii, (iname, idate, idx) in enumerate(zip(names, dates, idxs)):
    idy = np.argmax(probs[idx][0:4])
    color = colors[idy]

    tt = Time(idate, format='datetime')
    if idx == 40:
        lab = tt + TimeDelta(5*u.day)
    elif idx == 42:
        lab = tt + TimeDelta(5*u.day)
    else:
        lab = tt + TimeDelta(-5*u.day)
    lab = lab.datetime

    if (idy == 0) or (idx==42):
        level = bns_levels[bns_count]
        bns_count = bns_count + 1
    else:
        level = nsbh_levels[nsbh_count]
        nsbh_count = nsbh_count + 1

    vert = 'top' if level < 0 else 'bottom'

    if iname in names2:
        idz = names2.index(iname)
        facecolor = colors2[idz]
        ax.scatter(idate, level, s=100, facecolor=facecolor,
                   edgecolor='k', zorder=9999)
    else:
        ax.scatter(idate, level, s=100, facecolor='w', edgecolor='k', zorder=9999)
    # Plot a line up to the text
    ax.plot((idate, idate), (0, level), c=color, alpha=.7)
    # Give the text a faint background and align it properly
    if idx  == 40:
        ax.text(lab, level, iname,
            horizontalalignment='left', verticalalignment=vert, fontsize=14,
            backgroundcolor=(1., 1., 1., .3))
    elif idx == 42:
        ax.text(lab, level, iname,
            horizontalalignment='left', verticalalignment=vert, fontsize=14,
            backgroundcolor=(1., 1., 1., .3))
    else:
        ax.text(lab, level, iname,
            horizontalalignment='right', verticalalignment=vert, fontsize=14,
            backgroundcolor=(1., 1., 1., .3))

props = dict(boxstyle='round', facecolor=color1, alpha=0.5)
# place a text box in upper left in axes coords
textstr = "Second confirmed \n BNS event"
ax.text(0.085, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

props = dict(boxstyle='round', facecolor=color2, alpha=0.5)
# place a text box in upper left in axes coords
textstr = "First high \n significance \n BHNS$^*$ event"
ax.text(0.20, 0.45, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

props = dict(boxstyle='round', facecolor=color3, alpha=0.5)
# place a text box in upper left in axes coords
textstr = "First high \n significance \n BHNS candidate \n with counterpart \n initially probable"
ax.text(0.57, 0.32, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

props = dict(boxstyle='round', facecolor=color4, alpha=0.5)
# place a text box in upper left in axes coords
textstr = "First high significance \n candidate with NS and \n object in the Mass Gap"
ax.text(0.52, 0.92, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

# Set the xticks formatting
# format xaxis with 3 month intervals
ax.get_xaxis().set_major_locator(mdates.MonthLocator(interval=2))
ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%b %Y"))
fig.autofmt_xdate()

# Remove components for a cleaner look
plt.setp((ax.get_yticklabels() + ax.get_yticklines() +
          list(ax.spines.values())), visible=False)
plt.show()

plotName = '../output/O3_timeline.pdf'
plt.savefig(plotName, bbox_inches='tight')
plt.close()

