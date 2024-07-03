import requests
import json
import matplotlib.pyplot as plt
import numpy as np
import datetime
import prepare_report as pr
import jd_to_date as jd2date
import sys,getopt,argparse
import query_tns


#parser=argparse.ArgumentParser(description='ZTF name')
#parser.add_argument('ztfstring',metavar='N',type=str,help='ZTF name eg. ZTF18aajpkwd')
#args=parser.parse_args()
#print args.ztfstring
#ZTF_name=args.ztfstring


def main(argv):
    
    
    ZTF_name=str(sys.argv[1]) 
    reporter=sys.argv[2]   
 
    for arg in sys.argv:
        print(arg)
   
    username = "toomarshal"
    password = "findgws" 
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_programs.cgi', auth=(username, password))
    programs = json.loads(r.text)
    programidx = -1
    print("Programs you are a member of:")
    for index, program in enumerate(programs):
       if program['name'] == "Electromagnetic Counterparts to Gravitational Waves":
          programidx = program['programidx']
          print(program['programidx'])
          print(program['name'])
    submittedAnnotations = []
    
    mag,mager,time,filt,lim_mag,mag_last,programid=[],[],[],[],[],[],[]
    
    if programidx > 0:
      r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi', auth=(username, password), data={'programidx' : str(programidx)})
      sources = json.loads(r.text)
      #print sources
      #print "Sources in ZTF Science Validation:"
    
      for source in sources:
          if source['name']==ZTF_name:
              r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/source_summary.cgi', auth=(username, password),data={'sourceid' : str(source['id'])})
              sourceDict = json.loads(r.text)
              #f = open('temp.json', 'w+')
              #f.write(r.text)
              #f.close()
              #dict_index = (r.text).find('{') 
              #sourceDict = json.loads(r.text[dict_index:])
              #print sourceDict
              for i in range (0,len(sourceDict['uploaded_photometry'])):
                  #if sourceDict['uploaded_photometry'][i]['magpsf']!=99.0:
                  mag.append(sourceDict['uploaded_photometry'][i]['magpsf'])
                  mager.append(sourceDict['uploaded_photometry'][i]['sigmamagpsf'])
                  time.append(sourceDict['uploaded_photometry'][i]['jd'])
                  filt.append(sourceDict['uploaded_photometry'][i]['filter'])
                  lim_mag.append(sourceDict['uploaded_photometry'][i]['limmag'])
                  programid.append(sourceDict['uploaded_photometry'][i]['programid'])

    ind=np.argsort(time)
    time=np.asarray(time)[ind]
    mag=np.asarray(mag)[ind]
    mager=np.asarray(mager)[ind]
    filt=np.asarray(filt)[ind]
    lim_mag=np.asarray(lim_mag)[ind]
    programid=np.asarray(programid)[ind]

    print(programid)
    
    #pi_filt=(programid==2 or programid==1);
    
    
    #time=np.asarray(time)[pi_filt]
    #mag=np.asarray(mag)[pi_filt]
    #mager=np.asarray(mager)[pi_filt]
    #filt=np.asarray(filt)[pi_filt]
    #lim_mag=np.asarray(lim_mag)[pi_filt]
    #programid=np.asarray(programid)[pi_filt]
    
    print(programid)
  
    #print time
    
    up_ind=[]
    
    for i in range(0,len(mag)):
        if mag[i]<99:
            up_ind.append(i)
            
   

    if len(lim_mag)>0 and up_ind[0]>0:
        if time[up_ind[0]]-time[up_ind[0]-1] > 0.01:
           lim_mag=lim_mag[up_ind[0]-1]
           time_lim=time[up_ind[0]-1]
           filt_lim=filt[up_ind[0]-1]
        else:
           lim_mag=99
           time_lim=99
           filt_lim='r'
    else: 
        lim_mag=99
        time_lim=99
        filt_lim='r'
     
        
    #print up_ind

    
    #mag_unfiltered=mag
    #mager=np.asarray(mager)[ind]
    #filt=np.asarray(filt)[ind]
    #time_lim_mag=time
    
    
    #time=np.asarray(time)[ind]
    
    
    
    #ind2=np.argsort(time)
    #lim_mag_sorted=np.asarray(lim_mag)[ind]
    #time_lim_mag=np.asarray(time_lim_mag)[ind]
    
    '''
    if lim_mag_sorted!=-99:
        lim_mag=lim_mag_sorted
        
    if lim_mag_sorted==99:
        lim_mag='non detection'
    '''   
       
    #lim_mag=str("non_detection": {"archiveid": "0","archival_remarks": "Non existent in SDSS/PS1"},)
        
    
    
    #mager=np.asarray(mager)[ind]
    #filt=np.asarray(filt)[ind]
    #lim=np.asarray(lim)[ind]
    
    
                  
    
    
    
    #print sourceDict['ra'],sourceDict['dec']
    #print filt[ind],mag[ind],time[ind]
    
    mag_last=mag[len(mag)-1]
    magerr_last=mager[len(mag)-1]

    time_last=jd2date.jd_to_date(time[len(mag)-1])
    filt_last=filt[len(mag)-1]
    time_last=datetime.datetime(time_last[0],time_last[1],time_last[2],time_last[3],time_last[4],int(time_last[5]))

    obsdate = (jd2date.jd_to_date(time[up_ind[0]]))
    obsdate=datetime.datetime(obsdate[0],obsdate[1],obsdate[2],obsdate[3],obsdate[4],int(obsdate[5]))
    if time_lim == 99:
        time_lim=datetime.datetime(9999,9,9,9,9,9)
    else:
        time_lim=(jd2date.jd_to_date(time_lim))
        time_lim=datetime.datetime(time_lim[0],time_lim[1],time_lim[2],time_lim[3],time_lim[4],int(time_lim[5]))
    ra = sourceDict['ra']
    dec = sourceDict['dec']
    mag = mag[up_ind[0]]
    mag_err=mager[up_ind[0]]
    filter_name = filt[up_ind[0]]
    
    internal_name=ZTF_name
    #reporter=str('on behalf of the ZTF collaboration')
    print(internal_name, ra, dec, mag, mag_err,lim_mag,time_lim, filter_name,filt_lim, obsdate)
    
    #print internal_name+','+'2018XX'+','+str(ra)+','+str(dec)+','+str(obsdate)+','+str(mag)+' '+str(filter_name)+','+str(time_lim)+','+str(lim_mag)+'\n'
    #f = open('../notns_data.txt', 'w+')
    #f.write(internal_name+','+'2018XX'+','+str(ra)+','+str(dec)+','+str(obsdate)+','+str(mag)+','+str(filter_name)+','+str(time_lim)+','+str(lim_mag)+','+str(filt_lim)+','+str(time_last)+','+str(mag_last)+','+str(filt_last)+'\n')
    #f.close()



    tnsname=''
    query_tns.api_key="54916f1700966b3bd325fc1189763d86512bda1d"
    query_tns.url_tns_api="https://wis-tns.weizmann.ac.il/api/get"        
    query_tns.search_obj=[("ra",str(ra)), ("dec",str(dec)), ("radius","5"), ("units","arcsec")]
    response=query_tns.search(query_tns.url_tns_api,query_tns.search_obj)
    json_data=json.loads(response.text)
    if len(json_data['data']['reply']) > 0:
        tnsname=(json_data['data']['reply'][0]['prefix']+json_data['data']['reply'][0]['objname'])    
    #print json_data      
    query_tns.get_obj=[("objname",tnsname[2:len(tnsname)]), ("photometry","1"), ("spectra","0")]  
    response=query_tns.get(query_tns.url_tns_api,query_tns.get_obj)
    json_data2=json.loads(response.text)
    print(json_data2)
    str2=str(json_data2);
    ztfcheck=str2.find('ZTF');

    if ztfcheck == -1:
        pr.report(internal_name, ra, dec, mag, mag_err,lim_mag,time_lim, filter_name,filt_lim, obsdate, time_last, mag_last, filt_last, magerr_last, reporter)
    #    print('kek')
    else:
        print('Already reported by ZTF')

    #object_flag, object_name = pr.report(internal_name, ra, dec, mag, mag_err, filter_name, obsdate)
    #pr.report(internal_name, ra, dec, mag, mag_err,lim_mag,time_lim, filter_name,filt_lim, obsdate, time_last, mag_last, filt_last, magerr_last)


if __name__=='__main__':
    main(sys.argv[1:])
