#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 12:02:37 2017

Developed and tested in : 
- Python version 3.6.3
- Linux Ubuntu version 16.04 LTS (64-bit)

@author: Nikola Knezevic
"""

import os
import requests
import json
from collections import OrderedDict
import time    
import sys

############################# PARAMETERS #############################
api_key="54916f1700966b3bd325fc1189763d86512bda1d"                   #
list_of_filenames="Here put your list of filenames for uploading."   #
report_filename="Here put your report filename."                     #
# report type can only be "tsv" or "json"                            #
report_type="Here put the type of your report."                      #
id_report="Here put your report ID for getting report's reply."      #
######################################################################

#############################    URL-s   #############################
# url of TNS and TNS-sandbox api                                     #
url_tns_api="https://wis-tns.weizmann.ac.il/api"                     #
url_tns_sandbox_api="https://sandbox-tns.weizmann.ac.il/api"         #
######################################################################

############################# DIRECTORIES ############################
# current working directory                                          #
cwd=os.getcwd()                                                      #
# folder containing files for uploading                              #
upload_folder=os.path.join(cwd, 'files_for_uploading')               #
# folder containing tsv reports for sending                          #
tsv_reports_folder=os.path.join(cwd, 'tsv_reports_for_sending')      #
# folder containing json reports for sending                         #
json_reports_folder=os.path.join(cwd, '')    #
######################################################################

########################## API FUNCTIONS #############################
# function for changing data to json format                          #
def format_to_json(source):                                          #
    # change data to json format and return                          #
    parsed=json.loads(source,object_pairs_hook=OrderedDict)          #
    result=json.dumps(parsed,indent=4)                               #
    return result                                                    #
#--------------------------------------------------------------------#
# function for uploading files trough api                            #
def upload_files(url,list_of_files):                                 #
  try:                                                               #
    # url for uploading files                                        #
    upload_url=url+'/file-upload'                                    #
    # construct the list of (key,value) pairs                        #
    files=[('api_key',(None, api_key))]                              #
    # go trough every file in the list_of_files and add that         #
    # to the list of (key,value) pairs for uploading                 #
    for i in range(len(list_of_files)):                              #
        file_name=list_of_files[i]                                   #
        file_path=os.path.join(upload_folder,file_name)              #
        key='files['+str(i)+']'                                      #
        value=(file_name,open(file_path,'rb'))                       #
        files.append((key,value))                                    #
    # upload all files using request module                          #
    response=requests.post(upload_url, files=files)                  #
    # return response                                                #
    return response                                                  #
  except Exception as e:                                             #
    return [None,'Error message : \n'+str(e)]                        #
#--------------------------------------------------------------------#
# function for sending tsv reports (AT or Classification)            #
def send_tsv_report(url,tsv_report):                                 #
  try:                                                               #
    # url for sending tsv reports                                    #
    tsv_url=url+'/csv-report'                                        #
    # tsv report file path                                           #
    tsv_file_path=os.path.join(tsv_reports_folder,tsv_report)        #
    # construct list of (key,value) pairs                            #
    tsv_data=[('api_key',(None, api_key)),                           #
              ('csv',(tsv_report,open(tsv_file_path,'rb')))]         #
    # send tsv report using request module                           #
    response=requests.post(tsv_url, files=tsv_data)                  #
    # return response                                                #
    return response                                                  #
  except Exception as e:                                             #
    return [None,'Error message : \n'+str(e)]                        #
#--------------------------------------------------------------------#
# function for sending json reports (AT or Classification)           #
def send_json_report(url,json_report):                               #
  try:                                                               #
    # url for sending json reports                                   #
    json_url=url+'/bulk-report'                                      #
    # json report file path                                          #
    json_file_path=os.path.join(json_reports_folder,json_report)     #
    # read json data from file                                       #
    json_read=format_to_json(open(json_file_path).read())            #
    # construct list of (key,value) pairs                            #
    json_data=[('api_key',(None, api_key)),                          #
               ('data',(None,json_read))]                            #
    # send json report using request module                          #
    response=requests.post(json_url, files=json_data)                #
    # return response                                                #
    return response                                                  #
  except Exception as e:                                             #
    return [None,'Error message : \n'+str(e)]                        #
#--------------------------------------------------------------------#
# function for getting reply from report                             #
def reply(url, report_id):                                           #
  try:                                                               #
    # url for getting report reply                                   #
    reply_url=url+'/bulk-report-reply'                               #
    # construct list of (key,value) pairs                            #    
    reply_data=[('api_key',(None, api_key)),                         #
                ('report_id',(None,report_id))]                      #
    # send report ID using request module                            #
    response=requests.post(reply_url, files=reply_data)              #
    # return response                                                #
    return response                                                  #
  except Exception as e:                                             #
    return [None,'Error message : \n'+str(e)]                        #
######################################################################

################### FUNCTIONS FOR HANDLING OUTPUT #################### 
# http errors                                                        #
http_errors = {                                                      #                          
304: 'Error 304: Not Modified: There was no new data to return.',    #
400: 'Error 400: Bad Request: The request was invalid. '\
     'An accompanying error message will explain why.',              #
403: 'Error 403: Forbidden: The request is understood, but it has '\
     'been refused. An accompanying error message will explain why.',#
404: 'Error 404: Not Found: The URI requested is invalid or the '\
     'resource requested, such as a category, does not exists.',     #
500: 'Error 500: Internal Server Error: Something is broken.',       #
503: 'Error 503: Service Unavailable.'                               #
}                                                                    #                                                                
#--------------------------------------------------------------------#
# function that checks response and                                  #
# returns True if everything went OK                                 #
# or returns False if something went wrong                           #
def check_response(response):                                        #
    # if response exists                                             #
    if None not in response:                                         #
        # take status code of that response                          #
        status_code=int(response.status_code)                        #
        if status_code==200:                                         #
            # response as json data                                  #
            json_data=response.json()                                #
            # id code                                                #
            id_code=str(json_data['id_code'])                        #
            # id message                                             #
            id_message=str(json_data['id_message'])                  #
            # print id code and id message                           #
            print ("ID code = "+id_code)                             #
            print ("ID message = "+id_message)                       #
            # check if id code is 200 and id message OK              #
            if (id_code=="200" and id_message=="OK"):                #
                return True                                          #
            #special case                                            #
            elif (id_code=="400" and id_message=="Bad request"):     #
                return None                                          #
            else:                                                    #
                return False                                         #
        else:                                                        #
            # if status code is not 200, check if it exists in       #
            # http errors                                            #
            if status_code in list(http_errors.keys()):              #
                print (list(http_errors.values())                    #
                       [list(http_errors.keys()).index(status_code)])#
            else:                                                    #
                print ("Undocumented error.")                        #
            return False                                             #
    else:                                                            #
        # response doesn't exists, print error                       #
        print (response[1])                                          #
        return False                                                 #
#--------------------------------------------------------------------#
# find all occurrences of a specified key in json data               #
# and return all values for that key                                 # 
def find_keys(key, json_data):                                       #
    if isinstance(json_data, list):                                  #
        for i in json_data:                                          #
            for x in find_keys(key, i):                              #
               yield x                                               #
    elif isinstance(json_data, dict):                                #
        if key in json_data:                                         #
            yield json_data[key]                                     #
        for j in list(json_data.values()):                           #
            for x in find_keys(key, j):                              #
                yield x                                              #
#--------------------------------------------------------------------#
# print feedback                                                     #
def print_feedback(json_feedback):                                   #
    # find all message id-s in feedback                              #
    message_id=list(find_keys('message_id',json_feedback))           #
    # find all messages in feedback                                  #
    message=list(find_keys('message',json_feedback))                 #
    # find all obj names in feedback                                 #
    objname=list(find_keys('objname',json_feedback))                 #
    # find all new obj types in feedback                             #
    new_object_type=list(find_keys('new_object_type',json_feedback)) #
    # find all new obj names in feedback                             #
    new_object_name=list(find_keys('new_object_name',json_feedback)) #
    # find all new redshifts in feedback                             #
    new_redshift=list(find_keys('new_redshift',json_feedback))       #
    # index counters for objname, new_object_type, new_object_name   #
    # and new_redshift lists                                         #
    n_o=0                                                            #
    n_not=0                                                          #
    n_non=0                                                          #
    n_nr=0                                                           #
    # messages to print                                              #
    msg=[]                                                           #
    # go trough every message and print                              #
    for j in range(len(message)):                                    #
        m=str(message[j])                                            #
        m_id=str(message_id[j])                                      #
        if m_id not in ['102','103','110']:                          #
            if m.endswith('.')==False:                               #
                m=m+'.'                                              #
            if m_id=='100' or  m_id=='101':                          #
                m="Message = "+m+" Object name = "+str(objname[n_o]) #
                n_o=n_o+1                                            #
            elif m_id=='120':                                        #  
                m="Message = "+m+" New object type = "+str(          # 
                                            new_object_type[n_not])  #
                n_not=n_not+1                                        #
            elif m_id=='121':                                        #
                m="Message = "+m+" New object name = "+str(          #
                                             new_object_name[n_non]) #
                n_non=n_non+1                                        #
            elif m_id=='122' or  m_id=='123':                        #
                m="Message = "+m+" New redshift = "+str(             #
                                                 new_redshift[n_nr]) #
                n_nr=n_nr+1                                          #
            else:                                                    #
                m="Message = "+m                                     #
            msg.append(["Message ID = "+m_id,m])                     #
    # return messages                                                #        
    return msg                                                       #
#--------------------------------------------------------------------#
# sending report id to get reply of the report                       #
# and printing that reply                                            #
def print_reply(url,report_id):                                      #
    # sending reply using report id and checking response            #
    print ("Sending reply for the report id "+report_id+" ...")      #
    reply_res=reply(url, report_id)                                  #
    reply_res_check=check_response(reply_res)                        #
    # if reply is sent                                               #
    if reply_res_check==True:                                        #
        print ("The report was successfully processed on the TNS.\n")#
        # reply response as json data                                #
        json_data=reply_res.json()  
        # feedback of the response                                   #
        feedback=list(find_keys('feedback',json_data))               #
        # check if feedback is dict or list                          #
        if type(feedback[0])==type([]):                              #
            feedback=feedback[0]                                     #
        # go trough feedback                                         #
        for i in range(len(feedback)):                               #
            # feedback as json data                                  #
            json_f=feedback[i]                                       #
            # feedback keys                                          #
            feedback_keys=list(json_f.keys())                        #
            # messages for printing                                  #
            msg=[]                                                   #
            # go trough feedback keys                                #
            for j in range(len(feedback_keys)):                      #
                key=feedback_keys[j]                                 #
                json_feed=json_f[key]                                #
                msg=msg+print_feedback(json_feed)                    # 
            if msg!=[]:                                              #
                print ("-----------------------------------"\
                       "-----------------------------------" )       #
                for k in range(len(msg)):                            #
                    print (msg[k][0])                                #
                    print (msg[k][1])                                #
                print ("-----------------------------------"\
                       "-----------------------------------\n")      # 
    else:                                                            #
        if (reply_res_check!=None):                                  #
            print ("The report doesn't exist on the TNS.")           #
        else:                                                        #
            print ("The report was not processed on the TNS "\
                   "because of the bad request(s).")                 #
#--------------------------------------------------------------------#
# how many second to sleep                                           #
SLEEP_SEC = 1                                                        #
# max number of time to check response                               #
LOOP_COUNTER = 60                                                    #
# keeping sys.stdout                                                 #
old_stdout = sys.stdout                                              #
# Disable print                                                      #
def blockPrint():                                                    #
    sys.stdout = open(os.devnull, 'w')                               #
# Restore print                                                      #
def enablePrint():                                                   #
    sys.stdout.close()                                               #
    sys.stdout = old_stdout                                          #
#--------------------------------------------------------------------#
# sending tsv or json report (at or class) and printing reply        #
def send_report(url, report, type_of_report):                        #
    # sending report and checking response                           #
    print ("Sending "+report+" to the TNS...")                       #
    # choose which function to call                                  #
    if type_of_report=="tsv":                                        #
        response=send_tsv_report(url,report)                         #
    else:                                                            #
        response=send_json_report(url,report)                        #
    response_check=check_response(response)                          #
    # if report is sent                                              #
    if response_check==True:                                         #
        print ("The report was sent to the TNS.")                    # 
        # report response as json data                               #
        json_data=response.json()                                    #
        # taking report id                                           #
        report_id=str(json_data['data']['report_id'])                #
        print ("Report ID = "+report_id)                             #
        print ("")                                                   #
        # sending report id to get reply of the report               #
        # and printing that reply                                    #
        # waiting for report to arrive before sending reply          #
        # for report id                                              #
        blockPrint()                                                 #
        counter = 0                                                  #
        while True:                                                  #
            time.sleep(SLEEP_SEC)                                    #
            reply_response=reply(url,report_id)                      #
            reply_res_check=check_response(reply_response)           #
            if reply_res_check!=False or counter >= LOOP_COUNTER:    #
                break                                                #
            counter += 1                                             #
        enablePrint()                                                #
        print_reply(url,report_id)                                   #
    else:                                                            #
        print ("The report was not sent to the TNS.")                #
#--------------------------------------------------------------------#
# uploading files and printing reply                                 #
def upload(url, list_of_files):                                      #
    # upload files and checking response                             #
    print ("Uploading files on the TNS...")                          #
    response=upload_files(url,list_of_files)                         #
    response_check=check_response(response)                          #
    # if files are uploaded                                          #
    if response_check==True:                                         #
        print ("The following files are uploaded on the TNS : ")     #
        # response as json data                                      #
        json_data=response.json()                                    #
        # list of uploaded files                                     #
        uploaded_files=json_data['data']                             #
        for i in range(len(uploaded_files)):                         #
            print ("filename : "+str(uploaded_files[i]))             #
    else:                                                            #
        print ("Files are not uploaded on the TNS.")                 #
######################################################################

# EXAMPLE

# API key of your Bot:
api_key="54916f1700966b3bd325fc1189763d86512bda1d"

report_filename="data.txt"
report_type="json"
send_report(url_tns_api,report_filename,report_type)


# Comment/Uncomment sections for testing the various examples:

"""

# ---------------------------------------------------
# upload files
list_of_filenames=["rel_file_1.png","rel_file_2.jpg",
                   "spectra_2016jhs.asci.txt","spectra_2016jhs.fits"]
upload(url_tns_sandbox_api,list_of_filenames)


# ---------------------------------------------------
# send AT report
report_filename="json_at_report.txt"
report_type="json" 
# OR 
report_filename="bulk_tsv_at_report.txt"
report_type="tsv"
send_report(url_tns_sandbox_api,report_filename,report_type)


# ---------------------------------------------------
# send Classification report
report_filename="json_classification_report.txt"
report_type="json"
# OR 
report_filename="bulk_tsv_classification_report.txt"
report_type="tsv"
send_report(url_tns_sandbox_api,report_filename,report_type)

# ---------------------------------------------------
# reply
id_report="5775"
print_reply(url_tns_sandbox_api,id_report)

"""



