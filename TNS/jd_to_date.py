#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import math
import datetime


def jd_to_date(jd):
    """
    Convert Julian Day to date.

    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
        
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    
    """
    jd = jd + 0.5
    
    F, I = math.modf(jd)
    I = int(I)
    
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    
    D = math.trunc((C - 122.1) / 365.25)
    
    E = math.trunc(365.25 * D)
    
    G = math.trunc((C - E) / 30.6001)
    
    day = C - E + F - math.trunc(30.6001 * G)
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
    
    number_dec = float(str(day-int(day))[1:])   
    
    
    t=86400*number_dec
    
    hour=t/3600.
    hour_round=int(hour)
                   
    minute=(hour-int(t/3600))*60
    minute_round=int(minute)
    
    second=(minute-int(minute))*60
    
    
    #print t.date 
    #t2=datetime.datetime(t)
    #print t,type(t)
    #print t2
    
    #print t.days
    #print t.seconds
    
    
    
    #hour=t.hour
    #minute=t.minute
    #second=t.second
    
    #int(t[0:2])
    #minute=int(t[3:5])
    #second=float(t[6:-1])
    
    #hour, minute,second
        
    return year, month, int(day) ,hour_round, minute_round,second


