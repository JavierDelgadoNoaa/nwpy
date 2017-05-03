import time

'''
This module contains function to convert different time representations
'''

def yyyymmddHHMMSS_to_epoch(date_string):
    ''' Convert time string in yyyymmddHHMMSS to seconds since epoch'''
    if len(date_string) != 14: raise Exception('Invalid input to yyyymmddHHMMSS')
    return int( time.mktime( time.strptime(date_string, '%Y%m%d%H%M%S') ) )

def yyyymmddHHMM_to_epoch(date_string):
    ''' Convert time string in yyyymmddHHMM to seconds since epoch'''
    if len(date_string) != 12: raise Exception('Invalid input to yyyymmddHHMM')
    return int( time.mktime( time.strptime(date_string, '%Y%m%d%H%M') ) )

def yyyymmddHH_to_epoch(date_string):
    ''' Convert time string in yyyymmddHH to seconds since epoch'''
    s = '%s00' %date_string
    return yyyymmddHHMM_to_epoch(s)

def epoch_to_pretty_time_string(epoch_date):
    ''' Convert time in seconds since epoch to a date in MM/DD/YYYY @ hh:mm '''
    t = time.gmtime(epoch_date)
    return time.strftime('%m/%d/%Y @ %H:%M' , t)

def epoch_to_yyyymmddHHMM(epoch_date):
    t = time.gmtime(epoch_date)
    return time.strftime('%Y%m%d%H%M' , t)

