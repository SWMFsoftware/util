#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-

def download_ADAPT_magnetogram():
    '''
    This routine reads the PARAM.in and download the corresponding ADAPT GONG
    magnetogram from #STARTTIME
    '''
    
    from ftplib import FTP
    import gzip
    import shutil
    import math
    import sys

    # Ensure that we are using a version of python >= 3
    if sys.version_info < (3,0):
        print('ERROR: Python version must be >= 3')
        print('Current version: '+sys.version)
        exit(-1)

    yyyy = int(input('Enter year: ' ))
    mm   = int(input('Enter month: '))
    dd   = int(input('Enter day: '  ))
    hh   = int(input('Enter hour: ' ))

    StrTypeMap = input('Type of ADAPT maps: fixed or central?  ')

    if StrTypeMap == 'fixed':
        TypeMap = 0
    elif StrTypeMap == 'central':
        TypeMap = 1
    else:
        print('Unrecognized type of ADAPT map')
        return(-0)

    # ADAPT maps only contains the hours for even numbers
    if hh%2 != 0:
        hh = math.floor(hh/2)*2
        print(' Warning: Hour must be an even number.'\
                  +' The entered hour value is changed to ', hh)

    StringTime = str(yyyy).zfill(4)+'-'+str(mm).zfill(2)+'-' \
        +str(dd).zfill(2)+'T'+str(hh).zfill(2)

    # Go to the the ADAPT ftp server
    ftp=FTP('gong2.nso.edu')
    ftp.login()
    
    # Only ADAPT GONG is considered
    ftp.cwd('adapt/maps/gong')

    # Go to the specific year
    try:
        ftp.cwd(str(yyyy))
    except:
        print('Cannot go to the specific year directory')
        return(-1)

    print('')

    # Only consider the public (4) Carrington Fixed (0) GONG (3) ADAPT maps
    patten = 'adapt4'+str(TypeMap)+'3*' + str(yyyy).zfill(4) + \
        str(mm).zfill(2) + str(dd).zfill(2)  + str(hh).zfill(2) + '*'
    
    print('Trying to download the', StrTypeMap, 'ADAPT map', \
              ' for date:',StringTime)
    # print('The patten is:', patten)
    
    filenames = ftp.nlst(patten)
    
    if len(filenames) < 1:
        print('Could not find a file that matches the patten')
        return(-2)
    
    for filename in filenames:
        # open the file locally
        fhandle=open(filename, 'wb')
        
        # try to download the magnetogram
        try:
            ftp.retrbinary('RETR '+ filename, fhandle.write)
        except:
            print('Cannot download ', filename)
            return(-3)
        
        # close the file
        fhandle.close()
        print('Downloaded:',filename)

        #unzip the file
        if '.gz' in filename:
            print('Unzip',filename)
            filename_unzip = filename.replace('.gz','')
            with gzip.open(filename, 'rb') as s_file, \
                open(filename_unzip, 'wb') as d_file:
                    shutil.copyfileobj(s_file, d_file, 65536)
    
    ftp.quit()

####### Script:

download_ADAPT_magnetogram()
