#!/usr/bin/env python

# this magnetogram remapping can either be run as a script from the unix command line
# or imported into Python. See description of use in the README file
#         by  Richard A. Frazin July 2014 - February 2015

# April 2020: 
# script updated to work for both reading & remapping fits map & use in EEGGL
# uses pyfits instead of astropy

import pyfits as fits
from scipy import interpolate
from scipy import integrate
import numpy as np
import sys
import time
import argparse
import pdb

def remap(inputfile, outputfile, nlat = -1, nlong = -1, out_grid = 'unspecified', i=-1, nSmooth=1):
    """
    Flux-conserving magnetogram remapping tool.
    inputfile - FITS file containing original magnetogram (include path)
    outputfile - contains result in customized output format (include path)
    nlat  (opitonal) = desired number of latitude  points in output,
         if not specified, will be set to the same as the input file
    nlong (optional) = desired number of longitude points in output
        if not specified, will be set to the same as the input file
    out_grid (optional), choices are 'sin(lat)' or 'uniform'
        if not specified, the output grid will be the same type
    If nlat, nlong and out_grid are ALL left out, no remapping is done,
        and the code simply reformats.
    Note that ADAPT files may have multiple maps.  Only the 1st is utilized unless specified otherwise.
    My MATLAB code remap_mag.m uses the same algorithm but runs much faster.
       by Richard Frazin, July 2014 - Feb 2015
    """
    pi = 3.141592653589793
    if ( (out_grid != 'sin(lat)') and (out_grid != 'uniform') and (out_grid != 'unspecified') ):
        print ("Unknown output grid type.  Choices are blank, 'unspecified', 'uniform' and 'sin(lat)' ")
        return(-1)
    
    cc =  FITS_RECOGNIZE(inputfile)

    if cc == -1:
        print ("Input file not recognized.")
        return(-1)
    else:
        magtype = cc[0]
        grid_type = cc[1]
        map_data = cc[2]
        input_grid = grid_type

    #what kind of transformation are we doing?        
    if out_grid == 'unspecified':
        out_grid = grid_type
    
    if grid_type == out_grid:
        transformation = 'rebin' #no change in grid type, so just rebin
    elif ( (grid_type == 'uniform') and (out_grid == 'sin(lat)') ):
        transformation = 'reg2sin'
    elif ( (grid_type == 'sin(lat)') and (out_grid == 'uniform') ):
        transformation = 'sin2reg'
    else:
        print ("Unknown transformation type.")
        return(-1)

    # read the data
    g = fits.open(inputfile)
    d = g[0].data

    nlo = g[0].header['NAXIS1'] # number of longitude points
    nla = g[0].header['NAXIS2'] #           latitude
    
    if nlat == -1:
        nlat = nla
    elif nlat < 1:
        print ("nlat has to be -1 or a positive integer.")
        return(-1)

    if nlong == -1:
        nlong = nlo
    elif nlong < 1:
        print ("nlong has to be -1 or a positive integer.")
        return(-1)
        
    if magtype == 'ADAPT Synchronic':
        nim = g[0].header['NAXIS3'] # number of images
        imdex = i  #which of the 12 maps do you want?
        print ('This file contains ', str(nim),' images. Writing out file ', i)
        if nim > 1:  #just keep one of them for now
            d = d[imdex,:,:]        

    # Conservative smoothing. Boundary condition:
    # Periodic in Longitude, reflect in Latitude.
    if (nSmooth>2):
        d=smooth(nlong,nlat,nSmooth,d)
   
    g[0].header['NAXIS1'] = nlong # new number of longitude points
    g[0].header['NAXIS2'] = nlat  #               latitude

    try:
        g[0].header['GRID'] = out_grid #change the existing header value
    except KeyError as er:
        g[0].header.set('GRID',out_grid) #create FITS header keyword

    try:
        g[0].header['CTYPE2'] = out_grid #change the existing header value
    except KeyError as er:
        g[0].header.set('CTYPE2',out_grid) #create FITS header keyword

    if out_grid == 'sin(lat)':
        newlat = (180/pi)*np.arcsin(np.linspace(-1. + 1./2/nlat,1. - 1./2/nlat,nlat)) # in deg
    elif out_grid == 'uniform':
        newlat = (180/pi)*np.linspace(-pi/2 + pi/2/nlat,pi/2 - pi/2/nlat,nlat) # in deg
    else:
        print ("out_grid incorrectly set.")
        return(-1)

    if ( (nlo == nlong) and (nla == nlat) and (grid_type == out_grid) ):
        newmap = d  #no remapping
    else:
        #first make a hybrid map that is (nla X nlong) by using the rebin
        # and add alg. in the longitude direction. If nlo = nlong, hybrid --> d
        hybrid = np.zeros([nla,nlong]) 
        crap = np.arange(nlo+1)
        for pf in crap[nlo+1:0:-1]: #pf will be the greatest common factor of nlong and nlo
            if ( (np.mod(nlo,pf) == 0) and (np.mod(nlong,pf) == 0)): #common factor test
                nlo_fac   = nlo/pf
                nlong_fac = nlong/pf
                break
        for k in np.arange(nla):
            w = np.kron(d[k,:],np.ones(int(nlong_fac))) #this array has length pf*nlo_fac*nlong_fac
            for l in np.arange(nlong):  #take the average over nlo_fac bins of w
                hybrid[k,l] = np.sum(w[int(l*nlo_fac):(l+1)*int(nlo_fac)])/nlo_fac

        newmap = np.zeros([nlat,nlong]) #output map                           
        if transformation == 'rebin':  #do rebin and add in the latitude direction, if nlo = nlat, newmap --> d
            crap = np.arange(nla+1)
            for pf in crap[nla+1:0:-1]: #pf will be the greatest common factor of nla and nlat
                if ( (np.mod(nla,pf) == 0) and (np.mod(nlat,pf) == 0) ): # common factor test
                    nla_fac  = nla/pf
                    nlat_fac = nlat/pf
                    break
                
            for k in np.arange(nlong):
                w = np.kron(hybrid[:,k].T,np.ones(int(nlat_fac)))#length is pf*nla_fac*nlat_fac
                for l in np.arange(nlat):
                    newmap[l,k] = np.sum(w[int(l*nla_fac):(l+1)*int(nla_fac)])/nla_fac

        elif transformation == 'reg2sin':
            oldlat =  np.linspace(-pi/2 + pi/2/nla,pi/2 - pi/2/nla,nla) #old latitude grid
            oldlat = np.hstack((-pi/2-1.e-9,oldlat,pi/2+1.e-9)) #for the interpolator
            bin_boundary = np.arcsin(np.linspace(-1.,1.,nlat+1))   #boundaries of new sin(latitude) grid
            for k in np.arange(nlong):   #the magnetic field value assigned is the flux divided by the area.  
                u = np.hstack((hybrid[0,k],hybrid[:,k],hybrid[nla-1,k]))
                crap = interpolate.interp1d(oldlat,u,kind='linear') #magnetic field interpolator
                fcn = lambda x : crap(x)*np.cos(x)  #this is B(theta)*cos(theta)
                for l in np.arange(nlat):
                    result = integrate.quad(fcn,bin_boundary[l],bin_boundary[l+1],epsabs=1.e-3,epsrel=1.e-3)/(np.sin(bin_boundary[l+1]) - np.sin(bin_boundary[l]))
                    newmap[l,k] = result[0]

        elif transformation == 'sin2reg':
            oldlat = np.arcsin(np.linspace(-1. + 1./2/nla,1. - 1./2/nla,nla)) #arcsin(old sin(latitude) grid)
            oldlat = np.hstack((-pi/2-1.e-9,oldlat,pi/2+1.e-9)) #for the interpolator
            bin_boundary = np.linspace(-pi/2,pi/2,nlat+1) #boundaries of new latitude grid
            #pdb.set_trace()
            for k in np.arange(nlong):   #the magnetic field value assigned is the flux divided by the area.  
                u = np.hstack((hybrid[0,k],hybrid[:,k],hybrid[nla-1,k]))
                crap = interpolate.interp1d(oldlat,u,kind='linear') #magnetic field interpolator
                fcn = lambda x : crap(x)*np.cos(x)  #this is B(theta)*cos(theta)
                for l in np.arange(nlat):
                    result = integrate.quad(fcn,bin_boundary[l],bin_boundary[l+1],epsabs=1.e-3,epsrel=1.e-3)/(np.sin(bin_boundary[l+1]) - np.sin(bin_boundary[l]))
                    newmap[l,k] = result[0]

        else:
            print ("Unknown transformation type.")
            return(-1)

    #test for flux conservation in the transformation        
    test_flux = False 
    if test_flux:
        if grid_type == 'uniform':
            latt =  np.cos(np.linspace(-pi/2 + pi/2/nla,pi/2 - pi/2/nla,nla))
            cosgrid = np.kron(latt,np.ones((nlo,1))).T
            oldflux = np.sum(np.multiply(cosgrid,d))*2.*pi*pi/nlo/nla
        elif grid_type == 'sin(lat)':
            oldflux = np.sum(d)*4.*pi/nlo/nla
        else:
            print ("Bad grid_type.")
            return(-1)
        if out_grid == 'uniform':
            latt =  np.cos(np.linspace(-pi/2 + pi/2/nlat,pi/2 - pi/2/nlat,nlat))
            cosgrid = np.kron(latt,np.ones((nlong,1))).T
            newflux = np.sum(np.multiply(cosgrid,newmap))*2.*pi*pi/nlong/nlat
        elif out_grid == 'sin(lat)':
            newflux = np.sum(newmap)*4.*pi/nlong/nlat
        else:
            print ("Bad out_grid.")
            return(-1)
        print ("original flux =",str(oldflux),", new flux =",str(newflux))
    
    #try to get some context information from the FITS file to include in output
    try:
        CRnumber = str(g[0].header['MAPCR'])  #works for ADAPT - Center of CM 
        CR = str(int(g[0].header['MAPCR']))  #works for ADAPT
    except KeyError as er:
        CRnumber = '0'

    if CRnumber == '0':
        try :
            CRnumber = str(g[0].header['CRCENTER']) #works on GONG and MDI
        except KeyError as er:
            CRnumber = '0'
        try :
            CR = str(g[0].header['CAR_ROT']) #works on GONG and MDI
        except KeyError as er:
            CR ='0'
    if CRnumber == '0':
        try :
            CRnumber = str(g[0].header['CAR_ROT'])
        except KeyError as er:
            CRnumber = '0'

    if magtype.find('GONG') > -1:
        mapdate = g[0].header['DATE'] #works for GONG
    else:                          
        try:
            mapdate = g[0].header['MAPTIME']  #works for ADAPT
        except KeyError as er:
            mapdate = '0000-00-00T00:00:00'

        if mapdate == '0000-00-00T00:00:00':    
            try :
                mapdate = g[0].header['MAP_DATE']  #works for Hathaway
            except KeyError as er:
                mapdate = '0000-00-00T00:00:00'

        if mapdate == '0000-00-00T00:00:00':    
            try :
                mapdate = g[0].header['T_OBS']  #works for MDI, HMI
            except KeyError as er:
                mapdate = '0000-00-00T00:00:00'
    try:
        bunit = g[0].header['BUNIT']  #works on GONG, MDI
    except KeyError as er:   #Hathaway and ADAPT don't list units
        bunit = 'Gauss'  #assume it's Gauss if you don't know
                
    try :
        long0 = g[0].header['LONG0'] #works on GONG 
    except KeyError as er:
        long0 = 0
    if long0 == 0:
        try :
            long0 = g[0].header['CRLNGEDG']   #works on ADAPT
        except KeyError as er:
            long0 = 0

    #ascii output file, Gabor format, the first line is arbitary

    fid = open(outputfile,'w')
    if magtype == 'ADAPT Synchronic':
        line0 = 'MagnetogramType = '+magtype+'; ADAPTRealization = ' \
            +str(imdex+1)+ '; InstrumentName = '+map_data+'; InputLatGrid = '\
            +input_grid+'; OutputLatGrid = '+out_grid+'; MagUnit = ' \
            +bunit+'; InputMapResolution = '+str(nlo)+','+str(nla)+'; MagnetogramDate = '\
            +mapdate+'; CentralMeridianLong = '+CRnumber+'; InputFile = ' \
            +str(inputfile)+'; ASCIIFileCreationDate = '+time.strftime("%Y-%m-%dT%H:%M:%S")+'\n' 
    else:
        line0 = 'MagnetogramType = '+magtype+'; InstrumentName = '+map_data+'; InputLatGrid = '\
            +input_grid+'; OutputLatGrid = '+out_grid+'; MagUnit = '+bunit+'; InputMapResolution = '\
            +str(nlo)+','+str(nla)+'; MagnetogramDate = '+mapdate+'; CentralMeridianLong = '\
            +CRnumber+'; InputFile = "'+str(inputfile)+'"; ASCIIFileCreationDate = '\
            +time.strftime("%Y-%m-%dT%H:%M:%S")+'\n' 

    fid.write(line0)
    line0 = '       0      0.00000       2       2       1 \n'
    fid.write(line0)
    fid.write('      '+str(nlong)+'     '+str(nlat)+'\n')
    # Only adding the original longitude shift as read from the FITS file and the CM of the CR
    fid.write(str(long0) +'  '+str(CRnumber)+' \n') #longitude shift (important for GONG Hourly)
    fid.write('Longitude Latitude Br LongitudeShift CarringtonRotation \n')
    
    for k in np.arange(nlat):
        for l in np.arange(nlong):
            # Phi grid is also written out as cell centered similar to the latitude grid
            line0 = str((l+0.5)*360./nlong) + ' ' + str(newlat[k]) + ' ' + str(newmap[k,l]) + ' \n'
            fid.write(line0)
    g.close()
    fid.close()
    nParam = 2
    Param_I = np.zeros(nParam)
    Param_I[0] = long0
#    Param_I[1] = -1 # Long Earth
    Param_I[1] = CRnumber # Long Earth
    # Long and Lat in radians passed to GLSetup.py
    Long_I = (pi/180.) * np.linspace(0.5*360./nlong, 359.5*360./nlong, nlong) # in radians
    Lat_I = newlat * (pi/180.) # radians
#    Bmax = 1900
    return(nlong, nlat, nParam, Param_I, Long_I, Lat_I, newmap)

###############CONSERVATIVE (ON SIN(THETA) UNIFORM GRID########                 
def smooth(nLong, nLat, nSmooth, Br_C):
    nSmooth2 = nSmooth//2
    Coef    = 1./(nSmooth*nSmooth)
    BrOrig_G = np.zeros([nLat,nLong+2*nSmooth2])
    for iLat in np.arange(nLat):
        BrOrig_G[iLat,:] = np.hstack((
                Br_C[iLat,nLong-nSmooth2:nLong],
                Br_C[iLat,:],Br_C[iLat,0:nSmooth2]))
    Br_C=np.zeros([nLat,nLong])
    for iLat in np.arange(nLat):
        for iLong in np.arange(nLong):
            for iSubLat in np.arange(nSmooth):
                iLatExt  = iLat  + iSubLat  - nSmooth2
                iLatExt  = max([-iLatExt-1,min(
                            [iLatExt, 2*nLat-1-iLatExt])])
                Br_C[iLat,iLong] += np.sum(
                    BrOrig_G[iLatExt,iLong:iLong+nSmooth])
            Br_C[iLat,iLong]  *= Coef
    return(Br_C)
    
def FITS_RECOGNIZE(inputfile):
    """
    This function opens inputfile and tries to determine what type of magnetogram
    it is as well as the type of grid on which the datatype is represented.  The
    magnetogram types and grid types are: 
      Hathaway Synchronic, regular
      ADAPT Synchronic, regular 
      GONG Synoptic, sin(lat)
      GONG Hourly updated, sin(lat)
      MDI Synoptic, sin(lat)
    This function returns a tuple with this information.  The output tuple is:
      (magnetogram_type,grid_type,InstrumentType)
    """

    magnetogram_type = 'unknown'
    grid_type = 'unknown'
    map_data = 'unknown'
    g = fits.open(inputfile)
    g.info()
    header0 = g[0].header
    # Print out the headers
#    print("====================================================\n")
#    print("Primary Extension Header:\n")
#    print(repr(header0))
#    print("====================================================\n")

    try:
        telescope = g[0].header['TELESCOP'] #works for MDI, GONG, HMI
    except KeyError as er:
        telescope = 'unknown'
        
    try:
        inst = g[0].header['INSTRUME'] #works for MDI, HMI
    except KeyError as er:
        inst = 'unknown'

    try:
        ctyp = g[0].header['CTYPE2'] #works for MDI, GONG, HMI
    except KeyError as er:
        ctyp = 'unknown'

    try:
        cunit = g[0].header['CUNIT2'] #works for MDI, GONG, HMI
    except KeyError as er:
       cunit = 'unknown'

    try:
        model = g[0].header['MODEL'] #works for ADAPT
    except KeyError as err:
        model = 'unknown'

    try:
        sft = g[0].header['SFT_TYP'] #works for Hathaway
    except KeyError as er:
        sft = 'unknown'

    nlo = g[0].header['NAXIS1'] # number of longitude points
    nla = g[0].header['NAXIS2'] #           latitude
        
    g.close()

    if telescope.find('NSO-GONG') > -1 :
        magnetogram_type = 'NSO-GONG Synoptic'
        try:
            long0 = g[0].header['LONG0']
            if float(long0) > 0.:
                magnetogram_type = 'NSO-GONG Hourly'
        except KeyError as er:
            long0 = - 1
        if ctyp.find('CRLT-CEA') > -1:
            grid_type = 'sin(lat)'            
            map_data = 'GONG'
        else:
            print ("unknown NSO-GONG magnetogram type")
            return(-1)

    if telescope.find('SOHO') > -1:
        if ( (inst.find('MDI') > -1) & (ctyp.find('Sine Latitude') > -1) ):
            magnetogram_type = 'MDI Synoptic'
            grid_type = 'sin(lat)'
            map_data = 'MDI'
        else:
            print ("unknown SOHO magnetogram type")
            return(-1)

    if telescope.find('SDO/HMI') > -1:
        if ( (ctyp.find('CRLT-CEA') > -1) & (cunit.find('Sine Latitude') > -1) ):
            magnetogram_type = 'HMI Synoptic'
            grid_type = 'sin(lat)'
            map_data = 'HMI'
        else:
            print ("unknown SDO magnetogram type")
            return(-1)

    if model.find('ADAPT') > -1:
        magnetogram_type ='ADAPT Synchronic'
        try:
            adapt_grid = g[0].header['GRID']
        except KeyError as er:
            adapt_grid = -1.
        if adapt_grid == 1.:
            grid_type = 'uniform'
        else:
            print ("unknown ADAPT magnetogram type")
            return(-1)
        try:
            map_data = str(g[0].header['MAPDATA']) #works for ADAPT
        except KeyError as er:
            map_data = ' '


    if sft.find('Baseline / Assimilation') > -1:
        magnetogram_type = 'Hathaway Synchronic'
        grid_type = 'uniform'

    if  ( (magnetogram_type == 'unknown') or (grid_type == 'unknown') ):
        print ("I don't recognize the type of this magnetogram.")
        return(-1)
                
    print ("I think this is a",magnetogram_type,"magnetogram on a",str(nla),"X",str(nlo),grid_type,"grid.")
    return( (magnetogram_type, grid_type, map_data) )

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="""
    remap_magnetogram.py pre-processes the FITS format magnetograms 
    into ASCII files that can by read by FDIPS.exe, BATSRUS.exe and SWMF.exe
    and IDL macros. The script can read the following types of magnetograms:

       Hathaway Synchronic
       ADAPT Synchronic
       GONG Synoptic
       GONG Hourly updated
       MDI Synoptic

    The code opens the .fits file and automatically recognizes the type of
    map it is, which determines whether it is on a sin(latitude) or regular
    spherical grid.  The output can be any desired resolution, on either a
    sin(latitude) or regular spherical grid.  If the output grid type is not
    specified, it will be the same as the original .fits file.  If the
    resolution is not specified, it will be the same as the original .fits
    file. The calling syntax from the command line is shown above. Some examples:

    ./remap_magnetogram.py test.fits test.out
    ./remap_magnetogram.py test.fits test.out 180 360
    ./remap_magnetogram.py test.fits test.out -grid=uniform
    ./remap_magnetogram.py test.fits -istart 1 -iend 12

    Within Python, the remapping is done with the remap function contained
    in this file.
    
    The script uses the scipy and astropy packages that can be installed, 
    for example, with MacPorts.
    """)
    parser.add_argument('inputfile', help='Input FITS file name including path')
    parser.add_argument('outputfile', nargs='?', default="map", help='Output file name including path but without the map index and .out extension. Default name is "map".')
    parser.add_argument('nlat', nargs='?', type=int, default=-1, help='Number of latitude points in output. Default is same as input.')
    parser.add_argument('nlon', nargs='?', type=int, default=-1, help='Number of longitude points in output. Default is same as input.')
    parser.add_argument('-grid',choices=['uniform','sinlat'],help="type of latitude grid in the output. Default is same as input.")
    parser.add_argument('-istart', type=int, help='Initial map index. Default is the first map.', default=1)
    parser.add_argument('-iend', type=int, help='Final map number. Default is the first map.', default=1)

    args = parser.parse_args()

    if args.nlat < -1:
        print ("nlat must be -1 or a postive integer.  No output.")
        quit()
    if args.nlon < -1:
        print ("nlon must be -1 or a postive integer.  No output.")
        quit()

    grid_type = 'unspecified'
    if args.grid == 'sinlat':
        grid_type = 'sin(lat)'
    elif args.grid == 'uniform':
        grid_type = 'uniform'

    if args.istart == None:
            args.istart = 0
    if args.iend < args.istart:
        args.iend = args.istart

    for i in range(args.istart,args.iend+1):
        out=args.outputfile+'_'+str(i)+'.out' 
        remap(args.inputfile, out, args.nlat, args.nlon, grid_type, i-1)
