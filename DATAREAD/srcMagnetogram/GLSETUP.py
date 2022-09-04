#!/usr/bin/env python
####################### MASTER SCRIPT FOR EEGGL or SWMF_GLSETUP#######
####################### CAN BE IMPLEMENTED IN ANY SCRIPT LANGUAGE#####
####April 2020: Updated script for different types of magnetograms
# make FRM in SWMF/util/EMPIRICAL/srcEE/

import subprocess
import os
import fnmatch
import remap_magnetogram as rmag
import numpy as np
import argparse
import GLSETUPAlg as GL
from swmfpy.web import download_magnetogram_hmi as hmi_map
import datetime as dt

BMax = 1900.0
cPi = np.pi
Rad2Deg = 180/cPi
Deg2Rad = 1/Rad2Deg
IsPositionInput = 0

if __name__ == '__main__':
   parser = argparse.ArgumentParser(
      formatter_class=argparse.RawTextHelpFormatter)
   parser.add_argument('NameFile', help='Input FITS file name including path')
   parser.add_argument('nlat', nargs='?', type=int, default=180, help=
                       'Number of latitude points in output. Default is 180.')
   parser.add_argument('nlon', nargs='?', type=int, default=360, help= 
                       'Number of longitude points in output. Default is 360.')
   parser.add_argument('-outgrid',choices=['uniform','sinlat'],help=
                       'type of latitude grid in the output. Default is same as input.')
   parser.add_argument('-index', type=int, default=1, help=
                       'Initial map index. Default is the first map.')
   parser.add_argument('-nSmooth',type=int, default=5, help=
                       'If nSmooth is ODD integer larger than 1, apply boxcar smoothing on the magnetic field. This can help finding the PIL')
   parser.add_argument('-CMESpeed',type=float, default=-1.0, help=
                       'CME speed in km/s, recovered from observations')
   parser.add_argument('--CMEGrid',action='store_true', help=
                       'Output parameters of the refined CME grid')
   parser.add_argument('-Orientation',type=float,default=-1.0, help=
                       'Use specified orientation for GL_Orientation')
   parser.add_argument('-Helicity',type=int,default=0, help=
                       'Use specified helicity')
   parser.add_argument('-LonPosIn',type=float,default=999.0, help=
                       'Longitude for positive spot center (Deg)')
   parser.add_argument('-LatPosIn',type=float,default=999.0, help=
                       'Latitude for positive spot center (Deg)')
   parser.add_argument('-LonNegIn',type=float,default=999.0, help=
                       'Longitude for negative spot center (Deg)')
   parser.add_argument('-LatNegIn',type=float,default=999.0, help=
                       'Latitude for negative spot center (Deg)')
   parser.add_argument('-GLRadius',type=float, default=-1., help=
                       'Radius of the flux-rope before stretching.')
   parser.add_argument('-Stretch',type=float, default=-1., help=
                       'Stretching parameter of the flux-rope.')
   parser.add_argument('-Distance',type=float, default=-1., help=
                       'Distance parameter of the flux-rope.')
   parser.add_argument('--UsePNDist',action='store_true', help=
                       'Use the angular distance between positive and negative spot centers to calculate AR size and GL Radius')
   parser.add_argument('--UseARArea',action='store_true',help=
                       'Use Active Region area to calculate AR size and GL Radius')
   parser.add_argument('--DoScaling',action='store_true',help=
                       'Scale Distance, Stretch, Radius by a factor Alpha(a). Distance->a*Distance, Stretch->a*Stretch, GLRadius->a*GLRadius')
   parser.add_argument('-SizeFactor',type=float, default=1.0, help=
                       'Calculated from the angular distance between the  positive and negative spot centers. If specified, then the Angular width of the flux rope on the solar surface is SizeFactor * Distance between the spot centers. Default SizeFactor is 1.0.')
   parser.add_argument('--GLRadiusRange',nargs='+',type=float, 
                       default=[0.2,0.8], help=
                       'GLRadius is limited by GLRadiusRange = 2-elements array. Default is [0.2,2.0].')
   parser.add_argument('--DoHMI',action='store_true', help=
                       'Use HMI map for helicity determination')
   parser.add_argument('--UseBATS',action='store_true', help=
                       'Reading magnetogram in the ModPlotFile format')
   args = parser.parse_args()
   ##################OPTIONAL INPUT PARAMETERS######
   NameFile    = args.NameFile
   CMESpeed    = args.CMESpeed
   GLRadius    = args.GLRadius
   SizeFactor  = args.SizeFactor
   GLRadiusRange_I = args.GLRadiusRange
   UseCMEGrid  = args.CMEGrid
   nlat        = args.nlat
   nlon        = args.nlon
   i           = args.index
   nSmooth     = args.nSmooth
   Orientation = args.Orientation
   Helicity    = args.Helicity
   Stretch     = args.Stretch
   Distance    = args.Distance
   DoHMI       = args.DoHMI # default is False
   xPositive   = args.LonPosIn
   yPositive   = args.LatPosIn
   xNegative   = args.LonNegIn
   yNegative   = args.LatNegIn
   UsePNDist   = args.UsePNDist
   UseARArea   = args.UseARArea
   DoScaling   = args.DoScaling
   UseBATS     = args.UseBATS

   IdlFile = 'fitsfile.out'
   if UseBATS==False:
      # Check if the file extension is .out
      SplitName = NameFile.split('.')
      if  SplitName[-1]=='out':
         print('\n File name '+NameFile+
               ' has extension .out, is treated as ASCII converted file')
         UseBATS = True

   if not UsePNDist and not UseARArea and GLRadius <=0. :
      print('\n WARNING: User did not specify how to calculate GLRadius.'+
            ' Default is using Active Region area to estimate GLRadius.')
      UseARArea = True
  
   if (xPositive !=999. and yPositive !=999. and yNegative !=999.
       and xNegative !=999.):
      IsPositionInput = 1
      print('User input the x,y positions for Positive and Negative centers')
      print('Input Weighted centers :',xPositive,yPositive,xNegative,yNegative)
   else:
      IsPositionInput = 0
   
   # setting ouput grid if remapping is required
   grid_type = 'unspecified' # assumes the grid of the map
   if args.outgrid == 'sinlat':
      grid_type = 'sin(lat)'
   elif args.outgrid == 'uniform':
      grid_type = 'uniform'

   ##################END OF PARSER#####################
   #################SERVER SIDE, PYTHON################
   #################PROCESS MAGNETOGRAM###
   ##READ AND SMOOTH, IF DESIRED########################
   if UseBATS:
      cc =  rmag.read_bats(NameFile)
      nIndex_I     = cc[0]
      nLong        = nIndex_I[0]
      nLat         = nIndex_I[1]
      nVar         = cc[1]
      nParam       = cc[2]
      Param_I      = cc[3]
      Long0        = Param_I[0] # Longitude of left edge
      Time         = cc[7]
      LongEarth    = Param_I[1] + Time # CR number plus CRFraction
      Long_I       = cc[4]*Deg2Rad      # in radians
      Lat_I        = cc[5]*Deg2Rad      # in radians
      data         = cc[6]
      if nVar ==1:
         Br_C = data
      else:
         Br_C = data[:,:,0]
      if nSmooth > 2:
         Br_C = rmag.smooth(nLong,  nLat,  nSmooth, Br_C)
      IdlFile = NameFile
   else:
      # fits magnetogram is read, remapped (if required) using
      # remap_magnetogram.py to fitsfile.out
      cc = rmag.remap(NameFile, IdlFile, nlat, nlon, grid_type,
                      i-1, nSmooth,BMax)
      nLong        = cc[0]
      nLat         = cc[1]
      nParam       = cc[2]
      Param_I      = cc[3]
      Long0        = Param_I[0] # Longitude of left edge
      LongEarth    = Param_I[1] # CR number of central meridian
      Long_I       = cc[4]      # in radians
      Lat_I        = cc[5]      # in radians
      Br_C         = cc[6]
      Time         = LongEarth - int(LongEarth)
      if DoHMI:
         date         = cc[8]
         hmi_yymm = date.split("-")
         hmi_dd = hmi_yymm[2].split("T")
         hmi_hh = hmi_dd[1].split(":")
         hmi_yyyy     = int(hmi_yymm[0])
         hmi_mm       = int(hmi_yymm[1])
         hmi_dd       = int(hmi_dd[0])
         hmi_hh       = int(hmi_hh[0])
         cwd = os.getcwd()
         time_mag = dt.datetime(hmi_yyyy, hmi_mm, hmi_dd, hmi_hh)
         hmi_file = hmi_map(mag_time=time_mag, hmi_map='hmi.b_synoptic_small',
                         download_dir=cwd)

   #Info to the idl session is passed via the fitsfile.out file####
   ############END OF PYTHON FIRST SESSION##########
   ###IDL SESSION IN THE SWMF_GLSETUP/BROWSER SESSION IN EEGGL##

   if CMESpeed<= 0.0: 
      CMESpeed = float(raw_input(
            'Please Input the Observed CME Speed (km/s): '))

   if IsPositionInput == 0:
      print('Select the CME Source Region (POSITIVE) with the left button')
      print('Then select negative region with the right button')

      FileId=open('runidl1','w')
      FileId.write(';\n;\n')
      FileId.write(
         "      GLSETUP1,file='"+IdlFile+"',/UseBATS ")
      FileId.close()
   ########SHOW MAGNETOGRAM##########################
   # GLSETUP1.pro is run, it reads the magnetogram(fitsfile.out)
   # reads the cursor x,y indices for neg and pos. AR.
      ls = subprocess.Popen(["idl", "runidl1"],stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,text=True)
   #################PROCESSING STDOUT################
      stdout,stderr=ls.communicate()
      b=stdout[stdout.index('===')+4:len(stdout)]
      a=b.split() # x,y coordinates 
   ###### TAKE TWO COORDINATES FROM TWO CLICKS#######
      xPositive = float(a[0])
      yPositive = float(a[1])
      xNegative = float(a[2])
      yNegative = float(a[3])
   ##########SHAPE INPUTS FOR THE SECOND SERVER-SIDE SESSION####
   nParam  = 6
   Param_I = np.zeros(nParam)
   Param_I[0] = Long0
   Param_I[1] = LongEarth
   # Below the x,y positions are equal to location of clicks OR 
   # the location of weighted centers as input by the user, IsPositionInput =1
   # These are passed to GLSETUPALg.py and weighted centers are calculated
   Param_I[2] = xPositive
   Param_I[3] = yPositive
   Param_I[4] = xNegative
   Param_I[5] = yNegative

   ##SECOND SERVER-SIDE SESSION (PYTHON)#######################
   CC=GL.Alg(nLong,nLat,nParam,Param_I,Long_I,Lat_I,Br_C,
             CMESpeed,GLRadius,SizeFactor,
             GLRadiusRange_I, UseCMEGrid, Orientation,
             Stretch, Distance, Helicity, DoHMI, IsPositionInput,
             UsePNDist, UseARArea, DoScaling, Time)
   ##SHAPE INPUT PARAMETERS FOR THE CONCLUDING SESSION#########
   # nLong   =    CC[0] 
   # nLat    =  CC[1] 
   # nParam  =  CC[2]
   # Param_I =  CC[3]
   # Long0     =  Param_I[0]
   # LongEarth =  Param_I[1]
   # xPositive =  Param_I[2]
   # yPositive =  Param_I[3]
   # xNegative =  Param_I[4]
   # yNegative =  Param_I[5]
   # XyARCenter_D = Param_I[6:8]
   # nPIL        = (nParam - 8)//2
   # xPIL_I      =  Param_I[8:8+nPIL]
   # yPIL_I      =  Param_I[8+nPIL:nParam]
   # Long_I      = CC[4]
   # Lat_I       = CC[5]
   # Br_C        =  CC[6]
   # PSizeMap_C  =  CC[7]
   # NSizeMap_C  =  CC[8]
   # occPos      = CC[9]
   # occNeg      = CC[10]

   FileId=open('runidl','w')
   FileId.write(';\n;\n')
   FileId.write("GLSETUP2, file='AfterGLSETUP.out',/UseBATS \n")
   FileId.close()
   ###FINAL SESSION: SHOW MAGNETOGRAM AND BIPOLAR STRUCTURE OF AR
   subprocess.call(['idl','runidl'])
   ###IF THE MASTER SCRIPT IS IN PYTHON, AND A CHILD PROCESS IS IN IDL
   #(1) THE TIME OF IDL SESSION SHOULD BE LIMITED (30 seconds or so) 
   #(2) WINDOWS SHOULD BE CLOSED 
   #(3) FINAL EXIT COMMAND MUST BE PRESENT IN THE IDL SCRIPT######
   #print 'GLSETUP Session is closed. Bye!!!'
##############################################
