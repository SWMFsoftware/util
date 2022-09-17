#!/usr/bin/env python
import numpy as np
import remap_magnetogram as rmag
import GLSETUPAlg as GL
import os
import fnmatch
import math
import subprocess
BMax = 1900.0
cPi  = np.pi
Rad2Deg = 180/cPi
Deg2Rad = cPi/180

def get_weighted_center(X,Y,Br_C,BrThreshold,nLat,nLon,Lat_I,Lon_I,\
                           IsUniformLat):
   LonIndex = GL.round_my(X)
   LatIndex = GL.round_my(Y)
   print('\n Chosen Longitude, Latitude =',Lon_I[LonIndex]*Rad2Deg,
         Lat_I[LatIndex]*Rad2Deg)
   #Occcupancy matrix
   occ = np.zeros([nLat,nLon])
   occ[LatIndex,LonIndex] = 1
   #Occupancy level
   occ_level = 0
   #Occupancy level check
   occ_check = 1
   while occ_check > 0:
      occ_level = occ_level + 1
      [row,col] = np.where(occ == occ_level)
      n = row.size
      # row= lat, col = long
      for i in np.arange(n):
         rowp = row[i]
         colp = col[i]
         if BrThreshold > 0. :
            if rowp-1 > 0:
               if (Br_C[rowp-1,colp] > BrThreshold and occ[rowp-1,colp] == 0):
                  occ[rowp-1,colp] = occ_level + 1
            if rowp+1 < nLat:
               if (Br_C[rowp+1,colp] > BrThreshold and occ[rowp+1,colp] == 0):
                  occ[rowp+1,colp] = occ_level + 1
            if colp-1 > 0:
               if (Br_C[rowp,colp-1] > BrThreshold and occ[rowp,colp-1] == 0):
                  occ[rowp,colp-1] = occ_level + 1
            if colp+1 < nLon:
               if (Br_C[rowp,colp+1] > BrThreshold and occ[rowp,colp+1] == 0):
                  occ[rowp,colp+1] = occ_level + 1
         elif BrThreshold < 0.:
            if rowp-1 > 0:
               if (Br_C[rowp-1,colp] < BrThreshold and occ[rowp-1,colp] == 0):
                  occ[rowp-1,colp] = occ_level + 1
            if rowp+1 < nLat:
               if (Br_C[rowp+1,colp] < BrThreshold and occ[rowp+1,colp] == 0):
                  occ[rowp+1,colp] = occ_level + 1
            if colp-1 > 0:
               if (Br_C[rowp,colp-1] < BrThreshold and occ[rowp,colp-1] == 0):
                  occ[rowp,colp-1] = occ_level + 1
            if colp+1 < nLon:
               if (Br_C[rowp,colp+1] < BrThreshold and occ[rowp,colp+1] == 0):
                  occ[rowp,colp+1] = occ_level + 1
      occ_check = n
   #end whileloop
   SizeMap_C = np.zeros([nLat,nLon])
   #Calculate weighted center
   [LatOcc,LonOcc] = np.where(occ>0)
   nSize = LatOcc.size
   LatCenter=0.
   LonCenter=0.
   Flux=0.
   Area = 0.
   #flux = SUM(area * Br)
   dLon    = 2.*cPi/nLon
   dLat    = cPi/nLat
   dSinLat = 2.0/nLat
   for i in np.arange(nSize):
      iLon  = LonOcc[i]
      iLat  = LatOcc[i]
      if IsUniformLat : # uniform in lat
         dArea = np.cos(Lat_I[iLat]) * dLat * dLon   # in radians^2
      else:
         dArea = dSinLat * dLon   # in radians^2
      dFlux = Br_C[iLat,iLon] * dArea
      LonCenter += Lon_I[iLon] * dFlux
      LatCenter +=  Lat_I[iLat] * dFlux
      Flux += dFlux
      Area += dArea
      SizeMap_C[iLat,iLon]=Br_C[iLat,iLon]
   LonCenter /= Flux  # in Radians
   LatCenter /= Flux  # in Radians
   if BrThreshold > 0. :
      print("\n Flux from the positive spot={:6.2f} [Gs*Rs**2]".format(
         Flux))
   else:
      print("\n Flux from the negative spot={:6.2f} [Gs*Rs**2]".format(
         Flux))
   # return the longitude and latitude of the weighted center in radians, 
   # area in radians^2, and the occupancy matrix for plotting
   return(LatCenter,LonCenter,SizeMap_C,Flux, LatOcc, LonOcc)

def TwoPointsOnSph(Lon1, Lat1, Lon2, Lat2):
    SinLon1 = np.sin(Lon1)
    CosLon1 = np.cos(Lon1)
    SinLat1 = np.sin(Lat1)
    CosLat1 = np.cos(Lat1)
    Dir1_D = np.array([CosLat1*CosLon1, CosLat1*SinLon1, SinLat1])

    SinLon2 = np.sin(Lon2)
    CosLon2 = np.cos(Lon2)
    SinLat2 = np.sin(Lat2)
    CosLat2 = np.cos(Lat2)
    Dir2_D = np.array([CosLat2*CosLon2, CosLat2*SinLon2, SinLat2])

    HalfDist2  = np.sin(
       0.50*(Lat2 - Lat1))**2 + CosLat1*CosLat2*np.sin(
          0.50*(Lon2 - Lon1))**2
    rMidPoint = np.sqrt(1 - HalfDist2)

    HalfDist = np.sqrt(HalfDist2)
    # Depth = 1 - rMidPoint
    Depth = HalfDist2/(1 + rMidPoint)

    # Unit vector toward the mid point:
    DirMid_D = Dir2_D + Dir1_D
    # Normalization coefficient
    Coeff1 = 0.50/rMidPoint
    DirMid_D = DirMid_D*Coeff1
    SinLat = DirMid_D[2]
    CosLat = np.sqrt(1 - SinLat**2)
    CosLon = DirMid_D[0]/CosLat
    SinLon = DirMid_D[1]/CosLat
    Lat = np.arcsin(SinLat)
    Lon = np.arccos(CosLon)
    if SinLon < 0.0 :
       Lon = 2*cPi - Lon
    # Direction  vectors for parallel and meridian:
    DirPar_D = np.array([-SinLon, CosLon, 0.0])
    DirMer_D = np.array([-SinLat*CosLon, -SinLat*SinLon, CosLat])
    # Direction unit vector from point 1 to point 2:
    Dir12_D = Dir2_D - Dir1_D
    # Normalization coefficient
    Coeff1  = 0.50/HalfDist
    # Normalized vector from 1 to 2
    Dir12_D = Dir12_D*Coeff1
    # Its projections on the parallel and  meridian  directions:
    CosOrientation  = Dir12_D[0]*DirPar_D[0]
    CosOrientation += Dir12_D[1]*DirPar_D[1]
    CosOrientation += Dir12_D[2]*DirPar_D[2]
    SinOrientation  = Dir12_D[0]*DirMer_D[0]
    SinOrientation += Dir12_D[1]*DirMer_D[1]
    SinOrientation += Dir12_D[2]*DirMer_D[2]
    Orientation  = np.arccos(CosOrientation)
    if SinOrientation < 0.0 :
       Orientation = 2*cPi - Orientation

    return(Lon,Lat,Orientation,HalfDist)
   
def Alg(nLon, nLat, nParam, Param_I, Lon_I, Lat_I, Br_C, UseCMEGrid, 
        Helicity, Time):
   Lon0     = Param_I[0]
   LonEarth = Param_I[1]
   xPositive = Param_I[2]
   yPositive = Param_I[3]
   xNegative = Param_I[4]
   yNegative = Param_I[5]
   dLon    = 2.*cPi/nLon
   dLat    = cPi/nLat
   dSinLat = 2.0/nLat
   DsLat_C = np.zeros([nLat,nLon])
   Ds2_C   = np.zeros([nLat,nLon])
   # Check if the latitude grid is uniform
   if abs(Lat_I[2]-2*Lat_I[1]+Lat_I[0])<1.0e-5:
      IsUniformLat = True
      print('Uniform in Latitude grid')
      for k in np.arange(nLat):
         for l in np.arange(nLon):
            DsLat_C[k,l] = dLat
            Ds2_C[k,l]   = dLat**2 + (np.cos(Lat_I[k])*dLon)**2
   else:
      IsUniformLat = False
      print('Uniform in Sin(Latitude) grid')
      for k in np.arange(nLat):
         for l in np.arange(nLon):
            DsLat_C[k,l] = dSinLat/np.cos(Lat_I[k])
            Ds2_C[k,l]   = DsLat_I[k,l]**2 + (np.cos(Lat_I[k])*dLon)**2
   # Pass the x, y indices of the clicks to calculate weighted center
   # and their indices

   # get weighted centers(Lon,Lat), occupancy matrix, Area of AR for
   # positive and negative regions
   [LatPos,LonPos,PSizeMap_C,FluxP, LatP_I, LonP_I] = \
       get_weighted_center(xPositive,yPositive,Br_C,20.,\
                              nLat,nLon,Lat_I,Lon_I,IsUniformLat)
   LonPMin =  min(LonP_I)
   LonPMax =  max(LonP_I)
   print('\n Positive spot: minimum and maximum longitude  indexes',\
         LonPMin,LonPMax)
   LatPMin =  min(LatP_I)
   LatPMax =  max(LatP_I)
   print('\n Positive spot: minimum and maximum latitude  indexes',\
         LatPMin,LatPMax)
   
   LonPosIndex = GL.calculate_index(LonPos,Lon_I,nLon)
   LatPosIndex = GL.calculate_index(LatPos,Lat_I, nLat)
   print('\n Positive Weighted Center indexes (lon,lat) =',\
          LonPosIndex, LatPosIndex)
   print(
      '\n Positive Weighted Center Lon={0:6.2f}, Lat={1:6.2f} [deg]'.format(
         LonPos*Rad2Deg, LatPos*Rad2Deg))

   [LatNeg,LonNeg, NSizeMap_C,FluxN, LatN_I, LonN_I] = \
       get_weighted_center(xNegative,yNegative,Br_C,-20.,\
                              nLat,nLon,Lat_I,Lon_I,IsUniformLat)
   LonNMin =  min(LonN_I)
   LonNMax =  max(LonN_I)
   print('\n Negative spot: minimum and maximum longitude  indexes',\
         LonNMin,LonNMax)
   LatNMin =  min(LatN_I)
   LatNMax =  max(LatN_I)
   print('\n Negative spot: minimum and maximum latitude  indexes',\
         LatNMin,LatNMax)
   
   LonNegIndex = GL.calculate_index(LonNeg,Lon_I,nLon)
   LatNegIndex = GL.calculate_index(LatNeg,Lat_I, nLat)
   print(
      '\n Negative Weighted Center indexes (lon,lat) =',
      LonNegIndex, LatNegIndex)
   print(
      '\n Negative Weighted Center Lon={0:6.2f}, Lat={1:6.2f} [deg]'.format(
         LonNeg*Rad2Deg,LatNeg*Rad2Deg))

   [LonAR,LatAR,OrientationAR,HalfDist] = TwoPointsOnSph(
      LonPos,LatPos,LonNeg,LatNeg)
   # Major radius estimated from the distance between the spot centers
   RadiusAR = HalfDist*1.4
   print ("Mid point: Longitude: {0:5.2f} Latitude:{1:5.2f}".format(
         LonAR*Rad2Deg, LatAR*Rad2Deg))
   print ("Mid point: Orientation: {0:5.2f} Major Radius: {1:5.2f}".format(
         OrientationAR*Rad2Deg, RadiusAR))

   # Find center of the active region as the point on the line
   # connecting the positive and negative center at which the MF is minimal
   # (i.e as intersection of this with PIL,
   # herewith PIL=Polarity Inversion Line
   nProfile = max([abs(LonPosIndex - LonNegIndex),
                   abs(LatPosIndex - LatNegIndex)]) + 1
   BTmp = BMax + 1.0
   for i in np.arange(nProfile):
      LonProfile = LonPos+(
         LonNeg - LonPos)*i/(nProfile - 1)
      LatProfile = LatPos+(
         LatNeg - LatPos)*i/(nProfile - 1)
      IndexLon = GL.calculate_index(LonProfile,Lon_I,nLon)
      IndexLat = GL.calculate_index(LatProfile,Lat_I,nLat)
      AbsBr = abs(Br_C[IndexLat,IndexLon])
      if (AbsBr < BTmp):
         BTmp = AbsBr
         iLonAR =IndexLon
         iLatAR =IndexLat
   LonAR  = Lon_I[iLonAR]  # in radians
   LatAR  = Lat_I[iLatAR]  # in radians
   print (
      "Center for Active region: Lon={0:4.1f}, Lat={1:4.1f} [deg]:".format(
         LonAR*Rad2Deg,LatAR*Rad2Deg))
   # Rectangular box  for active region
   LonARMin=min([LonNMin,LonPMin])
   LonARMin=max([LonARMin-2,0])
   LonARMax=max([LonNMax,LonPMax])
   LonARMax=min([LonARMax+2,nLon-1])
   LatARMin =min([LatNMin, LatPMin])
   LatARMin =max([LatARMin-2,0])
   LatARMax =max([LatNMax, LatPMax])
   LatARMax =min([LatARMax+2,nLat-1])
   print('\n Box for AR: minimum and maximum longitude  indexes',\
         LonARMin,LonARMax)
   print('\n Box for AR: minimum and maximum longitude  indexes',\
         LatARMin,LatARMax)

   # Distance to negative spot
   nSizeN = LatN_I.size
   Dist2N_I=np.zeros(nSizeN)
   nSizeP=LatP_I.size
   Dist2P_I=np.zeros(nSizeP)
   nLonShort = LonARMax + 1 - LonARMin
   nLatShort = LatARMax + 1 - LatARMin
   Dist2Min_C=np.zeros([nLatShort,nLonShort])
   for k in range(LatARMin , LatARMax+1):
      CosLat=np.cos(Lat_I[k])
      for l in range(LonARMin, LonARMax+1):
         for  i in np.arange(nSizeN):
            Dist2N_I[i] = (Lat_I[LatN_I[i]]-Lat_I[k])**2+CosLat**2*(
               Lon_I[LonN_I[i]]-Lon_I[l])**2
         for  i in np.arange(nSizeP):
            Dist2P_I[i] = (Lat_I[LatP_I[i]]-Lat_I[k])**2+CosLat**2*(
               Lon_I[LonP_I[i]]-Lon_I[l])**2
         Dist2Min =max([min(Dist2N_I), min(Dist2P_I)])
         if Dist2Min<1.5*Ds2_C[k,l]:
            Dist2Min_C[k-LatARMin,l-LonARMin]=Ds2_C[k,l]/Dist2Min

   nParam = 8
   Param_I = np.array(
      [Lon0+LonARMin, LonEarth, LonPosIndex-LonARMin,LatPosIndex-LatARMin,
       LonNegIndex-LonARMin,LatNegIndex-LatARMin,iLonAR-LonARMin,iLatAR-LatARMin])
   nVar=4
   Data_IV=np.zeros([nLatShort,nLonShort,nVar])
   NameVar='Longitude Latitude Br MapP MapN  PIL PIL Lon0 LonEarth xP yP xN yN xC yC'
   for k in np.arange(nLatShort):
      for l in np.arange(nLonShort):
         Data_IV[k,l,0]=max([-BMax,min([BMax,Br_C[k+LatARMin,l+LonARMin]])])
         Data_IV[k,l,1]=PSizeMap_C[k+LatARMin,l+LonARMin]
         Data_IV[k,l,2]=NSizeMap_C[k+LatARMin,l+LonARMin]
         Data_IV[k,l,3]=Dist2Min_C[k,l]
   FinalFile=rmag.save_bats('AfterGLSETUP.out', 'After GLSETUP: Br[Gauss]', 
                            NameVar, [nLonShort,nLatShort], nVar, nParam,
                            Param_I, Rad2Deg*Lon_I[LonARMin:LonARMax+1],
                            Rad2Deg*Lat_I[LatARMin:LatARMax+1], Data_IV, Time)
   print('Select the CME Source Region (POSITIVE) with the left button')
   print('Then select negative region with the right button')

   FileId=open('runidl2','w')
   FileId.write(';\n;\n')
   FileId.write(
      "      TDSETUP1,file='"+FinalFile+"' \n")
   FileId.close()
   ########SHOW MAGNETOGRAM##########################
   # GLSETUP1.pro is run, it reads the magnetogram(fitsfile.out)
   # reads the cursor x,y indices for neg and pos. AR.
   ls = subprocess.Popen(["idl", "runidl2"],stdout=subprocess.PIPE,
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
   exit()
   GL_Latitude  = LatAR * Rad2Deg
   GL_Longitude = LonAR * Rad2Deg
   if Lon0>0:
      GL_Longitude +=Lon0
      if GL_Longitude>=360:
         GL_Longitude-=360
   print ("GL_Longitude: {0:4.1f} GL_Latitude:{1:4.1f}".format(
         GL_Longitude, GL_Latitude))
   PointN_I=[LonNeg,LatNeg] # negative spot
   PointP_I=[LonPos,LatPos] # positive spot

   AngularDistance = GL.get_angular_dist(PointN_I,PointP_I)


   # GL_Orientation calculation
   # Calculate the GL flux rope orientation from the two weighted points.
   #r1=[LonNegIndex-LonPosIndex,LatNegIndex-LatPosIndex] - incorrect
   r1 = [PointN_I[0] - PointP_I[0], PointN_I[1] - PointP_I[1]]
   r1[0] *= np.cos(LonAR)
   r1 /= np.sqrt(r1[0]**2+r1[1]**2)
   r2=[1.0,0.0]
   GL_Orientation=np.arccos(r1[0]*r2[0]+r1[1]*r2[1])*Rad2Deg
   if r1[1] < 0:
      # If sine of Orientation is negative
      GL_Orientation=360-GL_Orientation  
   if IsUniformLat :
      grid_type = 'uniform'
   else:
      grid_type = 'sin(lat)'

   if Helicity != 0 : # helicity input by user
      iHelicity = Helicity
      print('Using user input helicity = ', Helicity)
   else:
      # based on hemisphere
      iHelicity = 1
      if GL_Latitude > 0: 
         iHelicity = -1
      print('Helicity based on hemisphere: ',iHelicity)

      iHelicity = -1
      Depth=0
   #Recommended GL flux rope parameters
   ### TEMPORARY !!!!!
   AngularDistance = 0.2
   print ('========================================')
   print ('The Recommended GL FLux Rope Parameters')
   print ('========================================')
   print ('#CME')
   print ('                Latitude: %6.2f'%(GL_Latitude))
   print ('               Longitude: %6.2f'%(GL_Longitude))
   print ('             Orientation: %6.2f'%(GL_Orientation))
   print ('      Angular size [deg]: %6.2f'%(AngularDistance*Rad2Deg))
   print (' Poloidal flux: positive: %6.2f'%(FluxP))
   print (' Poloidal flux: negative: %6.2f'%(FluxN))
   print ('-----------------------------------------')
   FileId=open('CME.in','a')
   FileId.write("#CME \n")
   FileId.write("T                   UseCme \n")
   FileId.write("T                   DoAddFluxRope \n")
   FileId.write("%-10.2f          LongitudeCme \n"% GL_Longitude)
   FileId.write("%-10.2f          LatitudeCme \n"% GL_Latitude)
   FileId.write("%-10.2f          OrientationCme \n"% GL_Orientation)
   FileId.write("TD22                  TypeCme \n")
   FileId.write("%-+d                  iHelicity \n"% iHelicity)
   FileId.write("%-10.2f	  RadiusMajor \n"%(AngularDistance))
   FileId.write("%-10.2f	  RadiusMinor \n"%(0.35*AngularDistance))
   FileId.write("%-10.2f	  Depth \n"%Depth)
   FileId.write("1.0e-3                 PlasmaBeta \n")
   FileId.write("5.0e5                 EjectaTemperature \n")
   FileId.write("readbstrap            TypeBStrap \n")
   FileId.write("5.0                 BStrappingDim \n")
   FileId.write("none            TypeCharge \n")
   FileId.write(" \n")
   FileId.write("#END \n")
   FileId.write("\n")
   FileId.write("Angular Size            = %5.2f\n"%(AngularDistance
                                                     *Rad2Deg))
   print (' Poloidal flux: positive: %6.2f'%(FluxP))
   print (' Poloidal flux: negative: %6.2f'%(FluxN))
   FileId.close() 

   if UseCMEGrid:
      #Calculate the CME grid refinement parameters based on the flux rope
      #location and size.                                                
      
      print ('==========================================')
      print ('The Recommended Grid Refinement Parameters')
      print ('==========================================')
      print ('              R_Start: %6.2f'% (CMEbox_Start[0]))
      print ('                R_End: %6.2f'% (CMEbox_End[0]))
      print ('      Longitude_Start: %6.2f'% ( CMEbox_Start[1]))
      print ('        Longitude_End: %6.2f'% ( CMEbox_End[1]))
      print ('       Latitude_Start: %6.2f'% ( CMEbox_Start[2]))
      print ('         Latitude_End: %6.2f'% ( CMEbox_End[2]))
      print ('-----------------------------------------')
      FileId=open('CME_AMR.in','w')
      FileId.write("#AMRREGION \n")
      FileId.write("CMEbox              NameRegion \n")
      FileId.write(" \n")
      FileId.write("#END \n")
      FileId.close()
   #For comparison, make magnetogram of a flux rope field
   FileId=open('RunFRM','w')
   FileId.write('%-3d \n'%Lon0)
   if IsUniformLat :
      FileId.write('uniform latitude \n')
   else:
       FileId.write('sin(latitude) \n')
   FileId.close()
   FileId=open('RunFRM','r')
   subprocess.call('./FRMAGNETOGRAM.exe',stdin=FileId)
   FileId.close()
   
   nParam = 8
   Param_I = np.zeros(nParam)
   Param_I[0:8] = [Lon0,LonEarth,LonPosIndex,LatPosIndex,LonNegIndex,
                   LatNegIndex,iLonAR,iLatAR]
   FileId = open('AfterGLSETUP.out','w')
    
   FileId.write('After GLSETUP: Br[Gauss]'+'\n')
   FileId.write(
      '       0     '+str(Time)+'     2      %2d       3 \n'% nParam)
   FileId.write('      '+str(nLon)+'     '+str(nLat)+'\n')
   FileId.write(
      ' {0:5.1f} {1:5.1f} {2:5.1f} {3:5.1f} {4:5.1f} {5:5.1f} {6:5.1f} {7:5.1f}'.format(
         Lon0,LonEarth,LonPosIndex,LatPosIndex,LonNegIndex,LatNegIndex,iLonAR,iLatAR))

   FileId.write('\n')
   FileId.write(
      'Longitude Latitude Br PMap NMap Lon0 LonEarth xP yP xN yN xC yC \n')
   
   for k in np.arange(nLat):
      for l in np.arange(nLon):
         FileId.write("{0:6.1f} {1:6.1f} {2:14.6e} {3:14.6e} {4:14.6e} \n".format((180./cPi)*Lon_I[l],(180./cPi)*Lat_I[k],max([-BMax,min([BMax,Br_C[k,l]])]),PSizeMap_C[k,l],NSizeMap_C[k,l]))
    
   FileId.close()

   return(nLon,nLat,nParam,Param_I,Lon_I,Lat_I,Br_C,PSizeMap_C,NSizeMap_C)

