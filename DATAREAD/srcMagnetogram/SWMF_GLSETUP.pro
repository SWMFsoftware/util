pro SWMF_GLSETUP, DemoMode=DemoMode, PlotRadius=PlotRadius, $
                  USEPIL=USEPIL, CMEGrid=CMEGrid, ARMag=ARMag, $
                  GLRadius=GLRadius, SizeFactor=SizeFactor, $
                  nSmooth=nSmooth, FILE=FILE, CMEspeed=CMEspeed,$
                  Help=Help

;-----------------------------------------------------------------------
; NAME:
;   SWMF_GLSETUP
; PURPOSE:
;   Determine the Gibson-Low flux rope parameters from the input
;   magnetogram and observed CME speed.
;
; INPUT PARAMETERS:
;   The observed CME speed, the input magnetic field of SWMF, the
;   location of the CME source region(Interactive Selection).
;
; OUTPUTS:
;   Recommended GL flux rope parameters.
;
; KEYWORDS:
;   DemoMode = If set, a pre-saved magnetogram data will be loaded.
;
;   PlotRadius = Set up the layer of the magnetogram for 3D input
;   Cannot be used 2D file or with DemoMode. Default is 1.03.
;
;   UsePIL = If set, the orientation of the flux rope will be
;   calculated according to the PIL direction.
;
;   CMEGrid = If set, the grid refinement parameters for CME will be
;   calculated.
;
;   ARMag = 1 or 2, if 1, the AR magnetic parameter is calculated
;   based on the average Br around the weighted center; if 2, the
;   AR magnetic parameter is calculated based on the average total field
;   along the PIL. The corresponding empirical relationships will be
;   used accordingly. The default is 1.
;
;   GLRadius = Sets the GL flux rope radius (before shift). No default.
;
;   SizeFactor = If GLRadius is not set, the flux rope size is set to
;                the PIL length divided by SizeFactor. Default is 26.25.
;
;   nSmooth = If nSmooth is larger than 1, apply boxcar smoothing on the 
;             magnetic field. This can help finding the PIL.
;   
;   Help = If set, print all available keywords
;
; RESTRICTIONS:
;   - PlotRadius can only be 0.015 increament from 1.0. Please use
;     default value (1.03) for GL setup since the empirical
;     relationship is derived from that layer.
;   - For the active regions with mixed polarities, try to use 
;     smoothing or for 3D data use higher layers of the magnetogram.
;     to reduce complexity.
;
; MODIFICATION HISTORY:
;   Originally coded by Meng Jin@ AOSS, University of Michigan
;   v0.1 06/02/2014 Demo Version.
;   v0.2 06/05/2014 Added option to read binary data, remove the ASCII
;    template file dependance; added DemoMode and PlotRadius keywords.
;   v0.3 06/06/2014 Added option to calculate the GL orientation
;    according to the PIL.         
;   v0.4 11/24/2014 Added option to calculate the CME grid
;   v0.5 02/28/2015 Calculate the total magnetic field along the PIL
;   v0.6 03/10/2015 Added option for determining the AR parameter 
;    information.
;   v0.7 03/16/2015 Added option to manually set the radius of the
;   flux rope
;   v0.8 03/17/2015 Added two empirical equations to determine the
;   poloroidal flux of the GL flux rope.
;   v0.9 03/20/2015 Improved the coefficients of the empirical
;   equations.
;   v1.0 07/31/2015 Fixed reading new format FDIPS output. Fixed
;   PlotRadius issue. FDIPS output resolution can be flexible now.
;   2D Input option. Run without SSW option.
;   08/04/2015 G.Toth added optional nSmooth parameter for smoothing.
;      Added optional CMESpeed parameter to set the CME speed.
;      Added optional FILE parameter to set the name of the input
;      file. Show centerpoint of CME Source region with green.
;      Initialize !MOUSE.button before starting the cursor loop.
;      Removed DELVARX calls. Failed without solarsoft and not needed.
;------------------------------------------------------------------------

;Setup the color mode and a better IDL font.
  device,decomposed=0
  !p.font=1

;Setup common block for BATSRUS/Idl
  common file_head, $
     headline, it, time, gencoord, ndim, neqpar, nw, nx, eqpar, variables, $
     rbody

;Setup help option
  if not keyword_set(help) then help=0
  if help then begin
       PRINT,'Use:  IDL> SWMF_GLSETUP, [/DemoMode, PlotRadius=PlotRadius, /usepil, $ '
       PRINT,'                         /CMEGrid, ARMag=ARMag, GLRadius=GLRadius, $ '
       PRINT,'                         SizeFactor=SizeFactor, nSmooth=nSmooth, $ '
       PRINT,'                         File=File, CMEspeed=CMEspeed, /help] '
       RETURN
  endif

;Turn on/off demo mode. With the demo mode on, 
;the pre-saved 2D data will be read instead
;of reading from 3D data which is much more time consuming
  if not keyword_set(DemoMode) then DemoMode=0

;set default for nSmooth
  if not keyword_set(nSmooth) then nSmooth=0

;set default option for ARMag
  if not keyword_set(ARMag) then ARMag=1
  if ARMag ne 1 and ARMag ne 2 then begin
     print, 'ARMag can only be 1 or 2!'
     print, 'Use Default option: ARMag = 1'
     ARMag=1
  endif

;Read Observed CME speed.
  if not keyword_set(CMESpeed) then begin
     CMESpeed=0.0
     read,prompt='Please Input the Observed CME Speed (km/s): ',CMESpeed
  endif

;If GLRadius is not given then it is set to PIL_length/SizeFactor
  if not keyword_set(SizeFactor) then SizeFactor = 26.25

;Set DemoMode
  if keyword_set(DemoMode) and keyword_set(PlotRadius) then begin
     print,'DemoMode and PlotRadius Cannot be Used at the same time!'
     print,'Play DemoMode...'
     PlotRadius=1.03
  endif

;Read the SWMF input magnetic field
  if not DemoMode then begin
     if not keyword_set(file) then begin
        file=''
        read, prompt='Input file name (containing magnetic field data): ',file
     endif

     gettype, file, filetype, npictinfile
     openfile, 10, file, filetype
     get_pict, 10, file, filetype, 1, x, var, error 
     close, 10
     
     if gencoord then begin
        print, 'file '+file+' should contain a regular grid'
        retall
     endif

     case ndim of
        2:begin
           if not keyword_set(PlotRadius) then  PlotRadius=1.00
           if variables(0) ne "Longitude" or variables(1) ne "Latitude" or $
              variables(2) ne 'Br' then begin
              print, 'variables should be Longitude Latitude Br!'
              retall
           endif
           longitude = x(*,*,0)*!dtor
           latitude  = x(*,*,1)*!dtor
           br_field = var(*,*,0)
           bt_field = abs(br_field)
           nlon = nx[0]
           nlat = nx[1]
           bphi_field = fltarr(nlon, nlat)
           btheta_field = fltarr(nlon, nlat)
        end

        3:begin
           if not keyword_set(PlotRadius) then  PlotRadius=1.03
           if variables(0) ne "Radius" or variables(1) ne "Longitude" or $
              variables(2) ne "Latitude" or variables(3) ne 'Br' then begin
              print, 'variables should be Radius Longitude Latitude Br!'
              retall
           endif
           radius    = x(*,0,0,0)    
           longitude = x(0,*,*,1)
           latitude  = x(0,*,*,2)
           nlon = nx[1] - 1
           nlat = nx[2]
                                ; find index for the cut
           d = abs(radius - PlotRadius)
           icut = where( d eq min(d) )
           br_field     = var(icut,*,*,0)
           bphi_field   = var(icut,*,*,1)
           btheta_field = var(icut,*,*,2)
           bt_field = sqrt(br_field^2  + bphi_field^2 + btheta_field^2)

           br_field = reform(br_field[0,0:nlon-1,*])
           bphi_field = reform(bphi_field[0,0:nlon-1,*])
           btheta_field = reform(btheta_field[0,0:nlon-1,*])
           bt_field = reform(bt_field[0,0:nlon-1,*])
           longitude = reform(longitude[0,0:nlon-1,*])
           latitude = reform(latitude[0,0:nlon-1,*])
        end
        else: begin
           print, 'ndim=', ndim, ' should be 2 or 3'
           retall
        end
     endcase
  endif
    
;Restore pre-saved 2D data for the demo mode
  if DemoMode then begin
     restore,'magneticfield_1.0.sav'
  endif

  if nSmooth gt 1 then begin
     br_field     = smooth(br_field, nsmooth)
     btheta_field = smooth(btheta_field, nsmooth)
     bphi_field   = smooth(bphi_field, nsmooth)
  endif

  bt_field = sqrt(bphi_field^2+btheta_field^2+br_field^2)

;Display the magnetogram and let user interactively select the CME source region. The
;procedure to select is:
; 1. Click the CME source region of positive polarity with 'left' button of mouse
; 2. Click the CME source region of negative polarity with 'right' button of mouse
;
;Note that the user can click anywhere inside the active region. However, click closer
;to the center of the positive/negative patterns is recommended.
;
;Note the solar latitude is expressed in pixel due to the non-uniform spacing. The latitude
;is uniform in sin(latitude). This will be changed in the future to degree. 
  
  br_field_show=br_field
  index=where(br_field lt -20)
  br_field_show[index]=-20
  index=where(br_field gt 20)
  br_field_show[index]=20

  window,2,xs=1200,ys=1200.*float(nlat)/float(nlon)*4./3.
  loadct,0
  contour,br_field_show,min=-20,max=20,charsize=3,$
          title='SWMF Input Magnetogram (R ='$
          +strtrim(PlotRadius,2)+' Rs)',xtitle='Solar Longitude (Pixel)',$
          ytitle='Solar Latitude (Pixel)',/fill,nlevels=60,/iso,xstyle=1,ystyle=1
  
  print,'Please Select the CME Source Region (POSITIVE) with the left button'
  loadct,39
  !MOUSE.button = 0
  while(!MOUSE.button ne 1) do begin
     cursor,xPositiveSelect,yPositiveSelect,/data,/down
     if br_field[xPositiveSelect,yPositiveSelect] lt 0 then begin
        print,'Negative Polarity! Please Select POSITIVE Polarity!'   
        !MOUSE.button=0
     endif else begin
        plots,xPositiveSelect,yPositiveSelect,/data,psym=-2,color=250 
     endelse
  endwhile
  print,'Positive Source Region Selected: ',round(xPositiveSelect),round(yPositiveSelect)
  print,'Please Select the CME Source Region (NEGATIVE) with the right button'
  while(!MOUSE.button ne 4) do begin
     cursor,xNegativeSelect,yNegativeSelect,/data,/down
     if br_field[xNegativeSelect,yNegativeSelect] gt 0 then begin
        print,'Positive Polarity! Please Select NEGATIVE Polarity!'   
        !MOUSE.button=0
     endif else begin
        plots,xNegativeSelect,yNegativeSelect,/data,psym=-2,color=50
     endelse
  endwhile
  print,'Negative Source Region Selected: ',round(xNegativeSelect),round(yNegativeSelect)

;Set the box size. This box size is used to search the weight center around the
;selected positive/negative points in the interactive selection. It maybe increased
;for larger active regions.
  boxsize = round(long(16)*nlon/360)

;Calculate the weighted center for the positive polarity.
  xPositiveWeight=0.
  yPositiveWeight=0.
  TotalPositiveFlux=0.
  for i=xPositiveSelect-boxsize/2,xPositiveSelect+boxsize/2 do begin
     for j=yPositiveSelect-boxsize/2,yPositiveSelect+boxsize/2 do begin
        if br_field[i,j] gt 0 then begin
           xPositiveWeight=xPositiveWeight+br_field[i,j]*i
           yPositiveWeight=yPositiveWeight+br_field[i,j]*j
           TotalPositiveFlux=TotalPositiveFlux+br_field[i,j]
        endif
     endfor
  endfor
  xPositiveWeight=xPositiveWeight/TotalPositiveFlux
  yPositiveWeight=yPositiveWeight/TotalPositiveFlux

;Calculate the weighted center for the negative polarity.
  xNegativeWeight=0.
  yNegativeWeight=0.
  TotalNegativeFlux=0.
  for i=xNegativeSelect-boxsize/2,xNegativeSelect+boxsize/2 do begin
     for j=yNegativeSelect-boxsize/2,yNegativeSelect+boxsize/2 do begin
        if br_field[i,j] lt 0 then begin
           xNegativeWeight=xNegativeWeight+br_field[i,j]*i
           yNegativeWeight=yNegativeWeight+br_field[i,j]*j
           TotalNegativeFlux=TotalNegativeFlux+br_field[i,j]
        endif
     endfor
  endfor
  xNegativeWeight=xNegativeWeight/TotalNegativeFlux
  yNegativeWeight=yNegativeWeight/TotalNegativeFlux

;Plot the weighted centers on the magnetogram.
  
  device,decomposed=0
  loadct,0
  contour,br_field_show,min=-20,max=20,charsize=3,title='SWMF Input Magnetogram (R ='$
          +strtrim(PlotRadius,2)+' Rs)',xtitle='Solar Longitude (Pixel)',$
          ytitle='Solar Latitude (Pixel)',/fill,nlevels=60,/iso,xstyle=1,ystyle=1

  loadct,39
  plots,xPositiveWeight,yPositiveWeight,/data,psym=-2,color=250
  plots,xNegativeWeight,yNegativeWeight,/data,psym=-2,color=50

;Calculate the GL flux rope orientation from the two weighted points.
  r1=[xPositiveWeight-xNegativeWeight,yPositiveWeight-yNegativeWeight]
  r1=r1/sqrt(r1[0]^2+r1[1]^2)
  r2=[-1.0,0.0]
  GL_Orientation=acos(r1[0]*r2[0]+r1[1]*r2[1])*180/3.1415926
  if r1[1] lt 0 then begin
     GL_Orientation=360-GL_Orientation
  endif

;Extract the profile along the two weighted centers in order to determine the 
;center of the flux rope.
  aa=(yPositiveWeight-yNegativeWeight)/(xPositiveWeight-xNegativeWeight)
  bb=yPositiveWeight-aa*xPositiveWeight
  if abs(xPositiveWeight-xNegativeWeight) gt abs(yPositiveWeight-yNegativeWeight) then begin
     xProfile=min([xPositiveWeight,xNegativeWeight])+indgen(round(abs(xPositiveWeight-xNegativeWeight))+1)
     yProfile=round(aa*xProfile+bb)
     xProfile=round(xProfile)
  endif else begin
     yProfile=min([yPositiveWeight,yNegativeWeight])+indgen(round(abs(yPositiveWeight-yNegativeWeight))+1)
     xProfile=round((yProfile-bb)/aa)
     yProfile=round(yProfile)
  endelse

  nProfile=n_elements(xProfile)
  magProfile=fltarr(nProfile)

  for i=0,nProfile-1 do begin
     magProfile[i]=br_field[xProfile[i],yProfile[i]]
  endfor

  temp=min(abs(magProfile),index)
  plots,xProfile[index],yProfile[index],/data,psym=-2,color=150
  GL_Latitude=Latitude[xProfile[index],yProfile[index]]
  GL_Longitude=Longitude[xProfile[index],yProfile[index]]
  GL_Latitude=GL_Latitude*180./3.1415926
  GL_Longitude=GL_Longitude*180./3.1415926
  if ndim eq 2 then GL_Longitude = GL_Longitude + param 
  ar_center=[xProfile[index],yProfile[index]]


;Calculate the gradient of the Br field
  ddx=(shift(br_field,-1,0)-shift(br_field,1,0))/2.
  ddy=(shift(br_field,0,-1)-shift(br_field,0,1))/2.
  ddx[0,*]=br_field[1,*]-br_field[0,*]
  ddx[nlon-1,*]=br_field[nlon-1,*]-br_field[nlon-2,*]
  ddy[*,0]=br_field[*,1]-br_field[*,0]
  ddy[*,nlat-1]=br_field[*,nlat-1]-br_field[*,nlat-2]
  br_field_gradient=sqrt(ddx^2+ddy^2)

;Cell size is used to divide the magnetogram to sub regions in order to determine
;the PIL. 
  cell_size=1

;Setup the threshold for selecting cells near the PIL.
  flux_threshold=1.0

;Calculate the Bitmap (1/0) for determining the PIL.
  M=nlon/cell_size
  N=nlat/cell_size
  bitmap=fltarr(nlon,nlat)
  bitmap[*,*]=0.0
  for i=0,M-2 do begin
     for j=0,N-2 do begin
        index1=where(br_field[i*cell_size:(i+1)*cell_size,j*cell_size:(j+1)*cell_size] $
                     lt -flux_threshold)
        index2=where(br_field[i*cell_size:(i+1)*cell_size,j*cell_size:(j+1)*cell_size] $
                     gt flux_threshold)
        if index1[0] ne -1 and index2[0] ne -1 then begin
           bitmap[i*cell_size:(i+1)*cell_size,j*cell_size:(j+1)*cell_size]=1.0
        endif
     endfor
  endfor  

;Setup Bitmap for magnetic gradient
  bitmap_gradient=fltarr(nlon,nlat)
  bitmap_gradient[*,*]=1.0
  bitmap_gradient[where(br_field_gradient lt 0.5)]=0.0

;Distance cut-off for determining the PIL. 
  Dis_threshold=round(long(8)*nlon/360)
  DisCenter=fltarr(nlon,nlat)
  for i=0,nlon-1 do begin
     for j=0,nlat-1 do begin
        DisCenter[i,j]=sqrt((i-xProfile[index])^2+(j-yProfile[index])^2)
     endfor
  endfor

  dismap=fltarr(nlon,nlat)
  dismap[where(Discenter le Dis_threshold)]=1

;The final weighted map showing the PIL of the CME source region.
  wmap=bitmap*br_field*bitmap_gradient*dismap

;At this moment, the PIL length is represented by degree and does not 
;take into account the effect of different latitude. It will be improved later. 
  PIL_Length=(n_elements(where(wmap lt 0))+n_elements(where(wmap gt 0)))/2.*360./nlon

;Calculate the AR magnetic field strength in order to determine the
;right poloidal flux needed.
  showpoints=where(wmap ne 0)
  NN=n_elements(showpoints)
  pillines=intarr(NN,2)
  bt_pil=0.0
  RegionSize_ARMag=round(long(8)*nlon/360)
  for i=0,NN-1 do begin
     y_show=floor(showpoints[i]/nlon)
     x_show=showpoints[i]-(y_show*nlon)
     pillines[i,*]=[x_show,y_show]
     bt_pil=bt_pil+bt_field[x_show,y_show]
  endfor

  bt_pil=bt_pil/nn
  br_ar=mean(abs(br_field[ar_center[0]-RegionSize_ARMag/2:ar_center[0]+RegionSize_ARMag/2,$
                          ar_center[1]-RegionSize_ARMag/2:ar_center[1]+RegionSize_ARMag/2]))

;Showing the PIL
  showpoints=where(wmap gt 0)
  NN=n_elements(showpoints)
  for i=0,NN-1 do begin
     y_show=floor(showpoints[i]/nlon)
     x_show=showpoints[i]-(y_show*nlon)
     plots,x_show,y_show,psym=-1,color=200
  endfor

;Calculate the orientation of the flux rope according to PIL 
;(make it vertial to PIL).
  if keyword_set(USEPIL) then begin
     Dis_threshold_s=4
     dismap[where(Discenter gt Dis_threshold_s)]=0
     wmap_s=bitmap*br_field*bitmap_gradient*dismap
     PILpoints=where(wmap_s gt 0)
     MM=n_elements(PILpoints)
     PIL_x=fltarr(MM)
     PIL_y=fltarr(MM)
     for i=0,MM-1 do begin
        PIL_y[i]=floor(PILpoints[i]/nlon)
        PIL_x[i]=PILpoints[i]-(PIL_y[i]*nlon)
     endfor
     
     PIL_xx=PIL_x[sort(PIL_x)]
     PIL_yy=PIL_y[sort(PIL_x)]
     PIL_fit=ladfit(PIL_xx,PIL_yy,/double)  
     aa_PIL=-1./PIL_fit[1]
     if r1[0] lt 0 then begin 
        r3=[-1.,-aa_PIL]
     endif else begin
        r3=[1.,aa_PIL]
     endelse
     r3=r3/sqrt(r3[0]^2+r3[1]^2)
     GL_Orientation_s=acos(r3[0]*r2[0]+r3[1]*r2[1])*180/3.1415926
     if r3[1] lt 0 then begin
        GL_Orientation_s=360-GL_Orientation_s
     endif
     GL_Orientation=GL_Orientation_s
  endif

;Calculate the poloidal flux needed for the observed CME velocity.
  if ARMag eq 1 then begin
     GL_poloidal=(CMESpeed-655.00575+34.468279*br_ar-862.)/111.19675   
  endif else begin
     GL_poloidal=((CMESpeed+478.)/(33.2875/bt_pil)^0.494-655.00575)/111.19675
  endelse

;Relationship between the PIL length and the GL flux rope Radius.   
;This factor is now based on the 2011 March 7 CME. More tests  
;are needed in order to get a more precise value.  
  if not keyword_set(GLRadius) then GLRadius=PIL_Length/SizeFactor

;Relationship between the GL Poloidal flux and GL Bstrength.
  GL_Bstrength=(GL_poloidal-0.073579605)/(21.457435*GLRadius^4)

;Calculate the CME grid refinement parameters based on the flux rope
;location and size.                                                
  if keyword_set(CMEGrid) then begin
     CMEbox_Start=[1.1,GL_Longitude-40.*GLRadius,GL_Latitude-20.*GLRadius]
     CMEbox_End=[20.0,GL_Longitude+40.*GLRadius,GL_Latitude+20.*GLRadius]
  endif

;Recommended GL flux rope parameters
  Distance = 1.8
  Stretch  = 0.6

  print,'========================================'
  print,'The Recommended GL FLux Rope Parameters'
  print,'========================================'
  print,FORMAT='(A20,5X,F6.2)','Latitude: ',GL_Latitude
  print,FORMAT='(A20,5X,F6.2)','Longitude: ',GL_Longitude
  print,FORMAT='(A20,5X,F6.2)','Orientation: ',GL_Orientation
  print,FORMAT='(A20,5X,F6.2)','Radius: ', GLRadius
  print,FORMAT='(A20,5X,F6.2)','Bstrength: ',GL_Bstrength
  print,FORMAT='(A20,5X,F6.2)','Stretch (FIXED): ',Stretch
  print,FORMAT='(A20,5X,F6.2)','Distance (FIXED): ',Distance
  print,FORMAT='(A20,5X,F6.2)','Height [Rs]: ', $
        GLRadius + Distance - Stretch - 1.0
  print,FORMAT='(A20,5X,F6.2)','Angular size [deg]: ', $
        2*GLRadius/Distance/!dtor
  print,FORMAT='(A20,5X,F6.2)','Poloidal flux: ', GL_poloidal
  print,'-----------------------------------------'

  if keyword_set(CMEGrid) then begin
     print,'=========================================='
     print,'The Recommended Grid Refinement Parameters'
     print,'=========================================='
     print,FORMAT='(A20,5X,F6.2)','R_Start: ', CMEbox_start[0]
     print,FORMAT='(A20,5X,F6.2)','R_End: ', CMEbox_end[0]
     print,FORMAT='(A20,5X,F6.2)','Longitude_Start: ', CMEbox_start[1]
     print,FORMAT='(A20,5X,F6.2)','Longitude_End: ', CMEbox_end[1]
     print,FORMAT='(A20,5X,F6.2)','Latitude_Start: ', CMEbox_start[2]
     print,FORMAT='(A20,5X,F6.2)','Latitude_End: ', CMEbox_end[2]
     print,'-----------------------------------------'
  endif

;The region size is used to cover the whole area of active region in
;order to show a zoom-in image. Shorter RegionSize for near-Limb
;regions when needed.

  RegionSize=round(long(50)*nlon/360)

;Display the zoom-in image of the active region with weighted centers and PIL.
  window,3,xs=800,ys=800

  device,decomposed=0
  loadct,0
  contour,br_field_show[xProfile[index]-RegionSize/2:xProfile[index]+RegionSize/2,$
                        yProfile[index]-RegionSize/2:yProfile[index]+RegionSize/2],$
          min=-20,max=20,charsize=3,title='CME Source Region (R ='$
          +strtrim(PlotRadius,2)+' Rs)',xtitle='Solar Longitude (Degree)',$
          ytitle='Solar Latitude (Pixel)',/fill,nlevels=60,/iso,xstyle=1,ystyle=1

  loadct,39
  plots,xPositiveWeight-(xProfile[index]-RegionSize/2),yPositiveWeight-(yProfile[index]-RegionSize/2),$
        /data,psym=-2,color=250,symsize=3,thick=3
  plots,xNegativeWeight-(xProfile[index]-RegionSize/2),yNegativeWeight-(yProfile[index]-RegionSize/2),$
        /data,psym=-2,color=50,symsize=3,thick=3
  plots,RegionSize/2,RegionSize/2,/data,psym=-2,color=150,symsize=3,thick=3
  for i=0,NN-1 do begin
     y_show=floor(showpoints[i]/nlon)
     x_show=showpoints[i]-(y_show*nlon)
     y_show=y_show-(yProfile[index]-RegionSize/2)
     x_show=x_show-(xProfile[index]-RegionSize/2)
     plots,x_show,y_show,psym=-1,color=200,symsize=3,thick=3
  endfor

  !mouse.button=0
end
