C  Program for interpolation of NO2 AMF from look-up tables
C  Version 1.0 created on 13/01/2010 by BIRA-IASB
C  Contact: François Hendrick (franch@oma.be)
	PROGRAM interpolate_NO2_AMF

      INTEGER n_wv,n_lat,n_month,month_nbr,flag_am_pm,n_sza_max
	INTEGER n_sza_int,n_alb,albedo_flag,screen_flag,year
	PARAMETER (n_wv=6,n_lat=18,n_month=12,n_sza=14,n_sza_max=500000)
	PARAMETER (n_alb=2)
	REAL*8 wv_vec(n_wv),lat_vec(n_lat),sza_vec(n_sza),alb_vec(n_alb)
	REAL*8 no2_amftbl_sr(n_alb,n_wv,n_lat,n_month,n_sza)
	REAL*8 no2_amftbl_ss(n_alb,n_wv,n_lat,n_month,n_sza)
	REAL*8 no2_amftbl(n_alb,n_wv,n_lat,n_month,n_sza)
	REAL*8 daynb_vec(n_month),DT(3),UT,LT
	REAL*8 wvl_int,lat_int,daynb_int(n_sza_max),sza_int(n_sza_max)
	REAL*8 ut_int(n_sza_max),year_int(n_sza_max)
	REAL*8 dec_day_int(n_sza_max)
	REAL*8 no2_amftbl_alb(n_wv,n_lat,n_month,n_sza)
	REAL*8 no2_amftbl_dn(n_wv,n_lat,n_sza),no2_amftbl_lat(n_wv,n_sza)
	REAL*8 no2_amftbl_wv(n_sza),no2_amftbl_sza_int(500),long_int
	REAL*8 long_int_out,alb_int,no2_amf
	REAL*8 daynb_int_out,lat_int_out,YWORK(n_sza),Y2(n_sza),ycalc,hh
	REAL*8 alb_int_1,alb_int_2,alb_int_3,alb_int_4,alb_int_5
	REAL*8 alb_int_6,alb_int_7,alb_int_8,alb_int_9,alb_int_10
	REAL*8 alb_int_11,alb_int_12,long_alb_vec(360),lat_alb_vec(180)
	REAL*8 alb_mat_lat(360),alb_mat(360,180)
	CHARACTER*13 month_str
	CHARACTER*30 sza_file

C  Wavelength vector
	DATA wv_vec /
     $ 350.d0, 390.d0, 430.d0, 470.d0, 510.d0, 550.d0/
C  Latitude vector
	DATA lat_vec /
     $ -85.d0, -75.d0, -65.d0, -55.d0, -45.d0, -35.d0, -25.d0, 
     $ -15.d0, -5.d0, 5.d0, 15.d0, 25.d0, 35.d0, 45.d0, 55.d0, 65.d0,
     $  75.d0, 85.d0/
C  SZA vector
	DATA sza_vec /
     $ 10.d0, 30.d0, 50.d0, 70.d0, 80.d0, 82.5d0, 85.d0, 
     $ 86.d0, 87.d0, 88.d0, 89.d0, 90.d0, 91.d0, 92.d0/
C  Albedo vector
	DATA alb_vec /
     $ 0.d0, 1.d0/
C  Day number vector
      DATA daynb_vec /
     $ 15.d0, 46.d0, 74.d0, 105.d0, 135.d0, 166.d0, 196.d0, 
     $ 227.d0, 258.d0, 288.d0, 319.d0, 349.d0/


      
C  Read NO2 AMF interpolation input file
      OPEN(1,file='input_file_no2_amf.dat')
	DO i=1,3
		READ(1,*)
	ENDDO
      READ(1,*) wvl_int
	READ(1,*)
	READ(1,*) lat_int
	READ(1,*)
	READ(1,*) long_int
	READ(1,*)
	READ(1,*) albedo_flag
	READ(1,*) 
	READ(1,*) alb_int
	READ(1,*) 
      READ(1,100) sza_file
	READ(1,*) 
      READ(1,*) screen_flag
	CLOSE(1)
  100 format(a30)



C read SZA file

      OPEN(1,file=sza_file)
	DO i=1, n_sza_max
        READ(1,*,END=19) year_int(i),daynb_int(i),sza_int(i)
	  ut_int(i)= (daynb_int(i)-INT(daynb_int(i)))*24
	  dec_day_int(i)=daynb_int(i)
	  daynb_int(i)=INT(daynb_int(i))
c	write(*,*) year_int(i),daynb_int(i),sza_int(i), ut_int(i)
c	stop
	ENDDO
   19  n_sza_int = i-1
      CLOSE(1)

	

C Check input parameters values


	IF(wvl_int.LT.350.d0.OR.wvl_int.GT.550.d0) THEN
	WRITE(*,*)' !! PLEASE CHOOSE WAVELENGTH IN THE 350-550 NM RANGE 
     $!!'
	STOP
	ENDIF

      IF(lat_int.LT.-90.d0.OR.lat_int.GT.90.d0) THEN
	WRITE(*,*)' !! PLEASE CHOOSE LATITUDE IN THE -90deg - +90deg 
     $RANGE !!'
	STOP
	ENDIF

	IF(long_int.LT.-180.d0.OR.long_int.GT.180.d0) THEN
	WRITE(*,*)' !! PLEASE CHOOSE LONGITUDE IN THE -180deg - +180deg 
     $RANGE !!'
	STOP
	ENDIF

C Latitude check

      lat_int_out=lat_int

      IF(lat_int.LT.-85.d0) THEN
		lat_int=-84.99d0
	ELSEIF(lat_int.GT.85.d0) THEN
	    lat_int=84.99d0
	ENDIF

C Longitude check
      long_int_out=long_int
	IF(long_int.LT.-179.5d0) THEN
		long_int=-179.5d0
	ELSEIF(long_int.GT.179.5d0) THEN
	    long_int=179.5d0
	ENDIF

	

C  Read NO2 AMF look-up tables
      OPEN(1,file='no2_amf_sr.dat')
	DO i=1,16
			READ(1,*)
	ENDDO
      DO i_alb=1,n_alb
	DO i_wv=1,n_wv
		DO i_lat=1,n_lat
			DO i=1,2
				READ(1,*)
			ENDDO
			DO i_month=1,12
				READ(1,*) month_str, month_nbr,
     $			(no2_amftbl_sr(i_alb,i_wv,i_lat,i_month,j),j=1,14,1)
c	   write(*,*) i_month,i_lat,i_wv
c	   write(*,*) (no2_amftbl_sr(i_alb,i_wv,i_lat,i_month,j),j=1,14,1)
			ENDDO
		ENDDO
	ENDDO
	ENDDO
	CLOSE(1)

	OPEN(1,file='no2_amf_ss.dat')
	DO i=1,16
			READ(1,*)
	ENDDO
      DO i_alb=1,n_alb
	DO i_wv=1,n_wv
		DO i_lat=1,n_lat
			DO i=1,2
				READ(1,*)
			ENDDO
			DO i_month=1,12
				READ(1,*) month_str, month_nbr,
     $			(no2_amftbl_ss(i_alb,i_wv,i_lat,i_month,j),j=1,14,1)
c	   write(*,*) i_month,i_lat,i_wv
c	   write(*,*) (no2_amftbl_ss(i_alb,i_wv,i_lat,i_month,j),j=1,14,1)
			ENDDO
		ENDDO
	ENDDO
	ENDDO
	CLOSE(1)

C Surface albedo
C Koelemeijer et al, JGR, 108, D2, 4070, 2003
C
C Define lat_alb and long_alt vectors

      DO i=1,180
	    lat_alb_vec(i)=-89.5+1*(i-1)
	ENDDO

	DO i=1,360
	    long_alb_vec(i)=-179.5+1*(i-1)
	ENDDO

	
C load albedo data files and interpolation at the latitude and longitude of the station
      IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_01_380nm.dat')
	ELSE
          OPEN(2,file='albedo_01_494nm.dat')
	ENDIF
		DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_1=ycalc

	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_02_380nm.dat')
	ELSE
          OPEN(2,file='albedo_02_494nm.dat')
	ENDIF
	    DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_2=ycalc
	
	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_03_380nm.dat')
	ELSE
          OPEN(2,file='albedo_03_494nm.dat')
	ENDIF
	    DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_3=ycalc

	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_04_380nm.dat')
	ELSE
          OPEN(2,file='albedo_04_494nm.dat')
	ENDIF
   	    DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_4=ycalc
	
	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_05_380nm.dat')
	ELSE
          OPEN(2,file='albedo_05_494nm.dat')
	ENDIF
	    DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_5=ycalc
	
	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_06_380nm.dat')
	ELSE
          OPEN(2,file='albedo_06_494nm.dat')
	ENDIF
	    DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_6=ycalc
	
	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_07_380nm.dat')
	ELSE
          OPEN(2,file='albedo_07_494nm.dat')
	ENDIF
	    DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_7=ycalc
	    
	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_08_380nm.dat')
	ELSE
          OPEN(2,file='albedo_08_494nm.dat')
	ENDIF
	    DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_8=ycalc
	
	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_09_380nm.dat')
	ELSE
          OPEN(2,file='albedo_09_494nm.dat')
	ENDIF
	    DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_9=ycalc
	
	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_10_380nm.dat')
	ELSE
          OPEN(2,file='albedo_10_494nm.dat')
	ENDIF
	    DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_10=ycalc
	
	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_11_380nm.dat')
	ELSE
          OPEN(2,file='albedo_11_494nm.dat')
	ENDIF
	    DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_11=ycalc
	
	IF(wvl_int.LT.450) THEN
          OPEN(2,file='albedo_12_380nm.dat')
	ELSE
          OPEN(2,file='albedo_12_494nm.dat')
	ENDIF
      	DO i=1,360
			READ(2,*) alb_mat(i,:)
		ENDDO
		CLOSE(2)
		DO i=1,360
			CALL inter(180,180,2,lat_int,lat_alb_vec,
     $		alb_mat(i,:),ycalc,hh)
			alb_mat_lat(i)=ycalc
		ENDDO
			CALL inter(360,360,2,long_int,long_alb_vec,
     $		alb_mat_lat(:),ycalc,hh)
			alb_int_12=ycalc

C  Write parameters values in the output file
      OPEN(1,file='no2_amf_output.dat')
		WRITE(1,201) '% NO2 AMF'
	    WRITE(1,202) '% Latitude (°)                    : ',
     $     lat_int_out
          WRITE(1,202) '% Longitude (°)                   : ',
     $     long_int_out
		WRITE(1,202) '% Wavelength (nm)                 : ', wvl_int
	    WRITE(1,205) '% Year  Fractional day   SZA (°)    NO2 AMF'
	

C day number loop
      DO dd=1,n_sza_int
      daynb_int_out=daynb_int(dd)
	

      IF(daynb_int(dd).LT.15.d0) THEN
		daynb_int(dd)=15.d0
	ELSEIF(daynb_int(dd).GT.349.d0) THEN
	    daynb_int(dd)=349.d0
	ENDIF

	IF(albedo_flag.EQ.1) THEN

C Choose the albedo value corresponding to the day number in the input file
      IF(daynb_int(dd).GE.1.and.daynb_int(dd).LT.32) THEN
	    alb_int=alb_int_1
	ELSEIF(daynb_int(dd).GE.32.and.daynb_int(dd).LT.60) THEN
	    alb_int=alb_int_2
	ELSEIF(daynb_int(dd).GE.60.and.daynb_int(dd).LT.91) THEN
	    alb_int=alb_int_3
	ELSEIF(daynb_int(dd).GE.91.and.daynb_int(dd).LT.121) THEN
	    alb_int=alb_int_4
	ELSEIF(daynb_int(dd).GE.121.and.daynb_int(dd).LT.152) THEN
	    alb_int=alb_int_5
	ELSEIF(daynb_int(dd).GE.152.and.daynb_int(dd).LT.182) THEN
	    alb_int=alb_int_6
	ELSEIF(daynb_int(dd).GE.182.and.daynb_int(dd).LT.213) THEN
	    alb_int=alb_int_7
	ELSEIF(daynb_int(dd).GE.213.and.daynb_int(dd).LT.244) THEN
	    alb_int=alb_int_8
	ELSEIF(daynb_int(dd).GE.244.and.daynb_int(dd).LT.274) THEN
	    alb_int=alb_int_9
	ELSEIF(daynb_int(dd).GE.274.and.daynb_int(dd).LT.305) THEN
	    alb_int=alb_int_10
	ELSEIF(daynb_int(dd).GE.305.and.daynb_int(dd).LT.335) THEN
	    alb_int=alb_int_11
	ELSEIF(daynb_int(dd).GE.335) THEN
	    alb_int=alb_int_12
	ENDIF

	ENDIF

C     Calculate local time from fractional day and choose the right LUT
	CALL jjdat(daynb_int(dd),year_int(dd),DT)
      UT=ut_int(dd)
	CALL loctime_fct(DT,-long_int,UT,LT)
      
	IF(LT.LE.12.d0) THEN
		no2_amftbl=no2_amftbl_sr
	ELSE 
          no2_amftbl=no2_amftbl_ss
	ENDIF

C     Interpolation on the albedo
      DO i_wn=1,n_wv
		DO i_lat=1,n_lat
			DO i_sza=1,n_sza
	           DO i_month=1,12
				 CALL inter(n_month,n_month,2,alb_int,alb_vec,
     $				no2_amftbl(:,i_wn,i_lat,i_month,i_sza),ycalc,hh)
				    no2_amftbl_alb(i_wn,i_lat,i_month,i_sza)=ycalc
	           ENDDO
			ENDDO
		ENDDO
	ENDDO
	

C     Interpolation on the day number
      DO i_alb=1,n_alb
      DO i_wn=1,n_wv
		DO i_lat=1,n_lat
			DO i_sza=1,n_sza
			CALL SPLINE(daynb_vec,no2_amftbl_alb(i_wn,i_lat,:,i_sza),
     $      n_month,YWORK,Y2)
	        CALL SPLINT(daynb_vec,no2_amftbl_alb(i_wn,i_lat,:,i_sza),
     $      Y2,n_month,daynb_int(dd),ycalc)
				    no2_amftbl_dn(i_wn,i_lat,i_sza)=ycalc
			ENDDO
		ENDDO
	ENDDO
	ENDDO


C    Interpolation on the latitude
      DO i_wn=1,n_wv
		DO i_sza=1,n_sza
			CALL inter(n_lat,n_lat,2,lat_int,lat_vec,
     $		    no2_amftbl_dn(i_wn,:,i_sza),ycalc,hh)			
		no2_amftbl_lat(i_wn,i_sza)=ycalc
		ENDDO
	ENDDO

C    Interpolation on the wavelength
	DO i_sza=1,n_sza
		CALL SPLINE(wv_vec,no2_amftbl_lat(:,i_sza),n_wv,YWORK,Y2)
     	    CALL SPLINT(wv_vec,no2_amftbl_lat(:,i_sza),Y2,n_wv,wvl_int,
     $     ycalc)
		no2_amftbl_wv(i_sza)=ycalc
	ENDDO

C    Interpolation on the SZA
	    CALL SPLINE(sza_vec,no2_amftbl_wv(:),n_sza,YWORK,Y2)
	    CALL SPLINT(sza_vec,no2_amftbl_wv(:),Y2,n_sza,
     $      			sza_int(dd),ycalc)
		no2_amf=ycalc

	year=INT(year_int(dd))
	WRITE(1,206) year, dec_day_int(dd), sza_int(dd),no2_amf
	
      IF(screen_flag.EQ.1) THEN
	WRITE(*,207) year, dec_day_int(dd), sza_int(dd),no2_amf
	WRITE(*,*) ''
	ENDIF
	

      ENDDO

	CLOSE(1)

201   FORMAT(a9)
202   FORMAT(a36,F6.2)
205   FORMAT(a43)
206   FORMAT(I6,x,F12.5,x,F11.2,x,F10.2)
207   FORMAT(I6,x,F12.5,x,F9.2,x,F8.2)




	END
C  End of interplotate_NO2_AMF
*******************************************************************************
*******************************************************************************
       SUBROUTINE inter( dim, npoints, itype, arg,
     $     xarr, yarr, ynew, hh)
*
*     Interpolates at the x-point arg from x-value array xarr and
*     y-value array yarr. xarr and yarr are expected to have
*     ascending arguments
*     Input variables:
*     dim       Array dimension of xarr and yarr
*     npoints   No. points in arrays xarr and yarr
*     itype     Interpolation type
*     arg       Interpolation argument
*     xarr      array of x values
*     yarr      array of y values
*
*     Output variables:
*     ynew      Interpolated function value at arg
*     hh        gradient or scale height value  
*
      INTEGER dim
      REAL*8  arg, hh, xarr(dim), yarr(dim), ynew
*
C       write(*,*) xarr(:)
C	 write(*,*) yarr(:)
C	write(*,*) npoints
C	write(*,*) dim
C	write(*,*) arg
	 IF ( arg .LE. xarr(npoints) .AND. arg .GE. xarr(1)) THEN
         DO 10 iq = 1 , npoints-1
           IF ( arg .LT. xarr(iq+1) .AND. arg .GE. xarr(iq)) ip = iq
   10    CONTINUE
C        write(*,*) ip
         IF ( arg .EQ. xarr(npoints)) ip = npoints-1
       ELSEIF ( arg .GT. xarr(npoints)) THEN
         ip = npoints - 1
       ELSEIF ( arg .LT. xarr(1)) THEN
         ip = 1
       ENDIF

C	write(*,*) ip 
*
*     Interpolate function value at arg from data points ip to ip+1
*
*
*     exponential interpolation
*
C       IF ( itype .EQ. 1 ) THEN
C          IF ( yarr(ip+1) .EQ. yarr(ip) ) THEN
C             hh   = 0.0
C             ynew = yarr(ip)
C          ELSE
C             hh   = -( xarr(ip+1) - xarr(ip) ) / 
C     $            LOG( yarr(ip+1) / yarr(ip))
C             ynew = yarr(ip) * EXP(- ( arg - xarr(ip) ) / hh )
C          ENDIF
*
*     linear interpolation
*
C       ELSEIF ( itype .EQ. 2 ) THEN
         hh = ( yarr(ip+1) - yarr(ip) ) / ( xarr(ip+1) - xarr(ip) )
         ynew = yarr(ip) + hh*( arg - xarr(ip))
*     
C      ENDIF
      RETURN
      END       
C******************************************************************************
      SUBROUTINE SPLINE(X,Y,N,WORK,Y2)
C******************************************************************************
C     Routine to calculate 2.nd derivatives of tabulated function
C     Y(i)=Y(Xi), to be used for cubic spline calculation.
C
      IMPLICIT NONE
      INTEGER N,I
      REAL(KIND=8)  X(N),Y(N),WORK(N),Y2(N)
      REAL(KIND=8)  SIG,P,QN,UN,YP1,YPN
C
      YP1=(Y(2)-Y(1))/(X(2)-X(1))
      YPN=(Y(N)-Y(N-1))/(X(N)-X(N-1))
      IF(YP1.GT.99.0E+30) THEN
          Y2(1)=0.0
          WORK(1)=0.0
      ELSE
          Y2(1)=-0.5E0
          WORK(1)=(3.0E0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO I=2,N-1
          SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
          P=SIG*Y2(I-1)+2.0E0
          Y2(I)=(SIG-1.0E0)/P
          WORK(I)=(6.0E0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     +             /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*WORK(I-1))/P
      ENDDO
      IF(YPN.GT.99.0E+30) THEN
          QN=0.0
          UN=0.0
      ELSE
          QN=0.5E0
          UN=(3.0E0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*WORK(N-1))/(QN*Y2(N-1)+1.0E0)
      DO I=N-1,1,-1
          Y2(I)=Y2(I)*Y2(I+1)+WORK(I)
      ENDDO
C
      RETURN
      END
C
C******************************************************************************
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
C******************************************************************************
C     Cubic spline calculation
C
      IMPLICIT NONE
      INTEGER KLO,KHI,N,K
      REAL(KIND=8)  XA(N),YA(N),Y2A(N)
      REAL(KIND=8)  X,Y,H,A,B
C
      KLO=1
      KHI=N
 1    IF(KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).GT.X) THEN
              KHI=K
          ELSE
              KLO=K
          ENDIF
          GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     +        ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.0E0
C
      RETURN
      END



C******************************************************************************
C Subroutine calculating apparent local time given
C date, longitude, and UT
C LAT=loctime_fct(DT,LONG,UT) with
C	DT=date in the format [dd mm yyyy]
C	LONG=longitude (+W, -E) in decimal format
C	UT=universal time in the format decimal
C******************************************************************************

	SUBROUTINE loctime_fct(DT,LONG,UT,LT)

	REAL*8 DT(3),DT2(3),LT_HHMMSS(3)
	REAL*8 LONG,HR,JD1,JD2,DN,GMA,ET,LT,UT
	PARAMETER(PI=3.14159)

	DT2(1)=31
	DT2(2)=12
	DT2(3)=DT(3)-1

      HR=UT
	HR=HR-(LONG/15)

c	write(*,*) HR

	CALL dfact(DT,JD1)
	CALL dfact(DT2,JD2)

	DN=JD1-JD2


	GMA=(2*PI*(DN-1))/365
      ET=(0.000075+0.001868*COS(GMA)-0.032077*SIN(GMA)-0.014615
     $     *COS(2*GMA)-0.04089*SIN(2*GMA))*229.18
c	write(*,*) ET/60

	LT=HR+(ET/60)

	LT_HHMMSS(1)=INT(LT)
	LT_HHMMSS(2)=INT((LT-LT_HHMMSS(1))*60)
      LT_HHMMSS(3)=INT((LT-LT_HHMMSS(1)-LT_HHMMSS(2)/60)*3600)

c	write(*,*) LT_HHMMSS

	RETURN
	END


	SUBROUTINE dfact(DT,JD)

	REAL*8 DT(3)
	REAL*8 JD

    
      IF(DT(2).EQ.1.OR.DT(2).EQ.2) THEN
		JD=365*DT(3)+DT(1)+31*(DT(2)-1)+INT((DT(3)-1)/4)-INT(0.75*
     $      (INT(DT(3)-1)/100)+1)
	ELSE
          JD=365*DT(3)+DT(1)+31*(DT(2)-1)-INT(0.4*DT(2)+2.3)+INT(DT(3)/4
     $      )-INT(0.75*(INT(DT(3)/100)+1))
	ENDIF
      
	RETURN
	END

	SUBROUTINE jjdat(jul,year,DT)

	REAl*8 DT(3),Init(12),Initj(12)
	REAL*8 jul,year,remainder

	DATA Init /
     $ 1.d0, 32.d0, 60.d0, 91.d0, 121.d0, 152.d0, 182.d0, 213.d0,
     $ 244.d0, 274.d0, 305.d0, 335.d0/

	Initj = Init
	remainder=MOD(year,4.d0)
	IF(remainder.EQ.0) THEN
      Initj(3:12) = Init(3:12) + 1
	ELSE
      Initj(3:12) = Init(3:12)
	ENDIF

	DO i=1,12
		IF(Initj(i).GT.jul) THEN 
			GOTO 20
	    ENDIF
	ENDDO
   20 i=i-1

	DT(1)=jul-Initj(i)+1
	DT(2)=i
	DT(3)=year

	RETURN
	END