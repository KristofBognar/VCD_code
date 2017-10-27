C  Program for interpolating O3 AMF from look-up tables
C  Version 1.0 created on 14/10/2008 by BIRA-IASB
C  Contact: François Hendrick (franch@oma.be)
	PROGRAM interpolate_O3_AMF

      INTEGER n_wv,n_lat,n_month,month_nbr,flag_am_pm,n_sza_max
	INTEGER n_sza_int,n_o3,n_alb,n_alt,ind_lat,cc,inter_flag
      PARAMETER (n_wv=5,n_lat=18,n_month=12,n_sza=13,n_sza_max=1000)
	PARAMETER (n_o3=10,n_alb=2,n_alt=2)
	PARAMETER (n_ind_low=10,n_ind_up=10)
	INTEGER ind_low(n_ind_low),ind_up(n_ind_low),albedo_flag,month_nb
	REAL*8 wv_vec(n_wv),lat_vec(n_lat),sza_vec(n_sza),alb_vec(n_alb)
	REAL*8 alt_vec(n_alt),o3_scd_effective(n_o3),o3_amf(n_sza_max)
	REAL*8 o3_amftbl(n_alt,n_alb,n_wv,n_lat,n_month,n_sza,n_o3)
	REAL*8 o3_scdtbl(n_alt,n_alb,n_wv,n_lat,n_month,n_sza,n_o3)
	REAL*8 daynb_vec(n_month),o3_vcd_effective(n_o3)
	REAL*8 wvl_int,lat_int,daynb_int,sza_int(n_sza_max),ycalc,hh
	REAL*8 scd_int(n_sza_max),alb_int,alt_int,o3_col_int
	REAL*8 vcd_int(n_sza_max),dycalc,long_int
	REAL*8 o3_vcd(n_o3),daynb_int_out,lat_int_out,o3_vcd_vec(n_o3)
	REAL*8 o3_amftbl_dn(n_alt,n_alb,n_wv,n_lat,n_sza,n_o3)
	REAL*8 o3_amftbl_alt(n_alb,n_wv,n_lat,n_sza,n_o3)
	REAL*8 o3_amftbl_alb(n_wv,n_lat,n_sza,n_o3)
	REAL*8 o3_amftbl_wv(n_lat,n_sza,n_o3)
	REAL*8 o3_amftbl_lat(n_sza,n_o3),long_int_out
	REAL*8 o3_scdtbl_dn(n_alt,n_alb,n_wv,n_lat,n_sza,n_o3)
	REAL*8 o3_scdtbl_alt(n_alb,n_wv,n_lat,n_sza,n_o3)
	REAL*8 o3_scdtbl_alb(n_wv,n_lat,n_sza,n_o3),alb_mat_lat(360)
	REAL*8 o3_scdtbl_wv(n_lat,n_sza,n_o3),alb_mat(360,180)
	REAL*8 o3_scdtbl_lat(n_sza,n_o3),o3_amf_effective(n_o3)
	REAL*8 o3_amftbl_sza_int(n_sza_max,n_o3),long_alb_vec(360)
	REAL*8 o3_scdtbl_sza_int(n_sza_max,n_o3),lat_alb_vec(180)
	REAL*8 YWORK(100),Y2(100),lat_vec_work(2)
	CHARACTER*12 str_buffer
	CHARACTER*30 sza_file

C  Wavelength vector
	DATA wv_vec /
     $ 440.d0, 475.d0, 510.d0, 545.d0, 580.d0/
C  Latitude vector
	DATA lat_vec /
     $ -85.d0, -75.d0, -65.d0, -55.d0, -45.d0, -35.d0, -25.d0, 
     $ -15.d0, -5.d0, 5.d0, 15.d0, 25.d0, 35.d0, 45.d0, 55.d0, 65.d0,
     $  75.d0, 85.d0/
C  Albedo vector
	DATA alb_vec /
     $ 0.d0, 1.d0/
C  Altitude vector
	DATA alt_vec /
     $ 0.d0, 4.d0/
C  SZA vector
	DATA sza_vec /
     $ 30.d0, 50.d0, 70.d0, 80.d0, 82.5d0, 85.d0, 86.d0, 87.d0, 
     $ 88.d0, 89.d0, 90.d0, 91.d0, 92.d0/
C  Day number vector
      DATA daynb_vec /
     $ 15.d0, 46.d0, 74.d0, 105.d0, 135.d0, 166.d0, 196.d0, 
     $ 227.d0, 258.d0, 288.d0, 319.d0, 349.d0/
C  O3 VCD vector
      DATA o3_vcd_vec /
     $ 125.d0, 175.d0, 225.d0, 275.d0, 325.d0, 375.d0, 425.d0, 
     $ 475.d0, 525.d0, 575.d0/

      
C  Read O3 AMF interpolation input file
      OPEN(1,file='input_file_o3_amf.dat')
	DO i=1,3
		READ(1,*)
	ENDDO
      READ(1,*) wvl_int
	READ(1,*)
	READ(1,*) daynb_int
	READ(1,*)
	READ(1,*) lat_int
	READ(1,*) 
	READ(1,*) long_int
	READ(1,*) 
	READ(1,*) albedo_flag
	READ(1,*) 
	READ(1,*) alb_int
	READ(1,*) 
	READ(1,*) alt_int
	READ(1,*) 
      READ(1,100) sza_file
	READ(1,*) 
      READ(1,*) inter_flag
	CLOSE(1)
  100 format(a30)


C read SZA file

      OPEN(1,file=sza_file)
	DO i=1, n_sza_max
	IF(inter_flag.EQ.2) THEN
        READ(1,*,END=19) sza_int(i),scd_int(i)
	ELSE
        READ(1,*,END=19) sza_int(i),vcd_int(i)
	ENDIF

	ENDDO	
c     number of sza for interpolation  
  19  n_sza_int = i-1
      CLOSE(1)
      
	o3_vcd=o3_vcd_vec	

C Check input parameters values


	IF(wvl_int.LT.440.d0.OR.wvl_int.GT.580.d0) THEN
	WRITE(*,*)' !! PLEASE CHOOSE WAVELENGTH IN THE 440-580 NM RANGE 
     $!!'
	STOP
	ENDIF

      IF(daynb_int.LT.1.d0.OR.daynb_int.GT.366.d0) THEN
	WRITE(*,*)' !! PLEASE CHOOSE DAY NUMBER IN THE 1-365 RANGE
     $ OR IN THE 1-366 RANGE FOR LEAP YEAR !!'
	STOP
	ENDIF

      IF(lat_int.LT.-90.d0.OR.lat_int.GT.90.d0) THEN
	WRITE(*,*)' !! PLEASE CHOOSE LATITUDE IN THE -90deg - +90deg 
     $RANGE !!'
	STOP
	ENDIF

	IF(sza_int(1).LT.30.d0.OR.sza_int(n_sza_int).GT.92.d0) THEN
	WRITE(*,*)' !! PLEASE CHOOSE SZA IN THE 30deg - 92deg 
     $RANGE !!'
	STOP
	ENDIF

	IF((alb_int.LT.0.d0.OR.alb_int.GT.1.d0).and.albedo_flag.eq.2) THEN
	WRITE(*,*)' !! PLEASE CHOOSE GROUND ALBEDO IN THE 0 - 1 
     $RANGE !!'
	STOP
	ENDIF

	IF(alt_int.LT.0.d0.OR.alt_int.GT.4.d0) THEN
	WRITE(*,*)' !! PLEASE CHOOSE ALTITUDE IN THE 0 - 4 km 
     $RANGE !!'
	STOP
	ENDIF

C Day number check
      daynb_int_out=daynb_int

      IF(daynb_int.LT.15.d0) THEN
		daynb_int=15.d0
	ELSEIF(daynb_int.GT.349.d0) THEN
	    daynb_int=349.d0
	ENDIF

C Latitude check
      lat_int_out=lat_int

      IF(lat_int.LT.-85.d0) THEN
		lat_int=-84.99d0
	ELSEIF(lat_int.GE.85.d0) THEN
	    lat_int=84.99d0
	ENDIF

	DO i=1,n_lat-1
		IF(lat_int.GE.lat_vec(i).AND.lat_int.LT.lat_vec(i+1)) THEN
			ind_lat=i
		ENDIF
	ENDDO

	IF(lat_int.GE.lat_vec(n_lat)) THEN
		ind_lat=18
	ENDIF

C Longitude check
      long_int_out=long_int
	IF(long_int.LT.-179.5d0) THEN
		long_int=-179.5d0
	ELSEIF(long_int.GT.179.5d0) THEN
	    long_int=179.5d0
	ENDIF

	IF(albedo_flag.EQ.1) THEN

C interpolate albedo value if Koelemeijer database is used
C Koelemeijer et al, JGR, 108, D2, 4070, 2003
C
C Define lat_alb and long_alt vectors

      DO i=1,180
	    lat_alb_vec(i)=-89.5+1*(i-1)
	ENDDO

	DO i=1,360
	    long_alb_vec(i)=-179.5+1*(i-1)
	ENDDO

	
C Convert day number in month number
      If(daynb_int.GE.1.and.daynb_int.LT.32) THEN
	    OPEN(1,file='albedo_01_494nm.dat')
	ELSEIF(daynb_int.GE.32.and.daynb_int.LT.60) THEN
	    OPEN(1,file='albedo_02_494nm.dat')
	ELSEIF(daynb_int.GE.60.and.daynb_int.LT.91) THEN
	    OPEN(1,file='albedo_03_494nm.dat')
	ELSEIF(daynb_int.GE.91.and.daynb_int.LT.121) THEN
	    OPEN(1,file='albedo_04_494nm.dat')
	ELSEIF(daynb_int.GE.121.and.daynb_int.LT.152) THEN
	    OPEN(1,file='albedo_05_494nm.dat')
	ELSEIF(daynb_int.GE.152.and.daynb_int.LT.182) THEN
	    OPEN(1,file='albedo_06_494nm.dat')
	ELSEIF(daynb_int.GE.182.and.daynb_int.LT.213) THEN
	    OPEN(1,file='albedo_07_494nm.dat')
	ELSEIF(daynb_int.GE.213.and.daynb_int.LT.244) THEN
	    OPEN(1,file='albedo_08_494nm.dat')
	ELSEIF(daynb_int.GE.244.and.daynb_int.LT.274) THEN
	    OPEN(1,file='albedo_09_494nm.dat')
	ELSEIF(daynb_int.GE.274.and.daynb_int.LT.305) THEN
	    OPEN(1,file='albedo_10_494nm.dat')
	ELSEIF(daynb_int.GE.305.and.daynb_int.LT.335) THEN
	    OPEN(1,file='albedo_11_494nm.dat')
	ELSEIF(daynb_int.GE.335) THEN
	    OPEN(1,file='albedo_12_494nm.dat')
	ENDIF

	DO i=1,360
			READ(1,*) alb_mat(i,:)
	ENDDO

	CLOSE(1)


	DO i=1,360
	    CALL inter(180,180,2,lat_int,lat_alb_vec,
     $	alb_mat(i,:),ycalc,hh)
	    alb_mat_lat(i)=ycalc
	ENDDO

	
	CALL inter(360,360,2,long_int,long_alb_vec,
     $	alb_mat_lat(:),ycalc,hh)
	    alb_int=ycalc


	ENDIF
      


C  Write parameters values on the screen
	    WRITE(*,*) ''
	    WRITE(*,202) 'Latitude                        : ',
     $     lat_int_out
	    WRITE(*,202) 'Longitude                       : ',
     $     long_int_out
	    WRITE(*,203) 'Day number                      : ', 
     $     daynb_int_out
		WRITE(*,202) 'Wavelength (nm)                 : ', wvl_int
	    WRITE(*,204) 'Albedo                          : ', alb_int
	    WRITE(*,204) 'Altitude of the station (km)    : ', alt_int
          WRITE(*,*) ''
	

C  Read O3 AMF look-up table
      IF(ind_lat.EQ.1) THEN
		OPEN(1,file='o3_amf_lat_01.dat')
	    OPEN(2,file='o3_amf_lat_02.dat')
	ELSEIF (ind_lat.EQ.2) THEN
	    OPEN(1,file='o3_amf_lat_02.dat')
	    OPEN(2,file='o3_amf_lat_03.dat')
	ELSEIF (ind_lat.EQ.3) THEN
	    OPEN(1,file='o3_amf_lat_03.dat')
	    OPEN(2,file='o3_amf_lat_04.dat')
	ELSEIF (ind_lat.EQ.4) THEN
	    OPEN(1,file='o3_amf_lat_04.dat')
	    OPEN(2,file='o3_amf_lat_05.dat')
	ELSEIF (ind_lat.EQ.5) THEN
	    OPEN(1,file='o3_amf_lat_05.dat')
	    OPEN(2,file='o3_amf_lat_06.dat')
	ELSEIF (ind_lat.EQ.6) THEN
	    OPEN(1,file='o3_amf_lat_06.dat')
	    OPEN(2,file='o3_amf_lat_07.dat')
	ELSEIF (ind_lat.EQ.7) THEN
	    OPEN(1,file='o3_amf_lat_07.dat')
	    OPEN(2,file='o3_amf_lat_08.dat')
	ELSEIF (ind_lat.EQ.8) THEN
	    OPEN(1,file='o3_amf_lat_08.dat')
	    OPEN(2,file='o3_amf_lat_09.dat')
	ELSEIF (ind_lat.EQ.9) THEN
	    OPEN(1,file='o3_amf_lat_09.dat')
	    OPEN(2,file='o3_amf_lat_10.dat')
	ELSEIF (ind_lat.EQ.10) THEN
	    OPEN(1,file='o3_amf_lat_10.dat')
	    OPEN(2,file='o3_amf_lat_11.dat')
	ELSEIF (ind_lat.EQ.11) THEN
	    OPEN(1,file='o3_amf_lat_11.dat')
	    OPEN(2,file='o3_amf_lat_12.dat')
	ELSEIF (ind_lat.EQ.12) THEN
	    OPEN(1,file='o3_amf_lat_12.dat')
	    OPEN(2,file='o3_amf_lat_13.dat')
	ELSEIF (ind_lat.EQ.13) THEN
	    OPEN(1,file='o3_amf_lat_13.dat')
	    OPEN(2,file='o3_amf_lat_14.dat')
	ELSEIF (ind_lat.EQ.14) THEN
	    OPEN(1,file='o3_amf_lat_14.dat')
	    OPEN(2,file='o3_amf_lat_15.dat')
	ELSEIF (ind_lat.EQ.15) THEN
	    OPEN(1,file='o3_amf_lat_15.dat')
	    OPEN(2,file='o3_amf_lat_16.dat')
	ELSEIF (ind_lat.EQ.16) THEN
	    OPEN(1,file='o3_amf_lat_16.dat')
	    OPEN(2,file='o3_amf_lat_17.dat')
	ELSEIF (ind_lat.GE.17) THEN
	    OPEN(1,file='o3_amf_lat_17.dat')
	    OPEN(2,file='o3_amf_lat_18.dat')
	ENDIF
      
		DO i_alt=1,n_alt
			DO i_alb=1,n_alb
				DO i_wv=1,n_wv
					DO i_month=1,n_month
						DO i_sza=1,n_sza
			
						  READ(1,*) 
			
						  DO i_o3=1,n_o3

						   READ(1,*) 		
     $            o3_scdtbl(i_alt,i_alb,i_wv,1,i_month,i_sza,i_o3),
     $            o3_amftbl(i_alt,i_alb,i_wv,1,i_month,i_sza,i_o3) 
						  ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	 
	    CLOSE(1)
      
		DO i_alt=1,n_alt
			DO i_alb=1,n_alb
				DO i_wv=1,n_wv
					DO i_month=1,n_month
						DO i_sza=1,n_sza
			
						  READ(2,*) 
			
						  DO i_o3=1,n_o3

						   READ(2,*) 		
     $            o3_scdtbl(i_alt,i_alb,i_wv,2,i_month,i_sza,i_o3),
     $            o3_amftbl(i_alt,i_alb,i_wv,2,i_month,i_sza,i_o3) 
                 
						  ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	 
	    CLOSE(2)

	    IF(ind_lat.LT.17) THEN
			lat_vec_work(1)=lat_vec(ind_lat)
              lat_vec_work(2)=lat_vec(ind_lat+1)
	    ELSE
	        lat_vec_work(1)=lat_vec(17)
              lat_vec_work(2)=lat_vec(18)
	    ENDIF

C     Interpolation on the day number
      DO i_alt=1,n_alt
		DO i_alb=1,n_alb
			DO i_wv=1,n_wv
				DO i_lat=1,2
					DO i_sza=1,n_sza
						DO i_o3=1,n_o3
          CALL SPLINE(daynb_vec,o3_amftbl(i_alt,i_alb,i_wv,i_lat,:,
     $	     i_sza,i_o3),n_month,YWORK,Y2)
	    CALL SPLINT(daynb_vec,o3_amftbl(i_alt,i_alb,i_wv,i_lat,:,
     $		i_sza,i_o3),Y2,n_month,daynb_int,ycalc)
		o3_amftbl_dn(i_alt,i_alb,i_wv,i_lat,i_sza,i_o3)=ycalc
          CALL SPLINE(daynb_vec,o3_scdtbl(i_alt,i_alb,i_wv,i_lat,:,
     $	     i_sza,i_o3),n_month,YWORK,Y2)
	    CALL SPLINT(daynb_vec,o3_scdtbl(i_alt,i_alb,i_wv,i_lat,:,
     $		i_sza,i_o3),Y2,n_month,daynb_int,ycalc)
		o3_scdtbl_dn(i_alt,i_alb,i_wv,i_lat,i_sza,i_o3)=ycalc
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDDO
	 

c     Interpolation on the altitude
	DO i_alb=1,n_alb
		DO i_wv=1,n_wv
			DO i_lat=1,2
				DO i_sza=1,n_sza
					DO i_o3=1,n_o3
				CALL inter(n_alt,n_alt,2,alt_int,alt_vec,
     $	o3_amftbl_dn(:,i_alb,i_wv,i_lat,i_sza,i_o3),ycalc,hh)
		o3_amftbl_alt(i_alb,i_wv,i_lat,i_sza,i_o3)=ycalc
	            CALL inter(n_alt,n_alt,2,alt_int,alt_vec,
     $	o3_scdtbl_dn(:,i_alb,i_wv,i_lat,i_sza,i_o3),ycalc,hh)
		o3_scdtbl_alt(i_alb,i_wv,i_lat,i_sza,i_o3)=ycalc
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDDO


C     Interpolation on the albedo
	DO i_wv=1,n_wv
		DO i_lat=1,2
			DO i_sza=1,n_sza
				DO i_o3=1,n_o3
					CALL inter(n_alb,n_alb,2,alb_int,alb_vec,
     $	o3_amftbl_alt(:,i_wv,i_lat,i_sza,i_o3),ycalc,hh)
		o3_amftbl_alb(i_wv,i_lat,i_sza,i_o3)=ycalc
					CALL inter(n_alb,n_alb,2,alb_int,alb_vec,
     $	o3_scdtbl_alt(:,i_wv,i_lat,i_sza,i_o3),ycalc,hh)
		o3_scdtbl_alb(i_wv,i_lat,i_sza,i_o3)=ycalc
				ENDDO
			ENDDO
		ENDDO
	ENDDO
	

C     Interpolation on the wavelength
	DO i_lat=1,2
		DO i_sza=1,n_sza
			DO i_o3=1,n_o3
          CALL SPLINE(wv_vec,o3_amftbl_alb(:,i_lat,i_sza,i_o3),n_wv,
     $		YWORK,Y2)
	    CALL SPLINT(wv_vec,o3_amftbl_alb(:,i_lat,i_sza,i_o3),Y2,n_wv,
     $      			wvl_int,ycalc)
		o3_amftbl_wv(i_lat,i_sza,i_o3)=ycalc
          CALL SPLINE(wv_vec,o3_scdtbl_alb(:,i_lat,i_sza,i_o3),n_wv,
     $		YWORK,Y2)
	    CALL SPLINT(wv_vec,o3_scdtbl_alb(:,i_lat,i_sza,i_o3),Y2,n_wv,
     $      			wvl_int,ycalc)
		o3_scdtbl_wv(i_lat,i_sza,i_o3)=ycalc
			ENDDO
		ENDDO
	ENDDO
      
C     Interpolation on the latitude	
	DO i=1,n_o3
		IF (o3_scdtbl_wv(1,1,i).LT.1.d17) THEN
			ind_low(i)=1
	    ELSE
			ind_low(i)=0
	    ENDIF
	ENDDO


	DO i=1,n_o3
		IF (o3_scdtbl_wv(2,1,i).LT.1.d17) THEN
			ind_up(i)=1
		ELSE 
              ind_up(i)=0
	    ENDIF
	ENDDO

	DO i_sza=1,n_sza
		DO i_o3=1,n_o3
			CALL inter(2,2,2,lat_int,lat_vec_work,
     $	o3_amftbl_wv(:,i_sza,i_o3),ycalc,hh)
		o3_amftbl_lat(i_sza,i_o3)=ycalc
			CALL inter(2,2,2,lat_int,lat_vec_work,
     $	o3_scdtbl_wv(:,i_sza,i_o3),ycalc,hh)
		o3_scdtbl_lat(i_sza,i_o3)=ycalc

          IF(lat_int.GE.0) THEN
			IF(ABS(lat_int-lat_vec(2)).LE.5) THEN
				IF(ind_low(i_o3).EQ.1.AND.ind_up(i_o3).EQ.1) THEN
					o3_amftbl_lat(i_sza,i_o3)=999.d0
					o3_scdtbl_lat(i_sza,i_o3)=999.d0
				ELSEIF (ind_low(i_o3).EQ.1.AND.ind_up(i_o3).NE.1) THEN
	      o3_amftbl_lat(i_sza,i_o3)=o3_amftbl_wv(2,i_sza,i_o3)
		  o3_scdtbl_lat(i_sza,i_o3)=o3_scdtbl_wv(2,i_sza,i_o3)
				ENDIF
			ELSE
				IF(ind_low(i_o3).EQ.1) THEN
					o3_amftbl_lat(i_sza,i_o3)=999.d0
					o3_scdtbl_lat(i_sza,i_o3)=999.d0
				ENDIF
			ENDIF
		ELSE
			IF(ABS(lat_int-lat_vec(2)).LE.5) THEN
				IF(ind_up(i_o3).EQ.1) THEN
					o3_amftbl_lat(i_sza,i_o3)=999.d0
					o3_scdtbl_lat(i_sza,i_o3)=999.d0
				ENDIF
			ELSE	
				IF(ind_low(i_o3).EQ.1.AND.ind_up(i_o3).EQ.1) THEN
					o3_amftbl_lat(i_sza,i_o3)=999.d0
					o3_scdtbl_lat(i_sza,i_o3)=999.d0
				ELSEIF (ind_low(i_o3).NE.1.AND.ind_up(i_o3).EQ.1) THEN
	      o3_amftbl_lat(i_sza,i_o3)=o3_amftbl_wv(1,i_sza,i_o3)
		  o3_scdtbl_lat(i_sza,i_o3)=o3_scdtbl_wv(1,i_sza,i_o3)
				ENDIF
			ENDIF

	    ENDIF

		ENDDO
	ENDDO    

c     Interpolation on the SZA
	DO i_sza_int=1,n_sza_int
		DO i_o3=1,n_o3
	        CALL SPLINE(sza_vec,o3_amftbl_lat(:,i_o3),n_sza,YWORK,Y2)
			CALL SPLINT(sza_vec,o3_amftbl_lat(:,i_o3),Y2,n_sza,
     $      			sza_int(i_sza_int),ycalc)
	        o3_amftbl_sza_int(i_sza_int,i_o3)=ycalc
              CALL SPLINE(sza_vec,o3_scdtbl_lat(:,i_o3),n_sza,YWORK,Y2)
			CALL SPLINT(sza_vec,o3_scdtbl_lat(:,i_o3),Y2,n_sza,
     $      			sza_int(i_sza_int),ycalc)
			o3_scdtbl_sza_int(i_sza_int,i_o3)=ycalc
		ENDDO
	ENDDO

c     Interpolation on the O3 slant column
      IF(inter_flag.EQ.2) THEN
c Interpolation on O3 SCD in molec/cm2
		DO i_sza_int=1,n_sza_int
			cc=1
			DO i_o3=1,n_o3
				IF(o3_amftbl_sza_int(i_sza_int,i_o3).LT.997.d0) THEN
				o3_scd_effective(cc)=o3_scdtbl_sza_int(i_sza_int,i_o3)
	            o3_amf_effective(cc)=o3_amftbl_sza_int(i_sza_int,i_o3)
	            cc=cc+1
				ENDIF
			ENDDO

              CALL inter(cc-1,cc-1,2,scd_int(i_sza_int),
     $	    o3_scd_effective(1:cc-1),
     $			o3_amf_effective(1:cc-1),ycalc,hh)
	        o3_amf(i_sza_int)=ycalc
	        

			IF(scd_int(i_sza_int).GE.o3_scd_effective(1).AND.
     $        scd_int(i_sza_int).LE.o3_scd_effective(cc-1)) THEN 
	
	WRITE(*,27) '   SZA=',sza_int(i_sza_int),'deg: O3 AMF interpolated' 
	       
			ELSE
	write(*,28)	'   SZA=',sza_int(i_sza_int),'deg: O3 SCD OUT OF TOMS 
     $O3 SCD RANGE --> O3 AMF extrapolated' 
	        
			ENDIF
		ENDDO
	ELSE
c	Interpolation on O3 VCD in DU
		DO i_sza_int=1,n_sza_int
			cc=1
			DO i_o3=1,n_o3
                  IF(o3_amftbl_sza_int(i_sza_int,i_o3).LT.997.d0) THEN
				o3_vcd_effective(cc)=o3_vcd(i_o3)
	            o3_amf_effective(cc)=o3_amftbl_sza_int(i_sza_int,i_o3)
	            cc=cc+1
				ENDIF
			ENDDO

	       CALL inter(cc-1,cc-1,2,vcd_int(i_sza_int),
     $	   o3_vcd_effective(1:cc-1),o3_amf_effective(1:cc-1),ycalc,hh)
             o3_amf(i_sza_int)=ycalc
			IF(vcd_int(i_sza_int).GE.o3_vcd_effective(1).AND.
     $        vcd_int(i_sza_int).LE.o3_vcd_effective(cc-1)) THEN
	WRITE(*,27) '   SZA=',sza_int(i_sza_int),'deg: O3 AMF interpolated'
			
			ELSE
      write(*,28)	'   SZA=',sza_int(i_sza_int),'deg: O3 VCD OUT OF TOMS 
     $O3 VCD RANGE --> O3 AMF extrapolated' 
			ENDIF
		ENDDO
	ENDIF

27    FORMAT(a7,F5.2,a24)
28    FORMAT(a7,F5.2,a60)

C    Write AMF values in the output file
      OPEN(1,file='o3_amf_output.dat')
		WRITE(1,201) '% O3 AMF'
	    WRITE(1,202) '% Latitude (°)                    : ',
     $     lat_int_out
          WRITE(1,202) '% Longitude (°)                   : ',
     $     long_int_out
	    WRITE(1,203) '% Day number                      : ', 
     $     daynb_int_out
		WRITE(1,202) '% Wavelength (nm)                 : ', wvl_int
	    WRITE(1,204) '% Albedo                          : ', alb_int
	    WRITE(1,204) '% Altitude of the station (km)    : ', alt_int
	IF(inter_flag.EQ.2) THEN
	    WRITE(1,205) '% First column: SZA (°), second column: O3 SCD
     $ (molec/cm2), third column: O3 AMF'
	ELSE
          WRITE(1,208) '% First column: SZA (°), second column: O3 VCD
     $ (DU), third column: O3 AMF'
	ENDIF
	
      DO i_sza_int=1,n_sza_int
	IF(inter_flag.EQ.2) THEN
	    WRITE(1,206) sza_int(i_sza_int),scd_int(i_sza_int),
     $         o3_amf(i_sza_int)
	ELSE
          WRITE(1,207) sza_int(i_sza_int),vcd_int(i_sza_int),
     $         o3_amf(i_sza_int)
	ENDIF
	ENDDO
      
	CLOSE(1)

	

201   FORMAT(a8)
202   FORMAT(a36,F6.2)
203   FORMAT(a36,F5.1)
204   FORMAT(a36,F4.2)
205   FORMAT(a80)
206   FORMAT(F7.2,x,E10.4,x,F7.3)
207   FORMAT(F7.2,x,F7.1,x,F7.3)
208   FORMAT(a73)

	END
C  End of interplotate_O3_AMF
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
*     itype     Interpolation type (1: exponential; 2:linear)
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

	 IF ( arg .LE. xarr(npoints) .AND. arg .GE. xarr(1)) THEN
         DO 10 iq = 1 , npoints-1
           IF ( arg .LT. xarr(iq+1) .AND. arg .GE. xarr(iq)) ip = iq
   10    CONTINUE
         IF ( arg .EQ. xarr(npoints)) ip = npoints-1
       ELSEIF ( arg .GT. xarr(npoints)) THEN
         ip = npoints - 1
       ELSEIF ( arg .LT. xarr(1)) THEN
         ip = 1
       ENDIF

*
*     Interpolate function value at arg from data points ip to ip+1
*
*
*     exponential interpolation
*
       IF ( itype .EQ. 1 ) THEN
          IF ( yarr(ip+1) .EQ. yarr(ip) ) THEN
             hh   = 0.0
             ynew = yarr(ip)
          ELSE
             hh   = -( xarr(ip+1) - xarr(ip) ) / 
     $            LOG( yarr(ip+1) / yarr(ip))
             ynew = yarr(ip) * EXP(- ( arg - xarr(ip) ) / hh )
          ENDIF
*
*     linear interpolation
*
       ELSEIF ( itype .EQ. 2 ) THEN
         hh = ( yarr(ip+1) - yarr(ip) ) / ( xarr(ip+1) - xarr(ip) )
         ynew = yarr(ip) + hh*( arg - xarr(ip))
*     
      ENDIF
      RETURN
      END       
C******************************************************************************
      SUBROUTINE SPLINE(X,Y,N,WORK,Y2)
C******************************************************************************
C     Routine to calculate 2.nd derivatives of tabulated function
C     Y(i)=Y(Xi), to be used for cubic spline calculation.
C     From numerical recipes in fortran, Cambridge University Press, 1992
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
C     From numerical recipes in fortran, Cambridge University Press, 1992
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
