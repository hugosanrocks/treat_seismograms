!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +                                                                    +
!     + This program computes the first P and S arrival times in a given   +
!     + stations array from sources describing in a cartesian grid by mean +
!     + of finite difference approach developed by:                        +
!     +                                                                    +
!     + P.PODVIN, Geophysique, Ecole des Mines de Paris, Fontainebleau.    +
!     + e-mail: Pascal.Podvin@ensmp.fr    Tel: 33-(1) 64 69 49 25.         +
!     +                                                                    +
!     + Contact: Victor M. CRUZ-ATIENZA (cruz@geofisica.unam.mx)           +
!     +                                                                    +
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      PROGRAM arritime

      USE s_hs

      IMPLICIT NONE

      INTEGER(KIND=4) :: i,j,k
      INTEGER(KIND=4) :: ix,iy,iz
      INTEGER(KIND=4) :: nsta,nsrc,ista,isrc,nhypo
      INTEGER(KIND=4) :: ixsta,iysta
      INTEGER(KIND=4) :: msg,indx,nl0
      INTEGER(KIND=4) :: itemps,time_3d
      INTEGER(KIND=4) :: ilay
      ! Grids
      INTEGER(KIND=4) :: n1,n2,n3,ntot
      INTEGER(KIND=4) :: nx,ny,nz

      REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: stacoor,srcoor,arvsfs,arvpfs
      REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:,:) :: ts,tp,vitp,vits
      REAL(KIND=4), ALLOCATABLE, DIMENSION(:) :: TH,AL0,BE0
      ! Virtual Source
      REAL(KIND=4) :: xs,ys,zs
      REAL(KIND=4) :: hm
      REAL(KIND=4) :: eps
      REAL(KIND=4) :: vs,vp
      REAL(KIND=4) :: xmin,xmax,ymin,ymax
      REAL(KIND=4) :: xsize,ysize,zsize
      REAL(KIND=4) :: marge,hypdep
      REAL(KIND=4) :: ixd,iyd,izd
      REAL(KIND=4) :: trash,cdep,idep

      CHARACTER (LEN = 20) :: chypo
      CHARACTER (LEN = 30) :: vsfl,vpfl
      CHARACTER(LEN = 4), DIMENSION(:), ALLOCATABLE :: sname

      OPEN(10,FILE='arrtime.in',STATUS='old')
      OPEN(12,FILE='stations-utm.dat',STATUS='old')
      !OPEN(13,FILE='times.out',STATUS='unknown')
      !OPEN(14,FILE='slowness.out',STATUS='unknown')

      WRITE(*,*) 
      WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++'
      WRITE(*,*) '           PROGRAM: ARRTIME'
      WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++'
      WRITE(*,*) 

      ! Set constants
      msg=0
      eps=1e-3

      ! Read Input Parameters
      READ(10,*) nsta
      READ(10,*) hypdep
      READ(10,*) hm
      READ(10,*) marge
      READ(10,*) 
      READ(10,*) nl0
      ! Allocate properties of velocity model
      ALLOCATE ( TH(nl0),AL0(nl0),BE0(nl0) )
      DO i=1,nl0
        READ(10,*) TH(i),AL0(i),BE0(i)
      END DO

      ! Convert grid size to meters
      hm=hm*1e3
      marge=marge*1e3
      hypdep=hypdep*1e3
      TH(:)=TH(:)*1e3
      AL0(:)=AL0(:)*1e3
      BE0(:)=BE0(:)*1e3

      ! Allocate stations array
      ALLOCATE ( stacoor(nsta,2), sname(nsta) )

      ! Read Stations Coordinates (UTM, km)
      DO ista=1,nsta
        READ(12,FMT=104) sname(ista),stacoor(ista,1),stacoor(ista,2)
        stacoor(ista,1)=stacoor(ista,1)*1e3
        stacoor(ista,2)=stacoor(ista,2)*1e3
        !WRITE(*,*) sname(ista)
      END DO

      ! Determine the hypocenter grid size from stations' array (as done for NVTs location)
      xmin=MINVAL(stacoor(:,1))
      xmin=MIN(0.0,xmin)
      xmax=MAXVAL(stacoor(:,1))
      xmax=MAX(0.0,xmax)
      ymin=MINVAL(stacoor(:,2))
      ymin=MIN(0.0,ymin)
      ymax=MAXVAL(stacoor(:,2))
      ymax=MAX(0.0,ymax)

      xsize=xmax-xmin+2.0*marge
      ysize=ymax-ymin+2.0*marge
      zsize=hypdep+marge
      nsrc=1

      ! Allocate source array
      ALLOCATE ( srcoor(nsrc,3),arvsfs(nsta,nsrc), arvpfs(nsta,nsrc) ) 

      ! Hypocenter Location
      nhypo=1
      srcoor(nhypo,1)=0.0
      srcoor(nhypo,2)=0.0
      srcoor(nhypo,3)=hypdep

      ! Translate epicenter and stations to get positive coordinates
      srcoor(nhypo,1)=srcoor(nhypo,1)+ABS(xmin)+marge
      srcoor(nhypo,2)=srcoor(nhypo,2)+ABS(ymin)+marge
      stacoor(1:nsta,1)=stacoor(1:nsta,1)+ABS(xmin)+marge
      stacoor(1:nsta,2)=stacoor(1:nsta,2)+ABS(ymin)+marge
      !WRITE(*,*) srcoor(nhypo,1)/1e3,srcoor(nhypo,2)/1e3,srcoor(nhypo,3)/1e3

      ! Determine the size of the velocity model grid
      n1=NINT(xsize/hm)+1
      n2=NINT(ysize/hm)+1
      n3=NINT(zsize/hm)+1
      ntot=n1*n2*n3
      
      !WRITE(*,FMT=100) nx,ny,nz
      WRITE(*,FMT=101) n1,n2,n3
      WRITE(*,*) 

      ! Allocate arrival times and velocity models
      ALLOCATE ( ts(n1,n2,n3),tp(n1,n2,n3),vits(n1,n2,n3),vitp(n1,n2,n3) )

      ! Generate the Velocity Model Grid
      ! Homogeneous Model
!      vs=3300.0 ! Vs (m/s)
!      vp=5716.0 ! Vp (m/s)
!      indx=0          
!      DO i=1,n1
!        DO j=1,n2 
!          DO k=1,n3      
!             indx=indx+1
!             !v(indx)=vit(i,j,k)
!             vits(i,j,k)=vs
!             vitp(i,j,k)=vp
!          END DO               
!        END DO                
!      END DO  

      ! Heterogeneous Layered Model
       ilay=1
       cdep=0.0          
       idep=TH(ilay)
       vs=BE0(ilay)
       vp=AL0(ilay)
       DO k=1,n3      
         cdep=(k-1)*hm
         IF (cdep .GT. idep) THEN
           ilay=ilay+1
           idep=idep+TH(ilay)
           vs=BE0(ilay)
           vp=AL0(ilay)
         END IF
         DO j=1,n2 
           DO i=1,n1
              vits(i,j,k)=vs
              vitp(i,j,k)=vp
           END DO               
         END DO                
       END DO  

      ! Compute slowness models
      CALL sub_hs(vits,n1,n2,n3,hm)
      CALL sub_hs(vitp,n1,n2,n3,hm)

      ! Compute first P and S arrival times
      DO isrc=1,nsrc
      !DO isrc=1,10
      !DO isrc=1,1
        WRITE(*,FMT=102) nsrc-(isrc-1)
        ! Source coordinates
        xs=srcoor(isrc,1)/hm+1
        ys=srcoor(isrc,2)/hm+1
        zs=srcoor(isrc,3)/hm+1
        ! Compute arrival times
        itemps=time_3d(vits,ts,n1,n2,n3,xs,ys,zs,eps,msg)
        itemps=time_3d(vitp,tp,n1,n2,n3,xs,ys,zs,eps,msg)
        DO ista=1,nsta
          ixsta=NINT(stacoor(ista,1)/hm)+1
          iysta=NINT(stacoor(ista,2)/hm)+1
          ! Extract stations' arrival time
          arvsfs(ista,isrc) = ts(ixsta,iysta,1)
          arvpfs(ista,isrc) = tp(ixsta,iysta,1)
        END DO
      END DO

      WRITE(*,*) 
      WRITE(*,*) ' Writing Arrival Times...'
      DO ista=1,nsta
        OPEN(42,FILE='output.dat',STATUS='unknown',POSITION='append')
        DO isrc=1,nsrc
          WRITE(42,FMT=104) sname(ista),arvpfs(ista,isrc),arvsfs(ista,isrc)
        END DO
        CLOSE(40)
        CLOSE(41)
      END DO

      WRITE(*,*) 
      WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++'
      WRITE(*,*) 

!100   FORMAT('  Sources grid: nx=',I4,2X,'ny=',I4,2X,'nz=',I4)
101   FORMAT('    Model grid: n1=',I4,2X,'n2=',I4,2X,'n3=',I4)
102   FORMAT('    Remaining Sources:',I4)
103   FORMAT('    Writing station: ',A4)
104   FORMAT(A4,F12.5,F12.5)
105   FORMAT(2F12.5)

      ! To check
      !DO i=1,n1
      !  DO j=1,n2 
      !    WRITE(13,*) ts(i,j,1)
      !    DO k=1,n3      
      !      indx=indx+1
      !      !v(indx)=vits(i,j,k)
      !      !WRITE(13,*) ts(i,j,k)
      !    END DO               
      !  END DO                
      !END DO  

      END PROGRAM
