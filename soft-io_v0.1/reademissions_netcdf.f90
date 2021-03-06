subroutine reademissions_netcdf(nemission,outlon0,outlat0,&
  numxgrid,numygrid,dxout,dyout,nxmax,nymax,nzmax,area,heightnn,&
  emissions,index_cont,nbtime,date_flex)
  
  use netcdf

  implicit none

  INCLUDE 'netcdf.inc'

  integer nemission,numxgrid,numygrid,nxmax,nymax,nzmax,nbtime
  real outlon0,outlat0,dxout,dyout

  integer m,i,j
  real area(nxmax,nymax)
  real heightnn(nxmax,nymax,0:nzmax)

  real emissions(nxmax,nymax,nbtime)
  double precision date_flex(nbtime),juldate,julian_nc,julian_nc0
  integer index_cont(nxmax,nymax),indtime(nbtime)
  real pi,r_earth,pih,xl,xr,xm,yl,ym,yr
  integer maxemissions,nxem,nyem,maxmegaemissions,it
  integer maxpoint,indx,indy,ix,jy,n,neur,yeardate,timeday,year

  integer yeardate1,yeardate2,date_emis_beg,date_emis_end,date_emis_nc

  parameter(pi=3.14159265,r_earth=6.371e6,pih=pi/180.)
  parameter(maxemissions=4,nxem=360,nyem=180)
  parameter(maxmegaemissions=4)
  real emegcell(12),maxemis,xdiff,ydiff
  real emission_tami(360,180)
  integer indcont(nxem,nyem)
  parameter(maxpoint=1000000)
  real xlonl(maxpoint),xlonr(maxpoint),ylatl(maxpoint)
  real ylatr(maxpoint),zl(maxpoint),zr(maxpoint)

  character*256 maccfile,emislist,contlist,contfile,pixelareafile

  real weightmolar(0:maxemissions)
  data weightmolar/28.97,28.,46.,64.,12.01/

  integer ufile_nc,lon_dimid,lat_dimid,time_dimid,time_varid
  integer lon_varid,lat_varid,emis_varid,area_varid
  integer nblon_nc,nblat_nc,nbtime_nc

  integer ufile_cont,cont_varid

  real,allocatable :: lon_nc(:),lat_nc(:),emis_nc(:,:,:)
  real,allocatable :: cont_nc(:,:,:)
  real,allocatable :: area_nc(:,:)
  integer,allocatable :: time_nc(:)

  character*150 charr
  integer date_emis_anthro
  character*256 lab_anthropiques

  character*250 pref
  pref='         [reademis_netcdf]:'


!c Read emissions inventory
!**************************
  
  open(10,file='flexpart_contribution_levels.nml')
  read(10,'(a)') lab_anthropiques
  read(10,'(a)') 
  read(10,'(a)') maccfile
  read(10,'(a)') 
  read(10,'(a)') contfile
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') charr
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') pixelareafile
  close(10)
  read(charr,'(i8)')date_emis_anthro

  write(*,*) trim(pref),'flexpart_contribution_levels.nml:'
  write(*,*) trim(pref),' | lab_anthropiques:',trim(lab_anthropiques)
  write(*,*) trim(pref),' | maccfile        :',trim(maccfile)
  write(*,*) trim(pref),' | contfile        :',trim(contfile)
  write(*,*) trim(pref),' | date_emis_anthro:',date_emis_anthro


  write(*,*) 'Get emis data...'
  call check(nf_open(trim(maccfile),NF_NOWRITE,ufile_nc))

  call check(nf_inq_dimid(ufile_nc,'lon',lon_dimid))
  call check(nf_inq_dimlen(ufile_nc,lon_dimid,nblon_nc))
  call check(nf_inq_dimid(ufile_nc,'lat',lat_dimid))
  call check(nf_inq_dimlen(ufile_nc,lat_dimid,nblat_nc))
  call check(nf_inq_dimid(ufile_nc,'date',time_dimid)) ! should be only one timestep
  call check(nf_inq_dimlen(ufile_nc,time_dimid,nbtime_nc))

  allocate(lon_nc(nblon_nc),lat_nc(nblat_nc))
  allocate(time_nc(nbtime_nc))
  allocate(emis_nc(nblon_nc,nblat_nc,nbtime_nc)) ! kg/m2/s

  call check(nf_inq_varid(ufile_nc,'lon',lon_varid))
  call check(nf_get_var_real(ufile_nc,lon_varid,lon_nc))
  call check(nf_inq_varid(ufile_nc,'lat',lat_varid))
  call check(nf_get_var_real(ufile_nc,lat_varid,lat_nc))
  call check(nf_inq_varid(ufile_nc,'date',time_varid))
  call check(nf_get_var_int(ufile_nc,time_varid,time_nc))
  
  if(trim(lab_anthropiques).eq.'macc')then
     call check(nf_inq_varid(ufile_nc,'MACCity',emis_varid))
  elseif(trim(lab_anthropiques).eq.'edgar')then
     call check(nf_inq_varid(ufile_nc,'EDGARv4',emis_varid))
  elseif(trim(lab_anthropiques).eq.'reas')then
     call check(nf_inq_varid(ufile_nc,'REAS',emis_varid))
  endif
  call check(nf_get_var_real(ufile_nc,emis_varid,emis_nc))
  
  call check(nf_close(ufile_nc))
  
  !Pixel area
  write(*,*) 'Get pixel area...'
  call check(nf_open(trim(pixelareafile),NF_NOWRITE,ufile_nc))
  allocate(area_nc(nblon_nc,nblat_nc)) ! km2
  call check(nf_inq_varid(ufile_nc,'Pixel_area',area_varid))
  call check(nf_get_var_real(ufile_nc,area_varid,area_nc))
  call check(nf_close(ufile_nc))

!C Give indices to the different continents
!******************************************
!c  Water          0
!c  N America      1
!c  S America      2
!c  Europe         3
!c  Africa         4
!c  Asia           5
!c  Oceania        6

  write(*,*) trim(pref),trim(contfile)
  call check(nf_open(trim(contfile),NF_NOWRITE,ufile_cont))
  allocate(cont_nc(nblon_nc,nblat_nc,1)) !should have same dimensions as in emissions file
  call check(nf_inq_varid(ufile_cont,'GFED_regions_0_5deg',cont_varid))
  call check(nf_get_var_real(ufile_cont,cont_varid,cont_nc))
  call check(nf_close(ufile_cont))
  write(*,*) trim(pref),' >get',maxval(cont_nc),'regions codes in file'

  write(*,*) 'Check lon range...'
  do ix=1,nblon_nc
     if (lon_nc(ix).gt.(outlon0+numxgrid*dxout)) then
        lon_nc(ix)=lon_nc(ix)-360.
     endif
     if (lon_nc(ix).lt.outlon0) then
        lon_nc(ix)=lon_nc(ix)+360.
     endif
  end do

  call caldate(date_flex(1),yeardate1,timeday)      
  call caldate(date_flex(nbtime),yeardate2,timeday)
  julian_nc0=juldate(date_emis_anthro,0)
  call caldate(julian_nc0+time_nc(1),date_emis_beg,timeday)
  call caldate(julian_nc0+time_nc(nbtime_nc),date_emis_end,timeday)
  write(*,*) trim(pref),' >FLEXPART dates range from',yeardate1,'to',yeardate2
  write(*,*) trim(pref),' >EMIS dates range from',date_emis_beg,'to',date_emis_end
  write(*,*) trim(pref),date_emis_anthro,'from +',time_nc(1),'to +',time_nc(nbtime_nc),'days'

  call caldate(date_flex(nbtime),yeardate,timeday)      ! Now we search which netCDF time step
!  year=yeardate/10000                                   ! is to be used for every flexpart time step.
  do it=1,nbtime
     xdiff=999999.
     indx=1
     do n=1,nbtime_nc        
        julian_nc=juldate(date_emis_anthro,0)+time_nc(n)        
        if (abs(julian_nc-date_flex(it)).le.xdiff.and.julian_nc.le.date_flex(it)) then     
           xdiff=abs(julian_nc-date_flex(it))
           indx=n
        endif
     end do
     call caldate(date_flex(it),yeardate,timeday)
     call caldate(julian_nc0+time_nc(indx),date_emis_nc,timeday)
     write(*,*) trim(pref),'>>',it,yeardate,date_emis_nc,indx
     indtime(it)=indx                                    ! for each it, take time number indtime(it) in the netcdf
  end do

  write(*,*) trim(pref),'End reading anthropogenic emission data'
!  write(*,*) trim(pref),'CHECK emis index :emis_nc(363,82,1)          = ',emis_nc(363,82,1)
!  write(*,*) trim(pref),'CHECK emis index :emis_nc(363,82,indtime(1)) = ',emis_nc(363,82,indtime(1))





!*********************************************
!C Now attribute emissions to model grid cells
!*********************************************

  maxemis=0

  do ix=1,numxgrid ! start with empty emission grid
     do jy=1,numygrid
        do it=1,nbtime
           emissions(ix,jy,it)=0.
        end do
        index_cont(ix,jy)=0
     enddo !ix
  enddo !jy
  do ix=1,nblon_nc
     xdiff=360.
     indx=1
     do n=1,numxgrid ! find position of lon_nc(ix) in the flexpart grid
        xl=outlon0+float(n-1)*dxout
        xr=outlon0+float(n)*dxout
        xm=(xl+xr)/2.
        if (abs(lon_nc(ix)-xm).le.xdiff) then
           xdiff=abs(lon_nc(ix)-xm)
           indx=n
        endif
     enddo !n

     do jy=1,nblat_nc
        ydiff=360.
        indy=1
        do n=1,numygrid ! same for lat_nx(jy)
           yl=outlat0+float(n-1)*dyout
           yr=outlat0+float(n)*dyout
           ym=(yl+yr)/2.
           if (abs(lat_nc(jy)-ym).le.ydiff) then
              ydiff=abs(lat_nc(jy)-ym)
              indy=n
           endif
        enddo !n
        
        !c Now, we know that emis_nc(ix,jy) goes inside emissions(indx,indy)
        do it=1,nbtime
           emissions(indx,indy,it)=emissions(indx,indy,it)+&
                emis_nc(ix,jy,indtime(it))*&
                area_nc(ix,jy)*1e6/area(indx,indy)
        enddo
        index_cont(indx,indy)=cont_nc(ix,jy,1)
     enddo
  enddo

!NB:emis_nc(151,101,108)=emis(150,100,108) dans le netcdf


  !c Convert to kg/m2/s (injection height is applied in the main program)
  !C Account for molar weight to convert to volume mixing ratio
  maxemis=0
  do ix=1,numxgrid
     do jy=1,numygrid
        do it=1,nbtime
           emissions(ix,jy,it)=emissions(ix,jy,it)*weightmolar(0)/&
                weightmolar(nemission)
           if (emissions(ix,jy,it).gt.maxemis)maxemis=emissions(ix,jy,it)
        enddo
     enddo
  enddo
  write(*,*) trim(pref),'maxemis anthro',maxemis
  return
end subroutine reademissions_netcdf
