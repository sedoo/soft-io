subroutine reademissions_feux(nemission,outlon0,outlat0,&
     numxgrid,numygrid,dxout,dyout,nxmax,nymax,nzmax,area,heightnn,&
     emifeux,index_feux,nbtime,date_flex)

  use netcdf

  implicit none

  INCLUDE 'netcdf.inc'

  integer nemission,numxgrid,numygrid,nxmax,nymax,nzmax,nbtime
  real outlon0,outlat0,dxout,dyout

  integer m,i,j
  real area(nxmax,nymax)
  real heightnn(nxmax,nymax,0:nzmax)

  real emifeux(nxmax,nymax,nbtime)
  integer index_feux(nxmax,nymax),indtime(nbtime)
  double precision date_flex(nbtime),juldate,julian_nc,julian_nc0
  real pi,r_earth,pih,xl,xr,xm,yl,ym,yr
  integer maxemissions,maxmegaemissions,it
  !integer nxem,nyem
  integer maxpoint,indx,indy,ix,jy,n,neur,yeardate,yeardate2,timeday,year
  parameter(pi=3.14159265,r_earth=6.371e6,pih=pi/180.)
  parameter(maxemissions=4)
  !parameter(nxem=360,nyem=180)
  parameter(maxmegaemissions=4)
  real emegcell(12),maxemis,xdiff,ydiff
  !real emission_tami(360,180)
  !integer indcont(nxem,nyem)
  parameter(maxpoint=1000000)
  real xlonl(maxpoint),xlonr(maxpoint),ylatl(maxpoint)
  real ylatr(maxpoint),zl(maxpoint),zr(maxpoint)

  character*256 bbfile,emislist,contlist,contfile,pixelareafile
  character*32 lonname,latname,datename,coname

  real weightmolar(0:maxemissions)
  data weightmolar/28.97,28.,46.,64.,12.01/

  integer ufile_nc,lon_dimid,lat_dimid,time_dimid
  integer lon_varid,lat_varid,emis_varid,area_varid
  integer nblon_nc,nblat_nc,nbtime_nc
  integer nblon_cont,nblat_cont,icont,jcont
  integer nblon_pix,nblat_pix

  integer ufile_cont,cont_varid,time_varid

  real,allocatable :: lon_nc(:),lat_nc(:)
  integer,allocatable :: indlon_nc(:),indlat_nc(:)
  real,allocatable :: emis_nc(:,:,:)
  real,allocatable :: lon_cont(:),lat_cont(:)
  integer,allocatable :: indlon_cont(:),indlat_cont(:)
  integer,allocatable :: cont_nc(:,:,:)
  real,allocatable :: area_nc(:,:)
  integer,allocatable :: time_nc(:)

  character*150 charr
  integer date_emis_feux,date_emis_beg,date_emis_end,date_emis_nc
  character*256 lab_feuxbiomasse

  real meanemis
  character*250 pref
  pref='         [reademis_feux]:'

!c Read emissions inventory
!**************************

  open(10,file='flexpart_contribution_levels.nml')
  read(10,'(a)') 
  read(10,'(a)') lab_feuxbiomasse
  read(10,'(a)') 
  read(10,'(a)') bbfile
  read(10,'(a)') 
  read(10,'(a)') contfile
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') charr
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') pixelareafile
  close(10)
  read(charr,'(i8)')date_emis_feux

  write(*,*) trim(pref),'flexpart_contribution_levels.nml:'
  write(*,*) trim(pref),' >lab_feuxbiomasse: ',trim(lab_feuxbiomasse)
  write(*,*) trim(pref),' >bbfile          : ',trim(bbfile)
  write(*,*) trim(pref),' >contfile        : ',trim(contfile)
  write(*,*) trim(pref),' >date_emis_feux  : ',date_emis_feux   

  if(trim(lab_feuxbiomasse).eq.'gfas')then
     lonname='lon'
     latname='lat'
     datename='date'
     coname='GFAS'
  elseif(trim(lab_feuxbiomasse).eq.'gfas1.2_juelich')then
     lonname='lon'
     latname='lat'
     datename='date'
     coname='cofire'
  elseif(trim(lab_feuxbiomasse).eq.'gfas1.2')then
     lonname='lon'
     latname='lat'
     datename='date'
     coname='cofire'
  elseif(trim(lab_feuxbiomasse).eq.'gfas0.1')then
     lonname='longitude'
     latname='latitude'
     datename='time'
     coname='cofire'
  elseif(trim(lab_feuxbiomasse).eq.'gfed4')then
     lonname='longitude'
     latname='latitude'
     datename='Time'
     coname='cofire'
  elseif(trim(lab_feuxbiomasse).eq.'gfed3')then
     lonname='lon'
     latname='lat'
     datename='date'
     coname='GFED'
  endif 


  call check(nf_open(trim(bbfile),NF_NOWRITE,ufile_nc))

  call check(nf_inq_dimid(ufile_nc,trim(lonname),lon_dimid))
  call check(nf_inq_dimlen(ufile_nc,lon_dimid,nblon_nc))
!  write(*,*) trim(pref),' >check lon : dim =',nblon_nc

  call check(nf_inq_dimid(ufile_nc,trim(latname),lat_dimid))
  call check(nf_inq_dimlen(ufile_nc,lat_dimid,nblat_nc))
!  write(*,*) trim(pref),' >check lat : dim =',nblat_nc

  call check(nf_inq_dimid(ufile_nc,trim(datename),time_dimid)) ! should be only one timestep
  call check(nf_inq_dimlen(ufile_nc,time_dimid,nbtime_nc))
!  write(*,*) trim(pref),' >check date : dim =',nbtime_nc
  write(*,*) trim(pref),' >get dim :',nblon_nc,'lon',nblat_nc,'lat',nbtime_nc,'dates'

  allocate(lon_nc(nblon_nc),lat_nc(nblat_nc))
  allocate(indlon_nc(nblon_nc),indlat_nc(nblat_nc))
  allocate(indlon_cont(nblon_nc),indlat_cont(nblat_nc))
  allocate(time_nc(nbtime_nc))
  allocate(emis_nc(nblon_nc,nblat_nc,nbtime_nc)) ! kg/m2/s

  call check(nf_inq_varid(ufile_nc,trim(lonname),lon_varid))
  call check(nf_get_var_real(ufile_nc,lon_varid,lon_nc))
  call check(nf_inq_varid(ufile_nc,trim(latname),lat_varid))
  call check(nf_get_var_real(ufile_nc,lat_varid,lat_nc))
  call check(nf_inq_varid(ufile_nc,trim(datename),time_varid))
  call check(nf_get_var_int(ufile_nc,time_varid,time_nc))
  !write(*,*) trim(pref),' >get var : OK'

  call check(nf_inq_varid(ufile_nc,trim(coname),emis_varid))
  call check(nf_get_var_real(ufile_nc,emis_varid,emis_nc))
  call check(nf_close(ufile_nc))

  !Replace "epsilon" value to 0 for GFAS 1.2
  if(trim(lab_feuxbiomasse).eq.'gfas1.2')then
     do ix=1,nblon_nc
        do jy=1,nblat_nc
           do n=1,nbtime_nc              
              if(emis_nc(ix,jy,n).lt.1e-20)emis_nc(ix,jy,n)=0.
           enddo
        enddo
     enddo    
  endif
  !Replace missing or fill values to 0 for GFAS 1.2 @ 0.1°
  if(trim(lab_feuxbiomasse).eq.'gfas0.1')then
     !emis_nc=8.60349708996174e-11*emis_nc+2.81902185649686e-06
     do ix=1,nblon_nc
        do jy=1,nblat_nc
           do n=1,nbtime_nc              
              if(emis_nc(ix,jy,n).lt.1e-20)emis_nc(ix,jy,n)=0.
           enddo
        enddo
     enddo        
  endif

  write(*,*) trim(pref),' >get emissions : OK'
  write(*,*) trim(pref),' >range:',minval(emis_nc),maxval(emis_nc)

  !Inverse the latitude index of EMIS.nc to REGIONS.nc to fit -- ????? 

  !Pixel area
  write(*,*) trim(pref),' >get pixel area from:',trim(pixelareafile)
  call check(nf_open(trim(pixelareafile),NF_NOWRITE,ufile_nc))
  call check(nf_inq_dimid(ufile_nc,'lon',lon_dimid))
  call check(nf_inq_dimlen(ufile_nc,lon_dimid,nblon_pix))
  call check(nf_inq_dimid(ufile_nc,'lat',lat_dimid))
  call check(nf_inq_dimlen(ufile_nc,lat_dimid,nblat_pix))
  write(*,*) trim(pref),' >find',nblon_pix,'lon',nblat_pix,'lat'
!TODO check pixel area dim
  allocate(area_nc(nblon_nc,nblat_nc)) ! km2
  call check(nf_inq_varid(ufile_nc,'Pixel_area',area_varid))
  call check(nf_get_var_real(ufile_nc,area_varid,area_nc))
  call check(nf_close(ufile_nc))
  write(*,*) trim(pref),' >get pixel area: OK'

!C Give indices to the different continents
!******************************************
  write(*,*) trim(pref),' >Load Region codes from:',trim(contfile)
  call check(nf_open(trim(contfile),NF_NOWRITE,ufile_cont))
!TODO update region methodo
  !if (trim(contfile).eq.'GFED_fire_regions.nc')  
   call check(nf_inq_dimid(ufile_cont,'lon',lon_dimid))
   call check(nf_inq_dimlen(ufile_cont,lon_dimid,nblon_cont))
   call check(nf_inq_dimid(ufile_cont,'lat',lat_dimid))
   call check(nf_inq_dimlen(ufile_cont,lat_dimid,nblat_cont))
   write(*,*) trim(pref),' >find',nblon_cont,'longitudes &',nblat_cont,'latitudes'
   allocate(lon_cont(nblon_cont),lat_cont(nblat_cont))
   allocate(cont_nc(nblon_cont,nblat_cont,1))
   !get longitude and latitude arrays
   call check(nf_inq_varid(ufile_cont,'lon',lon_varid))
   call check(nf_get_var_real(ufile_cont,lon_varid,lon_cont))
   call check(nf_inq_varid(ufile_cont,'lat',lat_varid))
   call check(nf_get_var_real(ufile_cont,lat_varid,lat_cont))
   !get regions codes
   call check(nf_inq_varid(ufile_cont,'GFED_regions_0_5deg',cont_varid))
   call check(nf_get_var_int(ufile_cont,cont_varid,cont_nc))
  !endif
   call check(nf_close(ufile_cont))
   write(*,*) trim(pref),' >get',maxval(cont_nc),'regions codes in file'

  !Adapt the lon of EMIS.nc to the lon range of FLEXPART output
  write(*,*) trim(pref),' >adapt and indexing longitude range of emission to FLEXPART output'
  write(*,*) trim(pref),lon_nc(1),lon_nc(2),lon_nc(int(nblon_nc/2)),lon_nc(nblon_nc-1),lon_nc(nblon_nc)
  do ix=1,nblon_nc
   if (lon_nc(ix).gt.(outlon0+numxgrid*dxout)) then
    lon_nc(ix)=lon_nc(ix)-360.000000000
   endif
   if (lon_nc(ix).lt.outlon0) then
    lon_nc(ix)=lon_nc(ix)+360.000000000
   endif
  end do
  write(*,*) trim(pref),lon_nc(1),lon_nc(2),lon_nc(int(nblon_nc/2)),lon_nc(nblon_nc-1),lon_nc(nblon_nc)

  do ix=1,nblon_nc
   indx=1
   xl=outlon0
   xr=outlon0 + dxout
   do while ((lon_nc(ix).lt.xl).or.(lon_nc(ix).ge.xr).and.(indx.le.numxgrid))
    indx=indx+1
    xl=outlon0+float(indx-1)*dxout
    xr=outlon0+float(indx)*dxout
   enddo
   indlon_nc(ix)=indx
   !write(*,*) ix,lon_nc(ix),indx,xl,xr
  enddo
  write(*,*) trim(pref),indlon_nc(1),indlon_nc(2),indlon_nc(int(nblon_nc/2)),indlon_nc(nblon_nc-1),indlon_nc(nblon_nc)

  !Indexing the lat of EMIS.nc to the lat range of FLEXPART output
  write(*,*) trim(pref),' >indexing latitude range of emission to FLEXPART output'
  do jy=1,nblat_nc
   indy=1
   yl=outlat0
   yr=outlat0 + dxout
   do while ((lat_nc(jy).lt.yl).or.(lat_nc(jy).ge.yr).and.(indy.le.numygrid))
    indy=indy+1
    yl=outlat0+float(indy-1)*dyout
    yr=outlat0+float(indy)*dyout
   enddo
   indlat_nc(jy)=indy
  enddo
  write(*,*) trim(pref),lat_nc(1),lat_nc(2),lat_nc(int(nblat_nc/2)),lat_nc(nblat_nc-1),lat_nc(nblat_nc)
  write(*,*) trim(pref),indlat_nc(1),indlat_nc(2),indlat_nc(int(nblat_nc/2)),indlat_nc(nblat_nc-1),indlat_nc(nblat_nc)

  !indexing region to EMIS.nc dim
  write(*,*) trim(pref),' >indexing longitude range of regions codes to EMIS.nc'
  !write(*,*) trim(pref),lon_cont(1),lon_cont(2),lon_cont(int(nblon_cont/2)),lon_cont(nblon_cont-1),lon_cont(nblon_cont)
  do ix=1,nblon_nc
   indx=1
   xl=lon_cont(indx)
   xr=lon_cont(indx+1)
   do while (((lon_nc(ix).lt.xl).or.(lon_nc(ix).ge.xr)).and.(indx.lt.(nblon_cont-1)))
    indx=indx+1
    xl=lon_cont(indx)
    xr=lon_cont(indx+1)
   enddo
   indlon_cont(ix)=indx
   !write(*,*) ix,lon_nc(ix),indlon_cont(ix),lon_cont(indx),lon_cont(indx+1)
  enddo
  !write(*,*) trim(pref),indlon_cont(1),indlon_cont(2),indlon_cont(int(nblon_nc/2)),indlon_cont(nblon_nc-1),indlon_cont(nblon_nc)

  write(*,*) trim(pref),' >indexing latitude range of regions codes to EMIS.nc'
  write(*,*) trim(pref),lat_cont(1),lat_cont(2),lat_cont(int(nblat_cont/2)),lat_cont(nblat_cont-1),lat_cont(nblat_cont)
  do jy=1,nblat_nc
   indx=1
   xl=lat_cont(indx)-0.25
   xr=lat_cont(indx)+0.25
   do while (((lat_nc(jy).le.xl).or.(lat_nc(jy).gt.xr)).and.(indx.lt.(nblat_cont)))
    indx=indx+1
    xl=lat_cont(indx)-0.25
    xr=lat_cont(indx)+0.25
   enddo
   indlat_cont(jy)=indx
   !write(*,*) jy,lat_nc(jy),indlat_cont(jy),xr,xl
  enddo
  !write(*,*) trim(pref),indlon_cont(1),indlon_cont(2),indlon_cont(int(nblon_nc/2)),indlon_cont(nblon_nc-1),indlon_cont(nblon_nc)

  !CHECK date in GFAS 1.2 @ 0.1°
  if(trim(lab_feuxbiomasse).eq.'gfas1.0'.or.trim(lab_feuxbiomasse).eq.'gfas1.2')then
     do n=1,nbtime_nc
        time_nc(n)=n-1
     enddo
  endif

!  write(*,*) trim(pref),'1:',date_flex(1)
!  write(*,*) trim(pref),'2:',date_flex(2)
!  write(*,*) trim(pref),'3:',date_flex(3)
!  write(*,*) trim(pref),'4:',date_flex(4)
!  write(*,*) trim(pref),'nbtime:',date_flex(nbtime)

  call caldate(date_flex(1),yeardate,timeday)      
  call caldate(date_flex(nbtime),yeardate2,timeday)
  julian_nc0=juldate(date_emis_feux,0)
  call caldate(julian_nc0+time_nc(1),date_emis_beg,timeday)
  call caldate(julian_nc0+time_nc(nbtime_nc),date_emis_end,timeday)
  write(*,*) trim(pref),' >FLEXPART dates range from',yeardate,'to',yeardate2
  write(*,*) trim(pref),' >EMIS dates range from',date_emis_beg,'to',date_emis_end
  write(*,*) trim(pref),date_emis_feux,'from +',time_nc(1),'to +',time_nc(nbtime_nc),'days'

  call caldate(date_flex(nbtime),yeardate,timeday)      ! Now we search which netCDF time step
  do it=1,nbtime
     xdiff=999999.
     indx=1
     do n=1,nbtime_nc
        julian_nc=julian_nc0+time_nc(n)        
        if (abs(julian_nc-date_flex(it)).le.xdiff.and.julian_nc.le.date_flex(it)) then              
           xdiff=abs(julian_nc-date_flex(it))
           indx=n
        endif
     end do
     call caldate(date_flex(it),yeardate,timeday)
     call caldate(julian_nc0+time_nc(indx),date_emis_nc,timeday)
     indtime(it)=indx ! for each it, take time number indtime(it) in the netcdf  
     write(*,*) trim(pref),'>>',it,nbtime,yeardate,date_emis_nc,indtime(it),sum(emis_nc(:,:,indx))
  end do

  write(*,*) trim(pref),' >End reading fire emission data'


  !List all points correponding to islands with non assigned region
!  write(*,*)'----------------------------- Islands with non assigned region:'
!  do ix=1,nblon_nc
!     do jy=1,nblat_nc
!        if (cont_nc(ix,jy,1).eq.0.and.maxval(emis_nc(ix,jy,:)).ne.0)then
!           write(*,*) ix-1,'(+1)',jy-1,'(+1) Island (no assign region) : emis_nc(1)=',maxval(emis_nc(ix,jy,:))
!        endif
!     enddo
!  enddo

!*********************************************
!C Now attribute emissions to model grid cells
!*********************************************
  write(*,*) trim(pref),' >collect emission data and fire region index...'
  !initialize emission and fire index
  emifeux=0.
  index_feux=0

!old  do ix=1,numxgrid ! start with empty emission grid
!old     do jy=1,numygrid
!old        do it=1,nbtime
!old           emifeux(ix,jy,1)=0.
!old        end do
!old        index_feux(ix,jy)=0
!old     enddo !ix
!old  enddo !jy

  do ix=1,nblon_nc
!old     xdiff=360.
!old     indx=1
!old     do n=1,numxgrid ! find position of lon_nc(ix) in the flexpart grid
!old        xl=outlon0+float(n-1)*dxout
!old        xr=outlon0+float(n)*dxout
!old        xm=(xl+xr)/2.
!old        if (abs(lon_nc(ix)-xm).le.xdiff) then
!old           xdiff=abs(lon_nc(ix)-xm)
!old           indx=n
!old        endif
!old     enddo !n

     indx=indlon_nc(ix)
     icont=indlon_cont(ix)
     xl=outlon0+float(indx-1)*dxout
     xr=outlon0+float(indx)*dxout

     do jy=1,nblat_nc
!old        ydiff=360.
!old        indy=1
!old        do n=1,numygrid ! same for lat_nx(jy)
!old           yl=outlat0+float(n-1)*dyout
!old           yr=outlat0+float(n)*dyout
!old           ym=(yl+yr)/2.
!old       if (abs(lat_nc(jy)-ym).le.ydiff) then
!old              ydiff=abs(lat_nc(jy)-ym)
!old              indy=n
!old           endif
!old        enddo !n
        indy=indlat_nc(jy)
        jcont=indlat_cont(jy)
        yl=outlat0+float(indy-1)*dyout
        yr=outlat0+float(indy)*dyout

        
!c Now, we know that emis_nc(ix,jy) goes inside emissions(indx,indy)
        do it=1,nbtime
           emifeux(indx,indy,it)=emifeux(indx,indy,it)+&
                emis_nc(ix,jy,indtime(it))*&
                area_nc(ix,jy)*1e6/area(indx,indy)
        enddo
        if(cont_nc(icont,jcont,1).ne.0)then
           index_feux(indx,indy)=cont_nc(icont,jcont,1)
           !write(*,*) ix,lon_nc(ix),xl,xr,lon_cont(icont),&
           !jy,lat_nc(jy),yl,yr,lat_cont(jcont),cont_nc(icont,jcont,1)
        endif
     enddo
  enddo
  

!c Convert to kg/m2/s (injection height is applied in the main program)
!C Account for molar weight to convert to volume mixing ratio
  maxemis=0
  meanemis=0
  n=0

  do ix=1,numxgrid
     do jy=1,numygrid
        do it=1,nbtime
           emifeux(ix,jy,it)=emifeux(ix,jy,it)*weightmolar(0)/weightmolar(nemission)          
           if (emifeux(ix,jy,it).gt.maxemis)maxemis=emifeux(ix,jy,it)
           meanemis=meanemis+emifeux(ix,jy,it)
           n=n+1
        enddo
     enddo
  enddo
  write(*,*)trim(pref),'max  emifeux: ',maxemis,maxval(emifeux)
  write(*,*)trim(pref),'mean emifeux: ',meanemis/n

  return
end subroutine reademissions_feux
