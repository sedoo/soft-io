      subroutine reademissions_solene(nemission,outlon0,outlat0,
     +numxgrid,numygrid,dxout,dyout,nxmax,nymax,nzmax,area,heightnn,
     +emifeux,index_feux,nbtime,date_flex)

* Reads biomass burning emissions by S. Turquety (2007) for the ICARTT campaign

c      use netcdf

      implicit none

      INCLUDE 'netcdf.inc'

      integer nemission,numxgrid,numygrid,nxmax,nymax,nzmax,nbtime
      real outlon0,outlat0,dxout,dyout

      integer m,i,j
      real area(nxmax,nymax)
      real heightnn(nxmax,nymax,0:nzmax)

      real emifeux(nxmax,nymax,nbtime)
      integer index_feux(nxmax,nymax),indtime(nbtime)
      double precision date_flex(nbtime),juldate,julian_nc
      real pi,r_earth,pih,xl,xr,xm,yl,ym,yr
      integer maxemissions,maxmegaemissions,it
      integer maxpoint,indx,indy,ix,jy,n,neur,yeardate,timeday,year
      parameter(pi=3.14159265,r_earth=6.371e6,pih=pi/180.)
      parameter(maxemissions=4)
      parameter(maxmegaemissions=4)
      real emegcell(12),maxemis,xdiff,ydiff
      real emission_tami(360,180)
      parameter(maxpoint=1000000)
      real xlonl(maxpoint),xlonr(maxpoint),ylatl(maxpoint)
      real ylatr(maxpoint),zl(maxpoint),zr(maxpoint)

      character*256 emislist,contlist,contfile

      real weightmolar(0:maxemissions)
      data weightmolar/28.97,28.,46.,64.,12.01/

      integer lon_dimid,lat_dimid,time_dimid
      integer lon_varid,lat_varid
      integer nblon_nc,nblat_nc,nbtime_nc

      integer ufile_cont,cont_varid,time_varid

      real,allocatable :: lon_nc(:),lat_nc(:),emis_nc(:,:,:)
      real,allocatable :: cont_nc(:,:,:)

      integer nblon_sol,nblat_sol
      parameter (nblon_sol=360,nblat_sol=180)
      real lon_sol(nblon_sol),lat_sol(nblat_sol)
      real temp_sol(nblon_sol,nblat_sol)
      real emis_sol(nblon_sol,nblat_sol,nbtime)
      real*8 time_sol,ref_sol
      character*256 file_sol

*****************************
*** Read continents index ***
*****************************

!      contlist='/home/auba/couple_emissions/continents_file_feux.txt'
      contlist='continents_file_feux.txt'
      open(44,file=contlist)
      read(44,'(a)') contfile
      close(44)

      call check(nf_open(trim(contfile),NF_NOWRITE,ufile_cont))

      call check(nf_inq_dimid(ufile_cont,'lon',lon_dimid))
      call check(nf_inq_dimlen(ufile_cont,lon_dimid,nblon_nc))
      call check(nf_inq_dimid(ufile_cont,'lat',lat_dimid))
      call check(nf_inq_dimlen(ufile_cont,lat_dimid,nblat_nc))
      call check(nf_inq_dimid(ufile_cont,'date',time_dimid))
      call check(nf_inq_dimlen(ufile_cont,time_dimid,nbtime_nc))

      allocate(cont_nc(nblon_nc,nblat_nc,1))
      allocate(lon_nc(nblon_nc),lat_nc(nblat_nc))

      call check(nf_inq_varid(ufile_cont,'lon',lon_varid))
      call check(nf_get_var_real(ufile_cont,lon_varid,lon_nc))
      call check(nf_inq_varid(ufile_cont,'lat',lat_varid))
      call check(nf_get_var_real(ufile_cont,lat_varid,lat_nc))

      do ix=1,nblon_nc
       if (lon_nc(ix).gt.(outlon0+numxgrid*dxout)) then
        lon_nc(ix)=lon_nc(ix)-360.
       endif
       if (lon_nc(ix).lt.outlon0) then
        lon_nc(ix)=lon_nc(ix)+360.
       endif
      end do

      call check(nf_inq_varid(ufile_cont,'GFED_regions_0_5deg',
     +cont_varid))
      call check(nf_get_var_real(ufile_cont,cont_varid,cont_nc))
      call check(nf_close(ufile_cont))

********************************
*** Read emissions inventory ***
********************************
      
      file_sol='/home/auba/couple_emissions/emissions_test/'//
     +'OUT_emiss_BB2004/bioburnCO_peat_varlin_turetsky.dailyICARTT.'//
     +'generic.1x1.2004.bpch'
      do ix=1,nblon_sol
       lon_sol(ix)=-180.5+ix
       if (lon_sol(ix).gt.(outlon0+numxgrid*dxout)) then
        lon_sol(ix)=lon_sol(ix)-360.
       endif
       if (lon_sol(ix).lt.outlon0) then
        lon_sol(ix)=lon_sol(ix)+360.
       endif
      end do ! ix
      do jy=1,nblat_sol
       lat_sol(jy)=-90.5+jy
      end do ! jy

      do it=1,nbtime
       ref_sol=juldate(19850101,0)
       time_sol=aint(date_flex(it)-ref_sol)*24. ! looking for hours since 1st january 1985
!       print*,'Entering read_bpch2 for time',time_sol
       call read_bpch2(trim(file_sol),'BIOBSRCE',4,time_sol,nblon_sol,
     +  nblat_sol,1,temp_sol) ! temp_sol in molecules/cm2/s, need kg/m2/s
       do ix=1,nblon_sol
        do jy=1,nblat_sol
         emis_sol(ix,jy,it)=temp_sol(ix,jy)*weightmolar(nemission)*1e-3*
     +    1e4/6.022e23 ! g/mol / 1000 (kg) * 1e4 (m2) / molecules/mol
        end do ! ix
       end do ! jy
      end do ! it

***************************************************
*** Now attribute emissions to model grid cells ***
***************************************************

      maxemis=0

      do ix=1,numxgrid ! start with empty emission grid
       do jy=1,numygrid
        do it=1,nbtime
         emifeux(ix,jy,1)=0.
        end do
        index_feux(ix,jy)=0
       enddo !ix
      enddo !jy

      do 10 ix=1,nblon_nc
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
       enddo ! n

       do 10 jy=1,nblat_nc
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

c Now, we know that cont_nc(ix,jy) goes inside index_feux(indx,indy)
c        do it=1,nbtime
c         emifeux(indx,indy,it)=emifeux(indx,indy,it)+
c     +    emis_nc(ix,jy,indtime(it))*
c     +    area_nc(ix,jy,indtime(it))*1e6/area(indx,indy)
c        enddo
        index_feux(indx,indy)=cont_nc(ix,jy,1)
10    continue

      do ix=1,numxgrid
       xdiff=360.
       indx=1
       xm=outlon0+(float(ix)-0.5)*dxout
       do n=1,nblon_sol
        if (abs(lon_sol(n)-xm).le.xdiff) then
         xdiff=abs(lon_sol(n)-xm)
         indx=n
        end if
       end do ! n
       do jy=1,numygrid
        ydiff=360.
        indy=1
        ym=outlat0+(float(jy)-0.5)*dyout
        do n=1,nblat_sol
         if (abs(lat_sol(n)-ym).le.ydiff) then
          ydiff=abs(lat_sol(n)-ym)
          indy=n
         endif
        end do ! n
        do it=1,nbtime
         emifeux(ix,jy,it)=emis_sol(indx,indy,it)
        end do ! it
       end do ! jy
      end do ! ix


c Convert to kg/m2/s (injection height is applied in the main program)
C Account for molar weight to convert to volume mixing ratio
      maxemis=0
      do 61 ix=1,numxgrid
       do 61 jy=1,numygrid
        do 61 it=1,nbtime
         emifeux(ix,jy,it)=emifeux(ix,jy,it)*weightmolar(0)/
     +      weightmolar(nemission)
61        if (emifeux(ix,jy,it).gt.maxemis)
     +      maxemis=emifeux(ix,jy,it)
      print*,'maxemis feux',maxemis
      return
      end
