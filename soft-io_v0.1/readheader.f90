subroutine readheader(filename,nxmax,numxgrid,nymax,numygrid,&
     nzmax,numzgrid,outlon0,outlat0,dxout,dyout,outheight,ibdate,&
     ibtime,loutstep,maxspec,nspec,maxageclass,nageclass,lage,&
     ireleasestart,ireleaseend,maxpoint,numpoint,xpoint,ypoint,&
     zpoint1,zpoint2,heightnn,area,compoint)

  implicit none

  real pi,r_earth,pih,outlon0,outlat0,dxout,dyout
  real cosfact,cosfactm,cosfactp,gridarea,hzone,xp2,yp2
  real ylata,ylatm,ylatp
  integer nxmax,numxgrid,nymax,numygrid,nzmax,numzgrid
  integer ibdate,ibtime,loutstep,maxspec,nspec,maxageclass
  integer nageclass,maxpoint,numpoint,jjjjmmdd,ihmmss
  integer ix,j,jy,loutaver,loutsample,ltopo,method,n,i
  parameter(pi=3.14159265,r_earth=6.371e6,pih=pi/180.)

  character filename*150,compoint(maxpoint)*45,species(maxspec)*7
  real outheight(nzmax),heightnn(nxmax,nymax,0:nzmax)
  integer ireleasestart(maxpoint),ireleaseend(maxpoint)
  integer npart(maxpoint),kind(maxpoint),lage(0:maxageclass)
  real xpoint(maxpoint),ypoint(maxpoint),zpoint1(maxpoint)
  real zpoint2(maxpoint),xmass(maxpoint,maxspec)
  real oro(nxmax,nymax),area(nxmax,nymax)

  character*250 pref
  pref='         [readheader]:'

  write(*,*) trim(pref),'readheader'
  open(10,file=filename,form='unformatted',status='old')
  read(10) ibdate,ibtime
  write(*,*) trim(pref),'Simulation starts ',ibdate,ibtime
  read(10) loutstep,loutaver,loutsample
  read(10) outlon0,outlat0,numxgrid,numygrid,dxout,dyout
  write(*,*) trim(pref),'Lon start, delta, numpoints : ',outlon0,dxout,numxgrid
  write(*,*) trim(pref),'Lat start, delta, numpoints : ',outlat0,dyout,numygrid

  read(10) numzgrid,(outheight(i),i=1,numzgrid)
  write(*,*)trim(pref),'numzgrid:',numzgrid
  write(*,*)trim(pref),'outheight:',outheight

  read(10) jjjjmmdd,ihmmss
  read(10) nspec
  write(*,*)trim(pref),'nspec(1):',nspec
  nspec=nspec/3
  write(*,*)trim(pref),'nspec(2):',nspec
  do n=1,nspec
     read(10) numzgrid,species(n)
     read(10) numzgrid,species(n)
     read(10) numzgrid,species(n)
  enddo

  read(10) numpoint
  do i=1,numpoint
     read(10) ireleasestart(i),ireleaseend(i)

     read(10) xpoint(i),ypoint(i),xp2,yp2,zpoint1(i),zpoint2(i)
     read(10) npart(i),kind(i)
     read(10) compoint(i)
     do j=1,nspec
        read(10)
        read(10)
        read(10) xmass(i,j)
     enddo
  enddo
  read(10) method
  read(10) nageclass,(lage(i),i=1,nageclass)
  lage(0)=0
  write(*,*) trim(pref),(lage(i),i=0,nageclass)

  close(10)
  write(*,*) trim(pref),'Levels:',(outheight(i),i=1,numzgrid)

  if (loutstep.lt.0) nspec=numpoint
  write(*,*) trim(pref),'nspec(3):',nspec


  ! Calculate height, which is outheight plus topography
  !*****************************************************

  do ix=1,numxgrid
     do jy=1,numygrid
        if (ltopo.eq.1) then
           heightnn (ix,jy,0) = oro(ix,jy)
        else
           heightnn (ix,jy,0) = 0.
        endif
        do i=1,numzgrid
           if (ltopo.eq.1) then
              heightnn (ix,jy,i) = outheight(i) + oro(ix,jy)
           else
              heightnn (ix,jy,i) = outheight(i)
           endif
        enddo
     enddo
  enddo



  ! Determine area of each output grid box
  !***************************************

  do jy=1,numygrid
     ylata=outlat0+(float(jy-1)+0.5)*dyout
     ylatp=ylata+0.5*dyout
     ylatm=ylata-0.5*dyout
     if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
        hzone=dyout*r_earth*pih
     else
        cosfact=cos(ylata*pih)*r_earth
        cosfactp=cos(ylatp*pih)*r_earth
        cosfactm=cos(ylatm*pih)*r_earth
        if (cosfactp.lt.cosfactm) then
           hzone=sqrt(r_earth**2-cosfactp**2)-&
                 sqrt(r_earth**2-cosfactm**2)
        else
           hzone=sqrt(r_earth**2-cosfactm**2)-&
                 sqrt(r_earth**2-cosfactp**2)
        endif
     endif
     gridarea=2.*pi*r_earth*hzone*dxout/360.
     do ix=1,numxgrid
        area(ix,jy)=gridarea
     enddo
  enddo

end subroutine readheader
