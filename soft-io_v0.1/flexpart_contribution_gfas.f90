program flexpartplot_latlon
  
  use netcdf

  implicit none
  
  INCLUDE 'netcdf.inc'

  integer nxmax,nymax,nzmax,maxtime,levinteg
  integer maxpoint,maxspec,maxageclass
  integer nbtime

  integer i,j,k,iage,ibd,ibt,ibdate,ibtime,ilev
  integer jjjjmmdd1,ihmmss1,jjjjmmdd2,ihmmss2
  integer inday,intime,it,itime,ix,jy,kz,is1,js1
  integer len,len_fl,loutstep,lspec1,lspec2,nspec
  integer mother_or_nest,nageclass,nemission,nphys_spec
  integer nspec2,ntra,numpoint,numxgrid,numygrid,numzgrid
  real scaleconc

  real outlon0,outlat0

  parameter(nxmax=720,nymax=360,nzmax=13,maxtime=30) 
  parameter(maxpoint=1000,maxspec=1000,maxageclass=1)

  double precision jul,julstart,juldate
  real outheight(nzmax),decayconstant,dxout,dyout,xm,ym
  integer ireleasestart(maxpoint),ireleaseend(maxpoint)
  integer lage(0:maxageclass)
  real xpoint(maxpoint),ypoint(maxpoint),zpoint1(maxpoint)
  real zpoint2(maxpoint),maxcontrib,realdum(5)

  real,allocatable :: grid(:,:,:,:,:)
  double precision, allocatable :: date_flex(:)
  real area(nxmax,nymax),heightnn(nxmax,nymax,0:nzmax)

  character file_overview*120,strdum*120,strdum2*120
  character file_output*120,file_feux*120
  character dirname*120,datesname*120,flexdir*120
  character filename*150,aday*8,atime*6
  character*4 aspec
  character*250 filegrid,fmt*23
  character compoint(maxpoint)*45,anthrostr*64,firestr*64

  real,allocatable :: emissions(:,:,:),emifeux(:,:,:),frpfire(:,:,:),hfire(:,:,:)
  real frp,fre,firetop,inj
  integer index_cont(nxmax,nymax),index_feux(nxmax,nymax)

  integer nb_cont_anthro,nb_cont_fire
  parameter (nb_cont_anthro=14) ! number of geographical areas to separate contributions
  parameter (nb_cont_fire=14) ! same, for fire emissions
  double precision,allocatable :: agespectrum(:,:,:,:)
  double precision,allocatable :: firespectrum(:,:,:,:)

  integer maxtra
  parameter(maxtra=1000)
  real xtra(maxtra,maxtime),ytra(maxtra,maxtime)
  real ztra(maxtra,maxtime)
  real ablfract(maxtra,maxtime)
  real pvfract(maxtra,maxtime)
  real trfract(maxtra,maxtime),stratfract(maxtra,maxtime)
  real stratfract_sum(maxtra)

  real xclust(5,maxtra,maxtime),yclust(5,maxtra,maxtime)
  real zclust(5,maxtra,maxtime),fclust(5,maxtra,maxtime)
  real rmsclust(5,maxtra,maxtime)
  integer is(maxtra,maxtime),numtime(maxtra)

  double precision traj_jultrastart(maxtra)
  double precision traj_jultraend(maxtra)
  double precision traj_jultra(maxtra,maxtime)

  double precision,allocatable :: anthro_total(:,:)
  double precision,allocatable :: fire_total(:,:)

  double precision,allocatable :: dente_spec(:,:,:),dente_total(:)
  double precision dente_temp
  real,allocatable :: dente_hinj(:,:)
  character dentefile*150

  double precision,allocatable :: prm_spec(:,:,:),prm_total(:)
  double precision prm_temp
  real,allocatable :: prm_hinj(:,:)
  character prmfile*150

  double precision,allocatable :: h_spec(:,:,:),h_total(:)
  double precision h_temp
  

  double precision contribution,d1,d2,d3,d4,d5,d6,d7,d8,gridsum
  double precision contribfire,gridfire

  integer nbntra,nbflexdir,iflexdir,io
  
  integer anth1_both2_bbur3,flag300m,numzgrid_grid
  character outdir*150,flight*150

  character*250 pref

  pref='      [flexpart_contrib_files]:'
  flag300m=0

  !-------------------------------------------------


  open(10,file='flexpart_contribution_levels.nml')
  read(10,'(a)') anthrostr
  read(10,'(a)') firestr
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(a)') dentefile
  read(10,'(a)') prmfile
  read(10,'(a)') outdir
  read(10,'(a)') 
  read(10,'(a)') 
  read(10,'(i1)') anth1_both2_bbur3 
  read(10,'(a)') flexdir
  close(10)
  outdir= trim(outdir)//'/'

  write(*,*)trim(pref),'outdir ',outdir
  write(*,*)trim(pref),'anth1_both2_bbur3:',anth1_both2_bbur3

  nbflexdir=1
  write(*,*) flexdir
!  file_overview='run_over_following_files.txt'
!  open(400,file=file_overview)
!  nbflexdir=0
!  do
!     read(400,'(a)',iostat=io) flexdir
!     if (io.ne.0) exit
!     nbflexdir=nbflexdir+1
!  end do
!  rewind(400)





  do iflexdir=1,nbflexdir
     !read(400,'(a)') flexdir
     len_fl=index(flexdir,' ')-1
     dirname= flexdir(1:len_fl)//'/'
     len=index(dirname,' ')-1
     datesname=dirname(1:len)//'dates'

     write(*,*) trim(pref),'dirname:',trim(dirname),'   ',len,len_fl
     do k=1,(len-1)
        if(dirname(k:k).eq.'/')then 
           flight=trim(dirname(k+1:len))
        endif
     enddo


     write(*,*)
     write(*,*) trim(pref),'*****************************************************'
     write(*,*) trim(pref),'* Coupling FLEXPART residence times with emissions. *'
     write(*,*) trim(pref),'*****************************************************'
     write(*,*)
     write(*,*) trim(pref),'Input : ',trim(datesname)

     mother_or_nest=1 ! only mother grid, no nests

     !**********************************
     ! Read FLEXPART header information
     !**********************************
     write(*,*)
     write(*,*) trim(pref),'************* Read FLEXPART header information'
     write(*,*) trim(pref),trim(outdir)

     filename=dirname(1:len)//'header'
     write(*,*) trim(pref),trim(outdir)
     call readheader(filename,nxmax,numxgrid,nymax,numygrid,nzmax,&
          numzgrid,outlon0,outlat0,dxout,dyout,outheight,ibdate,ibtime,&
          loutstep,maxspec,nspec,maxageclass,nageclass,lage,ireleasestart,&
          ireleaseend,maxpoint,numpoint,xpoint,ypoint,zpoint1,zpoint2,&
          heightnn,area,compoint)    
     write(*,*) trim(pref),'Exit from readheader'
     write(*,*)


     if(numzgrid.eq.12)then
        flag300m=0
        numzgrid_grid=12
     else
        flag300m=1
        numzgrid=12
        numzgrid_grid=13
     endif

     allocate(grid(numxgrid,numygrid,numzgrid_grid,nspec,nageclass)) !rvrv
     
     if(anth1_both2_bbur3.le.2)then
        allocate(agespectrum(nspec,maxtime,nb_cont_anthro+1,numzgrid))
        allocate(anthro_total(nb_cont_anthro+1,numzgrid))
     endif
     if(anth1_both2_bbur3.ge.2)then
        allocate(firespectrum(nspec,maxtime,nb_cont_fire+1,numzgrid))
        allocate(fire_total(nb_cont_fire+1,numzgrid))
        allocate(dente_spec(nspec,maxtime,nb_cont_fire+1))
        allocate(dente_total(nb_cont_fire+1),dente_hinj(numzgrid,4))
        allocate(prm_spec(nspec,maxtime,nb_cont_fire+1))
        allocate(prm_total(nb_cont_fire+1),prm_hinj(numzgrid,6))
        allocate(h_spec(nspec,maxtime,nb_cont_fire+1),h_total(nb_cont_fire+1))
     endif




     julstart=juldate(ibdate,ibtime)
     lspec2=nspec
     lspec1=1

     !*********************************
     ! Read dentener injection heights
     !*********************************
     if(anth1_both2_bbur3.ge.2)then
        write(*,*)
        write(*,*) trim(pref),'************* Read dentener injection heights'
        !dentefile='hinj_dentener.txt'
        open(19,file=dentefile,form='formatted',status='old')
        read(19,*)
        do kz=1,numzgrid
           read(19,*) realdum(1),realdum(2),realdum(3),realdum(4),realdum(5)
           dente_hinj(kz,1)=realdum(2)
           dente_hinj(kz,2)=realdum(3)
           dente_hinj(kz,3)=realdum(4)
           dente_hinj(kz,4)=realdum(5)
        end do
        close(19)
     endif

     !****************************
     ! Read PRM injection heights
     !****************************
     if(anth1_both2_bbur3.ge.2)then
        write(*,*)
        write(*,*) trim(pref),'************* Read PRM injection heights'
        !prmfile='hinj_prmfrp.txt'
        open(199,file=prmfile,form='formatted',status='old')
        read(199,*)
        do kz=1,numzgrid
           read(199,*) realdum(1),realdum(2),realdum(3),realdum(4)
           prm_hinj(kz,1)=dente_hinj(kz,1)
           prm_hinj(kz,2)=dente_hinj(kz,2)
           prm_hinj(kz,3)=dente_hinj(kz,2)
           prm_hinj(kz,4)=realdum(2) ! boreal small fire
           prm_hinj(kz,5)=realdum(3) ! boreal medium fire
           prm_hinj(kz,6)=realdum(4) ! boreal big fire
        end do
        close(199)
     endif

     !***********************************
     ! Read all the emission information
     !***********************************
        write(*,*)
        write(*,*) trim(pref),'************* Read all the emission information'
        open(20,file=trim(outdir)//trim(flight)//'dates',form='formatted',status='old')
        nbtime=0
        do
           read(20,'(i8,i6)',iostat=io) inday,intime
           if (io.ne.0) exit
           nbtime=nbtime+1
        end do
        rewind(20)
        write(*,*) trim(pref),'nbtime:',nbtime

        allocate(date_flex(nbtime))
        
        do it=1,nbtime
           read(20,'(i8,i6)') inday,intime
           date_flex(it)=juldate(inday,intime)
        enddo
        rewind(20)

        nemission=1

        if(anth1_both2_bbur3.le.2)then
           allocate(emissions(nxmax,nymax,nbtime))
           call reademissions_netcdf(nemission,outlon0,outlat0,numxgrid,&
                numygrid,dxout,dyout,nxmax,nymax,nzmax,area,heightnn,&
                emissions,index_cont,nbtime,date_flex)
           write(*,*) trim(pref),''
           write(*,*) trim(pref),'Exit from reademissions anthro'
           write(*,*)
        endif

        if(anth1_both2_bbur3.ge.2)then
           allocate(emifeux(nxmax,nymax,nbtime))
           allocate(frpfire(nxmax,nymax,nbtime),hfire(nxmax,nymax,nbtime))
           call reademissions_feux(nemission,outlon0,outlat0,numxgrid,&
                numygrid,dxout,dyout,nxmax,nymax,nzmax,area,heightnn,&
                emifeux,index_feux,nbtime,date_flex,frpfire,hfire)
           write(*,*) trim(pref),'Exit from reademissions feux'
           write(*,*)
        endif

     !********************************************
     ! Initialize plotting fields and age spectra
     !********************************************
     write(*,*)
     write(*,*) trim(pref),'************* Initialize plotting fields'
     write (*,*) trim(pref),'lspe: ',lspec1,lspec2
     if(anth1_both2_bbur3.le.2)then
        do k=1,nspec
           do iage=1,maxtime              
              do kz=1,numzgrid
                 do i=1,nb_cont_anthro+1
                    agespectrum(k,iage,i,kz)=0.
                 end do
              enddo
           enddo
        enddo
     endif
     if(anth1_both2_bbur3.ge.2)then
        firespectrum=0.
        dente_spec=0.
        prm_spec=0.
        h_spec=0.

!        do k=1,nspec
!           do iage=1,maxtime
!              do i=1,nb_cont_fire+1
!                 dente_spec(k,iage,i)=0.
!                 prm_spec(k,iage,i)=0.
!              end do
!              do kz=1,numzgrid                 
!                 do i=1,nb_cont_fire+1
!                    firespectrum(k,iage,i,kz)=0.
!                 end do
!              enddo
!           enddo
!        enddo
     endif



     !**************************
     ! Read sorted trajectories
     !**************************
     write(*,*)
     write(*,*) trim(pref),'************* Read sorted trajectories'
     write(*,*) trim(pref),trim(flight)
     open(3010,file=trim(outdir)//trim(flight)//'sorted_trajectories.txt')
     write(*,*)
     write(*,*) trim(pref),'Reading sorted traj :',trim(outdir)//trim(flight),'sorted_trajectories.txt'
     write(*,*)

     nbntra=0
     do
        read(3010,'(i8,i6)',iostat=io) inday,intime
        nbntra=nbntra+1
        if (io/=0) exit
     end do
     rewind(3010)
     do ntra=1,nbntra

        stratfract_sum(ntra)=0
        read(3010,*) jjjjmmdd1,ihmmss1,jjjjmmdd2,ihmmss2,numtime(ntra)
        traj_jultrastart(ntra)=juldate(jjjjmmdd1,ihmmss1)
        traj_jultraend(ntra)=juldate(jjjjmmdd2,ihmmss2)
        do i=1,numtime(ntra)
           read(3010,98) is(ntra,i),xtra(ntra,i),ytra(ntra,i),&
                ztra(ntra,i),d1,d2,d3,d4,d5,d6,d7,d8,ablfract(ntra,i),&
                pvfract(ntra,i),trfract(ntra,i),&
                (xclust(j,ntra,i),yclust(j,ntra,i),zclust(j,ntra,i),&
                fclust(j,ntra,i),rmsclust(j,ntra,i),j=1,5)
           traj_jultra(ntra,i)=traj_jultrastart(ntra)+dble(float(is(ntra,i))/86400)

           if (abs(ytra(ntra,i)).gt.30.) then   ! use PV definition
              stratfract(ntra,i)=100.-pvfract(ntra,i)
           else                              ! use thermal definition
              stratfract(ntra,i)=100.-trfract(ntra,i)
           endif
           ! stratfract for 1st = -99 and for 21st 0 - leave them out      
           if ((i.gt.1).and.(i.lt.21))  then
              stratfract_sum(ntra)=stratfract_sum(ntra)+stratfract(ntra,i)
           endif
        enddo
98      format(i8,2f9.4,4f8.1,f8.2,4f8.1,3f6.1,5(2f8.3,f7.0,f6.1,f8.1))
     enddo
     ntra=ntra-1

     !*******************************************
     ! Loop about all times, given in file dates
     !*******************************************
     write(*,*)
     write(*,*) trim(pref),'************* Loop about all times, given in file'
 
     !test frp
     file_output=trim(outdir)//trim(flight)//'frp_file.txt'
     open(250,file=trim(outdir)//trim(flight)//'frp_file.txt')


     !      open(20,file=datesname,form='formatted',status='old')

     rewind(20)
     do it=1,nbtime    !(20+2 days)
     !do it=1,1
        read(20,'(i8,i6)') inday,intime
        write(aday,'(i8.8)') inday
        write(atime,'(i6.6)') intime

        aspec(1:1)='_'
        nphys_spec=1
        scaleconc=1
        write(aspec(2:4),'(i3.3)') nphys_spec
        filegrid=filename(1:len)//'/grid_time_'//aday//atime//aspec
        !write(*,*) trim(pref),'readgrid ',trim(filegrid),inday

        !write(*,*) 'lspec2',lspec2
        call readgrid(filegrid,itime,nxmax,numxgrid,nymax,numygrid,&
             nzmax,numzgrid,numzgrid_grid,maxspec,nspec,nspec2,maxageclass,nageclass,grid,&
             lage,scaleconc,decayconstant)
        !write(*,*) '!!>>>',grid(100,134,1,1,1),grid(100,134,1,2,1)
        !write(*,*) '!!>>>',grid(100,134,1,:,1) 




        !************************************************************************
        ! Add contributions from this time step to age spectra and gridded field
        !************************************************************************
        write(*,*) trim(pref),'inday (lspec1,lspec2)',inday,'(',lspec1,lspec2,')'

        

        do k=lspec1,lspec2    !(release boxes)
        !do k=1,30
           iage=abs(itime-ireleaseend(k))/abs(loutstep)+1
           !write(*,*) '///',it,k,iage,itime,ireleaseend(k),loutstep
           !if(it.eq.2.and.k.eq.3)stop


           if (iage.le.30) then
              do ix=1,numxgrid
                 do jy=1,numygrid
                    ym=outlat0+(float(jy)-0.5)*dyout ! for dentener contribution
                    do kz=1,numzgrid-1 ! then, all levels                               

                       
                       if(flag300m.eq.0)then 
                          if (kz.eq.1) then
                             if(anth1_both2_bbur3.le.2)&
                                  contribution=grid(ix,jy,kz,k,nageclass)*emissions(ix,jy,it)/outheight(kz)
                             if(anth1_both2_bbur3.ge.2)&
                                  contribfire=grid(ix,jy,kz,k,nageclass)*emifeux(ix,jy,it)/outheight(kz)
                          else
                             if(anth1_both2_bbur3.le.2)&
                                  contribution=grid(ix,jy,kz,k,nageclass)*emissions(ix,jy,it)/&
                                  (outheight(kz)-outheight(kz-1))
                             if(anth1_both2_bbur3.ge.2)&
                                  contribfire=grid(ix,jy,kz,k,nageclass)&
                                  *emifeux(ix,jy,it)/(outheight(kz)-outheight(kz-1))
                          endif
                       else ! add the two first layers for bioburn, and shift the other ones
                          if (kz.eq.1) then
                             if(anth1_both2_bbur3.le.2)&
                                  contribution=grid(ix,jy,kz,k,nageclass)*emissions(ix,jy,it)/outheight(kz)
                             if(anth1_both2_bbur3.ge.2)&
                                  contribfire=(grid(ix,jy,1,k,nageclass)+grid(ix,jy,2,k,nageclass))*&
                                  emifeux(ix,jy,it)/outheight(2)
                          else
                             if(anth1_both2_bbur3.le.2)&
                                  contribution=grid(ix,jy,kz,k,nageclass)&
                                  *emissions(ix,jy,it)/(outheight(kz)-outheight(kz-1))
                             if(anth1_both2_bbur3.ge.2)&
                                  contribfire=grid(ix,jy,kz+1,k,nageclass)&
                                  *emifeux(ix,jy,it)/(outheight(kz+1)-outheight(kz-1+1))
                          endif
                       endif                
                       if(anth1_both2_bbur3.le.2)&
                            agespectrum(k,iage,nb_cont_anthro+1,kz)=&
                            agespectrum(k,iage,nb_cont_anthro+1,kz)+contribution
                       !if(anth1_both2_bbur3.ge.2 .and. contribfire.gt.0.)&
                            !total firespectrum in last column
                       !     firespectrum(k,iage,nb_cont_fire+1,kz)=&
                       !     firespectrum(k,iage,nb_cont_fire+1,kz)+contribfire
                       ! Determine continental contributions for anthropogenic
                       if ((index_cont(ix,jy).ge.1).and.(index_cont(ix,jy).le.nb_cont_anthro))&
                            agespectrum(k,iage,index_cont(ix,jy),kz)=&
                            agespectrum(k,iage,index_cont(ix,jy),kz)+contribution


                       !if(k.le.2.and.kz.eq.1.and.index_cont(ix,jy).eq.1.and.contribution.ne.0)then
                       !   write(*,*) it,k,sum(agespectrum(k,:,index_cont(ix,jy),kz)),&
                       !        iage,nageclass,&
                       !        agespectrum(k,iage,index_cont(ix,jy),kz),ix,jy,contribution,&
                       !        grid(ix,jy,kz,k,nageclass),emissions(ix,jy,it),outheight(kz)
                       !endif


                       if(anth1_both2_bbur3.ge.2 .and. contribfire.gt.0.)then
                          ! Determine continental contributions for fire
                          !if ((index_feux(ix,jy).ge.1).and.(index_feux(ix,jy).le.nb_cont_fire))&
                          !     firespectrum(k,iage,index_feux(ix,jy),kz)=&
                          !     firespectrum(k,iage,index_feux(ix,jy),kz)+contribfire
                          ! Apply injection profile for fire emissions
                          dente_temp=0.
                          !if (abs(ym).ge.60.) then
                          !   if (index_cont(ix,jy).eq.1) then ! Canada
                          !      dente_temp=contribfire*dente_hinj(kz,4)
                          !   else ! Other boreal
                          !      dente_temp=contribfire*dente_hinj(kz,3)
                          !   endif
                          !else if (abs(ym).ge.30.) then ! Temperate
                          !   dente_temp=contribfire*dente_hinj(kz,2)
                          !else ! Tropical
                          !   dente_temp=contribfire*dente_hinj(kz,1)
                          !endif
                          !dente_spec(k,iage,nb_cont_fire+1)=&
                          !     dente_spec(k,iage,nb_cont_fire+1)+dente_temp
                          !if ((index_feux(ix,jy).ge.1).and.(index_feux(ix,jy).le.nb_cont_fire))&
                          !     dente_spec(k,iage,index_feux(ix,jy))=&
                          !     dente_spec(k,iage,index_feux(ix,jy))+dente_temp

                          prm_temp=0.
                          if (abs(ym).ge.60.) then
                             ! area range : 0.108E+9 - 12.36E+9 m2 (readheader.f90)
                             ! conversion factor beta_l = 1.55 kgDM.MJ-1 for extratropical forest with organic soil (Kaiser et al., 2012, Tab2)
                             ! CO emission factor = 106 gCO.kgDM-1 (Kaiser et al., 2012, Tab3)
                             ! frp compared to 10 or 100 TJ.day-1
                             fre=emifeux(ix,jy,it)*area(ix,jy)*24.*3600./(106.*1.55) !(GJ.day-1)
                             frp=frpfire(ix,jy,it)*24.*3600.
                             !write(900,*) ix,jy,fre,frp,hfire(ix,jy,it)
                             if (frp.gt.1e4) then
                                if (frp.gt.1e5) then
                                   prm_temp=contribfire*prm_hinj(kz,6)
                                else
                                   prm_temp=contribfire*prm_hinj(kz,5)
                                endif
                             else
                                prm_temp=contribfire*prm_hinj(kz,4)
                             endif
                          else if (abs(ym).ge.30.) then ! Temperate
                             prm_temp=contribfire*prm_hinj(kz,2)
                          else ! Tropical
                             prm_temp=contribfire*prm_hinj(kz,1)
                          endif
                          prm_spec(k,iage,nb_cont_fire+1)=&
                               prm_spec(k,iage,nb_cont_fire+1)+prm_temp
                          if ((index_feux(ix,jy).ge.1).and.(index_feux(ix,jy).le.nb_cont_fire))&
                               prm_spec(k,iage,index_feux(ix,jy))=&
                               prm_spec(k,iage,index_feux(ix,jy))+prm_temp

                          !report contribution with mean plumes altitude
                          h_temp=0.
                          inj=0.
                          firetop=hfire(ix,jy,it)
                           if (kz.eq.1) then
                            if (firetop.le.outheight(1)) then
                             !contribfire=contribfire*outheight(1)/firetop
                             inj=1.0

                             !write(*,'(2(i4,1x),2(i3,1x),2(f8.2,1x),2(e10.3,1x))')&
                             !ix,jy,kz,it,outheight(kz),firetop,inj,contribfire

                            else
                             inj=outheight(1)/firetop
                            endif
                           else
                            if (firetop.ge.outheight(kz)) then
                             inj=(outheight(kz)-outheight(kz-1))/firetop
                            else
                              if (firetop.gt.outheight(kz-1)) then
                               !contribfire=contribfire*(outheight(kz)-outheight(kz-1))/(firetop-outheight(kz-1))
                               inj=(firetop-outheight(kz-1))/firetop
                               !write(*,'(2(i4,1x),2(i3,1x),2(f8.2,1x),2(e10.3,1x))')&
                               !ix,jy,kz,it,outheight(kz),firetop,inj,contribfire
                              else
                              inj=0.0
                             endif
                            endif
                           endif

                          h_temp=contribfire*inj
                          h_spec(k,iage,nb_cont_fire+1)=&
                           h_spec(k,iage,nb_cont_fire+1)+h_temp
                          if ((index_feux(ix,jy).ge.1).and.(index_feux(ix,jy).le.nb_cont_fire))&
                           h_spec(k,iage,index_feux(ix,jy))=&
                            h_spec(k,iage,index_feux(ix,jy))+h_temp

                          !if (inj.gt.0.) then
                           !frp=frpfire(ix,jy,it)*24.*3600.
                           !write(250,'(2(i4,1x),2(i3,1x),2(f8.2,1x),2(e10.3,1x))')&
                            !ix,jy,kz,it,outheight(kz),firetop,inj,contribfire
                          !endif

                       endif !test fire flag and contribfire>0
                    enddo !(kz)
                 enddo !(jy)
              enddo !(ix)
           endif !(iage.le.30)
        enddo !(k)
     enddo !(it)
     close(20)



     !**********************************************
     ! Sum contributions over time and write output
     !**********************************************
     write(*,*) trim(pref),'************* Sum contributions over time and write'


     strdum=''
     write(strdum,'(i5)') floor(outheight(1))
     
     if(anth1_both2_bbur3.le.2)then
        write(*,*) '!!! n2: ',sum(agespectrum(1,:,1,1))
        write(*,*) '!!! n2: ',sum(agespectrum(2,:,1,1))

        write(*,*) trim(pref),'************************** ANTHROPIC...'
        file_output=trim(outdir)//trim(flight)//'co_'//trim(anthrostr)//'_0-'//trim(adjustl(strdum))//'m.txt'
        open(101,file=file_output)
        write(101,*) 0.
        write(101,*) outheight(1)
        write(*,*) trim(pref),'file_output: ',file_output
        do kz=2,numzgrid
           strdum=''
           strdum2=''
           write(strdum,'(i5)') floor(outheight(kz-1))
           write(strdum2,'(i5)') floor(outheight(kz))
           file_output=trim(outdir)//trim(flight)//'co_'//trim(anthrostr)//'_'//&
                trim(adjustl(strdum))//'-'//trim(adjustl(strdum2))//'m.txt'
           open(100+kz,file=file_output)
           write(100+kz,*) outheight(kz-1)
           write(100+kz,*) outheight(kz)
           write(*,*) trim(pref),'file_output: ',file_output           
        end do                       
        do k=lspec1,lspec2 ! loop on all release points
           do kz=1,numzgrid ! then, for all levels
              do iage=2,30
                 do i=1,nb_cont_anthro+1
                    agespectrum(k,iage,i,kz)=agespectrum(k,iage,i,kz)+  agespectrum(k,iage-1,i,kz)
                 enddo                 
              enddo ! iage loop
              do i=1,nb_cont_anthro+1
                 anthro_total(i,kz)=agespectrum(k,30,i,kz)
              end do            
              len=index(compoint(k),' ')-1
              strdum2=''
              write(strdum2,*) nb_cont_anthro+1
              strdum= '(a,1x,'//trim(adjustl(strdum2))//'(e12.4,1x),f10.5)'
              write(100+kz,strdum) trim(compoint(k)(1:len)),(anthro_total(i,kz),i=1,nb_cont_anthro+1),stratfract_sum(k)              
           enddo
        enddo
        do kz=1,numzgrid
           close(100+kz)
        end do     
        deallocate(agespectrum,anthro_total)
        deallocate(emissions)
     endif


     if(anth1_both2_bbur3.ge.2)then
        write(*,*) trim(pref),'************************** BIOMASS BURNING...'
        !strdum=''
        !write(strdum,'(i5)') floor(outheight(1))
        !file_feux=trim(outdir)//trim(flight)//'co_'//trim(firestr)//'_0-'//trim(adjustl(strdum))//'m.txt'
        !open(151,file=file_feux)
        !write(151,*) 0.
        !write(151,*) outheight(1)
        !write(*,*) trim(pref),'file_feux: ',file_feux
        !do kz=2,numzgrid
        !   strdum=''
        !   strdum2=''
        !   write(strdum,'(i5)') floor(outheight(kz-1))
        !   write(strdum2,'(i5)') floor(outheight(kz))           
        !   file_feux=trim(outdir)//trim(flight)//'co_'//trim(firestr)//'_'//&
        !        trim(adjustl(strdum))//'-'//trim(adjustl(strdum2))//'m.txt'
        !   open(150+kz,file=file_feux)
        !   write(150+kz,*) outheight(kz-1)
        !   write(150+kz,*) outheight(kz)
        !   write(*,*) trim(pref),'file_feux: ',file_feux
        !end do
        !open(200,file=trim(outdir)//trim(flight)//'co_'//trim(firestr)//'_dentener.txt')
        !write(200,*) 0.
        !write(200,*) 0.
        !write(*,*) trim(pref),trim(outdir)//trim(flight)//'co_'//trim(firestr)//'_dentener.txt'
        open(210,file=trim(outdir)//trim(flight)//'co_'//trim(firestr)//'_prmfrp.txt')
        write(210,*) 0.
        write(210,*) 0.
        write(*,*) trim(pref),trim(outdir)//trim(flight)//'co_'//trim(firestr)//'_prmfrp.txt'
 
        !hfile
        open(220,file=trim(outdir)//trim(flight)//'co_'//trim(firestr)//'_hmean.txt')
        write(220,*) 0.
        write(220,*) 0.
        write(*,*) trim(pref),trim(outdir)//trim(flight)//'co_'//trim(firestr)//'_hmean.txt'

        
        do k=lspec1,lspec2 ! loop on all release points

           !do i=1,nb_cont_fire+1 ! first,dentener
           !   do iage=2,30
           !      dente_spec(k,iage,i)= dente_spec(k,iage,i)+dente_spec(k,iage-1,i)
           !   end do
           !   dente_total(i)=dente_spec(k,30,i)
           !end do
           !strdum2=''
           !write(strdum2,*) nb_cont_fire+1
           !strdum= '(a,1x,'//trim(adjustl(strdum2))//'(e12.4,1x),f10.5)'
           !len=index(compoint(k),' ')-1
           !write(200,strdum) trim(compoint(k)(1:len)), (dente_total(i),i=1,nb_cont_fire+1),stratfract_sum(k)



           !write(*,'(a,1x,15e18.6,"///",2e18.6)') trim(compoint(k)(1:len)),(dente_total(i),i=1,nb_cont_fire+1),&
           !    & sum(dente_total(1:nb_cont_fire)),&
           !     & dente_total(nb_cont_fire+1)-sum(dente_total(1:nb_cont_fire))



           do i=1,nb_cont_fire+1 ! then, PRM
              do iage=2,30
                 prm_spec(k,iage,i)=  prm_spec(k,iage,i)+prm_spec(k,iage-1,i)
              end do
              prm_total(i)=prm_spec(k,30,i)
           end do
           strdum2=''
           write(strdum2,*) nb_cont_fire+1
           strdum= '(a,1x,'//trim(adjustl(strdum2))//'(e12.4,1x),f10.5)'
           len=index(compoint(k),' ')-1
           write(210,strdum) trim(compoint(k)(1:len)), (prm_total(i),i=1,nb_cont_fire+1),stratfract_sum(k)

           do i=1,nb_cont_fire+1 ! then, hmean
              do iage=2,30
                 h_spec(k,iage,i)=  h_spec(k,iage,i)+h_spec(k,iage-1,i)
              end do
              h_total(i)=h_spec(k,30,i)
           end do
           strdum2=''
           write(strdum2,*) nb_cont_fire+1
           strdum= '(a,1x,'//trim(adjustl(strdum2))//'(e12.4,1x),f10.5)'
           len=index(compoint(k),' ')-1
           write(220,strdum) trim(compoint(k)(1:len)), (h_total(i),i=1,nb_cont_fire+1),stratfract_sum(k)


           !do kz=1,numzgrid ! then, for all levels
           !   do iage=2,30                 
           !      do i=1,nb_cont_fire+1
           !         firespectrum(k,iage,i,kz)=firespectrum(k,iage,i,kz)+firespectrum(k,iage-1,i,kz)
           !      enddo
           !   enddo ! iage loop
           !   do i=1,nb_cont_fire+1
           !      fire_total(i,kz)=firespectrum(k,30,i,kz)
           !   end do
           !   len=index(compoint(k),' ')-1
           !   strdum2=''
           !   write(strdum2,*) nb_cont_fire+1
           !   strdum= '(a,1x,'//trim(adjustl(strdum2))//'(e12.4,1x),f10.5)'
           !   write(150+kz,strdum) trim(compoint(k)(1:len)),(fire_total(i,kz),i=1,nb_cont_fire+1),stratfract_sum(k)
           !enddo
        enddo

        !do kz=1,numzgrid
        !   close(150+kz)
        !end do
        !close(200)
        close(210)
        close(220)
        close(250)
        deallocate(firespectrum,fire_total)
        deallocate(dente_spec,dente_total,dente_hinj)
        deallocate(prm_spec,prm_total,prm_hinj)
        deallocate(h_spec,h_total)
        deallocate(emifeux,frpfire,hfire)
     endif
     deallocate(grid,date_flex)

  enddo
  close(400)

  write(*,*) 'END COUPLEMIS!!!'

end program
