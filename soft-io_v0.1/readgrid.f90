subroutine readgrid(filegrid,itime,nxmax,numxgrid,nymax,numygrid,&
  nzmax,numzgrid,numzgrid_grid,maxspec,nspec,nspec2,maxageclass,&
  nageclass,grid,lage,scaleconc,decayconstant)

  implicit none

  integer itime,nxmax,numxgrid,nymax,numygrid,nzmax,numzgrid,ix,jy,k
  integer numzgrid_grid
  integer maxspec,nspec,nspec2,maxageclass,nageclass,kz,n,lspec2
  integer nage
  real scaleconc,decayconstant

  real grid(numxgrid,numygrid,numzgrid_grid,nspec,nageclass)
  integer lage(0:maxageclass)

  !***************** new sparse output
  real smallnum
  integer fact,ii,ir,nline,iline,io
  real sparse_dump_r(nxmax*nymax*nzmax)
  integer sparse_dump_i(nxmax*nymax*nzmax)
  integer sp_count_i,sp_count_r
  !***************** new sparse output

  character*250 filegrid,charr
  character*250 pref
  pref='         [readgrid]:'
  
  open(10,file=filegrid,form='unformatted',status='old')
  
  nline = 0
  do
     read(10,*,iostat=io) charr
     if (io.ne.0) exit
     nline = nline + 1
  end do
  rewind(10)

  nline=1
  do iline=1,nline
     read(10) itime 

     ! Loop about all species
     !************************
     do k=1,nspec

        ! Loop about all age classes
        !****************************
        do nage=1,nageclass

           ! Initialize age class fields
           !*****************************
           do ix=1,numxgrid
              do jy=1,numygrid
                 do kz=1,numzgrid
                    grid(ix,jy,kz,k,nage)=0.
                 enddo
              enddo
           enddo

           ! Read wet deposition
           !*********************           
           fact=1
           read(10) sp_count_i
           read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
           read(10) sp_count_r
           read(10) (sparse_dump_r(ix),ix=1,sp_count_r)

           ! Read dry deposition
           !*********************           
           fact=1
           read(10) sp_count_i
           read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
           read(10) sp_count_r
           read(10) (sparse_dump_r(ix),ix=1,sp_count_r)

           ! Read concentrations
           !*********************
           fact=1
           read(10) sp_count_i
           read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
           read(10) sp_count_r
           read(10) (sparse_dump_r(ix),ix=1,sp_count_r)           
           ii=0
           do ir=1,sp_count_r
              if ((sparse_dump_r(ir)*fact).gt.smallnum) then
                 ii=ii+1
                 n=sparse_dump_i(ii)
                 fact=fact*(-1.)
              else
                 n=n+1
              endif
              kz=n/(numxgrid*numygrid)
              jy=(n-kz*numxgrid*numygrid)/numxgrid
              ix=n-numxgrid*numygrid*kz-numxgrid*jy
              !print*,'kz=',kz
              grid(ix+1,jy+1,kz,k,nage)=abs(sparse_dump_r(ir))
           enddo


           ! Sum up all age classes to total fields and scale values as defined by user
           !****************************************************************************
           do ix=1,numxgrid
              do jy=1,numygrid
                 do kz=1,numzgrid
                    grid(ix,jy,kz,k,nage)=grid(ix,jy,kz,k,nage)*scaleconc                    
                 enddo
              enddo
           enddo
 
           ! End species loop, end age class loop
           !**************************************
           if (k.eq.lspec2) write(*,*) trim(pref),'k=lspec2',k,lspec2
           if (k.eq.lspec2) exit
        enddo !(nage)
     enddo !(k)
  enddo !(iline)
  close(10)

end subroutine readgrid
