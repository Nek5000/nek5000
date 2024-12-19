      subroutine octuple_split_multirounds(n_rounds)
      use SIZE
      use OCTUPLESIZE
      integer n_rounds,i_round,iflag
      integer flag

       write(6,*)  'export to separate files of a big mesh'
       write(6,*)  'witch file ?'
       write(6,*)  ' (1 vertices, 2 edges, 3 bcs, 0 zero file)'
       read (5,*) iflag

      do i_round = 1,n_rounds
       if (i_round.eq.n_rounds) then
       call octuple_split1(iflag) ! no allocate arrays, and exporting on-run
	   else
       call octuple_split0()  ! allocate arrays, no exporting on-run
       endif
      enddo

      return
      end
!-----------------------------------------------------------------
      subroutine octuple_split1(iflag)
!
! convert nek hex20 element into 8 hex20 elements.
!  no allocate big arrays
!  export to re2 file on-run.
! 1. to export vertices on differnt files
! 2. to export mid 
! 3. to export sidesets.
!
      use SIZE
      use OCTUPLESIZE

      integer iflag

      real*8 hexver1(3,27) ! hex coordinates
      real*8 hexver2(3,27,8) ! hex coordinate
      integer hexss1(6),hexss2(6,8)
      integer hexss3(6),hexss4(6,8)
	  
      real xm3(3,3,3),ym3(3,3,3),zm3(3,3,3) ! hex coordinates
      real bc3(5,6)
      real curve3(6)
      character(3) cbc3(6)
      
! needded by open re2
      character(80) hdr
      integer nBCre2
      integer*8 eftot2

      real*4 test
      data   test  / 6.54321 /

! needed by write xyz
      real*8  xx(8), yy(8), zz(8)
      real*8  rgroup, buf2(30)
      integer igr 

      integer isym2pre(8)   ! Symmetric-to-prenek vertex ordering
      data    isym2pre    / 1 , 2 , 4 , 3 , 5 , 6 , 8 , 7 /

      integer edge_mid(12)
      save    edge_mid
      data    edge_mid /2,6,8,4,20,24,26,22,10,12,18,16/

      real    x3(27),y3(27),z3(27),xyz(3,3)
      real    mide(3),dist0,dist1

      integer e3(3,12)
      save    e3
      data    e3 /  1, 2, 3,    3, 6, 9,    9, 8, 7,    7, 4, 1, &
                   19,20,21,   21,24,27,   27,26,25,   25,22,19, &
                    1,10,19,    3,12,21,    9,18,27,    7,16,25  /

      character(1) cc
	
      integer*8 nbc,nbct,ncurve
      real*8 rbc,rbct,rcurve

      integer*8 iel,iedge,nedge,iel_new,iel_old,ivert,ifc,ihex,inekvert
      integer i,j,k,l,li,ierr

      character(3) ch3


      nxs = 3-1  ! always nx1=ny1=nz1=3 here
      nys = 3-1
      nzs = 3-1

      igr    = 0
      rgroup = igr

      write(6,'(A)') 'Starting octuple splitting and exporting'

!==================================================================
! open re2 file for writing
      !write(6,*) 'please give re2 file name:'
      !call read_re2_name
	  
      if (iflag.eq.1) then
      re2name = 'vertices.re2'
      else if (iflag.eq.2) then
      re2name = 'edges.re2'
      else if (iflag.eq.3) then
      re2name = 'bcs.re2'
      else if (iflag.eq.0) then
      re2name = 'zero.re2'
      endif

      call byte_open(re2name,ierr)


      num_elem = etot*8
      eftot2 = eftot*8
      !num_dim = 3
      nBCre2 = 1
      if (num_elem.ne.eftot2)  nBCre2  = 2

      if (iflag.eq.1) then
      write(6,*) 'writing file header'
!  Write the header
      call blank   (hdr,80)    
      write(hdr,1) num_elem, num_dim, eftot2, nBCre2
!      write(hdr,1) num_elem, num_dim, eftot
    1 format('#v004',i16,i3,i16,i4,' hdr')
!    1 format('#v002',i16,i3,i16,i4,' hdr')
!    1 FORMAT('#v002',i9,i3,i9,' this is the hdr')
      call byte_write(hdr,20,ierr)         
      call byte_write(test,1,ierr)     ! write the endian discriminator
! done write header

!============================================================
! write xyz
      write(6,*) 'splitting and writing xyz'
      !call rzero8(buf2,30)

      do iel_old=1,etot
       do ivert = 1,27
       hexver1(1,ivert) = xm1(ivert,1,1,iel_old)
       hexver1(2,ivert) = ym1(ivert,1,1,iel_old)
       hexver1(3,ivert) = zm1(ivert,1,1,iel_old)
       enddo
  
       do ifc = 1,6
       hexss1(ifc) = bc(5,ifc,iel_old)
       hexss3(ifc) = bc(4,ifc,iel_old) 
       enddo
 
! quadratic splitting hex20 elements..
! 
       call hex_octuple_splitting(hexver1,hexver2,hexss1,hexss2,hexss3,hexss4) 

       do ihex = 1,8

          ! export to new re2 files
          ! write to vertice, midpoint, and then sideset.. to seperate files
          do inekvert = 1,27
             xm3(inekvert,1,1) = hexver2(1,inekvert,ihex)
             ym3(inekvert,1,1) = hexver2(2,inekvert,ihex)
             zm3(inekvert,1,1) = hexver2(3,inekvert,ihex)
          enddo

          l = 0
          do k=0,1
          do j=0,1
          do i=0,1
            l=l+1
            li=isym2pre(l)

            xx(li) = xm3(1+i*nxs,1+j*nys,1+k*nzs)
            yy(li) = ym3(1+i*nxs,1+j*nys,1+k*nzs)
            zz(li) = zm3(1+i*nxs,1+j*nys,1+k*nzs)

          enddo
          enddo
          enddo

          call byte_write  (rgroup, 2,ierr)
          call copy        (buf2(1) ,xx,8)
          call copy        (buf2(9) ,yy,8)
          call copy        (buf2(17),zz,8)
          call byte_write  (buf2(1) ,16, ierr)
          call byte_write  (buf2(9) ,16, ierr)
          call byte_write  (buf2(17),16, ierr)

        enddo  !        do ihex = 1,8

      enddo !       do iel_old=1,etot

      endif !       if (iflag.eq.1) then !

! done write header and vertices

!=============================
! write curve
      if (iflag.eq.2) then
      write(6,*) 'splitting and writing curve'
      !call rzero8(buf2,30)

      nedge = 12
      ncurve = 0

      iel_new = 0
      do iel_old=1,etot
       do ivert = 1,27
       hexver1(1,ivert) = xm1(ivert,1,1,iel_old)
       hexver1(2,ivert) = ym1(ivert,1,1,iel_old)
       hexver1(3,ivert) = zm1(ivert,1,1,iel_old)
       enddo
  
       do ifc = 1,6
       hexss1(ifc) = bc(5,ifc,iel_old)
       hexss3(ifc) = bc(4,ifc,iel_old) 
       enddo
 
       call hex_octuple_splitting(hexver1,hexver2,hexss1,hexss2,hexss3,hexss4) 

       do ihex = 1,8
          iel_new = iel_new + 1
          ! export to new re2 files
          ! write to vertice, midpoint, and then sideset.. to seperate files
          do inekvert = 1,27
             xm3(inekvert,1,1) = hexver2(1,inekvert,ihex)
             ym3(inekvert,1,1) = hexver2(2,inekvert,ihex)
             zm3(inekvert,1,1) = hexver2(3,inekvert,ihex)
          enddo

         !call map2reg(x3,3,xm3(1,1,1),1)  ! Map to 3x3x3 array
         !call map2reg(y3,3,ym3(1,1,1),1)
         !call map2reg(z3,3,zm3(1,1,1),1)

         !call map2reg_3di_e(x3,3,xm3(1,1,1),3)
         !call map2reg_3di_e(y3,3,ym3(1,1,1),3)
         !call map2reg_3di_e(z3,3,zm3(1,1,1),3)

         !call rzero(curve3,5*12)

         do iedge=1,nedge

            do j=1,3
               xyz(1,j) = xm3(e3(j,iedge),1,1)
               xyz(2,j) = ym3(e3(j,iedge),1,1)
               xyz(3,j) = zm3(e3(j,iedge),1,1)
            enddo

            call average2vec(mide(1), xyz(1,1), xyz(1,3))  

            call distance(xyz(1,1),xyz(1,3),dist0)
            call distance(mide(1),xyz(1,2),dist1)

            if (dist1.gt.(1e-3*dist0)) then

            ncurve = ncurve + 1
            
			endif

         enddo

       enddo  !        do ihex = 1,8

      enddo !       do iel_old=1,etot

      !ncurve = num_elem*nedge
      rcurve = ncurve
      call byte_write(rcurve,2, ierr)
      write(6,*) 'curves: ',ncurve

      iel_new = 0
      do iel_old=1,etot
       do ivert = 1,27
       hexver1(1,ivert) = xm1(ivert,1,1,iel_old)
       hexver1(2,ivert) = ym1(ivert,1,1,iel_old)
       hexver1(3,ivert) = zm1(ivert,1,1,iel_old)
       enddo
  
       do ifc = 1,6
       hexss1(ifc) = bc(5,ifc,iel_old)
       hexss3(ifc) = bc(4,ifc,iel_old) 
       enddo
 
       call hex_octuple_splitting(hexver1,hexver2,hexss1,hexss2,hexss3,hexss4) 

       do ihex = 1,8
          iel_new = iel_new + 1
          ! export to new re2 files
          ! write to vertice, midpoint, and then sideset.. to seperate files
          do inekvert = 1,27
             xm3(inekvert,1,1) = hexver2(1,inekvert,ihex)
             ym3(inekvert,1,1) = hexver2(2,inekvert,ihex)
             zm3(inekvert,1,1) = hexver2(3,inekvert,ihex)
          enddo

         !call map2reg(x3,3,xm3(1,1,1),1)  ! Map to 3x3x3 array
         !call map2reg(y3,3,ym3(1,1,1),1)
         !call map2reg(z3,3,zm3(1,1,1),1)

         !call map2reg_3di_e(x3,3,xm3(1,1,1),3)
         !call map2reg_3di_e(y3,3,ym3(1,1,1),3)
         !call map2reg_3di_e(z3,3,zm3(1,1,1),3)

         !call rzero(curve3,5*12)

         do iedge=1,nedge

            do j=1,3
               xyz(1,j) = xm3(e3(j,iedge),1,1)
               xyz(2,j) = ym3(e3(j,iedge),1,1)
               xyz(3,j) = zm3(e3(j,iedge),1,1)
            enddo

            call average2vec(mide(1), xyz(1,1), xyz(1,3))  

            call distance(xyz(1,1),xyz(1,3),dist0)
            call distance(mide(1),xyz(1,2),dist1)

            if (dist1.gt.(1e-3*dist0)) then

            curve3(1) = xm3(edge_mid(iedge),1,1)
            curve3(2) = ym3(edge_mid(iedge),1,1)
            curve3(3) = zm3(edge_mid(iedge),1,1)

            cc='m'
            buf2(1) = iel_new
            buf2(2) = iedge
            call copy       (buf2(3),curve3(1),5)
            call blank      (buf2(8),8)
            call chcopy     (buf2(8),cc,1)
            call byte_write (buf2,16,ierr)

			endif

         enddo

       enddo  !        do ihex = 1,8

      enddo !       do iel_old=1,etot

      endif !       if (iflag.eq.2) then  

! done writing edges	

!=============================
! write bc
      if (iflag.eq.3) then

      write(6,*) 'splitting and writing vbc'
      !call rzero8(buf2,30)

      nbc   = 0
      nface = 2*num_dim
      do iel=1,eftot
        do ifc=1,nface
          if (cbc(ifc,iel).ne.'   ')  nbc = nbc + 1
        enddo
      enddo
      nbc = nbc*4
      rbc = nbc
      call byte_write(rbc,2, ierr)
      write(6,*) 'velocity boundary faces: ',nbc

      iel_new = 0
      do iel_old=1,eftot
       do ivert = 1,27
       hexver1(1,ivert) = xm1(ivert,1,1,iel_old)
       hexver1(2,ivert) = ym1(ivert,1,1,iel_old)
       hexver1(3,ivert) = zm1(ivert,1,1,iel_old)
       enddo
  
       do ifc = 1,6
       hexss1(ifc) = bc(5,ifc,iel_old)
       hexss3(ifc) = bc(4,ifc,iel_old) 
       enddo
 
! 
       call hex_octuple_splitting(hexver1,hexver2,hexss1,hexss2,hexss3,hexss4) 

       do ihex = 1,8
          iel_new = iel_new + 1

          !call rzero(bc3(1,1),30)

          do ifc=1,6
             cbc3(ifc) = '   '
             if(hexss2(ifc,ihex).gt.0) then
              bc3(5,ifc) = hexss2(ifc,ihex)
              bc3(4,ifc) = hexss4(ifc,ihex)
              cbc3(ifc)  = 'EXO' ! dummy boundary condition
             endif
          enddo

          do ifc = 1,2*num_dim
          ch3 = cbc3(ifc) 
          if (ch3.ne.'   ') then
            buf2(1)=iel_new
            buf2(2)=ifc
            call copy   (buf2(3),bc3(1,ifc),5)
            call blank  (buf2(8),8)
            call chcopy (buf2(8),ch3,3)
            if (eftot2.ge.1000000) then
              ibc     = bc3(1,ifc)
              buf2(3) = ibc
            endif
            call byte_write (buf2,16,ierr)
          endif
        enddo


       enddo  !        do ihex = 1,8

      enddo !       do iel_old=1,etot


      if(etot.ne.eftot) then

      write(6,*) 'splitting and writing tbc'
      !call rzero8(buf2,30)

! writing thermal bc to all elements
      nbct   = 0
      nface = 2*num_dim
      do iel=1,etot
        do ifc=1,nface
          if (cbc(ifc,iel).ne.'   ')  nbct = nbct + 1
        enddo
      enddo
      nbct = nbct*4
      rbct = nbct
      call byte_write (rbct,2, ierr)
      write(6,*) 'thermal boundary faces: ',nbct

      iel_new = 0
      do iel_old=1,etot
       do ivert = 1,27
       hexver1(1,ivert) = xm1(ivert,1,1,iel_old)
       hexver1(2,ivert) = ym1(ivert,1,1,iel_old)
       hexver1(3,ivert) = zm1(ivert,1,1,iel_old)
       enddo
  
       do ifc = 1,6
       hexss1(ifc) = bc(5,ifc,iel_old)
       hexss3(ifc) = bc(4,ifc,iel_old) 
       enddo
 
! 
       call hex_octuple_splitting(hexver1,hexver2,hexss1,hexss2,hexss3,hexss4) 

       do ihex = 1,8
          iel_new = iel_new + 1
          !call rzero(bc3(1,1),30)
		  
          do ifc=1,6
             cbc3(ifc) = '   '
             if(hexss2(ifc,ihex).gt.0) then
              bc3(5,ifc) = hexss2(ifc,ihex)
              bc3(4,ifc) = hexss4(ifc,ihex)
              cbc3(ifc)   = 'EXO' ! dummy boundary condition
             endif
          enddo

          do ifc = 1,2*num_dim
          ch3 = cbc3(ifc) 
          if (ch3.ne.'   ') then
            buf2(1)=iel_new
            buf2(2)=ifc
            call copy   (buf2(3),bc3(1,ifc),5)
            call blank  (buf2(8),8)
            call chcopy (buf2(8),ch3,3)
            if (num_elem.ge.1000000) then
              ibc     = bc3(1,ifc)
              buf2(3) = ibc
            endif
            call byte_write (buf2,16,ierr)
          endif
        enddo


       enddo  !        do ihex = 1,8

      enddo !       do iel_old=1,etot

     endif  ! if(etot.ne.eftot) then

     endif  ! if (iflag.eq.3) then
! done writign bcs file
!==============================================================
      if (iflag.eq.0) then
      rcurve = 0.0
      call byte_write(rcurve,2, ierr)
      endif

      write(6,'(A)') 'Done: octuple splitting'

! deallocate orgininal mesh info

      write(6,*) 'closing re2 file'
      call byte_close (ierr)
      write(6,*) 'Done: closing re2 file'


      deallocate (xm1,ym1,zm1)
      deallocate (ccurve,curve)
      deallocate (cbc,bc)

      return
      end
!--------------------------------------------------------------------  

!-----------------------------------------------------------------
      subroutine octuple_split0()
!
! convert nek hex20 element into 8 hex20 elements.
!  
! 
      use SIZE
      use OCTUPLESIZE

      real*8 hexver1(3,27) ! hex coordinates
      real*8 hexver2(3,27,8) ! hex coordinate
      integer hexss1(6),hexss2(6,8)
      integer hexss3(6),hexss4(6,8)

      write(6,'(A)') 'Starting octuple splitting'

      etot_new = etot*8

      allocate ( xm2                (3,3,3,etot_new)      )
      allocate ( ym2                (3,3,3,etot_new)      )
      allocate ( zm2                (3,3,3,etot_new)      )

      call rzero(xm2,3*3*3*etot_new)
      call rzero(ym2,3*3*3*etot_new)
      call rzero(zm2,3*3*3*etot_new)

!      allocate   (ccurve2 (4+8*(num_dim-2),etot_new) )
!      allocate   (curve2  (2*num_dim,12,   etot_new) )
!      call rzero (curve2,2*num_dim*12*etot_new)
!      call blank (ccurve2,(4+8*(num_dim-2))*etot_new)

      allocate   (cbc2    (2*num_dim,      etot_new) )
      allocate   (bc2     (5,2*num_dim,    etot_new) )
      call rzero (bc2,5*2*num_dim*etot_new)
      call blank (cbc2,3*2*num_dim*etot_new)

      iel_new = 0 
      do iel_old=1,etot
! split  one hex20 into 8 hexes...

! pack ehexver, and ehexss
       do ivert = 1,27
       hexver1(1,ivert) = xm1(ivert,1,1,iel_old)
       hexver1(2,ivert) = ym1(ivert,1,1,iel_old)
       hexver1(3,ivert) = zm1(ivert,1,1,iel_old)
       enddo
  
       do ifc = 1,6
       hexss1(ifc) = bc(5,ifc,iel_old)
       hexss3(ifc) = bc(4,ifc,iel_old) 
       enddo
 
! quadratic splitting hex20 elements..
! 
       call hex_octuple_splitting(hexver1,hexver2,hexss1,hexss2,hexss3,hexss4) 

       do ihex = 1,8
          iel_new = iel_new + 1
          do inekvert = 1,27
             xm2(inekvert,1,1,iel_new) = hexver2(1,inekvert,ihex)
             ym2(inekvert,1,1,iel_new) = hexver2(2,inekvert,ihex)
             zm2(inekvert,1,1,iel_new) = hexver2(3,inekvert,ihex)
          enddo
          do ifc=1,6
             cbc2(ifc,iel_new) = '   '
             if(hexss2(ifc,ihex).gt.0) then
              bc2(5,ifc,iel_new) = hexss2(ifc,ihex)
              bc2(4,ifc,iel_new) = hexss4(ifc,ihex)
              cbc2(ifc,iel_new)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

      enddo

      etot = etot*8
      eftot = eftot*8
      num_elem = etot
      write(6,'(A)') 'Done: octuple splitting'


! deallocate orgininal mesh info

      deallocate (xm1,ym1,zm1)
      deallocate (ccurve,curve)
      deallocate (cbc,bc)

! reallocate use new tot
      allocate ( xm1                (3,3,3,etot_new)      )
      allocate ( ym1                (3,3,3,etot_new)      )
      allocate ( zm1                (3,3,3,etot_new)      )

      call rzero(xm1,3*3*3*etot_new)
      call rzero(ym1,3*3*3*etot_new)
      call rzero(zm1,3*3*3*etot_new)

      allocate   (ccurve (4+8*(num_dim-2),etot_new) )
      allocate   (curve  (2*num_dim,12,   etot_new) )
      call rzero (curve,2*num_dim*12*etot_new)
      call blank (ccurve,(4+8*(num_dim-2))*etot_new)

      allocate   (cbc    (2*num_dim,      etot_new) )
      allocate   (bc     (5,2*num_dim,    etot_new) )
      call rzero (bc,5*2*num_dim*etot_new)
      call blank (cbc,3*2*num_dim*etot_new)


! copy new mesh to previous mesh array

      do iel =1,etot

       do ivert = 1,27
       xm1(ivert,1,1,iel) = xm2(ivert,1,1,iel)
       ym1(ivert,1,1,iel) = ym2(ivert,1,1,iel)
       zm1(ivert,1,1,iel) = zm2(ivert,1,1,iel)
       enddo
       
	   do ifc=1,6
       bc(5,ifc,iel) = bc2(5,ifc,iel)
       bc(4,ifc,iel) = bc2(4,ifc,iel)
       cbc(ifc,iel) = cbc2(ifc,iel) 
       enddo

      enddo

      deallocate (xm2,ym2,zm2)
!      deallocate (ccurve2,curve2)
      deallocate (cbc2,bc2)

      return
      end
!--------------------------------------------------------------------  
      subroutine hex_octuple_splitting(hexver1,hexver2,hexss1,hexss2,hexss3,hexss4) 
!
! splitting one hex into 8 hexes
! but quadratically
      integer exo_to_nek_vert3D_new(20)                   ! hex27 to nek numbering
      data    exo_to_nek_vert3D_new                     &
            / 1,3,9,7,19,21,27,25,2,6,8,4,10,12,18,16,20,24,26,22/

      real*8 hexver1(3,27)
      real*8 hexver2(3,27,8) 
      integer hexss1(6),hexss2(6,8)
      integer hexss3(6),hexss4(6,8)

      real*8 ehexver(3,20)
      real*8 hexface(3,6)
      real*8 hexcen(3) ! tet vol center
      real*8 subhexver(3,20)

      real*8 quad(3,8),fc(3),ver(3),fc2(3),tvec1(3),tvec2(3)

! convert nek hex27 to exo hex20
      do ivert =1,20
         jvert = exo_to_nek_vert3D_new(ivert)
         call assignvec(ehexver(1,ivert),hexver1(1,jvert))
      enddo

      do i=1,6*8
      hexss2(i,1)=0
      hexss4(i,1)=0
      enddo

! get face center, but quadratically
      call assignvec(quad(1,1),ehexver(1,1))
      call assignvec(quad(1,2),ehexver(1,2))
      call assignvec(quad(1,3),ehexver(1,6))
      call assignvec(quad(1,4),ehexver(1,5))
      call assignvec(quad(1,5),ehexver(1,9))
      call assignvec(quad(1,6),ehexver(1,14))
      call assignvec(quad(1,7),ehexver(1,17))
      call assignvec(quad(1,8),ehexver(1,13))
      call transfinite_quad2(quad(1,1),fc(1))
	 
      call average4vec(tvec1,quad(1,1),quad(1,2),&
      quad(1,3),quad(1,4))
      call average4vec(tvec2,quad(1,5),quad(1,6),&
      quad(1,7),quad(1,8))
      call average2vec(fc2,tvec1,tvec2)  
      call average2vec(hexface(1,1),fc,fc2)
      !call assignvec(hexface(1,1),fc(1))

      call assignvec(quad(1,1),ehexver(1,2))
      call assignvec(quad(1,2),ehexver(1,3))
      call assignvec(quad(1,3),ehexver(1,7))
      call assignvec(quad(1,4),ehexver(1,6))
      call assignvec(quad(1,5),ehexver(1,10))
      call assignvec(quad(1,6),ehexver(1,15))
      call assignvec(quad(1,7),ehexver(1,18))
      call assignvec(quad(1,8),ehexver(1,14))
      call transfinite_quad2(quad(1,1),fc(1))
	  
      call average4vec(tvec1,quad(1,1),quad(1,2),&
      quad(1,3),quad(1,4))
      call average4vec(tvec2,quad(1,5),quad(1,6),&
      quad(1,7),quad(1,8))
      call average2vec(fc2,tvec1,tvec2)  
      call average2vec(hexface(1,2),fc,fc2)
      !call assignvec(hexface(1,2),fc(1))

      call assignvec(quad(1,1),ehexver(1,3))
      call assignvec(quad(1,2),ehexver(1,4))
      call assignvec(quad(1,3),ehexver(1,8))
      call assignvec(quad(1,4),ehexver(1,7))
      call assignvec(quad(1,5),ehexver(1,11))
      call assignvec(quad(1,6),ehexver(1,16))
      call assignvec(quad(1,7),ehexver(1,19))
      call assignvec(quad(1,8),ehexver(1,15))
      call transfinite_quad2(quad(1,1),fc(1))
      call average4vec(tvec1,quad(1,1),quad(1,2),&
      quad(1,3),quad(1,4))
      call average4vec(tvec2,quad(1,5),quad(1,6),&
      quad(1,7),quad(1,8))
      call average2vec(fc2,tvec1,tvec2)  
      call average2vec(hexface(1,3),fc,fc2)
	  !call assignvec(hexface(1,3),fc(1))

      call assignvec(quad(1,1),ehexver(1,4))
      call assignvec(quad(1,2),ehexver(1,1))
      call assignvec(quad(1,3),ehexver(1,5))
      call assignvec(quad(1,4),ehexver(1,8))
      call assignvec(quad(1,5),ehexver(1,12))
      call assignvec(quad(1,6),ehexver(1,13))
      call assignvec(quad(1,7),ehexver(1,20))
      call assignvec(quad(1,8),ehexver(1,16))
      call transfinite_quad2(quad(1,1),fc(1))
      call average4vec(tvec1,quad(1,1),quad(1,2),&
      quad(1,3),quad(1,4))
      call average4vec(tvec2,quad(1,5),quad(1,6),&
      quad(1,7),quad(1,8))
      call average2vec(fc2,tvec1,tvec2)  
      call average2vec(hexface(1,4),fc,fc2)
      !call assignvec(hexface(1,4),fc(1))

      call assignvec(quad(1,1),ehexver(1,1))
      call assignvec(quad(1,2),ehexver(1,2))
      call assignvec(quad(1,3),ehexver(1,3))
      call assignvec(quad(1,4),ehexver(1,4))
      call assignvec(quad(1,5),ehexver(1,9))
      call assignvec(quad(1,6),ehexver(1,10))
      call assignvec(quad(1,7),ehexver(1,11))
      call assignvec(quad(1,8),ehexver(1,12))
      call transfinite_quad2(quad(1,1),fc(1))
      call average4vec(tvec1,quad(1,1),quad(1,2),&
      quad(1,3),quad(1,4))
      call average4vec(tvec2,quad(1,5),quad(1,6),&
      quad(1,7),quad(1,8))
      call average2vec(fc2,tvec1,tvec2)  
      call average2vec(hexface(1,5),fc,fc2)
      !call assignvec(hexface(1,5),fc(1))

      call assignvec(quad(1,1),ehexver(1,5))
      call assignvec(quad(1,2),ehexver(1,6))
      call assignvec(quad(1,3),ehexver(1,7))
      call assignvec(quad(1,4),ehexver(1,8))
      call assignvec(quad(1,5),ehexver(1,17))
      call assignvec(quad(1,6),ehexver(1,18))
      call assignvec(quad(1,7),ehexver(1,19))
      call assignvec(quad(1,8),ehexver(1,20))
      call transfinite_quad2(quad(1,1),fc(1))
      call average4vec(tvec1,quad(1,1),quad(1,2),&
      quad(1,3),quad(1,4))
      call average4vec(tvec2,quad(1,5),quad(1,6),&
      quad(1,7),quad(1,8))
      call average2vec(fc2,tvec1,tvec2)  
      call average2vec(hexface(1,6),fc,fc2)
      !call assignvec(hexface(1,6),fc(1))

! get hex center
      call average6vec(hexcen(1),hexface(1,1),hexface(1,2),hexface(1,3),hexface(1,4),hexface(1,5),hexface(1,6))

!  assign coordinates to 8 hex.
!  hex 1
      hexss2(1,1) = hexss1(1)
      hexss2(4,1) = hexss1(4)
      hexss2(5,1) = hexss1(5)
      hexss4(1,1) = hexss3(1)
      hexss4(4,1) = hexss3(4)
      hexss4(5,1) = hexss3(5)

      call assignvec(subhexver(1,1),ehexver(1,1)) 
      call assignvec(subhexver(1,2),ehexver(1,9)) 
      call assignvec(subhexver(1,3),hexface(1,5)) 
      call assignvec(subhexver(1,4),ehexver(1,12))
      call assignvec(subhexver(1,5),ehexver(1,13)) 
      call assignvec(subhexver(1,6),hexface(1,1)) 
      call assignvec(subhexver(1,7),hexcen(1))
      call assignvec(subhexver(1,8),hexface(1,4))

      call quick_quater_point_from_curve(ver(1),ehexver(1,1),ehexver(1,9),ehexver(1,2))
      call assignvec(subhexver(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,9),hexface(1,5),ehexver(1,11))
      call assignvec(subhexver(1,10),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,12),hexface(1,5),ehexver(1,10))
      call assignvec(subhexver(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,1),ehexver(1,12),ehexver(1,4))
      call assignvec(subhexver(1,12),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,1),ehexver(1,13),ehexver(1,5))
      call assignvec(subhexver(1,13),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,9),hexface(1,1),ehexver(1,17))
      call assignvec(subhexver(1,14),ver(1))
      call quick_quater_point_from_curve(ver(1),hexface(1,5),hexcen(1),hexface(1,6))
      call assignvec(subhexver(1,15),ver(1)) 
     ! call quick_quater_point_from_curve(ver(1),ehexver(1,2),hexface(1,4),ehexver(1,20)) ! error !
      call quick_quater_point_from_curve(ver(1),ehexver(1,12),hexface(1,4),ehexver(1,20))

      call assignvec(subhexver(1,16),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,13),hexface(1,1),ehexver(1,14))
      call assignvec(subhexver(1,17),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,1),hexcen(1),hexface(1,3))
      call assignvec(subhexver(1,18),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,4),hexcen(1),hexface(1,2))
      call assignvec(subhexver(1,19),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,13),hexface(1,4),ehexver(1,16))
      call assignvec(subhexver(1,20),ver(1)) 

      call hex20tohex27_2(hexver2(1,1,1),subhexver(1,1))	 

!  hex 2
      hexss2(1,2) = hexss1(1)
      hexss2(2,2) = hexss1(2)
      hexss2(5,2) = hexss1(5)
      hexss4(1,2) = hexss3(1)
      hexss4(2,2) = hexss3(2)
      hexss4(5,2) = hexss3(5)

      call assignvec(subhexver(1,1),ehexver(1,9)) 
      call assignvec(subhexver(1,2),ehexver(1,2)) 
      call assignvec(subhexver(1,3),ehexver(1,10)) 
      call assignvec(subhexver(1,4),hexface(1,5))
      call assignvec(subhexver(1,5),hexface(1,1)) 
      call assignvec(subhexver(1,6),ehexver(1,14)) 
      call assignvec(subhexver(1,7),hexface(1,2)) 
      call assignvec(subhexver(1,8),hexcen(1))

      call quick_quater_point_from_curve(ver(1),ehexver(1,2),ehexver(1,9),ehexver(1,1))
      call assignvec(subhexver(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,2),ehexver(1,10),ehexver(1,3))
      call assignvec(subhexver(1,10),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,10),hexface(1,5),ehexver(1,12))
      call assignvec(subhexver(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,9),hexface(1,5),ehexver(1,11))
      call assignvec(subhexver(1,12),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,9),hexface(1,1),ehexver(1,17))
      call assignvec(subhexver(1,13),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,2),ehexver(1,14),ehexver(1,6))
      call assignvec(subhexver(1,14),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,10),hexface(1,2),ehexver(1,18))
      call assignvec(subhexver(1,15),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,5),hexcen(1),hexface(1,6))
      call assignvec(subhexver(1,16),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,14),hexface(1,1),ehexver(1,13))
      call assignvec(subhexver(1,17),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,14),hexface(1,2),ehexver(1,15))
      call assignvec(subhexver(1,18),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,2),hexcen(1),hexface(1,4))
      call assignvec(subhexver(1,19),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,1),hexcen(1),hexface(1,3))
      call assignvec(subhexver(1,20),ver(1)) 

      call hex20tohex27_2(hexver2(1,1,2),subhexver(1,1))	

!  hex 3
      hexss2(2,3) = hexss1(2)
      hexss2(3,3) = hexss1(3)
      hexss2(5,3) = hexss1(5)
      hexss4(2,3) = hexss3(2)
      hexss4(3,3) = hexss3(3)
      hexss4(5,3) = hexss3(5)

      call assignvec(subhexver(1,1),hexface(1,5)) 
      call assignvec(subhexver(1,2),ehexver(1,10)) 
      call assignvec(subhexver(1,3),ehexver(1,3)) 
      call assignvec(subhexver(1,4),ehexver(1,11))
      call assignvec(subhexver(1,5),hexcen(1)) 
      call assignvec(subhexver(1,6),hexface(1,2)) 
      call assignvec(subhexver(1,7),ehexver(1,15))
      call assignvec(subhexver(1,8),hexface(1,3)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,10),hexface(1,5),ehexver(1,12))
      call assignvec(subhexver(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,3),ehexver(1,10),ehexver(1,2))
      call assignvec(subhexver(1,10),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,3),ehexver(1,11),ehexver(1,4))
      call assignvec(subhexver(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,11),hexface(1,5),ehexver(1,9))
      call assignvec(subhexver(1,12),ver(1)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,5),hexcen(1),hexface(1,6))
      call assignvec(subhexver(1,13),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,10),hexface(1,2),ehexver(1,18))
      call assignvec(subhexver(1,14),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,3),ehexver(1,15),ehexver(1,7))
      call assignvec(subhexver(1,15),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,11),hexface(1,3),ehexver(1,19))
      call assignvec(subhexver(1,16),ver(1)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,2),hexcen(1),hexface(1,4))
      call assignvec(subhexver(1,17),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,15),hexface(1,2),ehexver(1,14))
      call assignvec(subhexver(1,18),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,15),hexface(1,3),ehexver(1,16))
      call assignvec(subhexver(1,19),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,3),hexcen(1),hexface(1,1))
      call assignvec(subhexver(1,20),ver(1)) 

      call hex20tohex27_2(hexver2(1,1,3),subhexver(1,1))	

!  hex 4

      hexss2(3,4) = hexss1(3)
      hexss2(4,4) = hexss1(4)
      hexss2(5,4) = hexss1(5)
      hexss4(3,4) = hexss3(3)
      hexss4(4,4) = hexss3(4)
      hexss4(5,4) = hexss3(5)

      call assignvec(subhexver(1,1),ehexver(1,12)) 
      call assignvec(subhexver(1,2),hexface(1,5)) 
      call assignvec(subhexver(1,3),ehexver(1,11)) 
      call assignvec(subhexver(1,4),ehexver(1,4))
      call assignvec(subhexver(1,5),hexface(1,4)) 
      call assignvec(subhexver(1,6),hexcen(1)) 
      call assignvec(subhexver(1,7),hexface(1,3)) 
      call assignvec(subhexver(1,8),ehexver(1,16)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,12),hexface(1,5),ehexver(1,10))
      call assignvec(subhexver(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,11),hexface(1,5),ehexver(1,9))
      call assignvec(subhexver(1,10),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,4),ehexver(1,11),ehexver(1,3))
      call assignvec(subhexver(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,4),ehexver(1,12),ehexver(1,1))
      call assignvec(subhexver(1,12),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,12),hexface(1,4),ehexver(1,20))
      call assignvec(subhexver(1,13),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,5),hexcen(1),hexface(1,6))
      call assignvec(subhexver(1,14),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,11),hexface(1,3),ehexver(1,19))
      call assignvec(subhexver(1,15),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,4),ehexver(1,16),ehexver(1,8))
      call assignvec(subhexver(1,16),ver(1)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,4),hexcen(1),hexface(1,2))
      call assignvec(subhexver(1,17),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,3),hexcen(1),hexface(1,1))
      call assignvec(subhexver(1,18),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,16),hexface(1,3),ehexver(1,15))
      call assignvec(subhexver(1,19),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,16),hexface(1,4),ehexver(1,13))
      call assignvec(subhexver(1,20),ver(1)) 

      call hex20tohex27_2(hexver2(1,1,4),subhexver(1,1))	

!  hex 5
      hexss2(1,5) = hexss1(1)
      hexss2(4,5) = hexss1(4)
      hexss2(6,5) = hexss1(6)
      hexss4(1,5) = hexss3(1)
      hexss4(4,5) = hexss3(4)
      hexss4(6,5) = hexss3(6)

      call assignvec(subhexver(1,1),ehexver(1,13)) 
      call assignvec(subhexver(1,2),hexface(1,1)) 
      call assignvec(subhexver(1,3),hexcen(1)) 
      call assignvec(subhexver(1,4),hexface(1,4))
      call assignvec(subhexver(1,5),ehexver(1,5)) 
      call assignvec(subhexver(1,6),ehexver(1,17)) 
      call assignvec(subhexver(1,7),hexface(1,6)) 
      call assignvec(subhexver(1,8),ehexver(1,20)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,13),hexface(1,1),ehexver(1,14))
      call assignvec(subhexver(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,1),hexcen(1),hexface(1,3))
      call assignvec(subhexver(1,10),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,4),hexcen(1),hexface(1,2))
      call assignvec(subhexver(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,13),hexface(1,4),ehexver(1,16))
      call assignvec(subhexver(1,12),ver(1))

      call quick_quater_point_from_curve(ver(1),ehexver(1,5),ehexver(1,13),ehexver(1,1))
      call assignvec(subhexver(1,13),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,17),hexface(1,1),ehexver(1,9))
      call assignvec(subhexver(1,14),ver(1))
      call quick_quater_point_from_curve(ver(1),hexface(1,6),hexcen(1),hexface(1,5))
      call assignvec(subhexver(1,15),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,20),hexface(1,4),ehexver(1,12))
      call assignvec(subhexver(1,16),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,5),ehexver(1,17),ehexver(1,6))
      call assignvec(subhexver(1,17),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,17),hexface(1,6),ehexver(1,19))
      call assignvec(subhexver(1,18),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,20),hexface(1,6),ehexver(1,18))
      call assignvec(subhexver(1,19),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,5),ehexver(1,20),ehexver(1,8))
      call assignvec(subhexver(1,20),ver(1)) 

      call hex20tohex27_2(hexver2(1,1,5),subhexver(1,1))	

!  hex 6
      hexss2(1,6) = hexss1(1)
      hexss2(2,6) = hexss1(2)
      hexss2(6,6) = hexss1(6)
      hexss4(1,6) = hexss3(1)
      hexss4(2,6) = hexss3(2)
      hexss4(6,6) = hexss3(6)

      call assignvec(subhexver(1,1),hexface(1,1)) 
      call assignvec(subhexver(1,2),ehexver(1,14)) 
      call assignvec(subhexver(1,3),hexface(1,2)) 
      call assignvec(subhexver(1,4),hexcen(1))
      call assignvec(subhexver(1,5),ehexver(1,17)) 
      call assignvec(subhexver(1,6),ehexver(1,6)) 
      call assignvec(subhexver(1,7),ehexver(1,18)) 
      call assignvec(subhexver(1,8),hexface(1,6)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,14),hexface(1,1),ehexver(1,13))
      call assignvec(subhexver(1,9),ver(1)) 
      !call quick_quater_point_from_curve(ver(1),ehexver(1,14),hexface(1,1),ehexver(1,15)) !! error here!
      call quick_quater_point_from_curve(ver(1),ehexver(1,14),hexface(1,2),ehexver(1,15))  ! this is the correct version!
      call assignvec(subhexver(1,10),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,2),hexcen(1),hexface(1,4))
      call assignvec(subhexver(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,1),hexcen(1),hexface(1,3))
      call assignvec(subhexver(1,12),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,17),hexface(1,1),ehexver(1,9))
      call assignvec(subhexver(1,13),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,6),ehexver(1,14),ehexver(1,2))
      call assignvec(subhexver(1,14),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,18),hexface(1,2),ehexver(1,10))
      call assignvec(subhexver(1,15),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,6),hexcen(1),hexface(1,5))
      call assignvec(subhexver(1,16),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,6),ehexver(1,17),ehexver(1,5))
      call assignvec(subhexver(1,17),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,6),ehexver(1,18),ehexver(1,7))
      call assignvec(subhexver(1,18),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,18),hexface(1,6),ehexver(1,20))
      call assignvec(subhexver(1,19),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,17),hexface(1,6),ehexver(1,19))
      call assignvec(subhexver(1,20),ver(1)) 

      call hex20tohex27_2(hexver2(1,1,6),subhexver(1,1))	

!  hex 7
      hexss2(2,7) = hexss1(2)
      hexss2(3,7) = hexss1(3)
      hexss2(6,7) = hexss1(6)
      hexss4(2,7) = hexss3(2)
      hexss4(3,7) = hexss3(3)
      hexss4(6,7) = hexss3(6)

      call assignvec(subhexver(1,1),hexcen(1)) 
      call assignvec(subhexver(1,2),hexface(1,2))
      call assignvec(subhexver(1,3),ehexver(1,15)) 
      call assignvec(subhexver(1,4),hexface(1,3))
      call assignvec(subhexver(1,5),hexface(1,6))
      call assignvec(subhexver(1,6),ehexver(1,18)) 
      call assignvec(subhexver(1,7),ehexver(1,7))
      call assignvec(subhexver(1,8),ehexver(1,19)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,2),hexcen(1),hexface(1,4))
      call assignvec(subhexver(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,15),hexface(1,2),ehexver(1,14))
      call assignvec(subhexver(1,10),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,15),hexface(1,3),ehexver(1,16))
      call assignvec(subhexver(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,3),hexcen(1),hexface(1,1))
      call assignvec(subhexver(1,12),ver(1)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,6),hexcen(1),hexface(1,5))
      call assignvec(subhexver(1,13),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,18),hexface(1,2),ehexver(1,10))
      call assignvec(subhexver(1,14),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,7),ehexver(1,15),ehexver(1,3))
      call assignvec(subhexver(1,15),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,19),hexface(1,3),ehexver(1,11))
      call assignvec(subhexver(1,16),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,18),hexface(1,6),ehexver(1,20))
      call assignvec(subhexver(1,17),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,7),ehexver(1,18),ehexver(1,6))
      call assignvec(subhexver(1,18),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,7),ehexver(1,19),ehexver(1,8))
      call assignvec(subhexver(1,19),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,19),hexface(1,6),ehexver(1,17))
      call assignvec(subhexver(1,20),ver(1)) 

      call hex20tohex27_2(hexver2(1,1,7),subhexver(1,1))	

!  hex 8
      hexss2(3,8) = hexss1(3)
      hexss2(4,8) = hexss1(4)
      hexss2(6,8) = hexss1(6)
      hexss4(3,8) = hexss3(3)
      hexss4(4,8) = hexss3(4)
      hexss4(6,8) = hexss3(6)


      call assignvec(subhexver(1,1),hexface(1,4))
      call assignvec(subhexver(1,2),hexcen(1))
      call assignvec(subhexver(1,3),hexface(1,3)) 
      call assignvec(subhexver(1,4),ehexver(1,16)) 
      call assignvec(subhexver(1,5),ehexver(1,20))
      call assignvec(subhexver(1,6),hexface(1,6))
      call assignvec(subhexver(1,7),ehexver(1,19))
      call assignvec(subhexver(1,8),ehexver(1,8)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,4),hexcen(1),hexface(1,2))
      call assignvec(subhexver(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,3),hexcen(1),hexface(1,1))
      call assignvec(subhexver(1,10),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,16),hexface(1,3),ehexver(1,15))
      call assignvec(subhexver(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,16),hexface(1,4),ehexver(1,13))
      call assignvec(subhexver(1,12),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,20),hexface(1,4),ehexver(1,12))
      call assignvec(subhexver(1,13),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,6),hexcen(1),hexface(1,5))
      call assignvec(subhexver(1,14),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,19),hexface(1,3),ehexver(1,11))
      call assignvec(subhexver(1,15),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,8),ehexver(1,16),ehexver(1,4))
      call assignvec(subhexver(1,16),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,20),hexface(1,6),ehexver(1,18))
      call assignvec(subhexver(1,17),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,19),hexface(1,6),ehexver(1,17))
      call assignvec(subhexver(1,18),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,8),ehexver(1,19),ehexver(1,7))
      call assignvec(subhexver(1,19),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,8),ehexver(1,20),ehexver(1,5))
      call assignvec(subhexver(1,20),ver(1)) 

      call hex20tohex27_2(hexver2(1,1,8),subhexver(1,1))	

      return
      end
!--------------------------------------------------------------------  
!--------------------------------------------------------------------
      subroutine hex20tohex27_2(hex27,hex20)
!  convert hex20 coordinates exo order to hex27 coordinates in nek.

      integer exo_to_nek_vert3D_new(20)                   ! hex27 to nek numbering
      data    exo_to_nek_vert3D_new                     &
            / 1,3,9,7,19,21,27,25,2,6,8,4,10,12,18,16,20,24,26,22/

      real*8 hex20(3,20)
      real*8 hex27(3,27)

! convert nek hex27 to exo hex20
      do ivert =1,20
         jvert = exo_to_nek_vert3D_new(ivert)
         call assignvec(hex27(1,jvert),hex20(1,ivert))
      enddo

      return
      end
! -------------------------------------------------------------------
! -------------------------------------------------------------------
      subroutine average6vec(avg,a,b,c,d,e,f)
      real*8  a(3),b(3),c(3),d(3),e(3),f(3),avg(3)       
      avg(1) = (a(1)+b(1)+c(1)+d(1)+e(1)+f(1))/6.0
      avg(2) = (a(2)+b(2)+c(2)+d(2)+e(2)+f(2))/6.0
      avg(3) = (a(3)+b(3)+c(3)+d(3)+e(3)+f(3))/6.0
      return
      end
! -------------------------------------------------------------------
      subroutine quick_quater_point_from_curve(v4,v1,v2,v3)
      real*8 v1(3),v2(3),v3(3),v4(3)
      real*8 x(3),cr1(3,3),beta

      x(1) = 0.0
      x(2) = 0.5
      x(3) = 1.0

      call assignvec(cr1(1,1),v1(1))
      call assignvec(cr1(1,2),v2(1))
      call assignvec(cr1(1,3),v3(1))

      beta = 0.25
      call curve(cr1,x,beta,v4(1))

      !call average2vec(v4(1),v1(1),v2(1))  

      return
      end
! -------------------------------------------------------------------

