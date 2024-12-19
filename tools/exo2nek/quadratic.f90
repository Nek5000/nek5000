!-------------------------------------------------
      subroutine split_convert1_quadratic()
! only accepting tet10+wed15 mesh
! !
! 
! quadratic_option = 1:
! pre-check element, for wedge elements only
! if non-right-hand elements are generated.
! adjust wedge nodes on the triangle surface that are with sideset>0
!
! quadratic_option = 2:
! pre-check element, for tet/wedge/hex
! if non-right-hand elements are generated.
! adjust mid-edge nodes
!
! quadratic_option = 3:
! actual splitting elements
! 
 
      use SIZE
      include 'exodusII.inc'

      integer hex8_to_hex27_vertex(8)
      data hex8_to_hex27_vertex /1,3,9,7,19,21,27,25/ ! for nek non-right-hand element check

!  in fortran arrays higher order is at back.
      real*8 hexver(3,27,8) ! hex coordinates
      real*8 tetver(3,10),wedgever(3,15),wedgever2(3,15)
      real*8 ehexver(3,20)
      real*8 vert1(3),vert2(3)
      real*8 wedThickness,dist

!  exoss is used to store sideset information for all exo elements.
!  tet,hex,wedge
      integer exoss(6,num_elem) 

      integer tetss(4),wedgess(5),ehexss(6),hexss(6,8)
      integer*8 ehexnumber,tetnumber,wedgenumber,bctot
      integer*8 vert_index_exo,ivert_index_exo
      integer*8 iel_nek,iel_exo,ifc_exo,iel_exo_g
	  
      real*8 XYZorg(3,8),nv(3),delta_offset

	  integer iblk,iblk1,icorrection,n
      integer iss

      logical hasnrh,ifnrh_in_wedge,ifnonrighthand

      save ehexnumber,tetnumber,wedgenumber,bctot
      save vert_index_exo
      save iel_nek,iel_exo,ifc_exo,iel_exo_g
      eacc_old = eacc
      
      call rzero_int(exoss,6*num_elem)
! store sideset information
      if (num_side_sets.ne.0) then
      write(6,'(a)') ''
      write(6,'(a)') 'Store SideSet information from EXO file'
        do iss=1,num_side_sets   ! loop over ss 
        !write(6,'(a,i2,a)') 'Sideset ',idss(iss), ' ...'
        write(6,*) 'Sideset ',idss(iss), ' ...'
           do i=1,num_sides_in_set(iss) ! loop over sides in ss
             iel_exo = elem_list(i,iss) ! exo element number
             ifc_exo = side_list(i,iss) ! exo element face number
             exoss(ifc_exo,iel_exo) = idss(iss)
           enddo
        enddo
      endif

      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '

      !call rzero_int(hexss,6*8) 
      ! write(6,*) 'flag1'
      vert_index_exo = 0
      iel_nek =  eacc_old !0
	  ! write(6,*) 'flag2'

      do iblk = 1,num_elem_blk ! deal with multiple blocks.
	  
        write(6,*) 'Converting elements in block ',idblk(iblk)
        ! detect element type of this block
        nvert=num_nodes_per_elem(iblk)
        write(6,*) 'nvert, ', nvert
		

        if (nvert.eq.10) then ! splitting tet10
		
          do iel_exo = 1,num_elem_in_block(iblk)
           !write(6,*) 'iel_exo, ',iel_exo

         iel_exo_g = iel_exo ! iel_exo is element number in this block
		 !   iel_exo_g is global exo element number
         if (iblk.gt.1) then
            do iblk1 = 1,iblk-1
             iel_exo_g = iel_exo_g + num_elem_in_block(iblk1)
            enddo
         endif
		 
          if (quadratic_option.eq.1) then
! only read tet10, but doing nothing
       do ivert = 1,10
       vert_index_exo = vert_index_exo + 1  ! advance vertice index
       tetver(1,ivert) = x_exo(connect(vert_index_exo))
       tetver(2,ivert) = y_exo(connect(vert_index_exo))
       tetver(3,ivert) = z_exo(connect(vert_index_exo))
       enddo 
          else if  (quadratic_option.eq.2) then

! linearize tet elements if non-right hand elements are created
       n = 4
       call rzero_int2(tetss,n) 
       if (num_side_sets.ne.0) then
        do ifc_exo=1,4
         tetss(ifc_exo) = 0
         tetss(ifc_exo) = exoss(ifc_exo,iel_exo_g)
        enddo
       endif

       icorrection = 0

       do while (.True.)
       hasnrh = .false.
	   
       do ivert = 1, 10
       ivert_index_exo = vert_index_exo + ivert  
       tetver(1,ivert) = x_exo(connect(ivert_index_exo))
       tetver(2,ivert) = y_exo(connect(ivert_index_exo))
       tetver(3,ivert) = z_exo(connect(ivert_index_exo))
       enddo

       call tettohex_quadratic(hexver,tetver,hexss,tetss,hasnrh)
	   
!       do i = 1,4
!       iel_nek = iel_nek + 1
!       call if_in_nj_list(iel_nek,ifnj)
!       hasnrh = (hasnrh.or.ifnj) ! 
!       enddo

       if (hasnrh) then
	
          icorrection = icorrection + 1
          if (icorrection.gt.10) then
           write(6,*) 'stuck in non-right-hand correction of tet10'
           write(6,*) 'at ', tetver(1,1),tetver(2,1),tetver(3,1)
           exit
          endif
	   
          ! adjust mid-edge node position if non-right-hand elements detected
          call average2vec(vert1(1),tetver(1,1),tetver(1,2))
          call average2vec(vert2(1),vert1(1),tetver(1,5))
          call assignvec(tetver(1,5),vert2(1))
		  
          call average2vec(vert1(1),tetver(1,2),tetver(1,3))
          call average2vec(vert2(1),vert1(1),tetver(1,6))
          call assignvec(tetver(1,6),vert2(1))
		  
          call average2vec(vert1(1),tetver(1,1),tetver(1,3))
          call average2vec(vert2(1),vert1(1),tetver(1,7))
          call assignvec(tetver(1,7),vert2(1))
		  
          call average2vec(vert1(1),tetver(1,1),tetver(1,4))
          call average2vec(vert2(1),vert1(1),tetver(1,8))
          call assignvec(tetver(1,8),vert2(1))
		  
          call average2vec(vert1(1),tetver(1,2),tetver(1,4))
          call average2vec(vert2(1),vert1(1),tetver(1,9))
          call assignvec(tetver(1,9),vert2(1))
		  
          call average2vec(vert1(1),tetver(1,3),tetver(1,4))
          call average2vec(vert2(1),vert1(1),tetver(1,10))
          call assignvec(tetver(1,10),vert2(1))

!          call average2vec(tetver(1,6),tetver(1,2),tetver(1,3))
!          call average2vec(tetver(1,7),tetver(1,1),tetver(1,3))
!          call average2vec(tetver(1,8),tetver(1,1),tetver(1,4))
!          call average2vec(tetver(1,9),tetver(1,2),tetver(1,4))
!          call average2vec(tetver(1,10),tetver(1,3),tetver(1,4))

          do ivert = 1, 10
          ivert_index_exo = vert_index_exo + ivert  
          x_exo(connect(ivert_index_exo)) = tetver(1,ivert) ! return adjusted tet
          y_exo(connect(ivert_index_exo)) = tetver(2,ivert)
          z_exo(connect(ivert_index_exo)) = tetver(3,ivert)
          enddo

       else
        exit    ! break if no non-right-hand elements
       endif

      enddo
	  
       vert_index_exo = vert_index_exo + 10


          else if (quadratic_option.eq.3) then
! actual splitting.
       n = 4
       call rzero_int2(tetss,n)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,4
         tetss(ifc_exo) = 0
         tetss(ifc_exo) = exoss(ifc_exo,iel_exo_g)
        enddo
       endif

!  read tet 10
       do ivert = 1, 10
       vert_index_exo = vert_index_exo + 1  
       tetver(1,ivert) = x_exo(connect(vert_index_exo))
       tetver(2,ivert) = y_exo(connect(vert_index_exo))
       tetver(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

!  given tet10 vertices, and return you four hex27 coords.
!       call tettohex(hexver,tetver,hexss,tetss) 

! quadratic tet-to-hex conversion

       call tettohex_quadratic(hexver,tetver,hexss,tetss,hasnrh)

       do ihex = 1, 4
          iel_nek = iel_nek + 1
          eacc = eacc + 1
          !write(6,*)  iel_nek
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)& 
       + shiftvector(1)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)&
       + shiftvector(2)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)&
       + shiftvector(3)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              bc(4,ifc,iel_nek) = iexo
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo		
	   
          endif

          enddo

        else if (nvert.eq.15) then ! wedge15
          do iel_exo = 1,num_elem_in_block(iblk)

         iel_exo_g = iel_exo
         if (iblk.gt.1) then
            do iblk1 = 1,iblk-1
             iel_exo_g = iel_exo_g + num_elem_in_block(iblk1)
            enddo
         endif

          if (quadratic_option.eq.1) then

! fix invalid wedge element if it produce non-right-hand hex elements.
! 
       icorrection  = 0
       do while (.True.)

       do ivert = 1,6
       ivert_index_exo = vert_index_exo + ivert  
       wedgever(1,ivert) = x_exo(connect(ivert_index_exo))
       wedgever(2,ivert) = y_exo(connect(ivert_index_exo))
       wedgever(3,ivert) = z_exo(connect(ivert_index_exo))
       enddo

       do ivert = 7,15
       ivert_index_exo = vert_index_exo + ivert  
       wedgever2(1,ivert) = x_exo(connect(ivert_index_exo))
       wedgever2(2,ivert) = y_exo(connect(ivert_index_exo))
       wedgever2(3,ivert) = z_exo(connect(ivert_index_exo))
       enddo

       call average2vec(wedgever(1,7),wedgever(1,1), &
       wedgever(1,2))
       call average2vec(wedgever(1,8),wedgever(1,2), &
       wedgever(1,3))
       call average2vec(wedgever(1,9),wedgever(1,1), &
       wedgever(1,3))
       call average2vec(wedgever(1,10),wedgever(1,1), &
       wedgever(1,4))
       call average2vec(wedgever(1,11),wedgever(1,2), &
       wedgever(1,5))
       call average2vec(wedgever(1,12),wedgever(1,3), &
      wedgever(1,6))
       call average2vec(wedgever(1,13),wedgever(1,4), &
      wedgever(1,5))
       call average2vec(wedgever(1,14),wedgever(1,5), &
      wedgever(1,6))
       call average2vec(wedgever(1,15),wedgever(1,4), &
      wedgever(1,6))


!  assign sideset to wedge element
       n = 5
       call rzero_int2(wedgess,n)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,5
          wedgess(ifc_exo) = exoss(ifc_exo,iel_exo_g)
        enddo
       endif


!  given wedge15 vertices, and return you 3 hex coords.
! this is done in linear way,
! just to check if non-right-hand..
       call wedgetohex(hexver,wedgever,hexss,wedgess)

! evaluate non-right-hand ness
! if there is non-rigth-hand elements, then explicit fix this wedge elements

       ifnrh_in_wedge = .False.
       ifnonrighthand = .False.
       do ihex = 1,3
	     do ivert= 1,8
         call assignvec(XYZorg(1,ivert),hexver(1,hex8_to_hex27_vertex(ivert),ihex))
         enddo
         call nek_check_non_right_hand_per_element(XYZorg,ifnonrighthand)

         ifnrh_in_wedge = ( ifnrh_in_wedge.or.ifnonrighthand)
       enddo

       if (ifnrh_in_wedge) then
! fix fix this wedge element
         write(6,*) 'found non-right hand elements from wedge-to-hex '
         write(6,*) 'will fix this wedge element', iel_exo_g
		 
          icorrection = icorrection + 1
          if (icorrection.gt.10) then
           write(6,*) 'stuck in wedge fix'
           write(6,*) 'at ', wedgever(1,1),wedgever(2,1),wedgever(3,1)
           exit
          endif

!          delta_offset = 1e-4 ! adjustment thickness
! get wedge element thickness.
          wedThickness = 0.0
          do ivert = 1,3
          call distance(wedgever(1,ivert),wedgever(1,ivert+3),dist)
	      wedThickness = wedThickness + dist/3.0
          enddo
          delta_offset = wedThickness*0.01

! expand this wedge with delta_offset
! 1st layer elements are moved towards outside more, and inwards slightly.
! non-1st layer element, are moved in both directions. 

          if (wedgess(4).gt.0) then

          call cross_prod(nv,wedgever(1,1),wedgever(1,3),wedgever(1,2))

	      do ivert = 1,3
          wedgever(1,ivert) = wedgever(1,ivert) + nv(1)*delta_offset
          wedgever(2,ivert) = wedgever(2,ivert) + nv(2)*delta_offset
          wedgever(3,ivert) = wedgever(3,ivert) + nv(3)*delta_offset
          enddo
		  
	      do ivert = 7,9
          wedgever2(1,ivert) = wedgever2(1,ivert) + nv(1)*delta_offset
          wedgever2(2,ivert) = wedgever2(2,ivert) + nv(2)*delta_offset
          wedgever2(3,ivert) = wedgever2(3,ivert) + nv(3)*delta_offset
          enddo

	      do ivert = 4,6
          wedgever(1,ivert) = wedgever(1,ivert) - nv(1)*delta_offset*0.01
          wedgever(2,ivert) = wedgever(2,ivert) - nv(2)*delta_offset*0.01
          wedgever(3,ivert) = wedgever(3,ivert) - nv(3)*delta_offset*0.01
          enddo
		  
	      do ivert = 13,15
          wedgever2(1,ivert) = wedgever2(1,ivert) - nv(1)*delta_offset*0.01
          wedgever2(2,ivert) = wedgever2(2,ivert) - nv(2)*delta_offset*0.01
          wedgever2(3,ivert) = wedgever2(3,ivert) - nv(3)*delta_offset*0.01
          enddo

          elseif (wedgess(5).gt.0)  then

          call cross_prod(nv,wedgever(1,4),wedgever(1,5),wedgever(1,6))

	      do ivert = 4,6
          wedgever(1,ivert) = wedgever(1,ivert) + nv(1)*delta_offset
          wedgever(2,ivert) = wedgever(2,ivert) + nv(2)*delta_offset
          wedgever(3,ivert) = wedgever(3,ivert) + nv(3)*delta_offset
          enddo

	      do ivert = 13,15
          wedgever2(1,ivert) = wedgever2(1,ivert) + nv(1)*delta_offset
          wedgever2(2,ivert) = wedgever2(2,ivert) + nv(2)*delta_offset
          wedgever2(3,ivert) = wedgever2(3,ivert) + nv(3)*delta_offset
          enddo
		  
	      do ivert = 1,3
          wedgever(1,ivert) = wedgever(1,ivert) - nv(1)*delta_offset*0.01
          wedgever(2,ivert) = wedgever(2,ivert) - nv(2)*delta_offset*0.01
          wedgever(3,ivert) = wedgever(3,ivert) - nv(3)*delta_offset*0.01
          enddo
		  
	      do ivert = 7,9
          wedgever2(1,ivert) = wedgever2(1,ivert) - nv(1)*delta_offset*0.01
          wedgever2(2,ivert) = wedgever2(2,ivert) - nv(2)*delta_offset*0.01
          wedgever2(3,ivert) = wedgever2(3,ivert) - nv(3)*delta_offset*0.01
          enddo

          else  ! if non wedgess(4) nor wedgess(5) > 0, expand in both direction slightly..


          call cross_prod(nv,wedgever(1,1),wedgever(1,3),wedgever(1,2))

	      do ivert = 1,3
          wedgever(1,ivert) = wedgever(1,ivert) + nv(1)*delta_offset*0.01
          wedgever(2,ivert) = wedgever(2,ivert) + nv(2)*delta_offset*0.01
          wedgever(3,ivert) = wedgever(3,ivert) + nv(3)*delta_offset*0.01
          enddo
		  
	      do ivert = 7,9
          wedgever2(1,ivert) = wedgever2(1,ivert) + nv(1)*delta_offset*0.01
          wedgever2(2,ivert) = wedgever2(2,ivert) + nv(2)*delta_offset*0.01
          wedgever2(3,ivert) = wedgever2(3,ivert) + nv(3)*delta_offset*0.01
          enddo


          call cross_prod(nv,wedgever(1,4),wedgever(1,5),wedgever(1,6))

	      do ivert = 4,6
          wedgever(1,ivert) = wedgever(1,ivert) + nv(1)*delta_offset*0.01
          wedgever(2,ivert) = wedgever(2,ivert) + nv(2)*delta_offset*0.01
          wedgever(3,ivert) = wedgever(3,ivert) + nv(3)*delta_offset*0.01
          enddo

	      do ivert = 13,15
          wedgever2(1,ivert) = wedgever2(1,ivert) + nv(1)*delta_offset*0.01
          wedgever2(2,ivert) = wedgever2(2,ivert) + nv(2)*delta_offset*0.01
          wedgever2(3,ivert) = wedgever2(3,ivert) + nv(3)*delta_offset*0.01
          enddo

	      endif

          do ivert = 1,6
          ivert_index_exo = vert_index_exo + ivert  
          x_exo(connect(ivert_index_exo)) = wedgever(1,ivert)
          y_exo(connect(ivert_index_exo)) = wedgever(2,ivert) 
          z_exo(connect(ivert_index_exo)) = wedgever(3,ivert)
          enddo
		  
          do ivert = 7,15
          ivert_index_exo = vert_index_exo + ivert  
          x_exo(connect(ivert_index_exo)) = wedgever2(1,ivert) 
          y_exo(connect(ivert_index_exo)) = wedgever2(2,ivert)
          z_exo(connect(ivert_index_exo)) = wedgever2(3,ivert) 
          enddo
		  

       else 
        exit    ! break if no non-right-hand elements
       endif

       enddo

       vert_index_exo = vert_index_exo + 15		  

		  
          else if  (quadratic_option.eq.2) then
		 
!  given wedge15 vertices, and return you 3 hex coords.
!       call wedgetohex(hexver,wedgever,hexss,wedgess)
!       call wedgetohex_quadratic(hexver,wedgever,hexss,wedgess)
    
! loop until no non-right-hand elements are seem
 
       icorrection = 0

       do while (.True.)
	   
       hasnrh = .false.
	   
       do ivert = 1,15
       ivert_index_exo = vert_index_exo + ivert  
       wedgever(1,ivert) = x_exo(connect(ivert_index_exo))
       wedgever(2,ivert) = y_exo(connect(ivert_index_exo))
       wedgever(3,ivert) = z_exo(connect(ivert_index_exo))
       enddo	

       !call wedgetohex_quadratic(hexver,wedgever,hexss,wedgess,hasnrh)
       call wedgetohex_quadratic2(hexver,wedgever,hexss,wedgess,hasnrh)
!       do i = 1,3
!       iel_nek = iel_nek + 1
!       call if_in_nj_list(iel_nek,ifnj)
!       hasnrh = (hasnrh.or.ifnj) ! 
!       enddo

       if (hasnrh) then
	   
          icorrection = icorrection + 1
          if (icorrection.gt.10) then
           write(6,*) 'stuck in non-right-hand correction of wedge15'
           write(6,*) 'at ', wedgever(1,1),wedgever(2,1),wedgever(3,1)
           exit
          endif
	   
          ! adjust mid-edge node position if non-right-hand elements detected
          call average2vec(vert1(1),wedgever(1,1),wedgever(1,2))
          call average2vec(vert2(1),vert1(1),wedgever(1,7))
          call assignvec(wedgever(1,7),vert2(1))
		  
          call average2vec(vert1(1),wedgever(1,2),wedgever(1,3))
          call average2vec(vert2(1),vert1(1),wedgever(1,8))
          call assignvec(wedgever(1,8),vert2(1))
		  
          call average2vec(vert1(1),wedgever(1,1),wedgever(1,3))
          call average2vec(vert2(1),vert1(1),wedgever(1,9))
          call assignvec(wedgever(1,9),vert2(1))
		  
          call average2vec(vert1(1),wedgever(1,1),wedgever(1,4))
          call average2vec(vert2(1),vert1(1),wedgever(1,10))
          call assignvec(wedgever(1,10),vert2(1))
		  
          call average2vec(vert1(1),wedgever(1,2),wedgever(1,5))
          call average2vec(vert2(1),vert1(1),wedgever(1,11))
          call assignvec(wedgever(1,11),vert2(1))
		  
		  
          call average2vec(vert1(1),wedgever(1,3),wedgever(1,6))
          call average2vec(vert2(1),vert1(1),wedgever(1,12))
          call assignvec(wedgever(1,12),vert2(1))
		  
		  
          call average2vec(vert1(1),wedgever(1,4),wedgever(1,5))
          call average2vec(vert2(1),vert1(1),wedgever(1,13))
          call assignvec(wedgever(1,13),vert2(1))

          call average2vec(vert1(1),wedgever(1,5),wedgever(1,6))
          call average2vec(vert2(1),vert1(1),wedgever(1,14))
          call assignvec(wedgever(1,14),vert2(1))
		  
		  
          call average2vec(vert1(1),wedgever(1,4),wedgever(1,6))
          call average2vec(vert2(1),vert1(1),wedgever(1,15))
          call assignvec(wedgever(1,15),vert2(1))

          do ivert = 1, 15
          ivert_index_exo = vert_index_exo + ivert  
          x_exo(connect(ivert_index_exo)) = wedgever(1,ivert)
          y_exo(connect(ivert_index_exo)) = wedgever(2,ivert) 
          z_exo(connect(ivert_index_exo)) = wedgever(3,ivert)
          enddo

       else 
        exit    ! break if no non-right-hand elements
       endif

      enddo
	   
      vert_index_exo = vert_index_exo + 15

          else if (quadratic_option.eq.3) then

!  read wedge 15
       do ivert = 1,15
       vert_index_exo = vert_index_exo + 1  
       wedgever(1,ivert) = x_exo(connect(vert_index_exo))
       wedgever(2,ivert) = y_exo(connect(vert_index_exo))
       wedgever(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

!  assign sideset to wedge elements
       n = 5
       call rzero_int2(wedgess,n)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,5
          wedgess(ifc_exo) = exoss(ifc_exo,iel_exo_g)
        enddo
       endif

!  given wedge15 vertices, and return you 3 hex coords.
!       call wedgetohex(hexver,wedgever,hexss,wedgess)
!       call wedgetohex_quadratic(hexver,wedgever,hexss,wedgess,hasnrh)
       call wedgetohex_quadratic2(hexver,wedgever,hexss,wedgess,hasnrh)

	   
       do ihex = 1,6
          iel_nek = iel_nek + 1
          eacc = eacc + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)&
       + shiftvector(1)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)& 
       + shiftvector(2)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)&
       + shiftvector(3)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              bc(4,ifc,iel_nek) = iexo
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

          endif

          enddo  !           do iel_exo = 1,num_elem_in_block(iblk)
	
        else if (nvert.eq.20) then ! splitting hex20
		
          do iel_exo = 1,num_elem_in_block(iblk)
           !write(6,*) 'iel_exo, ',iel_exo

         iel_exo_g = iel_exo ! iel_exo is element number in this block
		 !   iel_exo_g is global exo element number
         if (iblk.gt.1) then
            do iblk1 = 1,iblk-1
             iel_exo_g = iel_exo_g + num_elem_in_block(iblk1)
            enddo
         endif
		 
          if (quadratic_option.eq.1) then
! only read hex20, but doing nothing
       do ivert = 1,20
       vert_index_exo = vert_index_exo + 1  ! advance vertice index
       ehexver(1,ivert) = x_exo(connect(vert_index_exo))
       ehexver(2,ivert) = y_exo(connect(vert_index_exo))
       ehexver(3,ivert) = z_exo(connect(vert_index_exo))
       enddo 

          else if  (quadratic_option.eq.2) then

! linearize tet elements if non-right hand elements are created
       n = 6
       call rzero_int2(ehexss,n) 
       if (num_side_sets.ne.0) then
        do ifc_exo=1,6
         ehexss(ifc_exo) = 0
         ehexss(ifc_exo) = exoss(ifc_exo,iel_exo_g)
        enddo
       endif

       icorrection = 0

       do while (.True.)
       hasnrh = .false.
	   
       do ivert = 1, 20
       ivert_index_exo = vert_index_exo + ivert  
       ehexver(1,ivert) = x_exo(connect(ivert_index_exo))
       ehexver(2,ivert) = y_exo(connect(ivert_index_exo))
       ehexver(3,ivert) = z_exo(connect(ivert_index_exo))
       enddo

       call ehexto8hex_quadratic(hexver,ehexver,hexss,ehexss,hasnrh) 

       if (hasnrh) then

          icorrection = icorrection + 1
          if (icorrection.gt.10) then
           write(6,*) 'stuck in non-right-hand correction of hex20'
           write(6,*) 'at ', ehexver(1,1),ehexver(2,1),ehexver(3,1)
           exit
          endif
	   
          ! adjust mid-edge node position if non-right-hand elements detected
          call average2vec(vert1(1),ehexver(1,1),ehexver(1,2))
          call average2vec(vert2(1),vert1(1),ehexver(1,9))
          call assignvec(ehexver(1,9),vert2(1))

          call average2vec(vert1(1),ehexver(1,2),ehexver(1,3))
          call average2vec(vert2(1),vert1(1),ehexver(1,10))
          call assignvec(ehexver(1,10),vert2(1))

          call average2vec(vert1(1),ehexver(1,3),ehexver(1,4))
          call average2vec(vert2(1),vert1(1),ehexver(1,11))
          call assignvec(ehexver(1,11),vert2(1))

          call average2vec(vert1(1),ehexver(1,1),ehexver(1,4))
          call average2vec(vert2(1),vert1(1),ehexver(1,12))
          call assignvec(ehexver(1,12),vert2(1))

          call average2vec(vert1(1),ehexver(1,1),ehexver(1,5))
          call average2vec(vert2(1),vert1(1),ehexver(1,13))
          call assignvec(ehexver(1,13),vert2(1))

          call average2vec(vert1(1),ehexver(1,2),ehexver(1,6))
          call average2vec(vert2(1),vert1(1),ehexver(1,14))
          call assignvec(ehexver(1,14),vert2(1))

          call average2vec(vert1(1),ehexver(1,3),ehexver(1,7))
          call average2vec(vert2(1),vert1(1),ehexver(1,15))
          call assignvec(ehexver(1,15),vert2(1))

          call average2vec(vert1(1),ehexver(1,4),ehexver(1,8))
          call average2vec(vert2(1),vert1(1),ehexver(1,16))
          call assignvec(ehexver(1,16),vert2(1))

          call average2vec(vert1(1),ehexver(1,5),ehexver(1,6))
          call average2vec(vert2(1),vert1(1),ehexver(1,17))
          call assignvec(ehexver(1,17),vert2(1))

          call average2vec(vert1(1),ehexver(1,6),ehexver(1,7))
          call average2vec(vert2(1),vert1(1),ehexver(1,18))
          call assignvec(ehexver(1,18),vert2(1))

          call average2vec(vert1(1),ehexver(1,7),ehexver(1,8))
          call average2vec(vert2(1),vert1(1),ehexver(1,19))
          call assignvec(ehexver(1,19),vert2(1))

          call average2vec(vert1(1),ehexver(1,5),ehexver(1,8))
          call average2vec(vert2(1),vert1(1),ehexver(1,20))
          call assignvec(ehexver(1,20),vert2(1))

          do ivert = 1, 20
          ivert_index_exo = vert_index_exo + ivert  
          x_exo(connect(ivert_index_exo)) = ehexver(1,ivert) ! return adjusted tet
          y_exo(connect(ivert_index_exo)) = ehexver(2,ivert)
          z_exo(connect(ivert_index_exo)) = ehexver(3,ivert)
          enddo

       else
        exit    ! break if no non-right-hand elements
       endif

      enddo ! do while (.True.)
	  
       vert_index_exo = vert_index_exo + 20


      else if (quadratic_option.eq.3) then
! actual splitting.
       n = 6
       call rzero_int2(ehexss,6)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,6
         ehexss(ifc_exo) = 0
         ehexss(ifc_exo) = exoss(ifc_exo,iel_exo_g)
        enddo
       endif

       do ivert = 1, 20
       vert_index_exo = vert_index_exo + 1  
       ehexver(1,ivert) = x_exo(connect(vert_index_exo))
       ehexver(2,ivert) = y_exo(connect(vert_index_exo))
       ehexver(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

       call ehexto8hex_quadratic(hexver,ehexver,hexss,ehexss,hasnrh) 

       do ihex = 1,8
          iel_nek = iel_nek + 1
          eacc = eacc + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)&
       + shiftvector(1)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)& 
       + shiftvector(2)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)& 
       + shiftvector(3)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              bc(4,ifc,iel_nek) = iexo
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

        endif

        enddo       !     do iel_exo = 1,num_elem_in_block(iblk)


	    else 

        write(6,*) 'ERROR, invalid element type'
        endif

      enddo

      if (quadratic_option.eq.3) then
      write(6,*) 'Converted elements in nek:',iel_nek-eacc_old
      write(6,'(A)') 'Done :: Converting elements '

      deallocate(x_exo,y_exo,z_exo,connect)
      deallocate (elem_list,side_list)
      deallocate ( idblk )
      deallocate ( num_nodes_per_elem )
      deallocate ( num_attr           )
      deallocate ( num_elem_in_block  )
      deallocate ( num_sides_in_set   )
      deallocate ( idss               )
      endif

      return
      end
!-------------------------------------------------
!-----------------------------------------------------------------------
      subroutine transfinite_tri(tri,fc,mp)
! transfinite interpolation for tri6
! get fc and mp9
!
      real*8 tri(3,6),fc(3),mp(3,9)
      real*8 cr1(3,3),cr2(3,3),cr3(3,3)
      real*8 kp(3,6),beta,x(3)
	  
      integer option
 
! pack curve infos

      x(1) = 0.0
      x(2) = 0.5
      x(3) = 1.0
	  
      call assignvec(cr1(1,1),tri(1,1))
      call assignvec(cr1(1,2),tri(1,4))
      call assignvec(cr1(1,3),tri(1,2))
	  
      call assignvec(cr2(1,1),tri(1,2))
      call assignvec(cr2(1,2),tri(1,5))
      call assignvec(cr2(1,3),tri(1,3))
	  
      call assignvec(cr3(1,1),tri(1,3))
      call assignvec(cr3(1,2),tri(1,6))
      call assignvec(cr3(1,3),tri(1,1))

      beta = 2.0/3.0
      call curve(cr3,x,beta,kp(1,1))
      beta = 1.0/3.0
      call curve(cr3,x,beta,kp(1,2))

      beta = 2.0/3.0
      call curve(cr2,x,beta,kp(1,3))
      beta = 1.0/3.0
      call curve(cr2,x,beta,kp(1,4))
	  
      beta = 2.0/3.0
      call curve(cr1,x,beta,kp(1,5))
      beta = 1.0/3.0
      call curve(cr1,x,beta,kp(1,6))
	
!      if (option.eq.1) then
! this is from transfinite interpolation
      do i = 1,3
      fc(i) = (1.0/3.0)*(kp(i,6)+kp(i,1)-tri(i,1)) + (1.0/3.0)*(kp(i,4)+kp(i,5)-tri(i,2)) + (1.0/3.0)*(kp(i,2)+kp(i,3)-tri(i,3))
      enddo
!      else if (option.eq.2) then
!      ! avoid non-right-hand elements
!      do i = 1,3
!      fc(i) = (1.0/3.0)*(tri(i,4)+tri(i,5)+tri(i,6))
!      enddo
!      endif	  

      beta = 0.25
      call curve(cr1,x,beta,mp(1,1))
      beta = 0.75
      call curve(cr1,x,beta,mp(1,2))
      beta = 0.25
      call curve(cr2,x,beta,mp(1,3))
      beta = 0.75
      call curve(cr2,x,beta,mp(1,4))
      beta = 0.25
      call curve(cr3,x,beta,mp(1,5))
      beta = 0.75
      call curve(cr3,x,beta,mp(1,6))
	  
      
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0
	  
      call assignvec(cr1(1,1),tri(1,3))
      call assignvec(cr1(1,2),fc(1))
      call assignvec(cr1(1,3),tri(1,4))
      beta = 5.0/6.0
      call curve(cr1,x,beta,mp(1,7))
	  
      call assignvec(cr1(1,1),tri(1,1))
      call assignvec(cr1(1,2),fc(1))
      call assignvec(cr1(1,3),tri(1,5))
      beta = 5.0/6.0
      call curve(cr1,x,beta,mp(1,8))
	  
      call assignvec(cr1(1,1),tri(1,2))
      call assignvec(cr1(1,2),fc(1))
      call assignvec(cr1(1,3),tri(1,6))
      beta = 5.0/6.0
      call curve(cr1,x,beta,mp(1,9)) 
	 
      return
      end
!----------------------------------------------------------------------
      subroutine transfinite_quad(quad,fc,mp)
! transfinite interpolation for quad8
! get fc, mp12
      real*8 quad(3,8),beta,fc(3),mp(3,12)
      real*8 cr1(3,3),cr2(3,3),cr3(3,3),cr4(3,3)
      real*8 x(3)

      x(1) = 0.0
      x(2) = 0.5
      x(3) = 1.0
	  
! pack curve infos
      call assignvec(cr1(1,1),quad(1,1))
      call assignvec(cr1(1,2),quad(1,5))
      call assignvec(cr1(1,3),quad(1,2))
	  
      call assignvec(cr2(1,1),quad(1,2))
      call assignvec(cr2(1,2),quad(1,6))
      call assignvec(cr2(1,3),quad(1,3))
	  
      call assignvec(cr3(1,1),quad(1,4))
      call assignvec(cr3(1,2),quad(1,7))
      call assignvec(cr3(1,3),quad(1,3))
	  
      call assignvec(cr4(1,1),quad(1,1))
      call assignvec(cr4(1,2),quad(1,8))
      call assignvec(cr4(1,3),quad(1,4))

	  
!      call curve(cr1,x,alpha(1),cr1a)
!      call curve(cr3,x,alpha(1),cr3a)	  
!      call curve(cr2,x,alpha(2),cr2a)
!      call curve(cr4,x,alpha(2),cr4a)

!      do i = 1,3
!      res(i) = (1-alpha(2))*cr1a +  alpha(2)*cr3a 
!     & + alpha(1)*cr2a + (1-alpha(1))*cr4a	 
!	 & -( (1-alpha(1))*(1-alpha(2))*quad(i,1) 
!     & + (alpha(1))*(1-alpha(2))*quad(i,2) 
!     & + (alpha(1))*(alpha(2))*quad(i,3) 
!     & + (1-alpha(1))*(alpha(2))*quad(i,4))
!      enddo
	 
! get fc first	,
! this is from transfinite interpolation
      do i = 1,3
      fc(i) = 0.5*(quad(i,5)+quad(i,6)+quad(i,7)+quad(i,8)) -0.25*(quad(i,1)+quad(i,2)+quad(i,3)+quad(i,4))
      enddo
	  
!      do i = 1,3
!      fc(i) = 0.25*(quad(i,1)+quad(i,2)+quad(i,3)+quad(i,4))
!      enddo
	 
! get mp points	 
      beta = 0.25
      call curve(cr1,x,beta,mp(1,1))
      beta = 0.75
      call curve(cr1,x,beta,mp(1,2))
      beta = 0.25
      call curve(cr2,x,beta,mp(1,3))
      beta = 0.75
      call curve(cr2,x,beta,mp(1,4))
	  
      beta = 0.75
      call curve(cr3,x,beta,mp(1,5))
      beta = 0.25
      call curve(cr3,x,beta,mp(1,6))
      beta = 0.75
      call curve(cr4,x,beta,mp(1,7))
      beta = 0.25
      call curve(cr4,x,beta,mp(1,8))
	
! get mp9 and mp11	
      call assignvec(cr1(1,1),quad(1,5))
      call assignvec(cr1(1,2),fc)
      call assignvec(cr1(1,3),quad(1,7))
      beta = 0.25
      call curve(cr1,x,beta,mp(1,9))
      beta = 0.75
      call curve(cr1,x,beta,mp(1,11))
	 
! get mp10 and mp12	
      call assignvec(cr1(1,1),quad(1,8))
      call assignvec(cr1(1,2),fc)
      call assignvec(cr1(1,3),quad(1,6))
      beta = 0.75
      call curve(cr1,x,beta,mp(1,10))
      beta = 0.25
      call curve(cr1,x,beta,mp(1,12))

      return
      end
!----------------------------------------------------------------------
      subroutine transfinite_quad2(quad,fc)
! transfinite interpolation for quad8
! get fc
      real*8 quad(3,8),fc(3)
	  
      do i = 1,3
      fc(i) = 0.5*(quad(i,5)+quad(i,6)+quad(i,7)+quad(i,8)) -0.25*(quad(i,1)+quad(i,2)+quad(i,3)+quad(i,4))
      enddo

      return
      end
!----------------------------------------------------------------------
      subroutine curve(cr,x,alpha,res)
! interpolation from curve defined by 3 points
      real*8 cr(3,3),alpha,res(3)
      real*8 x(3),y(3),z(3),r
	  
      y(1) = cr(1,1)
      y(2) = cr(1,2)
      y(3) = cr(1,3)
	  
      call lagrange_interpolation(x,y,alpha,r)
      res(1) = r

      y(1) = cr(2,1)
      y(2) = cr(2,2)
      y(3) = cr(2,3)
	  
      call lagrange_interpolation(x,y,alpha,r)
      res(2) = r

      y(1) = cr(3,1)
      y(2) = cr(3,2)
      y(3) = cr(3,3)
	  
      call lagrange_interpolation(x,y,alpha,r)
      res(3) = r
	  
      return
      end
!----------------------------------------------------------------------
      subroutine lagrange_interpolation(x,y,alpha,res)
! use lagrange_interpolation for 
      real*8 x(3),y(3),alpha,res,coeff
      integer i,j

      res = 0.0
	  
      do i = 1,3
      coeff = 1.0
       do j = 1,3
         if(j.ne.i) coeff = coeff*(alpha-x(j))/(x(i)-x(j))
       enddo	  
      res = res + coeff*y(i)
      enddo
	  
!      res  = (1.0-alpha)*y(1) + alpha*y(3)

      return
      end
!----------------------------------------------------------------------
!--------------------------------------------------------------------
      subroutine tettohex_quadratic(hexver,tetver,hexss,tetss,hasnrh)
! quadratic version of tet-to-hex
      real*8 tetver(3,10) ! tet vertices
      real*8 tetface(3,4) ! tet face center
      real*8 tetcen(3,1)  ! tet vol center
      real*8 mfc(3,4)     ! mid point of face center and vol center
      real*8 hex8(3,8),edge12(3,12) ! edge12 here follows nek definition of edge 
      real*8 hexver(3,27,4) ! four hex vertices in nek format.
      real*8 tempvec(3,3) ! temperary vector variable

      integer tetss(4),hexss(6,4)

      integer option
      logical hasnrh,ifnrh,ifnegjac


      real*8 tri(3,6,4),mp(3,9,4),fc(3,4)
       
! patch 1st tri6 info
      call assignvec(tri(1,1,1),tetver(1,1))
      call assignvec(tri(1,2,1),tetver(1,2))
      call assignvec(tri(1,3,1),tetver(1,4))
      call assignvec(tri(1,4,1),tetver(1,5))
      call assignvec(tri(1,5,1),tetver(1,9))
      call assignvec(tri(1,6,1),tetver(1,8))
      
      call transfinite_tri(tri(1,1,1),fc(1,1),mp(1,1,1))

! patch 2st tri6 info
      call assignvec(tri(1,1,2),tetver(1,2))
      call assignvec(tri(1,2,2),tetver(1,3))
      call assignvec(tri(1,3,2),tetver(1,4))
      call assignvec(tri(1,4,2),tetver(1,6))
      call assignvec(tri(1,5,2),tetver(1,10))
      call assignvec(tri(1,6,2),tetver(1,9))
      
      call transfinite_tri(tri(1,1,2),fc(1,2),mp(1,1,2))

! patch 3rd tri6 info
      call assignvec(tri(1,1,3),tetver(1,1))
      call assignvec(tri(1,2,3),tetver(1,3))
      call assignvec(tri(1,3,3),tetver(1,4))
      call assignvec(tri(1,4,3),tetver(1,7))
      call assignvec(tri(1,5,3),tetver(1,10))
      call assignvec(tri(1,6,3),tetver(1,8))
      
      call transfinite_tri(tri(1,1,3),fc(1,3),mp(1,1,3))

! patch 4th tri6 info
      call assignvec(tri(1,1,4),tetver(1,1))
      call assignvec(tri(1,2,4),tetver(1,2))
      call assignvec(tri(1,3,4),tetver(1,3))
      call assignvec(tri(1,4,4),tetver(1,5))
      call assignvec(tri(1,5,4),tetver(1,6))
      call assignvec(tri(1,6,4),tetver(1,7))

      call transfinite_tri(tri(1,1,4),fc(1,4),mp(1,1,4))

      do i=1,24
      hexss(i,1)=0
      enddo

!  get tet vol center from avg of face center
      call average4vec(tetcen(1,1),fc(1,1),fc(1,2),fc(1,3),fc(1,4))

      call assignvec(tetface(1,1),fc(1,1))
      call assignvec(tetface(1,2),fc(1,2))
      call assignvec(tetface(1,3),fc(1,3))
      call assignvec(tetface(1,4),fc(1,4))

      call average2vec(mfc(1,1),fc(1,1),tetcen(1,1))
      call average2vec(mfc(1,2),fc(1,2),tetcen(1,1))
      call average2vec(mfc(1,3),fc(1,3),tetcen(1,1))
      call average2vec(mfc(1,4),fc(1,4),tetcen(1,1))


! the rest of it is very complicated now:
! need to map projected points to hex20 elements now.
! but first, map to hex8 and edge12, and then convert to hex27 in nek

      hasnrh = .false.
      ifnrh = .false.
      ifnegjac = .false.

!  assign coordinates to four hex.
!  hex 1
      hexss(1,1) = tetss(1)
      hexss(4,1) = tetss(3)
      hexss(5,1) = tetss(4)

      ! pack hex8
      call assignvec(hex8(1,1),tetver(1,1))
      call assignvec(hex8(1,2),tetver(1,5))
      call assignvec(hex8(1,3),tetface(1,4))
      call assignvec(hex8(1,4),tetver(1,7))
      call assignvec(hex8(1,5),tetver(1,8))
      call assignvec(hex8(1,6),tetface(1,1))
      call assignvec(hex8(1,7),tetcen(1,1))
      call assignvec(hex8(1,8),tetface(1,3))
      
      ! pack edge12
      call assignvec(edge12(1,1),mp(1,1,1))
      call assignvec(edge12(1,2),mp(1,7,4))
      call assignvec(edge12(1,3),mp(1,9,4))
      call assignvec(edge12(1,4),mp(1,6,4))
      call assignvec(edge12(1,5),mp(1,9,1))
      call assignvec(edge12(1,6),mfc(1,1))
      call assignvec(edge12(1,7),mfc(1,3))
      call assignvec(edge12(1,8),mp(1,9,3))
      call assignvec(edge12(1,9),mp(1,6,1))
      call assignvec(edge12(1,10),mp(1,7,1))
      call assignvec(edge12(1,11),mfc(1,4))
      call assignvec(edge12(1,12),mp(1,7,3))
  
      call hex20tohex27(hexver(1,1,1),hex8(1,1),edge12(1,1))
      !call hex8tohex27(hexver(1,1,1),hex8(1,1))
	  
      ifnrh = .false.
      ifnegjac = .false.
  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
!      call chk_jac(hex8,edge12,ifnegjac)
!      hasnrh = (hasnrh.or.ifnegjac) ! 
       
!  hex 2
      hexss(1,2) = tetss(1)
      hexss(2,2) = tetss(2)
      hexss(5,2) = tetss(4)

      call assignvec(hex8(1,1),tetver(1,5))
      call assignvec(hex8(1,2),tetver(1,2))
      call assignvec(hex8(1,3),tetver(1,6))
      call assignvec(hex8(1,4),tetface(1,4))
      call assignvec(hex8(1,5),tetface(1,1))
      call assignvec(hex8(1,6),tetver(1,9))
      call assignvec(hex8(1,7),tetface(1,2))
      call assignvec(hex8(1,8),tetcen(1,1))

      ! pack edge12
      call assignvec(edge12(1,1),mp(1,2,1))
      call assignvec(edge12(1,2),mp(1,1,2))
      call assignvec(edge12(1,3),mp(1,8,4))
      call assignvec(edge12(1,4),mp(1,7,4))
      call assignvec(edge12(1,5),mp(1,8,1))
      call assignvec(edge12(1,6),mp(1,9,2))
      call assignvec(edge12(1,7),mfc(1,2))
      call assignvec(edge12(1,8),mfc(1,1))
      call assignvec(edge12(1,9),mp(1,7,1))
      call assignvec(edge12(1,10),mp(1,3,1))
      call assignvec(edge12(1,11),mp(1,7,2))
      call assignvec(edge12(1,12),mfc(1,4))
  
      call hex20tohex27(hexver(1,1,2),hex8(1,1),edge12(1,1))
      !call hex8tohex27(hexver(1,1,2),hex8(1,1))

      ifnrh = .false.
      ifnegjac = .false.
  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
!      call chk_jac(hex8,edge12,ifnegjac)
!      hasnrh = (hasnrh.or.ifnegjac) ! 

!  hex 3
      hexss(2,3) = tetss(2)
      hexss(3,3) = tetss(3)
      hexss(5,3) = tetss(4)

      call assignvec(hex8(1,1),tetface(1,4))
      call assignvec(hex8(1,2),tetver(1,6))
      call assignvec(hex8(1,3),tetver(1,3))
      call assignvec(hex8(1,4),tetver(1,7))
      call assignvec(hex8(1,5),tetcen(1,1))
      call assignvec(hex8(1,6),tetface(1,2))
      call assignvec(hex8(1,7),tetver(1,10))
      call assignvec(hex8(1,8),tetface(1,3))

      ! pack edge12
      call assignvec(edge12(1,1),mp(1,8,4))
      call assignvec(edge12(1,2),mp(1,4,4))
      call assignvec(edge12(1,3),mp(1,5,4))
      call assignvec(edge12(1,4),mp(1,9,4))
      call assignvec(edge12(1,5),mfc(1,2))
      call assignvec(edge12(1,6),mp(1,8,2))
      call assignvec(edge12(1,7),mp(1,8,3))
      call assignvec(edge12(1,8),mfc(1,3))
      call assignvec(edge12(1,9),mfc(1,4))
      call assignvec(edge12(1,10),mp(1,7,2))
      call assignvec(edge12(1,11),mp(1,3,2))
      call assignvec(edge12(1,12),mp(1,7,3))
  
      call hex20tohex27(hexver(1,1,3),hex8(1,1),edge12(1,1))
      !call hex8tohex27(hexver(1,1,3),hex8(1,1))

      ifnrh = .false.
      ifnegjac = .false.
  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
!      call chk_jac(hex8,edge12,ifnegjac)
!      hasnrh = (hasnrh.or.ifnegjac) ! 

!  hex 4
      hexss(2,4) = tetss(2)
      hexss(3,4) = tetss(3)
      hexss(6,4) = tetss(1)

      call assignvec(hex8(1,1),tetcen(1,1))
      call assignvec(hex8(1,2),tetface(1,2))
      call assignvec(hex8(1,3),tetver(1,10))
      call assignvec(hex8(1,4),tetface(1,3))
      call assignvec(hex8(1,5),tetface(1,1))
      call assignvec(hex8(1,6),tetver(1,9))
      call assignvec(hex8(1,7),tetver(1,4))
      call assignvec(hex8(1,8),tetver(1,8))

      ! pack edge12
      call assignvec(edge12(1,1),mfc(1,2))
      call assignvec(edge12(1,2),mp(1,8,2))
      call assignvec(edge12(1,3),mp(1,8,3))
      call assignvec(edge12(1,4),mfc(1,3))
      call assignvec(edge12(1,5),mp(1,8,1))
      call assignvec(edge12(1,6),mp(1,5,2))
      call assignvec(edge12(1,7),mp(1,5,3))
      call assignvec(edge12(1,8),mp(1,9,1))
      call assignvec(edge12(1,9),mfc(1,1))
      call assignvec(edge12(1,10),mp(1,9,2))
      call assignvec(edge12(1,11),mp(1,4,2))
      call assignvec(edge12(1,12),mp(1,9,3))
  
      call hex20tohex27(hexver(1,1,4),hex8(1,1),edge12(1,1))
      !call hex8tohex27(hexver(1,1,4),hex8(1,1))

      ifnrh = .false.
      ifnegjac = .false.
  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
!      call chk_jac(hex8,edge12,ifnegjac)
!      hasnrh = (hasnrh.or.ifnegjac) ! 
	  
      return
      end	  
!---------------------------------------------------------------------------------
!--------------------------------------------------------------------
      subroutine wedgetohex_quadratic(hexver,wedgever,hexss,wedgess,hasnrh)
! quadratic version of wedge-to-hex
      real*8 wedgever(3,15) ! tet vertices, actually only wedge15 is needed, but ICEM cannot dump wedge15...
      real*8 wedgeface(3,5) ! tet face center
      real*8 wedgecen(3,1) ! tet vol center
      real*8 hex8(3,8),edge12(3,12) ! edge12 here follows nek definition of edge 
      real*8 hexver(3,27,4) ! four hex vertices in nek format.
      real*8 tempvec(3,3) ! temperary vector variable

      integer wedgess(5),hexss(6,4)

      real*8 tri(3,6,2),tmp(3,9,2),tfc(3,2)
      real*8 quad(3,8,3),qmp(3,12,3),qfc(3,3)

      logical hasnrh,ifnrh,ifnegjac

! patch 1st tri6 info
      call assignvec(tri(1,1,1),wedgever(1,1))
      call assignvec(tri(1,2,1),wedgever(1,2))
      call assignvec(tri(1,3,1),wedgever(1,3))
      call assignvec(tri(1,4,1),wedgever(1,7))
      call assignvec(tri(1,5,1),wedgever(1,8))
      call assignvec(tri(1,6,1),wedgever(1,9))
      
      call transfinite_tri(tri(1,1,1),tfc(1,1),tmp(1,1,1))

! patch 2nd tri6 info
      call assignvec(tri(1,1,2),wedgever(1,4))
      call assignvec(tri(1,2,2),wedgever(1,5))
      call assignvec(tri(1,3,2),wedgever(1,6))
      call assignvec(tri(1,4,2),wedgever(1,13))
      call assignvec(tri(1,5,2),wedgever(1,14))
      call assignvec(tri(1,6,2),wedgever(1,15))
      
      call transfinite_tri(tri(1,1,2),tfc(1,2),tmp(1,1,2))

      call average2vec(wedgecen(1,1),tfc(1,1),tfc(1,2))

! patch 1st quad8 info
      call assignvec(quad(1,1,1),wedgever(1,1))
      call assignvec(quad(1,2,1),wedgever(1,2))
      call assignvec(quad(1,3,1),wedgever(1,5))
      call assignvec(quad(1,4,1),wedgever(1,4))
      call assignvec(quad(1,5,1),wedgever(1,7))
      call assignvec(quad(1,6,1),wedgever(1,11))
      call assignvec(quad(1,7,1),wedgever(1,13))
      call assignvec(quad(1,8,1),wedgever(1,10))
      
      call transfinite_quad(quad(1,1,1),qfc(1,1),qmp(1,1,1)) 

! patch 2nd quad8 info
      call assignvec(quad(1,1,2),wedgever(1,3))
      call assignvec(quad(1,2,2),wedgever(1,2))
      call assignvec(quad(1,3,2),wedgever(1,5))
      call assignvec(quad(1,4,2),wedgever(1,6))
      call assignvec(quad(1,5,2),wedgever(1,8))
      call assignvec(quad(1,6,2),wedgever(1,11))
      call assignvec(quad(1,7,2),wedgever(1,14))
      call assignvec(quad(1,8,2),wedgever(1,12))
      
      call transfinite_quad(quad(1,1,2),qfc(1,2),qmp(1,1,2))      
      
! patch 3nd quad8 info
      call assignvec(quad(1,1,3),wedgever(1,1))
      call assignvec(quad(1,2,3),wedgever(1,3))
      call assignvec(quad(1,3,3),wedgever(1,6))
      call assignvec(quad(1,4,3),wedgever(1,4))
      call assignvec(quad(1,5,3),wedgever(1,9))
      call assignvec(quad(1,6,3),wedgever(1,12))
      call assignvec(quad(1,7,3),wedgever(1,15))
      call assignvec(quad(1,8,3),wedgever(1,10))
      
      call transfinite_quad(quad(1,1,3),qfc(1,3),qmp(1,1,3))    

      do i=1,6*4
      hexss(i,1)=0
      enddo

!  assign coordinates to 3 hex.
!  hex 1
      hexss(1,1) = wedgess(1)
      hexss(4,1) = wedgess(3)
      hexss(5,1) = wedgess(4)
      hexss(6,1) = wedgess(5)

      call assignvec(hex8(1,1),wedgever(1,1))
      call assignvec(hex8(1,2),wedgever(1,7))
      call assignvec(hex8(1,3),tfc(1,1))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgever(1,4))
      call assignvec(hex8(1,6),wedgever(1,13))
      call assignvec(hex8(1,7),tfc(1,2))
      call assignvec(hex8(1,8),wedgever(1,15))
      
      ! pack edge12
      call assignvec(edge12(1,1),tmp(1,1,1))
      call assignvec(edge12(1,2),tmp(1,7,1))
      call assignvec(edge12(1,3),tmp(1,9,1))
      call assignvec(edge12(1,4),tmp(1,6,1))
      call assignvec(edge12(1,5),tmp(1,1,2))
      call assignvec(edge12(1,6),tmp(1,7,2))
      call assignvec(edge12(1,7),tmp(1,9,2))
      call assignvec(edge12(1,8),tmp(1,6,2))
      call assignvec(edge12(1,9),wedgever(1,10))
      call assignvec(edge12(1,10),qfc(1,1))
      call assignvec(edge12(1,11),wedgecen(1,1))
      call assignvec(edge12(1,12),qfc(1,3))
  
      call hex20tohex27(hexver(1,1,1),hex8(1,1),edge12(1,1))
      !call hex8tohex27(hexver(1,1,1),hex8(1,1))

      ifnrh = .false.
      ifnegjac = .false.
  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
!      call chk_jac(hex8,edge12,ifnegjac)
!      hasnrh = (hasnrh.or.ifnegjac) ! 
      
!  hex 2
      hexss(1,2) = wedgess(1)
      hexss(2,2) = wedgess(2)
      hexss(5,2) = wedgess(4)
      hexss(6,2) = wedgess(5)

      call assignvec(hex8(1,1),wedgever(1,7))
      call assignvec(hex8(1,2),wedgever(1,2))
      call assignvec(hex8(1,3),wedgever(1,8))
      call assignvec(hex8(1,4),tfc(1,1))
      call assignvec(hex8(1,5),wedgever(1,13))
      call assignvec(hex8(1,6),wedgever(1,5))
      call assignvec(hex8(1,7),wedgever(1,14))
      call assignvec(hex8(1,8),tfc(1,2))
      
      ! pack edge12
      call assignvec(edge12(1,1),tmp(1,2,1))
      call assignvec(edge12(1,2),tmp(1,3,1))
      call assignvec(edge12(1,3),tmp(1,8,1))
      call assignvec(edge12(1,4),tmp(1,7,1))
      call assignvec(edge12(1,5),tmp(1,2,2))
      call assignvec(edge12(1,6),tmp(1,3,2))
      call assignvec(edge12(1,7),tmp(1,8,2))
      call assignvec(edge12(1,8),tmp(1,7,2))
      call assignvec(edge12(1,9),qfc(1,1))
      call assignvec(edge12(1,10),wedgever(1,11))
      call assignvec(edge12(1,11),qfc(1,2))
      call assignvec(edge12(1,12),wedgecen(1,1))
  
      call hex20tohex27(hexver(1,1,2),hex8(1,1),edge12(1,1))
      !call hex8tohex27(hexver(1,1,2),hex8(1,1))

      ifnrh = .false.
      ifnegjac = .false.
  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
!      call chk_jac(hex8,edge12,ifnegjac)
!      hasnrh = (hasnrh.or.ifnegjac) ! 
      
      
!  hex 3
      hexss(2,3) = wedgess(2)
      hexss(3,3) = wedgess(3)
      hexss(5,3) = wedgess(4)
      hexss(6,3) = wedgess(5)

      call assignvec(hex8(1,1),tfc(1,1))
      call assignvec(hex8(1,2),wedgever(1,8))
      call assignvec(hex8(1,3),wedgever(1,3))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),tfc(1,2))
      call assignvec(hex8(1,6),wedgever(1,14))
      call assignvec(hex8(1,7),wedgever(1,6))
      call assignvec(hex8(1,8),wedgever(1,15))
      
      ! pack edge12
      call assignvec(edge12(1,1),tmp(1,8,1))
      call assignvec(edge12(1,2),tmp(1,4,1))
      call assignvec(edge12(1,3),tmp(1,5,1))
      call assignvec(edge12(1,4),tmp(1,9,1))
      call assignvec(edge12(1,5),tmp(1,8,2))
      call assignvec(edge12(1,6),tmp(1,4,2))
      call assignvec(edge12(1,7),tmp(1,5,2))
      call assignvec(edge12(1,8),tmp(1,9,2))
      call assignvec(edge12(1,9),wedgecen(1,1))
      call assignvec(edge12(1,10),qfc(1,2))
      call assignvec(edge12(1,11),wedgever(1,12))
      call assignvec(edge12(1,12),qfc(1,3))
  
      call hex20tohex27(hexver(1,1,3),hex8(1,1),edge12(1,1))
      !call hex8tohex27(hexver(1,1,3),hex8(1,1))

      ifnrh = .false.
      ifnegjac = .false.
  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
!      call chk_jac(hex8,edge12,ifnegjac)
!      hasnrh = (hasnrh.or.ifnegjac) ! 
	  
	  
      return
      end	 
!--------------------------------------------------------------------------s
!--------------------------------------------------------------------
      subroutine wedgetohex_quadratic2(hexver,wedgever,hexss,wedgess,hasnrh)
! quadratic version of wedge-to-hex
      real*8 wedgever(3,15) ! tet vertices, actually only wedge15 is needed, but ICEM cannot dump wedge15...
      real*8 wedgeface(3,5) ! tet face center
      real*8 wedgecen(3,1) ! tet vol center
      real*8 hex8(3,8),edge12(3,12) ! edge12 here follows nek definition of edge 
      real*8 hexver(3,27,6) ! four hex vertices in nek format.
      real*8 tempvec(3,3) ! temperary vector variable

      integer wedgess(5),hexss(6,6)

      real*8 tri(3,6,2),tmp(3,9,2),tfc(3,2)
      real*8 quad(3,8,3),qmp(3,12,3),qfc(3,3)
    
	  real*8 cr1(3,3),x(3),beta,ver(3)

      logical hasnrh,ifnrh,ifnegjac

! patch 1st tri6 info
      call assignvec(tri(1,1,1),wedgever(1,1))
      call assignvec(tri(1,2,1),wedgever(1,2))
      call assignvec(tri(1,3,1),wedgever(1,3))
      call assignvec(tri(1,4,1),wedgever(1,7))
      call assignvec(tri(1,5,1),wedgever(1,8))
      call assignvec(tri(1,6,1),wedgever(1,9))
      
      call transfinite_tri(tri(1,1,1),tfc(1,1),tmp(1,1,1))
      call assignvec(wedgeface(1,4),tfc(1,1))

! patch 2nd tri6 info
      call assignvec(tri(1,1,2),wedgever(1,4))
      call assignvec(tri(1,2,2),wedgever(1,5))
      call assignvec(tri(1,3,2),wedgever(1,6))
      call assignvec(tri(1,4,2),wedgever(1,13))
      call assignvec(tri(1,5,2),wedgever(1,14))
      call assignvec(tri(1,6,2),wedgever(1,15))
      
      call transfinite_tri(tri(1,1,2),tfc(1,2),tmp(1,1,2))
      call assignvec(wedgeface(1,5),tfc(1,2))

      call average2vec(wedgecen(1,1),tfc(1,1),tfc(1,2))

! patch 1st quad8 info
      call assignvec(quad(1,1,1),wedgever(1,1))
      call assignvec(quad(1,2,1),wedgever(1,2))
      call assignvec(quad(1,3,1),wedgever(1,5))
      call assignvec(quad(1,4,1),wedgever(1,4))
      call assignvec(quad(1,5,1),wedgever(1,7))
      call assignvec(quad(1,6,1),wedgever(1,11))
      call assignvec(quad(1,7,1),wedgever(1,13))
      call assignvec(quad(1,8,1),wedgever(1,10))
      
      call transfinite_quad(quad(1,1,1),qfc(1,1),qmp(1,1,1)) 
      call assignvec(wedgeface(1,1),qfc(1,1))

! patch 2nd quad8 info
      call assignvec(quad(1,1,2),wedgever(1,3))
      call assignvec(quad(1,2,2),wedgever(1,2))
      call assignvec(quad(1,3,2),wedgever(1,5))
      call assignvec(quad(1,4,2),wedgever(1,6))
      call assignvec(quad(1,5,2),wedgever(1,8))
      call assignvec(quad(1,6,2),wedgever(1,11))
      call assignvec(quad(1,7,2),wedgever(1,14))
      call assignvec(quad(1,8,2),wedgever(1,12))
      
      call transfinite_quad(quad(1,1,2),qfc(1,2),qmp(1,1,2))      
      call assignvec(wedgeface(1,2),qfc(1,2))
	
! patch 3nd quad8 info
      call assignvec(quad(1,1,3),wedgever(1,1))
      call assignvec(quad(1,2,3),wedgever(1,3))
      call assignvec(quad(1,3,3),wedgever(1,6))
      call assignvec(quad(1,4,3),wedgever(1,4))
      call assignvec(quad(1,5,3),wedgever(1,9))
      call assignvec(quad(1,6,3),wedgever(1,12))
      call assignvec(quad(1,7,3),wedgever(1,15))
      call assignvec(quad(1,8,3),wedgever(1,10))
      
      call transfinite_quad(quad(1,1,3),qfc(1,3),qmp(1,1,3))    
      call assignvec(wedgeface(1,3),qfc(1,3))

      do i=1,6*6
      hexss(i,1)=0
      enddo

!  assign coordinates to 6 hex.
!  hex 1
      hexss(1,1) = wedgess(1)
      hexss(4,1) = wedgess(3)
      hexss(5,1) = wedgess(4)

      call assignvec(hex8(1,1),wedgever(1,1))
      call assignvec(hex8(1,2),wedgever(1,7))
      call assignvec(hex8(1,3),wedgeface(1,4))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgever(1,10))
      call assignvec(hex8(1,6),wedgeface(1,1))
      call assignvec(hex8(1,7),wedgecen(1,1))
      call assignvec(hex8(1,8),wedgeface(1,3))


      call quick_quater_point_from_curve(ver(1),wedgever(1,1),wedgever(1,7),wedgever(1,2))
      call assignvec(edge12(1,1),ver(1)) 

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,3))
      call assignvec(cr1(1,2),wedgeface(1,4))
      call assignvec(cr1(1,3),wedgever(1,7)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,2))	  

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,2))
      call assignvec(cr1(1,2),wedgeface(1,4))
      call assignvec(cr1(1,3),wedgever(1,9)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,3))	  

      call quick_quater_point_from_curve(ver(1),wedgever(1,1),wedgever(1,9),wedgever(1,3))
      call assignvec(edge12(1,4),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgever(1,10),wedgeface(1,1),wedgever(1,11))
      call assignvec(edge12(1,5),ver(1)) 
 
 
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,12))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,1)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,6))	  
 
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,11))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,3)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,7))	  
  
      call quick_quater_point_from_curve(ver(1),wedgever(1,10),wedgeface(1,3),wedgever(1,12))
      call assignvec(edge12(1,8),ver(1)) 
   
   
      call quick_quater_point_from_curve(ver(1),wedgever(1,1),wedgever(1,10),wedgever(1,4))
      call assignvec(edge12(1,9),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgever(1,7),wedgeface(1,1),wedgever(1,13))
      call assignvec(edge12(1,10),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgeface(1,4),wedgecen(1,1),wedgeface(1,5))
      call assignvec(edge12(1,11),ver(1)) 
 
      call quick_quater_point_from_curve(ver(1),wedgever(1,9),wedgeface(1,3),wedgever(1,15))
      call assignvec(edge12(1,12),ver(1)) 
   
      call hex20tohex27(hexver(1,1,1),hex8(1,1),edge12(1,1))

      ifnrh = .false.
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 

!  hex 2
      hexss(1,2) = wedgess(1)
      hexss(2,2) = wedgess(2)
      hexss(5,2) = wedgess(4)

      call assignvec(hex8(1,1),wedgever(1,7))
      call assignvec(hex8(1,2),wedgever(1,2))
      call assignvec(hex8(1,3),wedgever(1,8))
      call assignvec(hex8(1,4),wedgeface(1,4))
      call assignvec(hex8(1,5),wedgeface(1,1))
      call assignvec(hex8(1,6),wedgever(1,11))
      call assignvec(hex8(1,7),wedgeface(1,2))
      call assignvec(hex8(1,8),wedgecen(1,1))
      
  
      call quick_quater_point_from_curve(ver(1),wedgever(1,2),wedgever(1,7),wedgever(1,1))
      call assignvec(edge12(1,1),ver(1)) 
      call quick_quater_point_from_curve(ver(1),wedgever(1,2),wedgever(1,8),wedgever(1,3))
      call assignvec(edge12(1,2),ver(1)) 
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,1))
      call assignvec(cr1(1,2),wedgeface(1,4))
      call assignvec(cr1(1,3),wedgever(1,8)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,3))	  
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,3))
      call assignvec(cr1(1,2),wedgeface(1,4))
      call assignvec(cr1(1,3),wedgever(1,7)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,4))	  

      call quick_quater_point_from_curve(ver(1),wedgever(1,11),wedgeface(1,1),wedgever(1,10))
      call assignvec(edge12(1,5),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgever(1,11),wedgeface(1,2),wedgever(1,12))
      call assignvec(edge12(1,6),ver(1)) 

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,10))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,2)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,7))	  

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,12))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,1)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,8))	  

      call quick_quater_point_from_curve(ver(1),wedgever(1,7),wedgeface(1,1),wedgever(1,13))
      call assignvec(edge12(1,9),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgever(1,2),wedgever(1,11),wedgever(1,5))
      call assignvec(edge12(1,10),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgever(1,8),wedgeface(1,2),wedgever(1,14))
      call assignvec(edge12(1,11),ver(1)) 
	  
      call quick_quater_point_from_curve(ver(1),wedgeface(1,4),wedgecen(1,1),wedgeface(1,5))
      call assignvec(edge12(1,12),ver(1)) 


      call hex20tohex27(hexver(1,1,2),hex8(1,1),edge12(1,1))
      !call hex8tohex27(hexver(1,1,1),hex8(1,1))

      ifnrh = .false.
  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 

!  hex 3
      hexss(2,3) = wedgess(2)
      hexss(3,3) = wedgess(3)
      hexss(5,3) = wedgess(4)

      call assignvec(hex8(1,1),wedgeface(1,4))
      call assignvec(hex8(1,2),wedgever(1,8))
      call assignvec(hex8(1,3),wedgever(1,3))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgecen(1,1))
      call assignvec(hex8(1,6),wedgeface(1,2))
      call assignvec(hex8(1,7),wedgever(1,12))
      call assignvec(hex8(1,8),wedgeface(1,3))

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,1))
      call assignvec(cr1(1,2),wedgeface(1,4))
      call assignvec(cr1(1,3),wedgever(1,8)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,1))	  

      call quick_quater_point_from_curve(ver(1),wedgever(1,3),wedgever(1,8),wedgever(1,2))
      call assignvec(edge12(1,2),ver(1)) 
      call quick_quater_point_from_curve(ver(1),wedgever(1,3),wedgever(1,9),wedgever(1,1))
      call assignvec(edge12(1,3),ver(1)) 

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,2))
      call assignvec(cr1(1,2),wedgeface(1,4))
      call assignvec(cr1(1,3),wedgever(1,9)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,4))	  

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,10))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,2)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,5))	  

      call quick_quater_point_from_curve(ver(1),wedgever(1,12),wedgeface(1,2),wedgever(1,11))
      call assignvec(edge12(1,6),ver(1)) 
      call quick_quater_point_from_curve(ver(1),wedgever(1,12),wedgeface(1,3),wedgever(1,10))
      call assignvec(edge12(1,7),ver(1)) 

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,11))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,3)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,8))	  


      call quick_quater_point_from_curve(ver(1),wedgeface(1,4),wedgecen(1,1),wedgeface(1,5))
      call assignvec(edge12(1,9),ver(1))
      call quick_quater_point_from_curve(ver(1),wedgever(1,8),wedgeface(1,2),wedgever(1,14))
      call assignvec(edge12(1,10),ver(1)) 
      call quick_quater_point_from_curve(ver(1),wedgever(1,3),wedgever(1,12),wedgever(1,6))
      call assignvec(edge12(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),wedgever(1,9),wedgeface(1,3),wedgever(1,15))
      call assignvec(edge12(1,12),ver(1)) 


      call hex20tohex27(hexver(1,1,3),hex8(1,1),edge12(1,1))
      !call hex8tohex27(hexver(1,1,2),hex8(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 

      
!  hex 4
      hexss(1,4) = wedgess(1)
      hexss(4,4) = wedgess(3)
      hexss(6,4) = wedgess(5)

      call assignvec(hex8(1,1),wedgever(1,10))
      call assignvec(hex8(1,2),wedgeface(1,1))
      call assignvec(hex8(1,3),wedgecen(1,1))
      call assignvec(hex8(1,4),wedgeface(1,3))
      call assignvec(hex8(1,5),wedgever(1,4))
      call assignvec(hex8(1,6),wedgever(1,13))
      call assignvec(hex8(1,7),wedgeface(1,5))
      call assignvec(hex8(1,8),wedgever(1,15))
	  
      call quick_quater_point_from_curve(ver(1),wedgever(1,10),wedgeface(1,1),wedgever(1,11))
      call assignvec(edge12(1,1),ver(1)) 
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,12))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,1)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,2))	  
 
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,11))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,3)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,3))	  
  
      call quick_quater_point_from_curve(ver(1),wedgever(1,10),wedgeface(1,3),wedgever(1,12))
      call assignvec(edge12(1,4),ver(1)) 

    
 
      call quick_quater_point_from_curve(ver(1),wedgever(1,4),wedgever(1,13),wedgever(1,5))
      call assignvec(edge12(1,5),ver(1)) 
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,6))
      call assignvec(cr1(1,2),wedgeface(1,5))
      call assignvec(cr1(1,3),wedgever(1,13)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,6))	  
 
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,5))
      call assignvec(cr1(1,2),wedgeface(1,5))
      call assignvec(cr1(1,3),wedgever(1,15)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,7))	  
  
      call quick_quater_point_from_curve(ver(1),wedgever(1,4),wedgever(1,15),wedgever(1,6))
      call assignvec(edge12(1,8),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgever(1,4),wedgever(1,10),wedgever(1,1))
      call assignvec(edge12(1,9),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgever(1,13),wedgeface(1,1),wedgever(1,7))
      call assignvec(edge12(1,10),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgeface(1,5),wedgecen(1,1),wedgeface(1,4))
      call assignvec(edge12(1,11),ver(1)) 
 
      call quick_quater_point_from_curve(ver(1),wedgever(1,15),wedgeface(1,3),wedgever(1,9))
      call assignvec(edge12(1,12),ver(1)) 


      call hex20tohex27(hexver(1,1,4),hex8(1,1),edge12(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
	  
!  hex 5
      hexss(1,5) = wedgess(1)
      hexss(2,5) = wedgess(2)
      hexss(6,5) = wedgess(5)

      call assignvec(hex8(1,1),wedgeface(1,1))
      call assignvec(hex8(1,2),wedgever(1,11))
      call assignvec(hex8(1,3),wedgeface(1,2))
      call assignvec(hex8(1,4),wedgecen(1,1))
      call assignvec(hex8(1,5),wedgever(1,13))
      call assignvec(hex8(1,6),wedgever(1,5))
      call assignvec(hex8(1,7),wedgever(1,14))
      call assignvec(hex8(1,8),wedgeface(1,5))
   
      call quick_quater_point_from_curve(ver(1),wedgever(1,11),wedgeface(1,1),wedgever(1,10))
      call assignvec(edge12(1,1),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgever(1,11),wedgeface(1,2),wedgever(1,12))
      call assignvec(edge12(1,2),ver(1)) 

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,10))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,2)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,3))	  

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,12))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,1)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,4))	

   
      call quick_quater_point_from_curve(ver(1),wedgever(1,5),wedgever(1,13),wedgever(1,4))
      call assignvec(edge12(1,5),ver(1)) 
      call quick_quater_point_from_curve(ver(1),wedgever(1,5),wedgever(1,14),wedgever(1,6))
      call assignvec(edge12(1,6),ver(1)) 
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,4))
      call assignvec(cr1(1,2),wedgeface(1,5))
      call assignvec(cr1(1,3),wedgever(1,14)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,7))	  
      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,6))
      call assignvec(cr1(1,2),wedgeface(1,5))
      call assignvec(cr1(1,3),wedgever(1,13)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,8))	  

      call quick_quater_point_from_curve(ver(1),wedgever(1,13),wedgeface(1,1),wedgever(1,7))
      call assignvec(edge12(1,9),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgever(1,5),wedgever(1,11),wedgever(1,2))
      call assignvec(edge12(1,10),ver(1)) 

      call quick_quater_point_from_curve(ver(1),wedgever(1,14),wedgeface(1,2),wedgever(1,8))
      call assignvec(edge12(1,11),ver(1)) 
	  
      call quick_quater_point_from_curve(ver(1),wedgeface(1,5),wedgecen(1,1),wedgeface(1,4))
      call assignvec(edge12(1,12),ver(1)) 

      call hex20tohex27(hexver(1,1,5),hex8(1,1),edge12(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
   
!  hex 6
      hexss(2,6) = wedgess(2)
      hexss(3,6) = wedgess(3)
      hexss(6,6) = wedgess(5)

      call assignvec(hex8(1,1),wedgecen(1,1))
      call assignvec(hex8(1,2),wedgeface(1,2))
      call assignvec(hex8(1,3),wedgever(1,12))
      call assignvec(hex8(1,4),wedgeface(1,3))
      call assignvec(hex8(1,5),wedgeface(1,5))
      call assignvec(hex8(1,6),wedgever(1,14))
      call assignvec(hex8(1,7),wedgever(1,6))
      call assignvec(hex8(1,8),wedgever(1,15))

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0
      call assignvec(cr1(1,1),wedgever(1,10))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,2))
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,1))	  

      call quick_quater_point_from_curve(ver(1),wedgever(1,12),wedgeface(1,2),wedgever(1,11))
      call assignvec(edge12(1,2),ver(1)) 
      call quick_quater_point_from_curve(ver(1),wedgever(1,12),wedgeface(1,3),wedgever(1,10))
      call assignvec(edge12(1,3),ver(1)) 

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,11))
      call assignvec(cr1(1,2),wedgecen(1,1))
      call assignvec(cr1(1,3),wedgeface(1,3)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,4))	  

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,4))
      call assignvec(cr1(1,2),wedgeface(1,5))
      call assignvec(cr1(1,3),wedgever(1,14)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,5))	  

      call quick_quater_point_from_curve(ver(1),wedgever(1,6),wedgever(1,14),wedgever(1,5))
      call assignvec(edge12(1,6),ver(1)) 
      call quick_quater_point_from_curve(ver(1),wedgever(1,6),wedgever(1,15),wedgever(1,4))
      call assignvec(edge12(1,7),ver(1)) 

      x(1) = 0.0
      x(2) = 2.0/3.0
      x(3) = 1.0	  
      call assignvec(cr1(1,1),wedgever(1,5))
      call assignvec(cr1(1,2),wedgeface(1,5))
      call assignvec(cr1(1,3),wedgever(1,15)) 
      beta = 5.0/6.0
      call curve(cr1,x,beta,edge12(1,8))	  

      call quick_quater_point_from_curve(ver(1),wedgeface(1,5),wedgecen(1,1),wedgeface(1,4))
      call assignvec(edge12(1,9),ver(1))
      call quick_quater_point_from_curve(ver(1),wedgever(1,14),wedgeface(1,2),wedgever(1,8))
      call assignvec(edge12(1,10),ver(1)) 
      call quick_quater_point_from_curve(ver(1),wedgever(1,6),wedgever(1,12),wedgever(1,3))
      call assignvec(edge12(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),wedgever(1,15),wedgeface(1,3),wedgever(1,9))
      call assignvec(edge12(1,12),ver(1)) 


      call hex20tohex27(hexver(1,1,6),hex8(1,1),edge12(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 

      return
      end	 
!--------------------------------------------------------------------------s

!--------------------------------------------------------------------  
      subroutine ehexto8hex_quadratic(hexver,ehexver,hexss,ehexss,hasnrh) 
!
      real*8 hexver(3,27,8)
      real*8 ehexver(3,20)
      real*8 hexface(3,6)
      real*8 hexcen(3) ! tet vol center
      real*8 subhexver(3,20)
      real*8 hex8(3,8),edge12(3,12) ! edge12 here follows nek definition of edge 

      integer tetss(4),wedgess(5),ehexss(6),hexss(6,8)

      real*8 quad(3,8),fc(3),ver(3),fc2(3),tvec1(3),tvec2(3)
 
      logical hasnrh,ifnrh

      do i=1,6*8
      hexss(i,1)=0
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
      !call average2vec(hexface(1,1),fc,fc2)
      call assignvec(hexface(1,1),fc(1))

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
      !call average2vec(hexface(1,2),fc,fc2)
      call assignvec(hexface(1,2),fc(1))

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
      !call average2vec(hexface(1,3),fc,fc2)
	  call assignvec(hexface(1,3),fc(1))

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
      !call average2vec(hexface(1,4),fc,fc2)
      call assignvec(hexface(1,4),fc(1))

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
      !call average2vec(hexface(1,5),fc,fc2)
      call assignvec(hexface(1,5),fc(1))

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
      !call average2vec(hexface(1,6),fc,fc2)
      call assignvec(hexface(1,6),fc(1))

! get hex center
      call average6vec(hexcen(1),hexface(1,1),hexface(1,2), &
      hexface(1,3),hexface(1,4),hexface(1,5),hexface(1,6))

!  hex 1
      hexss(1,1) = ehexss(1)
      hexss(4,1) = ehexss(4)
      hexss(5,1) = ehexss(5)

      call assignvec(hex8(1,1),ehexver(1,1)) 
      call assignvec(hex8(1,2),ehexver(1,9)) 
      call assignvec(hex8(1,3),hexface(1,5)) 
      call assignvec(hex8(1,4),ehexver(1,12))
      call assignvec(hex8(1,5),ehexver(1,13)) 
      call assignvec(hex8(1,6),hexface(1,1)) 
      call assignvec(hex8(1,7),hexcen(1))
      call assignvec(hex8(1,8),hexface(1,4))

      call quick_quater_point_from_curve(ver(1),ehexver(1,1),ehexver(1,9),ehexver(1,2))
      call assignvec(edge12(1,1),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,9),hexface(1,5),ehexver(1,11))
      call assignvec(edge12(1,2),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,12),hexface(1,5),ehexver(1,10))
      call assignvec(edge12(1,3),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,1),ehexver(1,12),ehexver(1,4))
      call assignvec(edge12(1,4),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,13),hexface(1,1),ehexver(1,14))
      call assignvec(edge12(1,5),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,1),hexcen(1),hexface(1,3))
      call assignvec(edge12(1,6),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,4),hexcen(1),hexface(1,2))
      call assignvec(edge12(1,7),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,13),hexface(1,4),ehexver(1,16))
      call assignvec(edge12(1,8),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,1),ehexver(1,13),ehexver(1,5))
      call assignvec(edge12(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,9),hexface(1,1),ehexver(1,17))
      call assignvec(edge12(1,10),ver(1))
      call quick_quater_point_from_curve(ver(1),hexface(1,5),hexcen(1),hexface(1,6))
      call assignvec(edge12(1,11),ver(1)) 
     ! call quick_quater_point_from_curve(ver(1),ehexver(1,2),hexface(1,4),ehexver(1,20)) ! error !
      call quick_quater_point_from_curve(ver(1),ehexver(1,12),hexface(1,4),ehexver(1,20))
      call assignvec(edge12(1,12),ver(1)) 

      call hex20tohex27(hexver(1,1,1),hex8(1,1),edge12(1,1))  ! edge should follow nek edge order, not exo edge order
      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 

!  hex 2
      hexss(1,2) = ehexss(1)
      hexss(2,2) = ehexss(2)
      hexss(5,2) = ehexss(5)

      call assignvec(hex8(1,1),ehexver(1,9)) 
      call assignvec(hex8(1,2),ehexver(1,2)) 
      call assignvec(hex8(1,3),ehexver(1,10)) 
      call assignvec(hex8(1,4),hexface(1,5))
      call assignvec(hex8(1,5),hexface(1,1)) 
      call assignvec(hex8(1,6),ehexver(1,14)) 
      call assignvec(hex8(1,7),hexface(1,2)) 
      call assignvec(hex8(1,8),hexcen(1))

      call quick_quater_point_from_curve(ver(1),ehexver(1,2),ehexver(1,9),ehexver(1,1))
      call assignvec(edge12(1,1),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,2),ehexver(1,10),ehexver(1,3))
      call assignvec(edge12(1,2),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,10),hexface(1,5),ehexver(1,12))
      call assignvec(edge12(1,3),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,9),hexface(1,5),ehexver(1,11))
      call assignvec(edge12(1,4),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,14),hexface(1,1),ehexver(1,13))
      call assignvec(edge12(1,5),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,14),hexface(1,2),ehexver(1,15))
      call assignvec(edge12(1,6),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,2),hexcen(1),hexface(1,4))
      call assignvec(edge12(1,7),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,1),hexcen(1),hexface(1,3))
      call assignvec(edge12(1,8),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,9),hexface(1,1),ehexver(1,17))
      call assignvec(edge12(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,2),ehexver(1,14),ehexver(1,6))
      call assignvec(edge12(1,10),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,10),hexface(1,2),ehexver(1,18))
      call assignvec(edge12(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,5),hexcen(1),hexface(1,6))
      call assignvec(edge12(1,12),ver(1)) 

      call hex20tohex27(hexver(1,1,2),hex8(1,1),edge12(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
!  hex 3
      hexss(2,3) = ehexss(2)
      hexss(3,3) = ehexss(3)
      hexss(5,3) = ehexss(5)

      call assignvec(hex8(1,1),hexface(1,5)) 
      call assignvec(hex8(1,2),ehexver(1,10)) 
      call assignvec(hex8(1,3),ehexver(1,3)) 
      call assignvec(hex8(1,4),ehexver(1,11))
      call assignvec(hex8(1,5),hexcen(1)) 
      call assignvec(hex8(1,6),hexface(1,2)) 
      call assignvec(hex8(1,7),ehexver(1,15))
      call assignvec(hex8(1,8),hexface(1,3)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,10),hexface(1,5),ehexver(1,12))
      call assignvec(edge12(1,1),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,3),ehexver(1,10),ehexver(1,2))
      call assignvec(edge12(1,2),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,3),ehexver(1,11),ehexver(1,4))
      call assignvec(edge12(1,3),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,11),hexface(1,5),ehexver(1,9))
      call assignvec(edge12(1,4),ver(1)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,2),hexcen(1),hexface(1,4))
      call assignvec(edge12(1,5),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,15),hexface(1,2),ehexver(1,14))
      call assignvec(edge12(1,6),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,15),hexface(1,3),ehexver(1,16))
      call assignvec(edge12(1,7),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,3),hexcen(1),hexface(1,1))
      call assignvec(edge12(1,8),ver(1)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,5),hexcen(1),hexface(1,6))
      call assignvec(edge12(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,10),hexface(1,2),ehexver(1,18))
      call assignvec(edge12(1,10),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,3),ehexver(1,15),ehexver(1,7))
      call assignvec(edge12(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,11),hexface(1,3),ehexver(1,19))
      call assignvec(edge12(1,12),ver(1)) 

      call hex20tohex27(hexver(1,1,3),hex8(1,1),edge12(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 


!  hex 4
      hexss(3,4) = ehexss(3)
      hexss(4,4) = ehexss(4)
      hexss(5,4) = ehexss(5)

      call assignvec(hex8(1,1),ehexver(1,12)) 
      call assignvec(hex8(1,2),hexface(1,5)) 
      call assignvec(hex8(1,3),ehexver(1,11)) 
      call assignvec(hex8(1,4),ehexver(1,4))
      call assignvec(hex8(1,5),hexface(1,4)) 
      call assignvec(hex8(1,6),hexcen(1)) 
      call assignvec(hex8(1,7),hexface(1,3)) 
      call assignvec(hex8(1,8),ehexver(1,16)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,12),hexface(1,5),ehexver(1,10))
      call assignvec(edge12(1,1),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,11),hexface(1,5),ehexver(1,9))
      call assignvec(edge12(1,2),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,4),ehexver(1,11),ehexver(1,3))
      call assignvec(edge12(1,3),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,4),ehexver(1,12),ehexver(1,1))
      call assignvec(edge12(1,4),ver(1)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,4),hexcen(1),hexface(1,2))
      call assignvec(edge12(1,5),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,3),hexcen(1),hexface(1,1))
      call assignvec(edge12(1,6),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,16),hexface(1,3),ehexver(1,15))
      call assignvec(edge12(1,7),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,16),hexface(1,4),ehexver(1,13))
      call assignvec(edge12(1,8),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,12),hexface(1,4),ehexver(1,20))
      call assignvec(edge12(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,5),hexcen(1),hexface(1,6))
      call assignvec(edge12(1,10),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,11),hexface(1,3),ehexver(1,19))
      call assignvec(edge12(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,4),ehexver(1,16),ehexver(1,8))
      call assignvec(edge12(1,12),ver(1)) 


      call hex20tohex27(hexver(1,1,4),hex8(1,1),edge12(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 
!  hex 5
      hexss(1,5) = ehexss(1)
      hexss(4,5) = ehexss(4)
      hexss(6,5) = ehexss(6)


      call assignvec(hex8(1,1),ehexver(1,13)) 
      call assignvec(hex8(1,2),hexface(1,1)) 
      call assignvec(hex8(1,3),hexcen(1)) 
      call assignvec(hex8(1,4),hexface(1,4))
      call assignvec(hex8(1,5),ehexver(1,5)) 
      call assignvec(hex8(1,6),ehexver(1,17)) 
      call assignvec(hex8(1,7),hexface(1,6)) 
      call assignvec(hex8(1,8),ehexver(1,20)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,13),hexface(1,1),ehexver(1,14))
      call assignvec(edge12(1,1),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,1),hexcen(1),hexface(1,3))
      call assignvec(edge12(1,2),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,4),hexcen(1),hexface(1,2))
      call assignvec(edge12(1,3),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,13),hexface(1,4),ehexver(1,16))
      call assignvec(edge12(1,4),ver(1))

      call quick_quater_point_from_curve(ver(1),ehexver(1,5),ehexver(1,17),ehexver(1,6))
      call assignvec(edge12(1,5),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,17),hexface(1,6),ehexver(1,19))
      call assignvec(edge12(1,6),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,20),hexface(1,6),ehexver(1,18))
      call assignvec(edge12(1,7),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,5),ehexver(1,20),ehexver(1,8))
      call assignvec(edge12(1,8),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,5),ehexver(1,13),ehexver(1,1))
      call assignvec(edge12(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,17),hexface(1,1),ehexver(1,9))
      call assignvec(edge12(1,10),ver(1))
      call quick_quater_point_from_curve(ver(1),hexface(1,6),hexcen(1),hexface(1,5))
      call assignvec(edge12(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,20),hexface(1,4),ehexver(1,12))
      call assignvec(edge12(1,12),ver(1)) 

      call hex20tohex27(hexver(1,1,5),hex8(1,1),edge12(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 

!  hex 6
      hexss(1,6) = ehexss(1)
      hexss(2,6) = ehexss(2)
      hexss(6,6) = ehexss(6)

      call assignvec(hex8(1,1),hexface(1,1)) 
      call assignvec(hex8(1,2),ehexver(1,14)) 
      call assignvec(hex8(1,3),hexface(1,2)) 
      call assignvec(hex8(1,4),hexcen(1))
      call assignvec(hex8(1,5),ehexver(1,17)) 
      call assignvec(hex8(1,6),ehexver(1,6)) 
      call assignvec(hex8(1,7),ehexver(1,18)) 
      call assignvec(hex8(1,8),hexface(1,6)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,14),hexface(1,1),ehexver(1,13))
      call assignvec(edge12(1,1),ver(1)) 
      !call quick_quater_point_from_curve(ver(1),ehexver(1,14),hexface(1,1),ehexver(1,15)) !! error here!
      call quick_quater_point_from_curve(ver(1),ehexver(1,14),hexface(1,2),ehexver(1,15))  ! this is the correct version!
      call assignvec(edge12(1,2),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,2),hexcen(1),hexface(1,4))
      call assignvec(edge12(1,3),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,1),hexcen(1),hexface(1,3))
      call assignvec(edge12(1,4),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,6),ehexver(1,17),ehexver(1,5))
      call assignvec(edge12(1,5),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,6),ehexver(1,18),ehexver(1,7))
      call assignvec(edge12(1,6),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,18),hexface(1,6),ehexver(1,20))
      call assignvec(edge12(1,7),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,17),hexface(1,6),ehexver(1,19))
      call assignvec(edge12(1,8),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,17),hexface(1,1),ehexver(1,9))
      call assignvec(edge12(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,6),ehexver(1,14),ehexver(1,2))
      call assignvec(edge12(1,10),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,18),hexface(1,2),ehexver(1,10))
      call assignvec(edge12(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,6),hexcen(1),hexface(1,5))
      call assignvec(edge12(1,12),ver(1)) 


      call hex20tohex27(hexver(1,1,6),hex8(1,1),edge12(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 

!  hex 7
      hexss(2,7) = ehexss(2)
      hexss(3,7) = ehexss(3)
      hexss(6,7) = ehexss(6)


      call assignvec(hex8(1,1),hexcen(1)) 
      call assignvec(hex8(1,2),hexface(1,2))
      call assignvec(hex8(1,3),ehexver(1,15)) 
      call assignvec(hex8(1,4),hexface(1,3))
      call assignvec(hex8(1,5),hexface(1,6))
      call assignvec(hex8(1,6),ehexver(1,18)) 
      call assignvec(hex8(1,7),ehexver(1,7))
      call assignvec(hex8(1,8),ehexver(1,19)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,2),hexcen(1),hexface(1,4))
      call assignvec(edge12(1,1),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,15),hexface(1,2),ehexver(1,14))
      call assignvec(edge12(1,2),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,15),hexface(1,3),ehexver(1,16))
      call assignvec(edge12(1,3),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,3),hexcen(1),hexface(1,1))
      call assignvec(edge12(1,4),ver(1)) 



      call quick_quater_point_from_curve(ver(1),ehexver(1,18),hexface(1,6),ehexver(1,20))
      call assignvec(edge12(1,5),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,7),ehexver(1,18),ehexver(1,6))
      call assignvec(edge12(1,6),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,7),ehexver(1,19),ehexver(1,8))
      call assignvec(edge12(1,7),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,19),hexface(1,6),ehexver(1,17))
      call assignvec(edge12(1,8),ver(1)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,6),hexcen(1),hexface(1,5))
      call assignvec(edge12(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,18),hexface(1,2),ehexver(1,10))
      call assignvec(edge12(1,10),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,7),ehexver(1,15),ehexver(1,3))
      call assignvec(edge12(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,19),hexface(1,3),ehexver(1,11))
      call assignvec(edge12(1,12),ver(1)) 

      call hex20tohex27(hexver(1,1,7),hex8(1,1),edge12(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 

!  hex 8
      hexss(3,8) = ehexss(3)
      hexss(4,8) = ehexss(4)
      hexss(6,8) = ehexss(6)


      call assignvec(hex8(1,1),hexface(1,4))
      call assignvec(hex8(1,2),hexcen(1))
      call assignvec(hex8(1,3),hexface(1,3)) 
      call assignvec(hex8(1,4),ehexver(1,16)) 
      call assignvec(hex8(1,5),ehexver(1,20))
      call assignvec(hex8(1,6),hexface(1,6))
      call assignvec(hex8(1,7),ehexver(1,19))
      call assignvec(hex8(1,8),ehexver(1,8)) 

      call quick_quater_point_from_curve(ver(1),hexface(1,4),hexcen(1),hexface(1,2))
      call assignvec(edge12(1,1),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,3),hexcen(1),hexface(1,1))
      call assignvec(edge12(1,2),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,16),hexface(1,3),ehexver(1,15))
      call assignvec(edge12(1,3),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,16),hexface(1,4),ehexver(1,13))
      call assignvec(edge12(1,4),ver(1)) 

      call quick_quater_point_from_curve(ver(1),ehexver(1,20),hexface(1,6),ehexver(1,18))
      call assignvec(edge12(1,5),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,19),hexface(1,6),ehexver(1,17))
      call assignvec(edge12(1,6),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,8),ehexver(1,19),ehexver(1,7))
      call assignvec(edge12(1,7),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,8),ehexver(1,20),ehexver(1,5))
      call assignvec(edge12(1,8),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,20),hexface(1,4),ehexver(1,12))
      call assignvec(edge12(1,9),ver(1)) 
      call quick_quater_point_from_curve(ver(1),hexface(1,6),hexcen(1),hexface(1,5))
      call assignvec(edge12(1,10),ver(1))
      call quick_quater_point_from_curve(ver(1),ehexver(1,19),hexface(1,3),ehexver(1,11))
      call assignvec(edge12(1,11),ver(1)) 
      call quick_quater_point_from_curve(ver(1),ehexver(1,8),ehexver(1,16),ehexver(1,4))
      call assignvec(edge12(1,12),ver(1)) 


      call hex20tohex27(hexver(1,1,8),hex8(1,1),edge12(1,1))

      ifnrh = .false.  
      call nek_check_non_right_hand_per_element(hex8,ifnrh)
      hasnrh = (hasnrh.or.ifnrh) ! 

      return
      end 
! -------------------------------------------------------------------
      subroutine hex20tohex27(hex27,hex8,edge12)
!  convert hex8 coordinates to hex27 coordinates in nek.
!  
      real*8 hex8(3,8),edge12(3,12)
      real*8 hex27(3,27)
      real*8 tempvec(3,3)

!   hex27 vert 1
      call assignvec(hex27(1,1),hex8(1,1))
!   hex27 vert 2
      call assignvec(hex27(1,2),edge12(1,1))
!   hex27 vert 3
      call assignvec(hex27(1,3),hex8(1,2))
!   hex27 vert 4
      call assignvec(hex27(1,4),edge12(1,4))
!   hex27 vert 5
      call average4vec(tempvec(1,1),hex8(1,1), &
      hex8(1,2),hex8(1,3),hex8(1,4))
      call assignvec(hex27(1,5),tempvec(1,1))
!   hex27 vert 6
      call assignvec(hex27(1,6),edge12(1,2))
!   hex27 vert 7
      call assignvec(hex27(1,7),hex8(1,4))
!   hex27 vert 8
      call assignvec(hex27(1,8),edge12(1,3))
!   hex27 vert 9
      call assignvec(hex27(1,9),hex8(1,3))

!   hex27 vert 10
      call assignvec(hex27(1,10),edge12(1,9))  
!   hex27 vert 11
      call average4vec(tempvec(1,1),hex8(1,1), &
      hex8(1,2),hex8(1,5),hex8(1,6))
      call assignvec(hex27(1,11),tempvec(1,1))
!   hex27 vert 12
      call assignvec(hex27(1,12),edge12(1,10))

!   hex27 vert 13
      call average4vec(tempvec(1,1),hex8(1,1), &
      hex8(1,4),hex8(1,5),hex8(1,8))
      call assignvec(hex27(1,13),tempvec(1,1))
!   hex27 vert 14
      call average4vec(tempvec(1,1),hex8(1,1), &
      hex8(1,4),hex8(1,5),hex8(1,8))
      call average4vec(tempvec(1,2),hex8(1,2), &
      hex8(1,3),hex8(1,6),hex8(1,7))
      call average2vec(tempvec(1,3),tempvec(1,1),tempvec(1,2))
      call assignvec(hex27(1,14),tempvec(1,3))
!   hex27 vert 15
      call average4vec(tempvec(1,1),hex8(1,2), & 
      hex8(1,3),hex8(1,6),hex8(1,7))
      call assignvec(hex27(1,15),tempvec(1,1))

!   hex27 vert 16
      call assignvec(hex27(1,16),edge12(1,12))
!   hex27 vert 17
      call average4vec(tempvec(1,1),hex8(1,3), &
      hex8(1,4),hex8(1,7),hex8(1,8))
      call assignvec(hex27(1,17),tempvec(1,1))
!   hex27 vert 18
      call assignvec(hex27(1,18),edge12(1,11))

!   hex27 vert 19
      call assignvec(hex27(1,19),hex8(1,5))
!   hex27 vert 20
      call assignvec(hex27(1,20),edge12(1,5))	  
!   hex27 vert 21
      call assignvec(hex27(1,21),hex8(1,6))

!   hex27 vert 22
      call assignvec(hex27(1,22),edge12(1,8))	  
!   hex27 vert 23
      call average4vec(tempvec(1,1),hex8(1,5), &
      hex8(1,6),hex8(1,7),hex8(1,8))
      call assignvec(hex27(1,23),tempvec(1,1))
!   hex27 vert 24
      call assignvec(hex27(1,24),edge12(1,6))
	  
!   hex27 vert 25
      call assignvec(hex27(1,25),hex8(1,8))
!   hex27 vert 26
      call assignvec(hex27(1,26),edge12(1,7))	  
!   hex27 vert 27
      call assignvec(hex27(1,27),hex8(1,7))

      return
      end
! -------------------------------------------------------------------