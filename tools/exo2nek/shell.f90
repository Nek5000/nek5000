!-----------------------------------------------------------------
      subroutine shell_wall()
!
! create shell solid elements
! 
      use SIZE
      use OCTUPLESIZE

	  real quad8p(3,8)
	  real quad8v(3,8)
      real hex27(3,27),hex27_new(3,27)
      real ss4(4),ss6(6)
      integer ifc,iel,ifc_new,ivert

      real thickness_ratio
      integer shell_ss,estot

      if (eftot.ne.etot) write(6,*) 'ERROR: this feature is only for fluid-only mesh'

      write(6,*) 'Please input sideset that you want to shell'
      read (5,*) shell_ss
      write(6,*) 'Please input thickness ratio: (0-1) '
      read (5,*) thickness_ratio

      estot = 0 

      do iel =1,etot
       do ifc = 1,6
           if(bc(5,ifc,iel).eq.shell_ss) then   
            estot = estot + 1
           endif
       enddo
      enddo

      etot_new = etot + estot
      etot = etot_new
      num_elem = etot_new 

      allocate ( xm2                (3,3,3,etot_new)      )
      allocate ( ym2                (3,3,3,etot_new)      )
      allocate ( zm2                (3,3,3,etot_new)      )

      call rzero(xm2,3*3*3*etot_new)
      call rzero(ym2,3*3*3*etot_new)
      call rzero(zm2,3*3*3*etot_new)

      allocate   (cbc2    (2*num_dim,      etot_new) )
      allocate   (bc2     (5,2*num_dim,    etot_new) )
      call rzero (bc2,5*2*num_dim*etot_new)
      call blank (cbc2,3*2*num_dim*etot_new)

! copy fluid elements
      do iel =1,eftot

       do ivert = 1,27
       xm2(ivert,1,1,iel) = xm1(ivert,1,1,iel)
       ym2(ivert,1,1,iel) = ym1(ivert,1,1,iel)
       zm2(ivert,1,1,iel) = zm1(ivert,1,1,iel)
       enddo
       
	   do ifc=1,6
       bc2(5,ifc,iel) = bc(5,ifc,iel)
       bc2(4,ifc,iel) = bc(4,ifc,iel)
       cbc2(ifc,iel) = cbc(ifc,iel) 
       enddo

      enddo

! add shell elements
      estot = 0 
      do iel =1,eftot
	  
       do ivert = 1,27
       hex27(1,ivert) = xm1(ivert,1,1,iel)
       hex27(2,ivert) = ym1(ivert,1,1,iel)
       hex27(3,ivert) = zm1(ivert,1,1,iel)
       enddo

       do ifc = 1,6
       ss6(ifc) = bc(5,ifc,iel)
       enddo

       do ifc = 1,6
           if(bc(5,ifc,iel).eq.shell_ss) then   
            estot = estot + 1

       !pack quad8 info
       call pack_quad8(ifc,hex27,quad8p,quad8v)

       !get edge_sideset()
       call edge_sideset(ifc,ss6,ss4)

       !shelling
       call shell_one_layer(quad8p,quad8v,hex27_new,thickness_ratio)

           iel_new = eftot+estot

           do ivert = 1,27
           xm2(ivert,1,1,iel_new) = hex27_new(1,ivert)
           ym2(ivert,1,1,iel_new) = hex27_new(2,ivert)
           zm2(ivert,1,1,iel_new) = hex27_new(3,ivert)
           enddo

           do ifc_new = 1,4
           !bc2(5,ifc_new,iel_new) = ss4(ifc_new)
		   if (ss4(ifc_new).gt.0.0) then
           cbc2(ifc_new,iel_new) = 'EXO'
           bc2(5,ifc_new,iel_new) = ss4(ifc_new)+100
           endif
           enddo

           cbc2(5,iel_new) = 'EXO'
           cbc2(6,iel_new) = 'EXO'
           bc2(5,5,iel_new) = shell_ss+100 ! bottom
           bc2(5,6,iel_new) = shell_ss+200 ! top

           endif
       enddo

      enddo

! copy new mesh info

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
      subroutine pack_quad8(ifc,hex27,quad8p,quad8v)

	  real quad8p(3,8),quad8p2(3,8)
	  real quad8v(3,8)
      real hex27(3,27)
      real mag
      integer ifc,iv


! node and face conversion (it works at least for cubit):
      integer face_to_quad(8,6)                   ! hex27 to nek numbering
      data    face_to_quad                     &
            / 1,3,21,19,2,12,20,10   &     
            , 3,9,27,21,6,18,24,12   &
            , 9,7,25,27,8,16,26,18   &
            , 7,1,19,25,4,10,22,16   &     
            , 1,7,9,3,4,8,6,2        &
            , 19,21,27,25,20,24,26,22/

      integer face_to_quad2(8,6)                   ! hex27 to nek numbering
      data    face_to_quad2                     &
            / 7,9,27,25,8,18,26,16     &     
            , 1,7,25,19,4,16,22,10     &
            , 3,1,19,21,2,10,20,12     &
            , 9,3,21,27,6,12,24,18     &     
            , 19,25,27,21,22,26,24,20  &
            , 1,3,9,7,2,6,8,4 /

       do iv = 1,8
        quad8p(1,iv) = hex27(1,face_to_quad(iv,ifc))
        quad8p(2,iv) = hex27(2,face_to_quad(iv,ifc))
        quad8p(3,iv) = hex27(3,face_to_quad(iv,ifc))
 
        quad8p2(1,iv) = hex27(1,face_to_quad2(iv,ifc))
        quad8p2(2,iv) = hex27(2,face_to_quad2(iv,ifc))
        quad8p2(3,iv) = hex27(3,face_to_quad2(iv,ifc))
 
        quad8v(1,iv) =  quad8p(1,iv) - quad8p2(1,iv)
        quad8v(2,iv) =  quad8p(2,iv) - quad8p2(2,iv)
        quad8v(3,iv) =  quad8p(3,iv) - quad8p2(3,iv)

!        mag = sqrt(quad8v(1,iv)**2.0+  &
!      quad8v(2,iv)**2.0+quad8v(3,iv)**2.0)
!
!        quad8v(1,iv) =  quad8v(1,iv)/mag
!        quad8v(2,iv) =  quad8v(2,iv)/mag
!        quad8v(3,iv) =  quad8v(3,iv)/mag

       enddo


      return
      end
!-------------------------------------------------------------------- 
!-------------------------------------------------------------------- 
      subroutine edge_sideset(ifc,ss6,ss4)
	  real ss6(6),ss4(4)
      integer ifc

      integer face_edge_sideset(4,6)                   ! hex27 to nek numbering
      data    face_edge_sideset                     &
            / 5,2,6,4 &     
            , 5,3,6,1 &
            , 5,4,6,2 &
            , 5,1,6,3 &      
            , 4,3,2,1 &
            , 1,2,3,4 /

      do ie = 1,4
       ss4(ie)= ss6(face_edge_sideset(ie,ifc))
      enddo

      return
      end
!-------------------------------------------------------------------- 
!-------------------------------------------------------------------- 
      subroutine shell_one_layer(quad8p,quad8v,hex27_new,ratio)

	  real quad8p(3,8),quad8p2(3,8),quad8p3(3,8)
	  real quad8v(3,8)
      real hex27_new(3,27)
      real ratio
      logical lefthand

      call rzero(hex27_new,3*27)

      do iv = 1,8
      quad8p2(1,iv) = quad8p(1,iv) + quad8v(1,iv)*ratio/2.0
      quad8p2(2,iv) = quad8p(2,iv) + quad8v(2,iv)*ratio/2.0
      quad8p2(3,iv) = quad8p(3,iv) + quad8v(3,iv)*ratio/2.0
      enddo

      do iv = 1,8
      quad8p3(1,iv) = quad8p(1,iv) + quad8v(1,iv)*ratio
      quad8p3(2,iv) = quad8p(2,iv) + quad8v(2,iv)*ratio
      quad8p3(3,iv) = quad8p(3,iv) + quad8v(3,iv)*ratio
      enddo

! assign hex27 points.
      call assignvec(hex27_new(1,1),quad8p(1,1))
      call assignvec(hex27_new(1,2),quad8p(1,5))
      call assignvec(hex27_new(1,3),quad8p(1,2))
      call assignvec(hex27_new(1,4),quad8p(1,8))
      !call assignvec(hex27_new(1,5),quad8p(1,1))
      call assignvec(hex27_new(1,6),quad8p(1,6))
      call assignvec(hex27_new(1,7),quad8p(1,4))
      call assignvec(hex27_new(1,8),quad8p(1,7))
      call assignvec(hex27_new(1,9),quad8p(1,3))

      call assignvec(hex27_new(1,10),quad8p2(1,1))
      !call assignvec(hex27_new(1,11),quad8p2(1,1))
      call assignvec(hex27_new(1,12),quad8p2(1,2))
      !call assignvec(hex27_new(1,13),quad8p2(1,1))
      !call assignvec(hex27_new(1,14),quad8p2(1,1))
      !call assignvec(hex27_new(1,15),quad8p2(1,1))
      call assignvec(hex27_new(1,16),quad8p2(1,4))
      !call assignvec(hex27_new(1,17),quad8p2(1,1))
      call assignvec(hex27_new(1,18),quad8p2(1,3))

      call assignvec(hex27_new(1,19),quad8p3(1,1))
      call assignvec(hex27_new(1,20),quad8p3(1,5))
      call assignvec(hex27_new(1,21),quad8p3(1,2))
      call assignvec(hex27_new(1,22),quad8p3(1,8))
      !call assignvec(hex27_new(1,23),quad8p3(1,1))
      call assignvec(hex27_new(1,24),quad8p3(1,6))
      call assignvec(hex27_new(1,25),quad8p3(1,4))
      call assignvec(hex27_new(1,26),quad8p3(1,7))
      call assignvec(hex27_new(1,27),quad8p3(1,3))

! check hex27 elements, if non-right-hand, report.
      lefthand = .false.
      call nek_check_non_right_hand_hex27(lefthand,hex27_new)
      if (lefthand) then 
       write(6,*) 'ERROR: left hand element is created when shelling'
       write(6,*) 'please reduce shell_thickness'
      endif

	  
      return
      end
!-------------------------------------------------------------------- 
