!--------------------------------------------------------------------
      subroutine fix_left_hand_elements_3d
! check if there is non-right hand elements (3D)
! because if mesh is from ICEM, and mirror operation is made in ICEM,
! the exported exo file will contain non-right hand elements.
! this subroutine will:
! 1. check right-hand
! 2. fix if not
      use SIZE
      integer iel
      logical lefthand
      character nek_check
	  
      write(6,*) 'do you want to fix left-hand element (y/n)'
      read (5,*) nek_check
 
      if (nek_check.eq.'y') then
        do iel=1,num_elem
          lefthand = .FALSE.
          call check_if_left_hand(lefthand,iel)
		  if (lefthand) call fix_left_hand(iel)
        enddo
      endif
      return
      end
!--------------------------------------------------------------------
      subroutine r_or_l_detect_for_quad(iel,rfflag)
! detect if iquad is right-hand or left-hand, elements
! borrowed from gmsh2nek
      use SIZE
      include 'exodusII.inc'

      integer iel, rfflag
      integer node(4)
      real vec12(3),vec14(3),cz
	  
      nvert = num_nodes_per_elem(1)
	  
      do inode = 1,4
         node(inode) = connect(nvert*(iel-1)+inode)
      enddo

      vec12(1) = x_exo(node(2)) - x_exo(node(1))
      vec12(2) = y_exo(node(2)) - y_exo(node(1))

      vec14(1) = x_exo(node(4)) - x_exo(node(1))
      vec14(2) = y_exo(node(4)) - y_exo(node(1))
	  
      cz = vec12(1)*vec14(2) - vec12(2)*vec14(1)

      if(cz.gt.0.0) rfflag = 0 ! right hand element
      if(cz.lt.0.0) rfflag = 1 ! left hand element

      return
      end
!--------------------------------------------------------------------
      subroutine right_hand_check(ne_nrh)
! check if there is non-right hand elements (3D)
! because if mesh is from ICEM, and mirror operation is made in ICEM,
! the exported exo file will contain non-right hand elements.
! this subroutine will:
! 1. check right-hand
! 2. fix if not
      use SIZE
      integer iel,ne_nrh
      logical ifnonrighthand
      character nek_check

      write(6,*) 'performing non-right-hand check'
	  
!      do iel=1,num_elem
!         !write(6,*) 'performing non-right-hand check on element ',iel
!         call check_if_non_right(ifnonrighthand,iel)
!         if (ifnonrighthand) call fix_if_non_right(iel)		 
!      enddo

!      write(6,*) 'done: non-right-hand check'
	  
!      write(6,*) 'using nek-method to do non-right-hand check? (y/n)'
!      read (5,*) nek_check
 
!      if (nek_check.eq.'y') then
!        do iel=1,num_elem
!          if (num_dim.eq.2) then
!          call nek_check_non_right_hand_2d(iel)
!          else 
!          call nek_check_non_right_hand(iel)
!          endif
!       enddo
!      endif

       ne_nrh = 0
        do iel=1,num_elem
          if (num_dim.eq.2) then
          call nek_check_non_right_hand_2d(iel,ne_nrh)
          else 
          call nek_check_non_right_hand(iel,ne_nrh)
          endif
       enddo
	  
      if (ne_nrh.gt.0) then
      write(6,*) 'WARNING: non-right-hand elements detected!'
      write(6,*) 'number of non-right-hand elements: ', ne_nrh
      endif
	  
      return 
      end
!--------------------------------------------------------------------
      subroutine nek_check_non_right_hand_2d(iel,ne_nrh)
      use SIZE
      logical ifnonrighthand
      integer iel,ne_nrh
      integer quad4_to_nek_quad9_vertex(4)
      data quad4_to_nek_quad9_vertex /1,3,7,9/ ! for nek non-right-hand element check
      real XYZ(2,4)
      real C1,C2,C3,C4

      do iver = 1,4
       XYZ(1,iver) = xm1(quad4_to_nek_quad9_vertex(iver),1,1,iel)
       XYZ(2,iver) = ym1(quad4_to_nek_quad9_vertex(iver),1,1,iel)
      enddo

!
!        CRSS2D(A,B,O) = (A-O) X (B-O)
!
         C1=CRSS2D(XYZ(1,2),XYZ(1,3),XYZ(1,1))
         C2=CRSS2D(XYZ(1,4),XYZ(1,1),XYZ(1,2))
         C3=CRSS2D(XYZ(1,1),XYZ(1,4),XYZ(1,3))
         C4=CRSS2D(XYZ(1,3),XYZ(1,2),XYZ(1,4))
!
      IF (C1.LE.0.0.OR.C2.LE.0.0.OR. &
            C3.LE.0.0.OR.C4.LE.0.0 ) THEN
       !write(6,*) 'WARNINGb: Detected non-right-handed element.'
       !write(6,*) 'at location:',XYZ(1,1),',',XYZ(2,1)
        ne_nrh = ne_nrh + 1
      ENDIF

      return
      end
!!--------------------------------------------------------------------
      subroutine nek_check_non_right_hand(iel,ne_nrh)
      use SIZE
      logical ifnonrighthand
      integer iel,ne_nrh
      integer hex8_to_hex27_vertex(8)
      data hex8_to_hex27_vertex /1,3,7,9,19,21,25,27/ ! for nek non-right-hand element check
      real XYZ(3,8)
      real V1,V2,V3,V4,V5,V6,V7,V8

      do iver = 1,8
       XYZ(1,iver) = xm1(hex8_to_hex27_vertex(iver),1,1,iel)
       XYZ(2,iver) = ym1(hex8_to_hex27_vertex(iver),1,1,iel)
       XYZ(3,iver) = zm1(hex8_to_hex27_vertex(iver),1,1,iel)       
      enddo
	  
      V1= VOLUM0(XYZ(1,2),XYZ(1,3),XYZ(1,5),XYZ(1,1))
      V2= VOLUM0(XYZ(1,4),XYZ(1,1),XYZ(1,6),XYZ(1,2))
      V3= VOLUM0(XYZ(1,1),XYZ(1,4),XYZ(1,7),XYZ(1,3))
      V4= VOLUM0(XYZ(1,3),XYZ(1,2),XYZ(1,8),XYZ(1,4))
      V5=-VOLUM0(XYZ(1,6),XYZ(1,7),XYZ(1,1),XYZ(1,5))
      V6=-VOLUM0(XYZ(1,8),XYZ(1,5),XYZ(1,2),XYZ(1,6))
      V7=-VOLUM0(XYZ(1,5),XYZ(1,8),XYZ(1,3),XYZ(1,7))
      V8=-VOLUM0(XYZ(1,7),XYZ(1,6),XYZ(1,4),XYZ(1,8))

      if ((V1.LE.0.0).OR.(V2.LE.0.0).OR. &
       (V3.LE.0.0).OR.(V4.LE.0.0).OR. &
       (V5.LE.0.0).OR.(V6.LE.0.0).OR. &
       (V7.LE.0.0).OR.(V8.LE.0.0)) then
   
      !write(6,*) 'WARNINGb: Detected non-right-handed element.'
      !write(6,*) 'at location:',XYZ(1,1),',',XYZ(2,1),',',XYZ(3,1)
       ne_nrh = ne_nrh + 1
      endif

      return
      end
!!--------------------------------------------------------------------
      subroutine nek_check_non_right_hand_hex27(ifleft,hex27)
      use SIZE
      logical ifnonrighthand,ifleft 
      real hex27(3,27)
      integer iel,ne_nrh
      integer hex8_to_hex27_vertex(8)
      data hex8_to_hex27_vertex /1,3,7,9,19,21,25,27/ ! for nek non-right-hand element check
      real XYZ(3,8)
      real V1,V2,V3,V4,V5,V6,V7,V8

      do iver = 1,8
       XYZ(1,iver) = hex27(1,hex8_to_hex27_vertex(iver))
       XYZ(2,iver) = hex27(2,hex8_to_hex27_vertex(iver))
       XYZ(3,iver) = hex27(3,hex8_to_hex27_vertex(iver))
      enddo
	  
      V1= VOLUM0(XYZ(1,2),XYZ(1,3),XYZ(1,5),XYZ(1,1))
      V2= VOLUM0(XYZ(1,4),XYZ(1,1),XYZ(1,6),XYZ(1,2))
      V3= VOLUM0(XYZ(1,1),XYZ(1,4),XYZ(1,7),XYZ(1,3))
      V4= VOLUM0(XYZ(1,3),XYZ(1,2),XYZ(1,8),XYZ(1,4))
      V5=-VOLUM0(XYZ(1,6),XYZ(1,7),XYZ(1,1),XYZ(1,5))
      V6=-VOLUM0(XYZ(1,8),XYZ(1,5),XYZ(1,2),XYZ(1,6))
      V7=-VOLUM0(XYZ(1,5),XYZ(1,8),XYZ(1,3),XYZ(1,7))
      V8=-VOLUM0(XYZ(1,7),XYZ(1,6),XYZ(1,4),XYZ(1,8))

      if ((V1.LE.0.0).OR.(V2.LE.0.0).OR. &
       (V3.LE.0.0).OR.(V4.LE.0.0).OR. &
       (V5.LE.0.0).OR.(V6.LE.0.0).OR. &
       (V7.LE.0.0).OR.(V8.LE.0.0)) then
   
      !write(6,*) 'WARNINGb: Detected non-right-handed element.'
      !write(6,*) 'at location:',XYZ(1,1),',',XYZ(2,1),',',XYZ(3,1)
      ! ne_nrh = ne_nrh + 1
    
	   ifleft = .true.

      endif

      return
      end
!------------------------------------------------------------------------
      subroutine nek_check_non_right_hand_per_element(XYZorg,ifnonrighthand)
      use SIZE
      logical ifnonrighthand
!      integer iel
!      integer hex8_to_hex27_vertex(8)
!      data hex8_to_hex27_vertex /1,3,7,9,19,21,25,27/ ! for nek non-right-hand element check
      real*8 XYZ(3,8),XYZorg(3,8)
      real V1,V2,V3,V4,V5,V6,V7,V8

      do iver = 1,8
       XYZ(1,iver) = XYZorg(1,iver)
       XYZ(2,iver) = XYZorg(2,iver)
       XYZ(3,iver) = XYZorg(3,iver)
      enddo
	  
      ! swap 3-4
       XYZ(1,3) = XYZorg(1,4)
       XYZ(2,3) = XYZorg(2,4)
       XYZ(3,3) = XYZorg(3,4)
       XYZ(1,4) = XYZorg(1,3)
       XYZ(2,4) = XYZorg(2,3)
       XYZ(3,4) = XYZorg(3,3)
	   
      ! swap 7-8
       XYZ(1,7) = XYZorg(1,8)
       XYZ(2,7) = XYZorg(2,8)
       XYZ(3,7) = XYZorg(3,8)
       XYZ(1,8) = XYZorg(1,7)
       XYZ(2,8) = XYZorg(2,7)
       XYZ(3,8) = XYZorg(3,7)

	  
      ifnonrighthand = .false.
	  
      V1= VOLUM0(XYZ(1,2),XYZ(1,3),XYZ(1,5),XYZ(1,1))
      V2= VOLUM0(XYZ(1,4),XYZ(1,1),XYZ(1,6),XYZ(1,2))
      V3= VOLUM0(XYZ(1,1),XYZ(1,4),XYZ(1,7),XYZ(1,3))
      V4= VOLUM0(XYZ(1,3),XYZ(1,2),XYZ(1,8),XYZ(1,4))
      V5=-VOLUM0(XYZ(1,6),XYZ(1,7),XYZ(1,1),XYZ(1,5))
      V6=-VOLUM0(XYZ(1,8),XYZ(1,5),XYZ(1,2),XYZ(1,6))
      V7=-VOLUM0(XYZ(1,5),XYZ(1,8),XYZ(1,3),XYZ(1,7))
      V8=-VOLUM0(XYZ(1,7),XYZ(1,6),XYZ(1,4),XYZ(1,8))

      if ((V1.LE.0.0).OR.(V2.LE.0.0).OR. &
       (V3.LE.0.0).OR.(V4.LE.0.0).OR. &
       (V5.LE.0.0).OR.(V6.LE.0.0).OR. &
       (V7.LE.0.0).OR.(V8.LE.0.0)) then
   
      !write(6,*) 'WARNINGb: Detected non-right-handed element.'
      !write(6,*) 'at location:',XYZ(1,1),',',XYZ(2,1),',',XYZ(3,1)
      ifnonrighthand = .true.
      endif

      return
      end
!--------------------------------------------------------------------
      subroutine check_if_left_hand(lefthand,iel)
      use SIZE
      logical lefthand
      integer iel
      integer hex8_to_hex27_vertex(8)
      data hex8_to_hex27_vertex /1,3,9,7,19,21,27,25/
      real*8 hex8_vertex(3,8),vec12(3),vec14(3),vec15(3)
      real*8 vec1(3),AA,dot_prod

      do iver = 1,8
       hex8_vertex(1,iver) = xm1(hex8_to_hex27_vertex(iver),1,1,iel)
       hex8_vertex(2,iver) = ym1(hex8_to_hex27_vertex(iver),1,1,iel)
       hex8_vertex(3,iver) = zm1(hex8_to_hex27_vertex(iver),1,1,iel)       
      enddo
	  
      vec12(1) = hex8_vertex(1,2) - hex8_vertex(1,1)
      vec12(2) = hex8_vertex(2,2) - hex8_vertex(2,1)
      vec12(3) = hex8_vertex(3,2) - hex8_vertex(3,1)

      vec14(1) = hex8_vertex(1,4) - hex8_vertex(1,1)
      vec14(2) = hex8_vertex(2,4) - hex8_vertex(2,1)
      vec14(3) = hex8_vertex(3,4) - hex8_vertex(3,1)

      vec15(1) = hex8_vertex(1,5) - hex8_vertex(1,1)
      vec15(2) = hex8_vertex(2,5) - hex8_vertex(2,1)
      vec15(3) = hex8_vertex(3,5) - hex8_vertex(3,1)

      call cross_product(vec12,vec14,vec1,AA) 
      dot_prod = vec1(1)*vec15(1) + vec1(2)*vec15(2) + vec1(3)*vec15(3)
  
      if(dot_prod.gt.0.0) then
       lefthand = .FALSE.
      else
       lefthand = .TRUE.
       !write(6,*) 'non-right hand element detected'
      endif  

      return
      end
!--------------------------------------------------------------------
      subroutine fix_left_hand(iel)
! fix fix_left_hand element
! 1 <->19, 3 <-> 21, 7<->25, 9 <->27 
      use SIZE
      integer iel,iver
      real*8 xm2(27),ym2(27),zm2(27)
      character(3) cbc5,cbc6
      real bc55,bc56

! swap vertex
      do iver = 1,27
       xm2(iver) = xm1(iver,1,1,iel)
       ym2(iver) = ym1(iver,1,1,iel)
       zm2(iver) = zm1(iver,1,1,iel)
      enddo

      do iver = 1,9
       xm1(iver,1,1,iel) = xm2(iver+18)
       ym1(iver,1,1,iel) = ym2(iver+18)
       zm1(iver,1,1,iel) = zm2(iver+18)
      enddo
     
      do iver = 19,27
       xm1(iver,1,1,iel) = xm2(iver-18)
       ym1(iver,1,1,iel) = ym2(iver-18)
       zm1(iver,1,1,iel) = zm2(iver-18)
      enddo

! swap face 5 <-> 6

      cbc5 = cbc(5,iel)
      bc55 = bc(5,5,iel)
 
      cbc6 = cbc(6,iel)
      bc56 = bc(5,6,iel)

      cbc(5,iel)   = cbc6
      bc (5,5,iel) = bc56

      cbc(6,iel)   = cbc5
      bc (5,6,iel) = bc55 

      return
      end
!--------------------------------------------------------------------
      FUNCTION VOLUM0(P1,P2,P3,P0)
!
!                           3
!     Given four points in R , (P1,P2,P3,P0), VOLUM0 returns
!     the volume enclosed by the parallelagram defined by the
!     vectors { (P1-P0),(P2-P0),(P3-P0) }.  This routine has
!     the nice feature that if the 3 vectors so defined are
!     not right-handed then the volume returned is negative.

      REAL*8 P1(3),P2(3),P3(3),P0(3)

         U1=P1(1)-P0(1)
         U2=P1(2)-P0(2)
         U3=P1(3)-P0(3)

         V1=P2(1)-P0(1)
         V2=P2(2)-P0(2)
         V3=P2(3)-P0(3)

         W1=P3(1)-P0(1)
         W2=P3(2)-P0(2)
         W3=P3(3)-P0(3)

         CROSS1 = U2*V3-U3*V2
         CROSS2 = U3*V1-U1*V3
         CROSS3 = U1*V2-U2*V1

         VOLUM0  = W1*CROSS1 + W2*CROSS2 + W3*CROSS3
         
      RETURN
      END	  
!-----------------------------------------------------------------------
      FUNCTION CRSS2D(XY1,XY2,XY0)
      REAL XY1(2),XY2(2),XY0(2)

         V1X=XY1(1)-XY0(1)
         V2X=XY2(1)-XY0(1)
         V1Y=XY1(2)-XY0(2)
         V2Y=XY2(2)-XY0(2)
         CRSS2D = V1X*V2Y - V1Y*V2X

      RETURN
      END