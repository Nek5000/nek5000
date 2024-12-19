      program exo2nek
! implement new interface
!
!
      use SIZE

      integer option,oct_input1,oct_input2,shell_input
      integer iexo1,flag, ne_nrh
      logical if_pre
!-----------------------------------------------------------

      etot_est = 0

      write(6,*) 'please input number of fluid exo files:'
      read (5,*) fnexo

      allocate ( fluidexo  (32,fnexo))
      allocate ( f_elem_exo  (fnexo))
	  
      do iexo = 1,fnexo

      write(6,*) 'please input exo file:'
      flag = 1
      call read_input_name(flag)
      call exodus_read_pre(flag)
 
      enddo

      write(6,*) 'please input number of solid exo files for CHT problem (input 0 for no solid mesh):'
      read (5,*) snexo

      if (snexo.gt.0) then
	  
      allocate ( solidexo  (32,snexo))
      allocate ( s_elem_exo  (snexo))

      do iexo = 1,snexo

      write(6,*) 'please input exo file:'
      flag = 2
      call read_input_name(flag)
      call exodus_read_pre(flag)

      enddo

      endif ! if (snexo.gt.0) then

      write(6,*) 'done pre-read exo files'
	  write(6,*) 'now converting to nek mesh'
	  
      num_elem = etot_est
	  
	  ! allocate Nek5000 array
      allocate ( xm1                (3,3,3,etot_est)      )
      allocate ( ym1                (3,3,3,etot_est)      )
      allocate ( zm1                (3,3,3,etot_est)      )

      call rzero(xm1,3*3*3*etot_est)
      call rzero(ym1,3*3*3*etot_est)
      call rzero(zm1,3*3*3*etot_est)

      allocate   (ccurve (4+8*(num_dim-2),etot_est) )
      allocate   (curve  (2*num_dim,12,   etot_est) )
      call rzero (curve,2*num_dim*12*etot_est)
      call blank (ccurve,(4+8*(num_dim-2))*etot_est)

      allocate   (cbc    (2*num_dim,      etot_est) )
      allocate   (bc     (5,2*num_dim,    etot_est) )
      call rzero (bc,5*2*num_dim*etot_est)
      call blank (cbc,3*2*num_dim*etot_est)

      allocate (r_or_l(num_elem))
	  call rzero_int(r_or_l,num_elem)


      eacc = 0
      eacc_old = 0

      do iexo = 1,fnexo
      flag = 1
      call trasnfer_exo_name(flag) ! copy fluidexo to exoname
      call exodus_read_new  ! actual read exodus file

      eacc_old = eacc
      if (num_dim.eq.2) then
         call convert_new
      else if (num_dim.eq.3) then
         if(converting_option.EQ.1) then
         write(6,*) 'Converting pure hex20 mesh'
         call convert_new

         elseif(converting_option.EQ.2) then
         write(6,*) 'Doing linear tet2hex conversion for hybrid (tet4+wedge6+hex8) mesh'
         call split_convert_new
		 
       elseif (converting_option.EQ.3) then ! for quadratic tet2hex
       write(6,*) 'Doing quadratic tet2hex conversion for hybrid (tet10+wedge15+hex20) mesh'

      ! this is to fix some wedge element, which maybe too thin
      quadratic_option = 1
      do i =1,4
      call split_convert1_quadratic
      enddo
      ! this is to linearize tet/wedge/hex to make sure no non-right-hand elements
      quadratic_option = 2
      do i = 1,2
      call split_convert1_quadratic
      enddo

      ! this is the actual splitting step
      quadratic_option = 3
      call split_convert1_quadratic

         endif
      endif

      f_elem_exo(iexo) = eacc - eacc_old

      !call offset_sideset()

      call checkXYZ_min_max()
      write(6,*) 'total element now is ',eacc
      write(6,*) 'fluid exo file ',iexo,' has elements ',f_elem_exo(iexo)
  
      enddo

      etot = eacc
      eftot = etot

      if (snexo.gt.0) then

      do iexo = 1,snexo
      flag = 2
      call trasnfer_exo_name(flag)  ! copy fluidexo to exoname
      call exodus_read_new
	  
      eacc_old = eacc

      if (num_dim.eq.2) then
         call convert_new
      else if (num_dim.eq.3) then
         if(converting_option.EQ.1) then
         write(6,*) 'Converting pure hex20 mesh'
         call convert_new
         elseif(converting_option.EQ.2) then

      write(6,*) 'Doing linear tet2hex conversion for hybrid (tet4+wedge6+hex8) mesh'

         call split_convert_new
		 
       elseif (converting_option.EQ.3) then ! for quadratic tet2hex

      write(6,*) 'Doing quadratic tet2hex conversion for hybrid (tet10+wedge15+hex20) mesh'

      ! this is to fix some wedge element, which maybe too thin
      quadratic_option = 1
      do i =1,4
      call split_convert1_quadratic
      enddo
      ! this is to linearize tet/wedge/hex to make sure no non-right-hand elements
      quadratic_option = 2
      do i = 1,2
      call split_convert1_quadratic
      enddo
      ! this is the actual splitting step
      quadratic_option = 3
      call split_convert1_quadratic

         endif
      endif

      s_elem_exo(iexo) = eacc - eacc_old
      !call offset_sideset()

      call checkXYZ_min_max()
      write(6,*) 'total element now is ',eacc
      write(6,*) 'solid exo file ',iexo,' has elements ',s_elem_exo(iexo)

      enddo

      endif ! if (snexo.gt.0) then

      etot = eacc
      num_elem = etot

!      write(6,*) 'EXPERIMENTAL feature: do you want to try shell feature? (1 for yes)'
!      write(6,*) 'this is only for fluid_only mesh 3D'
!      write(6,*) 'extrude a solid shell region with given thickness and on given sideset' 
!      read (5,*) shell_input
!	  if (shell_input.eq.1)  call shell_wall()

!      write(6,*) 'EXPERIMENTAL feature: do you want to try octuple splitting? (1 for yes)'
!      read (5,*) oct_input1
!	  if (oct_input1.eq.1)  then
!        write(6,*)  'how many rounds of octuple splitting? (2 most)'
!        read (5,*) oct_input2
!        call octuple_split_multirounds(oct_input2)
!        STOP
!        !call octuple_split0()  ! allocate arrays, no exporting on-run
!      endif

      ne_nrh = 0
      call right_hand_check(ne_nrh) ! check non-right-hand element here

      if ((num_dim.eq.3).and.(ne_nrh.gt.0)) then
      call fix_left_hand_elements_3d
      endif

      call gather_bc_info

      call set_periodicity(1)
      if (eftot.ne.num_elem)  call set_periodicity(2)
	  
      write(6,*) 'please give re2 file name:'
      call read_re2_name
      call gen_re2

      end 
!-----------------------------------------------------------------------
      subroutine checkXYZ_min_max
! return element bound for exo file iexo
      use SIZE
      real*8 xx,yy,zz
      integer ie,i

      maxxyz(1) = -1e6
      maxxyz(2) = -1e6
      maxxyz(3) = -1e6

      minxyz(1) = 1e6
      minxyz(2) = 1e6 
      minxyz(3) = 1e6 

      do ie = 1,eacc
        do i = 1,27
       
        xx = xm1(i,1,1,ie)        
        yy = ym1(i,1,1,ie)
        zz = zm1(i,1,1,ie)

        if(xx.gt.maxxyz(1)) maxxyz(1)  = xx
        if(yy.gt.maxxyz(2)) maxxyz(2)  = yy
        if(zz.gt.maxxyz(3)) maxxyz(3)  = zz

        if(xx.lt.minxyz(1)) minxyz(1)  = xx
        if(yy.lt.minxyz(2)) minxyz(2)  = yy
        if(zz.lt.minxyz(3)) minxyz(3)  = zz

        enddo
      enddo

      write(6,*) 'Domain max xyz:',maxxyz(1),maxxyz(2),maxxyz(3) 
      write(6,*) 'Domain min xyz:',minxyz(1),minxyz(2),minxyz(3) 

      return
      end
!------------------------------------------------------------------------
      subroutine read_input_name(flag)

      use SIZE

      character(1)  exonam1(32)
      character(32) fname
      integer flag
      logical e_file_exist,exo_file_exist
	  
      read (5,'(A32)') fname
      len = ltrunc(fname,32)
      
      call blank  (exonam1, 32)
      call chcopy (exonam1,fname,32)
      call chcopy (exonam1(len+1) ,'.exo',4)

      call blank  (exoname, 32)
      call chcopy (exoname,exonam1,len+4)
 
      inquire(FILE=exoname, EXIST=exo_file_exist)

      if ( .not. exo_file_exist) then

	  call blank  (exonam1, 32)
      call chcopy (exonam1,fname,32)
      call chcopy (exonam1(len+1) ,'.e',2)

      call blank  (exoname, 32)
      call chcopy (exoname,exonam1,len+4)
      inquire(FILE=exoname, EXIST=e_file_exist)

      if ( .not. e_file_exist) then
       write(6,*) "ERROR: input exodus file does not exist"	  
      endif
      endif
 
      if (flag.eq.1) then
      call blank  (fluidexo(1,iexo), 32)
      call chcopy (fluidexo(1,iexo),exonam1,len+4)
      else if (flag.eq.2) then
      call blank  (solidexo(1,iexo), 32)
      call chcopy (solidexo(1,iexo),exonam1,len+4)
      endif
 
      return 
      end
!------------------------------------------------------------------------
      subroutine trasnfer_exo_name(flag)
! copy fluidexo(1,iexo) into exoname

      use SIZE

      character(1)  exonam1(32)
      character(32) fname
      integer flag
	  
      if (flag.eq.1) then

      call blank  (exoname, 32)
      call chcopy (exoname,fluidexo(1,iexo),32)
      else if (flag.eq.2) then
      call blank  (exoname, 32)
      call chcopy (exoname,solidexo(1,iexo),32)
      endif

      return 
      end
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine read_re2_name

      use SIZE

      character(1)  re2nam1(80)
      character(32) fname

      read (5,'(A32)') fname
      len = ltrunc(fname,32)
      
      call blank  (re2nam1, 80)
      call chcopy (re2nam1,fname,80)
      call chcopy (re2nam1(len+1) ,'.re2',4)

      call blank  (re2name, 80)
      call chcopy (re2name,re2nam1,len+4)
 
      return 
      end
!-----------------------------------------------------------------------
      subroutine exodus_read_pre(flag)
!
! pre read exo file information
!
      use SIZE
      include 'exodusII.inc'

      integer exoid, cpu_ws, io_ws, flag

      character(MXSTLN) typ, qa_record(4,10)
      character(MXLNLN) titl
      character(1)      cdum
	  
      character(3)  typ3
      character(4)  typ4
      character(5)  typ5

!      integer,allocatable,dimension(:) :: idblk
!      integer,allocatable,dimension(:) :: num_nodes_per_elem
!      integer,allocatable,dimension(:) :: num_attr    !not used

!
! open EXODUS II file
!
      cpu_ws = 8 ! use real*8 to communicate with exodus
      io_ws  = 0
      exoid  = exopen (exoname, EXREAD, cpu_ws, io_ws, vers, ierr)
      if (ierr.lt.0) then
        write(6,'(2a)') "ERROR: cannot open file ", exoname 
        STOP
      endif
      write(6,*)
      write(6,'(a32,a,f4.2)') & 
            exoname," is an EXODUSII file; version ",vers
      write(6,'(a,i2)') "I/O word size", io_ws
!
! read database parameters
!
      call exgini (exoid, titl, num_dim, num_nodes, num_elem, &
                   num_elem_blk, num_node_sets, num_side_sets, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read exodusII parameters (exgini)"
        STOP
      endif
      write (6, '(/"database parameters:"/ /        &
                 "title         = ", a81 / /        &
                 "num_dim       = ", i8 /           &
                 "num_nodes     = ", i8 /           &
                 "num_elem      = ", i8 /           &
                 "num_elem_blk  = ", i8 /           &
                 "num_side_sets = ", i8)')          &
                 titl,num_dim, num_nodes, num_elem, &
                  num_elem_blk, num_side_sets
      write (6,*)
!
! allocate some arrays    
!
      ! EXODUS:
      allocate ( idblk              (num_elem_blk)        )
      allocate ( num_nodes_per_elem (num_elem_blk)        )
      allocate ( num_attr           (num_elem_blk)        )
      allocate ( num_elem_in_block  (num_elem_blk)        )

!
! read element block parameters
!
      call exgebi (exoid, idblk, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read block ids (exgebi)"
        STOP
      endif


      
      converting_option = 1


      do i = 1, num_elem_blk
        call exgelb (exoid, idblk(i), typ, num_elem_in_block(i), &
                     num_nodes_per_elem(i), num_attr(i), ierr)
        if (ierr.lt.0) then
          write(6,'(a,i3,a)') &
          "ERROR: cannot read parameters for block ",i," (exgelb)"
          STOP
        endif
        write (6, '("element block id   = ", i8,/       &
                   "element type       = ", 3x,a8,/     &
                   "num_elem_in_block  = ", i8,/        &
                   "num_nodes_per_elem = ", i8)')       &
                   idblk(i), typ, num_elem_in_block(i), &
                   num_nodes_per_elem(i)
        write(6,*)
      nvert=num_nodes_per_elem(i)

!----------------------------------------------------------------------
! check element type
!
      if (num_dim.eq.2) then
        ! if 2d mesh, only quad8 elements are allowed
        call chcopy(typ4,typ,4)
        if ((typ4.eq.'QUAD').and.(nvert.eq.8)) then
        ! valid mesh
		  write(6,*) "QUAD8 is the only valid element for 2D mesh"
           etot_est = etot_est + num_elem_in_block(i)
        else
          write(6,*) "ERROR: Only QUAD8 elements are allowed in a 2D mesh!"
          STOP
        endif
      endif

      if (num_dim.eq.3) then
	  
        call chcopy(typ3,typ,3) 
        call chcopy(typ5,typ,5) 

        if ((typ3.eq.'HEX').and.(nvert.eq.20)) then 
           write(6,*) "HEX20 is valid element in a 3D mesh."
           ! write(6,*) "in this scenario, element type should all be HEX20"
           converting_option = converting_option*1
        else if ((typ3.eq.'HEX').and.(nvert.eq.8)) then 
           write(6,*) "HEX8 is valid element in a 3D mesh."
           write(6,*) "assume linear hybrid mesh (tetra-hex-wedge)"
           write(6,*) "one HEX8 divide into 8 Nek hex elements"
           converting_option = 2
           !etot_est = etot_est + num_elem_in_block(i)*8
        else if ((typ5.eq.'TETRA').and.(nvert.eq.4))then 
           write(6,*) "TETRA4 is valid element in a 3D mesh."
           write(6,*) "assume linear hybrid mesh (tetra-hex-wedge)"
           write(6,*) "one TETRA4 divide into 4 Nek hex elements"
           converting_option = 2
           !etot_est = etot_est + num_elem_in_block(i)*4
        else if ((typ5.eq.'WEDGE').and.(nvert.eq.6)) then 
           write(6,*) "WEDGE6 is valid element in a 3D mesh."
           write(6,*) "assume linear hybrid mesh (tetra-hex-wedge)"
           write(6,*) "one WEDGE6 divide into 6 Nek hex elements"
           converting_option = 2
           !etot_est = etot_est + num_elem_in_block(i)*6
        else if ((typ5.eq.'TETRA').and.(nvert.eq.10))then 
           write(6,*) "TETRA10 is valid element in a 3D mesh."
           write(6,*) "assume quadratic hybrid mesh (tetra-wedge)"
           write(6,*) "one TETRA10 divide into 4 Nek hex elements"
           converting_option = 3
           !etot_est = etot_est + num_elem_in_block(i)*4
        else if ((typ5.eq.'WEDGE').and.(nvert.eq.15)) then 
           write(6,*) "WEDGE15 is valid element in a 3D mesh."
           write(6,*) "assume quadratic hybrid mesh (tetra-wedge)"
           write(6,*) "one WEDGE15 divide into 3 Nek hex elements"
           converting_option = 3
           !etot_est = etot_est + num_elem_in_block(i)*3
           !etot_est = etot_est + num_elem_in_block(i)*6
        else
          write(6,*) "ERROR: invalid element in a 3D mesh!"
          STOP
        endif
		

      endif
      enddo


        if ((flag.eq.1).and.(iexo.eq.1)) then
         converting_option_old = converting_option
        else
          if (converting_option.ne.converting_option_old) then
            write(6,*) "ERROR: inconsistent mesh type of exo file"
            STOP
          endif
        endif
		

!     etot_est = 0

      do i = 1, num_elem_blk
        call exgelb (exoid, idblk(i), typ, num_elem_in_block(i), &
                     num_nodes_per_elem(i), num_attr(i), ierr)
        if (ierr.lt.0) then
          write(6,'(a,i3,a)') &
          "ERROR: cannot read parameters for block ",i," (exgelb)"
          STOP
        endif
        write (6, '("element block id   = ", i8,/       &
                   "element type       = ", 3x,a8,/     &
                   "num_elem_in_block  = ", i8,/        &
                   "num_nodes_per_elem = ", i8)')       &
                   idblk(i), typ, num_elem_in_block(i), &
                   num_nodes_per_elem(i)
        write(6,*)
      nvert=num_nodes_per_elem(i)

!----------------------------------------------------------------------
! check element type
!
      if (num_dim.eq.2) then
        ! if 2d mesh, only quad8 elements are allowed
        call chcopy(typ4,typ,4)
        if ((typ4.eq.'QUAD').and.(nvert.eq.8)) then
        ! valid mesh
		  write(6,*) "QUAD8 is the only valid element for 2D mesh"
           etot_est = etot_est + num_elem_in_block(i)
        else
          write(6,*) "ERROR: Only QUAD8 elements are allowed in a 2D mesh!"
          STOP
        endif
      endif

      if (num_dim.eq.3) then
	  
        call chcopy(typ3,typ,3) 
        call chcopy(typ5,typ,5) 

        if (((typ3.eq.'HEX').and.(nvert.eq.20)).and.(converting_option.eq.1)) then 
           etot_est = etot_est + num_elem_in_block(i)
        else if (((typ3.eq.'HEX').and.(nvert.eq.20)).and.(converting_option.eq.3)) then 
           etot_est = etot_est + num_elem_in_block(i)*8
        else if ((typ3.eq.'HEX').and.(nvert.eq.8)) then 
           etot_est = etot_est + num_elem_in_block(i)*8
        else if ((typ5.eq.'TETRA').and.(nvert.eq.4))then 
           etot_est = etot_est + num_elem_in_block(i)*4
        else if ((typ5.eq.'WEDGE').and.(nvert.eq.6)) then 
           etot_est = etot_est + num_elem_in_block(i)*6
        else if ((typ5.eq.'TETRA').and.(nvert.eq.10))then 
           etot_est = etot_est + num_elem_in_block(i)*4
        else if ((typ5.eq.'WEDGE').and.(nvert.eq.15)) then 
           !etot_est = etot_est + num_elem_in_block(i)*3
           etot_est = etot_est + num_elem_in_block(i)*6
        else
          write(6,*) "ERROR: invalid element in a 3D mesh!"
          STOP
        endif

      endif
      enddo


          write(6,*) "converting_option:: ", converting_option
          write(6,*) "etot_est:: ", etot_est





      deallocate ( idblk )
      deallocate ( num_nodes_per_elem )
      deallocate ( num_attr           )
      deallocate ( num_elem_in_block  )
	  
      return
      end
! -------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine exodus_read_new
!
!  Subroutine to read an exodusII binary file containing a mesh.
!  It uses exodus fortran binding subroutines, which depend on
!  the netcdf library for low level data access.
!
      use SIZE
      include 'exodusII.inc'

      integer exoid, cpu_ws, io_ws

      character(MXSTLN) typ, qa_record(4,10)
      character(MXLNLN) titl
      character(1)      cdum

!      integer,allocatable,dimension(:) :: idblk
!      integer,allocatable,dimension(:) :: num_nodes_per_elem
!      integer,allocatable,dimension(:) :: num_attr    !not used

!
! open EXODUS II file
!
      cpu_ws = 8 ! use real*8 to communicate with exodus
      io_ws  = 0
      exoid  = exopen (exoname, EXREAD, cpu_ws, io_ws, vers, ierr)
!      if (ierr.lt.0) then
!        write(6,'(2a)') "ERROR: cannot open file ", exoname 
!        STOP
!      endif
!      write(6,*)
!      write(6,'(a32,a,f4.2)') & 
!            exoname," is an EXODUSII file; version ",vers
!      write(6,'(a,i2)') "I/O word size", io_ws
!
! read database parameters
!
      call exgini (exoid, titl, num_dim, num_nodes, num_elem, &
                   num_elem_blk, num_node_sets, num_side_sets, ierr)
!     if (ierr.lt.0) then
!       write(6,'(a)') "ERROR: cannot read exodusII parameters (exgini)"
!       STOP
!     endif
!     write (6, '(/"database parameters:"/ /        &
!                "title         = ", a81 / /        &
!                "num_dim       = ", i8 /           &
!                "num_nodes     = ", i8 /           &
!                "num_elem      = ", i8 /           &
!                "num_elem_blk  = ", i8 /           &
!                "num_side_sets = ", i8)')          &
!                titl,num_dim, num_nodes, num_elem, &
!                 num_elem_blk, num_side_sets
!     write (6,*)
!
! allocate some arrays    
!
      ! EXODUS:
      allocate ( idblk              (num_elem_blk)        )
      allocate ( num_nodes_per_elem (num_elem_blk)        )
      allocate ( num_attr           (num_elem_blk)        )
      allocate ( num_elem_in_block  (num_elem_blk)        )
      allocate ( num_sides_in_set   (num_side_sets)       )
      allocate ( idss               (num_side_sets)       )
      allocate ( connect            (3**num_dim*num_elem) )
      allocate ( x_exo              (3**num_dim*num_elem) )
      allocate ( y_exo              (3**num_dim*num_elem) )
      allocate ( z_exo              (3**num_dim*num_elem) )

      ! Nek5000:
!      allocate ( xm1                (3,3,3,num_elem*8)      )
!      allocate ( ym1                (3,3,3,num_elem*8)      )
!      allocate ( zm1                (3,3,3,num_elem*8)      )

!
! read element block parameters
!
      call exgebi (exoid, idblk, ierr)
!      if (ierr.lt.0) then
!        write(6,'(a)') "ERROR: cannot read block ids (exgebi)"
!        STOP
!      endif

      do i = 1, num_elem_blk
        call exgelb (exoid, idblk(i), typ, num_elem_in_block(i), &
                     num_nodes_per_elem(i), num_attr(i), ierr)
!        if (ierr.lt.0) then
!          write(6,'(a,i3,a)') &
!          "ERROR: cannot read parameters for block ",i," (exgelb)"
!          STOP
!        endif
!        write (6, '("element block id   = ", i8,/       &
!                   "element type       = ", 3x,a8,/     &
!                   "num_elem_in_block  = ", i8,/        &
!                   "num_nodes_per_elem = ", i8)')       &
!                   idblk(i), typ, num_elem_in_block(i), &
!                   num_nodes_per_elem(i)
!        write(6,*)

!        if (i.eq.1) then
!          nvert=num_nodes_per_elem(i)
!          if (num_dim.eq.2) then
!            if (nvert.ne.8) then
!              if (nvert.eq.9) then
!                write(6,*)
!                write(6,'(a)') &
!                "WARNING: QUAD9 elements are not officially supported"
!                write(6,'(a)') &
!                "as there is no exodus standard for this element type."
!                write(6,*)
!              else
!                write(6,*)
!                write(6,'(a)') &
!                "ERROR: Only QUAD8 elements are allowed in a 2D mesh!"
!                STOP
!              endif
!            endif      
!          elseif (num_dim.eq.3) then
!            if (nvert.ne.20) then
!              if (nvert.eq.27) then
!                write(6,*)
!                write(6,'(a)') &
!                "WARNING: HEX27 elements are not officially supported"
!                write(6,'(a)') &
!                "as there is no exodus standard for this element type."
!                write(6,*)
!              else
!                write(6,*)
!                write(6,'(a)') &
!                "ERROR: Only HEX20 elements are allowed in a 3D mesh!"
!                STOP
!              endif
!            endif      
!          else
!          write(6,'(a,i3)') &
!           "ERROR: Unknown number of dimensions! num_dim= ",num_dim
!          STOP
!          endif
!        endif
!
!        if (num_nodes_per_elem(i).ne.nvert) then
!          write(6,*)
!          write(6,'(a)') &
!           "ERROR: All blocks should contain elements of the same type!"
!          write(6,'(a,i3)') "num_nodes_per_elem= ",num_nodes_per_elem(i)
!          write(6,'(a,i3)') "num_nodes_per_elem of block 1 ",nvert
!          STOP
!        endif
      enddo
!
! read nodal coordinates values from database
!
      call exgcor (exoid, x_exo, y_exo, z_exo, ierr)
!      if (ierr.lt.0) then
!        write(6,'(a)') "ERROR: cannot read nodal coordinates (exgcor)"
!       STOP
!     endif
!
! read element connectivity
!
      iend = 0
      do 60 i = 1, num_elem_blk
        istart = iend + 1
        call exgelc (exoid, idblk(i), connect(istart), ierr)
        iend = iend+num_nodes_per_elem(i)*num_elem_in_block(i)
!        if (ierr.lt.0) then
!          write(6,'(a)') "ERROR: cannot read elm. connectivity (exgelc)"
!          STOP
!        endif
60    continue
!
! read individual side sets
!
      num_sides_tot = 0
      if (num_side_sets .gt. 0) then
        call exgssi (exoid, idss, ierr)
        if (ierr.lt.0) then
          write(6,'(a)') "ERROR: cannot read SideSet ids (exgssi)"
          STOP
        endif
        maxnss=0
        do i = 1, num_side_sets
          call exgsp (exoid,idss(i),num_sides_in_set(i),idum,ierr)
          if (ierr.lt.0) then
            write(6,'(a,i3,a)') &
            "ERROR: cannot read parameters for SideSet No.",i," (exgsp)"
            STOP
          endif
!          write (6, '("side set ", i2, " num_sides = ", i8)') &
!                 idss(i), num_sides_in_set(i)
          num_sides_tot = num_sides_tot + num_sides_in_set(i)
          maxnss        = max(maxnss,num_sides_in_set(i))
        enddo

        ! allocate sideset arrays
        allocate (elem_list(maxnss,num_side_sets) )
        allocate (side_list(maxnss,num_side_sets) )

        do i = 1, num_side_sets
          call exgss (exoid,idss(i),elem_list(1,i),side_list(1,i),ierr)
          if (ierr.lt.0) then
            write(6,'(a,i3,a)') &
            "ERROR: cannot read parameters for SideSet No.",i," (exgss)"
            STOP
          endif
        enddo
        write (6,*)
      else
        write(6,'(a)') "WARNING: No SideSets in exodus file!"
      endif
!
! read QA records
!
      call exinq (exoid, EXQA, num_qa_rec, fdum, cdum, ierr)
!      if (ierr.lt.0) then
!        write(6,'(a)') "ERROR: cannot read QA records (exinq QA) "
!        STOP
!      elseif (ierr.gt.0) then
!        write(6,'(a)') "INFO: file does not contain any QA records"
!      else
!        call exgqa (exoid, qa_record, ierr)
!        if (ierr.lt.0) then
!          write(6,'(a,i3)') "WARNING: cannot read QA records (exgqa)"
!        else
!          write (6, '("QA records = ")')
!          if (num_qa_rec.gt.10) then
!            write(6,'(2a)') &
!              'WARNING: Cannot handle more than 10 QA records', &
!              'Printing only the first 10...'
!          else
!            do i = 1, num_qa_rec
!              do j = 1, 4
!                write (6,'(a)') qa_record(j,i)
!              enddo
!            enddo
!          endif
!        endif
!      endif

      return
      end
! -------------------------------------------------------------------
!-----------------------------------------------------------------
      subroutine convert_new
!
!  Subroutine to convert an already read exodusII mesh to a nek
!  mesh. The idea is to fill each element's node coordinates of
!  size lx1**3 (3D) or lx1**2 (2D) (lx1=3) with the hex27/quad9
!  coordinates.
!
      use SIZE
      include 'exodusII.inc'
      integer ssID
      logical ifinterface

! node and face conversion (it works at least for cubit):
      integer exo_to_nek_vert3D(27)                   ! hex27 to nek numbering
      data    exo_to_nek_vert3D                     &
            / 19,  1,  7, 25, 21,  3,  9, 27, 10    &     
            ,  4, 16, 22, 20,  2,  8, 26, 12,  6    &
            , 18, 24, 14, 13, 15, 23,  5, 11, 17 /

      integer exo_to_nek_vert2D(9)                     ! quad9 to nek numbering
      data    exo_to_nek_vert2D   / 1, 3, 9, 7, 2, 6, 8, 4, 5  / 

      integer exo_to_nek_left(9)
      data    exo_to_nek_left /3,1,7,9,2,4,8,6,5/  ! LEFT HAND SIDE ELEMENT
 
      integer exo_to_nek_face3D(6)
      data    exo_to_nek_face3D  / 1, 5, 3, 6, 4, 2 /  ! symmetric face numbering

      integer exo_to_nek_face2D(4)
      data    exo_to_nek_face2D  / 1, 2, 3, 4 /        ! symmetric face numbering
      integer exo_to_nek_face2D_left(4)
      data    exo_to_nek_face2D_left  / 1,4,3,2 /        ! symmetric face numbering

      eacc_old = eacc
      nvert = num_nodes_per_elem(1)
      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '
      do iel=1,num_elem
        eacc = eacc + 1

        !if (eacc.gt.etot_est) write(6,*) 'ERROR: please increase estimate final total hex element number' 

        if (num_dim.eq.2) then
                      
          call r_or_l_detect_for_quad(iel,r_or_l(eacc))
	  
		  do ivert =1,nvert
            if (r_or_l(eacc).eq.0) then
            jvert = exo_to_nek_vert2D(ivert)
			xm1(jvert,1,1,eacc)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,eacc)=y_exo(connect(nvert*(iel-1)+ivert))
            else
            jvert = exo_to_nek_left(ivert)
			xm1(jvert,1,1,eacc)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,eacc)=y_exo(connect(nvert*(iel-1)+ivert))
            endif
          enddo
			
 
        else 

          do ivert =1,nvert
            jvert = exo_to_nek_vert3D(ivert)
            xm1(jvert,1,1,eacc)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,eacc)=y_exo(connect(nvert*(iel-1)+ivert))
            zm1(jvert,1,1,eacc)=z_exo(connect(nvert*(iel-1)+ivert))

            xm1(jvert,1,1,eacc)=xm1(jvert,1,1,eacc) + shiftvector(1)
            ym1(jvert,1,1,eacc)=ym1(jvert,1,1,eacc) + shiftvector(2)
            zm1(jvert,1,1,eacc)=zm1(jvert,1,1,eacc) + shiftvector(3)

          enddo

       endif

      enddo

      write(6,'(A)') 'done :: Converting elements '
!
! allocate and zero-out curve sides arrays
!
!      allocate   (ccurve (4+8*(num_dim-2),num_elem) )
!      allocate   (curve  (2*num_dim,12,   num_elem) )
!      call rzero (curve,2*num_dim*12*num_elem)
!      call blank (ccurve,(4+8*(num_dim-2))*num_elem)
!
! allocate and zero-out bc arrays only if sidesets are specified
!
!      if (num_side_sets.eq.0) return   

!      allocate   (cbc    (2*num_dim,      num_elem) )
!      allocate   (bc     (5,2*num_dim,    num_elem) ) 
!      call rzero (bc,5*2*num_dim*num_elem)
!      call blank (cbc,3*2*num_dim*num_elem)

!
! set bc's
!
      write(6,'(a)') ''
      write(6,'(a)') 'Converting SideSets ...'

      do iss=1,num_side_sets
        write(6,'(a)') ''
        write(6,'(a,i0,a)') 'Sideset ',idss(iss), ' ...'

        do i=1,num_sides_in_set(iss) 
          iel = elem_list(i,iss)
          ifc = side_list(i,iss)
          if (num_dim.eq.2) then 
            if (r_or_l(iel+eacc_old).eq.0) then
            jfc = exo_to_nek_face2D(ifc)
            else
            jfc = exo_to_nek_face2D_left(ifc)
            endif
          else
            jfc = exo_to_nek_face3D(ifc)
          endif
           cbc(jfc,iel+eacc_old)   = 'EXO' ! dummy exodus bc
           bc (5,jfc,iel+eacc_old) = idss(iss)
           bc (4,jfc,iel+eacc_old) = iexo
        enddo

        write(6,'(a,i0)') 'done :: Sideset ',idss(iss)
      enddo

      deallocate (x_exo,y_exo,z_exo,connect)
      deallocate (elem_list,side_list)
      deallocate ( idblk )
      deallocate ( num_nodes_per_elem )
      deallocate ( num_attr           )
      deallocate ( num_elem_in_block  )
      deallocate ( num_sides_in_set   )
      deallocate ( idss               )

      write(6,'(a)') ''
      write(6,'(a)') 'done :: Converting SideSets '

      return
      end
!--------------------------------------------------------------------
      subroutine cross_prod(nv,p1,p2,p3)
      real*8 nv(3),p1(3),p2(3),p3(3)
      real*8 v12(3),v13(3),mag
	  
      do i = 1,3
      v12(i) = p2(i) - p1(i)
      v13(i) = p3(i) - p1(i)
      enddo
	  
      nv(1) = v12(2)*v13(3) - v12(3)*v13(2)
      nv(2) = v12(3)*v13(1) - v12(1)*v13(3)
      nv(3) = v12(1)*v13(2) - v12(2)*v13(1)
	  
      mag = sqrt(nv(1)**2.0+nv(2)**2.0+nv(3)**2.0)

      do i = 1,3
      nv(i) = nv(i)/mag
      enddo

      return
      end
!--------------------------------------------------------------------
      subroutine distance(a,b,d)
      real*8 a(3),b(3),d
 
      d = 0.0

      do i=1,3
       d = d + (a(i)-b(i))**2.0
      enddo

      d = d**0.5

      return
      end
!------------------------------------
      subroutine assignvec(a,b)
      real*8 a(3),b(3)
      a(1) = b(1)
      a(2) = b(2)
      a(3) = b(3)
      return
      end
! -------------------------------------------------------------------
      subroutine average2vec(avg,a,b)
      real*8 a(3),b(3),avg(3)      
      avg(1) = (a(1)+b(1))*0.5
      avg(2) = (a(2)+b(2))*0.5
      avg(3) = (a(3)+b(3))*0.5
      return
      end
!--------------------------------------------------------------------
      subroutine average3vec(avg,a,b,c)
      real*8 a(3),b(3),c(3),avg(3)      
      avg(1) = (a(1)+b(1)+c(1))*0.33333333
      avg(2) = (a(2)+b(2)+c(2))*0.33333333
      avg(3) = (a(3)+b(3)+c(3))*0.33333333

      return
      end
! -------------------------------------------------------------------
      subroutine average4vec(avg,a,b,c,d)
      real*8  a(3),b(3),c(3),d(3),avg(3)       
      avg(1) = (a(1)+b(1)+c(1)+d(1))*0.25
      avg(2) = (a(2)+b(2)+c(2)+d(2))*0.25
      avg(3) = (a(3)+b(3)+c(3)+d(3))*0.25
      return
      end
! -------------------------------------------------------------------
!--------------------------------------------------------------------
      subroutine setbc_3d
      use SIZE
	 
      integer hex_face_node(4,6)
      data hex_face_node /1,3,21,19,3,9,27,21,7,9,27,25,1,7,25,19,1,7,9,3,19,21,27,25/
      
      character*3 ubc
      integer tags(2),ibc,nbc,io
      integer ip,np,ipe,ipe2,nipe(2)
      integer ptags(2)
      integer fnode(4)
      real pvec(3)
      real fpxyz(3,2)
      real AB_v(3),AD_v(3),farea,product_v(3)
      real dist,distMax,ptol

! boundary condition summary
      write(6,*) '******************************************************'
      write(6,*) 'Boundary info summary'
      write(6,*) 'sideSet ID'
      do ibc= 1,bcNumber
      write(6,*) bcID(ibc)
      enddo
      write(6,*) '******************************************************'
 
 
      write(6,*) 'Enter number of periodic boundary surface pairs:'
      read (5,*) nbc
	  
      if(nbc.ne.0) then
      write(6,*) 'Enter search tolerance (1 as default):'
      read (5,*) ptol
      endif
	  
      if(nbc.le.0) return
	  
      allocate ( parray (2,2,num_elem))

      do ibc = 1,nbc 
        write(6,*) 'input surface 1 and  surface 2  sideSet ID'
        read (5,*) ptags(1),ptags(2)
        write(6,*) 'input translation vector (surface 1 -> surface 2)'
        read (5,*) pvec(1),pvec(2),pvec(3)

          ipe = 0
          do ihex = 1, num_elem
            do iface = 1,6
               if(bc(5,iface,ihex).eq.ptags(1)) then   
                ipe = ipe + 1
                parray(1,1,ipe) = ihex
                parray(2,1,ipe) = iface
               endif
            enddo
          enddo
          nipe(1) = ipe
	  
          ipe = 0
          do ihex = 1, num_elem
            do iface = 1,6
               if(bc(5,iface,ihex).eq.ptags(2)) then
                ipe = ipe + 1
                parray(1,2,ipe) = ihex
                parray(2,2,ipe) = iface
               endif
            enddo
          enddo
          nipe(2) = ipe

          if(nipe(1).ne.nipe(2))  then
            write(6,*)'mapping sideset ',ptags(1),'with',nipe(1),'faces'
            write(6,*)'to sideset ',ptags(2),'with',nipe(2),'faces'
            write(6,*) 'EORROR, face numbers are not matching'
          endif

! 1st loop, loop faces on surface 1
      do ipe = 1,nipe(1)
         ihex = parray(1,1,ipe)
         iface = parray(2,1,ipe)
! get face center xyz
         call rzero(fpxyz(1,1),3)

         do ifnode = 1,4
             fnode(ifnode)=hex_face_node(ifnode,iface)
             fpxyz(1,1) = fpxyz(1,1)+xm1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(2,1) = fpxyz(2,1)+ym1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(3,1) = fpxyz(3,1)+zm1(fnode(ifnode),1,1,ihex)*0.25
         enddo

! 2nd loop over surface 2
         distMax = ptol
         do ipe2 = 1,nipe(2)
                ihex2 = parray(1,2,ipe2)
                iface2 = parray(2,2,ipe2)
! get face center xyz
               call rzero(fpxyz(1,2),3)
               do ifnode = 1,4
               fnode(ifnode)=hex_face_node(ifnode,iface2)
               fpxyz(1,2) = fpxyz(1,2)+xm1(fnode(ifnode),1,1,ihex2)*0.25
               fpxyz(2,2) = fpxyz(2,2)+ym1(fnode(ifnode),1,1,ihex2)*0.25
               fpxyz(3,2) = fpxyz(3,2)+zm1(fnode(ifnode),1,1,ihex2)*0.25
              enddo
 
       dist = sqrt((fpxyz(1,2)-fpxyz(1,1)-pvec(1))**2 &
      +(fpxyz(2,2)-fpxyz(2,1)-pvec(2))**2 &
      +(fpxyz(3,2)-fpxyz(3,1)-pvec(3))**2)

               if (dist.lt.distMax) then 
                  distMax = dist
                  bc(1,iface,ihex) = ihex2*1.0
                  bc(2,iface,ihex) = iface2*1.0
               endif
         enddo
      enddo

! change. only assign periodic face at the end of loop.
! this will assign the closest face for periodicity.

     do ipe = 1,nipe(1)
         ihex = parray(1,1,ipe)
         iface = parray(2,1,ipe)
         ihex2 = int(bc(1,iface,ihex))
         iface2 = int(bc(2,iface,ihex)) 
         bc(1,iface2,ihex2) = ihex*1.0
         bc(2,iface2,ihex2) = iface*1.0
         cbc(iface,ihex) = 'P  '
         cbc(iface2,ihex2) = 'P  '
     enddo

        nperror = 0
   
          write(6,*)'doing periodic check for surface',ptags(1)

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
             if (cbc(iface,ihex).ne.'P  ') then
                  nperror = nperror +1 
             endif
          enddo
          if (nperror.gt.0) write(6,*) 'ERROR,',nperror, ' faces did not for periodicity'

          nperror = 0

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
             ihex2 = int(bc(1,iface,ihex))
             iface2 = int(bc(2,iface,ihex))
             ihex3 = int(bc(1,iface2,ihex2))
             iface3 = int(bc(2,iface2,ihex2))
             if ((ihex.ne.ihex3).or.(iface.ne.iface3)) then
               nperror = nperror + 1
             endif
          enddo

          if (nperror.gt.0) then
          write(6,*) 'ERROR,',nperror,'faces are wrong out of total ',nipe(1),' faces'
          endif
  
          write(6,*)'doing periodic check for surface',ptags(2)

          nperror = 0

          do ipe = 1,nipe(2)
             ihex = parray(1,2,ipe)
             iface = parray(2,2,ipe)
             if (cbc(iface,ihex).ne.'P  ') then
                  nperror = nperror +1 
             endif
          enddo
          if (nperror.gt.0) write(6,*) 'ERROR,',nperror, ' faces did not map for periodicity'

          nperror = 0
          do ipe = 1,nipe(2)
             ihex = parray(1,2,ipe)
             iface = parray(2,2,ipe)
             ihex2 = int(bc(1,iface,ihex))
             iface2 = int(bc(2,iface,ihex))
             ihex3 = int(bc(1,iface2,ihex2))
             iface3 = int(bc(2,iface2,ihex2))
             if ((ihex.ne.ihex3).or.(iface.ne.iface3)) then
                nperror = nperror + 1
             endif
          enddo		  

          if (nperror.gt.0) then
          write(6,*) 'ERROR,',nperror,'faces are wrong out of total ',nipe(2),' faces'
          endif
 
      enddo

      deallocate(parray)
	  
      write(6,*) '******************************************************'
      write(6,*) 'Please set boundary conditions to all non-periodic boundaries'
      write(6,*) 'in .usr file usrdat2() subroutine'
      write(6,*) '******************************************************'

      return
      end
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine exodus_read
!
!  Subroutine to read an exodusII binary file containing a mesh.
!  It uses exodus fortran binding subroutines, which depend on
!  the netcdf library for low level data access.
!
      use SIZE
      include 'exodusII.inc'

      integer exoid, cpu_ws, io_ws

      character(MXSTLN) typ, qa_record(4,10)
      character(MXLNLN) titl
      character(1)      cdum

!      integer,allocatable,dimension(:) :: idblk
!      integer,allocatable,dimension(:) :: num_nodes_per_elem
!      integer,allocatable,dimension(:) :: num_attr    !not used

!
! open EXODUS II file
!
      cpu_ws = 8 ! use real*8 to communicate with exodus
      io_ws  = 0
      exoid  = exopen (exoname, EXREAD, cpu_ws, io_ws, vers, ierr)
      if (ierr.lt.0) then
        write(6,'(2a)') "ERROR: cannot open file ", exoname 
        STOP
      endif
      write(6,*)
      write(6,'(a32,a,f4.2)') & 
            exoname," is an EXODUSII file; version ",vers
      write(6,'(a,i2)') "I/O word size", io_ws
!
! read database parameters
!
      call exgini (exoid, titl, num_dim, num_nodes, num_elem, &
                   num_elem_blk, num_node_sets, num_side_sets, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read exodusII parameters (exgini)"
        STOP
      endif
      write (6, '(/"database parameters:"/ /        &
                 "title         = ", a81 / /        &
                 "num_dim       = ", i8 /           &
                 "num_nodes     = ", i8 /           &
                 "num_elem      = ", i8 /           &
                 "num_elem_blk  = ", i8 /           &
                 "num_side_sets = ", i8)')          &
                 titl,num_dim, num_nodes, num_elem, &
                  num_elem_blk, num_side_sets
      write (6,*)
!
! allocate some arrays    
!
      ! EXODUS:
      allocate ( idblk              (num_elem_blk)        )
      allocate ( num_nodes_per_elem (num_elem_blk)        )
      allocate ( num_attr           (num_elem_blk)        )
      allocate ( num_elem_in_block  (num_elem_blk)        )
      allocate ( num_sides_in_set   (num_side_sets)       )
      allocate ( idss               (num_side_sets)       )
      allocate ( connect            (3**num_dim*num_elem) )
      allocate ( x_exo              (3**num_dim*num_elem) )
      allocate ( y_exo              (3**num_dim*num_elem) )
      allocate ( z_exo              (3**num_dim*num_elem) )
      ! Nek5000:
      allocate ( xm1                (3,3,3,num_elem)      )
      allocate ( ym1                (3,3,3,num_elem)      )
      allocate ( zm1                (3,3,3,num_elem)      )
!
! read element block parameters
!
      call exgebi (exoid, idblk, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read block ids (exgebi)"
        STOP
      endif

      do i = 1, num_elem_blk
        call exgelb (exoid, idblk(i), typ, num_elem_in_block(i), &
                     num_nodes_per_elem(i), num_attr(i), ierr)
        if (ierr.lt.0) then
          write(6,'(a,i3,a)') &
          "ERROR: cannot read parameters for block ",i," (exgelb)"
          STOP
        endif
        write (6, '("element block id   = ", i8,/       &
                   "element type       = ", 3x,a8,/     &
                   "num_elem_in_block  = ", i8,/        &
                   "num_nodes_per_elem = ", i8)')       &
                   idblk(i), typ, num_elem_in_block(i), &
                   num_nodes_per_elem(i)
        write(6,*)

        if (i.eq.1) then
          nvert=num_nodes_per_elem(i)
          if (num_dim.eq.2) then
            if (nvert.ne.8) then
              if (nvert.eq.9) then
                write(6,*)
                write(6,'(a)') &
                "WARNING: QUAD9 elements are not officially supported"
                write(6,'(a)') &
                "as there is no exodus standard for this element type."
                write(6,*)
              else
                write(6,*)
                write(6,'(a)') &
                "ERROR: Only QUAD8 elements are allowed in a 2D mesh!"
                STOP
              endif
            endif      
          elseif (num_dim.eq.3) then
            if (nvert.ne.20) then
              if (nvert.eq.27) then
                write(6,*)
                write(6,'(a)') &
                "WARNING: HEX27 elements are not officially supported"
                write(6,'(a)') &
                "as there is no exodus standard for this element type."
                write(6,*)
              else
                write(6,*)
                write(6,'(a)') &
                "ERROR: Only HEX20 elements are allowed in a 3D mesh!"
                STOP
              endif
            endif      
          else
          write(6,'(a,i3)') &
           "ERROR: Unknown number of dimensions! num_dim= ",num_dim
          STOP
          endif
        endif

        if (num_nodes_per_elem(i).ne.nvert) then
          write(6,*)
          write(6,'(a)') &
           "ERROR: All blocks should contain elements of the same type!"
          write(6,'(a,i3)') "num_nodes_per_elem= ",num_nodes_per_elem(i)
          write(6,'(a,i3)') "num_nodes_per_elem of block 1 ",nvert
          STOP
        endif
      enddo
!
! read nodal coordinates values from database
!
      call exgcor (exoid, x_exo, y_exo, z_exo, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read nodal coordinates (exgcor)"
        STOP
      endif
!
! read element connectivity
!
      iend = 0
      do 60 i = 1, num_elem_blk
        istart = iend + 1
        call exgelc (exoid, idblk(i), connect(istart), ierr)
        iend = iend+num_nodes_per_elem(i)*num_elem_in_block(i)
        if (ierr.lt.0) then
          write(6,'(a)') "ERROR: cannot read elm. connectivity (exgelc)"
          STOP
        endif
60    continue
!
! read individual side sets
!
      num_sides_tot = 0
      if (num_side_sets .gt. 0) then
        call exgssi (exoid, idss, ierr)
        if (ierr.lt.0) then
          write(6,'(a)') "ERROR: cannot read SideSet ids (exgssi)"
          STOP
        endif
        maxnss=0
        do i = 1, num_side_sets
          call exgsp (exoid,idss(i),num_sides_in_set(i),idum,ierr)
          if (ierr.lt.0) then
            write(6,'(a,i3,a)') &
            "ERROR: cannot read parameters for SideSet No.",i," (exgsp)"
            STOP
          endif
          write (6, '("side set ", i2, " num_sides = ", i8)') &
                 idss(i), num_sides_in_set(i)
          num_sides_tot = num_sides_tot + num_sides_in_set(i)
          maxnss        = max(maxnss,num_sides_in_set(i))
        enddo

        ! allocate sideset arrays
        allocate (elem_list(maxnss,num_side_sets) )
        allocate (side_list(maxnss,num_side_sets) )

        do i = 1, num_side_sets
          call exgss (exoid,idss(i),elem_list(1,i),side_list(1,i),ierr)
          if (ierr.lt.0) then
            write(6,'(a,i3,a)') &
            "ERROR: cannot read parameters for SideSet No.",i," (exgss)"
            STOP
          endif
        enddo
        write (6,*)
      else
        write(6,'(a)') "WARNING: No SideSets in exodus file!"
      endif
!
! read QA records
!
      call exinq (exoid, EXQA, num_qa_rec, fdum, cdum, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read QA records (exinq QA) "
        STOP
      elseif (ierr.gt.0) then
        write(6,'(a)') "INFO: file does not contain any QA records"
      else
        call exgqa (exoid, qa_record, ierr)
        if (ierr.lt.0) then
          write(6,'(a,i3)') "WARNING: cannot read QA records (exgqa)"
        else
          write (6, '("QA records = ")')
          if (num_qa_rec.gt.10) then
            write(6,'(2a)') &
              'WARNING: Cannot handle more than 10 QA records', &
              'Printing only the first 10...'
          else
            do i = 1, num_qa_rec
              do j = 1, 4
                write (6,'(a)') qa_record(j,i)
              enddo
            enddo
          endif
        endif
      endif

      return
      end
!-----------------------------------------------------------------
      subroutine convert
!
!  Subroutine to convert an already read exodusII mesh to a nek
!  mesh. The idea is to fill each element's node coordinates of
!  size lx1**3 (3D) or lx1**2 (2D) (lx1=3) with the hex27/quad9
!  coordinates.
!
      use SIZE
      include 'exodusII.inc'

! node and face conversion (it works at least for cubit):
      integer exo_to_nek_vert3D(27)                   ! hex27 to nek numbering
      data    exo_to_nek_vert3D                     &
            / 19,  1,  7, 25, 21,  3,  9, 27, 10    &     
            ,  4, 16, 22, 20,  2,  8, 26, 12,  6    &
            , 18, 24, 14, 13, 15, 23,  5, 11, 17 /

      integer exo_to_nek_vert2D(9)                     ! quad9 to nek numbering
      data    exo_to_nek_vert2D   / 1, 3, 9, 7, 2, 6, 8, 4, 5  / 

      integer exo_to_nek_face3D(6)
      data    exo_to_nek_face3D  / 1, 5, 3, 6, 4, 2 /  ! symmetric face numbering

      integer exo_to_nek_face2D(4)
      data    exo_to_nek_face2D  / 1, 2, 3, 4 /        ! symmetric face numbering

      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '
      do iel=1,num_elem
        do ivert =1,nvert
          if (num_dim.eq.2) then
            jvert = exo_to_nek_vert2D(ivert)
            xm1(jvert,1,1,iel)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,iel)=y_exo(connect(nvert*(iel-1)+ivert))
          else
            jvert = exo_to_nek_vert3D(ivert)
            xm1(jvert,1,1,iel)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,iel)=y_exo(connect(nvert*(iel-1)+ivert))
            zm1(jvert,1,1,iel)=z_exo(connect(nvert*(iel-1)+ivert))
          endif
        enddo
      enddo
      deallocate(x_exo,y_exo,z_exo,connect)
      write(6,'(A)') 'done :: Converting elements '
!
! allocate and zero-out curve sides arrays
!
      allocate   (ccurve (4+8*(num_dim-2),num_elem) )
      allocate   (curve  (2*num_dim,12,   num_elem) )
      call rzero (curve,2*num_dim*12*num_elem)
      call blank (ccurve,(4+8*(num_dim-2))*num_elem)
!
! allocate and zero-out bc arrays only if sidesets are specified
!
      if (num_side_sets.eq.0) return   

      allocate   (cbc    (2*num_dim,      num_elem) )
      allocate   (bc     (5,2*num_dim,    num_elem) ) 
      call rzero (bc,5*2*num_dim*num_elem)
      call blank (cbc,3*2*num_dim*num_elem)
!
! set bc's
!
      write(6,'(a)') ''
      write(6,'(a)') 'Converting SideSets ...'
      do iss=1,num_side_sets
        write(6,'(a)') ''
        write(6,'(a,i2,a)') 'Sideset ',idss(iss), ' ...'
        do i=1,num_sides_in_set(iss) 
          iel = elem_list(i,iss)
          ifc = side_list(i,iss)
          if (num_dim.eq.2) then 
            jfc = exo_to_nek_face2D(ifc)
          else
            jfc = exo_to_nek_face3D(ifc)
          endif
          cbc(jfc,iel)   = 'EXO' ! dummy exodus bc 
          bc (5,jfc,iel) = idss(iss)
        enddo
      write(6,'(A,I2)') 'done :: Sideset ',idss(iss)
      enddo
      deallocate(elem_list,side_list)

      write(6,'(a)') ''
      write(6,'(a)') 'done :: Converting SideSets '

      return
      end
!--------------------------------------------------------------------
      subroutine offset_sideset()
! offset sideset number
! so sideset number is starting from 1 consecutively.
!
      use SIZE
	  
      integer ibc,ibc2,nbc
      logical newbc
      integer*8 iel,jfc
      integer bcID2

      allocate (bcID (100)) ! assuming there is no more than 100 sidesets in total
	  
      ibc = 0
      do iel= eacc_old+1,eacc
       do jfc =1,6

         if(ibc.eq.0) then
         if (bc(5,jfc,iel).ne.0) then
          ibc = ibc +1
          bcID(ibc) = bc(5,jfc,iel)
          nbc = ibc
         endif
         endif 
		 
         if(ibc.gt.0) then
         if (bc(5,jfc,iel).ne.0) then
          newbc = .TRUE.
          do ibc2 =1,nbc
           if (bc(5,jfc,iel).eq.bcID(ibc2)) then
            newbc = .FALSE.
           endif
          enddo 

          if(newbc) then		
          ibc = ibc +1		  
          bcID(ibc) = bc(5,jfc,iel)
          nbc = ibc
          endif

         endif
         endif 

       enddo
      enddo
	  
      ! sorting bcID array to ascend order 
	  
       do ibc2 =1,nbc-1
        do ibc =1,nbc-1
         if( bcID(ibc).gt.bcID(ibc+1)) then
         bcID2 = bcID(ibc)
         bcID(ibc) = bcID(ibc+1)
         bcID(ibc+1) = bcID2
         endif
        enddo
       enddo
	  
      ! reoder sideset number 
	  
      do ibc =1,nbc
       write(6,*) 'offset sideset number: ',bcID(ibc),'->',ibc
       do iel= eacc_old+1,eacc
       do jfc =1,6
         if (bc(5,jfc,iel).eq.bcID(ibc)) then  
           bc(5,jfc,iel) = ibc
         endif		 
       enddo
       enddo
      enddo
	  
      deallocate ( bcID )
	  
      return
      end
!--------------------------------------------------------------------
      subroutine gather_bc_info
      use SIZE

      integer ibc,ibc2,nbc
      logical newbc
      integer*8 iel,jfc
      integer bcID2
      
      allocate (bcID (num_side_sets))
!      allocate (bcID (1000))
	  
      ibc = 0
	  
      write(6,*) 'calling: gather_bc_info()'
	  
      ! gather all boundary information
	  
      do iel= 1,num_elem
       do jfc =1,2*num_dim

         if(ibc.eq.0) then
         if (bc(5,jfc,iel).ne.0) then
          ibc = ibc +1
          bcID(ibc) = bc(5,jfc,iel)
          nbc = ibc		  
         endif
         endif 
		 
         if(ibc.gt.0) then
         if (bc(5,jfc,iel).ne.0) then
          newbc = .TRUE.
          do ibc2 =1,nbc
           if (bc(5,jfc,iel).eq.bcID(ibc2)) then
            newbc = .FALSE.
           endif
          enddo 

          if(newbc) then		
          ibc = ibc +1		  
          bcID(ibc) = bc(5,jfc,iel)
          nbc = ibc
          endif

         endif
         endif 

       enddo
      enddo
	  
      ! sorting bcID array to ascend order 
	  
       do ibc2 =1,nbc-1
        do ibc =1,nbc-1
         if( bcID(ibc).gt.bcID(ibc+1)) then
         bcID2 = bcID(ibc)
         bcID(ibc) = bcID(ibc+1)
         bcID(ibc+1) = bcID2
         endif
        enddo
       enddo
      bcNumber = nbc

      write(6,*) 'done: gather_bc_info()'
	  
      return
      end
	  
!--------------------------------------------------------------------
      subroutine scale_mesh
      use SIZE
      real*8 xx,yy,zz,ss
      integer ie
      write(6,*) "please input scaling factor (1 for no scale):"
      read(5,*) ss
	  
      if (ss.eq.1) return
	  
      do ie = 1,num_elem
        do i = 1,27
        xx = xm1(i,1,1,ie)        
        yy = ym1(i,1,1,ie)
        zz = zm1(i,1,1,ie)
  
        xm1(i,1,1,ie) = xx*ss
        ym1(i,1,1,ie) = yy*ss
        zm1(i,1,1,ie) = zz*ss
        enddo
      enddo
	  
      call checkXYZ_min_max()

      return
      end
!--------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine cross_product(AB_v,AC_v,prod_v,area)
! calculate cross product of two vectors

      real*8 AB_v(3),AC_v(3),prod_v(3)
      real*8 area

      prod_v(1) = AB_v(2)*AC_v(3)-AB_v(3)*AC_v(2)
      prod_v(2) = AB_v(3)*AC_v(1)-AB_v(1)*AC_v(3)
      prod_v(3) = AB_v(1)*AC_v(2)-AB_v(2)*AC_v(1)

      area = prod_v(1)**2 + prod_v(2)**2+prod_v(3)**2
      area = sqrt(area)
      area = area/2.0
 
      return 
      end
!------------------------------------------------------------------------------------------
!--------------------------------------------------------------------
      subroutine gen_re2

      use SIZE

      write(6,*)
      write(6,'(A,A)') 'writing ', re2name

      call open_re2
      call write_xyz
      call write_curve
      call write_bc
      call close_re2

      return
      end
!--------------------------------------------------------------------
      subroutine open_re2

      use SIZE

      character(80) hdr
      integer nBCre2

      real*4 test
      data   test  / 6.54321 /

      call byte_open(re2name,ierr)
            
      num_elem = etot
      !num_dim = 3
      nBCre2 = 1
      if (num_elem.ne.eftot)  nBCre2  = 2

!  Write the header
      call blank   (hdr,80)    
      write(hdr,1) num_elem, num_dim, eftot, nBCre2
!      write(hdr,1) num_elem, num_dim, eftot
    1 format('#v004',i16,i3,i16,i4,' hdr')
!    1 format('#v002',i16,i3,i16,i4,' hdr')
!    1 FORMAT('#v002',i9,i3,i9,' this is the hdr')
      call byte_write(hdr,20,ierr)         
      call byte_write(test,1,ierr)     ! write the endian discriminator

      return
      end
!--------------------------------------------------------------------
      subroutine write_xyz

      use SIZE

      integer*8 iel
      real*8     xx(8), yy(8), zz(8)
      real*8   rgroup, buf2(30)

      integer isym2pre(8)   ! Symmetric-to-prenek vertex ordering
      data    isym2pre    / 1 , 2 , 4 , 3 , 5 , 6 , 8 , 7 /

      nxs = 3-1  ! always nx1=ny1=nz1=3 here
      nys = 3-1
      nzs = 3-1

      igr    = 0
      rgroup = igr

      do iel=1,num_elem
        l = 0
        if (num_dim.eq.3) then
          do k=0,1
          do j=0,1
          do i=0,1
            l=l+1
            li=isym2pre(l)
            xx(li) = xm1(1+i*nxs,1+j*nys,1+k*nzs,iel)
            yy(li) = ym1(1+i*nxs,1+j*nys,1+k*nzs,iel)
            zz(li) = zm1(1+i*nxs,1+j*nys,1+k*nzs,iel)
          enddo
          enddo
          enddo
        else
          do j=0,1
          do i=0,1
            l=l+1
            li=isym2pre(l)
            xx(li) = xm1(1+i*nxs,1+j*nys,1,iel)
            yy(li) = ym1(1+i*nxs,1+j*nys,1,iel)
          enddo
          enddo
        endif

        call byte_write(rgroup, 2,ierr)

        if (num_dim.eq.3) then
          call copy        (buf2(1) ,xx,8)
          call copy        (buf2(9) ,yy,8)
          call copy        (buf2(17),zz,8)
          call byte_write  (buf2(1) ,16, ierr)
          call byte_write  (buf2(9) ,16, ierr)
          call byte_write  (buf2(17),16, ierr)
        else
          call copy        (buf2(1),xx,4)
          call copy        (buf2(5),yy,4)
          call byte_write  (buf2(1),8, ierr)
          call byte_write  (buf2(5),8, ierr)
        endif
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine write_curve

      use SIZE
      integer*8 iel,iedge,ncurv,nedge
      real*8     buf2(30)
      real*8     rcurve

      character(1) cc

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


      do iel=1,num_elem
         call gen_rea_midside_e(iel)
      enddo

      nedge = 4 + 8*(num_dim-2)
      ncurv = 0
      do iel=1,num_elem
        do iedge=1,nedge
          if (ccurve(iedge,iel).ne.' ') ncurv = ncurv+1
        enddo
      enddo

      rcurve = ncurv
      call byte_write(rcurve,2, ierr)
      write(6,*) 'curves: ',ncurv

      do iel=1,num_elem
        do iedge=1,nedge
          if (ccurve(iedge,iel).ne.' ') then
            if (ccurve(iedge,iel).eq.'C') cc='C'
            if (ccurve(iedge,iel).eq.'s') cc='s'
            if (ccurve(iedge,iel).eq.'m') cc='m'
            buf2(1) = iel
            buf2(2) = iedge
            call copy       (buf2(3),curve(1,iedge,iel),5)
            call blank      (buf2(8),8)
            call chcopy     (buf2(8),cc,1)
            call byte_write (buf2,16,ierr)
          endif
        enddo
      enddo

!      nedge = 4 + 8*(num_dim-2)
!      ncurv = 0
!
!      do iel=1,num_elem
!        do iedge=1,nedge
!
!            do j=1,3
!               xyz(1,j) = xm1(e3(j,iedge),1,1,iel)
!               xyz(2,j) = ym1(e3(j,iedge),1,1,iel)
!               xyz(3,j) = zm1(e3(j,iedge),1,1,iel)
!            enddo
!
!            call average2vec(mide(1), xyz(1,1), xyz(1,3))  
!
!            call distance(xyz(1,1),xyz(1,3),dist0)
!            call distance(mide(1),xyz(1,2),dist1)
!
!            if (dist1.gt.(1e-3*dist0)) then
!
!            ncurv = ncurv + 1
!   
!            endif
!
!        enddo
!      enddo
!	  
!      rcurve = ncurv
!      call byte_write(rcurve,2, ierr)
!      write(6,*) 'curves: ',ncurv
!
!      do iel=1,num_elem
!        do iedge=1,nedge
!
!            do j=1,3
!               xyz(1,j) = xm1(e3(j,iedge),1,1,iel)
!               xyz(2,j) = ym1(e3(j,iedge),1,1,iel)
!               xyz(3,j) = zm1(e3(j,iedge),1,1,iel)
!            enddo
!
!            call average2vec(mide(1), xyz(1,1), xyz(1,3))  
!
!            call distance(xyz(1,1),xyz(1,3),dist0)
!            call distance(mide(1),xyz(1,2),dist1)
!
!            if (dist1.gt.(1e-3*dist0)) then
!
!            curve(1,iedge,iel) = xm1(edge_mid(iedge),1,1,iel)
!            curve(2,iedge,iel) = ym1(edge_mid(iedge),1,1,iel)
!            curve(3,iedge,iel) = zm1(edge_mid(iedge),1,1,iel)
!            cc='m'
!            buf2(1) = iel
!            buf2(2) = iedge
!            call copy       (buf2(3),curve(1,iedge,iel),5)
!            call blank      (buf2(8),8)
!            call chcopy     (buf2(8),cc,1)
!            call byte_write (buf2,16,ierr)
!   
!            endif
!
!        enddo
!      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine write_bc
      
      use SIZE
      integer*8 iel,ibc,ifc,nbc,nface
      real*8  rbc, buf2(30)

      character(3) ch3
      character(1) chdum
      data         chdum /' '/

      if (num_side_sets.eq.0) return

      !rbc = num_sides_tot
      !call byte_write (rbc,2,ierr)

      nbc   = 0
      nface = 2*num_dim
      do iel=1,eftot
        do ifc=1,nface
          if (cbc(ifc,iel).ne.'   ')  nbc = nbc + 1
        enddo
      enddo
      rbc = nbc
      call byte_write (rbc,2, ierr)
      write(6,*) 'velocity boundary faces: ',nbc


      do iel = 1,eftot
        do ifc = 1,2*num_dim
          ch3 = cbc(ifc,iel)
          if (ch3.ne.'   ') then
            buf2(1)=iel
            buf2(2)=ifc
            call copy   (buf2(3),bc(1,ifc,iel),5)
            call blank  (buf2(8),8)
            call chcopy (buf2(8),ch3,3)
            if (eftot.ge.1000000) then
              ibc     = bc(1,ifc,iel)
              buf2(3) = ibc
            endif
            call byte_write (buf2,16,ierr)
          endif
        enddo
      enddo

     if(num_elem.ne.eftot) then

! writing thermal bc to all elements
      nbc   = 0
      nface = 2*num_dim
      do iel=1,num_elem
        do ifc=1,nface
          if (cbc(ifc,iel).ne.'   ')  nbc = nbc + 1
        enddo
      enddo
      rbc = nbc
      call byte_write (rbc,2, ierr)
      write(6,*) 'thermal boundary faces: ',nbc

      do iel = 1,num_elem
        do ifc = 1,2*num_dim
          ch3 = cbc(ifc,iel)
          if (ch3.ne.'   ') then
            buf2(1)=iel
            buf2(2)=ifc
            call copy   (buf2(3),bc(1,ifc,iel),5)
            call blank  (buf2(8),8)
            call chcopy (buf2(8),ch3,3)
            if (num_elem.ge.1000000) then
              ibc     = bc(1,ifc,iel)
              buf2(3) = ibc
            endif
            call byte_write (buf2,16,ierr)
          endif
        enddo
      enddo

      endif

      return
      end 
!-----------------------------------------------------------------------
      subroutine close_re2

      write(6,*) 'calling close_re2'
      call byte_close (ierr)

      return
      end
!-----------------------------------------------------------------------
      subroutine gen_rea_midside_e(e)

      use SIZE

      real         len
      real         x3(27),y3(27),z3(27),xyz(3,3)
      character(1) ccrve(12)
      integer*8      e,edge

      integer e3(3,12)
      save    e3
      data    e3 /  1, 2, 3,    3, 6, 9,    9, 8, 7,    7, 4, 1, &
                   19,20,21,   21,24,27,   27,26,25,   25,22,19, &
                    1,10,19,    3,12,21,    9,18,27,    7,16,25  /

      call chcopy(ccrve,ccurve(1,e),12)

      call map2reg                   (x3,3,xm1(1,1,1,e),1)  ! Map to 3x3x3 array
      call map2reg                   (y3,3,ym1(1,1,1,e),1)
      if (num_dim.eq.3) call map2reg (z3,3,zm1(1,1,1,e),1)

!     Take care of spherical curved face defn
      if (ccurve(5,e).eq.'s') then
         call chcopy (ccrve(1),'ssss',4) ! face 5
         call chcopy (ccrve(5),' ',1)    ! face 5
      endif
      if (ccurve(6,e).eq.'s') then
         call chcopy (ccrve(5),'ssss',4) ! face 6
      endif

      tol   = 1.e-4
      tol2  = tol**2
      nedge = 4 + 8*(num_dim-2)

      do i=1,nedge
         if (ccrve(i).eq.' ') then
            do j=1,3
               xyz(1,j) = x3(e3(j,i))
               xyz(2,j) = y3(e3(j,i))
               xyz(3,j) = z3(e3(j,i))
            enddo
            len = 0.
            h   = 0.
            do j=1,num_dim
               xmid = .5*(xyz(j,1)+xyz(j,3))
               h    = h   + (xyz(j,2)-xmid)**2
               len  = len + (xyz(j,3)-xyz(j,1))**2
            enddo
            if (h.gt.tol2*len) ccurve(i,e) = 'm'
            if (h.gt.tol2*len) call copy (curve(1,i,e),xyz(1,2),num_dim)
         endif
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine map2reg(ur,n,u,nel)
!
!     Map scalar field u() to regular n x n x n array ur

      use SIZE

      real    ur(1), u(3*3*3,1)
      integer e

      ldr = n**num_dim

      k=1
      do e=1,nel
         if (num_dim.eq.2) call map2reg_2di_e(ur(k),n,u(1,e),3)
         if (num_dim.eq.3) call map2reg_3di_e(ur(k),n,u(1,e),3)
         k = k + ldr
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine map2reg_2di_e(uf,n,uc,m) ! Fine, uniform pt

      real uf(n,n),uc(m,m)

      parameter (l=50)
      common /cmap2d/ j(l*l),jt(l*l),w(l*l),z(l)

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no ) then

          call zwgll (z,w,m)
          call zuni  (w,n)

          call gen_int_gz(j,jt,w,n,z,m)

      endif

      call mxmf2(j,n,uc,m,w ,m)
      call mxmf2(w,n,jt,m,uf,n)

      return
      end
!-----------------------------------------------------------------------
      subroutine map2reg_3di_e(uf,n,uc,m) ! Fine, uniform pt

      real uf(n,n,n),uc(m,m,m)

      parameter (l=50)
      common /cmap3d/ j(l*l),jt(l*l),v(l*l*l),w(l*l*l),z(l)

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no ) then

          call zwgll (z,w,m)
          call zuni  (w,n)

          call gen_int_gz(j,jt,w,n,z,m)

      endif

      mm = m*m
      mn = m*n
      nn = n*n

      call mxmf2(j,n,uc,m,v ,mm)
      iv=1
      iw=1
      do k=1,m
         call mxmf2(v(iv),n,jt,m,w(iw),n)
         iv = iv+mn
         iw = iw+nn
      enddo
      call mxmf2(w,nn,jt,m,uf,n)

      return
      end
!-----------------------------------------------------------------------
      subroutine gen_int_gz(j,jt,g,n,z,m)

!     Generate interpolater from m z points to n g points

!        j   = interpolation matrix, mapping from z to g
!        jt  = transpose of interpolation matrix
!        m   = number of points on z grid
!        n   = number of points on g grid

      real j(n,m),jt(m,n),g(n),z(m)

      mpoly  = m-1
      do i=1,n
         call fd_weights_full(g(i),z,mpoly,0,jt(1,i))
      enddo

      call transpose(j,n,jt,m)

      return
      end
!-----------------------------------------------------------------------
      SUBROUTINE BLANK(A,N)
      INTEGER*8 N
      CHARACTER(1) A(1)
      CHARACTER(1) BLNK
      SAVE        BLNK
      DATA        BLNK /' '/
 
      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE
      RETURN
      END
!-----------------------------------------------------------------------
      subroutine exitt

      stop
      return
      end
!-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1),b(1)
      integer*8 n
      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine chcopy(a,b,n)
      CHARACTER(1) A(1), B(1)
      integer*8 n
 
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
!-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      INTEGER A(1), B(1)
      integer*8 n
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
!-----------------------------------------------------------------------
      subroutine fd_weights_full(xx,x,n,m,c)
!
!     This routine evaluates the derivative based on all points
!     in the stencils.  It is more memory efficient than "fd_weights"
!
!     This set of routines comes from the appendix of
!     A Practical Guide to Pseudospectral Methods, B. Fornberg
!     Cambridge Univ. Press, 1996.   (pff)
!
!     Input parameters:
!       xx -- point at wich the approximations are to be accurate
!       x  -- array of x-ordinates:   x(0:n)
!       n  -- polynomial degree of interpolant (# of points := n+1)
!       m  -- highest order of derivative to be approxxmated at xi
!
!     Output:
!       c  -- set of coefficients c(0:n,0:m).
!             c(j,k) is to be applied at x(j) when
!             the kth derivative is approxxmated by a
!             stencil extending over x(0),x(1),...x(n).
!
!
      real x(0:n),c(0:n,0:m)
 
      c1       = 1.
      c4       = x(0) - xx
 
      do k=0,m
      do j=0,n
         c(j,k) = 0.
      enddo
      enddo
      c(0,0) = 1.
 
      do i=1,n
         mn = min(i,m)
         c2 = 1.
         c5 = c4
         c4 = x(i)-xx
         do j=0,i-1
            c3 = x(i)-x(j)
            c2 = c2*c3
            do k=mn,1,-1
               c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
            enddo
            c(i,0) = -c1*c5*c(i-1,0)/c2
            do k=mn,1,-1
               c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
            enddo
            c(j,0) = c4*c(j,0)/c3
         enddo
         c1 = c2
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine transpose(a,lda,b,ldb)
      real a(lda,1),b(ldb,1)
 
      do j=1,ldb
         do i=1,lda
            a(i,j) = b(j,i)
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      function ltrunc(string,l)
      CHARACTER(1) STRING(L)
      CHARACTER(1)   BLNK
      DATA BLNK/' '/
 
      DO 100 I=L,1,-1
         L1=I
         IF (STRING(I).NE.BLNK) GOTO 200
  100 CONTINUE
      L1=0
  200 CONTINUE
      LTRUNC=L1
      return
      END
!-----------------------------------------------------------------------
      subroutine zuni(z,np)
!
!     Generate equaly spaced np points on the interval [-1:1]
!
      real z(1)

      dz = 2./(np-1)
      z(1) = -1.
      do i = 2,np-1
         z(i) = z(i-1) + dz
      enddo
      z(np) = 1.

      return
      end
!-----------------------------------------------------------------------
      subroutine rzero(A,N)
      integer*8 N,I
      real A(1)
      DO 100 I = 1, N
 100     A(I ) = 0.0
      return
      END
!-----------------------------------------------------------------------
      subroutine rzero8(A,N)
      integer*8 N,I
      real*8 A(1)
      DO 100 I = 1, N
 100     A(I ) = 0.0
      return
      END
!-----------------------------------------------------------------------
      subroutine rzero_int(A,N)
      integer*8 N,I
      integer A(1)
      DO 100 I = 1, N
 100     A(I) = 0
      return
      END
!-----------------------------------------------------------------------
      subroutine rzero_int2(A,N)
      integer N,I
      integer A(1)
      DO 100 I = 1, N
 100     A(I) = 0
      return
      END
