      subroutine iniproc
      include 'SIZE'
      include 'PARALLEL'
      include 'mpif.h'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

      logical flag

      call mpi_initialized(mpi_is_initialized, ierr) !  Initialize MPI
      if ( mpi_is_initialized .eq. 0 ) then
         call mpi_init (ierr)
      endif

      ! create communicator
      call init_nek_comm
      np  = np_
      nid = nid_

      if(nid.eq.0) call printHeader

      ! check upper tag size limit
      call mpi_attr_get(MPI_COMM_WORLD,MPI_TAG_UB,nval,flag,ierr)
      if (nval.lt.(10000+max(lp,lelg))) then
         if(nid.eq.0) write(6,*) 'ABORT: MPI_TAG_UB too small!'
         call exitt
      endif

      IF (NP.GT.LP) THEN
         WRITE(6,*) 
     $   'ERROR: Code compiled for a max of',LP,' processors.'
         WRITE(6,*) 
     $   'Recompile with LP =',NP,' or run with fewer processors.'
         WRITE(6,*) 
     $   'Aborting in routine INIPROC.'
         call exitt
      endif

      ! set word size for REAL
      wdsize=4
      eps=1.0e-12
      oneeps = 1.0+eps
      if (oneeps.ne.1.0) then
         wdsize=8
      else
         if(nid.eq.0) 
     &     write(6,*) 'ABORT: single precision mode not supported!'
         call exitt
      endif
      nekreal = mpi_real
      if (wdsize.eq.8) nekreal = mpi_double_precision

      ifdblas = .false.
      if (wdsize.eq.8) ifdblas = .true.

      ! set word size for INTEGER
      ! HARDCODED since there is no secure way to detect an int overflow
      isize = 4

      ! set word size for LOGICAL
      lsize = 4

      ! set word size for CHARACTER
      csize = 1
c
      PID = 0
      NULLPID=0
      NODE0=0
      NODE= NID+1

      if (nid.eq.0) then 
         write(6,*) 'Number of processors:',np
         WRITE(6,*) 'REAL    wdsize      :',WDSIZE
         WRITE(6,*) 'INTEGER wdsize      :',ISIZE
      endif

      call crystal_setup(cr_h,nekcomm,np)  ! set cr handle to new instance

      return
      end
c-----------------------------------------------------------------------
      subroutine init_nek_comm
      include 'mpif.h'
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
C
      call create_comm(nekcomm) ! set up nekton specific communicator
c
      nid_  = mynode()
      np_   = numnodes()
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gop( x, w, op, n)
c
c     Global vector commutative operation
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c
      real x(n), w(n)
      character*3 op
c
      if (op.eq.'+  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_prod,nekcomm,ierr)
      else
         write(6,*) nid,' OP ',op,' not supported.  ABORT in GOP.'
         call exitt
      endif

      call copy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine igop( x, w, op, n)
c
c     Global vector commutative operation
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      integer x(n), w(n)
      character*3 op

      if     (op.eq.'+  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_prod,nekcomm,ierr)
      else
        write(6,*) nid,' OP ',op,' not supported.  ABORT in igop.'
        call exitt
      endif

      call icopy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine i8gop( x, w, op, n)
c
c     Global vector commutative operation
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      integer*8 x(n), w(n)
      character*3 op

      if     (op.eq.'+  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_prod,nekcomm,ierr)
      else
        write(6,*) nid,' OP ',op,' not supported.  ABORT in igop.'
        call exitt
      endif

      call i8copy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine csend(mtype,buf,len,jnid,jpid)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real*4 buf(1)

      call mpi_send (buf,len,mpi_byte,jnid,mtype,nekcomm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine crecv(mtype,buf,lenm)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
C
      real*4 buf(1)
      len = lenm
      jnid = mpi_any_source

      call mpi_recv (buf,len,mpi_byte
     $              ,jnid,mtype,nekcomm,status,ierr)
c
      if (len.gt.lenm) then 
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine crecv3(mtype,buf,len,lenm)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
C
      real*4 buf(1)
      len = lenm
      jnid = mpi_any_source

      call mpi_recv (buf,len,mpi_byte
     $            ,jnid,mtype,nekcomm,status,ierr)
      call mpi_get_count (status,mpi_byte,len,ierr)
c
      if (len.gt.lenm) then 
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      integer function numnodes()
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      call mpi_comm_size (nekcomm, numnodes , ierr)

      return
      end
c-----------------------------------------------------------------------
      integer function mynode()
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer myid

      call mpi_comm_rank (nekcomm, myid, ierr)
      mynode = myid

      return
      end
c-----------------------------------------------------------------------
      real*8 function dnekclock()
      include 'mpif.h'
c
      dnekclock=mpi_wtime()
c
      return
      end
c-----------------------------------------------------------------------
      real*8 function dnekclock_sync()
      include 'mpif.h'
c
      call gsync()
      dnekclock_sync=mpi_wtime()
c
      return
      end
c-----------------------------------------------------------------------
      subroutine lbcast(ifif)
C
C     Broadcast logical variable to all processors.
C
      include 'SIZE'
      include 'PARALLEL'
      include 'mpif.h'

      logical ifif

      if (np.eq.1) return

      item=0
      if (ifif) item=1
      call bcast(item,isize)
      ifif=.false.
      if (item.eq.1) ifif=.true.

      return
      end
c-----------------------------------------------------------------------
      subroutine bcast(buf,len)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real*4 buf(1)

      call mpi_bcast (buf,len,mpi_byte,0,nekcomm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine create_comm(icomm)
      include 'mpif.h'

      call mpi_comm_group (mpi_comm_world,itmp,ierr)
      call mpi_comm_create (mpi_comm_world,itmp,icomm,ierr)
      call mpi_group_free (itmp,ierr)
c     write(6,*) 'nekcomm:',nekcomm

      return
      end
c-----------------------------------------------------------------------
      function isend(msgtag,x,len,jnid,jpid)
c
c     Note: len in bytes
c
      integer x(1)
C
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
C
      call mpi_isend (x,len,mpi_byte,jnid,msgtag
     $       ,nekcomm,imsg,ierr)
      isend = imsg
c     write(6,*) nid,' isend:',imsg,msgtag,len,jnid,(x(k),k=1,len/4)
c
      return
      end
c-----------------------------------------------------------------------
      function irecv(msgtag,x,len)
c
c     Note: len in bytes
c
      integer x(1)
C
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
C
      call mpi_irecv (x,len,mpi_byte,mpi_any_source,msgtag
     $       ,nekcomm,imsg,ierr)
      irecv = imsg
c     write(6,*) nid,' irecv:',imsg,msgtag,len
c
c
      return
      end
c-----------------------------------------------------------------------
      subroutine msgwait(imsg)
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
c
c     write(6,*) nid,' msgwait:',imsg
c
      call mpi_wait (imsg,status,ierr)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gsync()

      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      call mpi_barrier(nekcomm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine exitti(stringi,idata)
      character*1 stringi(132)
      character*1 stringo(132)
      character*11 s11
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      call blank(stringo,132)
      call chcopy(stringo,stringi,132)
      len = indx1(stringo,'$',1)
      write(s11,11) idata
   11 format(1x,i10)
      call chcopy(stringo(len),s11,11)

      if (nid.eq.0) write(6,1) (stringo(k),k=1,len+10)
    1 format('EXIT: ',132a1)

      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine err_chk(ierr,string)
      character*1 string(132)
      character*1 ostring(132)
      character*10 s10
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      ierr = iglsum(ierr,1)
      if(ierr.eq.0) return 

      len = indx1(string,'$',1)
      call blank(ostring,132)
      write(s10,11) ierr
   11 format(1x,' ierr=',i3)

      call chcopy(ostring,string,len-1)
      call chcopy(ostring(len),s10,10)

      if (nid.eq.0) write(6,1) (ostring(k),k=1,len+10)
    1 format('ERROR: ',132a1)

      call exitt

      return
      end
c


c-----------------------------------------------------------------------
      subroutine exitt
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'mpif.h'

      real*4 papi_mflops
      integer*8 papi_flops
c
      call gsync()

#ifdef PAPI
      call nek_flops(papi_flops,papi_mflops)
#endif

      tstop  = dnekclock()
      ttotal = tstop-etimes
      nxyz   = nx1*ny1*nz1

      if (nid.eq.0) then 
         dtmp1 = 0
         dtmp2 = 0
         dtmp3 = 0
         if(istep.gt.0) then
           dgp   = max(nvtot,1)
           dtmp1 = np*ttime/(dgp*max(istep,1))
           dtmp2 = ttime/max(istep,1)
           dtmp3 = 1.*papi_flops/1e6
         endif 
         write(6,*) ' '
         write(6,'(A)') 'call exitt: dying ...'
         write(6,*) ' '
         call print_stack()
         write(6,*) ' '
         write(6,'(4(A,1p1e13.5,A,/))') 
     &       'total elapsed time             : ',ttotal, ' sec'
     &      ,'total solver time incl. I/O    : ',ttime , ' sec'
     &      ,'time/timestep                  : ',dtmp2 , ' sec'
     &      ,'CPU seconds/timestep/gridpt    : ',dtmp1 , ' sec'
#ifdef PAPI
         write(6,'(2(A,1g13.5,/))') 
     &       'Gflops                         : ',dtmp3/1000.
     &      ,'Gflops/s                       : ',papi_mflops/1000.
#endif
      endif 
      call flush_io

      call mpi_finalize (ierr)
      call exit(0)

      return
      end
c-----------------------------------------------------------------------
      subroutine printHeader

      INCLUDE 'HEADER'

      return
      end
c-----------------------------------------------------------------------
      function igl_running_sum(in)
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
      integer x,w,r

      x = in  ! running sum
      w = in  ! working buff
      r = 0   ! recv buff

      call mpi_scan(x,r,1,mpi_integer,mpi_sum,nekcomm,ierr)
      igl_running_sum = r

      return
      end

   

