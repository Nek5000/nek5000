c-----------------------------------------------------------------------
      program re2torea

c     this program will take a nek all-ascii rea file and
c     extract geometry, bcs, and curve data to a new binary re2 file, plus
c     an ascii rea file for just the parameters

      
      character*80 fin
      character*1  fin1(80)
      equivalence (fin1,fin)

      character*80 fout
      character*1  fout1(80)
      equivalence (fout1,fout)

      character*80 fbin
      character*1  fbin1(80)
      equivalence (fbin1,fbin)

      character*80 string
      character*1  string1(80)
      equivalence (string,string1)
 
      character*80 sstring

      parameter(nelm=MAXNEL)
      real*8 x(8),y(8),z(8)
      real*8 bc(5,6,nelm)
      integer id, jd, wdsizi
      logical ifflow,ifheat

      character*80 hdr
      real*4 test
      character*5 version

      integer buf(30)

      real*8 rbc,rcurve
      character*3 cbc(6,nelm)
      logical ifbswap,if_byte_swap_test,ifexist,strcmp
      
c     for workstation:
      in = 5
 
c     Get file name
      write(6,*) 'Input .rea & .re2 name:'      
      call blank(fin,80)
      read(in,80) fin
      len = ltrunc(fin,80)
   80 format(a80)
 
c     Get file name
      write(6,*) 'Input .rea output name:'      
      call blank(fout,80)
      read(in,80) fout
      lou = ltrunc(fout,80)

      if (strcmp(fin,fout)) then
         write(6,*) 'Abort: filenames cannot be the same'
         stop
      endif

      fbin = fin
      call chcopy(fin1(len+1),'.rea',4)
      call chcopy(fout1(lou+1),'.rea',4)
      call chcopy(fbin1(len+1),'.re2'//CHAR(0),5)

      ifexist = .false.
      inquire(file=fin, exist=ifexist) 
      if(.not.ifexist)then
        write(6,*)'could not find file ',(fin1(i),i=1,len+4),', abort.'
        stop
      endif
      ifexist = .false.
      inquire(file=fbin, exist=ifexist) 
      if(.not.ifexist)then
        write(6,*)'could not find file ',(fbin1(i),i=1,len+4),', abort.'
        stop
      endif

      open(unit=10, file=fin)
      open(unit=11, file=fout)
      call byte_open(fbin,ierr)

      npscal = 0
      call scanout(sstring,'PARAMETERS FOLL',15,10,11)
      read(sstring,*) nparam
      write(11,80) sstring
      do ii=1,nparam
         read(10,80) sstring
         if(ii.eq.23) read(sstring,*) rnpscal
         if(ii.eq.32) read(sstring,*) rBCrea
         write(11,80) sstring
      enddo
      npscal=rnpscal
      nBCrea=rBCrea

      call scanout(sstring,'LOGICAL',7,10,11)
      read(sstring,*) nlogical
      write(11,80) sstring

      ifflow=.false.
      ifheat=.false. 
      do ii = 1,nlogical
         read(10,80) sstring
         if(indx1(sstring,'IFFLOW',6).ne.0) read(sstring,*) ifflow
         if(indx1(sstring,'IFHEAT',6).ne.0) read(sstring,*) ifheat
         write(11,80) sstring
      enddo
            

      call blank(hdr,80)
      call byte_read(hdr,20,ierr)

    1 format(a5)
    2 format(a5,i9,i3,i9)
    3 format(a5,i16,i3,i16,i4)
 
      read (hdr,1) version

      if (version.eq.'#v001') then
         wdsizi = 4
         read (hdr,2) version,nelt,ndim,nelv
      elseif (version.eq.'#v002') then
         wdsizi = 8
         read (hdr,2) version,nelt,ndim,nelv
      elseif (version.eq.'#v003') then
         wdsizi = 8
         read (hdr,2) version,nelt,ndim,nelv
      elseif (version.eq.'#v004') then
         wdsizi = 8
         read (hdr,3) version,nelt,ndim,nelv,nBCre2
         if (nBCrea.ne.nBCre2) then
            write(6,4) nBCrea,nBCre2
            call exitt
         endif
    4    format('Error: nBC mismatched! rea=',i3,' re2=',i3)
      else
         write(6,*) 'Error: Unsupported version number',version
         call exitt
      endif

      nelt=abs(nelt)

      write(11,*) 
     & ' 30.0000  20.0000  -6.00000  -10.0000 XFAC,YFAC,XZERO,YZERO'
      write(11,*) 
     & '**MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8'

      write(11,*) nelt,ndim,nelv

      call byte_read(test,1,ierr)
      ifbswap=if_byte_swap_test(test)

      print *, nelt, ndim, nelv


c mesh

      do ie=1,nelt

         nwds = (1 + ndim*(2**ndim))*(wdsizi/4)
         call byte_read(buf,nwds,ierr)
         call buf_to_xyz(buf,igroup,x,y,z,ifbswap,ndim,wdsizi)

         if(nelt.lt.100000) then
           write (11,12) ie, ie, 'a', igroup
         else
           write (11,'(A,I12,A,I1)') 
     &        '  ELEMENT ', ie, '  GROUP  ', igroup
         endif
         write (11,15)   (x(k),k=1,4)
         write (11,15)   (y(k),k=1,4)
 
         if(ndim.eq.3) then 
           write (11,15)   (z(k),k=1,4)
           write (11,15)   (x(k),k=5,8)
           write (11,15)   (y(k),k=5,8)
           write (11,15)   (z(k),k=5,8)
         endif
      enddo

 12   format(
     $     '            ELEMENT',i12,' [',i5,a1,']    GROUP     ',i1)

 15   format(4g14.6)

c curved sides

      if(wdsizi.eq.8) then
         call byte_read(rcurve,2,ierr)
         if(ifbswap) call byte_reverse8(rcurve,2,ierr)
         ncurve=rcurve
      else
         call byte_read(ncurve,1,ierr)
         if(ifbswap) call byte_reverse(ncurve,nwds,ierr)
      endif
      write(6,*) 'number of curved sides', ncurve

      write(11,13)
      write(11,14) ncurve
 13   format(' ***** CURVED SIDE DATA *****')
 14   format(i12
     $     ,' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE')


      nwds = (2 + 1 + 5)*(wdsizi/4) 

      IF (NCURVE.GT.0) THEN
         DO 50 ICURVE=1,NCURVE
            call byte_read(buf,nwds,ierr)
            call buf_to_curve(buf,nelt,ifbswap,wdsizi)
 50      CONTINUE
      endif


c boundary conditions
      nface=4
      if(ndim.eq.3) nface=6

      init_bc=2
      if(ifflow) init_bc=1
      nfld=1
      if(ifheat) nfld=2+npscal

      write(11,*) '***** BOUNDARY CONDITIONS *****'

      do ifld=init_bc,nfld
         if(ifld.eq.1) write(11,*) 
     &    '***** FLUID   BOUNDARY CONDITIONS *****'
         if(ifld.eq.2.and.init_bc.eq.2) write(11,*) 
     &    '***** NO FLUID   BOUNDARY CONDITIONS *****'
         if(ifld.eq.2) write(11,*) 
     &    '***** THERMAL  BOUNDARY CONDITIONS *****'
         if(ifld.gt.2) write(11,*) 
     &    '***** PASSIVE SCALAR  BOUNDARY CONDITIONS *****'

         nel=nelv
         if(ifld.ge.2) nel=nelt

         do ie = 1,nel
            do k=1,nface
               cbc(k,ie) = 'E  '
               do j=1,5
                 bc(j,k,ie) =0.0
               enddo
            enddo
         enddo

         if(wdsizi.eq.8) then
           call byte_read(rbc,2,ierr)
           if (ifbswap) call byte_reverse8(rbc,2,ierr) 
           nbc = rbc
         else
           call byte_read(nbc,1,ierr)
           if (ifbswap) call byte_reverse(nbc,1,ierr)
         endif

         write(6,*) ifld, nbc, 'number of bc'

         do ie = 1,nbc
            call byte_read(buf,nwds,ierr)
            call buf_to_bc(cbc,bc,buf,nelm,ifbswap,wdsizi)
         enddo

         do ie = 1,nel
            do k=1,nface

              if (nelt.lt.1 000) then
                  write (11,20) cbc(k,ie),ie,k,(bc(j,k,ie),j=1,5)
              elseif(nelt.lt.100 000) then
                  write (11,21) cbc(k,ie),ie,k,(bc(j,k,ie),j=1,5)
              elseif(nelt.lt.1 000 000) then
                  write (11,22) cbc(k,ie),ie,(bc(j,k,ie),j=1,5)
              else
                  write (11,23) cbc(k,ie),ie,(bc(j,k,ie),j=1,5)
              endif
            enddo
         enddo
      enddo

      if(nfld.eq.1.and.init_bc.eq.1) write(11,*) 
     &    '***** NO THERMAL BOUNDARY CONDITIONS *****'

   20 FORMAT(1x,A3,2I3,5G14.6)
   21 FORMAT(1x,A3,i5,i1,5G14.6)
   22 FORMAT(1x,A3,i6,5G14.7)
   23 FORMAT(1x,A3,i12,5G18.11)

      rewind(10) 
      call scanout(sstring,'PRESOLVE',8,10,99)
      write(11,80) sstring    
      call scanout(string,'endendend',9,10,11)

      close (unit=10)
      close (unit=11)
      call byte_close (fbin,ierr)

  999 continue
      stop
 9999 continue
      write(6,*) 'could not find file "indat". abort.'
      stop
      end
c-----------------------------------------------------------------------
      subroutine exitt

      stop 
 
      return
      end
c-----------------------------------------------------------------------
      logical function if_byte_swap_test(bytetest)

      real*4 bytetest,test2
      real*4 test_pattern
      save   test_pattern

      test_pattern = 6.54321
      eps          = 0.00020
      etest        = abs(test_pattern-bytetest)
      if_byte_swap_test = .true.
      if (etest.le.eps) if_byte_swap_test = .false.

      test2 = bytetest
      call byte_reverse(test2,1,ierr)
      write(6,*) 'byte swap:',if_byte_swap_test,bytetest,test2
      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_xyz(buf,igroup,x,y,z,ifbswap,ndim,wdsizi)

      logical ifbswap
      integer buf(0:49),wdsizi
      real*8 x(8),y(8),z(8)   

      nwds = (1 + ndim*(2**ndim))*(wdsizi/4) ! group + 2x4 for 2d, 3x8 for 3d

      if     (ifbswap.and.wdsizi.eq.8) then
          call byte_reverse8(buf,nwds,ierr)
      elseif (ifbswap.and.wdsizi.eq.4) then
          call byte_reverse (buf,nwds,ierr)
      endif


      if(wdsizi.eq.8) then
         call copyi4(igroup,buf(0),1) !0-1
         if (ndim.eq.3) then
            call copy8(x,buf( 2),8) !2 --17
            call copy8(y,buf(18),8) !18--33
            call copy8(z,buf(34),8) !34--49
         else
            call copy8(x,buf( 2),4) !2 --9
            call copy8(y,buf(10),4) !10--17
         endif
      else
         igroup = buf(0)
         if (ndim.eq.3) then
            call copy48(x,buf( 1),8)
            call copy48(y,buf( 9),8)
            call copy48(z,buf(17),8)
         else
            call copy48(x,buf( 1),4)
            call copy48(y,buf( 5),4)
         endif
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine buf_to_curve(buf,nel,ifbswap,wdsizi)

      integer e,f,wdsizi,buf(30)
      logical ifbswap
      character*1 ans
      real*8 r1,r2,r3,r4,r5

      nwds = (2 + 1 + 5)*(wdsizi/4) 

      if     (ifbswap.and.wdsizi.eq.8) then
          call byte_reverse8(buf,nwds-2,ierr)
      elseif (ifbswap.and.wdsizi.eq.4) then
          call byte_reverse (buf,nwds-1,ierr)
      endif

      if(wdsizi.eq.8) then
        call copyi4(ieg,buf(1),1) !1-2

        call copyi4(iedg,buf(3),1) !3-4

        call copy8(r1,buf(5) ,1)
        call copy8(r2,buf(7) ,1)
        call copy8(r3,buf(9) ,1)
        call copy8(r4,buf(11),1)
        call copy8(r5,buf(13),1) 
        call chcopy(ans,buf(15),1)
      else
        ieg   = buf(1)
        iedg  = buf(2)

        call copy48(r1,buf(3))
        call copy48(r2,buf(4))
        call copy48(r3,buf(5))
        call copy48(r4,buf(6))
        call copy48(r5,buf(7))
        call chcopy(ans,buf(8),1)
      endif
      IF (nel.LT.1000) THEN
        write(11,60) IEDG,IEG,R1,R2,R3,R4,R5,ANS
      ELSEIF (nel.LT.1 000 000) then
        write(11,61) IEDG,IEG,R1,R2,R3,R4,R5,ANS
      ELSE
        write(11,62) IEDG,IEG,R1,R2,R3,R4,R5,ANS
      ENDIF
 60   FORMAT(I3,I3 ,5G14.6,1X,A1)
 61   FORMAT(I2,I6 ,5G14.6,1X,A1)
 62   format(i2,i12,5g14.6,1x,a1)

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_bc(cbc,bc,buf,nelm,ifbswap,wdsizi)

      character*3 cbc(6,nelm)
      real*8       bc(5,6,nelm)

      integer e,f,buf(30),wdsizi
      logical ifbswap

      nwds = (2 + 1 + 5)*(wdsizi/4)      
      if     (ifbswap.and.wdsizi.eq.8) then
          call byte_reverse8(buf,nwds-2,ierr)
      elseif (ifbswap.and.wdsizi.eq.4) then
          call byte_reverse (buf,nwds-1,ierr)
      endif

      if(wdsizi.eq.8) then
        call copyi4(e,buf(1),1) !1-2

        call copyi4(f,buf(3),1) !3-4

        call copy8(bc(1,f,e),buf(5),5) !5--14
        call chcopy(cbc( f,e),buf(15),3)!15-16

      else
        e  = buf(1)
        f  = buf(2)

        call copy48 ( bc(1,f,e),buf(3),5)
        call chcopy (cbc(  f,e),buf(8),3)

      endif
      return
      end
c-----------------------------------------------------------------------
