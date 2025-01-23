c-----------------------------------------------------------------------
      subroutine scanout(string,input,len,infile,outfile)
c
c     scan through infile until input is found and output
c     all lines save that containing input to "outfile"
c
      character*80 string
      character*1 input(1)
      integer infile,outfile
c
      character*1 string1(80)
c
      do line=1,100000
         call blank(string,80)
         read (infile ,80,end=100,err=100) string
c        write(6,*) string
         if (indx1(string,input,len).ne.0) return
         call chcopy(string1,string,80)
         lout = ltrunc(string1,80)
         write (outfile,81) (string1(j),j=1,lout)
      enddo
   80 format(a80)
   81 format(80a1)
c
 100  continue
      return
      end

c-----------------------------------------------------------------------
      subroutine scanoutafter(string,input,len,infile,outfile)
c
c     scan through infile until input is found and output
c     all lines after save that containing input to "outfile"
c
      character*80 string
      character*1 input(1)
      integer infile,outfile
c
      character*1 string1(80)
c
      do line=1,100000
         call blank(string,80)
         read (infile ,80,end=100,err=100) string
         if (indx1(string,input,len).ne.0) return
         call chcopy(string1,string,80)
         lout = ltrunc(string1,80)
         write (outfile,81) (string1(j),j=1,lout)
      enddo
   80 format(a80)
   81 format(80a1)
c
 100  continue
      return
      end

c-----------------------------------------------------------------------
      subroutine copy8(a,b,n)
      real*8 a(1)
      real*8 b(1)
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1), b(1)
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine copy48(a,b,n)
      real*8 a(1)
      real*4 b(1)
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      integer a(1), b(1)
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine chcopy(a,b,n)
      character*1 a(1), b(1)
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine copyi4(a,b,n)
      integer a(1)
      real*8  b(1)
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine blank(s,n)
      character*1 s(1)
      do i=1,n
        s(i)=' '
      enddo
      return
      end
c-----------------------------------------------------------------------
c     comes from n2to3
      function ltrunc(s,n)
      character*1 s(1)
      ltrunc = 0
      do j=n,1,-1
         if (s(j).ne.' ') then
            ltrunc = j 
            return
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
c     comes from n2to3
      INTEGER FUNCTION INDX1(S1,S2,L2)
      CHARACTER*80 S1,S2
C
      N1=80-L2+1
      INDX1=0
      IF (N1.LT.1) RETURN
C
      DO 300 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX1=I
            RETURN
         ENDIF
300   CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
c     comes from n2to3
      subroutine rzero(x,n)
      real x(1)
      do i=1,n
         x(i) = 0.0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine izero(x,n)
      integer x(1)
      do i=1,n
         x(i) = 0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine lineout(io,s,n)
      character*1 s(n)
      m = ltrunc(s,n)
      write(io,1) (s(k),k=1,m)
    1 format(132a1)
      return
      end
c-----------------------------------------------------------------------
      subroutine fbyte_write(s,n)
      character*4 s(n)

      common /byte_flag/ ifbyte
      logical ifbyte

      if (ifbyte) then
         call byte_write(s,n)
      else
         call lineout(6,s,n)
      endif

      return
      end
c-----------------------------------------------------------------------
      logical function strcmp(str1in,str2in)

      character*80 str1in
      character*80 str1
      character*1  vec1(80)
      equivalence (vec1,str1)

      character*80 str2in
      character*80 str2
      character*1  vec2(80)
      equivalence (vec2,str2)

      strcmp = .true.

      call chcopy(str1,str1in,80)
      call chcopy(str2,str2in,80)

      n1 = ltrunc(str1,80)
      n2 = ltrunc(str2,80)
      if (n1.ne.n2) then
         strcmp = .false.
         return
      endif

      do i=1,n1
         if (vec1(i).ne.vec2(i)) then
            strcmp = .false.
            return
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine scan_all(num,file,input,len)
c
c     open file, scan and count how many input is found
c     close file and return the number
c
      character*1 input(1)
      integer num,infile
c
      character*80 file, string
c
      num = 0
      infile = 20
      open(unit=infile, file=file, err=100)
      do line=1,100000
         read (infile,80,end=100,err=101) string
         if (indx1(string,input,len).ne.0) then
            num = num + 1
         endif
      enddo
   80 format(a80)

 100  continue
      close (unit=infile)
      return

 101  continue
      write(6,*)'Error in scan_all',file

      return
      end
c-----------------------------------------------------------------------
