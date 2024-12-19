module OCTUPLESIZE
      integer*8  etot_new
      save     etot_new


      real,save,allocatable,dimension(:,:,:)   ::  bc2, curve2
      real,save,allocatable,dimension(:,:,:,:) ::  xm2, ym2, zm2

      character(1),save,allocatable,dimension(:,:) :: ccurve2
      character(3),save,allocatable,dimension(:,:) :: cbc2

!
end module OCTUPLESIZE