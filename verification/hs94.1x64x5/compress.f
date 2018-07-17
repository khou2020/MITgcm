      module compression
      implicit none

      contains

      subroutine CompressWr(unitid, R)
        integer :: unitid
        real*8, DIMENSION(:, :, :, :, :) :: R
        real*8, DIMENSION(getelement(shape(R), 1), :, :, :, :) :: T
        real*8, DIMENSION(SIZEOF(R)/8) :: F
        real*4, DIMENSION(SIZEOF(R)/8) :: C
        CALL TIMER_START('ArgCompress',myThid)
        F = reshape(R, shape(F))
        C = REAL(F, 4)
        CALL TIMER_STOP('ArgCompress',myThid)
        CALL TIMER_START('ArgStore',myThid)
        write(unit=unitid) C
        CALL TIMER_STOP('ArgStore',myThid)
        contains 
        function getelement(v, i)
          integer, DIMENSION(:) :: v
          integer :: i
          getelement = v(i)
          return 
        end function getelement
      END subroutine CompressWr

      subroutine CompressRd(unitid, D)
        integer :: unitid
        real*8, DIMENSION(:, :, :, :, :) :: D
        real*4, DIMENSION(SIZEOF(D)/8) :: C
        real*8, DIMENSION(SIZEOF(D)/8) :: F
        CALL TIMER_START('ArgRestore',myThid)
        read(unit=unitid) C
        CALL TIMER_STOP('ArgRestore',myThid)
        CALL TIMER_START('ArgDecompress',myThid)
        F = REAL(C, 8)
        D = reshape(F, shape(D))
        CALL TIMER_STOP('ArgDecompress',myThid)
      END subroutine CompressRd
      end module