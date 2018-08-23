




      module cpc
        implicit none
        contains
            
      subroutine CompressRdC_real_0(x)
        implicit none
        integer :: cpsize
        double precision :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_real(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_real_0(x)
        implicit none
        integer :: cpsize
        double precision :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_real(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_real_1(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_real(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_real_1(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_real(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_real_2(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_real(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_real_2(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_real(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_real_3(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_real(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_real_3(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_real(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_real_4(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_real(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_real_4(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_real(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_real_5(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_real(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_real_5(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_real(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_real_6(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:,:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_real(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_real_6(x)
        implicit none
        integer :: cpsize
        double precision, dimension(:,:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_real(x, cpsize) 
      end subroutine 


    
      subroutine CompressRdC_integer_0(x)
        implicit none
        integer :: cpsize
        integer :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_integer(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_integer_0(x)
        implicit none
        integer :: cpsize
        integer :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_integer(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_integer_1(x)
        implicit none
        integer :: cpsize
        integer, dimension(:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_integer(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_integer_1(x)
        implicit none
        integer :: cpsize
        integer, dimension(:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_integer(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_integer_2(x)
        implicit none
        integer :: cpsize
        integer, dimension(:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_integer(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_integer_2(x)
        implicit none
        integer :: cpsize
        integer, dimension(:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_integer(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_integer_3(x)
        implicit none
        integer :: cpsize
        integer, dimension(:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_integer(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_integer_3(x)
        implicit none
        integer :: cpsize
        integer, dimension(:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_integer(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_integer_4(x)
        implicit none
        integer :: cpsize
        integer, dimension(:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_integer(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_integer_4(x)
        implicit none
        integer :: cpsize
        integer, dimension(:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_integer(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_integer_5(x)
        implicit none
        integer :: cpsize
        integer, dimension(:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_integer(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_integer_5(x)
        implicit none
        integer :: cpsize
        integer, dimension(:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_integer(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_integer_6(x)
        implicit none
        integer :: cpsize
        integer, dimension(:,:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_integer(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_integer_6(x)
        implicit none
        integer :: cpsize
        integer, dimension(:,:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_integer(x, cpsize) 
      end subroutine 


    
      subroutine CompressRdC_bool_0(x)
        implicit none
        integer :: cpsize
        logical :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_bool(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_bool_0(x)
        implicit none
        integer :: cpsize
        logical :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_bool(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_bool_1(x)
        implicit none
        integer :: cpsize
        logical, dimension(:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_bool(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_bool_1(x)
        implicit none
        integer :: cpsize
        logical, dimension(:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_bool(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_bool_2(x)
        implicit none
        integer :: cpsize
        logical, dimension(:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_bool(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_bool_2(x)
        implicit none
        integer :: cpsize
        logical, dimension(:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_bool(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_bool_3(x)
        implicit none
        integer :: cpsize
        logical, dimension(:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_bool(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_bool_3(x)
        implicit none
        integer :: cpsize
        logical, dimension(:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_bool(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_bool_4(x)
        implicit none
        integer :: cpsize
        logical, dimension(:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_bool(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_bool_4(x)
        implicit none
        integer :: cpsize
        logical, dimension(:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_bool(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_bool_5(x)
        implicit none
        integer :: cpsize
        logical, dimension(:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_bool(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_bool_5(x)
        implicit none
        integer :: cpsize
        logical, dimension(:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_bool(x, cpsize) 
      end subroutine 


      subroutine CompressRdC_bool_6(x)
        implicit none
        integer :: cpsize
        logical, dimension(:,:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressRd_bool(x, cpsize) 
      end subroutine 

      subroutine CompressWrC_bool_6(x)
        implicit none
        integer :: cpsize
        logical, dimension(:,:,:,:,:,:) :: x
        cpsize = SIZEOF(x)
        CALL CompressWr_bool(x, cpsize) 
      end subroutine 



      end module