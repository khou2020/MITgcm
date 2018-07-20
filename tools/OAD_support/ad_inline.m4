dnl Process this m4 file to produce ]C] language file.
dnl
dnl If you see this line, you can ignore the next one.
C /* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
dnl
include(`foreach.m4')dnl
include(`forloop.m4')dnl
dnl
define(`CONCAT', `$1$2')dnl
dnl
define(`CP_ARG_STORE',changequote(`[', `]')dnl
[dnl
      subroutine cp_arg_store_$2_$3$4($8)
C $OpenAD$ INLINE DECLS
ifelse($4, `', `',`dnl
        use OAD_active
')dnl
        use OAD_cp
        implicit none
        CONCAT(ifelse($7, `', `$5', `type(active)'), ifelse($6, `0', `', ``, dimension(forloop(`i', `1', $6, `ifelse(i, `1', `:', `,:')'))'')) :: $8
C $OpenAD$ END DECLS
#ifdef OAD_DEBUG_CP
        !write(standardmessageunit,*)'OAD: cp write $8$7 ', CONCAT($8, ifelse($6, `0', `', `(forloop(`i', `1', $6, `ifelse(i, `1', `1', `,1')'))'))$7
#endif
C       cp_arg_store_$2_$3$4
        CALL TIMER_START('ArgStore',myThid)
ifelse(`$2', `bool',`dnl
        write(cp_io_unit) $8$7
', `dnl
        CALL ifelse($1, `store', `CompressWr_$2_$6', `CompressRd_$2_$6')(cp_io_unit, $8$7)
')dnl
        CALL TIMER_STOP('ArgStore',myThid)
      end subroutine 

]changequote([`], [']))dnl
dnl
define(`CP_ARG_RESTORE',changequote(`[', `]')dnl
[dnl
      subroutine cp_arg_restore_$2_$3$4($8)
C $OpenAD$ INLINE DECLS
ifelse($4, `', `',`dnl
        use OAD_active
')dnl
        use OAD_cp
        implicit none
        CONCAT(ifelse($7, `', `$5', `type(active)'), ifelse($6, `0', `', ``, dimension(forloop(`i', `1', $6, `ifelse(i, `1', `:', `,:')'))'')) :: $8
C $OpenAD$ END DECLS
C       cp_arg_restore_$2_$3$4
        CALL TIMER_START('ArgRestore',myThid)
ifelse(`$2', `bool',`dnl
        write(cp_io_unit, $8$7)
', `dnl
        CALL ifelse($1, `store', `CompressRd_$2_$6', `CompressRd_$2_$6')(cp_io_unit, $8$7)
')dnl
        CALL TIMER_STOP('ArgRestore',myThid)
#ifdef OAD_DEBUG_CP
        !write(standardmessageunit,*)'OAD: cp read $8$7 ', CONCAT($8, ifelse($6, `0', `', `(forloop(`i', `1', $6, `ifelse(i, `1', `1', `,1')'))'))$7
#endif
      end subroutine 

]changequote([`], [']))dnl
dnl
dnl
dnl
define(`CP_ARG_WR',
`dnl
CP_ARG_STORE(`store', $2, $3, $4, $5, $6, $7, $8)dnl
CP_ARG_RESTORE(`store', $2, $3, $4, $5, $6, $7, $8)dnl
')dnl
dnl
define(`CP_ARG_COMPOND',
`dnl
ifelse(`$2', `real',
`dnl
foreach(`dt', (`($1, $2, $3, `', $5, $6, `', $8)', dnl
               `($1, $2, $3, `_a', $5, $6, `%v', $8)', dnl
dnl               `($1, $2, $3, `_a_d', $5, $6, `%d', $8)', dnl
               ), `CP_ARG_WR(translit(dt, `()'))')dnl

', `dnl
ifelse(`$2', `realad', `CP_ARG_WR($1, `real', $3, `_a_d', $5, $6, `%d', $8)', `CP_ARG_WR($1, $2, $3, `', $5, $6, `', $8)')
')dnl
')dnl
dnl
define(`CP_ARG_DIM',
`dnl
C $2s -----------------------------------------------------
foreach(`dt', (`($1, $2, `scalar', $4, $5, `0', $7, $8)', dnl
               `($1, $2, `vector', $4, $5, `1', $7, $8)', dnl
               `($1, $2, `matrix', $4, $5, `2', $7, $8)', dnl
               `($1, $2, `three_tensor', $4, $5, `3', $7, $8)', dnl
               `($1, $2, `four_tensor', $4, $5, `4', $7, $8)', dnl
               `($1, $2, `five_tensor', $4, $5, `5', $7, $8)', dnl
               `($1, $2, `six_tensor', $4, $5, `6', $7, $8)', dnl
               ), `CP_ARG_COMPOND(translit(dt, `()'))')dnl
')dnl

define(`CP_ARG_TYPE',
`dnl
foreach(`dt', (`($1, `real', $3, $4, `double precision', $6, $7, `x')', dnl
               `($1, `integer', $3, $4, `integer', $6, $7, `i')', dnl
               `($1, `bool', $3, $4, `logical', $6, $7, `b')', dnl
               `($1, `realad', $3, $4, `double precision', $6, $7, `x')', dnl
               ), `CP_ARG_DIM(translit(dt, `()'))')dnl
')dnl

C taping --------------------------------------------


      subroutine push_s0(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x
C $OpenAD$ END DECLS
      if(oad_dt_sz .lt. oad_dt_ptr) call oad_dt_grow()
      oad_dt(oad_dt_ptr)=x; oad_dt_ptr=oad_dt_ptr+1
      end subroutine 

      subroutine pop_s0(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x
C $OpenAD$ END DECLS
      oad_dt_ptr=oad_dt_ptr-1
      x=oad_dt(oad_dt_ptr)
      end subroutine

      subroutine push_s1(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x(:)
C $OpenAD$ END DECLS
      oad_chunk_size=size(x,1)
      if(oad_dt_sz .lt. oad_dt_ptr+oad_chunk_size)
     + call oad_dt_grow()
      oad_dt(oad_dt_ptr:oad_dt_ptr+oad_chunk_size-1)=
     +x
      oad_dt_ptr=oad_dt_ptr+oad_chunk_size
      end subroutine 

      subroutine pop_s1(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x(:)
C $OpenAD$ END DECLS
      oad_chunk_size=size(x,1)
      oad_dt_ptr=oad_dt_ptr-oad_chunk_size
      x=oad_dt(oad_dt_ptr:oad_dt_ptr+oad_chunk_size-1)
      end subroutine

      subroutine push_s2(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x(:,:)
C $OpenAD$ END DECLS
      oad_chunk_size=size(x,1)*size(x,2)
      if(oad_dt_sz .lt. oad_dt_ptr+oad_chunk_size) 
     + call oad_dt_grow()
      oad_dt(oad_dt_ptr:oad_dt_ptr+oad_chunk_size-1)=
     +reshape(x,(/oad_chunk_size/))
      oad_dt_ptr=oad_dt_ptr+oad_chunk_size
      end subroutine 

      subroutine pop_s2(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      double precision :: x(:,:)
C $OpenAD$ END DECLS
      oad_chunk_size=size(x,1)*size(x,2)
      oad_dt_ptr=oad_dt_ptr-oad_chunk_size
        x=reshape(oad_dt(oad_dt_ptr:oad_dt_ptr+oad_chunk_size-1),
     +shape(x))
      end subroutine

      subroutine apush(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      use OAD_active
      implicit none
      type(active) :: x
C $OpenAD$ END DECLS
      if(oad_dt_sz .lt. oad_dt_ptr) call oad_dt_grow()
      oad_dt(oad_dt_ptr)=x%v; oad_dt_ptr=oad_dt_ptr+1
      end subroutine 

      subroutine apop(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      use OAD_active
      implicit none
      type(active) :: x
C $OpenAD$ END DECLS
      oad_dt_ptr=oad_dt_ptr-1
      x%v=oad_dt(oad_dt_ptr)
      end subroutine

      subroutine push_i_s0(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x
C $OpenAD$ END DECLS
      if(oad_it_sz .lt. oad_it_ptr) call oad_it_grow()
      oad_it(oad_it_ptr)=x; oad_it_ptr=oad_it_ptr+1
      end subroutine 

      subroutine pop_i_s0(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x
C $OpenAD$ END DECLS
      oad_it_ptr=oad_it_ptr-1
      x=oad_it(oad_it_ptr)
      end subroutine

      subroutine push_i_s1(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x(:)
C $OpenAD$ END DECLS
      oad_chunk_size=size(x,1)
      if(oad_it_sz .lt. oad_it_ptr+oad_chunk_size) 
     +call oad_it_grow()
      oad_it(oad_it_ptr:oad_it_ptr+oad_chunk_size-1)=
     +x 
      oad_it_ptr=oad_it_ptr+oad_chunk_size
      end subroutine 

      subroutine pop_i_s1(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x(:)
C $OpenAD$ END DECLS
      oad_chunk_size=size(x,1)
      oad_it_ptr=oad_it_ptr-oad_chunk_size
      x=oad_it(oad_it_ptr:oad_it_ptr+oad_chunk_size-1)
      end subroutine

      subroutine push_i_s2(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x(:,:)
C $OpenAD$ END DECLS
      oad_chunk_size=size(x,1)*size(x,2)
      if(oad_it_sz .lt. oad_it_ptr+oad_chunk_size) 
     + call oad_it_grow()
      oad_it(oad_it_ptr:oad_it_ptr+oad_chunk_size-1)=
     +reshape(x,(/oad_chunk_size/))
      oad_it_ptr=oad_it_ptr+oad_chunk_size
      end subroutine 

      subroutine pop_i_s2(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      integer :: x(:,:)
C $OpenAD$ END DECLS
      oad_chunk_size=size(x,1)*size(x,2)
      oad_it_ptr=oad_it_ptr-oad_chunk_size
        x=reshape(oad_it(oad_it_ptr:oad_it_ptr+oad_chunk_size-1),
     +shape(x))
      end subroutine

      subroutine push_b(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      logical :: x
C $OpenAD$ END DECLS
      if(oad_lt_sz .lt. oad_lt_ptr) call oad_lt_grow()
      oad_lt(oad_lt_ptr)=x; oad_lt_ptr=oad_lt_ptr+1
      end subroutine 

      subroutine pop_b(x)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      logical :: x
C $OpenAD$ END DECLS
      oad_lt_ptr=oad_lt_ptr-1
      x=oad_lt(oad_lt_ptr)
      end subroutine

      subroutine push_s(s)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      character*(80) :: s
C $OpenAD$ END DECLS
      if(oad_st_sz .lt. oad_st_ptr) call oad_st_grow()
      oad_st(oad_st_ptr)=s; oad_st_ptr=oad_st_ptr+1
      end subroutine 

      subroutine pop_s(s)
C $OpenAD$ INLINE DECLS
      use OAD_tape
      implicit none
      character*(80) :: s
C $OpenAD$ END DECLS
      oad_st_ptr=oad_st_ptr-1
      s=oad_st(oad_st_ptr)
      end subroutine

C ----------------------- Propagation -----------------------

      subroutine saxpy(a,x,y)
C $OpenAD$ INLINE DECLS
      use OAD_active
      implicit none
      double precision, intent(in) :: a
      type(active), intent(in) :: x
      type(active), intent(inout) :: y
C $OpenAD$ END DECLS
      y%d=y%d+x%d*(a)
      end subroutine

      subroutine zeroderiv(x)
C $OpenAD$ INLINE DECLS
      use OAD_active
      implicit none
      type(active), intent(out) :: x
C $OpenAD$ END DECLS
      x%d=0.0d0
      end subroutine

      subroutine setderiv(y,x)
C $OpenAD$ INLINE DECLS
      use OAD_active
      implicit none
      type(active), intent(out) :: x
      type(active), intent(in) :: y
C $OpenAD$ END DECLS
      x%d=y%d
      end subroutine

      subroutine incderiv(y,x)
C $OpenAD$ INLINE DECLS
      use OAD_active
      implicit none
      type(active), intent(out) :: x
      type(active), intent(in) :: y
C $OpenAD$ END DECLS
      x%d=x%d+y%d
      end subroutine

      subroutine decderiv(y,x)
C $OpenAD$ INLINE DECLS
      use OAD_active
      implicit none
      type(active), intent(out) :: x
      type(active), intent(in) :: y
C $OpenAD$ END DECLS
      x%d = x%d - y%d
      end subroutine decderiv

C Checkpointing stuff ---------------------------------------

CP_ARG_TYPE(1,2,3,4,5,6,7,8)

C strings  -----------------------------------------------------
      subroutine cp_arg_store_string_scalar(s)
C $OpenAD$ INLINE DECLS
      use OAD_cp
      implicit none
      character*(80) :: s
C $OpenAD$ END DECLS 
#ifdef OAD_DEBUG_CP
        !write(standardmessageunit,*)'OAD: cp write s ', s
#endif
C       cp_arg_store_string_scalar
        CALL TIMER_START('ArgStore',myThid)
        write(unit=cp_io_unit) s
        CALL TIMER_STOP('ArgStore',myThid)
      end subroutine 
      
      subroutine cp_arg_restore_string_scalar(s)
C $OpenAD$ INLINE DECLS
      use OAD_cp
      implicit none
      character*(80) :: s
C $OpenAD$ END DECLS
C       cp_arg_restore_real_scalar        
        CALL TIMER_START('ArgRestore',myThid)
        read (unit=cp_io_unit) s
        CALL TIMER_STOP('ArgRestore',myThid)
#ifdef OAD_DEBUG_CP
        !write(standardmessageunit,*)'OAD: cp read s ', s
#endif
      end subroutine 