module io
 implicit none
 save

 ! io directories
 character(100) :: dbhome = '/home/newton/database/casp8/nndb/40/db'

 ! databases (reference)
 character(100) :: vecfile
 character(100) :: rsafile
 character(100) :: inffile
 character(100) :: tabfile

 ! io length varibles for databases
 integer(4)   :: length_vec_db
 integer(4)   :: length_rsa_db
 integer(4)   :: length_inf_db

 ! io length varibles for neighbors
 integer(4)   :: length_indx
 integer(4)   :: length_dist
 integer(4)   :: length_zsco
 integer(4)   :: length_stat

 integer(4)   :: length_indx2
 integer(4)   :: length_dist2
 integer(4)   :: length_zsco2
 integer(4)   :: length_stat2

 ! io length variables for prediction
 integer(4)   :: length_s3
 integer(4)   :: length_s4
 integer(4)   :: length_s8
 integer(4)   :: length_a3
 integer(4)   :: length_a21
 integer(4)   :: length_a22
 integer(4)   :: length_a23
 integer(4)   :: length_a24
 integer(4)   :: length_rsa
 integer(4)   :: length_zz

 ! io units for databases
 integer(4), parameter :: io_vec_db =  10       ! for database.vec
 integer(4), parameter :: io_rsa_db =  12
 integer(4), parameter :: io_inf_db =  13
 integer(4), parameter :: io_dih_db =  14
 integer(4), parameter :: io_tab_db =  15

 ! io units for 2nd layer db
 integer(4), parameter :: io_v2c_db =  16       ! for database.vec


 ! results io
 integer(4), parameter :: iw_sec3 = 31          ! 3 state 
 integer(4), parameter :: iw_sec4 = 32          ! 4 state
 integer(4), parameter :: iw_sec8 = 33          ! 8 state
 integer(4), parameter :: iw_seav = 34          ! for average
 integer(4), parameter :: iw_sall = 35

 integer(4), parameter :: iw_st30 = 40
 integer(4), parameter :: iw_st21 = 41
 integer(4), parameter :: iw_st22 = 42
 integer(4), parameter :: iw_st23 = 43
 integer(4), parameter :: iw_st24 = 44
 integer(4), parameter :: iw_stav = 45

 integer(4), parameter :: io_log  = 55

 ! for predict
 integer(4), parameter :: in_prof = 60
 integer(4), parameter :: ss_prof = 61
 integer(4), parameter :: iw_s3   = 62
 integer(4), parameter :: iw_s4   = 63
 integer(4), parameter :: iw_s8   = 64
 integer(4), parameter :: iw_s3f  = 65
 integer(4), parameter :: iw_s4f  = 66
 integer(4), parameter :: iw_s8f  = 67
 integer(4), parameter :: iw_zz   = 68
 integer(4), parameter :: iw_zzf  = 69
 integer(4), parameter :: iw_dist = 50
 integer(4), parameter :: iw_cout = 51

 ! for sapred
 integer(4), parameter :: iw_a3   = 70
 integer(4), parameter :: iw_a21  = 71
 integer(4), parameter :: iw_a22  = 72
 integer(4), parameter :: iw_a23  = 73
 integer(4), parameter :: iw_a24  = 74
 integer(4), parameter :: iw_rsa  = 75

 integer(4), parameter :: iw_a3f  = 76
 integer(4), parameter :: iw_a21f = 77
 integer(4), parameter :: iw_a22f = 78
 integer(4), parameter :: iw_a23f = 79
 integer(4), parameter :: iw_a24f = 80
 integer(4), parameter :: iw_rsaf = 81
 integer(4), parameter :: iw_prsaf= 82

 ! for output
 logical :: out_a3   = .true.
 logical :: out_a21  = .true.
 logical :: out_a22  = .true.
 logical :: out_a23  = .true.
 logical :: out_a24  = .true.
 logical :: out_rsa  = .true.
 logical :: out_prsa = .true.
 logical :: out_zz   = .true.
 logical :: out_dist = .true.

 logical :: a3only = .false.

!interface get_length
!    module procedure get_length_database
!    module procedure get_length_rsa_database
!    module procedure get_length_int1
!    module procedure get_length_int2
!    module procedure get_length_int4
!    module procedure get_length_real4
!end interface

contains
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine outoff()
    out_a3   = .false.
    out_a21  = .false.
    out_a22  = .false.
    out_a23  = .false.
    out_a24  = .false.
    out_rsa  = .false.
    out_prsa = .false.
    out_zz   = .false.
    out_dist = .false.
 end subroutine outoff

 ! open database
 subroutine open_vec_database()
  implicit none
  open(io_vec_db,file=vecfile,form='unformatted',status='old',access='direct',recl=length_vec_db)
 end subroutine open_vec_database

 subroutine open_rsa_database()
  implicit none
  open(io_rsa_db,file=rsafile,form='unformatted',access='direct',recl=length_rsa_db)
 end subroutine open_rsa_database

 subroutine io_set()
 implicit none
 !character(100) :: dbhome

 !dbhome =trim(home)//'db/'

 ! databases (reference)
 vecfile=trim(dbhome)//'/database.vec'        ! unf direct
 rsafile=trim(dbhome)//'/database.rsa'        ! unf direct
 inffile=trim(dbhome)//'/database.inf'        ! unf direct
 tabfile=trim(dbhome)//'/database.tab'        ! text file

 end subroutine io_set

end module io
