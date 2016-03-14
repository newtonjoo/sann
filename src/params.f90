module params
 implicit none
 save

 integer(4) :: ncode
 integer(4) :: mxnaa
 integer(4) :: nvec

 !integer(4), parameter :: nvec=892369

 !integer(4), parameter :: knn=100               ! k-nearest neighbors
 integer(4), parameter :: knn=150               ! k-nearest neighbors
 integer(4), parameter :: mnn=300               ! max nearest neighbors

 integer(4), parameter :: hwin=7
 integer(4), parameter :: win=15
 integer(4), parameter :: vsize=20*win
 !integer(4), parameter :: vsize2=3*win

 integer(4) :: i, j
 private    :: i, j

 real(4)      :: w(vsize)  = (/ ( ( (8.-abs(8-i))**2,j=1,20),i=1,15 ) /)
 real(4)      :: zeros(20) = (/ (0.0, i=1,20) /)
 integer(4)   :: istart(-hwin:hwin)= (/ 1, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201, 221, 241, 261, 281 /)

 real(4)      :: vec(vsize)                        ! 20*win = 300
 !real(4)      :: vec2(vsize2)

 ! best for predict_zs
 !real(4), parameter    :: expo0=8.21
 !integer(4), parameter :: k0=34

 ! mapping vars
 !type :: invmap
 !    integer(4)   :: pid
 !    integer(4)   :: rid
 !end type invmap
 !integer(4)   :: map(mxnaa,ncode)               ! pid,rid -> idx
 !type (invmap):: inv(nvec)                      ! idx -> pid, rid (same with ii,jj)

 ! code variables (size=ncode)
 character(5),allocatable :: codes(:)
 integer(4),allocatable   :: nres(:)        ! sequence length of each code
 integer(4),allocatable   :: mres(:)        ! real length of each code


 ! database variables (size=nvec)
 integer(4),allocatable :: ii(:)           ! i'th proteins in the database
 integer(4),allocatable :: jj(:)           ! j'th residue of i'th proteins in the database
 integer(1),allocatable :: seqn(:)         ! residue number represent it's amino acid (0:22)
 !integer(1):: seqn(nvec)         ! residue number represent it's amino acid (0:22)
 integer(1),allocatable :: three(:)        ! 0: not, 1:coil, 2:helix, 3: strand
 integer(1),allocatable :: four(:)         ! 0: not, 1:coil, 2:helix, 3: strand, 4: turn
 integer(1),allocatable :: eight(:)        ! 0: not, COIL  = [1:'S', 2:'T', 3:' ' or '_']
                                  !         HELIX = [4:'H', 5:'G', 6:'I']
                                  !         STRAND= [7:'E', 8:'B']

 integer(2),allocatable :: acc(:)          ! -1: not, 0~   absolte accessibility
 real   (4),allocatable :: rsa(:)          ! -1: not, 0~   absolte accessibility
 real(4),allocatable    :: p(:,:)         ! feature vectors (vsize,nvec)
 !real   (4) :: rsa(nvec)          ! -1: not, 0~   absolte accessibility
 !real(4)      :: p(vsize,nvec)           ! feature vectors

 integer(1),allocatable :: state3(:)      ! (nvec)
 integer(1),allocatable :: state2(:,:)    ! (4,nvec)

 !integer(1) :: state3(nvec)      ! (nvec)
 !integer(1) :: state2(4,nvec)    ! (4,nvec)

 contains

 subroutine param_set(k0,e0)
 implicit none
 integer(4) :: k0
 real(4)    :: e0

 k0 = 62
 e0 = 6.48

 ncode  = 10021
 mxnaa  = 1733
 nvec   = 2251422

 end subroutine param_set

end module params
