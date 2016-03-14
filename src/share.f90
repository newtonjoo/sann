module share
 implicit none
 save
 character(20) :: amino20 = 'ACDEFGHIKLMNPQRSTVWY'           ! 'X':0
 character(23) :: amino22 = 'ARNDCQEGHILKMFPSTWYVBZ'         ! 'X':0
 character(1)  :: one20(0:20) = (/ 'X','A','C','D','E','F','G','H','I','K','L',   &
                                   'M','N','P','Q','R','S','T','V','W','Y' /)
 character(1)  :: one22(0:22) = (/ 'X','A','R','N','D','C','Q','E','G','H','I',   &
                                       'L','K','M','F','P','S','T','W','Y','V',   &
                                       'B','Z' /)
 character(1)  :: sss3(3)  = (/ 'C', 'H', 'E' /)
 character(1)  :: sss4(4)  = (/ 'C', 'H', 'E', 'T' /)
 character(1)  :: sss8(8)  = (/ 'S', 'T', '_', &
                                'H', 'G', 'I', &
                                'E', 'B' /)
 character(1)  :: saa3(3)  = (/ 'B', 'I', 'E' /)
 character(1)  :: saa2(2)  = (/ 'B', 'E' /)
contains

 subroutine seq2num(naa, seq, seqn)
  implicit none
  integer(4)   :: naa
  character(1) :: seq(naa)
  integer(1)   :: seqn(naa)
  integer(4)   :: i
  do i=1, naa
    seqn(i) = index(amino22,seq(i),BACK=.TRUE.)
  enddo
 end subroutine seq2num

end module share

