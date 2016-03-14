!
! Solvent Accessibility Parameters
!

module sacc
 implicit none
 save

 ! Solvent Accessibility
 !                                  X,  A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,
 ! Chotia                               L,  K,  M,  F,  P,  S,  T,  W,  Y,  V,  B,  Z 
 real(4)      :: mxacc(0:22) = (/ 180,115,225,160,150,135,180,190, 75,195,175, &
                                      170,200,185,210,145,115,140,255,230,155,155,185 /)
 character(1) :: amino(0:22) = (/'X','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z'/) 

 ! threshold from literature (Rost & Sander)
 real(4)      :: thr3(2) = (/ 0.09, 0.36 /)  ! three states (Buried,Intermediate,Exposed)
 real(4)      :: thr2(4) = (/ 0.0, 0.05, 0.16, 0.25 /)      ! 4 kind of two states threshold

 ! variables
 integer(2)   :: sa           ! absolute value
 real(4)      :: ra           ! relative value

 ! states variables
 integer(1)   :: st3, st2(4)

end module sacc
