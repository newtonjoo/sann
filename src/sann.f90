!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SAPRED OPENMP VERSION
! Solvent accessibility prediction based on nearest neighbor method
! Programmed by Keehyoung Joo at KIAS <newton@kias.re.kr>
! 2010. 02. OpenMP Version 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program sann
 use io
 use share
 use params
 use sacc
 implicit none

 ! input files
 character(200) :: inprof

 character(80) :: bout3, bout21, bout22, bout23, bout24
 character(80) :: fout3, fout21, fout22, fout23, fout24
 character(80) :: brout, bzout
 character(80) :: frout, fpout, fzout, fdout

 !character(1), allocatable :: seq(:)
 !real(4), allocatable      :: prf(:,:)
 !real(4), allocatable      :: q(:,:)
 character(1),allocatable :: seq(:)             ! mxnaa
 integer(1),  allocatable :: iseq(:)            ! seq index
 real(4),allocatable      :: prf(:,:)           ! (20,mxnaa)
 real(4),allocatable      :: q(:,:)             ! (vsize,mxnaa)
 integer(4) :: naa

 ! parameter
 integer(4) :: k0 = 0
 real(4)    :: e0 = 0.0

 ! for secondary structure
 character(1) :: ss
 integer(1)   :: s3,s4,s8

 ! for solvent accessibility                                    
 integer(2)   :: st30(knn), st21(knn), st22(knn), st23(knn), st24(knn)
 real(4),allocatable      :: rsaa(:,:)        ! (knn,mxnaa)
 real(4),allocatable      :: prsa(:)          ! prsa(mxnaa)
 real(4)      :: rtot

 ! residue number represent it's amino acid (0:22)
 integer(1)   :: iaa

 ! distance and indexing
 !real(4)      :: d(nvec)
 real(4), allocatable      :: d(:)
 !real(4)      :: z(nvec)
 real(8)      :: dsum, dsum2
 real(8)      :: dave, dsig
 real(4)      :: dmin, dmax
 !integer(4)   :: ndx(nvec)
 integer(4), allocatable   :: ndx(:)

 ! assessments
 real(8),allocatable      :: p2(:,:,:)      ! p2(2,4,mxnaa)
 real(8),allocatable      :: p3(:,:)        ! three state  p3(3,mxnaa)

 ! varibles
 integer(4)   :: ipdb, jres
 integer(4)   :: i, j, k
 integer(4)   :: idx, jdx

 ! cpu time varibles
 real(4)      :: time_begin, time_end, elapsed_time, total_time

 ! neighbor list
 real(4),allocatable      :: zsco(:,:)    ! zsco(mnn,mxnaa)
 real(4)      :: zz
 real(8)      :: zexp
 real(4),allocatable      :: dist(:,:)    ! dist(mnn,mxnaa) 
                                    !we save dist(1:knn,naa), dist(knn+1,:)=ave_dist, dist(knn+2,:)=std

 ! mpivars
 integer(4)   :: nload,ibegin,iend
 integer(4)   :: ierr, jerr

 ! arguments
 integer(4)     :: nargs, iargc
 character(100) :: string, str, nndb_home
 character(100) :: ck2dir
 character(100) :: code
 !character(100) :: dbhome
 !logical :: out_zz = .true.

 ! omp
 integer  omp_get_thread_num
 integer  :: np=8

 integer istatus, length

 !call omp_set_num_threads(4)

 code  = ''
 ck2dir='.'

 
 call get_environment_variable("NNDB_HOME", nndb_home, length, istatus)

 if (istatus ==0) dbhome = nndb_home

 nargs = iargc()
 select case (nargs)
 case (:1)
     print *, '% sann -np 8 -dbhome /home/newton/database  -i 7odcA'
     goto 1000
 case default
     do i=1, nargs
       call getarg(i, string)
       select case (string)
       case ('-i')
           call getarg(i+1, code)
       case ('-np')
           call getarg(i+1, str)
           read(str,*) np
       case ('-dbhome')
           call getarg(i+1, dbhome)
       end select
     enddo
 end select

 call param_set(k0,e0)
 call io_set()

 ! allocation
 allocate(seq(mxnaa))             ! mxnaa
 allocate(prf(20,mxnaa))          ! (20,mxnaa)
 allocate(q(vsize,mxnaa))         ! (vsize,mxnaa)

 allocate(rsaa(knn,mxnaa))        ! (knn,mxnaa)
 allocate(prsa(mxnaa))            ! prsa(mxnaa)

 allocate(p2(2,4,mxnaa))        ! p2(2,4,mxnaa)
 allocate(p3(3,mxnaa))          ! three state  p3(3,mxnaa)

 allocate(zsco(mnn,mxnaa))      ! zsco(mnn,mxnaa)
 allocate(dist(mnn,mxnaa))      ! dist(mnn,mxnaa) 

 call omp_set_num_threads(np)

 print *, 'Protein Solvent Accessibility Prediction Using Nearest Neighbor Method'
 print *, 'Keehyoung Joo at KIAS <newton@kias.re.kr>'
 print *, 'sapred with parameter (k,e)=',k0,e0
 print *, 'dbhome= ', trim(dbhome)
 print *, ''
 print *, 'target code= ', trim(code)

 !stop

 ! open input profiles
 inprof = trim(ck2dir)//'/'//trim(code)//'.ck2'
 open(in_prof,file=inprof,status='old')
 read(in_prof,*) naa
 !allocate (seq(naa))
 !allocate (prf(20,naa))
 !allocate (q(vsize,naa))
 read(in_prof,'(4000A1)') (seq(i),i=1,naa)
 print *, 'naa= ',naa
 print *, 'seq= ',seq(1:naa)
 do i=1, naa
   read(in_prof,*) (prf(k,i),k=1,20)
 enddo
 close(in_prof)
 
 ! mapping seq to iseq (integer number)
 allocate(iseq(naa))
 do i=1, naa
    do j=0,22
       if (amino(j).eq.seq(i)) iseq(i) = j
    enddo
 enddo

 call make_feature(naa,prf(:,1:naa),q(:,1:naa))

 !print *, q(:,1:naa)

 ! print weight
 !do i= 1, 281, 20
 !  write(6,'(20F4.0,1x)') (w(i:i+19))
 !enddo

 if (code .eq. '') then
     print *, '% sapred -i 7odcA'
     goto 1000
 endif

 print *, ''
 print *, 'number of processors used for calculation=',np
 !$omp parallel
 print *, "I am a processor ", omp_get_thread_num()
 !$omp end parallel

 if (a3only) then
    fout3 = trim(code)//'.sa2'
 else
    fout3 = trim(code)//'.a3'
 endif

 fout21= trim(code)//'.a21'
 fout22= trim(code)//'.a22'
 fout23= trim(code)//'.a23'
 fout24= trim(code)//'.a24'
 frout = trim(code)//'.rsa'
 fpout = trim(code)//'.prsa'
 fzout = trim(code)//'.zs'
 fdout = trim(code)//'.d'
 !print *, 'result 3   =', trim(fout3)
 !print *, 'result 2 1 =', trim(fout21)
 !print *, 'result 2 2 =', trim(fout22)
 !print *, 'result 2 3 =', trim(fout23)
 !print *, 'result 2 4 =', trim(fout24)
 !print *, 'result rsa =', trim(frout)
 !print *, 'result zz  =', trim(fzout)

 ! variable allocation
 allocate(p(vsize,nvec))
 !allocate(three(nvec))

 allocate(seqn(nvec))
 allocate(rsa(nvec))
 allocate(state3(nvec))
 allocate(state2(4,nvec))

 ! setup length of each data files

 ! open database

 ! read all profiles and dssp data from database
 call cpu_time( time_begin )
 inquire (IOLENGTH=length_vec_db) vec
 inquire (IOLENGTH=length_rsa_db) iaa, sa, ra, st3, st2(1:4)

 open(io_vec_db,file=vecfile,form='unformatted',action='read',status='old',access='direct',recl=length_vec_db)
 open(io_rsa_db,file=rsafile,form='unformatted',action='read',status='old',access='direct',recl=length_rsa_db)

 do idx=1, nvec
   read(io_vec_db, rec=idx) p(1:vsize,idx)
 enddo

 do idx=1, nvec
   read(io_rsa_db, rec=idx) seqn(idx), sa, rsa(idx), state3(idx), state2(1:4,idx)
 enddo

 call cpu_time( time_end )
 close(io_vec_db)
 close(io_rsa_db)
 print *, 'reading database is done...', time_end-time_begin, ' seconds'


 ! load balancing
 !call get_load(naa,nodes,me,nload,ibegin,iend)
 !write(*,'(A,I4,A,I,1x,I,A,I,1x,I)') 'me=', me,' my task range = (',ibegin, iend,') among ',nload,naa

!inquire (IOLENGTH=length_a3 ) st3,    p3(1:3), zz
!inquire (IOLENGTH=length_a21) st2(1), p2(1:2,1)
!inquire (IOLENGTH=length_a22) st2(2), p2(1:2,2)
!inquire (IOLENGTH=length_a23) st2(3), p2(1:2,3)
!inquire (IOLENGTH=length_a24) st2(4), p2(1:2,4)
!inquire (IOLENGTH=length_rsa) rsaa(1:knn)
!inquire (IOLENGTH=length_zz ) zsco(1:knn)
!open(iw_a3, file=bout3, access='direct',form='unformatted',status='unknown',recl=length_a3 )
!open(iw_a21,file=bout21,access='direct',form='unformatted',status='unknown',recl=length_a21)
!open(iw_a22,file=bout22,access='direct',form='unformatted',status='unknown',recl=length_a22)
!open(iw_a23,file=bout23,access='direct',form='unformatted',status='unknown',recl=length_a23)
!open(iw_a24,file=bout24,access='direct',form='unformatted',status='unknown',recl=length_a24)
!open(iw_rsa,file=brout, access='direct',form='unformatted',status='unknown',recl=length_rsa)

!if (out_zz) then
!  open(iw_zz,file=bzout,access='direct',form='unformatted',status='unknown',recl=length_zz)
!endif

 ! prediction
 !allocate(p3(3,naa))

 ! get neighbors
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 p2(:,:,1:naa) = 0.0d0
 p3(:,naa)     = 0.0d0

 call cpu_time( time_begin )
 !$omp parallel do private(idx,jdx,k,vec,d,ndx,st30,st21,st22,st23,st24,dsum,dsum2,dave,dsig,zexp) &
 !$omp shared(naa,w,p,q,state3,state2,k0,e0,p3,p2,rsa,rsaa,prsa,zsco,dist)
 do idx=1, naa
  
  if (.not. allocated(d)) allocate(d(nvec))
  if (.not. allocated(ndx)) allocate(ndx(nvec))

  do jdx=1, nvec
    d(jdx) = sum(w(:) * abs(q(:,idx) - p(:,jdx)))
  enddo
  
  !print *, 'idx=',idx, d(1), d(2),omp_get_thread_num()

  call qsortd(nvec,d,ndx)

  st30(1:knn) = state3(ndx(1:knn))
  st21(1:knn) = state2(1,ndx(1:knn))
  st22(1:knn) = state2(2,ndx(1:knn))
  st23(1:knn) = state2(3,ndx(1:knn))
  st24(1:knn) = state2(4,ndx(1:knn))

  ! statistics except ndx(1)
  !dsum = 0.d0
  !dsum2= 0.d0
  !do k=2, nvec
  !do k=1, nvec
  !  dsum  = dsum  + d(ndx(k))
  !  dsum2 = dsum2 + d(ndx(k))**2
  !enddo
  dsum = sum(d(1:nvec))
  dsum2 = dot_product(d(1:nvec),d(1:nvec))

  !print *, idx,'dsum=',dsum
  !print *, idx,'dsum2=',dsum2

  ! nvec-1 ???
  dave = dsum/dble(nvec)
  dsig = dsqrt( dsum2/dble(nvec) - dave**2 )

  !print *, idx,'dav=',dave
  !print *, idx,'dsig=',dsig

  !forall (k=1:k0)
  !forall (k=1:knn)
  zsco(1:knn,idx) = real( ( dave - dble(d(ndx(1:knn))) ) / dsig , 4)       ! -Z
  !end forall

  dist(1:knn,idx) = d(ndx(1:knn))
  dist(knn+1,idx) = dave
  dist(knn+2,idx) = dsig

  rsaa(1:knn,idx) = rsa(ndx(1:knn))

  !p2 = 0.0d0
  !p3 = 0.0d0
  prsa(idx) = 0.0 ; rtot = 0.0      ! for continuous prediction
  do k=1, k0
    zexp = zsco(k,idx) ** e0
    p3(st30(k),idx) = p3(st30(k),idx) + zexp
    p2(st21(k),1,idx) = p2(st21(k),1,idx) + zexp
    p2(st22(k),2,idx) = p2(st22(k),2,idx) + zexp
    p2(st23(k),3,idx) = p2(st23(k),3,idx) + zexp
    p2(st24(k),4,idx) = p2(st24(k),4,idx) + zexp
    prsa(idx) = prsa(idx) + zexp * rsaa(k,idx)
    rtot = rtot + zexp
  enddo
  prsa(idx) = prsa(idx) / rtot

  p3(:,idx) = p3(:,idx)/sum(p3(:,idx))
  forall (k=1:4)
    p2(:,k,idx) = p2(:,k,idx) / sum(p2(:,k,idx))
  end forall

  !deallocate(d)
  !deallocate(ndx)
 enddo
 !$omp end parallel do

 !do i=1,naa
 !  write(*, '(I4,1x,A1,1x,3(F6.3,1x),F6.3,2x,F6.3)') i, seq(i), p3(1:3,i), zsco(1,i), rsaa(1,i)
 !enddo

 call cpu_time( time_end )
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 elapsed_time = time_end - time_begin
 print *, 'Elapsed_time...', elapsed_time, 'seconds'
 print *, 'Done... Good luck! ^^'

 !call flush(iw_a3 )
 !call flush(iw_a21)
 !call flush(iw_a22)
 !call flush(iw_a23)
 !call flush(iw_a24)
 !call flush(iw_rsa)
 !if (out_zz) call flush(iw_zz)

 if (out_a3) then
    open(iw_a3f,file=fout3)
    write(iw_a3f, '(A,1x,I3,1x,F4.2,2x,I2)') '# SANN VFORMAT (SANN V1.0 by K. Joo)'
    write(iw_a3f,*)
 endif
 if (out_a21) then
    open(iw_a21f,file=fout21)
    write(iw_a21f,'(A,1x,I3,1x,F4.2,2x,I2)') '# SANN VFORMAT (SANN V1.0 by K. Joo)'
    write(iw_a21f,*)
 endif
 if (out_a22) then
    open(iw_a22f,file=fout22)
    write(iw_a22f,'(A,1x,I3,1x,F4.2,2x,I2)') '# SANN VFORMAT (SANN V1.0 by K. Joo)'
    write(iw_a22f,*)
 endif
 if (out_a23) then
    open(iw_a23f,file=fout23)
    write(iw_a23f,'(A,1x,I3,1x,F4.2,2x,I2)') '# SANN VFORMAT (SANN V1.0 by K. Joo)'
    write(iw_a23f,*)
 endif
 if (out_a24) then
    open(iw_a24f,file=fout24)
    write(iw_a24f,'(A,1x,I3,1x,F4.2,2x,I2)') '# SANN VFORMAT (SANN V1.0 by K. Joo)'
    write(iw_a24f,*)
 endif
 if (out_rsa) then
    open(iw_rsaf,file=frout)
    write(iw_rsaf,'(A,1x,I3,1x,F4.2,2x,I2,2x,I3)') '# SANN Nearest Relative Accessibilities (SANN V2.0 by K. Joo)'
    write(iw_rsaf,*)
 endif
 if (out_prsa) then
    open(iw_prsaf,file=fpout)
    write(iw_prsaf,'(A,1x,I3,1x,F4.2,2x,I2)') '# SANN Relative Accessibility (SANN V2.0 by K. Joo)'
    write(iw_prsaf,*)
 endif
 if (out_zz) then
    open(iw_zzf,file=fzout)
    write(iw_zzf,'(A,1x,I3,1x,F4.2,2x,I2)') '# Zscore from SANN (by K. Joo)'
    write(iw_zzf,*)
 endif
 if (out_dist) then
    open(iw_dist,file=fdout)
    write(iw_dist,'(A,1x,I3,1x,F4.2,2x,I2)') '# Dists from SANN (by K. Joo)'
    write(iw_dist,*)
 endif

 do idx=1, naa
   st3 = int(maxloc(p3(:,idx),1))
   !write(iw_a3f, 300) idx, seq(idx), saa3(st3), p3(1:3,idx), zsco(1,idx)
   if (out_a3) write(iw_a3f, 310) idx, seq(idx), saa3(st3), p3(1:3,idx), prsa(idx), prsa(idx)*mxacc(iseq(idx)), zsco(1,idx)

   forall (k=1:4)
     st2(k) = int(maxloc(p2(:,k,idx),1))
   end forall

   if (out_a21) write(iw_a21f,400) idx, seq(idx), saa2(st2(1)), p2(1:2,1,idx), zsco(1,idx)
   if (out_a21) write(iw_a22f,400) idx, seq(idx), saa2(st2(2)), p2(1:2,2,idx), zsco(1,idx)
   if (out_a21) write(iw_a23f,400) idx, seq(idx), saa2(st2(3)), p2(1:2,3,idx), zsco(1,idx)
   if (out_a21) write(iw_a24f,400) idx, seq(idx), saa2(st2(4)), p2(1:2,4,idx), zsco(1,idx)
   if (out_rsa) write(iw_rsaf,500) idx, seq(idx), saa3(st3), rsaa(1:knn,idx)
   if (out_zz ) write(iw_zzf, 600) idx, zsco(1:knn,idx)
   if (out_dist)write(iw_dist,650) idx, dist(knn+1,idx), dist(knn+2,idx), dist(1:knn,idx) ! index, ave, sig, dists

   if (prsa(idx) < 0.09) then
       st3 = 1
   elseif (prsa(idx) < 0.36) then
       st3 = 2
   else
       st3 = 3
   endif
   if (out_prsa)write(iw_prsaf,700) idx, seq(idx), saa3(st3), prsa(idx), prsa(idx)*mxacc(iseq(idx)), zsco(1,idx)

 enddo
 if (out_a3 ) close(iw_a3f )
 if (out_a21) close(iw_a21f)
 if (out_a21) close(iw_a22f)
 if (out_a21) close(iw_a23f)
 if (out_a21) close(iw_a24f)
 if (out_rsa) close(iw_rsaf)
 if (out_zz ) close(iw_zzf)
 if (out_dist)close(iw_dist)

 !close(iw_a3 )
 !close(iw_a21)
 !close(iw_a22)
 !close(iw_a23)
 !close(iw_a24)
 !close(iw_rsa)
 !close(iw_zz)

 !call system('rm -rf '//bout3)
 !call system('rm -rf '//bout21)
 !call system('rm -rf '//bout22)
 !call system('rm -rf '//bout23)
 !call system('rm -rf '//bout24)
 !call system('rm -rf '//brout)
 !if (out_zz) call system('rm -rf '//bzout)

 1000 continue

 100 format(I4, I3, I4, 1x, 3F7.3)
 300 format(I4,1x,A1,1x,A1,1x,3F7.3, 4x, F6.2)
 310 format(I4,1x,A1,1x,A1,1x,3F7.3,2x,F8.4,2x,F8.4,4x, F6.2)
 400 format(I4,1x,A1,1x,A1,1x,2F7.3, 4x, F6.2)
 500 format(I4,1x,A1,1x,A1,1x,300F6.3)
 600 format(I4,1x,300F7.3)
 650 format(I4,1x,F8.3,1x,F10.5,1x,300F8.3)
 700 format(I4,1x,A1,1x,A1,1x,F8.4,2x,F8.4,4x,F6.2)

end program sann

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_feature(naa,prf,f)
use params
implicit none
integer(4)  :: naa
real(4)     :: prf(20,naa)
real(4)     :: f(vsize,naa)
integer(4)  :: i,j,k, is,ie

do i=1, naa
 do k=-hwin, hwin
   j = i + k

   is = istart(k)
   ie = is+19

   if (j < 1) then
       f(is:ie,i) = zeros
   else if ( j > naa) then
       f(is:ie,i) = zeros
   else
       f(is:ie,i) = prf(1:20,j)
   endif

 enddo
enddo

end subroutine make_feature

