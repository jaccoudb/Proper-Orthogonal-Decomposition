!*< =================================================================
!*< Singular Value Decomposition (SVD) can be viewed as the extension of the
!*< eigenvalue decomposition for the case of non-square matrices. 
!*< For more detailed: Golub, G.H., Van Loan, C.F.: Matrix Computations. 
!*<                     The Johns Hopkins University Press, Baltimore/London (1993)
!*< In general, this decomposition states, that for any real rectangular
!*< matrix U(NxM), there exists orthogonal matrices, Vleft(NxN) and Vright(MxM)
!*< such that: (Vleft^T)UVright = S.
!*<
!*< This source code uses the dgesdd function from the lapack library
!*<
!*< Author: Bruno R. Jaccoud
!*< Version: 2.0
!*< Date: 2021-03-22
!*<
!*< =================================================================
module var_globais
  !* ... VARIABLES FOR COMPUTATIONAL TIME ...
  real(8) time_total,time_begin,time_end
  !* ... PARAMETERS AND VARIABLES TO POD ...
  integer :: nsamples,nsnaps
  integer :: truncEigenValues, epsilonEnergy
  real(8) :: tol
  real(8),allocatable,dimension(:,:) :: snapshot, matrixU
  real(8),allocatable,dimension(:,:) :: matrixVleft, matrixS, matrixVright
  real(8),allocatable,dimension(:) :: singular_values, eigenvalues, energy
  parameter(tol=1.d-6)
end module var_globais
!*
program pod_svd_program
  use var_globais
  implicit none
  integer :: i, j

  !*  ... READING AND ALLOCATE VARIABLES FOR SNAPSHOT MATRIX  ...
  open(999, file = 'snapshots_3.dat', status = 'old')
  read(999,*) nsamples,nsnaps

  ! Allocate snapshot matrix
  call allocation

  do i = 1, nsamples
    ! Snapshots Matrix POD
    read(999,*) (snapshot(i,j), j = 1, nsnaps)
  enddo
  ! Close snapshot file
  close(999)

  ! Initialize CPU time
  time_total=0
  call cpu_time(time_begin)

  ! Determine SVD results
  call SVD(nsamples,nsnaps,snapshot,matrixVleft,matrixS,matrixVright)

  ! Determine a eigenvalues from singular values of SVD.
  truncEigenValues = 0
  do i = 1, nsamples
    singular_values(i) = matrixS(i,i)
    eigenvalues(i) = (singular_values(i)) ** 2.d0
    if(eigenvalues(i).gt.tol) truncEigenValues = truncEigenValues + 1
  enddo

  ! Determine a relative energy of the snapshots captured by the k first POD basis vectors
  do i = 1, nsamples
    energy(i) = sum(eigenvalues(1:i)) / sum(eigenvalues)
    if((1.d0 - energy(i)).ge.1.0) epsilonEnergy = epsilonEnergy + 1
  enddo
  
  ! Output Solutions
  call output(nsamples,nsnaps,1,'matrixS',matrixS)
  call output(nsamples,nsamples,1,'matrixVleft',matrixVleft)
  call output(nsnaps,nsnaps,1,'matrixVright',matrixVright)
  call output(nsamples,1,1,'singular_values',singular_values)
  call output(nsamples,1,1,'eigenvalues',eigenvalues)
  call output(nsamples,1,1,'energy',energy)

  call cpu_time(time_end)
  time_total = time_end - time_begin 
  write(*,*) 
  write(*,*) 'timetotal: ',time_total

  !*
end program pod_svd_program
!*
!*< =================================================================
!*< ALLOCATE DEPENDENT VARIABLES
!*< =================================================================
subroutine allocation
  use var_globais
    !*
    allocate(snapshot(nsamples,nsnaps), matrixU(nsamples,nsnaps))
    allocate(matrixS(nsamples,nsnaps), matrixVleft(nsamples,nsamples), matrixVright(nsnaps,nsnaps))
    allocate(singular_values(nsamples), eigenvalues(nsamples), energy(nsamples))
    !*
end subroutine allocation
!*
!*< =================================================================
!*< ALLOCATE DEPENDENT VARIABLES
!*< =================================================================
subroutine deallocation
  use var_globais
    !*
    deallocate(snapshot, matrixU)
    deallocate(matrixS, matrixVleft, matrixVright)
    deallocate(singular_values, eigenvalues, energy)
    !*
end subroutine deallocation
!*< =================================================================
!*< SOLVED SVD PROBLEM WITH DGESDD LAPACK SUBROUTINE
!*< =================================================================
subroutine SVD(m,n,a,u,ms,v)
  implicit none
  external dgesdd
  !*------------------------------------------------
  integer :: i, m, n
  real(8) :: a(m,n),u(m,m),v(n,n),s(m),ms(m,n)
  !*
  real(8), dimension(:,:), allocatable :: vt
  real(8), dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: iwork
  integer :: ldu,lda,lwork
  integer :: ldvt,info,mx,mn,lwmax
  character(len=1) :: jobz
  !*------------------------------------------------
  !*
  mx = max(m,n)
  mn = min(m,n)
  !*
  jobz = 'a'
  lda = m
  ldu = m
  ldvt = max(1,n)
  !*
  lwork = -1
  lwmax = 3 * mn * mn + max(mx,4*mn*(mn+1)) + 500
  !*
  allocate(vt(n,n))
  allocate(work(lwmax))
  allocate(iwork(8*mn))
  !*  
  call dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info )
  lwork = min( lwmax, int( work( 1 ) ) )
  !write(*,*) lwork,lwmax,work(1)
  call dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info )
  !*
  v = transpose(vt)
  !*
  ms = 0.d0
  do i = 1, m
    ms(i,i) = s(i)
  enddo
  !*
end subroutine SVD
!*===================================================================
!*  PRINT OUTPUT RESULTS
!*===================================================================
subroutine output(m,n,k,name_variable,x)
implicit none
  !*------------------------------------------------
  integer :: i, j, k
  character(len=*):: name_variable 
  integer :: m, n
  real(8) :: x(m,n)
  !*------------------------------------------------
  !*
  open(unit = k, file = trim(name_variable)//".dat")
  do i = 1, m
    write(k,*) (x(i,j), j = 1, n)
  enddo
  close(k)
  !*
end subroutine output
!*