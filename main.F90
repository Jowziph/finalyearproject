program main

  
  use datastructures_mod
  use createhelmholtz_mod
  use applyhelmholtz_mod
  use jacobi_mod

  implicit none

  integer :: nx, ny ! Grid size in x- and y-direction
  real(kind=8) :: alpha, beta ! Model parameters
  real(kind=8) :: omega ! Jacobi relaxation parameter
  integer :: niter

  type(FieldType) :: u ! Numerical solution
  type(FieldType) :: u_exact ! Exact solution
  type(FieldType) :: b ! right hand side

  type(MatrixType) :: H ! Helmholtz matrix

  nx = 8
  ny = 8
  alpha = 1.0
  beta = 1.0
  omega = 1.0
  niter = 8

  ! Print out parameters
  write(*,'("nx    = ",I5)') nx
  write(*,'("ny    = ",I5)') ny
  write(*,'("alpha = ",F8.4)') alpha
  write(*,'("beta  = ",F8.4)') beta
  write(*,'("omega = ",F8.4)') omega
  write(*,'("niter = ",I5)') niter
  

  ! Create field with numerical solution
  u%nx = nx
  u%ny = ny
  u%pr = 1
  allocate(u%data_sgl(nx,ny))

  ! Create field with exact solution
  u_exact%nx = nx
  u_exact%ny = ny
  u_exact%pr = 1
  allocate(u_exact%data_sgl(nx,ny))
  
  ! Create field with RHS
  b%nx = nx
  b%ny = ny
  b%pr = 1
  allocate(b%data_sgl(nx,ny))

  ! Create Helmholtz matrix
  H%nx = nx
  H%ny = ny
  H%pr = 1
  allocate(H%data_sgl_D(nx,ny))
  allocate(H%data_sgl_N(nx,ny))
  allocate(H%data_sgl_S(nx,ny))
  allocate(H%data_sgl_E(nx,ny))
  allocate(H%data_sgl_W(nx,ny))

  ! Populate matrix
  call createH(H,alpha,beta)

  ! Fill u_exact with random values
  call random_number(u_exact%data_sgl(:,:))

  ! Apply matrix to get b = H.u_{exact}
  call applyH(H,u_exact,b)

  ! Initialise Jacobi
  call initialise_jacobi(u)
  
  ! Apply jacobi iteration
  call jacobi(H,b,u,omega,niter)

  ! Initialise Jacobi
  call finalise_jacobi()
  
  ! Deallocate memory
  deallocate(u%data_sgl)
  deallocate(b%data_sgl)
  deallocate(u_exact%data_sgl)
  deallocate(H%data_sgl_D)
  deallocate(H%data_sgl_N)
  deallocate(H%data_sgl_S)
  deallocate(H%data_sgl_E)
  deallocate(H%data_sgl_W)
  
end program main
