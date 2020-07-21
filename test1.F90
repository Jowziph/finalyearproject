program test1

  use datastructures_mod
  use createhelmholtz_mod
  use applyhelmholtz_mod
  use jacobi_mod
  
  implicit none

  integer :: nx, ny ! Grid size in x- and y-direction
  integer :: i, j
  real(kind=8) :: nxd, nyd, k1, k2, x, y, l1, l2
  real(kind=8) :: pi = 4*atan(1.0_8)
  real(kind=8) :: alpha, beta ! Model parameters
  real(kind=8) :: omega ! Jacobi relaxation parameter
  integer :: niter

  type(FieldType) :: u ! Numerical solution
  type(FieldType) :: v,c
  type(FieldType) :: u_exact ! Exact solution
  type(FieldType) :: b ! right hand side

  type(MatrixType) :: H ! Helmholtz matrix

  nx = 32
  ny = 32
  alpha = 1.0
  beta = 1.0
  omega = 1.0
  niter = 1

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
  !  call random_number(u_exact%data(:,:))

  v%nx = nx
  v%ny = ny
  v%pr = 1
  allocate(v%data_sgl(nx,ny))

  k1 = 2.0_8
  k2 = 2.0_8
  l1 = 0.0_8
  l2 = 2.0_8
  nxd = nx
  nyd = ny
  
  do i = 1,nx
     x = (i-0.5)/nxd
     do j = 1,ny
        y = (j-0.5)/nyd
        u_exact%data_sgl(i,j) = cos(k1*pi*x)*cos(k2*pi*y)
        v%data_sgl(i,j) = u_exact%data_sgl(i,j)*(k1*k1+k2*k2)*K(x,y)
        v%data_sgl(i,j) = v%data_sgl(i,j) + k1*l1*sin(k1*pi*x)*cos(k2*pi*y)*cos(l1*pi*x)*sin(l2*pi*y)
        v%data_sgl(i,j) = v%data_sgl(i,j) + k2*l2*cos(k1*pi*x)*sin(k2*pi*y)*sin(l1*pi*x)*cos(l2*pi*y)
        v%data_sgl(i,j) = alpha * pi*pi * v%data_sgl(i,j) + beta * u_exact%data_sgl(i,j)
     end do
  end do

  ! Apply matrix to get b = H.u_{exact}
  call applyH(H,u_exact,b)
  b%data_sgl(:,:) = b%data_sgl(:,:)

  !  print*,u_exact%data
  c%nx = nx
  c%ny = ny
  c%pr = 1
  allocate(c%data_sgl(nx,ny))
  c%data_sgl = v%data_sgl - b%data_sgl
  print*, two_norm(c)

  call initialise_jacobi(u)
  
  ! Apply jacobi iteration
  call jacobi(H,b,u,omega,niter)

  call finalise_jacobi()
  
  ! Deallocate memory
  deallocate(u%data_sgl)
  deallocate(b%data_sgl)
  deallocate(c%data_sgl)
  deallocate(u_exact%data_sgl)
  deallocate(H%data_sgl_D)
  deallocate(H%data_sgl_N)
  deallocate(H%data_sgl_S)
  deallocate(H%data_sgl_E)
  deallocate(H%data_sgl_W)
  deallocate(v%data_sgl)
  
end program test1
