program test2

  use datastructures_mod
  use createhelmholtz_mod
  use applyhelmholtz_mod
  use jacobi_mod
  
  implicit none

  integer :: nx, ny ! Grid size in x- and y-direction
  integer :: i, j
  real(kind=8) :: nxd, nyd, k1, k2, x, y, l1, l2, z, norm, norm_old
  real(kind=8) :: pi = 4*atan(1.0_8)
  real(kind=8) :: alpha, beta ! Model parameters
  real(kind=8) :: omega ! Jacobi relaxation parameter
  integer :: niter

  type(FieldType) :: u ! Numerical solution
  type(FieldType) :: v
  type(FieldType) :: u_exact ! Exact solution
  type(FieldType) :: b ! right hand side

  type(MatrixType) :: H ! Helmholtz matrix

  nx = 512
  ny = 512
  alpha = 1.0_8
  beta = 1.0_8
  omega = 1.0_8!(1.0_8/3.0_8)
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
  !  call random_number(u_exact%data(:,:))

  v%nx = nx
  v%ny = ny
  v%pr = 1
  allocate(v%data_sgl(nx,ny))

  l1 = 3.0_8
  l2 = 4.0_8
  nxd = nx
  nyd = ny
  !k1 = nxd/(2.0_8)
  !k2 = nxd/(2.0_8)
  k1 = 2.0_8
  k2 = 4.0_8
  
  do i = 1,nx
     x = (i-0.5)/nxd
     do j = 1,ny
        y = (j-0.5)/nyd
        u%data_sgl(i,j) = sin(k1*pi*x)*sin(k2*pi*y)
        v%data_sgl(i,j) = (1/(nxd*nyd))*beta+4*alpha*(sin(k1*pi/(2*nx))*sin(k1*pi/(2*nx)))
        v%data_sgl(i,j) = v%data_sgl(i,j)+4*alpha*sin(k2*pi/(2*ny))*sin(k2*pi/(2*ny))
        v%data_sgl(i,j) = v%data_sgl(i,j)*u_exact%data_sgl(i,j)
     end do
  end do
  
  u_exact%data_sgl(:,:) = 0
  
  ! Apply matrix to get b = H.u_{exact}
!  call applyH(H,u_exact,b)
  b%data_sgl(:,:) = 0

  call initialise_jacobi(u)
  
  z = (beta/(nxd*nyd))+4*alpha*(sin(pi*k1/(2*nxd))*sin(pi*k1/(2*nxd))+sin(pi*k2/(2*nyd))*sin(pi*k2/(2*nyd)))
  !z = (beta)+4*alpha*(sin(pi*k1/(2*nxd))*sin(pi*k1/(2*nxd))+sin(pi*k2/(2*nyd))*sin(pi*k2/(2*nyd)))
  z = 1 - omega* z/((beta/(nxd*nyd))+4*alpha)
  !z = 1 - omega* z/((beta)+4*alpha)
  print*,"true value",z

  ! Apply jacobi iteration
  norm_old = two_norm(u)
  do i =1,niter
!     print*,u%data
     call jacobi(H,b,u,omega,1)
     norm = two_norm(u)
     print*,"norm",i,(norm/norm_old)
     norm_old = norm
  end do
  print*,norm

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
  deallocate(v%data_sgl)
  
end program test2
