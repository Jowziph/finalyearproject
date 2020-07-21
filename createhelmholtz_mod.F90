module createhelmholtz_mod

  use datastructures_mod

  implicit none

contains

  ! The function K(x,y)
  function K(x,y) result(res)
    implicit none
    real(kind=8), intent(in) :: x
    real(kind=8), intent(in) :: y
    real(kind=8) :: res
    real(kind=8) :: pi = 4*atan(1.0_8)
    real(kind=8) :: l1,l2
    l1 = 2.0
    l2 = 2.0
!    res = 1.0 !2.0+sin(l1*pi*x)*sin(l2*pi*y)
    if ((0.2_8 < x) .and. (x < 0.8_8) .and. (0.2_8 < y) .and. (y < 0.8_8)) then
       res = 1.0
    else
       res = 1.E-6
    endif
  end function K

  ! Create entries of the helmholtz matrix. This matrix discretises the
  ! linear operator H such that
  !
  ! H(u) = -alpha*div(K*grad(u)) + beta*u
  !
  ! alpha and beta are positive parameters
  subroutine createH(H,alpha,beta)
    implicit none

    type(MatrixType), intent(inout) :: H ! Helmholtz matrix
    real(kind=8), intent(in) :: alpha, beta ! parameters
    integer :: nx, ny, i, j ! grid size in x- and y-direction
    real(kind=8) :: nxd, nyd, id, jd
    nx = H%nx
    ny = H%ny
    nxd = nx
    nyd = ny

    ! Helmholtz matrix constructed
    if (H%pr == 1) then
       do j = 1,ny
          do i = 1,nx
             H%data_sgl_N(i,j) = -alpha*(nyd/nxd)*K((i-0.5)/nxd,j/nyd)
             H%data_sgl_S(i,j) = -alpha*(nyd/nxd)*K((i-0.5)/nxd,(j-1)/nyd)
             H%data_sgl_E(i,j) = -alpha*(nxd/nyd)*K(i/nxd,(j-0.5)/nyd)
             H%data_sgl_W(i,j) = -alpha*(nxd/nyd)*K((i-1)/nxd,(j-0.5)/nyd)
             H%data_sgl_D(i,j) = -(H%data_sgl_N(i,j)+H%data_sgl_S(i,j)+H%data_sgl_E(i,j)+H%data_sgl_W(i,j))+beta/(nxd*nyd)
          end do
       end do
    else
       do i = 1,nx
          do j = 1,ny
             H%data_dbl_N(i,j) = -alpha*(nyd/nxd)*K((i-0.5)/nxd,j/nyd)
             H%data_dbl_S(i,j) = -alpha*(nyd/nxd)*K((i-0.5)/nxd,(j-1)/nyd)
             H%data_dbl_E(i,j) = -alpha*(nxd/nyd)*K(i/nxd,(j-0.5)/nyd)
             H%data_dbl_W(i,j) = -alpha*(nxd/nyd)*K((i-1)/nxd,(j-0.5)/nyd)
             H%data_dbl_D(i,j) = -(H%data_dbl_N(i,j)+H%data_dbl_S(i,j)+H%data_dbl_E(i,j)+H%data_dbl_W(i,j))+beta/(nxd*nyd)
          end do
       end do
    end if
    ! populate entries of the matrix her
   end subroutine createH

  
end module createhelmholtz_mod
