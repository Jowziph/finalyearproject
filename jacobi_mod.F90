module jacobi_mod
  
  use datastructures_mod
  use applyhelmholtz_mod
  
  implicit none

  type(FieldType) :: Hu

contains

  subroutine initialise_jacobi(u)
    implicit none
    type(FieldType), intent(in) :: u

    Hu%nx = u%nx
    Hu%ny = u%ny
    Hu%pr = u%pr

    if (Hu%pr == 1) then
       allocate(Hu%data_sgl(Hu%nx,Hu%ny))
    else
       allocate(Hu%data_dbl(Hu%nx,Hu%ny))
    end if
  end subroutine initialise_jacobi

  subroutine finalise_jacobi()
    implicit none
    if (allocated(Hu%data_sgl)) then
       deallocate(Hu%data_sgl)
    endif
    if (allocated(Hu%data_dbl)) then
       deallocate(Hu%data_dbl)
    endif
  end subroutine finalise_jacobi
  
  ! carry out niter iterations of the weighted Jacobi method:
  ! u <- u + omega*H_D^{-1}(b-H.u)
  ! This approximately solves the equations H.u = b
  subroutine jacobi(H,b,u,omega,niter)
    implicit none

    type(MatrixType), intent(in) :: H ! Helmholtz matrix
    type(FieldType), intent(in) :: b ! right hand side
    type(FieldType), intent(inout) :: u ! solution
    real(kind=8), intent(in) :: omega ! relaxation factor
    integer, intent(in) :: niter ! number of iterations
    integer :: nx, ny ! grid size in x- and y-direction
    integer :: k, i, j ! loop index
    nx = H%nx
    ny = H%ny

    ! Jacobi iteration
    
    if (H%pr == 1) then
       do k=1, niter
          call applyH(H,u,Hu)
          do j = 1,ny
             do i = 1,nx
                u%data_sgl(i,j) = u%data_sgl(i,j) + omega*(1./H%data_sgl_D(i,j))*(b%data_sgl(i,j)-Hu%data_sgl(i,j))
             end do
          end do
       end do
    else
       do k=1, niter
          call applyH(H,u,Hu)
          do j = 1,ny
             do i = 1,nx
                u%data_dbl(i,j) = u%data_dbl(i,j) + omega*(1./H%data_dbl_D(i,j))*(b%data_dbl(i,j)-Hu%data_dbl(i,j))
             end do
          end do
       end do
    end if
  end subroutine jacobi
  
  

end module jacobi_mod
