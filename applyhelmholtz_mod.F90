module applyhelmholtz_mod


  use datastructures_mod

  implicit none

contains

  ! Apply the Helmholtz operator, i.e. calculate v = H.u
  subroutine applyH(H,u,v) 
    implicit none

    type(MatrixType), intent(in) :: H ! Helmholtz operator
    type(FieldType), intent(in) :: u ! input field u
    type(FieldType), intent(inout) :: v ! output field v
    integer :: nx, ny, i, j ! grid size in x- and y-direction
    integer :: j_plus, j_minus
    nx = H%nx
    ny = H%ny

    ! Helmholtz matrix applied to vector u
    if (H%pr == 1) then
       do j = 1,ny
          j_plus = modulo(j+1-1,ny)+1
          j_minus = modulo(j-1-1,ny)+1
          i = 1
          v%data_sgl(i,j) = H%data_sgl_D(i,j)*u%data_sgl(i,j) &
               + H%data_sgl_N(i,j)*u%data_sgl(i,j_plus) &
               + H%data_sgl_S(i,j)*u%data_sgl(i,j_minus) &
               + H%data_sgl_E(i,j)*u%data_sgl(2,j) &
               + H%data_sgl_W(i,j)*u%data_sgl(nx,j)
          do i = 2,nx-1
             v%data_sgl(i,j) = H%data_sgl_D(i,j)*u%data_sgl(i,j) &
                  + H%data_sgl_N(i,j)*u%data_sgl(i,j_plus) &
                  + H%data_sgl_S(i,j)*u%data_sgl(i,j_minus) &
                  + H%data_sgl_E(i,j)*u%data_sgl(i+1,j) &
                  + H%data_sgl_W(i,j)*u%data_sgl(i-1,j)
          end do
          i = nx
          v%data_sgl(i,j) = H%data_sgl_D(i,j)*u%data_sgl(i,j) &
               + H%data_sgl_N(i,j)*u%data_sgl(i,j_plus) &
               + H%data_sgl_S(i,j)*u%data_sgl(i,j_minus) &
               + H%data_sgl_E(i,j)*u%data_sgl(1,j) &
               + H%data_sgl_W(i,j)*u%data_sgl(nx-1,j)
       end do
    else
       do j = 1,ny
          j_plus = modulo(j+1-1,ny)+1
          j_minus = modulo(j-1-1,ny)+1
          i = 1
          v%data_dbl(i,j) = H%data_dbl_D(i,j)*u%data_dbl(i,j) &
               + H%data_dbl_N(i,j)*u%data_dbl(i,j_plus) &
               + H%data_dbl_S(i,j)*u%data_dbl(i,j_minus) &
               + H%data_dbl_E(i,j)*u%data_dbl(2,j) &
               + H%data_dbl_W(i,j)*u%data_dbl(nx,j)
          do i = 2,nx-1
             v%data_dbl(i,j) = H%data_dbl_D(i,j)*u%data_dbl(i,j) &
                  + H%data_dbl_N(i,j)*u%data_dbl(i,j_plus) &
                  + H%data_dbl_S(i,j)*u%data_dbl(i,j_minus) &
                  + H%data_dbl_E(i,j)*u%data_dbl(i+1,j) &
                  + H%data_dbl_W(i,j)*u%data_dbl(i-1,j)
          end do
          i = nx
          v%data_dbl(i,j) = H%data_dbl_D(i,j)*u%data_dbl(i,j) &
               + H%data_dbl_N(i,j)*u%data_dbl(i,j_plus) &
               + H%data_dbl_S(i,j)*u%data_dbl(i,j_minus) &
               + H%data_dbl_E(i,j)*u%data_dbl(1,j) &
               + H%data_dbl_W(i,j)*u%data_dbl(nx-1,j)
       end do
    end if
  end subroutine applyH

  
end module applyhelmholtz_mod
