module datastructures_mod

  integer, parameter :: PREC_SINGLE = 1
  integer, parameter :: PREC_DOUBLE = 2
  
  integer, parameter :: sgl_kind = 4
  integer, parameter :: dbl_kind = 8

! type for storing discretised fields on a 2d grid
type FieldType
   integer :: nx ! number of grid cells in x-direction
   integer :: ny ! number of grid cells in y-direction
   integer :: pr ! 1 = single, 2 = double
   real(kind=sgl_kind), dimension(:,:), allocatable :: data_sgl
   real(kind=dbl_kind), dimension(:,:), allocatable :: data_dbl
end type FieldType

! type for storing the components of the Helmholtz matrix
type MatrixType
   integer :: nx ! number of grid cells in x-direction
   integer :: ny ! number of grid cells in y-direction
   integer :: pr ! 1 = single, 2 = double
   real(kind=sgl_kind), dimension(:,:), allocatable :: data_sgl_D ! diagonal entry
   real(kind=sgl_kind), dimension(:,:), allocatable :: data_sgl_N ! North
   real(kind=sgl_kind), dimension(:,:), allocatable :: data_sgl_S ! South
   real(kind=sgl_kind), dimension(:,:), allocatable :: data_sgl_E ! East
   real(kind=sgl_kind), dimension(:,:), allocatable :: data_sgl_W ! West
   real(kind=dbl_kind), dimension(:,:), allocatable :: data_dbl_D ! diagonal entry
   real(kind=dbl_kind), dimension(:,:), allocatable :: data_dbl_N ! North
   real(kind=dbl_kind), dimension(:,:), allocatable :: data_dbl_S ! South
   real(kind=dbl_kind), dimension(:,:), allocatable :: data_dbl_E ! East
   real(kind=dbl_kind), dimension(:,:), allocatable :: data_dbl_W ! West
end type MatrixType


contains

  ! Calculate the two-norm ||u||_2 of a discretised field
  function two_norm(u) result(nrm)
    implicit none
    type(FieldType), intent(in) :: u
    real(kind=8) :: nrm, nrm2
    integer :: i,j ! loop indices

    nrm2 = 0.0
    if (u%pr == 1) then
       do j=1,u%ny
          do i=1,u%nx
             nrm2 = nrm2 + u%data_sgl(i,j)**2
          end do
       end do
    else
       do j=1,u%ny
           do i=1,u%nx
             nrm2 = nrm2 + u%data_dbl(i,j)**2
          end do
       end do
    end if
    nrm = sqrt(nrm2)
  end function two_norm

  
  ! Calculate the dot-product of two discretised fields
  function dot(u,v) result(dot_uv)
    implicit none
    type(FieldType), intent(in) :: u
    type(FieldType), intent(in) :: v
    real(kind=8) :: dot_uv
    integer :: i,j ! loop indices

    dot_uv = 0.0
    if (u%pr == 1) then
       do j=1,u%ny
          do i=1,u%nx
             dot_uv = dot_uv + u%data_sgl(i,j)*v%data_sgl(i,j)
          end do
       end do
    else
       do j=1,u%ny
          do i=1,u%nx
             dot_uv = dot_uv + u%data_dbl(i,j)*v%data_dbl(i,j)
          end do
       end do
    end if
  end function dot

  
end module datastructures_mod
