program test3


  use datastructures_mod
  use applyhelmholtz_mod
  use createhelmholtz_mod
  use jacobi_mod
  use multigrid_mod

  implicit none

  type(FieldType) :: u_fine, u_exact
  type(FieldType) :: b_fine
  type(MatrixType) :: H_fine

!  type(FieldType) :: v
!  type(FieldType) :: d
!  real(kind=8) :: res2
  
  real(kind=8) :: nd, x, y, k1, k2
  real(kind=8) :: pi = 4*atan(1.0_8)
  integer :: prec

  call initialise_mg()
 
  nd = n

  k1 = 2.0_8
  k2 = 2.0_8
  prec = precision
  if (prec == 1) then
     allocate(H_fine%data_sgl_D(n,n),H_fine%data_sgl_N(n,n),H_fine%data_sgl_S(n,n))
     allocate(H_fine%data_sgl_W(n,n),H_fine%data_sgl_E(n,n))
     allocate(u_fine%data_sgl(n,n),b_fine%data_sgl(n,n),u_exact%data_sgl(n,n))
  else
     allocate(H_fine%data_dbl_D(n,n),H_fine%data_dbl_N(n,n),H_fine%data_dbl_S(n,n))
     allocate(H_fine%data_dbl_W(n,n),H_fine%data_dbl_E(n,n))
     allocate(u_fine%data_dbl(n,n),b_fine%data_dbl(n,n),u_exact%data_dbl(n,n))
  end if
  !allocate(v%data(n,n),d%data(n,n))

  H_fine%nx = n
  H_fine%ny = n
  H_fine%pr = prec
  
  u_fine%nx = n
  u_fine%ny = n
  u_fine%pr = prec

  u_exact%nx = n
  u_exact%ny = n
  u_exact%pr = prec

  b_fine%nx = n
  b_fine%ny = n
  b_fine%pr = prec
 
  !v%nx = n
  !v%ny = n
  !d%nx = n
  !d%ny = n
  
  call createH(H_fine,1.0_8,1.0_8)


  do i = 1,n
     x = (i-0.5)/nd
     do j = 1,n
        y = (j-0.5)/nd
        if (prec == 1) then
           u_exact%data_sgl(i,j) = sin(k1*pi*x)*sin(k2*pi*y)
        else
           u_exact%data_dbl(i,j) = sin(k1*pi*x)*sin(k2*pi*y)
        end if
     end do
  end do

  if (prec == 1) then
     u_fine%data_sgl(:,:) = 0.0_8
  else
     u_fine%data_dbl(:,:) = 0.0_8
  end if
  call applyH(H_fine,u_exact,b_fine)
  

  call mg_solve(u_fine,b_fine,H_fine)
  
  !call jacobi(H_fine,b_fine,v,0.5_8,100)
  !call applyH(H_fine,v,d)

!  print*,d%data
  call finalise_mg()

  if (prec == 1) then
     deallocate(H_fine%data_sgl_D,H_fine%data_sgl_N,H_fine%data_sgl_S)
     deallocate(H_fine%data_sgl_W,H_fine%data_sgl_E)
     deallocate(u_fine%data_sgl,b_fine%data_sgl,u_exact%data_sgl)
  else
     deallocate(H_fine%data_dbl_D,H_fine%data_dbl_N,H_fine%data_dbl_S)
     deallocate(H_fine%data_dbl_W,H_fine%data_dbl_E)
     deallocate(u_fine%data_dbl,b_fine%data_dbl,u_exact%data_dbl)
  end if
     !deallocate(v%data,d%data)
  
  
end program test3
