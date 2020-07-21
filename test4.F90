program test4


  
  use applyhelmholtz_mod
  use createhelmholtz_mod
  use jacobi_mod
  use multigrid_mod

  implicit none

  type(FieldType) :: u_fine
  type(FieldType) :: b_fine
  type(MatrixType) :: H_fine

  real(kind=8) :: nd, x, y, k1, k2
  real(kind=8) :: pi = 4*atan(1.0_8)

  call initialise_mg()

  nd = n

  k1 = 2.0_8
  k2 = 4.0_8

  if (precision == 1) then
     allocate(H_fine%data_sgl_D(n,n),H_fine%data_sgl_N(n,n),H_fine%data_sgl_S(n,n))
     allocate(H_fine%data_sgl_W(n,n),H_fine%data_sgl_E(n,n))
     allocate(u_fine%data_sgl(n,n),b_fine%data_sgl(n,n))
  else
     allocate(H_fine%data_dbl_D(n,n),H_fine%data_dbl_N(n,n),H_fine%data_dbl_S(n,n))
     allocate(H_fine%data_dbl_W(n,n),H_fine%data_dbl_E(n,n))
     allocate(u_fine%data_dbl(n,n),b_fine%data_dbl(n,n))
  end if
  H_fine%nx = n
  H_fine%ny = n
  H_fine%pr = precision
  u_fine%nx = n
  u_fine%ny = n
  u_fine%pr = precision
  b_fine%nx = n
  b_fine%ny = n
  b_fine%pr = precision
  
  
  call createH(H_fine,alpha,beta)

  if (precision == 1) then
     call random_number(u_fine%data_sgl(:,:))
     u_fine%data_sgl = 0.0_8
  else
     call random_number(u_fine%data_dbl(:,:))
     u_fine%data_dbl = 0.0_8
  end if

  call applyH(H_fine,u_fine,b_fine)

  call mg_solve(u_fine,b_fine,H_fine)


  call finalise_mg()

  if (precision == 1) then
     deallocate(H_fine%data_sgl_D,H_fine%data_sgl_N,H_fine%data_sgl_S,H_fine%data_sgl_W,H_fine%data_sgl_E)
     deallocate(u_fine%data_sgl,b_fine%data_sgl)
  else
     deallocate(H_fine%data_sgl_D,H_fine%data_sgl_N,H_fine%data_sgl_S,H_fine%data_sgl_W,H_fine%data_sgl_E)
     deallocate(u_fine%data_sgl,b_fine%data_sgl)
     
  end if
  
end program test4
