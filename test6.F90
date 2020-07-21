program test6

  
  use applyhelmholtz_mod
  use createhelmholtz_mod
  use jacobi_mod
  use multigrid_mod
  use solver_mod

  implicit none

  real :: start, finish, sum
  
  type(FieldType) :: u_fine
  type(FieldType) :: b_fine
  type(MatrixType) :: H_fine


  real(kind=8) :: nd, x, y, k1, k2, k3, k4, k5, k6

  real(kind=8) :: pi = 4*atan(1.0_8)


  call initialise_mg()


  nd = n

  k1 = 2
  k2 = 2
  k3 = 6
  k4 = 6
  k5 = 10
  k6 = 10



  allocate(H_fine%data_dbl_D(n,n),H_fine%data_dbl_N(n,n),H_fine%data_dbl_S(n,n))
  allocate(H_fine%data_dbl_W(n,n),H_fine%data_dbl_E(n,n))
  allocate(u_fine%data_dbl(n,n),b_fine%data_dbl(n,n))

  H_fine%nx = n
  H_fine%ny = n
  H_fine%pr = 2
  
  u_fine%nx = n
  u_fine%ny = n
  u_fine%pr = 2

  b_fine%nx = n
  b_fine%ny = n
  b_fine%pr = 2
  

  
  call createH(H_fine,alpha,beta)


  call applyH(H_fine,u_fine,b_fine)

  b_fine%data_dbl(:,:) = 0.0_8
  
  do i = 1,n
     x = (i-0.5)/nd
     do j = 1,n
        y = (j-0.5)/nd
        u_fine%data_dbl(i,j) = sin(k1*pi*x)*sin(k2*pi*y)! + sin(k3*pi*x)*sin(k4*pi*y) + sin(k5*pi*x)*sin(k6*pi*y)
     end do
  end do

  
  call cpu_time(start)
  call cg_solver(H_fine,b_fine,u_fine,maxiter,tolerance,precision)
  call cpu_time(finish)
  sum = sum + finish - start
  print '("Time = ",F6.3," seconds")', sum

  call finalise_mg()
  
  deallocate(H_fine%data_dbl_D,H_fine%data_dbl_N,H_fine%data_dbl_S)
  deallocate(H_fine%data_dbl_W,H_fine%data_dbl_E)
  deallocate(u_fine%data_dbl,b_fine%data_dbl)

  
end program test6
