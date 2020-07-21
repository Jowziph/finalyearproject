module multigrid_mod

  use createhelmholtz_mod
  use datastructures_mod
  use applyhelmholtz_mod
  use jacobi_mod
  implicit none

  
  integer :: npre, npost, ncoarse
  integer :: nlevel
  integer :: mg_iter, maxiter
  real(kind=8) :: omega, res, rho, res_old, res0, alpha, beta, tolerance
  integer :: i, j, n 
  integer :: dim, precision

  integer, parameter :: fileid = 10

  
  type(FieldType), dimension(:), allocatable :: U
  type(FieldType), dimension(:), allocatable :: b
  type(FieldType), dimension(:), allocatable :: r
  type(MatrixType), dimension(:), allocatable :: H


contains

  subroutine initialise_mg()

    implicit none
    
    namelist /parameters/ npre, npost, ncoarse, nlevel, mg_iter, omega, alpha, beta, dim, precision, tolerance, maxiter
    open(fileid,file='parameters.dat')
    read(fileid,nml=parameters)
    write(*,*) "npre      = ",npre
    write(*,*) "npost     = ",npost
    write(*,*) "ncoarse   = ",ncoarse
    write(*,*) "nlevel    = ",nlevel
    write(*,*) "mg_iter   = ",mg_iter
    write(*,*) "omega     = ",omega
    write(*,*) "alpha     = ",alpha
    write(*,*) "beta      = ",beta
    write(*,*) "dim       = ",dim
    write(*,*) "precison  = ",precision
    write(*,*) "tolerance = ",tolerance
    write(*,*) "maxiter   = ",maxiter
    close(fileid)


    allocate(u(nlevel))
    allocate(b(nlevel))
    allocate(r(nlevel))
    allocate(H(nlevel))

    !allocate and set dimensions of H,u,b,r at each level
    do i=1, nlevel
       n = dim * 2.0_8**(i-nlevel)
       if (precision == 1) then
          allocate(H(i)%data_sgl_N(n,n),H(i)%data_sgl_S(n,n),H(i)%data_sgl_E(n,n))
          allocate(H(i)%data_sgl_W(n,n),H(i)%data_sgl_D(n,n))
          allocate(u(i)%data_sgl(n,n),b(i)%data_sgl(n,n),r(i)%data_sgl(n,n))
       else
          allocate(H(i)%data_dbl_N(n,n),H(i)%data_dbl_S(n,n),H(i)%data_dbl_E(n,n))
          allocate(H(i)%data_dbl_W(n,n),H(i)%data_dbl_D(n,n))
          allocate(u(i)%data_dbl(n,n),b(i)%data_dbl(n,n),r(i)%data_dbl(n,n))   
       end if
       !H
       H(i)%nx = n; H(i)%ny = n; H(i)%pr = precision; call createH(H(i),alpha,beta)
       u(i)%nx = n; u(i)%ny = n; u(i)%pr = precision
       !b
       b(i)%nx = n; b(i)%ny = n; b(i)%pr = precision
       !r
       r(i)%nx = n; r(i)%ny = n; r(i)%pr = precision
    end do
    call initialise_jacobi(u(nlevel))
  end subroutine initialise_mg
  
  subroutine finalise_mg()
    ! Tidy up, i.e. deallocate data on all levels
    do i = 1, nlevel
       if (precision == 1) then
          deallocate(H(i)%data_sgl_N,H(i)%data_sgl_S,H(i)%data_sgl_E)
          deallocate(H(i)%data_sgl_W,H(i)%data_sgl_D)
          deallocate(u(i)%data_sgl,b(i)%data_sgl,r(i)%data_sgl)
       else
          deallocate(H(i)%data_dbl_N,H(i)%data_dbl_S,H(i)%data_dbl_E)
          deallocate(H(i)%data_dbl_W,H(i)%data_dbl_D)
          deallocate(u(i)%data_dbl,b(i)%data_dbl,r(i)%data_dbl)
       end if
    end do
    deallocate(H,u,b,r)
    call finalise_jacobi()
  end subroutine finalise_mg
  
  subroutine mg_solve(u_fine, b_fine, H_fine)
    ! Approximately solves H_fine.u_fine = b_fine using a couple of
    ! MG V-cycles
    
    Type(FieldType), intent(inout) :: u_fine
    Type(FieldType), intent(inout) :: b_fine
    Type(MatrixType), intent(in) :: H_fine

    !input data is assigned to larger arrays
    
    if (precision == 1) then
       u(nlevel)%data_sgl(:,:) = u_fine%data_sgl(:,:)
       b(nlevel)%data_sgl(:,:) = b_fine%data_sgl(:,:)
    else
       u(nlevel)%data_dbl(:,:) = u_fine%data_dbl(:,:)
       b(nlevel)%data_dbl(:,:) = b_fine%data_dbl(:,:)
    end if
 
    !residual from initial guess is calculated
    call applyH(H(nlevel),u(nlevel),r(nlevel))
    if (precision == 1) then
       r(nlevel)%data_sgl(:,:) = b(nlevel)%data_sgl(:,:) - r(nlevel)%data_sgl(:,:)
    else
       r(nlevel)%data_dbl(:,:) = b(nlevel)%data_dbl(:,:) - r(nlevel)%data_dbl(:,:)
    end if
    res0 = two_norm(r(nlevel))
    res_old = res0
    !write(*,'(I4," : ",E10.4," ",F6.3)') 0, res0, res0/res0

    !main loop applying mg v cycle
    do j=1, mg_iter

       call mg_vcycle(nlevel)


       ! Calculate fine-level residual r = b - H.u
       call applyH(H(nlevel),u(nlevel),r(nlevel))
       if (precision == 1) then
          r(nlevel)%data_sgl(:,:) = b(nlevel)%data_sgl(:,:) - r(nlevel)%data_sgl(:,:)
       else
          r(nlevel)%data_dbl(:,:) = b(nlevel)%data_dbl(:,:) - r(nlevel)%data_dbl(:,:)
       end if
       
       res = two_norm(r(nlevel))
       rho = res/res_old
       !write(*,'(I4," : ",E10.4," ",E10.4," ",F6.3)') j, res, res/res0, rho
       res_old = res
    end do
    !give final estimate then deallocate all
    if (precision == 1) then
       u_fine%data_sgl(:,:) = u(nlevel)%data_sgl(:,:)
    else
       u_fine%data_dbl(:,:) = u(nlevel)%data_dbl(:,:)  
    end if
    
  end subroutine mg_solve
  
  recursive subroutine mg_vcycle(L)

    implicit none

    integer, intent(in) :: L

    !set u to zero
    if (L /= nlevel) then
       if (precision == 1) then
          u(L)%data_sgl(:,:) = 0.0_4
       else
          u(L)%data_dbl(:,:) = 0.0_8
       end if
    end if
    !mg vcycle loop
    if (L /= 1) then
       call jacobi(H(L),b(L),u(L),omega,npre)
       call applyH(H(L),u(L),r(L))
       if (precision == 1) then
          r(L)%data_sgl(:,:) = b(L)%data_sgl(:,:) - r(L)%data_sgl(:,:)
       else
          r(L)%data_dbl(:,:) = b(L)%data_dbl(:,:) - r(L)%data_dbl(:,:)
       end if
       call restrict(r(L),b(L-1),L)
       call mg_vcycle(L-1)
       call prolong(u(L-1),r(L),L-1)
       if (precision == 1) then
          u(L)%data_sgl = u(L)%data_sgl + r(L)%data_sgl
       else
          u(L)%data_dbl = u(L)%data_dbl + r(L)%data_dbl
       end if
       call jacobi(H(L),b(L),u(L),omega,npost)
    else
       call jacobi(H(L),b(L),u(L),omega,ncoarse)
    end if
    
  end subroutine mg_vcycle

  
  subroutine restrict(Uin,Uout,L)

    implicit none

    type(FieldType), intent(in) :: Uin
    type(FieldType), intent(inout) :: Uout
    integer, intent(in) :: L
    integer :: nx, ny, i, j

    nx = Uout%nx
    ny = Uout%ny
    do j = 1,ny
       do i = 1,nx
          if (precision == 1) then
             Uout%data_sgl(i,j) = (Uin%data_sgl(2*i-1,2*j-1) + Uin%data_sgl(2*i-1,2*j) &
                  + Uin%data_sgl(2*i,2*j-1) + Uin%data_sgl(2*i,2*j))
          else
             Uout%data_dbl(i,j) = (Uin%data_dbl(2*i-1,2*j-1) + Uin%data_dbl(2*i-1,2*j) &
                  + Uin%data_dbl(2*i,2*j-1) + Uin%data_dbl(2*i,2*j))
          end if
          end do
    end do
  end subroutine restrict

  subroutine prolong(Uin,Uout,L)

    implicit none

    type(FieldType), intent(in) :: Uin
    type(FieldType), intent(inout) :: Uout
    integer, intent(in) :: L
    integer :: nx, ny, i, j

    nx = Uin%nx
    ny = Uin%ny
    do j = 1,ny
       do i = 1,nx
          if (precision == 1) then
             Uout%data_sgl(2*i-1,2*j-1) = Uin%data_sgl(i,j)
             Uout%data_sgl(2*i-1,2*j)   = Uin%data_sgl(i,j)
             Uout%data_sgl(2*i,2*j-1)   = Uin%data_sgl(i,j)
             Uout%data_sgl(2*i,2*j)     = Uin%data_sgl(i,j)
          else
             Uout%data_dbl(2*i-1,2*j-1) = Uin%data_dbl(i,j)
             Uout%data_dbl(2*i-1,2*j)   = Uin%data_dbl(i,j)
             Uout%data_dbl(2*i,2*j-1)   = Uin%data_dbl(i,j)
             Uout%data_dbl(2*i,2*j)     = Uin%data_dbl(i,j)
          end if
       end do
    end do
  end subroutine prolong

end module multigrid_mod
