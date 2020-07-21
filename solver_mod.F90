
module solver_mod


  use datastructures_mod
  use jacobi_mod
  use applyhelmholtz_mod
  use multigrid_mod
  
  implicit none

contains

  ! Apply preconditioner.
  ! Approximately solve H.u = b using a few iterations of Jacobi
  subroutine applyPrec(H,b,u)
    implicit none

    type(MatrixType), intent(in) :: H
    type(FieldType), intent(inout) :: b
    type(FieldType), intent(inout) :: u

    if (u%pr == PREC_DOUBLE) then
       u%data_dbl(:,:) = 0.0
    else
       u%data_sgl(:,:) = 0.0
    end if

    call mg_solve(u, b, H)
    
  end subroutine applyPrec
  
  ! Preconditioned conjugate gradient solver for problem
  ! H.u = b
  subroutine cg_solver(H,b,u,maxiter,tolerance,precision)

    implicit none

    type(MatrixType), intent(in) :: H
    type(FieldType), intent(in) :: b
    type(FieldType), intent(inout) :: u
    integer, intent(in) :: maxiter
    real(kind=8), intent(in) :: tolerance
    integer, intent(in) :: precision
    type(FieldType) :: r,z,p,q ! temporary fields
    type(FieldType) :: r_prec,z_prec
    type(MatrixType) :: H_prec
    real(kind=8) :: zdotr, zdotr_old
    real(kind=8) :: alpha, beta
    real(kind=8) :: r_two_norm, r0_two_norm, r_two_norm_old
    
    integer :: nx, ny, k

    real(kind=8) :: before, after, sum
    sum = 0
    
    nx = H%nx
    ny = H%ny

    ! Create temporary fields
    r%nx = nx; r%ny = ny; r%pr = 2; allocate(r%data_dbl(nx,ny))
    p%nx = nx; p%ny = ny; p%pr = 2; allocate(p%data_dbl(nx,ny))
    q%nx = nx; q%ny = ny; q%pr = 2; allocate(q%data_dbl(nx,ny))
    z%nx = nx; z%ny = ny; z%pr = 2; allocate(z%data_dbl(nx,ny))

    ! Create fields for single precision
    if (precision == PREC_SINGLE) then
       r_prec%nx = nx; r_prec%ny = ny; z_prec%pr = 1; allocate(r_prec%data_sgl(nx,ny))
       z_prec%nx = nx; z_prec%ny = ny; z_prec%pr = 1; allocate(z_prec%data_sgl(nx,ny)) 
       H_prec%nx = nx; H_prec%ny = ny; H_prec%pr = 1
       allocate(H_prec%data_sgl_N(nx,ny),H_prec%data_sgl_S(nx,ny),H_prec%data_sgl_E(nx,ny))
       allocate(H_prec%data_sgl_W(nx,ny),H_prec%data_sgl_D(nx,ny))
       H_prec%data_sgl_D(:,:) = H%data_dbl_D(:,:)
       H_prec%data_sgl_N(:,:) = H%data_dbl_N(:,:)
       H_prec%data_sgl_E(:,:) = H%data_dbl_E(:,:)
       H_prec%data_sgl_S(:,:) = H%data_dbl_S(:,:)
       H_prec%data_sgl_W(:,:) = H%data_dbl_W(:,:)
    end if
       
    ! Apply operator q = H.u
    call applyH(H,u,q)

    ! Calculate residual r = b - q = b - H.u
    r%data_dbl(:,:) = b%data_dbl(:,:) - q%data_dbl(:,:)
    ! Apply preconditioner
    call cpu_time(before)
    if (precision == PREC_SINGLE) then
       r_prec%data_sgl(:,:) = r%data_dbl(:,:)
       call applyPrec(H_prec,r_prec,z_prec)
       z%data_dbl(:,:) = z_prec%data_sgl(:,:)
    else
       call applyPrec(H,r,z)
    end if
    call cpu_time(after)
    sum = after-before+sum
    ! Write to file for analysis
    open(5, file = "results.txt")
    p%data_dbl(:,:) = z%data_dbl(:,:)
    ! Calculate <z,r>
    zdotr_old = dot(z,r)
    ! Calculate initial two-norm of residual
    r0_two_norm = two_norm(r)
    r_two_norm_old = r0_two_norm
    write (*, '(" initial residual : ",E12.4)') r0_two_norm
    write (*, '(" ",A4," : ",A12,"  ",A12,"  ",A6," ",A9)') &
         "iter","||r||","||r||/||r_0||","rho","prec time"
    do k=1,maxiter
       ! Apply operator q = H.p
       call applyH(H,p,q)
       ! Calculate alpha = <z_{old},r_{old}> / <p,q>
       alpha = zdotr_old/dot(p,q)
       ! axpy update of solution: u <- u + alpha*p
       ! axpy update of residual: r <- r - alpha*q
       u%data_dbl(:,:) = u%data_dbl(:,:) + alpha*p%data_dbl(:,:)
       r%data_dbl(:,:) = r%data_dbl(:,:) - alpha*q%data_dbl(:,:)
       ! Calculate norm of residual
       r_two_norm = two_norm(r)
       write (*, '(" ",I4," : ",E12.4,"  ",E12.4,"  ",F6.3," ",F6.3)') &
            k, r_two_norm, r_two_norm/r0_two_norm, r_two_norm/r_two_norm_old, after-before
       write(5,*) r_two_norm/r0_two_norm
       ! Exit loop if residual norm is sufficiently small
       if (r_two_norm/r0_two_norm < tolerance) exit
       ! Apply preconditioner
       call cpu_time(before)
       if (precision == 1) then
          r_prec%data_sgl(:,:) = r%data_dbl(:,:)
          call applyPrec(H_prec,r_prec,z_prec)
          z%data_dbl(:,:) = z_prec%data_sgl(:,:)
       else
          call applyPrec(H,r,z)
       end if
       call cpu_time(after)
       sum = after-before+sum
       ! Calculate <z,r>
       zdotr = dot(z,r)
       beta = zdotr/zdotr_old
       ! aypx update of p: p <- z + beta*p
       p%data_dbl(:,:) = z%data_dbl(:,:) + beta*p%data_dbl(:,:)
       zdotr_old = zdotr
       r_two_norm_old = r_two_norm
    end do
    print*,"prec time: ", sum
    ! Deallocate temporary vectors
    deallocate(r%data_dbl)
    deallocate(p%data_dbl)
    deallocate(q%data_dbl)
    deallocate(z%data_dbl)

    ! Deallocate single precision vectors
    if (precision == PREC_SINGLE) then
       deallocate(z_prec%data_sgl)
       deallocate(r_prec%data_sgl)
       deallocate(H_prec%data_sgl_D,H_prec%data_sgl_N,H_prec%data_sgl_W)
       deallocate(H_prec%data_sgl_S,H_prec%data_sgl_E)
    end if
    
  end subroutine cg_solver

  subroutine loop_solver(H,b,u,maxiter,tolerance,precision)

    type(MatrixType), intent(in) :: H
    type(FieldType), intent(in) :: b
    type(FieldType), intent(inout) :: u
    integer, intent(in) :: maxiter, precision
    real(kind=8), intent(in) :: tolerance

    type(FieldType) :: r,z,r_prec,z_prec
    type(MatrixType) :: H_prec
    integer :: nx, ny, k
    real(kind=8) :: r_two_norm, r0_two_norm, r_two_norm_old
    real(kind=8) :: after, before
    
    nx = H%nx
    ny = H%ny

    ! Create temporary fields
    r%nx = nx; r%ny = ny; r%pr = 2; allocate(r%data_dbl(nx,ny))
    z%nx = nx; z%ny = ny; z%pr = 2; allocate(z%data_dbl(nx,ny))

    ! Allocate fields for single precision
    if (precision == PREC_SINGLE) then
       r_prec%nx = nx; r_prec%ny = ny; z_prec%pr = 1; allocate(r_prec%data_sgl(nx,ny))
       z_prec%nx = nx; z_prec%ny = ny; z_prec%pr = 1; allocate(z_prec%data_sgl(nx,ny)) 
       H_prec%nx = nx; H_prec%ny = ny; H_prec%pr = 1
       allocate(H_prec%data_sgl_N(nx,ny),H_prec%data_sgl_S(nx,ny),H_prec%data_sgl_E(nx,ny))
       allocate(H_prec%data_sgl_W(nx,ny),H_prec%data_sgl_D(nx,ny))
       H_prec%data_sgl_D(:,:) = H%data_dbl_D(:,:)
       H_prec%data_sgl_N(:,:) = H%data_dbl_N(:,:)
       H_prec%data_sgl_E(:,:) = H%data_dbl_E(:,:)
       H_prec%data_sgl_S(:,:) = H%data_dbl_S(:,:)
       H_prec%data_sgl_W(:,:) = H%data_dbl_W(:,:)
    end if

    call applyH(H,u,r)
    r%data_dbl = b%data_dbl - r%data_dbl
    r0_two_norm = two_norm(r)
    r_two_norm = two_norm(r)
    r_two_norm_old = r0_two_norm
    write (*, '(" initial residual : ",E12.4)') r0_two_norm
    write (*, '(" ",A4," : ",A12,"  ",A12,"  ",A6," ",A9)') &
         "iter","||r||","||r||/||r_0||","rho","prec time"
    open(5, file = "results2.txt")
    
    do k=1,maxiter
       ! Exit loop if residual norm is sufficiently small
       if (r_two_norm/r0_two_norm < tolerance) exit
       call cpu_time(before)
       if (precision == 1) then
          r_prec%data_sgl(:,:) = r%data_dbl(:,:)
          call applyPrec(H_prec,r_prec,z_prec)
          z%data_dbl(:,:) = z_prec%data_sgl(:,:)
       else
          call applyPrec(H,r,z)
       end if
       call cpu_time(after)
       u%data_dbl = u%data_dbl + z%data_dbl
       r_two_norm_old = r_two_norm
       call applyH(H,u,r)
       r%data_dbl = b%data_dbl - r%data_dbl
       r_two_norm = two_norm(r)
       write (*, '(" ",I4," : ",E12.4,"  ",E12.4,"  ",F6.3," ",F6.3)') &
            k, r_two_norm, r_two_norm/r0_two_norm, r_two_norm/r_two_norm_old, after-before
       write(5,*) r_two_norm/r0_two_norm
    end do
    
    ! Deallocate temporary vectors
    deallocate(r%data_dbl)
    deallocate(z%data_dbl)
    
    ! Deallocate single precision vectors
    if (precision == PREC_SINGLE) then
       deallocate(z_prec%data_sgl)
       deallocate(r_prec%data_sgl)
       deallocate(H_prec%data_sgl_D,H_prec%data_sgl_N,H_prec%data_sgl_W)
       deallocate(H_prec%data_sgl_S,H_prec%data_sgl_E)
    end if
    
  end subroutine loop_solver
  
end module solver_mod
