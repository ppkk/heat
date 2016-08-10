program heat
  use mpi
  
  implicit none
  
  type mpi_info_t
    integer            :: my_rank, mpi_size, mpi_size_1D
    integer            :: i_pos, j_pos
  end type 

  integer, parameter :: num_iter = 100
  integer, parameter :: N_glob = 100 ! global number of inner points
  integer            :: N            ! local number of inner points
  integer            :: i, j, iter
  integer            :: ierr, neigh_rank
  integer            :: status(MPI_STATUS_SIZE)   
   
  real, allocatable  :: A (:,:), A_new(:,:)
  
  real               :: local_max(1), global_max(1)
  real :: h 
  
  type(mpi_info_t)  :: mpi_info

  
  
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpi_info%my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, mpi_info%mpi_size, ierr)

  mpi_info%mpi_size_1D = int(sqrt(real(mpi_info%mpi_size)))
  call assert(mpi_info%mpi_size_1D * mpi_info%mpi_size_1D == mpi_info%mpi_size)
  
  h = 1./ (N_glob + 1)
  N = N_glob / mpi_info%mpi_size_1D
  call assert(N*mpi_info%mpi_size_1D == N_glob)
 
  
  allocate( A(0:N+1, 0:N+1))
  allocate( A_new(1:N, 1:N))
  A = 0.0

  call rank_to_position(mpi_info%my_rank, mpi_info, mpi_info%i_pos, mpi_info%j_pos)    
    
      
  do iter = 1, num_iter
  
    ! left     
    call send(-1, 0, N, A(1, 1:N), mpi_info)
    call recv(-1, 0, N, A(0, 1:N), mpi_info)
    
    ! right
    call recv(+1, 0, N, A(N+1, 1:N), mpi_info)
    call send(+1, 0, N, A(N  , 1:N), mpi_info)
    
    ! top  
    call send(0, -1, N, A(1:N, 1), mpi_info)
    call recv(0, -1, N, A(1:N, 0), mpi_info)

    ! bottom
    call recv(0, +1, N, A(1:N, N+1), mpi_info)
    call send(0, +1, N, A(1:N, N  ), mpi_info)

    ! left top
    call send(-1, -1, 1, A(1,1), mpi_info)
    call recv(-1, -1, 1, A(0,0), mpi_info)
    
    ! right bottom
    call recv(+1, +1, 1, A(N+1, N+1), mpi_info)
    call send(+1, +1, 1, A(N,   N  ), mpi_info)
    
    ! right top 
    call send(+1, -1, 1, A(N, 1), mpi_info)
    call recv(+1, -1, 1, A(N+1,   0), mpi_info)
    
    ! left bottom
    call recv(-1, +1, 1, A(0, N+1), mpi_info)
    call send(-1, +1, 1, A(1, N  ), mpi_info)
    
!     call print_matrices(N+2, A, mpi_info)
    
    do i = 1, N 
      do j = 1, N               
        A_new(i, j) = 0.25 * (h*h + A(i-1, j) + A(i+1, j) + A(i, j-1) + A(i, j+1) )
      end do      
    end do
    A(1:N, 1:N) = A_new

!      call print_matrices(N+2, A, mpi_info)

  end do

  local_max = maxval(A)
  call MPI_Reduce(local_max, global_max, 1, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  if(mpi_info%my_rank == 0) then
    write(*,*) "maximum is ", global_max
  end if 
  
  call MPI_Finalize(ierr)
!   call neco()
  
  contains
  
  
  subroutine send(i_offs, j_offs, length, buf, mpi_info)
    implicit none
    integer, intent(in)           :: i_offs, j_offs, length
    real, intent(in)              :: buf(length)
    type(mpi_info_t), intent(in)  :: mpi_info
    
    integer        :: neigh_rank
    
    call position_to_rank(mpi_info%i_pos + i_offs, mpi_info%j_pos + j_offs, mpi_info, neigh_rank)
    if(neigh_rank /= -1) then
      call MPI_Send(buf, length, MPI_REAL, neigh_rank, 0, MPI_COMM_WORLD, ierr)
    end if

  end subroutine send
    
  subroutine recv(i_offs, j_offs, length, buf, mpi_info)
    implicit none
    integer, intent(in)           :: i_offs, j_offs, length
    real, intent(out)              :: buf(length)
    type(mpi_info_t), intent(in)  :: mpi_info
    
    integer        :: neigh_rank
    integer        :: status(MPI_STATUS_SIZE)   

    
    call position_to_rank(mpi_info%i_pos + i_offs, mpi_info%j_pos + j_offs, mpi_info, neigh_rank)
    if(neigh_rank /= -1) then
      call MPI_Recv(buf, length, MPI_REAL, neigh_rank, 0, MPI_COMM_WORLD, status, ierr)
    end if    
    
  end subroutine recv

  subroutine position_to_rank(i, j, mpi_info, rank)
    implicit none
    integer, intent(in)              :: i, j
    integer, intent(out)             :: rank
    type(mpi_info_t), intent(in)     :: mpi_info
    
    if( i < 0 .or. i >= mpi_info%mpi_size_1D .or. j < 0 .or. j >= mpi_info%mpi_size_1D) then
      rank = -1
    else
      rank = i * mpi_info%mpi_size_1D + j
    end if

  end subroutine

  subroutine rank_to_position(rank, mpi_info, i, j)
    implicit none
    integer, intent(in)             :: rank
    integer, intent(out)            :: i, j
    type(mpi_info_t), intent(in)    :: mpi_info

    i = rank / mpi_info%mpi_size_1D
    j = mod(rank, mpi_info%mpi_size_1D)
    
    call assert(i >= 0 .and. i < mpi_info%mpi_size_1D .and. j >= 0 .and. j <= mpi_info%mpi_size_1D)

  end subroutine
    
  subroutine print_matrices(length, mat, mpi_info)
    implicit none
    integer, intent(in)          :: length
    real, intent(in)             :: mat(length, length)
    type(mpi_info_t), intent(in) :: mpi_info
    
    integer         :: buf(1), ierr
    integer         :: status(MPI_STATUS_SIZE)   

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if(mpi_info%my_rank > 0) then
      call MPI_Recv(buf, 1, MPI_INT, mpi_info%my_rank - 1, 0, MPI_COMM_WORLD, status, ierr) 
    end if
    
    write(*,*) "RANK ", mpi_info%my_rank, ": "
    do, i=1,length
      write(*,"(10g13.4)") ( mat(i,j), j=1,length )
    end do
    write(*,*)

    if(mpi_info%my_rank < mpi_info%mpi_size - 1) then
      call MPI_Send(buf, 1, MPI_INT, mpi_info%my_rank + 1, 0, MPI_COMM_WORLD, ierr) 
    end if
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    !      if(mpi_info%my_rank == 0) then
!        write(*,*) mat
!      end if
    
  end subroutine 
    
  subroutine assert(expr)
    implicit none
    logical, intent(in)  :: expr
    
    if (expr .eqv. .false.) then
      write(*,*) "assertion failed "
      call exit(0)
    end if
  end subroutine assert

end program heat
