program Example 

    use :: ISO_Precisions 
    use :: CSR_type, only: CSR
    use :: Spinifex
    use :: mpi
    use :: omp_lib

    type(CSR(dp)) :: A

    ! Spinifex_SpMV
    integer(sp), dimension(:), allocatable :: partition_table
    integer(sp) :: row_offset
    complex(dp), dimension(:), allocatable :: vec, vec_local
    integer(sp), dimension(:,:), allocatable :: reconciled_sends
    logical :: new = .true.

    !MPI environment
    integer(sp) :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer(sp) :: rank
    integer(sp) :: flock
    integer, dimension(:), allocatable :: block_lens, disps
    integer :: MASTER = 0

    character(128) :: CSR_filename, vec_filename, output_filename

    real(dp) :: start_time, end_time

    call MPI_INIT(ierr)

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, flock, ierr)

    call get_command_argument(1, CSR_filename)
    call get_command_argument(2, vec_filename)
    call get_command_argument(3, output_filename)

    write(*,*) CSR_filename, vec_filename, output_filename

    call Spinifex_Import_Vector(vec, trim(vec_filename))

    write(*,*) "Rank", rank, "Vector Imported"

    call Spinifex_Sparse_Read_and_Distribute(A, row_offset, partition_table, trim(CSR_filename))

    write(*,*) "Rank", rank, "Matrix Imported"

    start_time = omp_get_wtime()

    if (rank == MASTER) then
        write(*,*) "Calculating (A^10).u = v"
    endif

    do i = 1, 5

        call Spinifex_SpMV(A, partition_table, row_offset, vec, reconciled_sends, new)

    enddo

    end_time = omp_get_wtime()

    if (rank == MASTER) then
        write(*,*) "Completed multiplication. Time: ", end_time - start_time
    endif

    allocate(block_lens(flock), disps(flock))

    do j = 1, flock

         block_lens(j) = partition_table(j+1)-partition_table(j)

    end do

    do j = 1, flock

      disps(j) = partition_table(j) - 1

    end do

    allocate(vec_local(size(vec)))

    call MPI_ALLGATHERV(vec&
            & (partition_table(rank+1):partition_table(rank+2)-1), A%rows, &
                & MPI_DOUBLE_COMPLEX,vec_local, block_lens, &
                    & disps, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

    write(*,*) "Rank", rank, "Result gathered"

    if (rank == MASTER) then
        call Nullarbor_Export_Vector(vec_local, trim(output_filename))
    endif

    if (rank == MASTER) then
        write(*,*) "Result saved."
    endif

    call MPI_FINALIZE(ierr)

end program Example 

