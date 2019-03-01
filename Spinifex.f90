module Spinifex

    use :: ISO_Precisions 
    use :: MPI
    use :: omp_lib
    use :: CSR_type
    use :: Nullarbor

    contains

    subroutine Spinifex_Import_Vector(vec, filename)

    !************************************************************************80
    !
    !! Spinifex_Import_Vector
    !
    ! Modified:
    !
    !   28 October 2018
    !
    ! Description:
    !
    !   Opens and partitions a 1D dense array via MPI-I/O.
    !
    !   'filename’ points to an unformatted binary array file with the
    !   following data structure:
    !
    !       integer (4) ::  number of elements
    !       complex(dp) :: elements
    !
    ! Assumptions:
    !
    !  Subroutine must be called in MPI instance.


        complex(dp), intent(out), dimension(:), allocatable :: vec
        character(len=*), intent(in) :: filename

        integer(sp) :: vec_size

        ! MPI environment.
        integer(sp) :: ierr
        integer(sp), dimension (MPI_STATUS_SIZE) :: status
        integer(kind = MPI_OFFSET_KIND) :: disp
        integer(sp) :: rank
        integer(sp) :: flock ! Number of MPI instances.
        integer(sp) :: vec_file_ref

        call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, flock, ierr )

        ! Open binary vec file.
        call MPI_FILE_OPEN( MPI_COMM_WORLD, filename, &
                    & MPI_MODE_RDONLY, MPI_INFO_NULL, vec_file_ref, ierr )

        ! Obtain vec length from header.
        call MPI_FILE_READ( vec_file_ref, vec_size, 1, MPI_INT, status, ierr )

        allocate (vec(vec_size))

        ! Read partitioned vec from disp, skipping header integer.
        disp = 4

        call MPI_FILE_READ_AT( vec_file_ref, disp, vec, &
                & vec_size, MPI_DOUBLE_COMPLEX, status, ierr )

        call MPI_FILE_CLOSE( vec_file_ref, ierr )

    end subroutine Spinifex_Import_Vector

    subroutine Spinifex_Sparse_Read_and_Distribute ( A_local, &
            & row_offset, partition_table, filename )

    !************************************************************************80
    !
    !! Spinifex_Sparse_Read_and_Distribute
    !
    ! Modified:
    !
    !   21 October 2018
    !
    ! Description:
    !
    !   Opens and partitions Sparse array via MPI. Returns Sparse ‘A_local’,
    !   its row offset ‘row_offset’ relative to its non-partitioned parent and
    !   an array 'partition table' containing the starting row index of each
    !   partition ordered by rank.
    !
    !   'filename’ points to an unformatted binary Sparse file with the
    !   following data structure:
    !
    !       integer ( 4 ) ::  header array ( nonzeros, columns, rows)
    !       integer ( 4 ) ::  row pointer array ‘ir’
    !       integer ( 4 ) :: column index array ‘jc’
    !       complex ( 8 ) :: complex numbers array 'num'
    !
    !   Sparse array is partitioned row-wise with the number of ('A_local')
    !   partitions equal to the number of MPI instances.
    !
    ! Assumptions:
    !
    !   Subroutine must be called in MPI instance.
    !

        type(CSR(dp)), intent(out) :: A_local
        integer(sp), intent(out) :: row_offset
        integer(sp), intent(out), dimension (:), allocatable :: partition_table
        character (len = *), intent(in) :: filename

        integer(sp), dimension(3) :: header
        integer(sp) :: partition_size
        integer(sp) :: jc_and_num_size
        integer(sp) :: realign_ir
        integer(sp) :: i

        ! MPI environment.
        integer(sp) :: ierr
        integer(sp), dimension ( MPI_STATUS_SIZE ) :: status
        integer(kind = MPI_OFFSET_KIND) :: disp
        integer(sp) :: rank
        integer(sp) :: flock ! Number of MPI instances.
        integer(sp) :: sparse_file_ref

        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, flock, ierr)

       if (allocated(A_local%row_start)) then

           deallocate (A_local%row_start, A_local%col_index, A_local%values)

        end if

        write(*,*) "Rank: ", rank, "Importing: ", filename

        ! Open binary Sparse file.
        call MPI_FILE_OPEN ( MPI_COMM_WORLD, filename, &
                    & MPI_MODE_RDONLY, MPI_INFO_NULL, sparse_file_ref, ierr )

        ! Number of nonzeros, columns and rows stored in array "header".
        call MPI_FILE_READ_ALL ( sparse_file_ref, header, 3, MPI_INT, status, ierr )


        ! Determine file offset for ir array in bytes.
        if ( rank .ne. (flock - 1) ) then
            partition_size = header(2)/flock
            row_offset = partition_size*rank
            disp = (3 + row_offset)*4
        else
            partition_size = header(2) - (header(2)/flock)*(flock - 1)
            row_offset = header(3) - partition_size
            disp = (3 + row_offset)*4
        end if

        allocate ( partition_table ( flock + 1) )

        do i = 1, flock
            partition_table(i) = 1 + (header(2)/flock)*(i - 1)
        end do

        partition_table(flock + 1) = header(2) + 1

        allocate ( A_local%row_start (partition_size + 1) )

        ! Read partitioned ir array from disp.
        call MPI_FILE_READ_AT ( sparse_file_ref, disp, A_local%row_start, &
                & partition_size + 1, MPI_INT, status, ierr )

        ! Determine size of partitioned jc and num arrays.
        jc_and_num_size =  A_local%row_start(partition_size + 1) - A_local%row_start(1)

        allocate ( A_local%col_index ( jc_and_num_size ), &
                & A_local%values ( jc_and_num_size ) )

        ! Determine file offset for partitioned jc array in bytes.
        disp = (3 + header(2) + A_local%row_start(1))*4

        ! Read partitioned jc array from disp.
        call MPI_FILE_READ_AT ( sparse_file_ref, disp, A_local%col_index, &
                & jc_and_num_size, MPI_INT, status, ierr )

        ! Determine file offset for partitioned num array in bytes.
        disp =  (3 + (header(3) + 1) + header(1))*4 + ((A_local%row_start(1) - 1)*2)*8

        call MPI_FILE_READ_AT ( sparse_file_ref, disp, A_local%values, &
                    & jc_and_num_size, MPI_DOUBLE_COMPLEX, status, &
                            & ierr )

        call MPI_FILE_CLOSE ( sparse_file_ref, ierr )

        ! Determine offset in ir array required to access the partitioned jc and
        ! ir arrays as per the standard Sparse format.
        realign_ir = A_local%row_start(1) - 1

        ! Offset partitioned ir array.
        A_local%row_start = A_local%row_start - realign_ir

        ! Partitioned array dimensions.
        A_local%columns = header(2)
        A_local%rows = partition_size

    end subroutine Spinifex_Sparse_Read_and_Distribute

    subroutine Spinifex_SpMV ( A_local, partition_table, row_offset, vec_local,&
                & reconciled_sends, new )

    !************************************************************************80
    !
    !! Spinifex_Sparse_SpMV
    !
    ! Modified:
    !
    !   28 October 2018
    !
    ! Description:
    !
    !    Performs a distributed multiplication of Sparse matrix 'A' with dense
    !    vector 'v'. For subsequent multiplications, only locally required
    !    updated values are received.
    !
    ! Assumptions:
    !
    !   Subroutine must be called in MPI instance.
    !   'A' has been imported using 'Spinifex_Sparse_Read_and_Distribute'.
    !   'v' has been imported using 'Spinifex_Import_Vector'.
    !
    !
    !


        ! In and Out
        type(CSR(dp)), intent(in) :: A_local
        integer(sp), dimension(:), intent(in) :: partition_table
        integer(sp), intent(in) :: row_offset
        complex(dp), dimension(:), intent(inout) :: vec_local

        integer(sp), dimension(:,:), allocatable, intent(inout) :: reconciled_sends
        logical, intent(inout) :: new

        ! Local MPI environment
        integer(sp), dimension(MPI_STATUS_SIZE) :: status
        integer(sp) :: ierr
        integer(sp) :: rank
        integer(sp) :: flock ! Number of nodes
        integer(sp) :: position
        integer(sp), dimension(:), allocatable :: requests_in, requests_out
        integer(sp) :: buffers_size
        character(len=1000), dimension(:), allocatable :: buffers_in, buffers_out

        ! Local
        integer(sp) :: number_nz_columns, max_number_nz_columns
        integer(sp), dimension(:,:), allocatable :: requested_indexes
        integer(sp), dimension(:), allocatable :: indexes_out
        integer(sp) :: index
        integer(sp) :: i, j
        integer(sp) :: in_length
        integer(sp), dimension(:), allocatable :: indexes_in
        complex(dp), dimension(:), allocatable :: incom_values
        complex(dp), dimension(:), allocatable :: values_out, values_in
        complex(dp), dimension(:), allocatable :: vec_local_temp
        logical, dimension(:), allocatable :: populated_columns

        !call MPI_INIT (ierr)

        call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, flock, ierr )

        if (new) then

            ! Find and count all non-zero columns in A_local.

            allocate(populated_columns(A_local%columns))

            populated_columns = .false.

            do i = 1, size(A_local%col_index)
                populated_columns(A_local%col_index(i)) = .true.
            end do

            number_nz_columns = 0

            do i = 1, A_local%columns
                if (populated_columns(i)) then
                    number_nz_columns = number_nz_columns + 1
                end if
            end do

            ! Get the system wide maximum of 'number_nz_columns'.

            call MPI_ALLREDUCE (number_nz_columns, max_number_nz_columns, 1, &
                        & MPI_INT, MPI_MAX, MPI_COMM_WORLD, ierr )

            ! max_number_nz_columns gives a size limit for 'requested_indexes'.
            allocate(requested_indexes(flock, max_number_nz_columns))

            requested_indexes = 0

            ! Populate row of 'requested_indexes' corresponding to rank + 1
            ! with all non-zero A_local columns.

            index = 1
            do i = 1, A_local%columns
                if (populated_columns(i)) then
                    requested_indexes(rank+1, index) = i
                    index = index + 1
                end if
            end do

            ! Fully populate 'requested_indexes' across all nodes.
            do i = 1, flock
                call MPI_BCAST ( requested_indexes(i,:), &
                    max_number_nz_columns, MPI_INT, i-1, MPI_COMM_WORLD, ierr )
            end do

            ! This will contain a lookup table of how many updated vec_local
            ! values should be sent to each rank and their indexes.
            allocate(reconciled_sends(flock,A_local%rows + 1))

            reconciled_sends = 0

            ! Scan 'requested_indexes', find which entries fall within the local
            ! matrix/vector partition, place in 'reconciled_sends'.
            do j = 1, flock
                index = 1
                do i = 1, max_number_nz_columns
                    if (requested_indexes(j,i).eq.0) then
                        exit
                    else if(requested_indexes(j,i).ge.partition_table(rank+2)) then
                        exit
                    else if (requested_indexes(j,i).ge.partition_table(rank+1)) then
                            reconciled_sends(j,index+1) = &
                                    & requested_indexes(j, i)
                           index = index + 1
                    end if
                end do
                ! The number of reconciled sends per node.
                reconciled_sends(j,1) = index - 1
            end do

            deallocate(requested_indexes)

            new = .false.

        endif

        allocate(vec_local_temp(size(vec_local)))
        vec_local_temp = 0


        ! Perform multiplication of A_local on vec_local.
        ! row offset is used to populate the correct vector index.
        !$omp parallel do schedule(guided) private(j) firstprivate(row_offset)
        do i = 1, A_local%rows
            do j = A_local%row_start(i), A_local%row_start(i+1) - 1
                vec_local_temp(i+row_offset) = vec_local_temp(i+row_offset) &
                        & + A_local%values(j)*vec_local(A_local%col_index(j))
            end do
        end do
        !$omp end parallel do

        vec_local = vec_local_temp

        deallocate(vec_local_temp)

        ! For each non-local rank, using 'reconciled_sends' allocate
        ! indexes_out and values_out. Populate them with data to send to
        ! rank = i - 1, then MPI pack data and send.
        buffers_size = 1000
        allocate(requests_out(flock),buffers_out(flock),buffers_in(flock), &
                    & requests_in(flock))
        do i = 1, flock

            if (rank.ne.i-1) then

                allocate(indexes_out(reconciled_sends(i,1)+1),&
                        & values_out(reconciled_sends(i,1)))

                do j = 1, reconciled_sends(i,1) + 1

                    indexes_out(j) = reconciled_sends(i,j)

                end do

                do j = 1, reconciled_sends(i,1)

                    values_out(j) = vec_local(reconciled_sends(i,j + 1))

                end do

                position = 0
                call MPI_PACK(indexes_out, indexes_out(1) + 1, MPI_INT, &
                        & buffers_out(i), 1000, position, &
                                & MPI_COMM_WORLD, ierr)
                call MPI_PACK(values_out, indexes_out(1), MPI_DOUBLE_COMPLEX, &
                        & buffers_out(i), 1000, position, &
                                & MPI_COMM_WORLD, ierr)
                call MPI_ISEND(buffers_out(i), 1000, MPI_PACKED, i-1, &
                        & (i-1)+rank*10, MPI_COMM_WORLD, requests_out(i), ierr)

                deallocate (indexes_out, values_out)
            end if
        end do

    do i = 1, flock

        if (rank.ne.i-1) then

            position = 0

            call MPI_RECV(buffers_in(i), 1000, MPI_PACKED, i-1, 10*(i-1)+rank, &
                & MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            call MPI_UNPACK(buffers_in(i),1000,position,in_length,1,MPI_INT,&
                & MPI_COMM_WORLD,ierr)

            allocate(indexes_in(in_length),values_in(in_length))

            values_in = 0

            call MPI_UNPACK(buffers_in(i),1000,position,indexes_in, &
                & in_length,MPI_INT,MPI_COMM_WORLD,ierr)

            call MPI_UNPACK(buffers_in(i),1000,position,values_in, &
                & in_length,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)

            do j = 1, in_length
                vec_local(indexes_in(j)) = values_in(j)
            end do

            deallocate(indexes_in, values_in)

        end if

    end do

    end subroutine Spinifex_SpMV

    function Spinifex_Sparse_OneNormEst(A_local, A_Dagger, A_rows, row_offset, &
                & partition_table, n)

    !   A streamlined version of the one norm estimation algorithm used in Nullarbor.mod

        real(dp) :: Spinifex_Sparse_OneNormEst
        type(CSR(dp)), intent(in) :: A_local, A_Dagger
        integer, intent(in) :: A_rows !No. of rows accross all partitions.
        integer, intent(in) :: row_offset
        integer, dimension(:), intent(in) :: partition_table
        integer, intent(in) :: n

        real(dp) :: est, est_old, Z_norm_max, Z_norm_temp
        real(dp), dimension(:), allocatable ::Z_norm_array
        complex(dp), dimension(:), allocatable :: Y, S, S_Old, Z, Y_temp
        logical, dimension(:), allocatable :: ind_i_hist
        logical :: not_sorted, all_done
        integer, dimension(:), allocatable :: ind_i, v
        integer :: ind
        integer :: ind_i_temp, itmax
        integer :: i, j, k
        integer :: exit_onenormest

        integer, dimension(:,:), allocatable :: reconciled_sends_1, reconciled_sends_2
        logical :: new_1, new_2

        ! MPI environment.
        integer :: ierr
        integer, dimension (MPI_STATUS_SIZE) :: status
        integer :: rank
        integer :: flock ! Number of MPI instances.
        integer :: MASTER = 0
        integer, dimension(:), allocatable :: block_lens, disps

        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        call MPI_COMM_SIZE( MPI_COMM_WORLD, flock, ierr )

        ! Allocate and populate arrays to gather Spinifex_SpMV result.
        allocate(block_lens(flock), disps(flock))

        do j = 1, flock

           block_lens(j) = partition_table(j+1)-partition_table(j)

        end do

        do j = 1, flock

            disps(j) = partition_table(j) - 1

        end do

        est_old = 0
        itmax = 15
        new_1 = .true.
        new_2 = .true.

        allocate(Z_norm_array(A_rows), &
                & Y(A_rows), S(A_rows), S_old(A_rows), &
                        & Z(A_rows), ind_i_hist(A_rows), &
                                & ind_i(A_rows), Y_temp(A_rows))

        ind_i_hist = .true.
        exit_onenormest = 1

        Y_temp = 1.d0/cmplx(A_rows,A_rows,8)  !Allocate chunk at each node length A_local

        k = 1

        do

            do i = 1, n
            call Spinifex_SpMV ( A_local, partition_table, row_offset, Y_temp,&
                  & reconciled_sends_1, new_1 )
            end do

             call MPI_ALLGATHERV(Y_temp(partition_table(rank+1):partition_table(rank+2)-1), &
                & A_local%ROWS, MPI_DOUBLE_COMPLEX, Y, block_lens, disps, &
                    & MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

            if(k.gt.itmax) exit

            Z = 0

            !$omp parallel do schedule(guided)
            do i = partition_table(rank+1), partition_table(rank+2)-1

                if(Y(i).eq.0) then

                    Z(i) = 1

                else

                    Z(i) = Y(i)/abs(Y(i))

                end if
            end do
            !$omp end parallel do

            call MPI_Allreduce( Z, S, A_rows, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)

            ! Performing multiplication of A transposed.
            ! Transposed multiplication is performed row-wise to generate sub-vectors which are
            ! summed to arrive at the final result.
            do i = 1, n

                Z = Spinifex_Sparse_Trans_VecMul(A_Dagger,S,row_offset)
                call MPI_Allreduce( Z, S, A_rows, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)

            end do

            ! START MASTER BLOCK
            if (rank.eq.MASTER) then

                !$omp parallel do
                do i = 1, A_rows

                    Z_norm_array(i) = abs(S(i))
                    ind_i(i) = i

                end do
                !$omp end parallel do

                Z_norm_max = Vector_MaxVal(Z_norm_array)

                ! Should implement a MPI-based merge sort.
                call Vector_Merge_Sort_Ordered_Pair_Real_Integer(Z_norm_array,ind_i,1, &
                        & A_rows)

            end if
            ! END MASTER BLOCK


            est = Z_norm_max

            call MPI_BCAST(est, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, ierr)

            if((k.gt.2).and.(est.le.est_old)) then

                est = est_old

                call MPI_BARRIER(MPI_COMM_WORLD, ierr)

                exit

            end if

            est_old = est

            ! START MASTER BLOCK
            if (rank.eq.MASTER) then

                if(k.gt.1) then

                    all_done = .true.

                    if(ind_i_hist(ind)) then

                        all_done = .false.
                        ind_i_hist(ind) = .false.

                    end if

                    if (all_done) then

                        exit_onenormest = 1

                    end if

                    do j = 1, A_rows !MPI OMP

                        if (ind_i_hist(j)) then

                            ind = j
                            exit_onenormest = 1

                        end if

                    end do

                else

                    ind = 1

                end if

            end if

            call MPI_BCAST(exit_onenormest, 1, MPI_INT, MASTER, MPI_COMM_WORLD, ierr)

            if (exit_onenormest.eq.1) then

                call MPI_BARRIER(MPI_COMM_WORLD, ierr)

                exit

            endif


            Y = 0
            Y(ind) = (1,0)

            k = k + 1

        end do

        Spinifex_Sparse_OneNormEst = est

    end function Spinifex_Sparse_OneNormEst

    function Spinifex_Sparse_VecMul(A_local,v,row_offset)

        ! Computes "A.v = b". Where v is a dense array of length A%COLUMNS.
        ! 'b' entires are offset bu an amount corresponding to the 'A_local'
        ! partition.

        type(CSR(dp)), intent(in) :: A_local
        integer, intent(in) :: row_offset
        complex(dp), dimension(A_local%COLUMNS) :: Spinifex_Sparse_VecMul

        complex(dp), dimension(:), intent(in) :: v

        integer :: i, j

        Spinifex_Sparse_VecMul = 0

        !$omp parallel do schedule(guided)
        do i = 1, A_local%ROWS
            do j = A_local%row_start(i), A_local%row_start(i+1) - 1
                Spinifex_Sparse_VecMul(i+row_offset) = Spinifex_Sparse_VecMul(i+row_offset) + &
                        & A_local%values(j)*v(A_local%col_index(j))
            end do
        end do
        !$omp end parallel do

    end function Spinifex_Sparse_VecMul

    function Spinifex_Sparse_Trans_VecMul(A_T_local,v,row_offset)

        ! Computes "A^T.v = b". Where v is a dense array of length A%COLUMNS.
        ! 'b' entires are offset by an amount corresponding to the 'A_local'
        ! partition.

        type(CSR(dp)), intent(in) :: A_T_local
        integer, intent(in) :: row_offset
        complex(dp), dimension(A_T_local%ROWS) :: Spinifex_Sparse_Trans_VecMul
        complex(dp), dimension(:), intent(in) :: v

        integer :: i, j

        Spinifex_Sparse_Trans_VecMul = 0

        !$omp parallel do schedule(guided)
        do i = 1, A_T_local%ROWS
            do j = A_T_local%row_start(i), A_T_local%row_start(i+1) - 1
                Spinifex_Sparse_Trans_VecMul(i) = Spinifex_Sparse_Trans_VecMul(i) + &
                        & A_T_local%values(j)*v(row_offset+A_T_local%col_index(j))
            end do
        end do
        !$omp end parallel do

    end function Spinifex_Sparse_Trans_VecMul

    function Spinifex_Vector_InfNorm(vec_local, partition_table)

        real(dp) :: Spinifex_Vector_InfNorm_local, Spinifex_Vector_InfNorm
        complex(dp), dimension(:), intent(in) :: vec_local
        integer, dimension(:), intent(in) :: partition_table

        real(dp) :: temp
        integer :: i

        ! MPI environment.
        integer :: ierr
        integer, dimension (MPI_STATUS_SIZE) :: status
        integer :: rank
        integer :: flock ! Number of MPI instances.
        integer :: MASTER

        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        call MPI_COMM_SIZE( MPI_COMM_WORLD, flock, ierr )

        Spinifex_Vector_InfNorm_local = 0

        do i = partition_table(rank+1), partition_table(rank+2)-1
            temp = abs(vec_local(i))
            if(temp.gt.Spinifex_Vector_InfNorm_local)then
                Spinifex_Vector_InfNorm_local = temp
            end if
        end do

        call MPI_Allreduce( Spinifex_Vector_InfNorm_local, Spinifex_Vector_InfNorm, 1, &
                & MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)


    end function Spinifex_Vector_InfNorm

end module Spinifex
