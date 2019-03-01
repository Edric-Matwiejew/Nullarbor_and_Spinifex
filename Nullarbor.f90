module Nullarbor

    use :: ISO_Precisions 
    use :: omp_lib, only: omp_get_max_threads, omp_get_thread_num, &
            & omp_get_num_threads
    use :: CSR_type

    implicit none

    private

  ! Public Subroutines.
    public :: Nullarbor_Export_Sparse_Matrix, Nullarbor_Import_Sparse_Matrix
    public :: Nullarbor_Import_Vector, Nullarbor_Export_Vector, Sparse_Allocate
    public :: Sparse_Assignment, Sparse_Deallocate, Vector_MaxValPos_Real
    public :: Vector_Merge_Ordered_Pair_Real_Integer
    public :: Vector_Merge_Sort_Ordered_Pair_Real_Integer
    public ::  Nullarbor_Import_Dense_Matrix
    public :: Sparse_Write, Nullarbor_Export_Unformatted, Vec_to_Binary

    ! Public functions.
    public :: Sparse, SparseDense_MatMul, Sparse_Conjugate, Sparse_OneNormEst
    public :: Sparse_ScaMul_Real, Sparse_Size, Sparse_Transpose, Sparse_VecMul
    public :: Vector_InfNorm, Vector_MaxVal


    ! Public Overloadings.
    public :: operator(*), operator(.dot.), assignment(=)


    ! Operator overloading for CSR type.

    interface operator(*)

        procedure :: Sparse_ScaMul_Real

    end interface operator (*)

    interface operator(.dot.)

        procedure :: Sparse_VecMul, SparseDense_MatMul

    end interface operator (.dot.)

    interface assignment(=)

        module procedure :: Sparse_Assignment

    end interface assignment (=)

    contains

    ! Subroutines.

    subroutine Sparse_Allocate(A, NONZEROS, ROWS, COLUMNS)

        type(CSR(dp)), intent(inout) :: A
        integer(sp), intent(in) :: NONZEROS, ROWS, COLUMNS

        allocate(A%values(NONZEROS), A%col_index(NONZEROS), A%row_start(ROWS+1))

        A%ROWS = ROWS
        A%COLUMNS = COLUMNS

    end subroutine Sparse_Allocate

    subroutine Sparse_Deallocate(A)

        type(CSR(dp)), intent(inout) :: A

        deallocate(A%values, A%col_index, A%row_start)

        A%ROWS = 0
        A%COLUMNS = 0

    end subroutine Sparse_Deallocate

    subroutine Nullarbor_Import_Dense_Matrix(A, filename)

        ! comment on format

        complex(8), dimension(:,:), allocatable, intent(out) :: A
        character(len=*), intent(in) :: filename

        real(dp), dimension(:,:), allocatable :: temp
        integer(sp) :: ROWS, COLUMNS, i, j
        integer(sp) :: dummy

        open(17, file = filename)

        read(17,*) ROWS
        read(17,*) COLUMNS

        allocate(temp(ROWS,COLUMNS*2), A(ROWS,COLUMNS))

        do i = 1, ROWS

            read(17,*) temp(i,:)

        end do

        close(17)

        do j = 1, COLUMNS
            do i = 1, ROWS
                A(i,j) = cmplx(temp(i,(2*j)-1), temp(i,(2*j)),8)
            end do
        end do

    end subroutine Nullarbor_Import_Dense_Matrix

    subroutine Nullarbor_Import_Sparse_Matrix(A, filename)

        type(CSR(dp)), intent(inout) :: A
        character(len=*), intent(in) :: filename

        integer(sp) :: NZ, ROWS, COLUMNS, i

        open(unit = 10, file = trim(filename), status = 'old')

        read(10,*) nz
        read(10,*) ROWS
        read(10,*) COLUMNS

        call Sparse_Allocate(A, NZ, ROWS, COLUMNS)

        do i = 1, size(A%values)
            read(10,*) A%values(i)
        end do

        do i = 1, size(A%col_index)
            read(10,*) A%col_index(i)
        end do

        do i = 1, size(A%row_start)
            read(10,*) A%row_start(i)
        end do

        close(10)

    end subroutine Nullarbor_Import_Sparse_Matrix

    subroutine Nullarbor_Import_Vector(v, filename)

        complex(8), dimension(:), allocatable, intent(out) :: v
        character(len=*), intent(in) :: filename

        integer(sp) :: length
        real(dp), dimension(:,:), allocatable :: temp
        integer(sp) :: i

        open(10,file=filename)

        read(10,*) length

        allocate(v(length),temp(length,2))

        do i = 1, length

            read(10,*) temp(i,:)

        end do

        do i = 1, length

                v(i) = cmplx(temp(i,1), temp(i,2),8)

            end do

    end subroutine Nullarbor_Import_Vector

    subroutine Nullarbor_Export_Sparse_Matrix(A, filename)

        type(CSR(dp)), intent(in) :: A
        character(len=*), intent(in) :: filename

        integer(sp) :: i

        open(unit = 10, file = trim(filename), status = 'replace')

        write(10,*) Sparse_Nonzero(A)

        write(10,*) A%ROWS
        write(10,*) A%COLUMNS

        do i = 1, size(A%values)
            write(10,*) A%values(i)
        end do

        do i = 1, size(A%col_index)
            write(10,*) A%col_index(i)
        end do

        do i = 1, size(A%row_start)
            write(10,*) A%row_start(i)
        end do

        close(10)

    end subroutine Nullarbor_Export_Sparse_Matrix

    subroutine Nullarbor_Export_Vector(v, filename)

        complex(8), dimension(:), intent(in) :: v
        character(len=*), intent(in) :: filename

        integer(sp) :: i

        open(unit = 10, file = trim(filename), status = 'replace')

        do i = 1, size(v)
            write(10,*) v(i)
        end do

        close(10)

    end subroutine Nullarbor_Export_Vector

        subroutine Nullarbor_Export_Unformatted ( A, filename )

        type(CSR(dp)), intent ( in ) :: A
        character ( len = * ), intent ( in ) :: filename
        integer(sp) :: lunit, i

        open ( newunit = lunit , file = filename, status = "replace", access="stream", &
        & form = "unformatted" )

        write ( lunit ) Sparse_Nonzero(A)
        write ( lunit ) A%columns
        write ( lunit ) A%rows

        do i = 1, size ( A%row_start )
            write ( lunit ) A%row_start ( i )
        end do

        do i = 1, size( A%col_index )
            write( lunit ) A%col_index ( i )
        end do

        do i = 1, size ( A%values )
                    write ( lunit ) A%values(i)
        end do

        close (lunit)

    end subroutine Nullarbor_Export_Unformatted

    subroutine Sparse_Assignment(B, A)

        type(CSR(dp)), intent(in) :: A
        type(CSR(dp)), allocatable, intent(out) :: B

        integer(sp) :: i, length

        if(allocated(B%row_start).eqv..true.) then
            call Sparse_Deallocate(B)
        end if

        length = size(A%col_index)

        call Sparse_Allocate(B, length, A%ROWS, A%COLUMNS)

        B%values = A%values
        B%col_index = A%col_index
        B%row_start = A%row_start

        B%COLUMNS = A%COLUMNS
        B%ROWS = A%ROWS

    end subroutine Sparse_Assignment

    subroutine Vector_MaxValPos_Real(v,maximum,location)

        real(dp), dimension(:), intent(in) :: v
        real(dp), intent(out) :: maximum
        integer(sp), intent(out) :: location

        real(dp), dimension(:), allocatable :: maxs
        real(dp) :: local_max, local_loc
        integer(sp), dimension(:), allocatable :: locs
        integer(sp) :: i, start, finish, nthrds, thrdID, length

        length = size(v)

        nthrds = omp_get_max_threads()

        allocate(maxs(nthrds),locs(nthrds))

        !$omp parallel private(thrdID,i,start,finish,local_max,local_loc) &
        !$omp   shared (nthrds,maxs,locs)
        thrdID = omp_get_thread_num()
        start = (thrdID*length/nthrds) + 1
        finish = (thrdID+1)*length/nthrds
        if((thrdID+1).eq.nthrds) then
            finish = length
        end if

        local_max = v(start)
        local_loc = start

        do i = start+1, finish
            if(v(i).gt.local_max) then
                local_max = v(i)
                local_loc = i
            end if
        end do

        maxs(thrdID+1) = local_max
        locs(thrdID+1) = local_loc
        !$omp end parallel

        maximum = maxs(1)
        location = locs(1)
        do i = 1, nthrds
            if(maxs(i).gt.maximum)then
                maximum = maxs(i)
                location = locs(i)
            end if
        end do

    end subroutine Vector_MaxValPos_Real

    subroutine Vector_Merge_Ordered_Pair_Real_Integer(A_1,A_2,p,q,r)

        real(dp), dimension(:), intent(inout) :: A_1
        integer(sp), dimension(:), intent(inout) :: A_2
        integer(sp), intent(in) :: p, q, r

        integer(sp) :: n1, n2
        real(dp), dimension(:), allocatable :: left_1, right_1, left_2, right_2

        integer(sp) :: i, j, k

        n1 = q - p + 1
        n2 = r - q

        allocate(left_1(n1+1), right_1(n2+1), left_2(n1+1), right_2(n2+1))

        do i = 1, n1
            left_1(i) = A_1(p +i-1)
            left_2(i) = A_2(p +i-1)
        end do

        do j = 1, n2
            right_1(j) = A_1(q+j)
            right_2(j) = A_2(q+j)
        end do

        left_1(n1+1) = huge(left_1(n1+1))
        left_2(n1+1) = huge(left_2(n1+1))
        right_1(n2+1) = huge(right_1(n2+1))
        right_2(n2+1) = huge(right_2(n2+1))

        i = 1
        j = 1

        do k = p, r
            if(left_1(i).le.right_1(j))then
                A_1(k) = left_1(i)
                A_2(k) = left_2(i)
                i = i + 1
            else
                A_1(k) = right_1(j)
                A_2(k) = right_2(j)
                j = j + 1
            end if
        end do

    end subroutine Vector_Merge_Ordered_Pair_Real_Integer

    recursive subroutine Vector_Merge_Sort_Ordered_Pair_Real_Integer(A_1,A_2,p,r)
        ! sorts with respect to array A_1.
        real(dp), dimension(:), intent(inout) :: A_1
        integer(sp), dimension(:), intent(inout) :: A_2
        integer(sp), intent(in) :: p, r

        integer(sp) :: q
        if(p.lt.r)then
            q = floor((real(p) + real(r))/2.d0)
            call Vector_Merge_Sort_Ordered_Pair_Real_Integer(A_1,A_2,p,q)
            call Vector_Merge_Sort_Ordered_Pair_Real_Integer(A_1,A_2,q+1,r)
            call Vector_Merge_Ordered_Pair_Real_Integer(A_1,A_2,p,q,r)

        end if

    end subroutine Vector_Merge_Sort_Ordered_Pair_Real_Integer

    subroutine Sparse_Write(A, row1, column1, row2, column2)

    ! Write CSR type to the terminal.

    type(CSR(dp)), intent(in) :: A
    integer(sp), intent(in) :: row1, row2, column1, column2

    integer(sp) :: i, j, k
    logical :: empty

    do i = row1, row2

        do k = column1, column2

            empty = .true.

            if((A%row_start(i+1)-A%row_start(i)).eq.0)then

                    write(*,'(A11,I1,A11)', advance='no')"", 0,""
                    empty = .false.

            else

                do j = A%row_start(i), A%row_start(i + 1) - 1

                        if(A%col_index(j).eq.k)then

                            write(*,'(A1,f8.3,A3,f8.3,A1,A2)',&
                                    & advance='no')"", real(A%values(j)), &
                                            & " +", aimag(A%values(j)), "i", ""

                            empty = .false.

                            exit

                        end if

                end do

            end if

            if(empty)then

                write(*,'(A11,I1,A11)', advance='no')"", 0,""

            end if

        end do

        write(*,*)

    end do

    end subroutine Sparse_Write

    subroutine Vec_to_Binary(v, filename)

        complex(dp), dimension(:), intent(in) :: v
        character(len=*), intent(in) ::filename

        integer(sp) :: i

        open( 12 , file = filename, status = "replace", &
             & access="stream", form = "unformatted" )

        write(12) size(v)

        do i = 1, size(v)
          write(12) v(i)
        end do

        close(12)

    end subroutine Vec_to_Binary


    ! Functions

    function Sparse(dense_in)

        !Populate CSR derived type.

        type(CSR(dp)) :: Sparse
        complex(8), dimension(:,:), intent(in) :: dense_in

        logical :: newline
        integer(sp) :: nz, ia_length, a_counter, ia_counter
        integer(sp) :: rows, COLUMNS
        integer(sp) :: i, j

        nz = 0
        ia_length = 1
        a_counter = 1
        ia_counter = 1

        ROWS = size(dense_in,1)
        COLUMNS = size(dense_in,2)


        do i = 1, rows
            newline = .true.
            do j = 1,COLUMNS
                if(dense_in(i,j).ne.0) then
                    nz = nz + 1
                        if(newline) then
                                ia_length = ia_length + 1
                                newline = .false.
                        end if
                end if
            end do
            if(newline) then
                ia_length = ia_length + 1
            end if
        end do

        if(allocated(Sparse%row_start).eqv..true.)then

            call Sparse_Deallocate(Sparse)

        end if

        call Sparse_Allocate(Sparse, nz, ROWS, COLUMNS)

        Sparse%values = 0
        Sparse%row_start = 0
        Sparse%col_index = 0

        do i = 1, rows
            newline = .true.
            do j = 1, COLUMNS
                if(dense_in(i,j).ne.0) then
                    Sparse%values(a_counter) = dense_in(i,j)
                    Sparse%col_index(a_counter) = j
                    if(newline) then
                        Sparse%row_start(ia_counter) = a_counter
                        ia_counter = ia_counter + 1
                        newline = .false.
                    end if
                a_counter = a_counter + 1
                end if
            end do
            if(newline) then
                Sparse%row_start(ia_counter) = a_counter
                ia_counter = ia_counter + 1
            end if
        end do

        Sparse%row_start(ia_counter) = nz + 1
        Sparse%ROWS = rows
        Sparse%COLUMNS = COLUMNS

    end function Sparse

    function Sparse_VecMul(A,v)

        ! Computes "A.v = B". Where v is a dense array of length A%COLUMNS.

        type(CSR(dp)), intent(in) :: A
        complex(8), dimension(A%ROWS) :: Sparse_VecMul

        complex(8), dimension(:), intent(in) :: v

        integer(sp) :: i, j

        Sparse_VecMul = 0

        !$omp parallel do schedule(guided)
        do i = 1, A%ROWS
            do j = A%row_start(i), A%row_start(i+1) - 1
                Sparse_VecMul(i) = Sparse_VecMul(i) + &
                        & A%values(j)*v(A%col_index(j))
            end do
        end do
        !$omp end parallel do

    end function Sparse_VecMul

    function Sparse_Size(A, dims)

        integer(sp) :: Sparse_Size
        type(CSR(dp)), intent(in) :: A
        integer(sp), intent(in) :: dims

        if(dims.eq.1) then

            Sparse_Size = A%ROWS

        else if(dims.eq.2) then

            Sparse_Size = A%COLUMNS

        end if

    end function Sparse_Size

    function Sparse_Nonzero(A)

        integer(sp) :: Sparse_Nonzero
        type(CSR(dp)), intent(in) :: A

        Sparse_Nonzero = size(A%col_index)

    end function Sparse_Nonzero

    function Sparse_ScaMul_Real(x,A)

        type(CSR(dp)) :: Sparse_ScaMul_Real
        real(dp), intent(in) :: x
        type(CSR(dp)), intent(in) :: A

        integer(sp) :: i, length

        length = size(A%values)

        call Sparse_Allocate(Sparse_ScaMul_Real, size(A%col_index), A%rows, A%columns)

        Sparse_ScaMul_Real%values = A%values
        Sparse_ScaMul_Real%col_index = A%col_index
        Sparse_ScaMul_Real%row_start = A%row_start

        Sparse_ScaMul_Real%values = x*Sparse_ScaMul_Real%values

    end function Sparse_ScaMul_Real

    function Vector_InfNorm(v)

        real(dp) :: Vector_InfNorm
        complex(8), dimension(:), intent(in) :: v

        real(dp) :: temp
        integer(sp) :: i

        Vector_InfNorm = 0

        do i = 1, size(v)
            temp = abs(v(i))
            if(temp.gt.Vector_InfNorm)then
                Vector_InfNorm = temp
            end if
        end do

    end function Vector_InfNorm

    function SparseDense_MatMul(A,B)

        complex(8), dimension(:,:), allocatable :: SparseDense_MatMul
        type(CSR(dp)), intent(in) :: A
        complex(8), dimension(:,:), intent(in) ::B

        integer(sp) :: i, j ,k

        allocate(SparseDense_MatMul(A%ROWS, size(B,2)))

        SparseDense_MatMul = 0

        !$omp parallel do schedule(guided)
        do i = 1, A%ROWS
            do j = 1, size(B,2)
                do k = A%row_start(i), A%row_start(i+1) - 1
                    SparseDense_MatMul(i,j) = SparseDense_MatMul(i,j) + &
                            & A%values(k)*B(A%col_index(k),j)
                end do
            end do
        end do
        !$omp end parallel do

    end function SparseDense_Matmul

    function Sparse_Conjugate(A)

        type(CSR(dp)) :: Sparse_Conjugate
        type(CSR(dp)) :: A
        integer(sp) :: i

        call Sparse_Allocate(Sparse_Conjugate, size(A%col_index), A%rows, A%columns)

        Sparse_Conjugate%values = A%values
        Sparse_Conjugate%col_index = A%col_index
        Sparse_Conjugate%row_start = A%row_start

        Do i = 1, size(A%col_index)
            Sparse_Conjugate%values(i) = conjg(Sparse_Conjugate%values(i))
        end do

    end function

    function Sparse_OneNormEst(A, A_Dagger, n,l)

        real(dp) :: Sparse_OneNormEst
        type(CSR(dp)), intent(in) :: A, A_Dagger
        integer(sp), intent(in) :: n, l

        real(dp) :: est, est_old, Z_norm_max, Z_norm_temp
        real(dp), dimension(:), allocatable :: Y_norm_array, Z_norm_array
        complex(8), dimension(:), allocatable :: w
        complex(8), dimension(:,:), allocatable :: Y, S, S_Old, X, Z
        logical, dimension(:), allocatable :: ind_i_hist
        logical :: not_sorted, all_done
        integer(sp), dimension(:), allocatable :: ind_l, ind_i, v
        integer(sp) :: ind_i_temp, ind_best, A_rows, itmax
        integer(sp) :: i, j, k

        A_rows = Sparse_Size(A,1)
        est_old = 0
        itmax = 20

        allocate(Y_norm_array(l), Z_norm_array(A_rows), w(A_rows), &
                & Y(A_rows,l), S(A_rows,l), S_old(A_rows, l), X(A_rows,l), &
                        & Z(A_rows,l), ind_i_hist(A_rows), ind_l(l), &
                                & ind_i(A_rows))

        ind_i_hist = .true.


        X(:,1) = 1.d0/real(A_rows,8)

        !$omp parallel do
        do i = 2, l
            do j = 1, A_rows
                X(j,i) = real(Pick_Random_Int((/-1,1/)),8)/real(A_rows,8)
            end do
        end do
        !$omp end parallel do

        k = 1

        do

            Y = A.dot.X

            do i = 1, n-1
                Y= A.dot.Y
            end do

            Y_norm_array = 0

            !Y 1-norm
            !$omp parallel do schedule(guided)
            do i = 1, l
                do j = 1, A_rows
                    Y_norm_array(i) = Y_norm_array(i) + abs(Y(j,i))
                end do
            end do
            !$omp end parallel do

            call Vector_MaxValPos_Real(Y_norm_array,est,ind_best)

            !$omp parallel do schedule(guided)
                do i = 1, A_rows
                    w(i) = Y(i,ind_best)
                end do
            !$omp end parallel do

            if((k.gt.2).and.(est.lt.est_old)) then
                est = est_old
                exit
            end if

            est_old = est
            S_old = S

            if(k.gt.itmax) exit

            !$omp parallel do schedule(guided)
            do i = 1, l
                do j = 1, A_rows
                    if(Y(j,i).eq.0) then
                        S(j,i) = 1
                    else
                        S(j,i) = Y(j,i)/abs(Y(j,i))
                    end if
                end do
            end do
            !$omp end parallel do

            Z = A_Dagger.dot.S

            do i = 1, n-1
                Z = A_Dagger.dot.Z
            end do

            Z_norm_array = 0

            !$omp parallel do
            do i = 1, A_rows
                do j = 1, l
                    Z_norm_array(i) = Z_norm_array(i) + abs(Z(i,l))
                end do
                ind_i(i) = i
            end do
            !$omp end parallel do

            Z_norm_max = Vector_MaxVal(Z_norm_array)

            if((k.gt.2).and.(Z_norm_max.eq.Z_norm_array(ind_best))) then
                exit
            end if

            call Vector_Merge_Sort_Ordered_Pair_Real_Integer(Z_norm_array,ind_i,1, &
                    & A_rows)

            if(k.gt.1) then
                all_done = .true.
                !$omp parallel do schedule(guided)
                do i = 1, l
                        if(ind_i_hist(ind_l(i))) then
                            all_done = .false.
                            ind_i_hist(ind_l(i)) = .false.
                        end if
                end do
                !$omp end parallel do
                if (all_done) exit
                ind_l = 0

                do i = 1, l
                    do j = 1, A_rows
                        if(ind_i_hist(j))then
                            ind_l(i) = j
                            exit
                        end if
                    end do
                end do
            else

            do i = 1, l
                ind_l(i) = i
            end do

            end if

            X = 0
            do i = 1, l
                X(ind_l(i), i) = (1,0)
            end do

            k = k + 1

        end do

        Sparse_OneNormEst = est

    end function Sparse_OneNormEst

    function Sparse_Transpose(A)

        type(CSR(dp)) :: Sparse_Transpose
        type(CSR(dp)), intent(in) :: A

        integer(sp) :: i, j, NZ, loc
        integer(sp), dimension(:), allocatable :: intra, csrRowIdx
        integer(sp), dimension(:,:), allocatable :: inter

        integer(sp) :: start, finish, thrdID, nthrds
        integer(sp), dimension(:), allocatable :: prefix_intermediate

        NZ = Sparse_Nonzero(A)

        allocate(intra(NZ),csrRowIdx(NZ),inter(omp_get_max_threads()+1,A%COLUMNS))

        call Sparse_Allocate(Sparse_Transpose,NZ,A%COLUMNS,A%ROWS)

        Sparse_Transpose%ROWS = A%COLUMNS
        Sparse_Transpose%COLUMNS = A%ROWS
        Sparse_Transpose%row_start = 0
        Sparse_Transpose%row_start(1) = 1

        intra = 0
        inter = 0

        nthrds = omp_get_max_threads()

        ! Parallel reigon #1

        !$omp parallel do schedule(guided)
        do i = 1, A%ROWS
            do j = A%row_start(i), A%row_start(i+1)-1
                csrRowIdx(j) = i
            end do
        end do
        !$omp end parallel do

        start = 1
        finish = NZ
        thrdID = 0

        ! Parallel reigon #2

        !$omp parallel private(thrdID,i,start,finish)
        thrdID = omp_get_thread_num()
        start = (thrdID*NZ/nthrds) + 1
        finish = (thrdID+1)*NZ/nthrds

        if((thrdID+1).eq.nthrds) then
            finish = NZ
        end if

        do i = start, finish
            inter(thrdID+2,A%col_index(i)) = inter(thrdID+2,A%col_index(i)) + 1
            if((inter(thrdID+2,A%col_index(i))).gt.1) then
                intra(i) = inter(thrdID+2,A%col_index(i)) - 1
            end if
        end do
        !$omp end parallel

        ! Parallel reigon #3

        !$omp parallel do schedule(guided)
        do i = 1, A%COLUMNS
            do j = 2, nthrds + 1
                inter(j,i) = inter(j,i) + inter(j-1,i)
            end do
        end do
        !$omp end parallel do
        ! Parallel reigon #4
        !$omp parallel do schedule(guided)
        do i = 1, A%COLUMNS
            Sparse_Transpose%row_start(i+1) = inter(nthrds + 1, i)
        end do
        !$omp end parallel do

        !parallel prefix sum

        start = 1
        finish = A%COLUMNS

        allocate(prefix_intermediate(omp_get_max_threads()))
        prefix_intermediate = 0

        !$omp parallel private(thrdID,i,start,finish) shared (prefix_intermediate)
        thrdID = omp_get_thread_num()
        start = thrdID*(A%COLUMNS+1)/nthrds + 1
        finish = (thrdID+1)*(A%COLUMNS+1)/nthrds

        if((thrdID+1).eq.nthrds) then
            finish = (A%COLUMNS+1)
        end if

        do i = start, finish-1
                Sparse_Transpose%row_start(i+1) = Sparse_Transpose%row_start(i+1) + Sparse_Transpose%row_start(i)
        end do

        if((thrdID+1).ne.nthrds) then
            prefix_intermediate(thrdID+2) = Sparse_Transpose%row_start(finish)
        end if

        !$omp barrier

        !$omp single
        do i = 2, nthrds
            prefix_intermediate(i) = prefix_intermediate(i-1) + prefix_intermediate(i)
        end do
        !$omp end single

        do i = start, finish
            Sparse_Transpose%row_start(i) = Sparse_Transpose%row_start(i) + &
                    prefix_intermediate(thrdID+1)
        end do
        !$omp end parallel

        !end parallel prefix sum

        ! Parallel reigon #5
        !$omp parallel private(thrdID,i,start,finish,loc)
        thrdID = omp_get_thread_num()
        start = thrdID*NZ/nthrds + 1
        finish = (thrdID+1)*NZ/nthrds

        if((thrdID+1).eq.nthrds) then
            finish = NZ
        end if

        do i = start, finish
            loc = Sparse_Transpose%row_start(A%col_index(i)) + inter(thrdID+1,A%col_index(i)) + intra(i)
            Sparse_Transpose%col_index(loc) = csrRowIdx(i)
            Sparse_Transpose%values(loc) = A%values(i)
        end do
        !$omp end parallel

    end function Sparse_Transpose

    function Vector_MaxVal(v)

        real(dp) :: Vector_MaxVal
        real(dp), dimension(:), intent(in) :: v

        real(dp), dimension(:), allocatable :: maxs
        real(dp) :: local_max
        integer(sp) :: i, start, finish, thrdID, length

        length = size(v)

        allocate(maxs(omp_get_max_threads()))

        !$omp parallel private(thrdID,i,start,finish,local_max) shared(maxs)
        thrdID = omp_get_thread_num()
        start = (thrdID*length/omp_get_num_threads()) + 1
        finish = (thrdID+1)*length/omp_get_num_threads()
        if((thrdID+1).eq.omp_get_num_threads()) then
            finish = length
        end if

        local_max = v(start)

        do i = start+1, finish
            if(v(i).gt.local_max) then
                local_max = v(i)
            end if
        end do

        maxs(thrdID+1) = local_max
        !$omp end parallel

        Vector_MaxVal = maxs(1)

        do i = 1, size(maxs)
            if(maxs(i).gt.Vector_MaxVal)then
                Vector_MaxVal = maxs(i)
            end if
        end do

    end function Vector_MaxVal

    function Pick_Random_Int(intlist)

        integer(sp) :: Pick_Random_Int
        integer(sp), dimension(:), intent(in) :: intlist

        real :: r

            call random_seed
            call random_number(r)
            Pick_Random_Int = intlist(int(r*size(intlist))+1)

    end function Pick_Random_Int

end module Nullarbor


