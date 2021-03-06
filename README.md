# Nullarbor_and_Spinifex
MPI + OpenMP sparse matrix libraries. 

# Overview #

Nullarbor.f90 and Spinifex.f90 contain modules designed for working with [compressed sparse row matrices](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)) in an MPI + OpenMP environment. They were developed as part of my coursework in semester 2 2018. </P>

The primary intent of this repository is to demonstrate use of MPI-I/O and a repeated matrix-dense-vector multiplication method to calculate <b>A</b><sup>n</sup>.<b>u</b> = <b>v</b>, where <b>A</b> is a sparse matrix, n is an integer and <b>u</b>,<b>v</b> are dense vectors. The approach used is appropriate for a scenario where MPI communication between networked computational nodes is a primary cost factor. The subroutines of interest are Spinifex_Import_Vector, Spinifex_Sparse_Read_and_Distribute and Spinifex_Sparse_SpMV. These are, unsurprisingly, contained within Spinifex.f90 with the process briefly detail by the code comments.

# Build Instructions #

Currently these modules require a Fortran compiler which has implemented support for parameterized types (part of the Fortran 2008 standard). The Intel Fortran compiler fits this criteria. </P>

The included makefile compiles Example.f90. To run simply type mpiexec -n \<number of nodes> ./Example FMO.bin statevector.bin \<output file name>.
