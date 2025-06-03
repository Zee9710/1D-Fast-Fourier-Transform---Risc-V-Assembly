#define STDOUT 0xd0580000

.section .text
.global _start
_start:
## START YOUR CODE HERE

la a0, matrix
lw a1, size
call printToLog

la a0, matrix
lw a1, size
call transpose

la a0, matrix
lw a1, size
call printToLog


j _finish

transpose:
    # Prologue
    addi sp, sp, -16
    sw ra, 12(sp)
    sw s0, 8(sp)
    sw s1, 4(sp)
    sw s2, 0(sp)

    mv s0, a0      # s0 = base address of matrix
    mv s1, a1      # s1 = size N

    li t0, 0       # i = 0 (outer loop counter)

outer_loop:
    bge t0, s1, end_transpose  # if i >= N, exit

    addi t1, t0, 1  # j = i + 1 (inner loop counter)

inner_loop:
    bge t1, s1, next_row  # if j >= N, go to next row

    # Compute offsets: (i * N + j) and (j * N + i)
    mul t2, t0, s1  # t2 = i * N
    add t2, t2, t1  # t2 = i * N + j
    slli t2, t2, 2  # t2 *= 4 (float size)

    mul t3, t1, s1  # t3 = j * N
    add t3, t3, t0  # t3 = j * N + i
    slli t3, t3, 2  # t3 *= 4

    add t2, t2, s0  # Address of matrix[i][j]
    add t3, t3, s0  # Address of matrix[j][i]

    # Swap elements
    flw ft0, 0(t2)  # Load matrix[i][j] into ft0
    flw ft1, 0(t3)  # Load matrix[j][i] into ft1
    fsw ft1, 0(t2)  # Store matrix[j][i] into matrix[i][j]
    fsw ft0, 0(t3)  # Store matrix[i][j] into matrix[j][i]

    addi t1, t1, 1  # j++
    j inner_loop

next_row:
    addi t0, t0, 1  # i++
    j outer_loop

end_transpose:
    # Epilogue
    lw ra, 12(sp)
    lw s0, 8(sp)
    lw s1, 4(sp)
    lw s2, 0(sp)
    addi sp, sp, 16
    ret
## END YOU CODE HERE

# Function to print a matrix for debugging purposes
# This function iterates over all elements of a matrix stored in memory.
# Instead of calculating the end address in each iteration, it precomputes 
# the end address (baseAddress + size^2 * 4) to optimize the loop.
# Input:
#   a0: Base address of the matrix
#   a1: Size of matrix
# Clobbers:
#   t0, t1, ft0
printToLog:
    li t0, 0x123                #  Identifiers used for python script to read logs
    li t0, 0x456
    mv a1, a1                   # moving size to get it from log 
    mv t0, a0                   # Copy the base address of the matrix to t0 to avoid modifying a0
    mul t1, a1, a1              # size^2 
    slli  t1, t1, 2             # size^2 * 4 (total size of the matrix in bytes)
    add t1, a0, t1              # Calculate the end address (base address + total size)

    printMatrixLoop:
        bge t0, t1, printMatrixLoopEnd 
        flw ft0, 0(t0)          # Load from array
        addi t0, t0, 4          # increment address by elem size
        j printMatrixLoop
    printMatrixLoopEnd:

    li t0, 0x123                #  Identifiers used for python script to read logs
    li t0, 0x456

    jr ra


# Function: _finish
# VeeR Related function which writes to to_host which stops the simulator
_finish:
    li x3, 0xd0580000
    addi x5, x0, 0xff
    sb x5, 0(x3)
    beq x0, x0, _finish

    .rept 100
        nop
    .endr


.data
## ALL DATA IS DEFINED HERE LIKE MATRIX, CONSTANTS ETC

## DATA DEFINE START
.equ MatrixSize, 5
matrix:
    .float -10.0, 13.0, 10.0, -3.0, 2.0
    .float 6.0, 15.0, 4.0, 13.0, 4.0
    .float 18.0, 2.0, 9.0, 8.0, -4.0
    .float 5.0, 4.0, 12.0, 17.0, 6.0
    .float -10.0, 7.0, 13.0, -3.0, 16.0
## DATA DEFINE END
size: .word MatrixSize