#define STDOUT 0xd0580000

.section .text
.global _start
_start:
## START YOUR CODE HERE

la a0, matrix
lw a1, size
call printToLogVectorized

la a0, matrix
lw a1, size
call transpose

la a0, matrix
lw a1, size
call printToLogVectorized


j _finish

transpose:
    # a0 = base address, a1 = n (matrix dimension)

    # Initialize outer loop index: i = 0
    addi    s0, zero, 0         # s0 holds i

outer_loop:
    bge     s0, a1, done         # if i >= n, we are done
    addi    s1, s0, 1           # initialize inner loop index: j = i + 1

inner_loop:
    bge     s1, a1, next_i       # if j >= n, advance to next row

    # Compute remaining elements count: t1 = n - j
    sub     t1, a1, s1

    # Set vector length (VL) for 32-bit floats; t0 = actual vl = min(n-j, VLEN)
    vsetvli t0, t1, e32, m1

    # --- Compute address for matrix[i][j] -----------------------
    # Element (i,j) address = a0 + ((i*n + j) * 4)
    mul     t2, s0, a1         # t2 = i * n
    slli    t2, t2, 2          # t2 = i * n * 4
    add     t3, a0, t2         # t3 = base address of row i
    slli    t4, s1, 2          # t4 = j * 4
    add     t5, t3, t4         # t5 = address of matrix[i][j]

    # Load vector of row elements starting at matrix[i][j]
    vle32.v v0, (t5)           # v0 <- matrix[i][j ... j+vl-1]

    # --- Compute index vector for loading column elements ---
    # For each k in 0..vl-1, we want the address for matrix[j+k][i]:
    # Address = a0 + [ (j+k)*n*4 + i*4 ]
    slli    t6, a1, 2          # t6 = n * 4 (stride for one row)
    slli    s7, s0, 2          # s7 = i * 4 (column offset)

    # Create a vector of indices: v2 = {0,1,...,vl-1}
    vid.v     v2                 # v2[i] = i for i=0...vl-1
    vadd.vx v2, v2, s1         # v2 = {j, j+1, ..., j+vl-1}
    vmul.vx v2, v2, t6         # v2 = {j*n*4, (j+1)*n*4, ...}
    vadd.vx v2, v2, s7         # v2 = {j*n*4 + i*4, (j+1)*n*4 + i*4, ...}

    # Use indexed load to fetch column elements: matrix[j...][i]
    vluxei32.v v1, (a0), v2       # v1 <- { matrix[j][i], matrix[j+1][i], ... }

    # --- Swap the two sets of elements -----------------------
    vse32.v v1, (t5)           # Store v1 into row i at matrix[i][j ...]
    vsuxei32.v v0, (a0), v2       # Store v0 into column i at matrix[j...][i]

    # Update inner loop index: j = j + vl (t0 holds vl)
    add     s8, s1, t0
    mv      s1, s8
    j       inner_loop

next_i:
    addi    s0, s0, 1          # i++
    j       outer_loop

done:
    ret
## END YOU CODE HERE

# Function: print
# Logs values from array in a0 into registers v1 for debugging and output.
# Inputs:
#   - a0: Base address of array
#   - a1: Size of array i.e. number of elements to log
# Clobbers: t0,t1, t2,t3 ft0, ft1.
printToLogVectorized:        
    addi sp, sp, -4
    sw a0, 0(sp)

    li t0, 0x123                 # Pattern for help in python script
    li t0, 0x456                 # Pattern for help in python script
    mv a1, a1                   # moving size to get it from log 
    mul a1, a1, a1              # sqaure matrix has n^2 elements 
	li t0, 0		                # load i = 0
    printloop:
        vsetvli t3, a1, e32           # Set VLEN based on a1
        slli t4, t3, 2                # Compute VLEN * 4 for address increment

        vle32.v v1, (a0)              # Load real[i] into v1
        add a0, a0, t4                # Increment pointer for real[] by VLEN * 4
        add t0, t0, t3                # Increment index

        bge t0, a1, endPrintLoop      # Exit loop if i >= size
        j printloop                   # Jump to start of loop
    endPrintLoop:
    li t0, 0x123                    # Pattern for help in python script
    li t0, 0x456                    # Pattern for help in python script
	
    lw a0, 0(sp)
    addi sp, sp, 4

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
.equ MatrixSize, 20
matrix:
    .float 333.75, 588.0, -957.75, -914.75, -150.0, -966.75, 744.25, 876.0, 83.25, -302.25, -49.0, 588.0, -6.25, 36.25, 118.5, 367.25, 343.25, 439.75, 20.5, 8.5
    .float 700.75, 141.5, -464.25, 698.25, 370.75, 835.75, 492.25, 424.75, 958.75, -293.25, 92.75, 716.5, -478.25, 791.5, -142.25, 761.0, 235.5, -650.0, 345.5, 946.25
    .float -311.25, 220.25, 94.25, -788.75, 485.5, -114.25, -365.25, -223.0, -250.25, 511.75, 397.0, 278.75, 139.0, -914.5, 550.5, 114.75, 418.25, 861.75, -241.25, -686.75
    .float 663.5, 77.25, 486.75, 125.75, 335.0, 918.25, 128.0, 442.5, 997.0, 32.5, -270.5, -497.75, -375.75, 803.25, 279.75, 991.5, 876.25, 929.75, 245.5, 917.0
    .float 961.25, -541.75, 465.25, 576.75, 674.25, 803.0, -180.25, 259.0, -240.75, -417.0, 453.0, 607.75, -428.25, 442.25, -616.0, -874.25, -473.0, 457.25, 491.0, 530.25
    .float 921.75, 600.0, 61.5, 7.5, -680.25, 784.75, 290.0, 482.75, -342.0, 117.75, 805.75, 349.75, -382.25, -397.0, 483.5, 757.5, -237.0, -33.0, -569.25, 767.0
    .float -658.75, -412.5, 419.0, -969.0, -289.25, 607.5, -555.25, -577.5, -819.0, 656.0, -852.75, -211.75, 291.75, -661.5, -790.25, 652.75, -320.75, -312.75, 207.5, -866.75
    .float 952.5, 155.0, -592.0, 325.0, 232.5, -455.5, -833.5, 91.25, 444.25, 4.0, 147.75, 595.25, -920.0, 856.5, 170.5, -68.25, -677.0, 914.25, 717.0, 427.75
    .float 815.25, -245.0, 41.5, 891.5, 259.0, -658.25, -6.0, -472.5, -745.75, -689.5, -730.0, 285.25, 206.75, 668.5, -279.5, 696.75, -259.0, 780.75, -296.5, 145.75
    .float 957.5, -672.5, -647.75, 699.25, -986.5, 87.5, 639.0, -343.0, -782.0, -560.75, 527.0, 516.0, -315.5, -360.5, -363.0, 825.25, 889.25, 646.0, -568.25, 518.5
    .float 53.25, -820.0, 405.0, -798.0, -674.5, 595.75, -877.75, 900.25, -878.75, 497.25, -181.0, 54.5, 100.25, 405.0, -762.0, 646.25, -492.0, 619.5, -223.75, -510.25
    .float 941.75, -474.25, 395.5, 813.0, 494.75, -182.75, 324.5, -586.5, -4.75, -784.75, 880.75, -971.75, -380.0, -645.5, -592.5, -726.75, -150.0, -405.5, -854.5, -388.5
    .float -961.75, -122.0, -734.75, -147.0, 741.25, -139.0, 630.0, 468.75, 62.5, -103.75, 875.75, 970.0, 477.75, -503.0, 663.0, -213.0, 797.0, 956.5, -739.5, -813.5
    .float -567.0, -935.0, -912.75, -779.0, 91.75, -584.5, 320.25, 494.75, 877.75, 268.5, -66.0, -589.25, 885.5, -839.75, -413.5, 411.5, 209.75, -829.25, 933.75, -698.25
    .float 159.25, 978.5, 81.75, 242.75, -63.5, -824.75, 200.25, -433.25, 622.5, 453.5, 905.25, -657.0, 180.25, 674.75, -310.5, -823.5, 767.75, 675.25, -962.5, -516.75
    .float -853.0, 930.25, -340.25, -413.5, 473.0, -597.25, -411.75, -313.0, -59.0, -295.75, -587.25, 294.5, -0.5, 244.75, -69.75, -982.25, 20.0, 27.5, -184.75, 158.0
    .float 34.5, -184.5, 491.25, -957.25, 410.5, 622.75, 338.25, 469.25, -816.0, 120.5, 845.0, -950.5, -678.75, -307.0, 667.75, -639.75, 538.25, -72.0, -555.75, -906.25
    .float -501.0, -651.75, -390.5, -88.5, 240.75, -307.75, -41.0, -16.25, 277.75, -690.25, -209.0, -793.5, 306.5, -400.0, 268.25, -83.0, -553.25, -660.0, 513.5, -167.0
    .float 464.75, 450.5, -606.5, -832.25, 679.5, -628.5, 542.5, -9.25, 591.75, -618.75, -845.25, 699.0, -23.5, 600.0, -506.0, -762.25, -472.0, 605.0, -846.75, -202.5
    .float -65.75, -333.25, 495.0, -268.0, 281.75, -817.75, 319.75, 990.75, 793.0, 640.0, -870.75, -792.5, 29.75, 122.75, -168.0, 802.25, 802.75, -269.75, 354.75, -402.5
## DATA DEFINE END
size: .word MatrixSize