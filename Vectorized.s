.section .text
    .global _start
    .option norvc         # Disable compressed instructions (use full 32-bit)

# Stack offset constants for local storage in butterfly_stages
.equ OFF_SAVED_REGS_SIZE, 16          # Space for saving return address and callee-saved registers

# 32 bytes = 8 floats (VLEN=8 for e32), used for each vector buffer
.equ VL_MAX_TIMES_4, 32

# Offsets for temporary scalar arrays (indices) used in vector memory operations
.equ OFF_offsets_a_real, 0           # Offsets for a_real vector access
.equ OFF_offsets_b_real, 32          # Offsets for b_real vector access
.equ OFF_offsets_a_imag, 64          # Offsets for a_imag vector access
.equ OFF_offsets_b_imag, 96          # Offsets for b_imag vector access

# Storage offsets for vector register contents during computation
.equ OFF_v0_data, 128                # Temporary store for v0 (a_real)
.equ OFF_v1_data, 160                # Temporary store for v1 (a_imag)
.equ OFF_v8_data, 192                # Temporary store for v8 (b_real)
.equ OFF_v9_data, 224                # Temporary store for v9 (b_imag)

# Temporary buffer for complex multiplication results (interleaved format)
.equ OFF_temp_cmult_real, 400        # 64 bytes if vl_max = 8 (8 complex floats)

# Output storage offsets for computed results (planar layout)
.equ OFF_new_a_real, 528             # Result: real part of new a
.equ OFF_new_a_imag, 560             # Result: imag part of new a
.equ OFF_new_b_real, 592             # Result: real part of new b
.equ OFF_new_b_imag, 624             # Result: imag part of new b

# Total space required for butterfly_stages local stack usage
.equ BUTTERFLY_STACK_LOCALS_SIZE, 736  # Minimum required size for locals
.equ BUTTERFLY_STACK_FRAME_SIZE, 896   # Full frame (locals + saved regs + alignment)


# ------------------------------------------------------------
#  RISC‑V RV32 1‑D FFT (N=16 complex samples)
# ------------------------------------------------------------
_start:
    # 1) Bit-reverse input_data → bitrev_output
    la      a0, input_data
    la      a1, bitrev_output
    li      a2, 16
    call    bit_reverse

    # --- Save bit-reversed output to file ---
    la      a0, bitrev_filename   # filename string
    la      a1, bitrev_output     # buffer address
    li      a2, 128               # bytes (16 complex × 8 bytes)
    call    write_to_file

    # 2) FFT butterfly using twiddles
    la      a0, bitrev_output
    la      a1, twiddles
    la      a2, fft_output
    li      a3, 16
    call    butterfly_stages

    # --- Save FFT output to file ---
    la      a0, fft_filename      # filename string
    la      a1, fft_output        # buffer address
    li      a2, 128               # bytes (16 complex × 8 bytes)
    call    write_to_file

    # Infinite loop (signal "done" to host)
_finish:
    lui     x3, %hi(STDOUT_ADDR)
    addi    x3, x3, %lo(STDOUT_ADDR)
    li      x5, 0xff
    sb      x5, 0(x3)
    beq     x0, x0, _finish

# ----------------------------------------------------
#  bit_reverse: a0=inPtr, a1=outPtr, a2=N
# ----------------------------------------------------
.globl bit_reverse
bit_reverse:
    li      t0, 0
    mv      t1, a2
bitrev_loop:
    beq     t0, t1, bitrev_done

    # For N=16, need a 4-bit reversal
    andi    t2, t0, 1       # bit 0
    slli    t2, t2, 3       # shift to position 3

    andi    t3, t0, 2       # bit 1
    slli    t3, t3, 1       # shift to position 2

    andi    t4, t0, 4       # bit 2
    srli    t4, t4, 1       # shift to position 1

    andi    t5, t0, 8       # bit 3
    srli    t5, t5, 3       # shift to position 0

    // combine the reversed
    or      t6, t2, t3
    or      t6, t6, t4
    or      t6, t6, t5

    slli    s7, t0, 3       # byte_offset = i * 8
    add     s7, a0, s7      # in_addr = inPtr + byte_offset
    flw     f0, 0(s7)       # load input_real
    flw     f1, 4(s7)       # load input_imag

    slli    s7, t6, 3       # byte_offset_rev = rev_i * 8
    add     s7, a1, s7      # out_addr = outPtr + byte_offset_rev
    fsw     f0, 0(s7)       # store input_real to output[rev_i]
    fsw     f1, 4(s7)       # store input_imag to output[rev_i]

    addi    t0, t0, 1
    j       bitrev_loop
bitrev_done:
    ret

# ----------------------------------------------------
#  butterfly_stages: a0=inPtr (current data), a1=twiddlePtr, a2=outPtr (same as inPtr after init_copy), a3=N
# ------------------------------------------------------------
.globl butterfly_stages
butterfly_stages:
    addi    sp, sp, -BUTTERFLY_STACK_FRAME_SIZE # Allocate stack frame
    sw      ra, (BUTTERFLY_STACK_FRAME_SIZE - 4)(sp)
    sw      s0, (BUTTERFLY_STACK_FRAME_SIZE - 8)(sp)
    sw      s1, (BUTTERFLY_STACK_FRAME_SIZE - 12)(sp)
    sw      s2, (BUTTERFLY_STACK_FRAME_SIZE - 16)(sp)
    # s3-s11 can be used as callee-saved if needed. We use s3, s4 in loops.
    # Ensure enough saved register space if more are used heavily across calls.
    # For now, assuming s0-s2 are primary callee-saved for this function's context.
    # Add s3, s4 to saved registers if they need to persist across calls made by this func.
    # However, current scalar loops are local and reuse s3, s4 without calls.

    # Copy input to output for first stage preparation (if inPtr and outPtr are different)
    # Current code implies a0 (bitrev_output) is copied to a2 (fft_output)
    # This is an in-place FFT, so a0 and a2 effectively point to the same buffer after this copy.
    li      t0, 0           # loop counter i
init_copy_loop:
    beq     t0, a3, init_copy_done # if i == N, done
    slli    t1, t0, 3       # t1 = i * 8 (byte offset for complex float)
    add     t2, a0, t1      # t2 = &inPtr[i]
    add     t3, a2, t1      # t3 = &outPtr[i] (fft_output)
    
    flw     f0, 0(t2)       # load real part from inPtr
    flw     f1, 4(t2)       # load imag part from inPtr
    fsw     f0, 0(t3)       # store real part to outPtr
    fsw     f1, 4(t3)       # store imag part to outPtr
    
    addi    t0, t0, 1       # i++
    j       init_copy_loop
init_copy_done:

    # Base address for FFT data is now a2 (fft_output)
    # Process FFT stages (len = 2, 4, 8, 16)
    li      t0, 2           # len = 2 (current stage length)
stage_loop:
    bgt     t0, a3, stage_done # if len > N, FFT done
    srli    t1, t0, 1       # half = len / 2
    
    # Process blocks within this stage
    li      t2, 0           # block_start = 0
block_loop:
    bge     t2, a3, block_done # if block_start >= N, done with this stage
    
    # Vector processing for butterflies within this block
    li      t3, 0           # j = 0 (index within a half-block)
butterfly_vector_loop:
    sub     t4, t1, t3      # remaining_in_half_block = half - j
    beqz    t4, butterfly_block_processed # if no more elements in this half_block pass, done

    vsetvli t5, t4, e32, m1 # t5 = vl (actual vector length for this iteration)
    
    # Generate indices for the first part of each butterfly: idx1 = block_start + j + {0..vl-1}
    vid.v   v2              # v2 = {0, 1, ..., vl-1}
    mv      s0, t3          # s0 = j (use a saved reg or reload for inner loops)
    vadd.vx v2, v2, s0      # v2 = {j, j+1, ..., j+vl-1}
    mv      s0, t2          # s0 = block_start
    vadd.vx v2, v2, s0      # v2 = {block_start+j, ..., block_start+j+vl-1} (these are idx1 values)
    
    # Generate indices for the second part: idx2 = idx1 + half
    vmv.v.v v3, v2          # v3 = v2 (copy of idx1 values)
    mv      s0, t1          # s0 = half
    vadd.vx v3, v3, s0      # v3 = {idx1+half, ..., idx1+half+vl-1} (these are idx2 values)

    # --- Scalar loop to prepare byte offsets for vector gather/scatter ---
    # v4 will store byte_offsets for idx1 elements (real)
    # v5 will store byte_offsets for idx2 elements (real)
    # v6 will store byte_offsets for imag_part of idx1 elements
    # v7 will store byte_offsets for imag_part of idx2 elements
    li      s0, 0           # loop_iter = 0 (re-using s0, careful if it held other value)
offset_loop_scalar:
    beq     s0, t5, offset_loop_scalar_done # if loop_iter == vl, done

    add     s1, t2, t3      # s1 = block_start + j
    add     s1, s1, s0      # s1 = idx1_scalar = block_start + j + loop_iter
    
    add     s2, s1, t1      # s2 = idx2_scalar = idx1_scalar + half

    slli    s4, s1, 3       # byte_offset_idx1_real = idx1_scalar * 8
    slli    s5, s2, 3       # byte_offset_idx2_real = idx2_scalar * 8

    add     s6, sp, OFF_offsets_a_real
    slli    s7, s0, 2       # offset_within_array = loop_iter * 4 (for float/int offsets)
    add     s6, s6, s7
    sw      s4, 0(s6)       # Store byte_offset_idx1_real

    add     s6, sp, OFF_offsets_b_real
    add     s6, s6, s7
    sw      s5, 0(s6)       # Store byte_offset_idx2_real

    addi    s4, s4, 4       # byte_offset_idx1_imag
    addi    s5, s5, 4       # byte_offset_idx2_imag

    add     s6, sp, OFF_offsets_a_imag
    add     s6, s6, s7
    sw      s4, 0(s6)       # Store byte_offset_idx1_imag

    add     s6, sp, OFF_offsets_b_imag
    add     s6, s6, s7
    sw      s5, 0(s6)       # Store byte_offset_idx2_imag

    addi    s0, s0, 1       # loop_iter++
    j       offset_loop_scalar
offset_loop_scalar_done:

    add     s0, sp, OFF_offsets_a_real
    vle32.v v4, (s0)
    add     s0, sp, OFF_offsets_b_real
    vle32.v v5, (s0)
    add     s0, sp, OFF_offsets_a_imag
    vle32.v v6, (s0)
    add     s0, sp, OFF_offsets_b_imag
    vle32.v v7, (s0)

    vluxei32.v v0, (a2), v4
    vluxei32.v v1, (a2), v6
    vluxei32.v v8, (a2), v5
    vluxei32.v v9, (a2), v7

    add     s0, sp, OFF_v0_data
    vse32.v v0, (s0)
    add     s0, sp, OFF_v1_data
    vse32.v v1, (s0)
    add     s0, sp, OFF_v8_data
    vse32.v v8, (s0)
    add     s0, sp, OFF_v9_data
    vse32.v v9, (s0)

    li      s0, 0           # loop_iter = 0
twiddle_multiply_loop:
    beq     s0, t5, twiddle_multiply_done

    add     s1, t3, s0                      # s1 = j + loop_iter = twiddle index
    div     s2, a3, t0                      # s2 = N / len (length of sub-DFT block)
    mul     s1, s1, s2                      # twiddle_index = (j + i) * (N / len)
    
    # Compute the twiddle index from j, loop_iter, N, and len.



slli    s1, s1, 3       # byte_offset = twiddle_index * 8
add     s1, a1, s1      # s1 = address of twiddles[twiddle_index]
flw     f4, 0(s1)       # load twiddle real part into f4
flw     f5, 4(s1)       # load twiddle imaginary part into f5


add     s1, sp, OFF_v8_data   # s1 = address of base of vector v8_data in stack OFF_v8_data--------> b's real number data 
slli    s2, s0, 2            # s2 = loop_iter * 4 (byte offset for float)
add     s1, s1, s2           # s1 = address of element v8_data[loop_iter]
flw     f8, 0(s1)            # load float value from v8_data[loop_iter] into f8

add     s1, sp, OFF_v9_data  # s1 = address of base of vector v9_data in stack  OFF_v9_data---------> b's imaginary number data
add     s1, s1, s2           # s1 = address of element v9_data[loop_iter]
flw     f9, 0(s1)            # load float value from v9_data[loop_iter] into f9


fmul.s  f10, f8, f4    # f10 = b.real * a.real  (f8 * f4)
fmul.s  f11, f9, f5    # f11 = b.imag * a.imag  (f9 * f5)
fsub.s  f6, f10, f11   # f6 = (b.real * a.real) - (b.imag * a.imag)  => c.real


fmul.s  f10, f8, f5    # f10 = b.real * a.imag  (f8 * f5)
fmul.s  f11, f9, f4    # f11 = b.imag * a.real  (f9 * f4)
fadd.s  f7, f10, f11   # f7 = (b.real * a.imag) + (b.imag * a.real) => c.imag


    add     s1, sp, OFF_temp_cmult_real   # Load base address of temp storage for complex results into s1
    slli    s2, s0, 3                     # Calculate byte offset: s0 * 8 (each complex number = 8 bytes)
    add     s1, s1, s2                    # Add offset to base address to get current element's address

    fsw     f6, 0(s1)                    # Store real part (f6) of complex multiplication at address s1
    fsw     f7, 4(s1)                    # Store imaginary part (f7) at address s1 + 4 bytes (right after real part)


    addi    s0, s0, 1
    j       twiddle_multiply_loop
twiddle_multiply_done:







   li      s0, 0           # Initialize loop counter s0 = 0

butterfly_op_scalar_loop:
    beq     s0, t5, butterfly_op_scalar_done  # If s0 == t5 (vector length), exit loop

    # Load first complex vector element's real part
    add     s1, sp, OFF_v0_data              # Load base address of v0_data into s1
    slli    s2, s0, 2                        # Calculate offset = s0 * 4 bytes (each float is 4 bytes)
    add     s1, s1, s2                       # Add offset to base to get address of element s0
    flw     f0, 0(s1)                        # Load float from memory into f0 (real part of vector 0)

    # Load first complex vector element's imaginary part
    add     s1, sp, OFF_v1_data              # Load base address of v1_data into s1
    add     s1, s1, s2                       # Add same offset to get element s0 address
    flw     f1, 0(s1)                        # Load float from memory into f1 (imaginary part of vector 0)

    # Load twiddle multiplied complex value (result of previous multiplication)
    add     s1, sp, OFF_temp_cmult_real      # Base address of temporary complex multiplied results
    slli    s3, s0, 3                        # Offset = s0 * 8 bytes (each complex number has real + imag, 8 bytes total)
    add     s1, s1, s3                       # Add offset for element s0
    flw     f10, 0(s1)                       # Load real part of twiddle multiplied complex number
    flw     f11, 4(s1)                       # Load imaginary part of twiddle multiplied complex number

    # Perform butterfly operations:
    fadd.s  f6, f0, f10                      # real part: f6 = f0 + f10
    fadd.s  f7, f1, f11                      # imag part: f7 = f1 + f11

    fsub.s  f8, f0, f10                      # real part: f8 = f0 - f10
    fsub.s  f9, f1, f11                      # imag part: f9 = f1 - f11



    # --- CORRECTED SECTION: Store new_a/b real/imag to PLANAR stack arrays ---
    # s2 (loop_iter * 4) is already calculated and suitable for float offset in planar array
    
    # Store new_a_real[loop_iter]
    add     s3, sp, OFF_new_a_real  # s3 now holds base address for this store
    add     s3, s3, s2              # Address = sp + OFF_new_a_real + loop_iter*4
    fsw     f6, 0(s3)               # f6 holds new_a_real

    # Store new_a_imag[loop_iter]
    add     s3, sp, OFF_new_a_imag
    add     s3, s3, s2              # Address = sp + OFF_new_a_imag + loop_iter*4
    fsw     f7, 0(s3)               # f7 holds new_a_imag

    # Store new_b_real[loop_iter]
    add     s3, sp, OFF_new_b_real
    add     s3, s3, s2              # Address = sp + OFF_new_b_real + loop_iter*4
    fsw     f8, 0(s3)               # f8 holds new_b_real

    # Store new_b_imag[loop_iter]
    add     s3, sp, OFF_new_b_imag
    add     s3, s3, s2              # Address = sp + OFF_new_b_imag + loop_iter*4
    fsw     f9, 0(s3)               # f9 holds new_b_imag
    # --- END OF CORRECTED SECTION ---

    addi    s0, s0, 1
    j       butterfly_op_scalar_loop






butterfly_op_scalar_done:

    # Load computed new_a and new_b values from PLANAR stack arrays into vector registers
    # This part of the code remains unchanged, but will now load correct planar data
    add     s0, sp, OFF_new_a_real
    vle32.v v12, (s0)
    add     s0, sp, OFF_new_a_imag
    vle32.v v13, (s0)

    add     s0, sp, OFF_new_b_real
    vle32.v v14, (s0)
    add     s0, sp, OFF_new_b_imag
    vle32.v v15, (s0)

    # Store results back to main buffer a2 using original byte_offset vectors v4,v5,v6,v7

vsuxei32.v v12, (a2), v4   # store new_a_real at offsets in v4
vsuxei32.v v13, (a2), v6   # store new_a_imag at offsets in v6
vsuxei32.v v14, (a2), v5   # store new_b_real at offsets in v5
vsuxei32.v v15, (a2), v7   # store new_b_imag at offsets in v7

    
    add     t3, t3, t5
    j       butterfly_vector_loop
    
butterfly_block_processed:
    add     t2, t2, t0
    j       block_loop
    
block_done:
    slli    t0, t0, 1
    j       stage_loop
    
stage_done:
    lw      s2, (BUTTERFLY_STACK_FRAME_SIZE - 16)(sp)
    lw      s1, (BUTTERFLY_STACK_FRAME_SIZE - 12)(sp)
    lw      s0, (BUTTERFLY_STACK_FRAME_SIZE - 8)(sp)
    lw      ra, (BUTTERFLY_STACK_FRAME_SIZE - 4)(sp)
    addi    sp, sp, BUTTERFLY_STACK_FRAME_SIZE
    ret

# ----------------------------------------------------
#  write_to_file: a0=filenamePtr, a1=bufferPtr, a2=lengthBytes
# ----------------------------------------------------
.globl write_to_file
write_to_file:
    addi    sp, sp, -16
    sw      ra,  0(sp)
    sw      s0,  4(sp)      # fd
    sw      s1,  8(sp)      # buffer
    sw      s2, 12(sp)      # length

    mv      s1, a1          # s1 = buffer
    mv      s2, a2          # s2 = length in bytes

    # open(filename, O_WRONLY|O_CREAT|O_TRUNC, 0666)
    # Using typical Linux/Spike syscall values, Venus might differ.
    # Assuming a linked 'open' function as per original structure.
    # Original flags: 0x601, mode: 0x1B6 (0666 octal)
    li      a1, 0x601       # Flags for open
    li      a2, 0x1B6       # Mode 0666
    call    open
    mv      s0, a0          # s0 = fd

    # write(fd, buffer, len)
    mv      a0, s0
    mv      a1, s1
    mv      a2, s2
    call    write

    # close(fd)
    mv      a0, s0
    call    close

    lw      s2, 12(sp)
    lw      s1,  8(sp)
    lw      s0,  4(sp)
    lw      ra,  0(sp)
    addi    sp, sp, 16
    ret

# ----------------------------------------------------
#  Data Section
# ----------------------------------------------------
    .section .data
    .align 2
input_data:
    .float 1.0, 0.0
    .float 2.0, 0.0
    .float 3.0, 0.0
    .float 4.0, 0.0
    .float 5.0, 0.0
    .float 6.0, 0.0
    .float 7.0, 0.0
    .float 8.0, 0.0
    .float 9.0, 0.0
    .float 10.0, 0.0
    .float 11.0, 0.0
    .float 12.0, 0.0
    .float 13.0, 0.0
    .float 14.0, 0.0
    .float 15.0, 0.0
    .float 16.0, 0.0

    .align 2
bitrev_output:
    .space 128

    .align 2
twiddles:
    .float 1.0,  0.0
    .float 0.9238795325112867, -0.3826834323650898
    .float 0.7071067811865475, -0.7071067811865475
    .float 0.3826834323650898, -0.9238795325112867
    .float 0.0, -1.0
    .float -0.3826834323650898, -0.9238795325112867
    .float -0.7071067811865475, -0.7071067811865475
    .float -0.9238795325112867, -0.3826834323650898
    .float -1.0, 0.0
    .float -0.9238795325112867, 0.3826834323650898
    .float -0.7071067811865475, 0.7071067811865475
    .float -0.3826834323650898, 0.9238795325112867
    .float 0.0, 1.0
    .float 0.3826834323650898, 0.9238795325112867
    .float 0.7071067811865475, 0.7071067811865475
    .float 0.9238795325112867, 0.3826834323650898

    .align 2
fft_output:
    .space 128

    .align 2
bitrev_filename:
    .string "bitrev_output.hex"
    .align 2
fft_filename:
    .string "fft_output.hex"

    .equ STDOUT_ADDR, 0xd0580000
