.section .text
    .global _start
    .option norvc          # no compressed instructions

# ------------------------------------------------------------
#  RISC‑V RV32 1‑D FFT (N=32 complex samples)
# ------------------------------------------------------------
_start:
    # 1) Bit-reverse input_data → bitrev_output
    la      a0, input_data
    la      a1, bitrev_output
    li      a2, 32
    call    bit_reverse

    # --- Save bit-reversed output to file ---
    la      a0, bitrev_filename    # filename string
    la      a1, bitrev_output      # buffer address
    li      a2, 256                # bytes (32 complex × 8 bytes)
    call    write_to_file

    # 2) FFT butterfly using twiddles
    la      a0, bitrev_output
    la      a1, twiddles
    la      a2, fft_output
    li      a3, 32
    call    butterfly_stages

    # --- Save FFT output to file ---
    la      a0, fft_filename        # filename string
    la      a1, fft_output          # buffer address
    li      a2, 256                 # bytes (32 complex × 8 bytes)
    call    write_to_file

    # Exit program properly
    li      a7, 93          # sys_exit
    li      a0, 0           # exit status
    ecall

# ----------------------------------------------------
#  System call implementations
# ----------------------------------------------------
.globl open
open:
    li      a7, 56          # sys_openat
    mv      a3, a2          # mode
    mv      a2, a1          # flags  
    mv      a1, a0          # filename
    li      a0, -100        # AT_FDCWD
    ecall
    ret

.globl write
write:
    li      a7, 64          # sys_write
    ecall
    ret

.globl close
close:
    li      a7, 57          # sys_close
    ecall
    ret

# ----------------------------------------------------
#  bit_reverse: a0=inPtr, a1=outPtr, a2=N
# ----------------------------------------------------
.globl bit_reverse
bit_reverse:
    li      t0, 0
    mv      t1, a2 
bitrev_loop:
    beq     t0, t1, bitrev_done

    # For N=32, need a 5-bit reversal
    # Extract individual bits
    andi    t2, t0, 1       # bit 0
    slli    t2, t2, 4       # shift to position 4
    
    andi    t3, t0, 2       # bit 1
    slli    t3, t3, 2       # shift to position 3
    
    andi    t4, t0, 4       # bit 2
    # bit 2 stays in position 2
    
    andi    t5, t0, 8       # bit 3
    srli    t5, t5, 2       # shift to position 1
    
    andi    t6, t0, 16      # bit 4
    srli    t6, t6, 4       # shift to position 0
    
    # Combine the reversed bits
    or      s7, t2, t3
    or      s7, s7, t4
    or      s7, s7, t5
    or      s7, s7, t6

    # load input
    slli    s0, t0, 3
    add     s0, a0, s0
    flw     f0, 0(s0)
    flw     f1, 4(s0)

    # store to bitrev_output[r]
    slli    s0, s7, 3
    add     s0, a1, s0
    fsw     f0, 0(s0)
    fsw     f1, 4(s0)

    addi    t0, t0, 1
    j       bitrev_loop
bitrev_done:
    ret

# ----------------------------------------------------
#  butterfly_stages: a0=inPtr, a1=twiddlePtr, a2=outPtr, a3=N
# ------------------------------------------------------------
.globl butterfly_stages
butterfly_stages:
    addi    sp, sp, -32
    sw      ra, 0(sp)
    sw      s0, 4(sp)
    sw      s1, 8(sp)
    sw      s2, 12(sp)
    sw      s3, 16(sp)
    sw      s4, 20(sp)
    sw      s5, 24(sp)
    sw      s6, 28(sp)
    
    # Copy input to output for first stage preparation
    li      t0, 0
init_copy:
    beq     t0, a3, init_copy_done
    slli    t1, t0, 3         # t1 = i * 8 (byte offset)
    add     t2, a0, t1        # t2 = input + offset
    add     t3, a2, t1        # t3 = output + offset
    
    flw     f0, 0(t2)         # load real part
    flw     f1, 4(t2)         # load imag part
    fsw     f0, 0(t3)         # store real part to output
    fsw     f1, 4(t3)         # store imag part to output
    
    addi    t0, t0, 1
    j       init_copy
init_copy_done:

    # Process FFT stages (len = 2, 4, 8, 16, 32)
    li      t0, 2             # len = 2 (first stage)
    
stage_loop:
    bgt     t0, a3, stage_done
    srli    t1, t0, 1         # half = len/2
    
    # Process blocks within this stage
    li      t2, 0             # block_start = 0
block_loop:
    bge     t2, a3, block_done
    
    # Process butterflies within this block
    li      t3, 0             # j = 0
    
butterfly_loop:
    bge     t3, t1, butterfly_done_block
    
    # Calculate indices
    add     s4, t2, t3        # s4 = block_start + j = a_idx
    add     s5, s4, t1        # s5 = a_idx + half = b_idx
    
    # Calculate twiddle factor index
    div     s6, a3, t0        # s6 = N/len
    mul     s6, t3, s6        # s6 = j * (N/len) = twiddle_idx
    
    # Load a (real and imaginary)
    slli    s0, s4, 3         # s0 = a_idx * 8 (byte offset)
    add     s0, a2, s0        # s0 = output + a_offset
    flw     f0, 0(s0)         # f0 = a_real
    flw     f1, 4(s0)         # f1 = a_imag
    
    # Load b (real and imaginary)
    slli    s1, s5, 3         # s1 = b_idx * 8 (byte offset)
    add     s1, a2, s1        # s1 = output + b_offset
    flw     f2, 0(s1)         # f2 = b_real
    flw     f3, 4(s1)         # f3 = b_imag
    
    # Load twiddle factor (real and imaginary)
    slli    s2, s6, 3         # s2 = twiddle_idx * 8 (byte offset)
    add     s2, a1, s2        # s2 = twiddles + twiddle_offset
    flw     f4, 0(s2)         # f4 = twiddle_real
    flw     f5, 4(s2)         # f5 = twiddle_imag
    
    # Complex multiplication: b * twiddle
    # Real part = b_real*tw_real - b_imag*tw_imag
    fmul.s  f6, f2, f4        # f6 = b_real * tw_real
    fmul.s  f7, f3, f5        # f7 = b_imag * tw_imag
    fsub.s  f8, f6, f7        # f8 = b_real*tw_real - b_imag*tw_imag
    
    # Imaginary part = b_real*tw_imag + b_imag*tw_real
    fmul.s  f9, f2, f5        # f9 = b_real * tw_imag
    fmul.s  f10, f3, f4       # f10 = b_imag * tw_real
    fadd.s  f11, f9, f10      # f11 = b_real*tw_imag + b_imag*tw_real
    
    # Butterfly operation
    # New a = a + b*twiddle
    fadd.s  f12, f0, f8       # f12 = a_real + (b*tw)_real = new_a_real
    fadd.s  f13, f1, f11      # f13 = a_imag + (b*tw)_imag = new_a_imag
    
    # New b = a - b*twiddle
    fsub.s  f14, f0, f8       # f14 = a_real - (b*tw)_real = new_b_real
    fsub.s  f15, f1, f11      # f15 = a_imag - (b*tw)_imag = new_b_imag
    
    # Store new a
    fsw     f12, 0(s0)        # Store new_a_real
    fsw     f13, 4(s0)        # Store new_a_imag
    
    # Store new b
    fsw     f14, 0(s1)        # Store new_b_real
    fsw     f15, 4(s1)        # Store new_b_imag
    
    addi    t3, t3, 1         # j++
    j       butterfly_loop
    
butterfly_done_block:
    add     t2, t2, t0        # block_start += len
    j       block_loop
    
block_done:
    slli    t0, t0, 1         # len *= 2
    j       stage_loop
    
stage_done:
    lw      ra, 0(sp)
    lw      s0, 4(sp)
    lw      s1, 8(sp)
    lw      s2, 12(sp)
    lw      s3, 16(sp)
    lw      s4, 20(sp)
    lw      s5, 24(sp)
    lw      s6, 28(sp)
    addi    sp, sp, 32
    ret

# ----------------------------------------------------
#  write_to_file: a0=filenamePtr, a1=bufferPtr, a2=lengthBytes
# ----------------------------------------------------
.globl write_to_file
write_to_file:
    addi    sp, sp, -16
    sw      ra,  0(sp)
    sw      s0,  4(sp)
    sw      s1,  8(sp)
    sw      s2, 12(sp)

    mv      s1, a1          # buffer
    mv      s2, a2          # length in bytes

    # open(filename, O_WRONLY|O_CREAT|O_TRUNC, 0666)
    li      a1, 0x601       # O_WRONLY|O_CREAT|O_TRUNC
    li      a2, 0x1B6       # 0666 permissions
    call    open
    mv      s0, a0          # fd

    # Check if file opened successfully
    bltz    s0, write_error

    # write(fd, buffer, len)
    mv      a0, s0
    mv      a1, s1
    mv      a2, s2
    call    write

    # close(fd)
    mv      a0, s0
    call    close

write_error:
    # restore
    lw      ra,  0(sp)
    lw      s0,  4(sp)
    lw      s1,  8(sp)
    lw      s2, 12(sp)
    addi    sp, sp, 16
    ret

# ----------------------------------------------------
#  Data Section
# ----------------------------------------------------
    .section .data
    .align 2
input_data:
    .float 1.0, 0.0    # Sample 0
    .float 2.0, 0.0    # Sample 1
    .float 3.0, 0.0    # Sample 2
    .float 4.0, 0.0    # Sample 3
    .float 5.0, 0.0    # Sample 4
    .float 6.0, 0.0    # Sample 5
    .float 7.0, 0.0    # Sample 6
    .float 8.0, 0.0    # Sample 7
    .float 9.0, 0.0    # Sample 8
    .float 10.0, 0.0   # Sample 9
    .float 11.0, 0.0   # Sample 10
    .float 12.0, 0.0   # Sample 11
    .float 13.0, 0.0   # Sample 12
    .float 14.0, 0.0   # Sample 13
    .float 15.0, 0.0   # Sample 14
    .float 16.0, 0.0   # Sample 15
    .float 17.0, 0.0   # Sample 16
    .float 18.0, 0.0   # Sample 17
    .float 19.0, 0.0   # Sample 18
    .float 20.0, 0.0   # Sample 19
    .float 21.0, 0.0   # Sample 20
    .float 22.0, 0.0   # Sample 21
    .float 23.0, 0.0   # Sample 22
    .float 24.0, 0.0   # Sample 23
    .float 25.0, 0.0   # Sample 24
    .float 26.0, 0.0   # Sample 25
    .float 27.0, 0.0   # Sample 26
    .float 28.0, 0.0   # Sample 27
    .float 29.0, 0.0   # Sample 28
    .float 30.0, 0.0   # Sample 29
    .float 31.0, 0.0   # Sample 30
    .float 32.0, 0.0   # Sample 31

    .align 2
bitrev_output:
    .space 256  # 32 complex numbers * 8 bytes

    .align 2
twiddles:
    # Twiddles for N=32: W_k = exp(-2πi * k / 32)
    # W_0 = exp(-2πi * 0/32) = 1.0 + 0.0i
    .float 1.0, 0.0
    # W_1 = exp(-2πi * 1/32) = 0.9808 - 0.1951i
    .float 0.9807852804032304, -0.19509032201612825
    # W_2 = exp(-2πi * 2/32) = 0.9239 - 0.3827i
    .float 0.9238795325112867, -0.3826834323650898
    # W_3 = exp(-2πi * 3/32) = 0.8315 - 0.5556i
    .float 0.8314696123025452, -0.5555702330196022
    # W_4 = exp(-2πi * 4/32) = 0.7071 - 0.7071i
    .float 0.7071067811865476, -0.7071067811865475
    # W_5 = exp(-2πi * 5/32) = 0.5556 - 0.8315i
    .float 0.5555702330196023, -0.8314696123025452
    # W_6 = exp(-2πi * 6/32) = 0.3827 - 0.9239i
    .float 0.3826834323650898, -0.9238795325112867
    # W_7 = exp(-2πi * 7/32) = 0.1951 - 0.9808i
    .float 0.19509032201612833, -0.9807852804032304
    # W_8 = exp(-2πi * 8/32) = 0.0 - 1.0i
    .float 0.0, -1.0
    # W_9 = exp(-2πi * 9/32) = -0.1951 - 0.9808i
    .float -0.19509032201612833, -0.9807852804032304
    # W_10 = exp(-2πi * 10/32) = -0.3827 - 0.9239i
    .float -0.3826834323650897, -0.9238795325112867
    # W_11 = exp(-2πi * 11/32) = -0.5556 - 0.8315i
    .float -0.5555702330196022, -0.8314696123025453
    # W_12 = exp(-2πi * 12/32) = -0.7071 - 0.7071i
    .float -0.7071067811865475, -0.7071067811865476
    # W_13 = exp(-2πi * 13/32) = -0.8315 - 0.5556i
    .float -0.8314696123025453, -0.5555702330196022
    # W_14 = exp(-2πi * 14/32) = -0.9239 - 0.3827i
    .float -0.9238795325112867, -0.3826834323650899
    # W_15 = exp(-2πi * 15/32) = -0.9808 - 0.1951i
    .float -0.9807852804032304, -0.19509032201612836
    # W_16 = exp(-2πi * 16/32) = -1.0 + 0.0i
    .float -1.0, 0.0
    # W_17 = exp(-2πi * 17/32) = -0.9808 + 0.1951i
    .float -0.9807852804032304, 0.1950903220161283
    # W_18 = exp(-2πi * 18/32) = -0.9239 + 0.3827i
    .float -0.9238795325112867, 0.3826834323650897
    # W_19 = exp(-2πi * 19/32) = -0.8315 + 0.5556i
    .float -0.8314696123025453, 0.5555702330196022
    # W_20 = exp(-2πi * 20/32) = -0.7071 + 0.7071i
    .float -0.7071067811865477, 0.7071067811865475
    # W_21 = exp(-2πi * 21/32) = -0.5556 + 0.8315i
    .float -0.5555702330196023, 0.8314696123025452
    # W_22 = exp(-2πi * 22/32) = -0.3827 + 0.9239i
    .float -0.3826834323650899, 0.9238795325112867
    # W_23 = exp(-2πi * 23/32) = -0.1951 + 0.9808i
    .float -0.19509032201612836, 0.9807852804032304
    # W_24 = exp(-2πi * 24/32) = 0.0 + 1.0i
    .float 0.0, 1.0
    # W_25 = exp(-2πi * 25/32) = 0.1951 + 0.9808i
    .float 0.1950903220161283, 0.9807852804032304
    # W_26 = exp(-2πi * 26/32) = 0.3827 + 0.9239i
    .float 0.3826834323650897, 0.9238795325112868
    # W_27 = exp(-2πi * 27/32) = 0.5556 + 0.8315i
    .float 0.5555702330196022, 0.8314696123025453
    # W_28 = exp(-2πi * 28/32) = 0.7071 + 0.7071i
    .float 0.7071067811865475, 0.7071067811865477
    # W_29 = exp(-2πi * 29/32) = 0.8315 + 0.5556i
    .float 0.8314696123025452, 0.5555702330196023
    # W_30 = exp(-2πi * 30/32) = 0.9239 + 0.3827i
    .float 0.9238795325112867, 0.3826834323650899
    # W_31 = exp(-2πi * 31/32) = 0.9808 + 0.1951i
    .float 0.9807852804032304, 0.19509032201612833

    .align 2
fft_output:
    .space 256  # 32 complex numbers * 8 bytes

    .align 2
bitrev_filename:
    .string "bitrev_output.hex"
    .align 2
fft_filename:
    .string "fft_output.hex"