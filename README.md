# 🎙️ PCM-Based Voice Codec (CIE 337 Project)

This repository contains the implementation of a **PCM-based Voice Codec** developed as part of **CIE 337 - Digital Communications Project (Spring 2025)**.  
The project demonstrates the design and testing of a voice codec system for transmitting audio over data networks using **sampling, quantization, encoding, and decoding** techniques.

---

## 📌 Project Overview
The project implements a digital voice codec using **Pulse Code Modulation (PCM)**.  
It includes the following system blocks:

1. **Sampler**  
   - Samples input audio at user-specified frequencies.  

2. **Quantizer**  
   - Supports **Mid-rise** and **Mid-tread** uniform quantizers.  
   - Configurable number of levels `L` and peak amplitude `mp`.  
   - Outputs quantized signals, Mean Square Quantization Error (MSQE), and bitstreams.  

3. **Encoder**  
   - Converts bitstreams into signals using:  
     - **Unipolar NRZ signaling**  
     - **Polar NRZ signaling**  
   - Configurable pulse amplitude and bit duration.  

4. **Decoder**  
   - Reconstructs quantized samples from encoded bitstreams.  
   - Outputs decoded audio files for analysis.

---

## 🧪 Testing

### **Test 1: Codec Performance (No Noise)**
- Import and process a **20s+ CD-quality audio file**.  
- Parameters:  
  - Sampling frequencies: `fs = 40, 20, 5 kHz`  
  - Quantizer levels: `L = 4, 8, 64`  
  - Encoder: Both signaling techniques (show first 10 bits).  
- Output: Decoded audio saved as `.wav`.

### **Test 2: Codec Performance (With AWGN Noise)**
- Encoder amplitude: `A = 2V`  
- Parameters:  
  - `fs = 40 kHz`, `L = 64`  
- Channel noise: `AWGN ∼ N(0, N0)` with `N0 = 1, 4, 16`.  
- Includes a **regenerative function** before decoding.  
- Output: Noisy decoded audio saved as `.wav`.

---

## 📂 Repository Structure


```bash
├── Test I/ # Test 1 (no noise)
│ ├── Test1.m # MATLAB script for Test 1
│ ├── decoded_L4.wav # Output with quantizer levels L=4
│ ├── decoded_L8.wav # Output with quantizer levels L=8
│ ├── decoded_L64.wav # Output with quantizer levels L=64
│ └── ... # Other audio results
├── Test II/ # Test 2 (with AWGN noise)
│ ├── Test2.m # MATLAB script for Test 2
│ ├── decoded_N0_1.wav # Output with noise variance N0=1
│ ├── decoded_N0_4.wav # Output with noise variance N0=4
│ ├── decoded_N0_16.wav # Output with noise variance N0=16
│ └── ... # Other audio results
├── Input Song.wav # Original input audio file (20+ seconds, CD quality)
├── Report.pdf # Final project report (documentation + results)
├── README.md # Project description and usage guide

