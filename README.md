# DoubleSideBand_Modulation_Detection_NoiseEffect
Double Side Band (DSB) modulation for both suppressed carrier and transmitted carrier, detection using envelope detector and synchronous detector in addition to showing the effect of noise on the signal with different SNR (signal to noise ratio) values.

* Procedures:
1. Use Matlab to read the attached audio file, which has a sampling frequency Fs= 48 KHz. Find the spectrum of this signal (the signal in frequency domain). [audioread, fft , fftshift , plot]
2. Using an ideal Filter, remove all frequencies greater than 4 KHz.
3. Obtain the filtered signal in time domain and frequency domain, this is a band limited signal of BW=4KHz. [ifftshift ,ifft]
4. sound the filtered audio signal (make sure that there is only a small error in the filtered signal) [sound] 
5. Modulate the carrier with the filtered signal you obtained, you are required to generate both types of modulation (DSB-TC and DSB-SC). Choose a carrier frequency of 100 KHz. For the DSB-TC take the DC bias added to message before modulation to be twice the maximum of the message (modulation index =0.5 in this case). Note: You will also need to increase the sampling frequency of the filtered audio signal, the sampling frequency must be at least 2 times the carrier frequency, In this simulation use Fs = 5 𝐹𝑐 [resample] You have to sketch the modulated signal of both DSB-TC & DSB-SC in frequency domain.
6. For both types of modulations (DSB-SC & DSB-TC), use envelop detector to receive the message (assume no noise). Note: to obtain the envelope you can use the following matlabcommand.𝑒𝑛𝑣𝑒𝑙𝑜𝑝𝑒 = 𝑎𝑏𝑠(ℎ𝑖𝑙𝑏𝑒𝑟𝑡(𝑚𝑜𝑑𝑢𝑙𝑎𝑡𝑒𝑑 𝑠𝑖𝑔𝑛𝑎𝑙))
7. After the reception of both modulation types using envelope detector, sketch the received signal in time domain, and Play the received signal back (Note: to sound signal after demodulation process you have to decrease the sampling frequency again). What observation can you make of this or which type of modulation the envelope detector can be used with? For DSB-SC, perform steps 8-10.
8. Use coherent detection to receive the modulated signal with SNR=0, 10, 30dB then sound the received signals and plot them in both time and frequency domain.
9. Repeat the coherent detection with frequency error, F=100.1 KHz instead of 100 KHz and find the error. Do you have a name for this phenomenon ? 
10. Repeat the coherent detection with phase error = 200.
