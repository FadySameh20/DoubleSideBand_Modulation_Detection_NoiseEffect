clear;
close all;

% Exp 1

[S, Fs] = audioread('eric.wav'); % read the attached audio file

% Fourier transform
L = length(S);
F = fftshift(fft(S));
f = Fs/2*linspace(-1,1,L);

% Plot
figure; plot(f, abs(F)/L); title('Original Signal Spectrum');

% Filter as 4kHz
W = 4000;
F(f>=W | f<=-W) = 0;
y =ifft(ifftshift(F));

% Plot after filter
L = length(y);
F = fftshift(fft(y));
f = Fs/2*linspace(-1,1,L);

figure; plot(f, abs(F)/L); title('Filtered Signal Spectrum');

%calculate time vector
tstart = 0;
tend = tstart + length(y) / Fs;
t1 = linspace(tstart, tend, length(y));
t1 = t1';
figure; plot(t1, y); title('Filtered Signal Time Domain');

% Hear the signal
sound(abs(y),Fs);


% Calculate constants
fm = W;
fc = 100000;
mu = 0.5;
Am = max(y);
Ac = Am/mu;

% resample at 5Fc
y = resample(y,5*fc,Fs);
Fs = 5*fc;

% DSBSC

% calculate time vector
tstart = 0;
tend = tstart + length(y) / Fs;
t1 = linspace(tstart, tend, length(y));
t1 = t1';

% DSBSC generation
carrier_signal = Ac .* cos(2*pi*fc*t1);
DSBSC = y.*carrier_signal; % y = modulating signal

% fourier transform
L = length(DSBSC);
F = fftshift(fft(DSBSC));
f = Fs/2*linspace(-1,1,L);

% Plot DSBSC frequency domain
figure; plot(f, abs(F) / L); title('DSBSC Frequency Domain');


% DSBTC
DSBTC=(1 + mu * y / Am) .* carrier_signal;

% fourier transform
L = length(DSBTC);
F = fftshift(fft(DSBTC));
f = Fs/2*linspace(-1,1,L);

% Plot DSBTC frequency  domain
figure; plot(f, abs(F) / L); title('DSBTC Frequency Domain');

% envelope detector DSBSC
envelopeDSBSC = abs(hilbert(DSBSC));

% Plot DSBSC envelope detector time domain
figure; plot(t1, DSBSC);
hold on;
plot(t1,-envelopeDSBSC,'r-',t1, envelopeDSBSC,'-r','Linewidth',1.5); % phase reversal occurs
hold off;
title('DSBSC time (blue) domain with envelope detector (red)');
ylim([-2 2])
xlim([2 2.5])

% ressample to sound the signal
envelopeDSBSC = resample(abs(envelopeDSBSC), Fs/5, Fs);
sound(abs(envelopeDSBSC), Fs/5) % not detected very well, distorted sound


% envelope detector DSBTC
envelopeDSBTC = abs(hilbert(DSBTC));

figure; plot(t1, DSBTC);
hold on;
plot(t1,-envelopeDSBTC,'r-',t1,envelopeDSBTC,'-r','Linewidth',1.5); % no phase reversal occurs
title('DSBTC time (blue) domain with envelope detector (red)');
hold off;
ylim([-5 5])
xlim([3 3.5])

% resample to sound the signal
envelopeDSBTC = resample(envelopeDSBTC, Fs/5, Fs);
sound(abs(envelopeDSBTC), Fs/5) % more accurately detected, less distortion


% Coherent Detection DSBSC

% 0 SNR
% generate signal+noise
noisy_DSBSC_0dB = awgn(DSBSC, 0);

% demodulate using coherent detector
demodulated = noisy_DSBSC_0dB.*cos(2*pi*fc*t1);

% fourier transform
demodulated_FFT = fftshift(fft(demodulated));

% LPF at Fm
demodulated_FFT(f>=W | f<=-W) = 0;

% inverse fourier transform to get demodulated signal in time domain
demodulated = ifft(ifftshift(demodulated_FFT));

% plot demodulated signal in time domain
figure; plot(t1, demodulated); title('0 SNR demodulated signal in time domain');

% fourier transform
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2*linspace(-1,1,L);

% plot demodulated signal in frequency domain
figure; plot(f, abs(F) / L); title('0 SNR demodulated signal in frequency domain');

% resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);



% 10 SNR
% generate signal+noise
noisy_DSBSC_10dB = awgn(DSBSC, 10);

% demodulate using coherent detector
demodulated = noisy_DSBSC_10dB.*cos(2*pi*fc*t1);

% fourier transform
demodulated_FFT = fftshift(fft(demodulated));

% LPF at Fm
demodulated_FFT(f>=W | f<=-W) = 0;

% inverse fourier transform to get demodulated signal in time domain
demodulated = ifft(ifftshift(demodulated_FFT));

% plot demodulated signal in time domain
figure; plot(t1, demodulated); title('10 SNR demodulated signal in time domain');

% fourier transform
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2*linspace(-1,1,L);

% plot demodulated signal in frequency domain
figure; plot(f, abs(F) / L); title('10 SNR demodulated signal in frequency domain');

% resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);


% 30 SNR
% generate signal+noise
noisy_DSBSC_30dB = awgn(DSBSC, 30);

% demodulate using coherent detector
demodulated = noisy_DSBSC_30dB.*cos(2*pi*fc*t1);

% fourier transform
demodulated_FFT = fftshift(fft(demodulated));

% LPF at Fm
demodulated_FFT(f>=W | f<=-W) = 0;

% inverse fourier transform to get demodulated signal in time domain
demodulated = ifft(ifftshift(demodulated_FFT));

% plot demodulated signal in time domain
figure; plot(t1, demodulated); title('30 SNR demodulated signal in time domain');

% fourier transform
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2*linspace(-1,1,L);

% plot demodulated signal in frequency domain
figure; plot(f, abs(F) / L); title('30 SNR demodulated signal in frequency domain');

% resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);


% Coherent Detection with frequency error
fc = 100100;

% 0 SNR
% generate signal+noise
noisy_DSBSC_0dB = awgn(DSBSC, 0);

% demodulate using coherent detector
demodulated = noisy_DSBSC_0dB.*cos(2*pi*fc*t1);

% fourier transform
demodulated_FFT = fftshift(fft(demodulated));

% LPF at Fm
demodulated_FFT(f>=W | f<=-W) = 0;

% inverse fourier transform to get demodulated signal in time domain
demodulated = ifft(ifftshift(demodulated_FFT));

% plot demodulated signal in time domain
figure; plot(t1, demodulated); title('0 SNR demodulated signal with frequency error in time domain');

% fourier transform
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2*linspace(-1,1,L);

% plot demodulated signal in frequency domain
figure; plot(f, abs(F) / L); title('0 SNR demodulated signal with frequency error in frequency domain');

% resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);

% 10 SNR
% generate signal+noise
noisy_DSBSC_10dB = awgn(DSBSC, 10);

% demodulate using coherent detector
demodulated = noisy_DSBSC_10dB.*cos(2*pi*fc*t1);

% fourier transform
demodulated_FFT = fftshift(fft(demodulated));

% LPF at Fm
demodulated_FFT(f>=W | f<=-W) = 0;

% inverse fourier transform to get demodulated signal in time domain
demodulated = ifft(ifftshift(demodulated_FFT));

% plot demodulated signal in time domain
figure; plot(t1, demodulated); title('10 SNR demodulated signal with frequency error in time domain');

% fourier transform
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2*linspace(-1,1,L);

% plot demodulated signal in frequency domain
figure; plot(f, abs(F) / L); title('10 SNR demodulated signal with frequency error in frequency domain');

% resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);


% 30 SNR
% generate signal+noise
noisy_DSBSC_30dB = awgn(DSBSC, 30);

% demodulate using coherent detector
demodulated = noisy_DSBSC_30dB.*cos(2*pi*fc*t1);

% fourier transform
demodulated_FFT = fftshift(fft(demodulated));

% LPF at Fm
demodulated_FFT(f>=W | f<=-W) = 0;

% inverse fourier transform to get demodulated signal in time domain
demodulated = ifft(ifftshift(demodulated_FFT));

% plot demodulated signal in time domain
figure; plot(t1, demodulated); title('30 SNR demodulated signal with frequency error in time domain');

% fourier transform
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2*linspace(-1,1,L);

% plot demodulated signal in frequency domain
figure; plot(f, abs(F) / L); title('30 SNR demodulated signal with frequency error in frequency domain');

% resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);


% Coherent Detection with phase error
fc = 100000;


% 0 SNR
% generate signal+noise
noisy_DSBSC_0dB = awgn(DSBSC, 0);

% demodulate using coherent detector
demodulated = noisy_DSBSC_0dB.*cos(2*pi*fc*t1 + pi/9);

% fourier transform
demodulated_FFT = fftshift(fft(demodulated));

% LPF at Fm
demodulated_FFT(f>=W | f<=-W) = 0;

% inverse fourier transform to get demodulated signal in time domain
demodulated = ifft(ifftshift(demodulated_FFT));

% plot demodulated signal in time domain
figure; plot(t1, demodulated); title('0 SNR demodulated signal with phase error in time domain');

% fourier transform
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2*linspace(-1,1,L);

% plot demodulated signal in frequency domain
figure; plot(f, abs(F) / L); title('0 SNR demodulated signal with phase error in frequency domain');

% resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);



% 10 SNR
% generate signal+noise
noisy_DSBSC_10dB = awgn(DSBSC, 10);

% demodulate using coherent detector
demodulated = noisy_DSBSC_10dB.*cos(2*pi*fc*t1 + pi/9);

% fourier transform
demodulated_FFT = fftshift(fft(demodulated));

% LPF at Fm
demodulated_FFT(f>=W | f<=-W) = 0;

% inverse fourier transform to get demodulated signal in time domain
demodulated = ifft(ifftshift(demodulated_FFT));

% plot demodulated signal in time domain
figure; plot(t1, demodulated); title('10 SNR demodulated signal with phase error in time domain');

% fourier transform
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2*linspace(-1,1,L);

% plot demodulated signal in frequency domain
figure; plot(f, abs(F) / L); title('10 SNR demodulated signal with phase error in frequency domain');

% resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);


% 30 SNR
% generate signal+noise
noisy_DSBSC_30dB = awgn(DSBSC, 30);

% demodulate using coherent detector
demodulated = noisy_DSBSC_30dB.*cos(2*pi*fc*t1 + pi/9);

% fourier transform
demodulated_FFT = fftshift(fft(demodulated));

% LPF at Fm
demodulated_FFT(f>=W | f<=-W) = 0;

% inverse fourier transform to get demodulated signal in time domain
demodulated = ifft(ifftshift(demodulated_FFT));

% plot demodulated signal in time domain
figure; plot(t1, demodulated); title('30 SNR demodulated signal with phase error in time domain');

% fourier transform
L = length(demodulated);
F = fftshift(fft(demodulated));
f = Fs/2*linspace(-1,1,L);

% plot demodulated signal in frequency domain
figure; plot(f, abs(F) / L); title('30 SNR demodulated signal with phase error in frequency domain');

% resample to sound the demodulated signal
demodulated = resample(abs(demodulated), Fs/5, Fs);
sound(abs(demodulated), Fs/5);
