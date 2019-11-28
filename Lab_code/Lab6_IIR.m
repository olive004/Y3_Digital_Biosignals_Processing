% Sixth tutorial.
close all; clear all; clc
 
load('EMG_6.mat'); % Load the EMG signal
L = length(EMG); % Duration of the signal in samples
 
Fs = 2500; % Sample frequency in Hz
t_ax = (0:L-1)/Fs; % Time axis of the signal in seconds
 
figure(1), plot(t_ax,EMG)
title('Intramuscular EMG signal')
xlabel('Time [s]')
ylabel('AU')
 
F_EMG = fftshift(fft(EMG)); % Find the DFT of the EMG 
f_ax = (-L/2:L/2-1)*Fs/L;
figure(2), plot(f_ax,abs(F_EMG));   % Plot of rectified EMG
xlabel('Frequency (Hz)'); ylabel('signal strenght (AU)');
title('Spectrum of the signal before the filtering','fontsize', 18);
legend('signal'); 
ylim([0 20*10^8]);
 
% Create band-stop filters to remove 60 Hz and harmonics
N = 3; % Filter order
band1 = [58 62];
[B1,A1] = butter(N,band1/(Fs/2),'stop'); % Generate filter coefficients
band2 = 2.*band1;
[B2,A2] = butter(N,band2/(Fs/2),'stop'); % Generate filter coefficients no. 2
band3 = 3.*band1;
[B3,A3] = butter(N,band3/(Fs/2),'stop'); % Generate filter coefficients no. 2
 
% Analyze the properties of the filter
H1 = tf(B1,A1,1/Fs,'variable','z^-1'); % Create transfer function object for filter 1
[z1,p1,k1] = tf2zp(B1,A1); % Calculate zeros and poles
figure(3), freqz(B1,A1);
title('Analysis in frequency domain of the filter that remove 60Hz-band frequency ')
figure(4), hold on, zplane(B1,A1);
title('z-plane to represent zeros and poles of the filter that remove 60Hz-band frequency')
legend('Zeros', 'Poles');

H2 = tf(B2,A2,1/Fs,'variable','z^-1'); % Create transfer function object
[z2,p2,k2] = tf2zp(B2,A2); % Calculate zeros and poles
figure(5), freqz(B2,A2);
title('Analysis in frequency domain of the filter that remove 120Hz-band frequency ')
figure(6), zplane(B2,A2);
title('z-plane to represent zeros and poles of the filter that remove 120Hz-band frequency')
legend('Zeros', 'Poles');



H3 = tf(B3,A3,1/Fs,'variable','z^-1'); % Create transfer function object
[z3,p3,k3] = tf2zp(B3,A3); % Calculate zeros and poles
figure(7), freqz(B3,A3);
title('Analysis in frequency domain of the filter that remove 180Hz-band frequency ')
figure(8), zplane(B3,A3);
title('z-plane to represent zeros and poles of the filter that remove 180Hz-band frequency')
legend('Zeros', 'Poles');




% Filter the signal
% 'Plot the magnitude of the DFT of the signal at each of 60Hz filter cascade'
figure(9), hold on, plot(t_ax,abs(F_EMG),'b'); %original signal
 %Filter 1
EMG_filt = filter(B1,A1,F_EMG);
figure(9), plot(t_ax,abs(EMG_filt),'r');
title('EMG filtered w/ 60Hz 1st harmonic','fontsize', 18);
legend('EMG', 'Filtered EMG');
xlabel('Frequency (Hz)'); ylabel('signal strenght (AU)');
% Filter 2
figure(10), hold on, plot(t_ax,abs(F_EMG),'b'); %original signal
EMG_filt = filter(B2,A2,EMG_filt);
figure(10), plot(t_ax,abs(EMG_filt),'r');
title('EMG filtered w/ 60Hz 1st and 2nd harmonic','fontsize', 18)
legend('EMG', 'Filtered EMG');
xlabel('Frequency (Hz)'); ylabel('signal strenght (AU)');
%Filter 3
figure(11), hold on, plot(t_ax,abs(F_EMG),'b'); %original signal
EMG_filt = filter(B3,A3,EMG_filt);
figure(11), plot(t_ax,abs(EMG_filt),'r');
title('EMG filtered w/ 60Hz 1st, 2nd, & 3rd harmonic','fontsize', 18)
legend('EMG', 'Filtered EMG');
xlabel('Frequency (Hz)'); ylabel('signal strenght (AU)');




