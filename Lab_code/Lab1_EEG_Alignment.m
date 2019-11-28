% Lab 1: EEG signal alignment

close all
clear all
fclose('all')

for i_M = 1
M_idx = [1,2,4,8];    
load('Signals_EMG.mat'); % Loading the recorded EMGs (two channels)
n_step = 100; % Number of steps in the loop for optimal alignment
stepCount = 0.2; % Size of the step of non-integer delay applied to one of the signals for alignment
ied=24; % Interelectrode distance of the recordings in mm
Fs = 2048; % Recording sampling frequency
Ts = 1/Fs; % Sampling interval
 
% Downsampling by an integer
M = M_idx(i_M); % Downsampling factor
channel1 = channel1(1:M:end); % Downsampling first channel
channel2 = channel2(1:M:end); % Downsampling second channel
Fs=Fs/M;
Ts=Ts*M;
% End downsampling
 
timeAxis=[1:length(channel1)].*Ts.*1000; % Definition of time axis in ms
freqAxis=fftshift([-0.5:1/(length(channel1)):0.5-1/(length(channel1))]); % Definition of discrete frequency axis
 
figure(1);plot(timeAxis,channel1,'k');hold on; plot(timeAxis,channel2-1000,'k'); % Plot of the two recordings after down-sampling
xlabel('Time (ms)');
ylabel('Signals (AU)');
 
channel1_ft = fft(channel1); % Fourier transform of the first channel
 
figure(2)
for uu = 1 : n_step
    channel1_dt = (channel1_ft).*exp(-i*2*pi*stepCount*uu*freqAxis); % complex exponential multiplication (delay in frequency domain)
    channel1_dt = real(ifft((channel1_dt))); % inverse transform to go back to the time domain 
    plot(timeAxis,channel1_dt,'r'); % Plot of the time-shifted signal
    hold on
    plot(timeAxis,channel2,'k');
    uu;
    MSE_vect(uu)= sum((channel1_dt - channel2).^ 2)./sum(channel2.^ 2).*100; % normalized mean square error between aligned signals
    delay(uu) = stepCount*uu; % Imposed delay in samples
end;
xlabel('Time (ms)')
ylabel('Signal amplitude (AU)')
 
% Identification of the optimal delay (minimum mean sqaure error)
[MSEopt, optDelay] = min(MSE_vect);

% Plot of optimal alignment
est_cond_vel = ied/(delay(optDelay)*Ts*1000);

channel1_opt_delay = (channel1_ft).*exp(-i*2*pi*stepCount*optDelay*freqAxis); % complex exponential multiplication (delay in frequency domain)
channel1_opt_delay = real(ifft((channel1_opt_delay))); % inverse transform to go back to the time domain 
figure(3)
plot(timeAxis, channel1_opt_delay, 'r')
hold on
plot(timeAxis, channel2, 'k')
title(['Optimal Alignment M=', num2str(M)])
xlabel('Time (ms)')
ylabel('Signal amplitude (AU)')
textbox = sprintf('Step size: %3.2f\nNum steps: %3.2f\nDownsampling factor:%3.2f\n\nOptimal Delay: %3.2fms\nEst. conduction velocity:%3.2fm/s', stepCount, n_step, M,delay(optDelay)*Ts*1000,est_cond_vel);
text(60, 400, textbox);
hold off


% Delay and conduction velocity estimate
fprintf('The optimal delay is %2.2f ms \n',delay(optDelay)*Ts*1000);
fprintf('The estimated conduction velocity is %2.2f m/s \n',est_cond_vel);
fprintf('Optimal MSE between signals: %2.2f %%\n',MSEopt);







% Plot of estimation error as a function of delay
figure(4)
color = ['r','b','g','k'];
plot(delay, MSE_vect, color(i_M))
hold on
title(['Estimation Error vs. Delay']);
xlabel('Delay in s')
legend(num2str(M_idx(1)),num2str(M_idx(2)),num2str(M_idx(3)),num2str(M_idx(4)))

end


