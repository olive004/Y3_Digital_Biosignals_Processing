% Third tutorial
clear all
close all
clc
 
% Signal loading
load('EEG.mat');
 
% Sampling frequency
fsamp = 512;
 
% Select duration to analyze
% Change Duration to change signal interval to analyze
Duration_s = 1; % Duration in seconds (max 15 seconds)
Duration = round(Duration_s*fsamp);
EEG = EEG(1:Duration); 
 
% Signal duration in samples
L = length(EEG);
 
% Plot the EEG signal
time_ax = [0:L-1]./fsamp;
figure(1)
plot(time_ax, EEG);
xlabel('Time (s)')
title(['EEG signal over ',num2str(Duration_s),'s'])
ylabel('Amplitude EEG (Arbitrary Units)')
 
% Compute DFT of the two channels
mean_EEG = mean(EEG);
X1 = fft( EEG - mean_EEG );
 
% Compute PSD (power spectral density) of the two channels and adjust
% frequency axis according to Matlab notation
PSD1 = fftshift(abs(X1).^2)/L;

% Build the frequency axis in radiants
freq_a_rad = [-pi+pi/L:2*pi/L:pi-pi/L];
% Convert the frequency axis in Hz
freq_a_Hz = freq_a_rad./(2*pi).*fsamp;
 
% Plot the PSD of the two channels with frequencies in radiants and Hz
figure(2), subplot(2,1,1), plot(freq_a_rad,PSD1);
xlabel('Frequency (radians)')
title('Power Spectral Density of EEG')
ylabel('PSD (Arbitrary Units)')

figure(2), subplot(2,1,2), plot(freq_a_Hz,PSD1);
xlabel('Frequency (Hz)')
title('Power Spectral Density of EEG')
ylabel('PSD (Arbitrary Units)')
 
% Compute the percentage of power in different subbands
% [Here please complete with instructions for computing the relative power in the frequency bands.]

% Percentage of power in freq bands

subbands = [0.5,4,8,13,30,42];  % alpha to theta subbands (Hz) 
% Get indices of frequencies correspoding to the PSD values
for i=1:length(subbands)
    subbands_idx(i) = find(abs(freq_a_Hz - subbands(i)) < 0.5);  
end

% Calculating percent power
sum_PSD1 = sum(PSD1);    % Speed up calculation
for k = 1:(length(subbands_idx)-1)
    perc_power(k) = sum(PSD1(subbands_idx(k):subbands_idx(k+1)))./sum_PSD1.*100;
end


fprintf('Results: Percent power for %2d s signal interval is \n %8.3f for delta band, \n %8.3f for theta band, \n %8.3f for alpha band,\n %8.3f for beta band,\n and %8.3f for gamma band\n',Duration_s,  perc_power)






