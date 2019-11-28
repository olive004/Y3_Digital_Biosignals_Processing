% Fifth Tutorial
clear all; close all; clc;
 
load 'emg'; % Load the EMG signal
fs = 1600; % Sampling frequency
L = length(emg); % Duration of the signal in samples
time_ax=[0:1/fs:(L-1)/fs]; % Time axis of the signal in seconds

figure(7); plot(time_ax, emg)
title('EMG signal'); xlabel('time (s)'); ylabel('signal strength');
 
fc = 2; % Filter cut-off frequency in (Hz)
wc = 2*fc/fs*pi; % Normalized cut-off frequency
 
% Generate Sinc function with a given sampling rate
t = -floor(L):floor(L); % Form the time axis of the Sinc function (theoretically it can be infinetly long)
sinc_func=wc*sinc((2*fc/fs)*t);

% 1. Plots of the Sinc function used for designing the FIR filter and of its DFT [10%]
figure(1); plot(t, sinc_func);
title('Sinc function w/ cut-off freq 2Hz', 'fontsize', 18);
legend('sinc'); xlabel('n'); ylabel('AU');
 
f_snc=fftshift(fft(sinc_func)); % Find the DFT of the sinc function 
f_ax =(-pi+pi/length(sinc_func):2*pi/length(sinc_func):pi-pi/length(sinc_func))./pi; % Frequency axis for the DFT of sinc function

% 1. Plots of the Sinc function used for designing the FIR filter and of its DFT [10%]
figure(2); plot(f_ax,abs(f_snc));
title('Magnitude of DFT of the Sinc', 'fontsize', 18);
legend('DFT(sinc)'); xlabel('Frequency (rad)'); ylabel('AU');
 

% 2. Plot obtained by two filters
% 3. hanning / rectwin
% Create the FIR filter by truncating the very long sinc function and by using
% one of the following windows: hanning; rectwin;
for filter_lengths = [10, 100, 1000]

    FIR_duration= filter_lengths;
    truncation_section=floor(length(sinc_func)/2-FIR_duration/2):floor(length(sinc_func)/2+FIR_duration/2);
    fir_filt = sinc_func(truncation_section).*rectwin(length(truncation_section))';


    % Create moving avarege filter with the lenght of MA_coef_num
    MA_coef_num= filter_lengths;
    MA = ones(1,MA_coef_num)/MA_coef_num; % Impulse response of the moving average filter (see slide 30)

    figure(3);freqz(fir_filt);     
    figure(4);freqz(MA); 


    % Rectified EMG
    rect_EMG = abs(emg);

    % Filter the rectified EMG using FIR filter
    emg_FIR= conv(rect_EMG, fir_filt);
    emg_FIR = emg_FIR(length(fir_filt)/2:end-length(fir_filt)/2);   % Selecting relevant part of filtered
    emg_FIR_norm = emg_FIR./(sqrt(sum(abs(emg_FIR.^2)) / length(emg_FIR)));  % Normalizing signal using Root-mean-square
    
    figure(5); plot(time_ax, emg_FIR_norm, 'Linewidth', 2); hold on
    title(['Normalized, rectified EMG filtered w FIR'], 'fontsize', 18); xlabel('time (s)'); ylabel('signal strength (AU)');
    lgd = legend(num2str(10), num2str(100), num2str(1000)); title(lgd, 'Filter length')
    
    % Filter the rectified EMG using moving avarage filter
    emg_MA = conv(rect_EMG, MA);
    emg_MA = emg_MA(MA_coef_num/2:end-MA_coef_num/2);
    emg_MA_norm = emg_MA./(sqrt(sum(abs(emg_MA.^2)) / length(emg_MA)));  % Normalizing signal using RMS
    
    
    figure(6); plot(time_ax, emg_MA_norm, 'Linewidth', 2); hold on
    title(['Normalized, rectified EMG filtered w MA'], 'fontsize', 18); xlabel('time (s)'); ylabel('signal strength (AU)');
    lgd = legend(num2str(10), num2str(100), num2str(1000)); title(lgd, 'Filter length')
    
    
    if filter_lengths == 100
        % FIR filter 100 points, hanning vs rectwin
        figure(8); plot(time_ax, emg_FIR_norm);
        title('Norm. rect. filtered EMG, FIR rectangle window', 'fontsize', 18); xlabel('time (s)'); ylabel('signal strength (AU)');
        legend('EMG signal')
    end;

end;

