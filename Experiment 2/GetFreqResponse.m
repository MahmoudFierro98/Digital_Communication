function [X, f_axis] = GetFreqResponse(x,fs)
%
% Inputs:
%   x:      Signal in time domain
%   fs:     Sampling frequency
% Outputs:
%   X:      Magnitude of the signal in frequency domain
%   f_axis: Frequency axis
%
% This function generates the magnitude of the signal in the frequency
% domain after performing FFT on the signal in time domain.

N = length(x);
X = fftshift(fft(x));

f_axis = -fs/2:fs/N:fs/2-1/N;