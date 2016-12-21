% IMPORTANT NOTE:
%
% This function was modified from the following source:

%%%%%%%%%%%%%%%%%%%%%%

%FileName - genSignalForSNR.m
 
%Author: Mathuranathan
 
%http://www.gaussianWaves.com
 
%Creative Commons - cc by-nc-sa

%%%%%%%%%%%%%%%%%%%%%
 
 
function [signalWithNoise,noise,measuredSNR]= genSignalForSNR(signal,SNR)
 
% Generate Poisson Noise independently for each pixel such that the mean
% level of distribution (lambda) is equal to that particular pixel value

[m,n]=size(signal);

noise = zeros(size(signal));

for i= 1:m
    for j= 1:n
        noise(i,j)=poissrnd(signal(i,j), 1);
    end
end

 
%Scale the input signal accordingly for the given SNR.
 
scaledSignal = std(noise)/std(signal)*(sqrt(10^(SNR/10)))*signal;
 
%calculate Signal power and noise power
 
%signalPower = (norm(scaledSignal)^2)/length(scaledSignal);
%noisePower = (norm(noise)^2)/length(noise);
 
%Alternative way of calculating Signal and noise power from their variance
 
signalPower = var(scaledSignal); 
noisePower = var(noise);
 
 
%Calculate Signal to noise ratio for the scaledSignal and generated Noise
 
SNRratio = signalPower/noisePower;
measuredSNR=10*log10(SNRratio);
 
%Add the scaled signal with the generated noise
 
signalWithNoise=scaledSignal+noise;
 
end