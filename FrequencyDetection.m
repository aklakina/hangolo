function [fundFreq] = FrequencyDetection(inputData, SampleFrequency)

n = 2^nextpow2(length(inputData)); %point of FFT

Freq = linspace(0, 1, fix(n/2))*SampleFrequency/2;
Iv = 1:length(Freq);
y0 = abs(fft(inputData,n))/n; % calculate magnitude spectrum of Signal
y0 = y0(Iv); % One-sided spectrum
%Harmonic Product Spectrum

hps1 = downsample(y0,1);
hps2 = downsample(y0,2);
hps3 = downsample(y0,3);
hps4 = downsample(y0,4);
hps5 = downsample(y0,5);
HPS = zeros(1,length(hps5));
for i = 1:length(hps5)
prod = hps1(i)*hps2(i)*hps3(i)*hps4(i)*hps5(i);
HPS(i) =[prod];
end

%Finding the fundamental frequency

[maxvalHPS,idx] = max(HPS);
fundFreq=Freq(idx);
end

