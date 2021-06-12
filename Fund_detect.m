function [FundFreq] = Fund_detect(data,SampleFreq)
[spectra,Freq]= pwelch(data,[],[],[],SampleFreq);
decibel = pow2db(spectra);
 [maxval,idx] = max(decibel);
Fundfreq=Freq(idx)
end

