function [fundFreq] = welchest(filename)
[sigsound, Fs] = audioread(filename);
[spectra,fb]= pwelch(sigsound,[],[],[],Fs);

decibel = pow2db(spectra);
 [maxval,idx] = max(decibel);
 fundFreq = fb(idx);
end

