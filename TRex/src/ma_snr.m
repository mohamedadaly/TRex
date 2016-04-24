function [snr] = ma_snr(gt, noise)
%MA_SNR Computes SNR between input gt and approximation.

% make sure they have the same size
assert(all(size(gt) == size(noise)));

% sum square of gt
num = sum(gt(:).^2);
% sum square of noise 
den = sum(noise(:).^2);

snr = 10 * log10(num / den);

end

