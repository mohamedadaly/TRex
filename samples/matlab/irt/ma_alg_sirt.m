function [rec, times, snrs, iters] = ma_alg_sirt(in_params, alg_params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[times, snrs, iters] = deal(zeros(alg_params.iter, 1));

% init in the transpose to make row-major
rec = zeros(size(in_params.gt_vol))';

% take transpose to make row major because the projection matrix assumes
% row-major sinogram, and gives row-major volume
sino = in_params.sino';
[ndet, nviews] = size(sino);

% sum of rows
sum_rows = sum(in_params.A, 2);
sum_rows(abs(sum_rows)<1e-16) = 1;
% sum over pixels 
sum_cols = sum(in_params.A, 1)';
sum_cols(abs(sum_cols) < 1e-16) = 1;

% precondition
if alg_params.precon
  [ny, nx] = size(rec);
  prec = zeros(size(rec));
  prec(floor(ny/2)+1, floor(nx/2)+1) = 1;
  prec = reshape(in_params.A' * (in_params.A * prec(:)), ny, nx);
  prec = abs(fft2(prec));
  % add constant
  prec = prec + alg_params.gamma * max(prec(:));
  
  % load
  xx = load('pre.mat');
  prec = (xx.p1 + 5e3 * xx.p2);
  prec = prec / max(prec(:));
end

% iterations
for it = 1:alg_params.iter
  % loop on views
  tic;
    
  % forward project current view
  Ax = in_params.A * rec(:);

  % difference
  diff = (sino(:) - Ax) ./ sum_rows;

  % backproject
  bp = (in_params.A' * diff) ./ sum_cols;
  alg_params.gamma = max(bp(:));

  % precondition
  if alg_params.precon
    bp = real(ifft2(fft2(reshape(bp,ny,nx)) ./ (prec')));
    bp = bp(:) / max(bp(:)) * alg_params.gamma; % / max(bp(:)) * .2; % * 1000;
  end

  % update
%   figure(1), imshow(rec',[]); colorbar;
  rec(:) = rec(:) + alg_params.alpha * bp; 
%   figure(2), imshow(rec', []); colorbar;
%   pause;

  if ~isempty(alg_params.min_val)
    rec = max(alg_params.min_val, rec);
  end
  
  times(it) = toc;
  snrs(it) = ma_snr(in_params.gt_vol, in_params.gt_vol - rec');
  iters(it) = it;
  
end

% transpose
rec = rec';

times = cumsum(times);
end

