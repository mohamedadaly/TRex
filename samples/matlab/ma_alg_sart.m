function [rec, times, snrs, iters] = ma_alg_sart(in_params, alg_params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[times, snrs, iters] = deal(zeros(alg_params.iter, 1));

% init in the transpose to make row-major
rec = zeros(size(in_params.gt_vol))';
[ny, nx] = size(rec);

% take transpose to make row major because the projection matrix assumes
% row-major sinogram, and gives row-major volume
sino = in_params.sino';
[ndet, nviews] = size(sino);

% precondition
if alg_params.precon
  prec = zeros(size(rec));
  prec(floor(ny/2)+1, floor(nx/2)+1) = 1;
  prec = reshape(in_params.A' * (in_params.A * prec(:)), ny, nx);
  prec = abs(fft2(prec));
  
  xx = load('pre.mat');
  prec = (xx.p1 + 5e3 * xx.p2);
  
end

% scale matrix and sinogram by sum_rows to solve a least-square instead of
% a weighted least-square
if alg_params.wls
  % sum of rows
  sum_rows = sum(in_params.A, 2);
  sum_rows(abs(sum_rows)<1e-16) = 1;
  % scale
  in_params.A = bsxfun(@(a,b)(a.*b), in_params.A, sum_rows);
  sino = reshape(sum_rows .* sino(:), ndet, nviews);
end

% initialize for proxisart
if alg_params.prox
  rec = in_params.prox_in;
  y = zeros(size(sino));
  in_params.A = sqrt(2 * alg_params.lambda) * in_params.A;
  sino = sqrt(2 * alg_params.lambda) * sino;
end;

% compute sum_cols
sum_cols = sum(in_params.A, 1)';
sum_cols(abs(sum_cols) < 1e-16) = 1;

% sum of rows
sum_rows = sum(in_params.A, 2);
sum_rows(abs(sum_rows)<1e-16) = 1;

% If BSSART, then initialize rec with one of the rows of A
if alg_params.BSSART
  rec = reshape(in_params.A(10,:), ny, nx);
end

% iterations
for it = 1:alg_params.iter
  % loop on views
  tic;
  pstart = 1;
  for v = 1:nviews
    % rows
    rows = pstart : pstart + ndet - 1;
    
    if ~alg_params.BSSART
      % sum over pixels for this view
      sum_cols = sum(in_params.A(rows, :), 1)';
      sum_cols(abs(sum_cols) < 1e-16) = 1;
    end
    
    % forward project current view
    Ax = in_params.A(rows, :) * rec(:);
    
    % difference
    if alg_params.prox
      diff = (sino(:,v) - Ax - y(:,v)) ./ (sum_rows(rows) + 1);
      y(:,v) = y(:,v) + alg_params.alpha * diff;
    else
      diff = (sino(:,v) - Ax) ./ sum_rows(rows);
    end
    
    % backproject
    bp = (in_params.A(rows,:)' * diff) ./ sum_cols;
    
    % precondition
    if alg_params.precon
      bp = real(ifft2(fft2(reshape(bp,ny,nx)) ./ (prec)));
      bp = bp(:) / max(bp(:)) * 0.2;
    end
    
    % update
%     figure(1), imshow(rec',[]); colorbar;
    rec(:) = rec(:) + alg_params.alpha * bp; 
%     figure(2), imshow(rec', []); colorbar;
%     pause;
    
    if ~isempty(alg_params.min_val)
      rec = max(alg_params.min_val, rec);
    end

    % update
    pstart = pstart + ndet;
  end;
  
  
  times(it) = toc;
  snrs(it) = ma_snr(in_params.gt_vol, in_params.gt_vol - rec');
  iters(it) = it;
  
end

% transpose
rec = rec';

times = cumsum(times);
end

