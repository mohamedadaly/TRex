function [rec, times, snrs, iters] = ma_alg_bicav(in_params, alg_params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[times, snrs, iters] = deal(zeros(alg_params.iter, 1));

% take the transpose
in_params.prox_in = in_params.prox_in';
in_params.fbp = in_params.fbp';

% init in the transpose to make row-major
rec = zeros(size(in_params.gt_vol))';
[ny, nx] = size(rec);
if alg_params.init_fbp
  assert(all(size(in_params.fbp) == size(rec)));
  rec = in_params.fbp;
end


% take transpose to make row major because the projection matrix assumes
% row-major sinogram, and gives row-major volume
sino = in_params.sino';
[ndet, nviews] = size(sino);

% WLS solution for Poisson noise: multiply both sinogram and A by wi
if alg_params.wls
%   in_params.wi = sqrt(in_params.wi');
  in_params.wi = (in_params.wi .^ (1/6))';
  in_params.A = bsxfun(@(x,y)(x .* y), in_params.A, in_params.wi(:));
  sino = in_params.wi .* sino;
end

% initialize for proxisart
if alg_params.prox
  rec = in_params.prox_in;
  if alg_params.init_fbp, rec = rec + in_params.fbp; end
  y = zeros(size(sino));
  in_params.A = sqrt(2 * alg_params.lambda) * in_params.A;
  sino = sqrt(2 * alg_params.lambda) * sino;
end;

% sum of rows squares
sum_rows = sum(in_params.A .^ 2, 2);
sum_rows(abs(sum_rows)<1e-16) = 1;

% imitating ordered-subsets with more views per inner iteration
if isempty(alg_params.nsubsets), alg_params.nsubsets = nviews; end;
% how many views per subset
subset_size =floor(nviews / alg_params.nsubsets); 
% how many inner iterations
inner_iter = alg_params.nsubsets;

% iterations
for it = 1:alg_params.iter
  % loop on views
  tic;
  pstart = 1;
  for v = 1:inner_iter
    % rows of A
    arows = pstart : pstart + ndet * subset_size - 1;
    % cols of sino
    scols = (v - 1) * subset_size + 1 : v * subset_size;
    
    % sum over pixels for this view, where 1 for non-zero
    sum_cols = full(sum(abs(in_params.A(arows, :)) > 1e-16, 1)');
    
    % forward project current view
    Ax = full(reshape(in_params.A(arows, :) * rec(:), ndet, []));
    
    % difference
    if alg_params.prox
      diff = (sino(:,scols) - Ax - y(:,scols)) ./ (sum_rows(arows) + 1);
      y(:,scols) = y(:,scols) + alg_params.alpha * diff;
      % backproject
      bp = (in_params.A(arows,:)' * diff(:)) ./ sum_cols;

    else
      diff = (sino(:,scols) - Ax) ./ reshape(sum_rows(arows), ndet, []);
      % backproject
      bp = (in_params.A(arows,:)' * diff(:)) ./ sum_cols;
    end
        
%     figure(1), imshow(rec',[]); colorbar;
      % update
      rec(:) = rec(:) + alg_params.alpha * bp;
%     figure(2), imshow(rec', []); colorbar;
%     pause;
    
    if ~isempty(alg_params.min_val)
      rec = max(alg_params.min_val, rec);
    end

    % update
    pstart = pstart + ndet * subset_size;
  end;
  
  
  times(it) = toc;
  snrs(it) = ma_snr(in_params.gt_vol, in_params.gt_vol - rec');
  iters(it) = it;
  
end

% transpose
rec = rec';

times = cumsum(times);
end

