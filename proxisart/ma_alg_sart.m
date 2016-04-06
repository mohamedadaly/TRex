function [rec, times, snrs, iters] = ma_alg_sart(in_params, alg_params)
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
  in_params.wi = sqrt(in_params.wi');
  in_params.A = bsxfun(@(x,y)(x .* y), in_params.A, in_params.wi(:));
  sino = in_params.wi .* sino;
end

% precondition
if alg_params.precon
  prec = zeros(size(rec));
  prec(floor(ny/2)+1, floor(nx/2)+1) = 1;
  prec = reshape(in_params.A' * (in_params.A * prec(:)), ny, nx);
  prec = abs(fft2(prec));
  
  xx = load('pre.mat');
  prec = (xx.p1 + 5e3 * xx.p2);
  
end

% % scale matrix and sinogram by sum_rows to solve a least-square instead of
% % a weighted least-square
% if alg_params.wls
%   % sum of rows
%   sum_rows = sum(in_params.A, 2);
%   sum_rows(abs(sum_rows)<1e-16) = 1;
%   % scale
%   in_params.A = bsxfun(@(a,b)(a.*b), in_params.A, sum_rows);
%   sino = reshape(sum_rows .* sino(:), ndet, nviews);
% end

% initialize for proxisart
switch alg_params.prox
  case 'sart'
    rec = in_params.prox_in;
    if alg_params.init_fbp, rec = rec + in_params.fbp; end
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
%   rec = reshape(in_params.A(floor(ndet * nviews/2),:), ny, nx);
end

% imitating ordered-subsets with more views per inner iteration
if isempty(alg_params.nsubsets), alg_params.nsubsets = nviews; end;
% how many views per subset
subset_size =floor(nviews / alg_params.nsubsets); 
% how many inner iterations
inner_iter = alg_params.nsubsets;

% OS-SQS
if alg_params.sqs 
  % set BSSART to 1 so that we don't update column sums
  alg_params.BSSART = 1;
  % set row sums to 1 to neutralize
  sum_rows(:) = 1;
  % compute the column sums which are sum(A' * A)
  sum_cols = in_params.A' * sum(in_params.A, 2);
  
  % prox operator?
  if strcmp(alg_params.prox, 'sqs')
    % scale sum_cols by 2 lambda and add 1
    sum_cols = sum_cols * 2 * alg_params.lambda + alg_params.nu;
  end
  
  % scale
  sum_cols = sum_cols / alg_params.nsubsets;
end

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
    
    if ~alg_params.BSSART
      % sum over pixels for this view
      sum_cols = full(sum(in_params.A(arows, :), 1)');
%       sum_cols(abs(sum_cols) < 1e-16) = 1;
    end
    
    % forward project current view
    Ax = full(reshape(in_params.A(arows, :) * rec(:), ndet, []));
    
    % difference
    switch (alg_params.prox)
      case 'sart'
        diff = (sino(:,scols) - Ax - y(:,scols)) ./ (sum_rows(arows) + 1);
        y(:,scols) = y(:,scols) + alg_params.alpha * diff;
        % backproject
        bp = (in_params.A(arows,:)' * diff(:)) ./ sum_cols;

      case 'sqs'
        diff = (sino(:,scols) - Ax) ./ reshape(sum_rows(arows), ndet, []);          
        % backproject and add
        bp = 2 * alg_params.lambda * (in_params.A(arows,:)' * diff(:));
        % subset from volume to update 
        vsubset = ceil(length(rec(:)) / inner_iter);
        vids = (vsubset * (v-1) + 1 : min(length(rec(:)),vsubset*v))';
        bp(vids) = bp(vids) + (in_params.prox_in(vids) - rec(vids)) * alg_params.nu;
        % divide
        bp = bp ./ sum_cols;

      otherwise
        diff = (sino(:,scols) - Ax) ./ reshape(sum_rows(arows), ndet, []);
        % backproject
        bp = (in_params.A(arows,:)' * diff(:)) ./ sum_cols;
    end
    
    
    % precondition
    if alg_params.precon
      bp = real(ifft2(fft2(reshape(bp,ny,nx)) ./ (prec)));
      bp = bp(:) / max(bp(:)) * 0.2;
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

