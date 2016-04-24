  function [xs, info] = pwls_sqs_os_mom(x, yi, wi, Ab, R, niter, pixmax, ...
		denom, aai, relax0, conv, mom, ig, roi, chat)
%|function [xs, info] = pwls_sqs_os_mom(x, yi, wi, Ab, R, niter, pixmax, ...
%|		denom, aai, relax0, conv, mom, ig, roi, chat)
%|
%| penalized weighted least squares estimation/reconstruction
%| using separable quadaratic surrogates algorithm with
%| (optionally relaxed) ordered subsets.  (relaxation ensures convergence.)
%| + Nesterov's fast gradient method
%| + Kim and Fessler's optimal gradient method [kim::ofo]
%|
%| cost(x) = (y-Gx)' W (y-Gx) / 2 + R(x)
%|
%| in
%|	x	[np 1]		initial estimate
%|	yi	[ns nt na]	measurements (noisy sinogram)
%|	wi	[ns nt na]	weighting sinogram (or [] for uniform)
%|	Ab	[nd np]		Gblock object, aij >= 0 required!
%|					or sparse matrix (implies nsubset=1)
%|	R			penalty object (see Robject.m)
%|	niter			# of iterations (including 0)
%|
%| optional
%|	pixmax	[1] or [2]	max pixel value, or [min max] (default [0 inf])
%|	denom	[np 1]		precomputed denominator
%|	aai	[ns nt na]	precomputed row sums of |Ab|
%|	relax0	[1] or [2]	relax0 or (relax0, relax_rate)
%|	conv	[np 1]		converged image
%|	mom	[1]		0 if no, 
%|				1 if Nesterov's, 
%|				2 if Kim and Fessler's,
%| 	ig			image_geom
%|	roi			iz roi index
%|	chat
%|
%| out
%|	xs	[np 1]		last iterate
%|	info	[niter 2]	time
%|
%| Copyright 2002-2-12, Jeff Fessler, University of Michigan
%| Copyright 2015-8-13, Donghwan Kim, University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

cpu etic
info = zeros(niter,1);

Ab = block_op(Ab, 'ensure'); % make it a block object (if not already)
nblock = block_op(Ab, 'n');
starts = subset_start(nblock);

if ~isvar('conv')       || isempty(conv),       comp_rmsd = 0;
else    comp_rmsd = 1; info = zeros(niter,2);
end
if ~isvar('mom')	|| isempty(mom),	mom = 0;	end

if ~isvar('niter')	|| isempty(niter),	niter = 1;	end
if ~isvar('pixmax')	|| isempty(pixmax),	pixmax = inf;	end
if ~isvar('chat')	|| isempty(chat),	chat = false;	end
if isempty(wi)
	wi = ones(size(yi));
end
if ~isvar('aai') || isempty(aai)
	aai = reshape(sum(Ab'), size(yi)); % a_i = sum_j |a_ij|
					% requires real a_ij and a_ij >= 0
end

if ~isvar('relax0') || isempty(relax0)
	relax0 = 1;
end
if length(relax0) == 1
	relax_rate = 0;
elseif length(relax0) == 2
	relax_rate = relax0(2);
	relax0 = relax0(1);
else
	error relax
end

if length(pixmax) == 2
	pixmin = pixmax(1);
	pixmax = pixmax(2);
elseif length(pixmax) == 1
	pixmin = 0;
else
	error pixmax
end

%
% likelihood denom, if not provided
%
if ~isvar('denom') || isempty(denom)
	denom = Ab' * col(aai .* wi);	% requires real a_ij and a_ij >= 0
end, clear aai

if ~isvar('R') || isempty(R)
	pgrad = 0;		% unregularized default
	Rdenom = 0;
end

[ns nt na] = size(yi);


%
% loop over iterations
%

xs = zeros(numel(x), 1); %niter);
x = max(x,pixmin);
x = min(x,pixmax);
%xs(:,1) = x;

if (comp_rmsd)
        iter = 1;
        xtmp = ig.embed(x);
        ctmp = ig.embed(conv);
        info(iter,2) = sqrt(mean(col(xtmp(:,:,roi) - ctmp(:,:,roi)).^2));
        disp(sprintf('rmsd = %0.5g', info(iter,2)));
end

% for mom
if (mom == 1)
	tprev = 1; xprev = x;
elseif (mom == 2)
	tprev = 1; xprev = x; xprev0 = x;
end

for iter = 2:niter
	ticker(mfilename, iter, niter)

	relax = relax0 / (1 + relax_rate * (iter-2));

	%
	% loop over subsets
	%
	for iset=1:nblock
		iblock = starts(iset);
		ia = iblock:nblock:na;

		li = Ab{iblock} * x;
		li = reshape(li, [ns nt length(ia)]);
		resid = wi(:,:,ia) .* (yi(:,:,ia) - li);
		grad = Ab{iblock}' * resid(:); % G' * W * (y - G*x)

		if ~isempty(R)
			pgrad = col(R.cgrad(R, ig.embed(x)));
			%Rdenom = R.denom(R, x);
			Rdenom = col(R.denom_sqs1(R, zeros(size(ig.embed(x)))));
		end

		num = nblock * grad - pgrad;
		den = denom + Rdenom;
		den(denom == 0) = 0; %		

		x = x + relax * div0(num, den);     % relaxed update
    x = max(x,pixmin);              % lower bound
    x = min(x,pixmax);              % upper bound
		xs = x; %%%%%%%%%%
		
		if (mom == 1) 
			t = 1/2*(1 + sqrt(1+4*tprev^2));
      xtmp = x;
      x = x + (tprev - 1)/t * (x - xprev);
      xprev = xtmp;
      tprev = t;
    elseif (mom == 2)
			t = 1/2*(1 + sqrt(1+4*tprev^2));
      xtmp = x;
      x = x + (tprev - 1)/t * (x - xprev) ...
        + tprev/t * (x - xprev0);
      xprev = xtmp;
			xprev0 = x;
      tprev = t;
		end
	        
		if (comp_rmsd)
                	xtmp = ig.embed(xs);
                	ctmp = ig.embed(conv);
                	rmsd = sqrt(mean(col(xtmp(:,:,roi) - ctmp(:,:,roi)).^2));
                	disp(sprintf('%d.%d of %d.%d: rmsd = %0.5g', ...
				iter, iset, niter, nblock, rmsd));
    end
	end

	if chat, printm('range(x) = %g %g', min(x), max(x)), end
	%xs(:,iter) = x;
	info(iter,1) = cpu('etoc');
	if (comp_rmsd)
		xtmp = ig.embed(xs);
		ctmp = ig.embed(conv);
                info(iter,2) = sqrt(mean(col(xtmp(:,:,roi) - ctmp(:,:,roi)).^2));
                disp(sprintf('rmsd = %0.5g', info(iter,2)));
        end
end

