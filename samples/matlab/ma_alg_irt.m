function [rec, times, snrs, iters] = ma_alg_irt(alg, in_params, alg_params)
% MA_ALG_IRT a wrapper around running IRT methods.
% 
% INPUT
% - alg: the name of the algorithm 'admm', 'sqs-os', 'sqs-os-mom', 'pcg'
% 
% - in_params: the input paramters
%     vol_geom: the volume geometry
%     proj_geom: the projection geometry
%     gt_vol: the ground truth volume 
%     sino: the input projections
%     A: the explicit projection matrix
%     wi: the weights for the input projections
%     fbp: the initial FBP of the sino
%     
% - alg_params: the parameters for the algorithm, depends on the kind of 
%               algorithm used
%
% OUTPUT
% - rec: the reconstructed volume of the last iteration
% - times: the running times
% - snr: the SNR corresponding to iterations
% - iters: the inner iterations 
%
% Adapted from IRT toolbox by Fessler.
%

addpath samples/matlab/irt/

% init.
params = struct();
init();

switch alg
  % Ramani ADMM
  case 'admm'
    params.AL.mu = alg_params.mu; %0.001; %medwi;
    params.CG.nCG = alg_params.nCG; % 2;
    params.lambda = alg_params.lambda; % 0.02
%     params.minx = -Inf;
    params.AL.nu1 = nuoptapprox / params.AL.nu1AL1factor;
    params.CG.precon = alg_params.precon; % Do Precondition CG-solver inverting (A'A + mu*R'R) in ADMM

    %disp(sprintf('nu = %f', params.AL.nu1));
    if ~isempty(alg_params.nu), params.AL.nu1 = alg_params.nu; end

    params.figno = 0;
    kWAL1 = (mxW + params.AL.mu)/(mnW + params.AL.mu);
    params.AL.iWmu = 1 ./ (in_params.wi + params.AL.mu);
    params.maxitr = alg_params.iter; % 25; % Max # of outer iterations
    params.dispitrnum = 0;

    printm(['Doing ADMM with mu = ' num2str(params.AL.mu) '; nu1 = ' num2str(params.AL.nu1) '...'])
    [rec CADMM TADMM l2DADMM snrs] = runADMM(in_params.sino, ...
      params.xini, params);
    times = cumsum(TADMM); times(1) = [];
    snrs(1) = [];
    % inner iterations
    iters = (1:alg_params.iter) * alg_params.nCG;
    
  % Ramani MFISTA
  case 'mfista'
    params.figno = 0;
    params.doMFISTA = 1;
    params.maxitr = 20; % Max # of outer iterations
    params.AWy = params.A' * params.Wy;
    params.dispitrnum = 0;
    
    params.lambda = .1; %.02
    params.MFISTA.nCG = alg_params.nCG; % 5;

    [rec CMFIS TFIS l2DFIS snrs] = runMFISTA(in_params.sino, ...
      params.AWy, params.xini, params);
    times = cumsum(TFIS); times(1) = [];
    snrs(1) = [];
    
    % inner iterations
    iters = (1:alg_params.iter) * alg_params.nCG;    
  
  case {'sqs-os', 'sqs-os-mom', 'pcg'}
    % ids to convert from row-major to column-major in both input and output
    cids = reshape(reshape(1:nx*ny, nx,ny)',[],1);
    rids = reshape(reshape(1:nxs*nys, nxs,nys)',[],1);
    params.Ab = Gmatrix(in_params.A(rids, cids));

    % Regularizer
    params.R = [];
    if alg_params.reg
      % fail potential  wpot(t) = (1 + a * |t/d|) / (1 + b * |t/d|)
      params.R = Reg1(ones(size(in_params.gt_vol)), 'type_denom', 'matlab', ...
        'pot_arg',alg_params.pot_arg, 'beta',alg_params.beta); %1
  %         'pot_arg',{'gf1', 1, [1, 1]}, 'beta', 1); %1
      %   'pot_arg',{'gf1', 1, [1, 1]}, 'beta',0.1); 
      %   'pot_arg', {'hyper3', 1}, 'beta', 0.1);
      %   'pot_arg', {'huber', 1}, 'beta', 10);  %1
    end

    if strcmpi(alg, 'pcg')
      [xos timos] = pwls_pcg1(params.xini(:), params.Ab, ...
        Gdiag(in_params.wi), in_params.sino(:),...
         params.R, 'niter',alg_params.iter, 'isave', 'all');
    else      
      nsubsets = alg_params.nsubsets; %10; %10 or 15
      % convert to blocks
      params.Ab = Gblock(params.Ab, nsubsets);
      params.relax = alg_params.relax; % [1 1e-5]
      
      % SQS-OS
      [xos timos] = pwls_sqs_os(in_params.fbp, params.Ab, in_params.sino, ...
        params.R, 'wi', in_params.wi, 'niter', alg_params.iter, ...
        'pixmax', [0 inf], 'isave', 'all', 'chat',0, ...
        'relax0',alg_params.relax, 'mom', alg_params.mom);
    end;
    
    % compute SNR from return values
    snrs = compute_snr(xos(:,end - alg_params.iter + 1:end));
    % times
    times = cumsum(timos);    
    % rec
    rec = reshape(xos(:,end), ny, nx);    
    
    % inner iterations
    iters = (1:alg_params.iter);
    
  otherwise
    error(['Unknown algorithm: ' alg]);
end;


% adjust output
snrs = snrs(:);
times = times(:);
iters = iters(:);

  % compute SNR from the reconstruction in each iteration
  function snrs = compute_snr(xs)
    snrs = zeros(size(xs, 2), 1);
    for i=1:length(snrs)
      snrs(i) = ma_snr(in_params.gt_vol(:), in_params.gt_vol(:) - xs(:,i));
    end;    
  end

  % Stuff common to all gorithms.
  function init()
    params.clim = [0, 1];
    params.ig.mask = true(size(in_params.gt_vol));
    params.xini = in_params.fbp;

    mxW = max(in_params.wi(:));
    mnW = min(in_params.wi(:));
    params.Wy = in_params.sino;
    params.W = in_params.wi;

    params.Nmask = sum(params.ig.mask(:)>0);
    params.img = in_params.gt_vol;

    nx = in_params.vol_geom.GridColCount;
    ny = in_params.vol_geom.GridRowCount;
    nxs = in_params.proj_geom.DetectorCount;
    nys = length(in_params.proj_geom.ProjectionAngles);
    N = [nx ny];
    farg.A = in_params.A; %here
    farg.ny = ny;
    farg.nx = nx;
    farg.nys = nys;
    farg.nxs = nxs;
    params.A = Fatrix(size(in_params.A), farg, 'caller', 'malaa', ...
        'forw', @(arg,x)(reshape(arg.A * reshape(x',[],1), arg.nxs, arg.nys)'), ...
        'back', @(arg,x)(reshape(arg.A' * reshape(x',[],1), arg.nx, arg.ny)'));
    params.zoomr = 1:ny;
    params.zoomc = 1:nx;
    params.scale = 1;

    % Kappa for space-varying weights
    % if ~isvar('kappa'), printm('kappas')
      kappa = sqrt( div0(params.A' * in_params.wi, ...
        params.A' * ones(size(in_params.wi))) );
      skap = sort(kappa(:));
      skap = skap(skap > 0);
      kappa(kappa==0) = skap(1); % Ensure kappa does not have any zeros
      rw = kappa .^ 2; % spatially varying regularization weights; usually kappa ^ 2, but can be adjusted for quality
      im(kappa)
    % end

    params.scale = 1;
    params.kappa = kappa;
    params.rw = rw;

    % Regularization parameter for statistical recon
    lambda = 0.01; % medwi * 20; % heuristic!
    params.lambda = lambda;

    % Parameters for Preconditioned CG-solver inside ADMM for "inverting" A'A + nu * R'R
    params.CG.precon = 1; % Do Precondition CG-solver inverting (A'A + mu*R'R) in ADMM
    params.CG.restol = 1e-8; % Tol for residue
    params.CG.nCG = 5; % # of CGs
    % 	params.MFISTA.nCG = 5; % maximum # of CG iterations inside MFISTA

    % x_infinity solution
    params.xinf = zeros(nx,ny); % to be populated when x_infinity ( minimizer of the cost ) is available
    params.xinfnorm = 1;

    % Parameters for cost decrease, etc
    params.dcosttol = 1e18;
    params.dxtol = 0;

    % Display parameters
    params.dispfig = 0; % turn on to display recon at every iteration
    params.dispitr = 1; % turn on to display cost, etc., at every iteration
    params.dispitrnum = 3; % display info once these many iterations

    % ADMM parameters
    params.AL.mu = 1; % To be populated later
    params.AL.nu1 = 1; % To be populated later
    params.AL.nu1AL1factor = 100; % divide nu in (CAA + nu * R'R) by this factor

    params.AL.iWmu = [];
    params.AL.iRnu2nu1 = [];

    params.formatstringC = '%0.5E';
    params.formatstringT = '%0.3E';
    params.formatstringO = '%0.2E';

    params.Constants.EPSL = 1.0E-10; % Some constants needed for Brent iterations for minimization without derivatives; used for condition number minimization
    params.Constants.CGOLD = 0.3819660;
    params.Constants.GOLD = 1.618034;
    params.Constants.GLIMIT = 100.0;
    params.Constants.Brenttol = 0.01;
    params.Constants.maxBitr = 100;

    % Recon Setup
    % Potential Function
%     params.PriorType = 'l1'; % Absolute function = |t|
    % params.PriorType = 'FP'; % Fair potential = alpha^2 * ( |t| / alpha - log(1 + |t| / alpha ) )
    % params.Prior.alpha = 1e-4; % smoothing parameter for Fair potential

    % Type of reg op
    % params.Operator = 'AFD'; % plain finite differences (bilateral form)
    % params.Operator = 'FD'; % gradient-norm of finite differences; **** for TV use l1-FD ****
    % params.Operator = 'W'; % Undecimanted (shift-invariant) wavelet transform
    
    switch alg
      case {'admm', 'mfista'}
        params.PriorType = alg_params.prior_type;
        params.Operator = alg_params.operator;
      otherwise
        params.PriorType = 'l1';
        params.Operator = 'W';
    end;

    % Wavelet Options
    params.dwtmode = 'per'; % Period boundaries for wavelet implementation
    params.redundancy = 'undecimated'; % Undecimated wavelet transform using wavelet filters corresponding to standard orthonormal wavelets
    params.wname = 'haar'; % wavelet name, see waveinfo.m
    params.nlev = 2; % # of levels in wavelet transform
    params.includeApprox = false; % do not include wavelet approximation level in the reg op

    % Wavelet filters
    [lod, hid, lor, hir] = wfilters(params.wname);

    % Normalize filters so as to avoid a product with 0.5 during inverse undecimated wavelet transform
    params.lod = lod/sqrt(2);
    params.hid = hid/sqrt(2);
    params.lor = lor/sqrt(2);
    params.hir = hir/sqrt(2);
    % end


    % Regularization object
    if ~isvar('params.R'), printm('R'); end
    switch params.Operator
    case {'W'}
      % With mask supporting only the object, used in the recon algo
      R = Gwave2('mask', params.ig.mask, 'dwtmode', 'per', ...
        'redundancy', 'undecimated', 'wname', 'haar', ...
        'nlevel', params.nlev, 'includeApprox', false, ...
        'scale_filters', 1/sqrt(2));

      % Without mask, used only to obtain the freq.resp of Rf'Rf for preconditioner
      Rf = Gwave2('mask', true(N), 'dwtmode', 'per', ...
        'redundancy', 'undecimated', 'wname', 'haar', ...
        'nlevel', params.nlev, 'includeApprox', false, ...
        'scale_filters', 1/sqrt(2));

    case {'AFD'}
      R = fatrix2('imask', params.ig.mask, 'odim',[size(params.ig.mask) 2], ...
        'arg',[], 'abs',[], 'power',[], ...
        'forw', @(arg,x)(fdiff(x)), ...
        'back', @(arg,x)(ndiv(x)));

      Rf = fatrix2('imask',true(N), 'odim',[1 N], ...
        'arg',[], 'abs',[], 'power',[], ...
        'forw', @(arg,x)(fdiff(x)), ...
        'back', @(arg,x)(ndiv(x)));
    case { 'FD'}
      fail('Fatrix for finite differences should be inserted here...')

    otherwise
      fail('Wrong choice for regularization operator')
    end
    params.R = R;

    % preconditioner related
    ei = zeros(nx,ny); ei(nx/2+1,ny/2+1) = 1; RRei = R' * ( R * ei );
    RR = abs(fft2(RRei)); % freq. response corresponding to circulant R'R (uses Rf without masks)
    mxRR = max(RR(:));
    mnRR = min(RR(:));

    params.RR = RR;
    params.mxRR = mxRR;
    params.mnRR = mnRR;

    warn 'no file?'
    params.eigtol = eps; % Matlab epsilon
    params.eigpowermaxitr = 10000;
    params.eigdispitr = 10;
    params.N = N;

    % eigenvalue related
    printm('Computing max eigvalue of AWA using Power method for MFISTA...');
    tic
    mEAWA = get_MaxEigenValue_1by1(params, 'AWA'); % Find maximum eigenvalues of the forward system and the regularization
    toc
    printm(['Max eigvalue of AWA =' num2str(mEAWA)]);

    params.mEAWA = mEAWA;
    params.eiginflate = 1.0025;
    params.mEval_adjust = mEAWA * params.eiginflate;

    % Circulant approximation to A'A by taking response of A'A to an impulse at the center
    % if ~isvar('CAA'), printm('ADMM: Minimizing the condition number of CAA + nu1*RR');
      ei = zeros(nx,ny); ei(nx/2+1,ny/2+1) = 1; 
      AAei = params.A' * params.A * ei;
    %   AAei = params.A' * params.A * ei(:);
    %   AAei = reshape(AAei, ny, nx);
    CAA = abs( real( fft2( fftshift( AAei ) ) ) ); % freq. resp corresponding to the circulant approximation to A'A
    params.CAA = CAA;
    params.BRbx = 1e-2; %% Center of first bracket for Brent optimization of condition number

    [minCondCAARRapprox, nuoptapprox] = Brent_Linmin(params); %% Input arguments unused except for "params"; Need to simplify code using varargin!

    params.nuoptapprox = nuoptapprox; % approximate nu that minimizes condition number of original system (A'A + nu * R'R)
    params.CondCAARRapprox = minCondCAARRapprox; % approximate min.cond.num
    printm(['ADMM: Approx min cond num = ' num2str(minCondCAARRapprox, '%0.4E'), '; approx optimal nu = ' num2str(nuoptapprox, '%0.4E')]);

    % SPS the same as SQS
    % [xos timos] = pwls_sps_os(params.xini, dd.sinogram, wi, sos.Ab,  ...
    %   sos.R, 10, [0 inf], [], [], [1 1e-5], 1);  %[1 1e-5]

  end

end

% computes the first forward difference in both horizontal and vertical
% direction with Neumann boundary conditions i.e. zeros at the end
% returns the horizontal difference and vertical difference on third
% dimension
function [dx] = fdiff(x)
sz = size(x);
assert(length(sz) == 2);

% horizontal difference
dxh = zeros(sz);
dxh(1:end, 1:end-1) = x(1:end, 2:end) - x(1:end, 1:end-1);

% vertical
dxv = zeros(sz);
dxv(1:end-1, 1:end) = x(2:end, 1:end) - x(1:end-1, 1:end);

% stack in third dimension
dx = cat(3,dxh, dxv);
end

% computes the negative divergence (transpose of fdiff)
% returns a matrix with same columns and rows
function [x] = ndiv(dx)
sz = size(dx);
assert(length(sz) == 3);

% return matrix
x = zeros(sz(1:2));

% contribution of left with +ve
x(:, 2:end) = dx(:, 1:end-1,1);
% up with +ve
x(2:end, :) = x(2:end, :) + dx(1:end-1, 1:end, 2);
% down and right with -ve
x =  x - dx(:, 1:end, 1) - dx(:, 1:end, 2);
end
