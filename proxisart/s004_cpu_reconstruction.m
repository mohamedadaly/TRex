% -----------------------------------------------------------------------
% This file is part of the ASTRA Toolbox
% 
% Copyright: 2010-2015, iMinds-Vision Lab, University of Antwerp
%            2014-2015, CWI, Amsterdam
% License: Open Source under GPLv3
% Contact: astra@uantwerpen.be
% Website: http://sf.net/projects/astra-toolbox
% -----------------------------------------------------------------------
% addpath matlab/tools/
% addpath bin/x64/Release/
colormap gray

vol_geom = astra_create_vol_geom(256, 256);
proj_geom = astra_create_proj_geom('parallel', 1.0, 384, ...
  linspace2(0,2*pi,30));
% vol_geom = astra_create_vol_geom(4,4);
% proj_geom = astra_create_proj_geom('parallel', 1.0, 10, linspace2(0,pi,4));

% For CPU-based algorithms, a "projector" object specifies the projection
% model used. In this case, we use the "strip" model.
proj_id = astra_create_projector('linear', proj_geom, vol_geom);
% assumes the volume is a vector in row-major order, and the output
% is also in row-major order
% to emulate projection, do the following:
%   sino = reshape(proj_mat * reshape(P',[],1), 384, 30)';
% which should be the same as sinogram
proj_mat_id = astra_mex_projector('matrix', proj_id);
proj_mat = astra_mex_matrix('get', proj_mat_id);

% save for input to IRT toolbox for comparison
% save sl.mat P proj_mat sinogram vol_geom proj_geom fbp

% Create a sinogram from a phantom
P = phantom(256);
% create gt astra data
P_id = astra_mex_data2d('create', '-vol', vol_geom, P);

% P = phantom(4);
[sinogram_id, sinogram] = astra_create_sino(P, proj_id);
% figure(1); imshow(P, []); axis image; axis off;
% figure(2); imshow(sinogram, []); axis image; axis off;
% imshow(sinogram, []);

% % compare full projection matrix
% sino2 = reshape(proj_mat * reshape(P',[],1), 384, 30)';
% disp(norm(sino2 - sinogram, 'fro'));

% add sinogram noise
rng('default') 
rng(123);
% sinogram = sinogram + randn(size(sinogram)) * 0.05 * mean(sinogram(:));
% sinogram = sinogram .* (randn(size(sinogram)) * .1 + 1);
      
% copied from Fessler's IRT
I0 = 1e5;  %1e5 % incident photons; decrease this for "low dose" scans
poissonfactor = 0.4; % for generating poissonn noise using rejection method
% scale the sinogram to have a max of 10, then scale again after
% sampling
scale_factor = 10 / max(sinogram(:));
yi = poisson(I0 * exp(-sinogram * scale_factor), 0, 'factor', poissonfactor); % poisson noise for transmission data:
if any(yi(:) == 0)
  warn('%d of %d values are 0 in sinogram!', sum(yi(:)==0), length(yi(:)));
end
sinogram = log(I0 ./ max(yi,1)) / scale_factor; % noisy fan-beam sinogram: form of the data in the quadratic approximation to the actual log-likelihood
% PWLS weights
wi = yi / max(yi(:)); %
% wi
wi_id = astra_mex_data2d('create', '-sino', proj_geom, wi);

astra_mex_data2d('delete', sinogram_id);

% We now re-create the sinogram data object as we would do when loading
% an external sinogram
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogram);

% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom);

% create a data object for the proximal input
prox_in = zeros(size(P));
prox_in_id = astra_mex_data2d('create', '-vol', vol_geom, prox_in);

% create object for metrics
metrics_id = astra_mex_data2d('create', '-vol', vol_geom);

%%

cfgi = astra_struct('FBP');
cfgi.ReconstructionDataId = rec_id;
cfgi.ProjectionDataId = sinogram_id;
cfgi.ProjectorId = proj_id;
alg_id = astra_mex_algorithm('create', cfgi);
astra_mex_algorithm('iterate', alg_id, 5);
astra_mex_algorithm('delete', alg_id);
% now initial value is stored in rec_id
fbp = astra_mex_data2d('get',rec_id);

%%
% addpath samples/matlab

in_params = struct('vol_geom',vol_geom, 'proj_geom',proj_geom, ...
  'gt_vol',P, 'sino',sinogram, 'A',proj_mat, ...
  'wi',wi, 'fbp',fbp);
%%
alg = 'admm'; %.0001&[]&2e3
alg_params = struct('iter',10, 'mu',.0001, 'nCG',2, 'lambda',.05, ... .001&.02
  'nu',2e3, 'operator','W', 'prior_type','l1', ...
  'init_fbp',1, 'precon',1, 'minx',-Inf); %nu=200
% alg_params = struct('iter',10, 'mu',0.0001, 'nCG',2, 'lambda',.1, ... .001&.01
%   'nu', 600, 'operator','AFD', 'prior_type','l1', 'precon',1, ... %nu=200
%   'init_fbp',1, 'minx',-Inf);

% profile on      
[rec, times, snrs, iters] = ma_alg_irt(alg, in_params, alg_params);
% profile viewer
[times snrs iters]
% figure(1), imshow(rec,[])

%%
alg = 'mfista';
alg_params = struct('iter',10, 'nCG',2, 'lambda',0.1, ...
  'prior_type','l1', 'operator','AFD', 'init_fbp',1, 'precon',0);  
%%
alg = 'pcg';
precon = spdiag(1 ./ sum(in_params.A .^ 2, 1));
alg_params = struct('iter',10, 'beta',0, 'pot_arg',{{'gf1',1,[1 1]}}, ...
  'reg', 1, 'init_fbp',0, 'pixmax',[-Inf Inf], 'precon',precon);  
%%
alg = 'sqs-os';
alg_params = struct('iter',10, 'beta',1, 'pot_arg',{{'gf1',1,[1 1]}}, ...
  'mom',0, 'relax',[1 0], 'nsubsets',30, 'reg',0, 'bit_reverse',0, ...
  'init_fbp',0);  %[1 1e-5]
%%
alg = 'sqs-os-mom';
alg_params = struct('iter',10, 'beta',.01, 'pot_arg',{{'huber', 1}}, ...
  'mom',2, 'relax',[1 1e-5], 'nsubsets',10, 'reg',0);  
%%
%     % Regularizer
%     % fair potential  wpot(t) = (1 + a * |t/d|) / (1 + b * |t/d|)
%     params.R = Reg1(ones(size(in_prams.gt_vol)), 'type_denom', 'matlab', ...
%       'pot_arg',alg_params.pot_arg, 'beta',alg_params.beta); %1
% %         'pot_arg',{'gf1', 1, [1, 1]}, 'beta', 1); %1
%     %   'pot_arg',{'gf1', 1, [1, 1]}, 'beta',0.1); 
%     %   'pot_arg', {'hyper3', 1}, 'beta', 0.1);
%     %   'pot_arg', {'huber', 1}, 'beta', 10);  %1
%     params.relax = alg_params.relax; % [1 1e-5]


[rec, times, snrs, iters] = ma_alg_irt(alg, in_params, alg_params);
[times snrs iters]

%%
in_params = struct('vol_geom',vol_geom, 'proj_geom',proj_geom, ...
  'gt_vol',P, 'sino',sinogram, 'proj_id',proj_id, ...
  'wi',ones(size(sinogram)), 'fbp',fbp, 'prox_in',zeros(size(P)));

alg = 'CGLS';
alg_params = struct('iter',40, ...
    'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
      'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
      'Alpha',2, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
      'ClearReconstruction',1));
      
[rec, times, snrs, iters] = ma_alg_astra(alg, in_params, alg_params);
[times snrs iters]

%% PSART
in_params = struct('vol_geom',vol_geom, 'proj_geom',proj_geom, ...
  'gt_vol',P, 'sino',sinogram, 'proj_id',proj_id, ...
  'wi',ones(size(sinogram)), 'fbp',fbp, 'A',proj_mat);
%%
alg = 'admm';
alg_params = struct('iter',5, 'sigma',25, 'rho',10, 'mu',[], ...1/(8*rho), ... 40&3 (no fbp) 100&5 (fbp)
  'theta',1, 'init_fbp',1, 'data','l2', 'prior','atv', ...
  'sigma_with_data', 1, 'prox','astra', ...
  'psart_params', struct('iter',2, ...
      'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
        'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
        'Alpha',1, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
        'ClearReconstruction',1)));
[rec, times, snrs, iters] = ma_alg_psart(alg, in_params, alg_params);
[times snrs iters]

%% ADMM with ma_alg_sart SART
alg = 'admm';
alg_params = struct('iter',2, 'sigma',.25, 'rho',10.0, 'mu',[], ...1/(8*rho), ... 40&3 (no fbp) 100&5 (fbp)
  'theta',1, 'init_fbp',1, 'data','l2', 'prior','atv', ...
  'sigma_with_data', 1, 'prox','matlab', ...
  'psart_params', struct('iter',2, 'alpha',1, 'min_val',0, 'precon',0, ...
    'wls',0, 'BSSART',0, 'prox','sart', 'lambda',1e7, 'nsubsets',30, ...
    'sqs',0, 'init_fbp',0));
% struct('iter',2, ...
%       'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
%         'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
%         'Alpha',2, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
%         'ClearReconstruction',1)));


[rec, times, snrs, iters] = ma_alg_psart(alg, in_params, alg_params);
[times snrs iters]

%% ADMM with ma_alg_sart SQS
alg = 'admm';
alg_params = struct('iter',2, 'sigma',.25, 'rho',.10, 'mu',[], ...1/(8*rho), ... 40&3 (no fbp) 100&5 (fbp)
  'theta',1, 'init_fbp',0, 'data','l2', 'prior','atv', ...
  'sigma_with_data', 1, 'prox','matlab', ...
  'psart_params', struct('iter',2, 'alpha',1, 'min_val',0, 'precon',0, ...
    'wls',0, 'BSSART',0, 'prox','sqs', 'lambda',1e7, 'nsubsets',30, ...
    'sqs',1, 'init_fbp',0, 'nu',1));
% struct('iter',2, ...
%       'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
%         'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
%         'Alpha',2, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
%         'ClearReconstruction',1)));


[rec, times, snrs, iters] = ma_alg_psart(alg, in_params, alg_params);
[times snrs iters]
      
%%
alg = 'cp';
alg_params = struct('iter',20, 'sigma',100, 'lambda',0.03, 'mu',[], ... 1/(8*.010), ...
  'theta',1, 'init_fbp',0, 'data','l2', 'prior','atv', ...
  'psart_params', struct('iter',2, ...
      'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
        'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
        'Alpha',2, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
        'ClearReconstruction',1)));

profile on      
[rec, times, snrs, iters] = ma_alg_psart(alg, in_params, alg_params);
profile viewer
[times snrs iters]


%% Matrix Experimental SART
in_params = struct('vol_geom',vol_geom, 'proj_geom',proj_geom, ...
  'gt_vol',P, 'sino',sinogram, 'proj_id',proj_id, ...
  'wi',wi, 'fbp',fbp, 'A',proj_mat, 'prox_in',fbp*0);

%%    SART
% --> debug sqs
alg_params = struct('iter',5, 'alpha',1, 'min_val',0, 'precon',0, ...
  'wls',1, 'BSSART',0, 'prox','sart', 'lambda',1e0, 'nsubsets',30, 'sqs',0, ...
  'init_fbp',0);
[rec, times, snrs, iters] = ma_alg_sart(in_params, alg_params);
[times snrs iters]

%% BICAV
alg_params = struct('iter',5, 'alpha',1, 'min_val',0, 'precon',0, ...
  'wls',1, 'prox',1, 'lambda',1e0, 'nsubsets',30, 'init_fbp',0);
[rec, times, snrs, iters] = ma_alg_bicav(in_params, alg_params);
[times snrs iters]

%%    SQS
% --> debug sqs
alg_params = struct('iter',5, 'alpha',1, 'min_val',0, 'precon',0, ...
  'wls',0, 'BSSART',0, 'prox','sqs', 'lambda',1e0, 'nsubsets',30, 'sqs',1, ...
  'init_fbp',0, 'nu',1);
[rec, times, snrs, iters] = ma_alg_sart(in_params, alg_params);
[times snrs iters]

%%    SIRT
alg_params = struct('iter',10, 'alpha',2, 'min_val',0, 'precon',0, ...
  'gamma', .08);
[rec, times, snrs, iters] = ma_alg_sirt(in_params, alg_params);
[times snrs iters]

% figure(1), imshow(rec, [])
%%

% % preconditioner
% prec = ones(size(P));
% prec(128,128)=1;
% prec = reshape(proj_mat' * (proj_mat * prec(:)), 256, 256);
% % prec(prec < 1e-16) = 1;
% % prec = 1;
% prec_id = astra_mex_data2d('create', '-vol', vol_geom, prec);

% Set up the parameters for a reconstruction algorithm using the CPU
% The main difference with the configuration of a GPU algorithm is the
% extra ProjectorId setting.
cfg = astra_struct('SART');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ProjectorId = proj_id;
cfg.ProxInputDataId = prox_in_id;
cfg.option.UseMinConstraint = 1;
cfg.option.MinConstraintValue = 0;
cfg.option.UseMaxConstraint = 0;
cfg.option.MaxConstraintValue = 1;
cfg.option.Alpha = 2;
cfg.option.Lambda = 1e3;
cfg.option.ComputeIterationMetrics = 1;
cfg.option.GTReconstructionId = P_id;
cfg.option.IterationMetricsId = metrics_id;
cfg.option.ClearReconstruction = 1;
cfg.option.PreconditionerId = -1; %prec_id;
cfg.option.UseJacobiPreconditioner = 1;
cfg.option.UseBSSART = 0;
cfg.option.ProjectionOrder = 'sequential'; 'random';
cfg.option.ProjectionOrderList = subset_start(30)-1;
astra_mex_data2d('set', wi_id, sqrt(wi)); %ones(size(wi)));
cfg.option.WlsWeightDataId = -1;wi_id;
% cfg.option.ProjectionOrder = 'random';

% Available algorithms:
% ART, SART, SIRT, CGLS, FBP
 
% initialize with FBP
if cfg.option.ClearReconstruction == 0
  cfgi = cfg;
  cfgi.type = 'FBP';
  alg_id = astra_mex_algorithm('create', cfgi);
  astra_mex_algorithm('iterate', alg_id, 5);
  astra_mex_algorithm('delete', alg_id);
  % now initial value is stored in rec_id
  fbp = astra_mex_data2d('get',rec_id);
end

% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% Run 20 iterations of the algorithm
% This will have a runtime in the order of 10 seconds.
tic;
% astra_mex_algorithm('iterate', alg_id, 30*30);
astra_mex_algorithm('iterate', alg_id, 5);
toc;

% get the metrics
metrics = astra_mex_data2d('get', metrics_id)

% Get the result
rec = astra_mex_data2d('get', rec_id);
astra_mex_data2d('get', rec_id)

% figure(3); imshow(rec, []); axis image; axis off; %, []);
fprintf('%s SNR=%f\n', cfg.type, snr(P, (P-rec)));

%%
% Clean up. 
astra_mex_projector('delete', proj_id);
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sinogram_id);
astra_mex_data2d('delete', prox_in_id);
astra_mex_data2d('delete', P_id);

