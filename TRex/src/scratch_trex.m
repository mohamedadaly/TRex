% -----------------------------------------------------------------------
% This file is part of the ASTRA Toolbox
% 
% Copyright: 2010-2015, iMinds-Vision Lab, University of Antwerp
%            2014-2015, CWI, Amsterdam
% License: Open Source under GPLv3
% Contact: astra@uantwerpen.be
% Website: http://sf.net/projects/astra-toolbox
% -----------------------------------------------------------------------
addpath matlab/tools/
addpath bin/x64/Release/
% colormap gray
N = 256;

vol_geom = astra_create_vol_geom(N, N);
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
P = phantom(N);
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

%

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
% extra ProjectorId setting.
cfg = astra_struct('TREX');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ProjectorId = proj_id;
% cfg.ProxInputDataId = prox_in_id;
cfg.option.UseMinConstraint = 1;
cfg.option.MinConstraintValue = 0;
cfg.option.UseMaxConstraint = 0;
cfg.option.MaxConstraintValue = 1;
cfg.option.Alpha = 1;
cfg.option.Lambda = 1e-3;
cfg.option.Rho = 25;
cfg.option.Mu = -1;
cfg.option.Sigma = 0.05;
cfg.option.WlsRoot = 1;
cfg.option.Data = 'WLS';
cfg.option.Prior = 'SAD'; 'ITV'; % 'ATV'
cfg.option.DataProx = 'SART-PROX';
cfg.option.InnerIter = 2;
cfg.option.ComputeIterationMetrics = 1;
cfg.option.GTReconstructionId = P_id;
cfg.option.IterationMetricsId = metrics_id;
cfg.option.ClearReconstruction = 1;
% cfg.option.PreconditionerId = -1; %prec_id;
% cfg.option.UseJacobiPreconditioner = 1;
% cfg.option.UseBSSART = 0;
% cfg.option.ProjectionOrder = 'sequential'; 'random';
% cfg.option.ProjectionOrderList = subset_start(30)-1;
astra_mex_data2d('set', wi_id, wi); %(wi.^(1/2))); %ones(size(wi)));
cfg.option.WlsWeightDataId =  wi_id;
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
astra_mex_algorithm('iterate', alg_id, 2);
toc;

% get the metrics
metrics = astra_mex_data2d('get', metrics_id)

% Get the result
rec = astra_mex_data2d('get', rec_id);
astra_mex_data2d('get', rec_id)

% figure(3); imshow(rec, []); axis image; axis off; %, []);
fprintf('%s SNR=%f\n', cfg.type, snr(P, (P-rec)));

%%
in_params = struct('vol_geom',vol_geom, 'proj_geom',proj_geom, ...
  'gt_vol',P, 'sino',sinogram, 'proj_id',proj_id, ...
  'wi',ones(size(sinogram)), 'fbp',fbp, 'A',proj_mat);
%
alg = 'admm';
alg_params = struct('iter',2, 'sigma',.05, 'rho',25, 'mu',[], ...1/(8*rho), ... 40&3 (no fbp) 100&5 (fbp)
  'theta',1, 'init_fbp',0, 'data','wls', 'prior','sad', ...
  'sigma_with_data', 0, 'prox','astra', 'prox_name','SART-PROX',...
  'wls_rt',1, ...
  'prox_params', struct('iter',2, ...
      'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
        'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
        'Alpha',1, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
        'ClearReconstruction',1)));
[rec, times, snrs, iters] = ma_alg_psart(alg, in_params, alg_params);
[times snrs iters]
