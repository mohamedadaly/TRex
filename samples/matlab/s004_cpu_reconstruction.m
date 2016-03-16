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
colormap gray

vol_geom = astra_create_vol_geom(256, 256);
proj_geom = astra_create_proj_geom('parallel', 1.0, 384, ...
  linspace2(0,pi,30));
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
figure(1); imshow(P, []); axis image; axis off;
figure(2); imshow(sinogram, []); axis image; axis off;
% imshow(sinogram, []);

% compare full projection matrix
sino2 = reshape(proj_mat * reshape(P',[],1), 384, 30)';
disp(norm(sino2 - sinogram, 'fro'));

% add sinogram noise
rng(123);
sinogram = sinogram + randn(size(sinogram)) * 0.5;

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

% Set up the parameters for a reconstruction algorithm using the CPU
% The main difference with the configuration of a GPU algorithm is the
% extra ProjectorId setting.
cfg = astra_struct('BICAV');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ProjectorId = proj_id;
cfg.ProxInputDataId = prox_in_id;
cfg.option.UseMinConstraint = 1;
cfg.option.MinConstraintValue = 0;
cfg.option.UseMaxConstraint = 0;
cfg.option.MaxConstraintValue = 1;
cfg.option.ClearRayLength = 1;
cfg.option.Alpha = 1;
cfg.option.Lambda = 1e3;
cfg.option.ComputeIterationMetrics = 1;
cfg.option.GTReconstructionId = P_id;
cfg.option.IterationMetricsId = metrics_id;
cfg.option.ClearReconstruction = 1;

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
end


% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% Run 20 iterations of the algorithm
% This will have a runtime in the order of 10 seconds.
tic;
% astra_mex_algorithm('iterate', alg_id, 30*30);
astra_mex_algorithm('iterate', alg_id, 20);
toc;

% get the metrics
metrics = astra_mex_data2d('get', metrics_id)

% Get the result
rec = astra_mex_data2d('get', rec_id);
figure(3); imshow(rec, []); axis image; axis off; %, []);
fprintf('%s SNR=%f\n', cfg.type, snr(P, (P-rec)));

% Clean up. 
astra_mex_projector('delete', proj_id);
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sinogram_id);
astra_mex_data2d('delete', prox_in_id);
astra_mex_data2d('delete', P_id);

