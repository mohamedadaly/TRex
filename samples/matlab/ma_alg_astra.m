function [rec, times, snrs, iters] = ma_alg_astra(alg, in_params, alg_params)
% MA_ALG_ASTRA a wrapper around running ASTRA methods.
%
% INPUT
% - alg: the name of the algorithm 'SART', 'ART', 'SIRT', 'FBP', 'BICAV',
%        'CGLS', 'PSART'
% 
% - in_params: the input paramters
%     vol_geom: the volume geometry
%     proj_geom: the projection geometry
%     gt_vol: the ground truth volume 
%     sino: the input projections
%     proj_id: the explicit projection matrix
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
%

% create GT id
gt_vol_id = astra_mex_data2d('create','-vol', in_params.vol_geom, ...
  in_params.gt_vol);
% create prox id
prox_in_id = astra_mex_data2d('create','-vol', in_params.vol_geom, ...
  in_params.prox_in);
% rec and copy in fbp for initialization (if available)
rec_id = astra_mex_data2d('create', '-vol', in_params.vol_geom, ...
  in_params.fbp);
% sino
sino_id = astra_mex_data2d('create', '-sino', in_params.proj_geom, ...
  in_params.sino);


% create object for metrics
metrics_id = astra_mex_data2d('create', '-vol', in_params.vol_geom);

% config
cfg = astra_struct(alg);

cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sino_id;
cfg.ProjectorId = in_params.proj_id;
cfg.ProxInputDataId = prox_in_id;

cfg.option = alg_params.option; 
cfg.option.GTReconstructionId = gt_vol_id;
cfg.option.IterationMetricsId = metrics_id;    

% Create the algorithm and run
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id, alg_params.iter);
astra_mex_algorithm('delete', alg_id);

% get the metrics
metrics = astra_mex_data2d('get', metrics_id);

% returns
times = metrics(:, 1);
snrs =  metrics(:, 2);
iters = (1:length(times))';

% Get rec
rec = astra_mex_data2d('get', rec_id);

% clear some volumes
astra_mex_data2d('delete', gt_vol_id);
astra_mex_data2d('delete', metrics_id);
astra_mex_data2d('delete', prox_in_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sino_id);

end