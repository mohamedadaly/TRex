function [ output_args ] = ma_fig_periter( num_proj )
% MA_FIG_PERITER Plots performance per iteration for different phantoms.

% phantom and size
phan = 'sl';
phan_size = 512;
phan_file = sprintf('%s-%d.png', phan, phan_size);

% amount and type of noise
noise_type = 'gauss';
noise_level = 0.05;

% num of projections
num_proj = 30;

% iterations
iter = 40;

% figure file
pat = sprintf('per_iter-ph_%s-%d_nt_%s-nl_%.2f-np_%d', ...
  phan, phan_size, noise_type, noise_level, num_proj);
fig_file = [pat '.pdf'];
mat_file = [pat '.mat'];

% geometries
vol_geom = astra_create_vol_geom(phan_size, phan_size);
proj_geom = astra_create_proj_geom('parallel', 1.0, ceil(1.25 * phan_size), ...
  linspace2(0,pi,num_proj));
% projector
proj_id = astra_create_projector('linear', proj_geom, vol_geom);
proj_mat_id = astra_mex_projector('matrix', proj_id);
proj_mat = astra_mex_matrix('get', proj_mat_id);

% load phantom
P = imread(phan_file);
P = im2double(P);
% create gt astra data
P_id = astra_mex_data2d('create', '-vol', vol_geom, P);
% generate projections
[sinogram_id, sinogram] = astra_create_sino(P, proj_id);

% add sinogram noise
rng('default') 
rng(123);
switch noise_type
case 'guass'
  sinogram = sinogram + ...
    randn(size(sinogram)) * nosie_level * mean(sinogram(:));
case 'poisson'
  error('poisson');
end
% delete and recreate
astra_mex_data2d('delete', sinogram_id);
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogram);

% FBP initialization
rec_id = astra_mex_data2d('create', '-vol', vol_geom);
cfgi = astra_struct('FBP');
cfgi.ReconstructionDataId = rec_id;
cfgi.ProjectionDataId = sinogram_id;
cfgi.ProjectorId = proj_id;
alg_id = astra_mex_algorithm('create', cfgi);
astra_mex_algorithm('iterate', alg_id, 5);
astra_mex_algorithm('delete', alg_id);
% now initial value is stored in rec_id
fbp = astra_mex_data2d('get',rec_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', rec_id);

% Algorithms to run
algs = {
  % ASTRA types
  %
  struct('name','SART', 'type','astra', 'clr','-kd', ...
    'alg','SART', ...
    'alg_params', struct('iter',iter, ...
      'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
      'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
      'Alpha',1, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
      'ClearReconstruction',1, 'UseJacobiPreconditioner',0, ...
      'UseBSSART',0)))

  struct('name','BSSART', 'type','astra', 'clr','-.ko', ...
    'alg','SART', ...
    'alg_params', struct('iter',iter, ...
      'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
      'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
      'Alpha',1, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
      'ClearReconstruction',1, 'UseJacobiPreconditioner',0, ...
      'UseBSSART',1)))

  struct('name','ART', 'type','astra', 'clr','-rs', ...
    'alg','ART', ...
    'alg_params', struct('iter',iter, ...
      'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
      'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
      'Alpha',1, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
      'ClearReconstruction',1, 'UseJacobiPreconditioner',0, ...
      'UseBSSART',0)))

  struct('name','SIRT', 'type','astra', 'clr','-b+', ...
    'alg','SIRT', ...
    'alg_params', struct('iter',iter, ...
      'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
      'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
      'Alpha',1, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
      'ClearReconstruction',1, 'UseJacobiPreconditioner',0, ...
      'UseBSSART',0)))

  struct('name','BICAV', 'type','astra', 'clr','-g*', ...
    'alg','BICAV', ...
    'alg_params', struct('iter',iter, ...
      'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
      'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
      'Alpha',1, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
      'ClearReconstruction',1, 'UseJacobiPreconditioner',0, ...
      'UseBSSART',0)))

  struct('name','OS-SQS', 'type','astra', 'clr','-m^', ...
    'alg','OS-SQS', ...
    'alg_params', struct('iter',iter, ...
      'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
      'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
      'Alpha',1, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
      'ClearReconstruction',1, 'UseJacobiPreconditioner',0, ...
      'UseBSSART',0)))

  struct('name','CGLS', 'type','astra', 'clr','-cv', ...
    'alg','CGLS', ...
    'alg_params', struct('iter',iter, ...
      'option',struct('UseMinConstraint',1, 'MinConstraintValue',0, ...
      'UseMaxConstraint',0, 'MaxConstraintValue',1, ...
      'Alpha',2, 'Lambda',1e3, 'ComputeIterationMetrics',1, ...
      'ClearReconstruction',1, 'UseJacobiPreconditioner',0, ...
      'UseBSSART',0)))      
};

fprintf('Fig %s\n', fig_file);

% loop and compute
results = cell(size(algs));
for a = 1:length(algs)
  alg = algs{a};
  
  switch alg.type
  case 'astra'
    % input params
    in_params = struct('vol_geom',vol_geom, 'proj_geom',proj_geom, ...
      'gt_vol',P, 'sino',sinogram, 'proj_id',proj_id, ...
      'wi',ones(size(sinogram)), 'fbp',fbp, 'prox_in',zeros(size(P)));
    % run
    fprintf('Alg %s\n', alg.name);
    [rec, tt, ss, it] = ma_alg_astra(alg.alg, in_params, alg.alg_params);
  case 'irt'
  end
  
  % put results
  results{a}.rec = rec;
  results{a}.times = tt;
  results{a}.snrs = ss;
  results{a}.iter = it;
end

% save
save(mat_file, 'results', 'algs'); 

% plot
hfig = figure('Name',fig_file, 'Position',[1, 1, 800, 800]);
hold on;
legends = cell(size(results));
for a = 1:length(algs)
  res = results{a};
  alg = algs{a};
  
  plot(res.iter, res.snrs, alg.clr, 'LineWidth',2);
  legends{a} = alg.name;

end
hold off;
title(sprintf('Phantom %s', phan));
xlabel('Iteration');
ylabel('SNR (db)');
xlim([1 iter]);
legend(legends, 'Location','SouthEast');

% save figure
save_fig(fig_file, hfig, 'pdf');
end

