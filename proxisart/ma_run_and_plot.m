function [ output_args ] = ma_run_and_plot(plt, arg, algs )
% MA_FIG_PERITER Plots performance per iteration for different phantoms.

[fig_file, mat_file, vol_geom, proj_geom, P, sinogram, proj_id, ...
  proj_mat, fbp] = deal([]);

init;
eval(plt);

  % initialize sinogram, projection, noise, ... etc
  function init

    % phantom and size
    phan = arg.phan;
    phan_size = arg.phan_size;
    phan_file = sprintf('phantoms/%s-%d.mat', phan, phan_size);

    % amount and type of noise
    noise_type = arg.noise_type;
    noise_level = arg.noise_level;

    % num of projections
    num_proj = arg.num_proj;

    % iterations
    iter = arg.iter;

    % figure file
    pat = sprintf('%s/%s%s-ph_%s-%d_nt_%s-nl_%.3f-np_%d-p_%s', ...
      arg.path, arg.prefix, plt, phan, phan_size, noise_type, ...
      noise_level, num_proj, arg.proj_type);
    fig_file = [pat '.pdf'];
    mat_file = [pat '.mat'];

    % geometries
    vol_geom = astra_create_vol_geom(phan_size, phan_size);
    % num_det = ceil(1.1 * phan_size);
    % fan_angle = 30 * pi/180;
    % det_size = 2;
    % voxel_size = 1;
    % sdd = 0.5 * det_size * num_det / tan(fan_angle);
    % sid = 0.5 * voxel_size * phan_size / tan(fan_angle) * 1.1;
    switch arg.proj_type
    case 'parallel'
      num_det = ceil(1.1 * phan_size);
      proj_geom = astra_create_proj_geom('parallel', 1.0, num_det, ...
        linspace2(0,2*pi,num_proj));

      % projector
      proj_id = astra_create_projector('linear', proj_geom, vol_geom);

    case 'fan'
      det_size = 1.0239;
      num_det = 888;
      sdd = 949.075;
      sid = sdd * (phan_size/2) / (num_det/2);
      proj_geom = astra_create_proj_geom('fanflat', det_size, num_det, ...
        linspace2(0,2*pi,num_proj), sid, sdd - sid); %2*pi

      % projector
      proj_id = astra_create_projector('line_fanflat', proj_geom, vol_geom);

    case 'mouse'
      det_size = 4/3; 0.16176;
      num_det = 512;
      sdd = 529.59 * (4/3) / .16176;
      sid = sdd * (phan_size/2) / (num_det/2 * det_size); %395.730011;
      % projections
      projs = (0:194)*pi/180 + pi/2;
      projs_ids = 1:floor(195/num_proj):195;
      proj_geom = astra_create_proj_geom('fanflat', det_size, num_det, ...
        projs(projs_ids), sid, sdd - sid); %2*pi

      % projector
      proj_id = astra_create_projector('line_fanflat', proj_geom, vol_geom);

    end
    proj_mat_id = astra_mex_projector('matrix', proj_id);
    proj_mat = astra_mex_matrix('get', proj_mat_id);
    astra_mex_matrix('delete', proj_mat_id);

    % load phantom
    switch phan
    case {'sl', 'mod-sl', 'ncat'}
      P = load(phan_file);
      P = P.P;

      % generate projections
      [sinogram_id, sinogram] = astra_create_sino(P, proj_id);
      astra_mex_data2d('delete', sinogram_id);

    case 'mouse'
      pp = load(phan_file);
      P = pp.P;
      sinogram = pp.sino;
      % choose the relevant projections
      sinogram = sinogram(projs_ids, :);
    end

    % create gt astra data
    P_id = astra_mex_data2d('create', '-vol', vol_geom, P);

    % add sinogram noise
    rng('default') 
    rng(123);
    switch noise_type
    case 'gauss'
      sinogram = sinogram + ...
        randn(size(sinogram)) * noise_level * max(sinogram(:));
    case 'poisson'
      error('poisson');
    otherwise
      error(['Unsupported nosie: ' noise_type]);
    end
    % add negative value to make all positive
    sinogram = sinogram - min(sinogram(:));

    % delete and recreate
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

    % clear
    astra_mex_data2d('delete', P_id);
    astra_mex_data2d('delete', rec_id);
    astra_mex_data2d('delete', sinogram_id);
    
  end

  % Plot SNR per iteration
  function per_iter
    fprintf('Fig %s\n', fig_file);

    % loop and compute
    if arg.recompute
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
          alg.alg_params.iter = arg.iter;

          [rec, tt, ss, it] = ma_alg_astra(alg.alg, in_params, ...
            alg.alg_params);
        case 'irt'
        end

        % put results
      %   results{a}.rec = rec;
        results{a}.times = tt;
        results{a}.snrs = ss;
        results{a}.iter = it;
      end

      % save
      save(mat_file, 'results', 'algs'); 
    else
      rr = load(mat_file);
      results = rr.results;
      algs = rr.algs;
    end

    % plot
    hfig = figure('Name',fig_file, 'Position',[1, 1, 800, 800]);
    hold on;
    legends = cell(size(results));
    for a = 1:length(algs)
      res = results{a};
      alg = algs{a};

      plot(res.iter, res.snrs, 'Color',alg.clr, ...
        'LineStyle',alg.lstyle, 'Marker', alg.marker, ...
        'LineWidth',2);
      legends{a} = alg.name;
    end
    hold off;
%     title(sprintf('Phantom %s', phan));
    xlabel('Iteration');
    ylabel('SNR (db)');
    xlim([1 arg.iter]);
    legend(legends, 'Location','SouthEast');

    % save figure
    save_fig(fig_file, hfig, 'pdf');
  end

  % Plot SNR per number of projections
  function per_num_proj
  end
end

