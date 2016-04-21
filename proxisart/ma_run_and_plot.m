function [ output_args ] = ma_run_and_plot(plt, arg, algs )
% MA_FIG_PERITER Plots performance per iteration for different phantoms.

[fig_file, mat_file, vol_geom, proj_geom, P, sinogram, proj_id, ...
  proj_mat, fbp, wi] = deal([]);

% init;
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
%       proj_id = astra_create_projector('strip_fanflat', proj_geom, vol_geom);

    case 'mouse'
      det_size = 4/3; 0.16176;
      num_det = 512;
      sdd = 529.59 * (4/3) / .16176;
      sid = sdd * (phan_size/2) / (num_det/2 * det_size); %395.730011;
      % projections
      projs = (0:194)*pi/180 + pi/2;
      projs_ids = 1:floor(195/num_proj):195;
      if isfield(arg, 'force_num_proj') && arg.force_num_proj
        projs_ids = projs_ids(1:num_proj);
      end
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
      wi = [];

    case 'mouse'
      pp = load(phan_file);
      P = pp.P;
      sinogram = pp.sino;
      % choose the relevant projections
      sinogram = sinogram(projs_ids, :);
      
      % get the weights and scale
      wi = pp.wi(projs_ids, :);
      wi = wi / max(wi(:));
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
      % copied from Fessler's IRT
      I0 = noise_level;  %1e5 % incident photons; decrease this for "low dose" scans
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
      wi = yi / max(yi(:)); % PWLS weights; gives 0 weight to any ray where yi=0!
      mxW = max(wi(:));
      mnW = min(wi(:));      
    otherwise
      error(['Unsupported nosie: ' noise_type]);
    end
    % add negative value to make all positive
%     sinogram = sinogram - min(sinogram(:));
%     sinogram(sinogram < 0) = 0;

    % default wi
    if isempty(wi)
      wi = ones(size(sinogram));
    end
    
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
    % plot per time or per iteration
    if ~isfield(arg,'per_time') 
      arg.per_time = 0;
    end
    
    % figure file
    pat = sprintf('%s/%s%s-ph_%s-%d_nt_%s-nl_%.3f-np_%d-p_%s', ...
      arg.path, arg.prefix, plt, arg.phan, arg.phan_size, arg.noise_type, ...
      arg.noise_level, arg.num_proj, arg.proj_type);
    fig_file = pat; if arg.per_time, fig_file = [fig_file '_t-1']; end
    mat_file = [pat '.mat'];
    
    % initialize
    init;
    
    % figure file
    snr_file = [fig_file '-snr'];
    resid_file = [fig_file '-resid'];
    rec_file = [fig_file '-rec'];
    
    fprintf('Fig %s\n', fig_file);

    % loop and compute
    if arg.recompute
      results = cell(size(algs));
      for a = 1:length(algs)
        alg = algs{a};

        % run
        [rec, tt, ss, it, rs] = run_alg(alg);

        % put results
      %   results{a}.rec = rec;
        results{a}.times = tt;
        results{a}.snrs = ss;
        results{a}.iter = it;
        results{a}.resid = rs;
      end

      % save
      if ~isfield(arg, 'save') || arg.save
        save(mat_file, 'results', 'algs'); 
      end
    else
      rr = load(mat_file);
      results = rr.results;
%       algs = rr.algs;
    end

    % SNR plot
    fig_size = [800, 800];
    if isfield(arg, 'fig_size'), fig_size = arg.fig_size; end
    snr_fig = figure('Name',snr_file, 'Position',[1, 1, fig_size]);
    hold on;
    
    legends = cell(size(results));
    handles = [];
    max_it = 0;
    for a = 1:length(algs)
      res = results{a};
      alg = algs{a};
      
      if arg.per_time
        xs = res.times;
      else
        xs = res.iter;
        max_it = max(max_it, xs(end));
      end
      handles(a) = plot(xs, res.snrs, 'Color',alg.clr, ...
        'LineStyle',alg.lstyle, 'Marker', alg.marker, ...
        'LineWidth',1);
            
      legends{a} = alg.name;
    end
%     title(sprintf('Phantom %s', phan));
    if arg.per_time
      xlabel('Time (sec)');
    else
      xlabel('Iteration');
      xlim([1 max_it]);
    end
    ylabel('SNR (db)');    
    % exclude certain plots from legend
    if isfield(arg, 'exclude_legend')
      legends(arg.exclude_legend) = [];
      handles(arg.exclude_legend) = [];
    end;
    hleg = legendflex(handles, legends, 'anchor',{'se','se'}, ...
      'ref',gca, 'buffer',[-10 2], 'ncol',arg.legend_cols, ...
      'box','on', 'padding',[2, 1, 6]);
%     hleg = columnlegend(arg.legend_cols, legends, 'boxon',...
%       'Location','SouthEast');
%     set(hleg, 'box','off');
%     legHdl = gridLegend(gca, arg.legend_cols, legends, ...
%       'Location','SouthEast'); 
    % title
    if isfield(arg,'title') 
      switch arg.title
        % number of projections
        case 'num_proj', title(sprintf('%d Projections', arg.num_proj));
      end
    end
    
    % save figure
    if ~isfield(arg, 'save') || arg.save
      save_fig(snr_file, snr_fig, 'pdf', '-fonts');
    end
    
    % Residual plot
    if arg.plot_resid
      resid_fig = figure('Name',resid_file, 'Position',[1, 1, 800, 800]);    
      hold on;
      legends = cell(size(results));
      for a = 1:length(algs)
        res = results{a};
        alg = algs{a};

        hold on;
        plot(res.iter, 20*log10(res.resid), 'Color',alg.clr, ...
          'LineStyle',alg.lstyle, 'Marker', alg.marker, ...
          'LineWidth',1);

        legends{a} = alg.name;
      end
      xlabel('Iteration');
      ylabel('Residual (db)');
      xlim([1 arg.iter]);
      columnlegend(arg.legend_cols, legends, 'Location','NorthEast');    

      if ~isfield(arg, 'save') || arg.save
        save_fig(resid_file, resid_fig, 'pdf');
      end
    end    
  end

  % Plot SNR per number of projections
  function per_num_proj
    % figure file
    pat = sprintf('%s/%s%s-ph_%s-%d_nt_%s-nl_%.3f-it_%d-p_%s', ...
      arg.path, arg.prefix, plt, arg.phan, arg.phan_size, arg.noise_type, ...
      arg.noise_level, arg.iter, arg.proj_type);
    fig_file = [pat];
    mat_file = [pat '.mat'];
    
    % figure file
    snr_file = [fig_file '-snr'];
    resid_file = [fig_file '-resid'];

    fprintf('Fig %s\n', fig_file);
    
    results = cell(size(algs));

    if arg.recompute
      % init structure
      for a = 1:length(algs)
        results{a} = struct('snrs',[], 'resid',[]);
      end

      % loop over num_proj
      for n=1:length(arg.num_projs)
        % put in arg
        arg.num_proj = arg.num_projs(n);
        % init
        init;
        
        % loop over algorithms
        for a = 1:length(algs)
          alg = algs{a};          

          % run
          [rec, tt, ss, it, rs] = run_alg(alg);

          % put results
          if arg.max_snr 
            results{a}.snrs(end+1) = max(ss);
          else
            results{a}.snrs(end+1) = ss(end);
          end
%           results{a}.resid(end+1) = rs(end);
        end
      end

      % save
      save(mat_file, 'results', 'algs'); 
    else
      rr = load(mat_file);
      results = rr.results;
      algs = rr.algs;
    end

    % SNR plot
    fig_size = [800, 800];
    if isfield(arg, 'fig_size'), fig_size = arg.fig_size; end
    snr_fig = figure('Name',snr_file, 'Position',[1, 1, fig_size]);
    hold on;

    legends = cell(size(results));
    for a = 1:length(algs)
      res = results{a};
      alg = algs{a};

      plot(arg.num_projs, res.snrs, 'Color',alg.clr, ...
        'LineStyle',alg.lstyle, 'Marker', alg.marker, ...
        'LineWidth',1);

      legends{a} = alg.name;
    end
%     title(sprintf('Phantom %s', phan));
    xlabel('Number of Projections');
    ylabel('SNR (db)');
%     xlim([1 arg.iter]);
%     legend(legends, 'Location','NorthWest');
    hleg = legendflex(legends, 'anchor',{'se','se'}, ...
      'ref',gca, 'buffer',[-2 2], 'ncol',arg.legend_cols, ...
      'box','on', 'padding',[2, 1, 6]);
    set(gca, 'xtick',arg.num_projs);

    % save figure
    save_fig(snr_file, snr_fig, 'pdf');
  end

  % Plot SNR per number of projections
  function per_noise_level
    % figure file
    pat = sprintf('%s/%s%s-ph_%s-%d_nt_%s-it_%d-np_%d-p_%s', ...
      arg.path, arg.prefix, plt, arg.phan, arg.phan_size, arg.noise_type, ...
      arg.iter, arg.num_proj, arg.proj_type);
    fig_file = [pat];
    mat_file = [pat '.mat'];
    
    % figure file
    snr_file = [fig_file '-snr'];
    resid_file = [fig_file '-resid'];

    fprintf('Fig %s\n', fig_file);
    
    results = cell(size(algs));

    if arg.recompute
      % init structure
      for a = 1:length(algs)
        results{a} = struct('snrs',[], 'resid',[]);
      end

      % loop over num_proj
      for n=1:length(arg.noise_levels)
        % put in arg
        arg.noise_level = arg.noise_levels(n);
        % init
        init;
        
        % loop over algorithms
        for a = 1:length(algs)
          alg = algs{a};          

          % run
          [rec, tt, ss, it, rs] = run_alg(alg);

          % put results
          results{a}.snrs(end+1) = ss(end);
          results{a}.resid(end+1) = rs(end);
        end
      end

      % save
      save(mat_file, 'results', 'algs'); 
    else
      rr = load(mat_file);
      results = rr.results;
      algs = rr.algs;
    end

    % SNR plot
    snr_fig = figure('Name',snr_file, 'Position',[1, 1, 800, 800]);
    hold on;

    legends = cell(size(results));
    for a = 1:length(algs)
      res = results{a};
      alg = algs{a};

      plot(arg.noise_levels * 100, res.snrs, 'Color',alg.clr, ...
        'LineStyle',alg.lstyle, 'Marker', alg.marker, ...
        'LineWidth',1);

      legends{a} = alg.name;
    end
%     title(sprintf('Phantom %s', phan));
    xlabel('Noise Level \sigma (%)');
    ylabel('SNR (db)');
%     xlim([1 arg.iter]);
%     legend(legends, 'Location','NorthWest');
    hleg = legendflex(legends, 'anchor',{'ne','ne'}, ...
      'ref',gca, 'buffer',[-2 -2], 'ncol',arg.legend_cols, ...
      'box','on', 'padding',[2, 1, 6]);
    set(gca, 'xtick',arg.noise_levels * 100);

    % save figure
    save_fig(snr_file, snr_fig, 'pdf');
  end

  % Plot Phantom image 
  function phantom
    % phantom and size
    phan = arg.phan;
    phan_size = arg.phan_size;
    phan_file = sprintf('phantoms/%s-%d.mat', phan, phan_size);

    % figure file
    pat = sprintf('%s/%s%s-ph_%s-%d', ...
      arg.path, arg.prefix, plt, arg.phan, arg.phan_size);
    fig_file = [pat];

    % load phantom
    switch phan
    case {'sl', 'mod-sl', 'ncat'}
      P = load(phan_file);
      P = P.P;

    case 'mouse'
      pp = load(phan_file);
      P = pp.P;
    end
    
    % plot
    figh = figure('Name',fig_file, 'Position',[1, 1, 800, 800]);
    imshow(P, arg.range); colorbar;
    
    if ~isfield(arg, 'save') || arg.save
      save_fig(fig_file, figh, 'pdf');
    end
    
  end

 % Plot Reconstruction
  function recon
   
    % figure file
    pat = sprintf('%s/%s%s-ph_%s-%d_nt_%s-nl_%.3f-np_%d-p_%s', ...
      arg.path, arg.prefix, plt, arg.phan, arg.phan_size, arg.noise_type, ...
      arg.noise_level, arg.num_proj, arg.proj_type);
    fig_file = pat;
    mat_file = [pat '.mat'];
    
    % initialize
    init;
    
    % figure file
    rec_file = [fig_file '-rec'];
    
    fprintf('Fig %s\n', fig_file);

    % loop and compute
    if arg.recompute
      results = cell(size(algs));
      for a = 1:length(algs)
        alg = algs{a};

        % run
        [rec, tt, ss, it, rs] = run_alg(alg);

        % put results
      %   results{a}.rec = rec;
        results{a}.rec = rec;
      end

      % save
      if ~isfield(arg, 'save') || arg.save
        save(mat_file, 'results', 'algs'); 
      end
    else
      rr = load(mat_file);
      results = rr.results;
%       algs = rr.algs;
    end

    % Reconstruction plot
    fig_size = [2400, 800];
    if isfield(arg, 'fig_size'), fig_size = arg.fig_size; end
    rec_fig = figure('Name',rec_file, 'Position',[1, 1, fig_size]);
    hold on;
    
    legends = cell(size(results));
    handles = [];
    for a = 1:length(algs)
      res = results{a};
      alg = algs{a};

      % subplot
      subplot(arg.subplot_r, arg.subplot_c, a);
      imshow(res.rec, arg.range); colorbar;
      title(alg.name);
    end
%     % title
%     if isfield(arg,'title') 
%       switch arg.title
%         % number of projections
%         case 'num_proj', title(sprintf('%d Projections', arg.num_proj));
%       end
%     end
    
    % save figure
    if ~isfield(arg, 'save') || arg.save
      save_fig(rec_file, rec_fig, 'pdf', '-fonts');
    end    
  end  

  % run the algorithm
  function [rec, tt, ss, it, rs] = run_alg(alg)
      
    fprintf('Alg %s\n', alg.name);
    
    switch alg.type
    % ASTRA algorithms
    case 'astra'
      % input params
      in_params = struct('vol_geom',vol_geom, 'proj_geom',proj_geom, ...
        'gt_vol',P, 'sino',sinogram, 'proj_id',proj_id, ...
        'wi',wi, 'fbp',fbp, 'prox_in',zeros(size(P)));
      % put fbp in prox input?
      if isfield(arg,'prox_fbp') && arg.prox_fbp
        in_params.prox_in = fbp;
      end
      alg.alg_params.iter = arg.iter;
      % residual
      alg.alg_params.option.ComputeIterationResidual = arg.plot_resid;          

%       if ~alg.alg_params.wls, in_params.wi = ones(size(wi)); end
      
      [rec, tt, ss, it, rs] = ma_alg_astra(alg.alg, in_params, ...
        alg.alg_params);
      
    % IRT toolbox
    case 'irt'
      in_params = struct('vol_geom',vol_geom, 'proj_geom',proj_geom, ...
        'gt_vol',P, 'sino',sinogram, 'A',proj_mat, ...
        'wi',wi, 'fbp',fbp);
      
      if ~alg.alg_params.wls, in_params.wi = ones(size(wi)); end
      
      % set iterations according to if they have inner iterations or not
      switch(alg.alg)
      case {'admm', 'mfista'}
        alg.alg_params.iter = arg.iter / alg.alg_params.nCG;
      otherwise
        alg.alg_params.iter = arg.iter ;
      end


      [rec, tt, ss, it] = ma_alg_irt(alg.alg, in_params, ...
        alg.alg_params);
      rs = 0 * tt;
      
    % Proximal algorithms
    case 'psart'
      % input params
      in_params = struct('vol_geom',vol_geom, 'proj_geom',proj_geom, ...
        'gt_vol',P, 'sino',sinogram, 'proj_id',proj_id, ...
        'wi',wi, 'fbp',fbp);

%       if ~alg.alg_params.wls, in_params.wi = ones(size(wi)); end
      
      alg.alg_params.iter = arg.iter / alg.alg_params.prox_params.iter;

      [rec, tt, ss, it] = ma_alg_psart(alg.alg, in_params, alg.alg_params);
      rs = 0 * tt;
      
    end
  end

end

