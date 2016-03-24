function [rec, times, snrs, iters] = ma_alg_psart(alg, in_params, alg_params)
% MA_ALG_PSART runs ProxiSART
% min (1/sigma) ||Ax-b||^2 + R(x)
% INPUT
% - alg: the name of the algorithm 'admm', 'cp'
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


% x = [1 2 3; 5 8 12; 15, 21, 30]
% dx = fdiff(x)
% xx = ndiv(dx)
[times snrs iters] = deal(zeros(alg_params.iter, 1));

switch alg
  case 'admm'
    rec = admm();
    
  case 'cp'
    rec = cp();
    
  otherwise
    error(['Unknown algorithm: ' alg]);
end

times = cumsum(times);
iters = cumsum(iters);

  function x = admm()
    rho = alg_params.rho;
    mu = alg_params.mu;
    if isempty(mu), mu = 1 / (rho * 8); end    
    sigma = alg_params.sigma;
    mu/sigma
    
    % initialize from FBP
    if alg_params.init_fbp 
      % read initial value from recid
      x = in_params.fbp;
    else
      x = zeros(size(in_params.fbp));
    end
    % inital values
    z = fdiff(x);
    u = z;

%     [ny, nx] = size(in_params.gt_vol);
%     [nys, nxs] = size(in_params.sino);
%     cids = reshape(reshape(1:nx*ny, nx,ny)',[],1);
%     rids = reshape(reshape(1:nxs*nys, nxs,nys)',[],1);
%     A = Gmatrix(in_params.A(rids, cids));
%     % solve with CG
%     AA = @(x)(2*mu/sigma * (A' * A) * x + x);
%     Ab = sqrt(2*mu/sigma) * A' * in_params.sino(:) ;
    
    % loop
    for it=1:alg_params.iter

      % x-step
    %   new_x = l2_data_prox(cfg, x - lambda * ndiv(z), lambda, 20);
%       x = l2_data_prox(cfg, x - rho*mu * ndiv(fdiff(x)-z+u), ...
%         mu/sigma, 2);
      % 1/sigma with data term
      if alg_params.sigma_with_data
        [x, itimes, isnrs, iiters] = l2_data_prox(...
          x - rho*mu * ndiv(fdiff(x)-z+u), mu / sigma, ... 1e-3
          alg_params.psart_params, in_params, alg_params.prox);
      else
        % sigma with regularizer
        [x, itimes, isnrs, iiters] = l2_data_prox(...
          x - rho*mu * ndiv(fdiff(x)-z+u), mu, ...
          alg_params.psart_params, in_params, alg_params.prox);
      end
      x = max(0, x);
      times(it) = itimes(end);
      snrs(it) = isnrs(end);
      iters(it) = iiters(end);
    
%       % ids to convert from row-major to column-major in both input and output
%       [x] = pcg(AA, Ab + reshape(x - rho*mu * ndiv(fdiff(x)-z+u),[],1), ...
%         1e-10, alg_params.psart_params.iter, [], [], x(:));
%       x = reshape(x, ny, nx);
%       x = max(0, min(1, x));
% %       max(x(:))
%       times(it) = 0;
%       snrs(it) = ma_snr(in_params.gt_vol, in_params.gt_vol - x);
%       iters(it) = 2;      

      % z-step
      tt = tic;
      if alg_params.sigma_with_data
        z = atv_prox(fdiff(x) + u, 1/rho, 1);
      else
        z = atv_prox(fdiff(x) + u, 1/rho, sigma);
      end
      
    %   z = z - (1/rho) * atv_conj_prox((fdiff(x) + u) * rho, rho, 1);

      % u-step
      u = u + fdiff(x) - z;
      
      times(it) = times(it) + toc(tt);
      
      % objective, can be increasing or decreasing...
      [obj1, obj2] = objective(x, z, sigma);
      
      % primal residual, should be decreasing per iteration...
      r = rho * norm(reshape(fdiff(x),[],1) - z(:));

      fprintf('it=%d\tobj1=%.2f\tobj2=%.2f\tsnr=%.2f\tres=%f\n', it, obj1, obj2, snrs(it), r);
    end
%   figure(3), imshow(x, []);
  end % admm

  function x = cp()
    % sigma
    sigma = alg_params.sigma;
    lambda = alg_params.lambda;  %0.01; %.001
    mu = alg_params.mu; %1 / (lambda * 8) % * sigma^2)
    if isempty(mu), mu = 1 / (lambda * 8); end
    theta = alg_params.theta; % 1;

    % initialize from FBP
    if alg_params.init_fbp 
      % read initial value from recid
      x = in_params.fbp;
    else
      x = zeros(size(in_params.fbp));
    end
    xbar = x;
    z = fdiff(x);


    % gt dual variable z (to measure dual gap)
    % gt_z = sign(fdiff(P));
    % gt_z(gt_z == 0) = 1;

    % loop
    for it=1:alg_params.iter
      % z-step
      tt = tic;
      z = atv_conj_prox(z + mu * fdiff(xbar), mu, 1);
      times(it) = times(it) + toc(tt);
    %   z = atv_conj_prox(z + mu * fdiff(xbar), mu, sigma);
    %   z = itv_conj_prox(z + mu * fdiff(xbar), mu, sigma);

      % x-step
      % new_x = l2_data_prox(cfg, x - lambda * ndiv(z), lambda, 20);
      % new_x = l2_data_prox(cfg, x - lambda * ndiv(z), lambda/sigma, 20);      
      [new_x, itimes, isnrs, iiters] = l2_data_prox(...
        x - lambda * ndiv(z), lambda / sigma, ...
        alg_params.psart_params, in_params);
      times(it) = itimes(end);
      snrs(it) = isnrs(end);
      iters(it) = iiters(end);

      % xbar
      tt = tic;
      xbar = new_x + theta * (new_x - x);
      x = new_x;
      times(it) = times(it) + toc(tt);

%       % objective
%       obj = objective(x, sigma, sinogram, cfg);
% 
%       % primal-dual gap
%       pdg = norm(P - x, 'fro')^2 / (2*lambda) + sum(reshape(gt_z - z, [],1).^2) / (2*mu);
% 
%       fprintf('it=%d\tobj=%.2f\tsnr=%.2f\tgap=%.2f\n', it, obj, snr(P, x-P), pdg);
    end;
  end

  % Computes the objective function
  function [val1, val2] = objective(x, z, sigma)
    % geom
    proj_id = in_params.proj_id;
    proj_geom = astra_mex_projector('projection_geometry', proj_id);
    vol_geom = astra_mex_projector('volume_geometry', proj_id);

    % create astra volume 
    volume_id = astra_mex_data2d('create','-vol', vol_geom, x);

    % crate sino
    sino_id = astra_mex_data2d('create','-sino', proj_geom, 0);

    % project
    cfg = astra_struct('FP');
    cfg.ProjectorId = proj_id;
    cfg.ProjectionDataId = sino_id;
    cfg.VolumeDataId = volume_id;
    alg_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('iterate', alg_id);
    astra_mex_algorithm('delete', alg_id);

    % compute objective
    sino = astra_mex_data2d('get',sino_id);
    if alg_params.sigma_with_data
      val1 = norm(sino - in_params.sino, 'fro')^2 / sigma;
      val2 = sum(abs(z(:))); %atv(x);
    else
      val1 = norm(sino - in_params.sino, 'fro')^2;
      val2 = sigma * sum(abs(z(:))); % atv(x);
    end

    astra_mex_data2d('delete', volume_id);
    astra_mex_data2d('delete', sino_id);
  end
% vol_geom = astra_create_vol_geom(256, 256);
% proj_geom = astra_create_proj_geom('parallel', 1.0, 384, ...
%   linspace2(0,pi,30));
% % vol_geom = astra_create_vol_geom(4,4);
% % proj_geom = astra_create_proj_geom('parallel', 1.0, 10, linspace2(0,pi,4));
% 
% % For CPU-based algorithms, a "projector" object specifies the projection
% % model used. In this case, we use the "strip" model.
% proj_id = astra_create_projector('linear', proj_geom, vol_geom);
% 
% % Create a sinogram from a phantom
% P = phantom(256);
% % P = phantom(4);
% [sinogram_id, sinogram] = astra_create_sino(P, proj_id);
% figure(1); imshow(P, []); axis image; axis off;
% figure(2); imshow(sinogram, []); axis image; axis off;
% % imshow(sinogram, []);
% 
% % add sinogram noise
% rng(123);
% sinogram = sinogram + randn(size(sinogram)) * 0.5;
% 
% astra_mex_data2d('delete', sinogram_id);
% 
% % We now re-create the sinogram data object as we would do when loading
% % an external sinogram
% sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogram);
% sinogram_id2 = astra_mex_data2d('create', '-sino', proj_geom, sinogram);
% 
% % Create a data object for the reconstruction
% rec_id = astra_mex_data2d('create', '-vol', vol_geom);
% 
% % create a data object for the proximal input
% prox_in = zeros(size(P));
% prox_in_id = astra_mex_data2d('create', '-vol', vol_geom, prox_in);
% 
% % Set up the parameters for a reconstruction algorithm using the CPU
% % The main difference with the configuration of a GPU algorithm is the
% % extra ProjectorId setting.
% cfg = astra_struct('PSART');
% cfg.ReconstructionDataId = rec_id;
% cfg.ProjectionDataId = sinogram_id;
% cfg.ProjectorId = proj_id;
% cfg.ProxInputDataId = prox_in_id;
% cfg.option.UseMinConstraint = 1;
% cfg.option.MinConstraintValue = 0;
% cfg.option.UseMaxConstraint = 0;
% cfg.option.MaxConstraintValue = 1;
% cfg.option.ClearRayLength = 1;
% cfg.option.Alpha = 2;
% cfg.option.Lambda = 1e8;
% 
% % initialize with FBP
% cfgi = cfg;
% cfgi.type = 'FBP';
% alg_id = astra_mex_algorithm('create', cfgi);
% astra_mex_algorithm('iterate', alg_id, 5);
% astra_mex_algorithm('delete', alg_id);
% % now initial value is stored in rec_id
% 
% init = true;
% 
% % CP
% % CP(cfg, P, sinogram);
% 
% % ADMM
% ADMM(cfg, P, sinogram, init);
% 
% % cleanup
% astra_mex_projector('delete', proj_id);
% astra_mex_data2d('delete', rec_id);
% astra_mex_data2d('delete', sinogram_id);
% astra_mex_data2d('delete', prox_in_id);

end

% function admm(cfg, P, sinogram, init)
% % sigma
% sigma = 100; %100
% 
% % ADMM algorithm parameters
% rho = 10; %10
% 
% mu = 1 / (8 * rho)
% % lambda * mu * 8
% theta = 1;
% 
% % init
% if init 
%   % read initial value from recid
%   x = astra_mex_data2d('get', cfg.ReconstructionDataId);
% else
%   x = zeros(size(P));
% end
% z = fdiff(x);
% u = z;
% 
% % loop
% for it=1:20
%   
%   % x-step
% %   new_x = l2_data_prox(cfg, x - lambda * ndiv(z), lambda, 20);
%   x = l2_data_prox(cfg, x - rho*mu * ndiv(fdiff(x)-z+u), ...
%     mu/sigma, 2);
% 
%   % z-step
%   z = atv_prox(fdiff(x) + u, 1/rho, 1);
% %   z = z - (1/rho) * atv_conj_prox((fdiff(x) + u) * rho, rho, 1);
% %   z = atv_conj_prox(z + mu * fdiff(xbar), mu, sigma);
% %   z = itv_conj_prox(z + mu * fdiff(xbar), mu, sigma);
%   
%   % u-step
%   u = u + fdiff(x) - z;
%   
% %   x = new_x;
%   
%   % objective
%   obj = objective(x, sigma, sinogram, cfg);
%     
%   fprintf('it=%d\tobj=%.2f\tsnr=%.2f\t\n', it, obj, snr(P, x-P));
% end;
% 
% figure(3), imshow(x, []);
% end

% function CP(cfg, P, sinogram)
% % sigma
% sigma = 20;
% 
% % CP algorithm parameters
% lambda = 0.01; %.001
% % mu = 10; %100
% mu = 1 / (lambda * 8) % * sigma^2)
% % mu = 0.1
% % lambda * mu * 8
% theta = 1;
% 
% % gamma = 0.7 / sigma/2;
% 
% % init
% x = zeros(size(P));
% xbar = x;
% z = fdiff(x);
% 
% % Power Iteration to estimate spectral norm of first diff matrix
% % xx = randn(size(P));
% % for i=1:1e3
% %   xx = ndiv(fdiff(xx));
% %   xx = xx / sqrt(sum(reshape(xx,[],1).^2));
% % end
% % val = sigma^2 * sum(reshape(fdiff(xx).^2,[],1)) / sum(reshape(xx.^2,[],1))
% % return;
% 
% % gt dual variable z 
% gt_z = sign(fdiff(P));
% % gt_z(gt_z == 0) = 1;
% 
% % loop
% for it=1:20
%   % z-step
%   z = atv_conj_prox(z + mu * fdiff(xbar), mu, 1);
% %   z = atv_conj_prox(z + mu * fdiff(xbar), mu, sigma);
% %   z = itv_conj_prox(z + mu * fdiff(xbar), mu, sigma);
%   
%   % x-step
% %   new_x = l2_data_prox(cfg, x - lambda * ndiv(z), lambda, 20);
%   new_x = l2_data_prox(cfg, x - lambda * ndiv(z), lambda/sigma, 20);
%   
% %   theta = 1 / sqrt(1 + 2 * gamma * lambda);
% %   lambda = theta * lambda;
% %   mu = mu / theta;
%   
%   % xbar
%   xbar = new_x + theta * (new_x - x);
%   x = new_x;
%   
%   % objective
%   obj = objective(x, sigma, sinogram, cfg);
%   
%   % primal-dual gap
%   pdg = norm(P - x, 'fro')^2 / (2*lambda) + sum(reshape(gt_z - z, [],1).^2) / (2*mu);
%   
%   fprintf('it=%d\tobj=%.2f\tsnr=%.2f\tgap=%.2f\n', it, obj, snr(P, x-P), pdg);
% end;
% 
% figure(3), imshow(x, []);
% end

% % Computes the objective function
% function val = objective(x, sigma, sinogram, psart_cfg)
% % geom
% proj_id = psart_cfg.ProjectorId;
% proj_geom = astra_mex_projector('projection_geometry', proj_id);
% vol_geom = astra_mex_projector('volume_geometry', proj_id);
% 
% % create astra volume 
% volume_id = astra_mex_data2d('create','-vol', vol_geom, x);
% 
% % crate sino
% sino_id = astra_mex_data2d('create','-sino', proj_geom, 0);
% 
% % project
% cfg = astra_struct('FP');
% cfg.ProjectorId = proj_id;
% cfg.ProjectionDataId = sino_id;
% cfg.VolumeDataId = volume_id;
% alg_id = astra_mex_algorithm('create', cfg);
% astra_mex_algorithm('iterate', alg_id);
% astra_mex_algorithm('delete', alg_id);
% 
% % compute objective
% sino = astra_mex_data2d('get',sino_id);
% % val = norm(sino - sinogram, 'fro')^2 + sigma * atv(x);
% val = norm(sino - sinogram, 'fro')^2 / sigma + atv(x);
% 
% astra_mex_data2d('delete', volume_id);
% astra_mex_data2d('delete', sino_id);
% end

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

% L2 data proximal operator using PSART
% function [xp] = l2_data_prox(cfg, prox_in, lambda, it)
function [xp, time, snr, iter] = l2_data_prox(prox_in, lambda, ...
  psart_params, in_params, which_prox)

% Using SartProx
switch which_prox
  case 'astra'
    % Set prox_in and Lambda
    in_params.prox_in = prox_in;
    psart_params.option.Lambda = lambda;
    % Run
    [xp, time, snr, iter] = ma_alg_astra('SART-PROX', ...
            in_params, psart_params);
  case 'matlab'
    % Set prox_in and Lambda
    in_params.prox_in = prox_in;
    psart_params.lambda = lambda;
    % Run  
    [xp, time, snr, iter] = ma_alg_sart(in_params, psart_params);  
end
      
% % store input
% astra_mex_data2d('store', cfg.ProxInputDataId, prox_in);
% % set lambda
% cfg.option.Lambda = lambda;
% 
% % Create the algorithm object from the configuration structure
% alg_id = astra_mex_algorithm('create', cfg);
% 
% % run
% tic;
% astra_mex_algorithm('iterate', alg_id, it);
% toc;
% 
% % Get the result
% xp = astra_mex_data2d('get', cfg.ReconstructionDataId);
% % figure(3); imshow(rec, []); axis image; axis off; %, []);
% % fprintf('%s SNR=%f\n', cfg.type, snr(P, (P-rec)));
% 
% % clear algorithm
% astra_mex_algorithm('delete', alg_id);
end

% Anisotropic TV convex conjugate proximal operator
% z is a "gradient" object i.e. 3-dim with horizontal diff in first channel
% and vertical diff in second channel
% This is equivalent to the proximal operator of the convex conjugate of
%   sigma * ||z||_1 i.e. prox_{mu conj(sigma * ||.||_1)}
% and mu does not affect the result because it's simply a projection
function zp = atv_conj_prox(z, mu, sigma)
sz = size(z);
assert(length(sz) == 3);

% project on L_inf ball with radius sigma
zp = max(-sigma, min(sigma, z));
end

% This is equivalent to the proximal operator of sigma * ||z||_1 i.e.
%   prox{mu * \sigma ||z||_1} and we get it via Moreau decomposition from
%   its convex conjugate
function zp = atv_prox(z, mu, sigma)
sz = size(z);
assert(length(sz) == 3);

% Use Moreau decomposition plus the relation of convex conjugate
% zp = z - atv_conj_prox(z, 1 / (mu * sigma), mu * sigma);

% Use normal formula
zp = sign(z) .* max(0, abs(z) - mu * sigma);
end

% Anisotropic TV objective
function g = atv(x)
z = fdiff(x);
g = sum(abs(z(:)));
end

% Isotropic TV convex conjugate proximal operator
% z is a "gradient" object i.e. 3-dim with horizontal diff in first channel
% and vertical diff in second channel
function zp = itv_conj_prox(z, mu, sigma)
sz = size(z);
assert(length(sz) == 3);

% compute magnitude of gradient at each voxel
mag = sqrt(sum(z.^2, 3));
assert(all(size(mag) == sz(1:2)));

% project on L_2 ball with radius sigma
zp = bsxfun(@rdivide, sigma*z, max(sigma, mag));
end

% Isotropic TV objective
function g = itv(x)
z = fdiff(x);
mag = sqrt(sum(z.^2, 3));
g = sum(mag(:));
end

%%
% P = zeros(256, 256);
% 
% % Wavelet decomposition operator
% R = Gwave2('mask', true(size(P)), 'dwtmode', 'per', ...
%   'redundancy', 'undecimated', 'wname', 'haar', ...
%   'nlevel', 2, 'includeApprox', false, ...
%   'scale_filters', 1/sqrt(2));
% 
% % Power Iteration to estimate spectral norm of first diff matrix
%   xx = randn(size(P));
%   for i=1:1e2
% %     xx = ndiv(fdiff(xx));
%     xx = R' * R * xx;
%     xx = xx / sqrt(sum(reshape(xx,[],1).^2));
%   end
% %   val = sum(reshape(fdiff(xx).^2,[],1)) / sum(reshape(xx.^2,[],1))
% val = sum(reshape(R* xx.^2,[],1)) / sum(reshape(xx.^2,[],1))
  
