function [] = CP_PSART()
% minimizes ||Ax-b||^2 / sigma + TV(x)

x = [1 2 3; 5 8 12; 15, 21, 30]
dx = fdiff(x)
xx = ndiv(dx)

vol_geom = astra_create_vol_geom(256, 256);
proj_geom = astra_create_proj_geom('parallel', 1.0, 384, ...
  linspace2(0,pi,30));
% vol_geom = astra_create_vol_geom(4,4);
% proj_geom = astra_create_proj_geom('parallel', 1.0, 10, linspace2(0,pi,4));

% For CPU-based algorithms, a "projector" object specifies the projection
% model used. In this case, we use the "strip" model.
proj_id = astra_create_projector('linear', proj_geom, vol_geom);

% Create a sinogram from a phantom
P = phantom(256);
% P = phantom(4);
[sinogram_id, sinogram] = astra_create_sino(P, proj_id);
figure(1); imshow(P, []); axis image; axis off;
figure(2); imshow(sinogram, []); axis image; axis off;
% imshow(sinogram, []);

% add sinogram noise
rng(123);
sinogram = sinogram + randn(size(sinogram)) * 0.5;

astra_mex_data2d('delete', sinogram_id);

% We now re-create the sinogram data object as we would do when loading
% an external sinogram
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogram);
sinogram_id2 = astra_mex_data2d('create', '-sino', proj_geom, sinogram);

% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom);

% create a data object for the proximal input
prox_in = zeros(size(P));
prox_in_id = astra_mex_data2d('create', '-vol', vol_geom, prox_in);

% Set up the parameters for a reconstruction algorithm using the CPU
% The main difference with the configuration of a GPU algorithm is the
% extra ProjectorId setting.
cfg = astra_struct('PSART');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ProjectorId = proj_id;
cfg.ProxInputDataId = prox_in_id;
cfg.option.UseMinConstraint = 1;
cfg.option.MinConstraintValue = 0;
cfg.option.UseMaxConstraint = 1;
cfg.option.MaxConstraintValue = 1;
cfg.option.ClearRayLength = 1;
cfg.option.Alpha = 2;
cfg.option.Lambda = 1e8;

% CP
% CP(cfg, P, sinogram);

% ADMM
ADMM(cfg, P, sinogram);

% cleanup
astra_mex_projector('delete', proj_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sinogram_id);
astra_mex_data2d('delete', prox_in_id);

end

function ADMM(cfg, P, sinogram)
% sigma
sigma = 1000;

% CP algorithm parameters
rho = 10; %.001

mu = 1 / (8 * rho)
% lambda * mu * 8
theta = 1;

% init
x = zeros(size(P));
z = fdiff(x);
u = z;

% loop
for it=1:20
  
  % x-step
%   new_x = l2_data_prox(cfg, x - lambda * ndiv(z), lambda, 20);
  x = l2_data_prox(cfg, x - rho*mu * ndiv(fdiff(x)-z+u), ...
    mu/sigma, 20);

  % z-step
  z = atv_prox(fdiff(x) + u, 1/rho, 1);
%   z = z - (1/rho) * atv_conj_prox((fdiff(x) + u) * rho, rho, 1);
%   z = atv_conj_prox(z + mu * fdiff(xbar), mu, sigma);
%   z = itv_conj_prox(z + mu * fdiff(xbar), mu, sigma);
  
  % u-step
  u = u + fdiff(x) - z;
  
%   x = new_x;
  
  % objective
  obj = objective(x, sigma, sinogram, cfg);
    
  fprintf('it=%d\tobj=%.2f\tsnr=%.2f\t\n', it, obj, snr(P, x-P));
end;

figure(3), imshow(x, []);
end

function CP(cfg, P, sinogram)
% sigma
sigma = 20;

% CP algorithm parameters
lambda = 0.01; %.001
% mu = 10; %100
mu = 1 / (lambda * 8) % * sigma^2)
% mu = 0.1
% lambda * mu * 8
theta = 1;

% gamma = 0.7 / sigma/2;

% init
x = zeros(size(P));
xbar = x;
z = fdiff(x);

% Power Iteration to estimate spectral norm of first diff matrix
% xx = randn(size(P));
% for i=1:1e3
%   xx = ndiv(fdiff(xx));
%   xx = xx / sqrt(sum(reshape(xx,[],1).^2));
% end
% val = sigma^2 * sum(reshape(fdiff(xx).^2,[],1)) / sum(reshape(xx.^2,[],1))
% return;

% gt dual variable z 
gt_z = sign(fdiff(P));
% gt_z(gt_z == 0) = 1;

% loop
for it=1:20
  % z-step
  z = atv_conj_prox(z + mu * fdiff(xbar), mu, 1);
%   z = atv_conj_prox(z + mu * fdiff(xbar), mu, sigma);
%   z = itv_conj_prox(z + mu * fdiff(xbar), mu, sigma);
  
  % x-step
%   new_x = l2_data_prox(cfg, x - lambda * ndiv(z), lambda, 20);
  new_x = l2_data_prox(cfg, x - lambda * ndiv(z), lambda/sigma, 20);
  
%   theta = 1 / sqrt(1 + 2 * gamma * lambda);
%   lambda = theta * lambda;
%   mu = mu / theta;
  
  % xbar
  xbar = new_x + theta * (new_x - x);
  x = new_x;
  
  % objective
  obj = objective(x, sigma, sinogram, cfg);
  
  % primal-dual gap
  pdg = norm(P - x, 'fro')^2 / (2*lambda) + sum(reshape(gt_z - z, [],1).^2) / (2*mu);
  
  fprintf('it=%d\tobj=%.2f\tsnr=%.2f\tgap=%.2f\n', it, obj, snr(P, x-P), pdg);
end;

figure(3), imshow(x, []);
end

% Computes the objective function
function val = objective(x, sigma, sinogram, psart_cfg)
% geom
proj_id = psart_cfg.ProjectorId;
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
% val = norm(sino - sinogram, 'fro')^2 + sigma * atv(x);
val = norm(sino - sinogram, 'fro')^2 / sigma + atv(x);

astra_mex_data2d('delete', volume_id);
astra_mex_data2d('delete', sino_id);
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

% L2 data proximal operator using PSART
function xp = l2_data_prox(cfg, prox_in, lambda, it)
% store input
astra_mex_data2d('store', cfg.ProxInputDataId, prox_in);
% set lambda
cfg.option.Lambda = lambda;

% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% run
tic;
astra_mex_algorithm('iterate', alg_id, it);
toc;

% Get the result
xp = astra_mex_data2d('get', cfg.ReconstructionDataId);
% figure(3); imshow(rec, []); axis image; axis off; %, []);
% fprintf('%s SNR=%f\n', cfg.type, snr(P, (P-rec)));

% clear algorithm
astra_mex_algorithm('delete', alg_id);
end

% Anisotropic TV convex conjugate proximal operator
% z is a "gradient" object i.e. 3-dim with horizontal diff in first channel
% and vertical diff in second channel
function zp = atv_conj_prox(z, mu, sigma)
sz = size(z);
assert(length(sz) == 3);

% project on L_inf ball with radius sigma
zp = max(-sigma, min(sigma, z));
end

function zp = atv_prox(z, mu, sigma)
sz = size(z);
assert(length(sz) == 3);

% Use Moreau decomposition
zp = z - mu * atv_conj_prox(z / mu, 1/mu, sigma);
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