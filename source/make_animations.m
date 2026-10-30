% Animations: 2-D k-Wave runs with PROTEUS transducer geometry, delays and
% pulses (PROTEUS itself runs in 3-D; 2-D keeps these animations small).
here = fileparts(mfilename('fullpath'));
repo = fullfile(here, '..', 'PROTEUS');
out  = fullfile(here, 'assets');
addpath(repo, fullfile(repo,'GUIfunctions'), fullfile(repo,'acoustic-module'), ...
    fullfile(repo,'external','k-wave','k-Wave'), fullfile(repo,'tests','helpers'));
set(groot,'defaultFigureVisible','off','defaultAxesFontSize',11);

P41 = reference_transducer(reset_transducer());
C16.Type = 'C1-6D'; C16 = get_transducer_properties(C16);

cfg(1) = struct('T', P41, 'focus', 30e-3, 'angle', 0, 'depth', 45e-3, 'width', 36e-3, ...
    'ppw', 5, 'tend_factor', 2.0, 'name', 'kwave_P41_focused', ...
    'title', 'P4-1, focused at 30 mm: transmit, scattering, echoes');
cfg(2) = struct('T', C16, 'focus', Inf, 'angle', 0, 'depth', 60e-3, 'width', 80e-3, ...
    'ppw', 3, 'tend_factor', 1.2, 'name', 'kwave_C16D_diverging', ...
    'title', 'C1-6D, diverging wave from the centre of curvature');

for k = 1:numel(cfg)
    run_case(cfg(k), out);
end
disp('animations done')

function run_case(cfg, out)
T = cfg.T; c0 = 1540;
Tx = reset_transmit(T);                     % pulse at Tx.SamplingRate
Tx.LateralFocus = cfg.focus; Tx.Angle = cfg.angle;
Tx = compute_delays(Tx, T, struct('SpeedOfSound', c0));
El = get_element_geometry(T);

lambda = c0/T.CenterFrequency;
dx = lambda/cfg.ppw;
x0 = -2e-3;                                  % grid starts behind the array
Nx = round((cfg.depth - x0)/dx); Ny = round(cfg.width/dx);
Nx = Nx + mod(Nx,2); Ny = Ny + mod(Ny,2);
kgrid = kWaveGrid(Nx, dx, Ny, dx);
xg = x0 + (0:Nx-1)*dx; yg = ((0:Ny-1) - Ny/2)*dx;

% Medium: weak speckle, coupling region behind curved arrays, 4 point
% scatterers (stand-ins for bubbles: in PROTEUS, bubbles are mass sources
% driven by the bubble ODE, not density contrasts).
rng(3);
[X, Y] = ndgrid(xg, yg);
bg = 1 + 0.005*randn(Nx, Ny);
if isfinite(El.Radius)
    bg((X - El.Center).^2 + Y.^2 <= El.Radius^2) = 1;
end
medium.sound_speed = c0*bg; medium.density = 1000*bg;
scat = [0.45 -0.15; 0.6 0.05; 0.75 0.2; 0.9 -0.05].*[cfg.depth cfg.width/2];
for s = 1:size(scat,1)
    m = hypot(X - scat(s,1), Y - scat(s,2)) < 0.25e-3;
    medium.density(m) = 1500; medium.sound_speed(m) = 2000;
end
medium.alpha_coeff = 0.5; medium.alpha_power = 1.5;
kgrid.makeTime(medium.sound_speed, 0.3, cfg.tend_factor*cfg.depth/c0);

% Source: element surface points (z = 0 row of the integration points),
% snapped to the grid, each driven by the pulse delayed by its element delay.
G.dy = dx; G.dz = dx; T.IntegrationDensity = 2;
T = get_transducer_integration_points(T, G);
P = T.integration_points;                    % N_el-by-N_int-by-3
zmid = abs(P(:,:,3)) < dx/2;
ix = []; iy = []; el = [];
for e = 1:T.NumberOfElements
    pts = squeeze(P(e, zmid(e,:), 1:2));
    if isvector(pts); pts = pts(:)'; end
    ix = [ix; round((pts(:,1) - x0)/dx) + 1]; %#ok<AGROW>
    iy = [iy; round(pts(:,2)/dx + Ny/2) + 1]; %#ok<AGROW>
    el = [el; e*ones(size(pts,1),1)];        %#ok<AGROW>
end
[lin, ia] = unique(sub2ind([Nx Ny], ix, iy)); el = el(ia);
source.p_mask = false(Nx, Ny); source.p_mask(lin) = true;
[~, order] = sort(lin); el = el(order);      % k-Wave orders sources by index
pulse = Tx.PressureSignal; tp = (0:numel(pulse)-1)/Tx.SamplingRate;
t = kgrid.t_array;
sig = zeros(numel(el), numel(t));
for i = 1:numel(el)
    sig(i,:) = interp1(tp + Tx.Delays(el(i)), pulse, t, 'linear', 0);
end
source.p = sig/1e3;                          % kPa scale
source.p_mode = 'additive';

% Sensor: every 2nd grid point, all time steps
step = 2;
sensor.mask = false(Nx, Ny); sensor.mask(1:step:end, 1:step:end) = true;
sensor.record = {'p'};
data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false, ...
    'PMLSize', 20, 'PlotSim', false, 'DataCast', 'single', 'Smooth', false);
nx = numel(1:step:Nx); ny = numel(1:step:Ny);
Pm = reshape(data.p, nx, ny, []);
xs = xg(1:step:end); ys = yg(1:step:end);

% Frames: compressed amplitude so that weak echoes stay visible
nf = 90; idx = round(linspace(1, size(Pm,3), nf));
front = xs > max(El.x) + 3e-3;                % exclude the source region
a = 0.08*max(abs(Pm(front,:,:)), [], 'all');  % one scale, saturating near focus
cmap = interp1([0 0.5 1], [0.13 0.40 0.69; 1 1 1; 0.80 0.10 0.10], linspace(0,1,256));
f = figure('Position', [0 0 520 round(520*(cfg.depth-x0)/cfg.width)+80]);
theme(f, 'light');
ax = axes(f, 'Position', [0.12 0.08 0.84 0.82]);
vw = VideoWriter(fullfile(out, [cfg.name '.mp4']), 'MPEG-4');
vw.FrameRate = 12; vw.Quality = 90; open(vw);
fp = f.Position; sz = 16*ceil([fp(4) fp(3)]/16);   % frame size (pixels)
Eg = [El.x El.y];
for k = 1:nf
    img = Pm(:,:,idx(k));
    img = sign(img).*min(abs(img)/a, 1).^0.5;   % same scale for all frames
    imagesc(ax, ys*1e3, xs*1e3, img); clim(ax, [-1 1]); colormap(ax, cmap); axis(ax, 'image')
    hold(ax, 'on')
    plot(ax, Eg(:,2)*1e3, Eg(:,1)*1e3, '-', 'Color', [0 0 0], 'LineWidth', 3)
    plot(ax, scat(:,2)*1e3, scat(:,1)*1e3, 'o', 'Color', [0.2 0.2 0.2], 'MarkerSize', 7)
    hold(ax, 'off')
    xlabel(ax, 'y  lateral (mm)'); ylabel(ax, 'x  depth (mm)')
    title(ax, sprintf('%s\nt = %.1f \\mus', cfg.title, t(idx(k))*1e6), 'FontSize', 10)
    drawnow
    im = frame2im(getframe(f));
    im = im(1:2:end, 1:2:end, :);            % screen captures are 2x
    fr = 255*ones([sz 3], 'uint8');          % fixed size, white padding
    h = min(sz(1), size(im,1)); w = min(sz(2), size(im,2));
    fr(1:h, 1:w, :) = im(1:h, 1:w, :);
    writeVideo(vw, fr);
end
close(vw)
close(f)
fprintf('%s: grid %dx%d, %d steps\n', cfg.name, Nx, Ny, kgrid.Nt);
end
