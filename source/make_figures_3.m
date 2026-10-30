% Figures 3: flow tracking on a synthetic bifurcation (real generate_streamlines),
% grid vs integration points, and the size of bubble-bubble interactions.
here = fileparts(mfilename('fullpath'));
repo = fullfile(here, '..', 'PROTEUS');
out  = fullfile(here, 'assets');
addpath(fullfile(repo,'tests'), fullfile(repo,'tests','helpers'));
set(groot,'defaultFigureVisible','off','defaultAxesFontSize',11,'defaultLineLineWidth',1.4);
C = [0.13 0.40 0.69; 0.85 0.37 0.01; 0.17 0.63 0.17; 0.58 0.40 0.74; 0.5 0.5 0.5; 0.8 0.1 0.3];

tc  = matlab.unittest.TestCase.forInteractiveUse;
env = tc.applyFixture(ProteusTestEnvironment);
addpath(fullfile(repo,'streamline-module'), fullfile(repo,'microbubble-simulator'), ...
    fullfile(repo,'microbubble-simulator','functions'));

%% A. Synthetic steady velocity field on a Cartesian lattice
dx = 25e-6; Umax = 0.05;                     % [m], [m/s]
xv = 0:dx:4e-3; yv = -1.25e-3:dx:1.25e-3; zv = (-1:1)*dx;
[X, Y] = ndgrid(xv, yv);
% Stream function of a symmetric bifurcation: incompressible by construction.
% psi = f(eta), eta = position across the half channel between the septum
% (or axis) and the wall; the profile blends from a tube centred on the axis
% to a channel with no slip on both walls.
sstep = @(t) 3*min(max(t,0),1).^2 - 2*min(max(t,0),1).^3;
cx  = @(x) 0.85e-3*sstep((x - 1.4e-3)/2.6e-3);          % daughter centreline
rx  = @(x) 0.30e-3 - 0.08e-3*sstep((x - 1.4e-3)/1.2e-3); % local radius
bx  = @(x) sstep((x - 1.5e-3)/1.0e-3);                   % profile blend
yUp = cx(X) + rx(X); yLo = max(0, cx(X) - rx(X));
ay  = abs(Y); eta = (ay - yLo)./(yUp - yLo);
inside = eta >= 0 & eta <= 1;
e = min(max(eta,0),1);
psi = (1 - bx(X)).*(3*e - e.^3)/2 + bx(X).*(3*e.^2 - 2*e.^3);
psi = psi.*sign(Y)*(2/3)*Umax*0.30e-3*1.5;               % flux scale
[dPdy, dPdx] = deal(zeros(size(psi)));
dPdx(2:end-1,:) = (psi(3:end,:) - psi(1:end-2,:))/(2*dx);
dPdy(:,2:end-1) = (psi(:,3:end) - psi(:,1:end-2))/(2*dx);
U = dPdy; V = -dPdx;
U(~inside) = 0; V(~inside) = 0;
% Cells on the wall, axis or edge of the array keep a small forward velocity:
inside = inside & hypot(U,V) > 0;

% Write it in the format of vtu2matlab (only what load_vessel_data reads):
[iX, iY] = find(inside); nc = numel(iX); nz = numel(zv);
pts = zeros(nc*nz, 3); vel = zeros(nc*nz, 3);
for k = 1:nz
    r = (1:nc) + (k-1)*nc;
    pts(r,:) = [xv(iX)' - dx/2, yv(iY)' - dx/2, zv(k)*ones(nc,1) - dx/2];
    sel = sub2ind(size(X), iX, iY);
    vel(r,:) = [U(sel), V(sel), zeros(nc,1)];
end
keep = vecnorm(vel,2,2) > 0; pts = pts(keep,:); vel = vel(keep,:);
vtuStruct.points = pts; vtuStruct.cells = (1:size(pts,1))';
vtuStruct.cellData.velocity_phy = vel;

P = env.PATHS; G.Folder = 'deck_bifurcation';
geo = fullfile(P.GeometriesPath, G.Folder); mkdir(geo);
save(fullfile(geo,'vtu.mat'),'vtuStruct');
vtuProperties = struct('lengthUnit',1,'velocityUnit',1,'velocityField','velocity_phy');
options = odeset('MaxStep', dx/Umax);
save(fullfile(geo,'GeometryProperties.mat'),'vtuProperties','options');
cc = pts + dx/2;                               % cell centres
in = cc(:,1) < 2*dx;
inlet = struct('points', cc(in,:), 'cellsize', [dx dx dx]);
save(fullfile(geo,'inlet.mat'),'inlet');

%% B. The real tracker
rng(4);
Mb = reset_microbubble(); Mb.Number = 10;
Acq = reset_acquisition(); Acq.FrameRate = 100; Acq.NumberOfFrames = 120; Acq.NumberOfPulses = 1;
evalc('generate_streamlines(G, Mb, Acq, P, ''deck_flow'', false)');
pos = nan(Acq.NumberOfFrames, Mb.Number, 3); sn = nan(Acq.NumberOfFrames, Mb.Number); rd = sn;
nd = num2str(length(num2str(Acq.NumberOfFrames)));
for m = 1:Acq.NumberOfFrames
    F = load(fullfile(P.GroundTruthPath,'deck_flow',['Frame_' num2str(m,['%0' nd 'i']) '.mat']));
    pos(m,:,:) = F.Frame.Pulse1.Points; sn(m,:) = F.Frame.Pulse1.StreamNumber; rd(m,:) = F.Frame.Pulse1.Radius;
end

% Figure: the field and the trajectories
speed = hypot(U,V); speed(~inside) = NaN;
cols = lines(Mb.Number);
f = figure('Position',[0 0 1180 420]); theme(f,'light');
tl = tiledlayout(f,1,2,'TileSpacing','compact');
ax = nexttile(tl); h = imagesc(ax, xv*1e3, yv*1e3, speed'*100); set(h,'AlphaData',~isnan(speed')); axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,turbo); hold(ax,'on')
st = 12; quiver(ax, X(1:st:end,1:st:end)*1e3, Y(1:st:end,1:st:end)*1e3, U(1:st:end,1:st:end).*inside(1:st:end,1:st:end), V(1:st:end,1:st:end).*inside(1:st:end,1:st:end), 1.1, 'k', 'LineWidth', 0.6);
cb = colorbar(ax); cb.Label.String = '|v| (cm/s)'; xlabel(ax,'mm'); ylabel(ax,'mm')
title(ax,'Input: one steady velocity field, known at lattice cells')
ax = nexttile(tl); im = imagesc(ax, xv*1e3, yv*1e3, double(inside')); set(im,'AlphaData',0.35*inside'); colormap(ax,[1 1 1; 0.55 0.75 0.55]); axis(ax,'image'); set(ax,'YDir','normal'); hold(ax,'on')
for n = 1:5
    for u = unique(sn(1:40,n))'
        r = sn(:,n) == u & (1:Acq.NumberOfFrames)' <= 40;
        plot(ax, pos(r,n,1)*1e3, pos(r,n,2)*1e3, '-', 'Color', cols(n,:), 'LineWidth', 0.8)
        plot(ax, pos(r,n,1)*1e3, pos(r,n,2)*1e3, '.', 'Color', cols(n,:), 'MarkerSize', 9)
    end
end
xlabel(ax,'mm'); title(ax,'Output of generate\_streamlines: positions at each frame (dots)')
theme(f,'light'); exportgraphics(f, fullfile(out,'flow_field_tracks.png'), 'Resolution', 170); close(f)

% Video of the bubbles moving
f = figure('Position',[0 0 720 460]); theme(f,'light'); ax = axes(f);
vw = VideoWriter(fullfile(out,'flow_tracks.mp4'),'MPEG-4'); vw.FrameRate = 12; vw.Quality = 90; open(vw);
sz = 16*ceil([f.Position(4) f.Position(3)]/16);
for m = 1:Acq.NumberOfFrames
    cla(ax); im = imagesc(ax, xv*1e3, yv*1e3, double(inside')); set(im,'AlphaData',0.35*inside'); colormap(ax,[1 1 1; 0.55 0.75 0.55]); axis(ax,'image'); set(ax,'YDir','normal'); hold(ax,'on')
    scatter(ax, pos(m,:,1)*1e3, pos(m,:,2)*1e3, 12*(rd(m,:)*1e6).^2, cols, 'filled')
    xlabel(ax,'mm'); ylabel(ax,'mm'); title(ax, sprintf('generate\\_streamlines output, frame %d (%.0f ms); marker size follows radius', m, (m-1)/Acq.FrameRate*1e3), 'FontSize', 10)
    drawnow; im2 = frame2im(getframe(f)); im2 = im2(1:2:end,1:2:end,:);
    fr = 255*ones([sz 3],'uint8'); h2 = min(sz(1),size(im2,1)); w2 = min(sz(2),size(im2,2)); fr(1:h2,1:w2,:) = im2(1:h2,1:w2,:);
    writeVideo(vw, fr);
end
close(vw); close(f)
fprintf('streams per bubble: %s\n', mat2str(arrayfun(@(n) numel(unique(sn(:,n))), 1:Mb.Number)));

%% C. Grid nodes, sensor mask and integration points
S = make_tiny_curved_settings(P);
S.SimulationParameters.GridSize = S.Medium.SpeedOfSound/(S.Transmit.CenterFrequency*S.SimulationParameters.PointsPerWavelength);
Grid = define_grid(S.SimulationParameters, S.Geometry);
T = S.Transducer; T.OnGrid = false; T.IntegrationDensity = 1;
T = get_transducer_integration_points(T, Grid);
sensor = define_sensor_transducer(T, Grid);
[~, kz] = min(abs(Grid.z));
[gy, gx] = meshgrid(Grid.y, Grid.x);
f = figure('Position',[0 0 1180 470]); theme(f,'light');
tl = tiledlayout(f,1,2,'TileSpacing','compact');
ax = nexttile(tl); hold(ax,'on')
m = sensor.mask(:,:,kz);
plot(ax, gy(:)*1e3, gx(:)*1e3, '.', 'Color', [0.72 0.72 0.72], 'MarkerSize', 6)
plot(ax, gy(m)*1e3, gx(m)*1e3, 's', 'MarkerFaceColor', [0.65 0.8 0.95], 'MarkerEdgeColor', 'none', 'MarkerSize', 6)
Pp = T.integration_points; dz = abs(Pp(:,:,3)) < Grid.dz/2;
for e = 1:T.NumberOfElements
    plot(ax, Pp(e,dz(e,:),2)*1e3, Pp(e,dz(e,:),1)*1e3, 'k+', 'MarkerSize', 6, 'LineWidth', 1, 'Color', C(1+mod(e,2),:).*(mod(e,2)==1))
end
set(ax,'YDir','reverse'); axis(ax,'equal'); xlim(ax,[-1.7 1.7]); ylim(ax,[-0.8 1.4]); grid(ax,'off')
xlabel(ax,'y  lateral (mm)'); ylabel(ax,'x  depth (mm)')
title(ax,'Grid nodes (grey), sensor mask (blue), integration points (+)')

ax = nexttile(tl); hold(ax,'on')
C16.Type = 'C1-6D'; C16 = get_transducer_properties(C16); El = get_element_geometry(C16);
rectangle(ax,'Position',[-40 -1 80 62],'EdgeColor',C(3,:),'LineWidth',2)
text(ax,-38,45,{'grid in PROTEUS','orientation (x depth,','y lateral): contains','the array'},'Color',C(3,:),'FontSize',9)
rectangle(ax,'Position',[-6 -1 12 62],'EdgeColor',C(6,:),'LineWidth',2,'LineStyle','--')
text(ax,7,30,{'grid with y = height (±6 mm):','most of the array is outside,','pointsOutsideGrid error'},'Color',C(6,:),'FontSize',9)
plot(ax, El.y*1e3, El.x*1e3, 'k-', 'LineWidth', 3)
set(ax,'YDir','reverse'); axis(ax,'equal'); xlim(ax,[-45 45]); ylim(ax,[-5 65]); grid(ax,'on')
xlabel(ax,'y (mm)'); ylabel(ax,'x (mm)'); title(ax,'Why the grid and the points must overlap (C1-6D, top view)')
exportgraphics(f, fullfile(out,'grid_vs_points.png'), 'Resolution', 170); close(f)
fprintf('tiny array: %d points in mask slice, %d mask nodes total, %d integration points\n', nnz(m), nnz(sensor.mask), numel(Pp)/3);

%% D. Size of bubble-bubble interactions
S2 = make_tiny_settings(P); Mb2 = S2.Microbubble; Mb2.BatchSize = 5; Mb2.UseParfor = 'off';
dt = 1/100e6; kg.dt = dt; kg.k_max = pi/(S2.Medium.SpeedOfSound/(S2.Transmit.CenterFrequency*8));
p = S2.Transmit.PressureSignal; t0 = (0:numel(p)-1)/S2.Transmit.SamplingRate;
t = 0:dt:(t0(end)+4e-6); p = interp1(t0,p,t,'linear',0);
q = compute_bubble_mass_source(p, 2e-6, kg, S2.Medium, Mb2, S2.Transmit);
avg.c = S2.Medium.SpeedOfSound; avg.a = S2.Medium.AttenuationA; avg.b = S2.Medium.AttenuationB;
kg.Nt = numel(t) + round(3e-6/dt);
r = [20 50 100 200 500 1000]*1e-6;
ps = calc_scatter_attenuated(q, avg, kg, r(:));
fprintf('peak neighbour pressure (2 um bubble, 200 kPa drive):\n');
for k = 1:numel(r); fprintf('  r = %5.0f um: %6.1f kPa\n', r(k)*1e6, max(abs(ps(k,:)))/1e3); end
delete(tc); disp('figures 3 done')

function [idx, d] = knnsearch_base(ref, q)
% Nearest point of ref for each row of q (base MATLAB).
idx = zeros(size(q,1),1); d = idx;
for i = 1:size(q,1)
    [d2, idx(i)] = min((ref(:,1)-q(i,1)).^2 + (ref(:,2)-q(i,2)).^2);
    d(i) = sqrt(d2);
end
end
