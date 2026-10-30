% Figures 1: geometry, transducers, transmit chain, delays.
% Run from this folder. Uses the PROTEUS repo next to it (../PROTEUS).
repo = fullfile(fileparts(mfilename('fullpath')), '..', 'PROTEUS');
out  = fullfile(fileparts(mfilename('fullpath')), 'assets');
if ~isfolder(out); mkdir(out); end
addpath(repo, fullfile(repo,'GUIfunctions'), fullfile(repo,'acoustic-module'), ...
    fullfile(repo,'delay-and-sum'), fullfile(repo,'tests','helpers'));
set(groot,'defaultFigureVisible','off','defaultAxesFontSize',11, ...
    'defaultFigureColor','w','defaultLineLineWidth',1.4);
C = [0.13 0.40 0.69; 0.85 0.37 0.01; 0.17 0.63 0.17; 0.58 0.40 0.74; 0.5 0.5 0.5];
save_fig = @(f, name) export_light(f, fullfile(out,[name '.png']));

PATHS.GeometriesPath = fullfile(repo,'geometry_data');
P41 = reset_transducer();
C16.Type = 'C1-6D'; C16 = get_transducer_properties(C16);
Medium = reset_medium();

%% 1. Domain: transducer, beam projection, geometry, simulation box
for T = {P41, C16}
    T = T{1};
    S.Transducer = T; S.Medium = Medium;
    S.Transmit = reset_transmit(T);
    S.Microbubble = reset_microbubble(); S.Acquisition = reset_acquisition();
    S.SimulationParameters = reset_simulation_parameters();
    G.Folder = 'renal_tree'; S.Geometry = reset_geometry(G, PATHS);
    if strcmp(T.Configuration,'Curved')
        El = get_element_geometry(T);
        S.Geometry.startDepth = max(S.Geometry.startDepth, El.Apex + 5e-3);
    end
    S = update_derived_settings(S);
    D = S.Geometry.Domain;
    f = figure('Position',[0 0 760 520]); ax = axes(f); hold(ax,'on')
    Sv = D.TransducerSurface*1e3; P = D.TransducerProjection*1e3;
    NP = size(Sv,1)/4; FT = reshape(1:4*NP,4,NP)';
    patch(ax,'Vertices',Sv,'Faces',FT,'FaceColor',C(1,:),'EdgeColor','none')
    V = [Sv; P]; c = [4*NP-3, 2, 3, 4*NP]; p = 4*NP + (1:4);
    FB = [c(1) c(2) p(2) p(1); c(2) c(3) p(3) p(2); c(3) c(4) p(4) p(3); c(4) c(1) p(1) p(4)];
    patch(ax,'Vertices',V,'Faces',FB,'FaceColor',C(1,:),'FaceAlpha',0.12,'EdgeColor',C(1,:),'EdgeAlpha',0.3)
    Vd = D.Vertices*1e3; FD = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
    patch(ax,'Vertices',Vd,'Faces',FD,'FaceColor','none','EdgeColor',[0 0 0],'LineStyle','--')
    % Geometry bounding box (rotated), centred at Geometry.Center
    Pu = [-1 -1 -1 -1 1 1 1 1; -1 -1 1 1 -1 -1 1 1; -1 1 -1 1 -1 1 -1 1]/2;
    BB = (S.Geometry.Rotation*(S.Geometry.BoundingBox.Diagonal.*Pu) + S.Geometry.Center)'*1e3;
    k = convhull(BB); trisurf(k, BB(:,1),BB(:,2),BB(:,3),'FaceColor',C(2,:),'FaceAlpha',0.18,'EdgeColor','none','Parent',ax)
    axis(ax,'equal'); grid(ax,'on'); view(ax,[-38 22])
    xlabel(ax,'x  axial / depth (mm)'); ylabel(ax,'y  lateral (mm)'); zlabel(ax,'z  elevation (mm)')
    title(ax, sprintf('%s: transducer, beam, geometry box, simulation domain', T.Type))
    set(ax,'ZDir','normal')
    save_fig(f, ['domain_' strrep(T.Type,'-','')]); close(f)
end

%% 2. Element geometry and integration points
Grid.dx = 0.15e-3; Grid.dy = 0.15e-3; Grid.dz = 0.15e-3;
f = figure('Position',[0 0 1450 430]);
tl = tiledlayout(f,1,3,'TileSpacing','compact');
for k = 1:2
    if k == 1, T = P41; else, T = C16; end
    T.IntegrationDensity = 1;
    T = get_transducer_integration_points(T, Grid);
    El = get_element_geometry(T);
    ax = nexttile(tl); hold(ax,'on')
    Pp = T.integration_points;
    if k == 1
        % Linear array: face view (lateral y vs elevation z) of 8 elements
        for e = 45:52
            plot(ax, squeeze(Pp(e,:,2))*1e3, squeeze(Pp(e,:,3))*1e3, '.', 'Color', C(1+mod(e,2),:), 'MarkerSize', 5)
            yl = El.y(e)*1e3 + [-1 1]*T.ElementWidth/2*1e3;
            rectangle(ax, 'Position', [yl(1), -T.ElementHeight/2*1e3, diff(yl), T.ElementHeight*1e3], 'EdgeColor', C(5,:))
        end
        grid(ax,'on'); xlim(ax,[-1.5 1.5]); xlabel(ax,'y  lateral (mm)'); ylabel(ax,'z  elevation (mm)')
        title(ax, sprintf('P4-1 face view (x = 0), 8 of %d elements, %d points each', ...
            T.NumberOfElements, size(Pp,2)))
    else
        for e = 1:T.NumberOfElements
            plot(ax, squeeze(Pp(e,:,2))*1e3, squeeze(Pp(e,:,1))*1e3, '.', 'Color', C(1+mod(e,2),:), 'MarkerSize', 4)
        end
        th = linspace(-El.ThetaMax*1.15, El.ThetaMax*1.15, 200);
        plot(ax, El.Radius*sin(th)*1e3, (El.Center + El.Radius*cos(th))*1e3, ':', 'Color', C(5,:))
        plot(ax, [0 0], [El.Center*1e3, El.Apex*1e3], '--', 'Color', C(5,:))
        plot(ax, 0, El.Center*1e3, 'k+', 'MarkerSize', 10, 'LineWidth', 2)
        yline(ax, 0, ':', 'x = 0 (outer edges)');
        text(ax, 1.5, El.Center*1e3, 'centre of curvature', 'FontSize', 9)
        text(ax, 1.5, El.Apex*1e3 - 2, sprintf('apex x = %.1f mm', El.Apex*1e3), 'FontSize', 9)
        set(ax,'YDir','reverse'); axis(ax,'equal'); grid(ax,'on')
        ylim(ax, [El.Center*1e3-3, 15])
        xlabel(ax,'y  lateral (mm)'); ylabel(ax,'x  depth (mm)')
        title(ax, sprintf('C1-6D top view (z = 0): %d elements along an arc', T.NumberOfElements))
    end
end
% C1-6D angled view: the same points in 3-D
ax = nexttile(tl); hold(ax,'on')
Pp = T.integration_points; stride = 2;
for e = 1:stride:T.NumberOfElements
    plot3(ax, squeeze(Pp(e,:,2))*1e3, squeeze(Pp(e,:,1))*1e3, squeeze(Pp(e,:,3))*1e3, '.', 'Color', C(1+mod((e-1)/stride,2),:), 'MarkerSize', 2)
end
set(ax,'YDir','reverse'); axis(ax,'equal'); grid(ax,'on'); view(ax,[-35 28])
xlabel(ax,'y  lateral (mm)'); ylabel(ax,'x  depth (mm)'); zlabel(ax,'z  elevation (mm)')
title(ax,'C1-6D angled view: curved surface (every 2nd element)')
save_fig(f, 'integration_points'); close(f)

% Zoom: one element in 3D, points in y and z
T = P41; T.IntegrationDensity = 1.5; T = get_transducer_integration_points(T, Grid);
f = figure('Position',[0 0 560 420]); ax = axes(f); hold(ax,'on')
for e = 47:50
    Pp = squeeze(T.integration_points(e,:,:));
    plot3(ax, Pp(:,1)*1e3, Pp(:,2)*1e3, Pp(:,3)*1e3, '.', 'Color', C(1+mod(e,2),:), 'MarkerSize', 7)
end
xlabel(ax,'x (mm)'); ylabel(ax,'y  lateral (mm)'); zlabel(ax,'z  elevation (mm)')
title(ax,'P4-1: integration points of 4 elements (density 1.5)'); grid(ax,'on'); view(ax,[-60 20])
save_fig(f, 'integration_zoom'); close(f)

%% 3. Off-grid point: band-limited delta on the grid
G2.x = (-10:10)*1e-4; G2.y = (-10:10)*1e-4; G2.z = (-10:10)*1e-4;
G2.dx = 1e-4; G2.dy = 1e-4; G2.dz = 1e-4; G2.Nx = 21; G2.Ny = 21; G2.Nz = 21;
G2.full_size = [21; 21; 21]; G2.sensor_on_grid = false;
pt = [0.33e-4 -0.21e-4 0];
[pt2, ~, idx] = voxelize_media_points(pt, G2);
sens.mask = false(21,21,21);
[sens, w] = update_sensor_fast(sens, pt2, idx, G2, false);
W = zeros(21,21,21); W(find(sens.mask)) = full(w); %#ok<FNDSB>
f = figure('Position',[0 0 980 400]); tl = tiledlayout(f,1,2,'TileSpacing','compact');
ax = nexttile(tl); imagesc(ax, G2.y*1e3, G2.x*1e3, W(:,:,11)); axis(ax,'image'); colorbar(ax)
hold(ax,'on'); plot(ax, pt(2)*1e3, pt(1)*1e3, 'k+', 'MarkerSize', 12, 'LineWidth', 2)
xlabel(ax,'y (mm)'); ylabel(ax,'x (mm)'); title(ax,'weights of one off-grid point (z = 0 slice)')
colormap(ax, redblue_map());
cl = max(abs(W(:))); clim(ax, [-cl cl]);
ax = nexttile(tl); stem(ax, G2.x*1e3, squeeze(W(:,11,11)), 'filled', 'Color', C(1,:)); hold(ax,'on')
xx = linspace(-1,1,400)*1e-3; plot(ax, xx*1e3, max(abs(squeeze(W(:,11,11))))*sinc((xx-pt(1))/1e-4), ':', 'Color', C(2,:))
xline(ax, pt(1)*1e3, 'k'); grid(ax,'on'); xlabel(ax,'x (mm)'); title(ax,'1-D cut: truncated sinc (±4 grid points)')
save_fig(f, 'offgrid_weights'); close(f)

%% 4. Transmit chain: voltage -> transmit impulse response -> pressure
Tx = reset_transmit(P41);
f = figure('Position',[0 0 1100 600]); tl = tiledlayout(f,2,3,'TileSpacing','compact');
tV = (0:numel(Tx.VoltageSignal)-1)/Tx.SamplingRate*1e6;
tI = (0:numel(P41.TransmitImpulseResponse)-1)/P41.SamplingRate*1e6;
tP = (0:numel(Tx.PressureSignal)-1)/Tx.SamplingRate*1e6;
ax = nexttile(tl); plot(ax, tV, Tx.VoltageSignal, 'Color', C(1,:)); title(ax,'1. Voltage (three-level)'); xlabel(ax,'t (\mus)'); grid(ax,'on')
ax = nexttile(tl); plot(ax, tI, P41.TransmitImpulseResponse/max(abs(P41.TransmitImpulseResponse)), 'Color', C(2,:)); title(ax,'2. Transmit impulse response'); xlabel(ax,'t (\mus)'); grid(ax,'on')
ax = nexttile(tl); plot(ax, tP, Tx.PressureSignal/1e3, 'Color', C(3,:)); title(ax,'3. Pressure = 1 * 2 (200 kPa)'); xlabel(ax,'t (\mus)'); ylabel(ax,'kPa'); grid(ax,'on')
sp = @(x,fs) deal((0:2^14-1)*fs/2^14/1e6, 20*log10(abs(fft(x,2^14))/max(abs(fft(x,2^14)))));
[fa, A] = sp(Tx.VoltageSignal, Tx.SamplingRate); ax = nexttile(tl); plot(ax, fa, A, 'Color', C(1,:)); xlim(ax,[0 8]); ylim(ax,[-40 0]); xlabel(ax,'f (MHz)'); ylabel(ax,'dB'); grid(ax,'on')
[fa, A] = sp(P41.TransmitImpulseResponse, P41.SamplingRate); ax = nexttile(tl); plot(ax, fa, A, 'Color', C(2,:)); xlim(ax,[0 8]); ylim(ax,[-40 0]); xlabel(ax,'f (MHz)'); grid(ax,'on')
xline(ax, [1.5 3.5], '--', {'1.5','3.5 MHz'})
[fa, A] = sp(Tx.PressureSignal, Tx.SamplingRate); ax = nexttile(tl); plot(ax, fa, A, 'Color', C(3,:)); xlim(ax,[0 8]); ylim(ax,[-40 0]); xlabel(ax,'f (MHz)'); grid(ax,'on')
save_fig(f, 'transmit_chain'); close(f)

%% 5. Transmit delays
M.SpeedOfSound = 1540;
cases = {P41, Inf, 0, 'plane wave'; P41, Inf, 15, 'steered 15°'; P41, 40e-3, 0, 'focused 40 mm'; ...
    P41, -20e-3, 0, 'diverging (virtual source -20 mm)'; C16, Inf, 0, 'curved: diverging from centre'; C16, 50e-3, 0, 'curved: focused 50 mm'};
f = figure('Position',[0 0 1100 560]); tl = tiledlayout(f,2,3,'TileSpacing','compact');
for k = 1:size(cases,1)
    T = cases{k,1}; X.LateralFocus = cases{k,2}; X.Angle = cases{k,3};
    X = compute_delays(X, T, M); El = get_element_geometry(T);
    X.Delays(abs(X.Delays) < 1e-15) = 0;
    ax = nexttile(tl); yyaxis(ax,'left'); plot(ax, El.y*1e3, X.Delays*1e6, '.-', 'Color', C(1,:)); ylabel(ax,'delay (\mus)')
    yyaxis(ax,'right'); plot(ax, El.y*1e3, El.x*1e3, '-', 'Color', C(5,:)); ylabel(ax,'element x (mm)'); set(ax,'YDir','reverse')
    xlabel(ax,'y  lateral (mm)'); title(ax, [T.Type ': ' cases{k,4}]); grid(ax,'on')
end
save_fig(f, 'delays'); close(f)

disp('figures 1 done')

function export_light(f, file)
theme(f, 'light');
exportgraphics(f, file, 'Resolution', 170);
end

function m = redblue_map()
n = 128; r = [linspace(0.13,1,n) ones(1,n)]'; g = [linspace(0.4,1,n) linspace(1,0.2,n)]'; b = [ones(1,n) linspace(1,0.13,n)]';
m = [r g b];
m(1:n,1) = linspace(0.13,1,n); m(1:n,2) = linspace(0.4,1,n); m(1:n,3) = 1;
m(n+1:end,1) = 1; m(n+1:end,2) = linspace(1,0.2,n); m(n+1:end,3) = linspace(1,0.13,n);
end
