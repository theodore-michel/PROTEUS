% Figures 2: microbubble response, propagation, DAS, end-to-end RF.
here = fileparts(mfilename('fullpath'));
repo = fullfile(here, '..', 'PROTEUS');
out  = fullfile(here, 'assets');
addpath(fullfile(repo,'tests'), fullfile(repo,'tests','helpers'));
set(groot,'defaultFigureVisible','off','defaultAxesFontSize',11,'defaultLineLineWidth',1.4);
C = [0.13 0.40 0.69; 0.85 0.37 0.01; 0.17 0.63 0.17; 0.58 0.40 0.74; 0.5 0.5 0.5; 0.8 0.1 0.3];

tc  = matlab.unittest.TestCase.forInteractiveUse;
env = tc.applyFixture(ProteusTestEnvironment);   % paths, k-Wave, temp folders
addpath(fullfile(repo,'microbubble-simulator'), fullfile(repo,'microbubble-simulator','functions'));
data = struct();                                 % numbers for the HTML animations

S = make_tiny_settings(env.PATHS);
Mb = S.Microbubble; Mb.BatchSize = 50; Mb.UseParfor = 'off';
dt = 1/100e6;
kgrid.dt = dt; kgrid.k_max = pi/(S.Medium.SpeedOfSound/(S.Transmit.CenterFrequency*8));
rho = S.Medium.Vessel.Density;

if ~isfile(fullfile(out,'das_point_scatterers.png')) || true
%% 6. Bubble response: radius and mass source for several sizes
p  = S.Transmit.PressureSignal;
t0 = (0:numel(p)-1)/S.Transmit.SamplingRate;
t  = 0:dt:(t0(end) + 4e-6);
p  = interp1(t0, p, t, 'linear', 0);
radii = [1.0 1.5 2.0 2.5 3.0]'*1e-6;
q = compute_bubble_mass_source(repmat(p, numel(radii), 1), radii, kgrid, S.Medium, Mb, S.Transmit);
V = 4/3*pi*radii.^3 + cumsum(q, 2)*dt/rho;          % q = rho dV/dt
R = (3*V/(4*pi)).^(1/3);
f = figure('Position',[0 0 1100 650]); tl = tiledlayout(f,3,1,'TileSpacing','compact');
ax = nexttile(tl); plot(ax, t*1e6, p/1e3, 'Color', C(5,:)); ylabel(ax,'p_{inc} (kPa)'); grid(ax,'on')
title(ax,'Incident pressure at the bubble (P4-1 pulse, 200 kPa)')
ax = nexttile(tl); hold(ax,'on')
for k = 1:numel(radii); plot(ax, t*1e6, R(k,:)*1e6, 'Color', C(k,:)); end
ylabel(ax,'R(t) (\mum)'); grid(ax,'on'); title(ax,'Radius from the Rayleigh-Plesset ODE with shell (Marmottant model, default)')
legend(ax, compose('R_0 = %.1f \\mum', radii*1e6), 'Location','eastoutside')
ax = nexttile(tl); hold(ax,'on')
for k = 1:numel(radii); plot(ax, t*1e6, q(k,:)*1e9, 'Color', C(k,:)); end
ylabel(ax,'q = 4\pi\rho R^2 dR/dt (\mug/s)'); xlabel(ax,'t (\mus)'); grid(ax,'on')
title(ax,'Mass source = what the bubble injects into the acoustic simulation')
legend(ax, compose('R_0 = %.1f \\mum', radii*1e6), 'Location','eastoutside')
theme(f,'light'); exportgraphics(f, fullfile(out,'bubble_response.png'), 'Resolution', 170); close(f)
data.bubble.t = t*1e6; data.bubble.p = p/1e3;
data.bubble.R = R(3,:)*1e6; data.bubble.R0 = radii(3)*1e6;

% Resonance sweep: radial excursion vs size
radii2 = linspace(0.75, 4, 27)'*1e-6;
q2 = compute_bubble_mass_source(repmat(p, numel(radii2), 1), radii2, kgrid, S.Medium, Mb, S.Transmit);
V2 = 4/3*pi*radii2.^3 + cumsum(q2, 2)*dt/rho; R2 = (3*V2/(4*pi)).^(1/3);
exc = (max(R2,[],2) - min(R2,[],2))./radii2/2;
f = figure('Position',[0 0 620 380]); ax = axes(f);
plot(ax, radii2*1e6, exc*100, 'o-', 'Color', C(1,:)); grid(ax,'on')
xlabel(ax,'R_0 (\mum)'); ylabel(ax,'(R_{max}-R_{min})/2R_0 (%)')
title(ax,'Bubble excursion vs size at 2.5 MHz, 200 kPa: resonance')
theme(f,'light'); exportgraphics(f, fullfile(out,'bubble_resonance.png'), 'Resolution', 170); close(f)

%% 7. Spherical spreading: mass source -> pressure at distance r
avg.c = S.Medium.SpeedOfSound; avg.a = S.Medium.AttenuationA; avg.b = S.Medium.AttenuationB;
kg.dt = dt; kg.Nt = numel(t) + round(25e-6/dt);
te = (0:kg.Nt-1)*dt;
r = [2 10 30]*1e-3;
ps = calc_scatter_attenuated(q(3,:), avg, kg, r(:));
f = figure('Position',[0 0 1000 360]); ax = axes(f); hold(ax,'on')
for k = 1:numel(r); plot(ax, te*1e6, ps(k,:)*r(k)/r(1), 'Color', C(k,:)); end
legend(ax, compose('r = %g mm (x r/2 mm)', r*1e3), 'Location','eastoutside'); grid(ax,'on')
xlabel(ax,'t (\mus)'); ylabel(ax,'p \cdot r/r_1 (Pa)')
title(ax,'Echo of a 2 \mum bubble: delay r/c, 1/r spreading (compensated), attenuation')
theme(f,'light'); exportgraphics(f, fullfile(out,'green_propagation.png'), 'Resolution', 170); close(f)

%% 8. Point-scatterer RF and DAS image (P4-1, plane wave), synthetic
T = reference_transducer(reset_transducer());
Tx = reset_transmit(T); Tx.LateralFocus = Inf; Tx.Angle = 0;
M.SpeedOfSound = 1540; Tx = compute_delays(Tx, T, M);
El = get_element_geometry(T);
Fs = 40e6; c = 1540;
twoway = conv(Tx.PressureSignal, T.ReceiveImpulseResponse)/T.SamplingRate;
tw = (0:numel(twoway)-1)/T.SamplingRate;
scat = [20 -6; 30 0; 40 5; 30 8]*1e-3;      % [axial lateral]
tt = (0:round(70e-6*Fs))/Fs;
RF = zeros(numel(tt), T.NumberOfElements);
for s = 1:size(scat,1)
    rr = hypot(El.x' - scat(s,1), El.y' - scat(s,2));
    t_tx = min(Tx.Delays + rr/c);                        % first arrival
    for e = 1:T.NumberOfElements
        RF(:,e) = RF(:,e) + interp1(tw, twoway, tt' - t_tx - rr(e)/c, 'linear', 0)/sqrt(rr(e));
    end
end
[~, I] = max(abs(hilbert(twoway)));  tpk = tw(I);
recon = struct('x', (5:0.1:50)*1e-3, 'y', (-14:0.1:14)*1e-3, 'c', c, 't', tt - tpk, ...
    'Fs', Fs, 'alpha_th', pi/3, 'apodization', 'cosine');
tr = struct('x_el', El.x', 'y_el', El.y', 'azimuth', El.theta', 'center', El.Center, 'radius', El.Radius);
Mdas = compute_das_matrix_general(recon, tr, struct('focus', Inf, 'steer', 0));
IMG = reshape(Mdas*reshape(hilbert(RF),[],1), numel(recon.x), numel(recon.y));
IMG = 20*log10(abs(IMG)/max(abs(IMG(:))));
f = figure('Position',[0 0 1100 470]); tl = tiledlayout(f,1,2,'TileSpacing','compact');
ax = nexttile(tl); imagesc(ax, 1:T.NumberOfElements, tt*1e6, RF/max(abs(RF(:)))); colormap(ax, gray); clim(ax,[-0.5 0.5])
xlabel(ax,'element'); ylabel(ax,'t (\mus)'); title(ax,'RF data: one column per element (hyperbolas)')
ax = nexttile(tl); imagesc(ax, recon.y*1e3, recon.x*1e3, IMG); clim(ax,[-40 0]); colormap(ax, gray); axis(ax,'image')
hold(ax,'on'); plot(ax, scat(:,2)*1e3, scat(:,1)*1e3, 'o', 'Color', C(6,:), 'MarkerSize', 12)
xlabel(ax,'y  lateral (mm)'); ylabel(ax,'x  depth (mm)'); title(ax,'Delay-and-sum image (dB): echoes summed back'); colorbar(ax)
theme(f,'light'); exportgraphics(f, fullfile(out,'das_point_scatterers.png'), 'Resolution', 170); close(f)

end
%% 9. End-to-end: tiny hybrid simulation RF (real pipeline)
St = make_tiny_settings(env.PATHS);
St.Microbubble.Number = 6; St.Acquisition.NumberOfFrames = 1;
St = update_derived_settings(St);
[o, resultDir, settingsFile] = run_rf_simulation(env, St, 'deck_tiny'); %#ok<ASGLU>
RFt = o.RF{1}{1};
tt2 = (0:size(RFt,2)-1)*o.dt*1e6;
f = figure('Position',[0 0 1100 460]); tl = tiledlayout(f,1,2,'TileSpacing','compact');
ax = nexttile(tl); hold(ax,'on'); sc = max(abs(RFt(:)));
for e = 1:size(RFt,1); plot(ax, tt2, e + 0.45*RFt(e,:)/sc, 'Color', C(1,:)); end
xlabel(ax,'t (\mus)'); ylabel(ax,'element'); grid(ax,'on')
title(ax,'Output of main\_RF (tiny 8-element run)')
ax = nexttile(tl); hold(ax,'on'); late = tt2 > 2.5; sc2 = max(abs(RFt(:,late)),[],'all');
for e = 1:size(RFt,1); plot(ax, tt2(late), e + 0.45*RFt(e,late)/sc2, 'Color', C(2,:)); end
xlabel(ax,'t (\mus)'); grid(ax,'on')
title(ax, sprintf('After the transmit pulse, x%.0f: bubble and tissue echoes', sc/sc2))
theme(f,'light'); exportgraphics(f, fullfile(out,'rf_tiny.png'), 'Resolution', 170); close(f)

save(fullfile(here,'deck_data.mat'),'data');
fid = fopen(fullfile(here,'deck_data.json'),'w'); fprintf(fid,'%s',jsonencode(data)); fclose(fid);
delete(tc)
disp('figures 2 done')
