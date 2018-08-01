clear;
%fprintf('\n\nrstor_peakdetect\n');

global debug_plot;
debug_plot = false;

%values = ([1.5 2.0 2.5]); % 1000um - 100um
values = 1.5:0.1:2.5;
%values = 2.3;
n = length(values);

values = fliplr(values);
n = length(values);

fig rmodl:kt2; clf;
ax = axes;

for k = 1:n

    fprintf('Iteration %d:\n', k);
    
    [Zel, f, ~, k2] = rstor_model('E33', values(k), 'tPiezo', 500e-9, 'tSi', 1e-9, 'tSiO', 1e-9, 'tBotE', 1e-9, 'tTopE', 1e-9, 'fstep', .1e6);
    [f, ~, i] = resample(f, log(abs(Zel)), 50e6);
    Zel = Zel(i);
    
    fig rmodl:kt2;
    hplot1 = plot_model(f, Zel);
    hplot1.DisplayName = sprintf('k_t^2 = %.2f%% (Piezo Only)', k2*100);
    %skiplegend; ax.ColorOrderIndex = ax.ColorOrderIndex -1;
    
%     [Zel, f, ~, k2] = rstor_model('E33', values(k)); %('tPiezo', values(k), 'tSi', 725-6, 'tBotE', 100e-9);
%     hplot = plot_model(f, Zel);
%     hplot.DisplayName = sprintf('k_t^2 = %.1f%%', k2*100);

    [ZelMax, ir] = max(Zel);
    [ZelMin, ia] = min(Zel);
    fr = f(ir);
    fa = f(ia);
    
    data(k) = bvdfit(f, Zel);
    
    kt2_ref(k)   = k2;
    kt2_extr1(k) = pi^2/4 * (fr - fa)/fr; % paul's book (fa / fr swapped)
    kt2_extr2(k) = pi/2 * fa/fr / tan(pi/2 * fa/fr); % finnish thesis
    kt2_extr3(k) = pi/2 * fa/fr * tan(pi/2 * (fr - fa)/ fa);
    kt2_extr4(k) = data(k).kt2;
    
    
    
    fR(k) = fr;
    fA(k) = fa;
    
%     fprintf('\nExtracted Values:\n');
%     fprintf('resonace peak:      fr = %.3f GHz\n', fr/1e9);
%     fprintf('anti-resonace peak: fa = %.3f GHz\n', fa/1e9);
%     fprintf('coupling factor:   kt2 = %.1f%%\n', kt2_extr*100);
    
    fprintf('kt2 reference: kt2 = %.1f%%\n', kt2_ref(k)*100);
    disp(data(k));

    drawnow;
    fprintf('done %d/%d\n', k, n);
end

%legend show;
%legend location best;

%%
fig rmodl:kt2extr; clf;
plot(kt2_ref*100, kt2_extr1*100, 'o-' , 'DisplayName', 'model 1');
plot(kt2_ref*100, kt2_extr2*100, 'o-' , 'DisplayName', 'model 2' );
plot(kt2_ref*100, kt2_extr3*100, 'o-' , 'DisplayName', 'model 3' );
plot(kt2_ref*100, kt2_extr4*100, 'o-' , 'DisplayName', 'BvD fit');

plot([5 20], [5 20], 'k--', 'DisplayName', 'Reference (x=y)'); skiplegend;
axis equal;
xlim([5 20]); ylim([5 20]); 
fillmarkers;

xscale lin; xlabel 'reference k_t^2 [%]';
yscale lin; ylabel 'extracted k_t^2 [%]';
title 'Extracted k_t^2 vs model k_t^2';

legend show;
legend location SouthEast;

% fig rmodl:resfreq; clf;
% plot(kt2_ref*100, fR/1e9, 'DisplayName', 'Resonance Frequency ');
% plot(kt2_ref*100, fA/1e9, 'DisplayName', 'Antiresonance Frequency');
% xscale lin; xlabel 'reference k_t^2 [%]';
% yscale lin; ylabel 'f [GHz]';
% title 'Resonance / Antiresonance Peaks vs k_t^2';

%%
fig rmodl:fitcoeffs; clf;
subplot(2, 3, 1);
plot(kt2_ref*100, [data.Q]);
title 'Q factor';
xscale lin;

%fig rmodl:rmech; clf;
subplot(2, 3, 2);
plot(kt2_ref*100, [data.R]);
title 'Damping R_m';
xscale lin;

%fig rmodl:cmech; clf;
subplot(2, 3, 3);
plot(kt2_ref*100, [data.Cm]*1e12);
title 'Compliance C_m';
xscale lin;

%fig rmodl:c0; clf;
subplot(2, 3, 4);
plot(kt2_ref*100, [data.C0]*1e12);
title 'C_0 [pF]';
xscale lin; yscale log;

subplot(2, 3, 5);
plot(kt2_ref*100, [data.Rs]*1e12);
title 'R_s';
xscale lin; yscale log;

subplot(2, 3, 6);
plot(kt2_ref*100, fR/1e9, 'DisplayName', 'Resonance Frequency ');
plot(kt2_ref*100, fA/1e9, 'DisplayName', 'Antiresonance Frequency');
xscale lin; xlabel 'reference k_t^2 [%]';
yscale lin; ylabel 'f [GHz]';
title 'f_r, f_a vs k_t^2';


fig rmodl:fitcurve; clf;
ax = axes;
for k = 1:n
    plot(data(k).f, abs(data(k).zraw), '-'); %skipcolor;
    %plot(data(k).f, abs(data(k).zfit), '--');
end
yscale log;

insetax NorthEast Normal;
a1 = gca;
for k = 1:n
    df = (fR(k) - fA(k)) / 20;
    istart = find(data(k).f < fR(k) - df, 1, 'last');
    istop  = find(data(k).f > fR(k) + df, 1, 'first');
    irange = istart:istop;
    
    plot(data(k).f(irange), abs(data(k).zraw(irange)), ':', 'LineWidth', 3); %skipcolor;
    plot(data(k).f(irange), abs(data(k).zfit(irange)), 'k', 'LineWidth', 1);
    yscale log;
    axis tight; ylim(ylim + diff(ylim)/10);
    set(gca, 'XTick', [], 'YTick', []);
end

axes(ax);
insetax SouthWest Normal;
for k = 1:n
    df = (fR(k) - fA(k)) / 5;
    istart = find(data(k).f < fA(k) - df, 1, 'last');
    istop  = find(data(k).f > fA(k) + df, 1, 'first');
    irange = istart:istop;
    
    plot(data(k).f(irange), abs(data(k).zraw(irange)), ':', 'LineWidth', 3); %skipcolor;
    plot(data(k).f(irange), abs(data(k).zfit(irange)), 'k', 'LineWidth', 1);
    yscale log;
    axis tight; ylim(ylim - diff(ylim)/10);
    set(gca, 'XTick', [], 'YTick', []);
end
axes(a1);