clear;
%fprintf('\n\nrstor_peakdetect\n');

global debug_plot;
debug_plot = false;


values = 1.5:0.1:2.5;
values = fliplr(values);
n = length(values);



fig cextr:kt2; clf;
ax = axes; 
xscale log;
yscale log;

for k = 1:n

    fprintf('Iteration %d:\n', k);
    
    [Zel, f, ~, k2] = rstor_model('E33', values(k), 'tPiezo', 500e-9, 'tSi', 725e-6, 'tSiO', 200e-9, 'tBotE', 200e-9, 'tTopE', 100e-9, 'fstep', .1e6);
    
    fig cextr:kt2;
    hplot1 = plot_model(f, Zel);
    hplot1.DisplayName = sprintf('k_t^2 = %.2f%% (Piezo Only)', k2*100);

    
    fi = linspace(min(f), max(f), 10000);
    [Zel, Zup, Zlo] = deripple(f, abs(Zel), fi); 
    f = fi;
    
    %[fr, fa, ir, ia] = resonance_find(f, Zel);
    [fr, ~, ir, ~] = resonance_find(f, Zlo);
    [~, fa, ~, ia] = resonance_find(f, Zup);
    
    data(k) = bvdfit(f, Zel);
    
    kt2_ref(k)   = k2;
    kt2_extr1(k) = pi^2/4 * (fa - fr)/fa; % paul's book (fa / fr swapped)
    kt2_extr2(k) = pi/2 * fr/fa / tan(pi/2 * fr/fa); % finnish thesis
    %kt2_extr3(k) = pi/2 * fr/fa * tan(pi/2 * (fa - fr)/ fa);
    kt2_extr4(k) = data(k).kt2;

    fR(k) = fr(1);
    fA(k) = fa(1);
    zR(k) = Zel(ir(1));
    zA(k) = Zel(ia(1));
    
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
fig cextr:kt2extr; clf;
plot(kt2_ref*100, kt2_extr1*100, 'o-' , 'DisplayName', 'model 1');
plot(kt2_ref*100, kt2_extr2*100, 'o-' , 'DisplayName', 'model 2' );
%plot(kt2_ref*100, kt2_extr3*100, 'o-' , 'DisplayName', 'model 3' );
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

% fig cextr:resfreq; clf;
% plot(kt2_ref*100, fR/1e9, 'DisplayName', 'Resonance Frequency ');
% plot(kt2_ref*100, fA/1e9, 'DisplayName', 'Antiresonance Frequency');
% xscale lin; xlabel 'reference k_t^2 [%]';
% yscale lin; ylabel 'f [GHz]';
% title 'Resonance / Antiresonance Peaks vs k_t^2';

%%
idat = [data.initial];

fig cextr:fitcoeffs; clf;
subplot(2, 3, 1);
plot(kt2_ref*100, [data.Q]);
plot(kt2_ref*100, [idat.Q], '--');
title 'Q factor';
xscale lin;

%fig cextr:rmech; clf;
subplot(2, 3, 2);
plot(kt2_ref*100, [data.R]);
plot(kt2_ref*100, [idat.R], '--');
title 'Damping R_m';
xscale lin;

%fig cextr:cmech; clf;
subplot(2, 3, 3);
plot(kt2_ref*100, [data.Cm]*1e12);
title 'Compliance C_m';
xscale lin;

%fig cextr:c0; clf;
subplot(2, 3, 4);
plot(kt2_ref*100, [data.C0]*1e12);
plot(kt2_ref*100, [idat.C0]*1e12, '--');
title 'C_0 [pF]';
xscale lin; yscale log;

subplot(2, 3, 5);
plot(kt2_ref*100, [data.Rs]);
title 'R_s';
xscale lin; yscale log;

subplot(2, 3, 6);
plot(kt2_ref*100, fR/1e9, 'DisplayName', 'Resonance Frequency ');
plot(kt2_ref*100, fA/1e9, 'DisplayName', 'Antiresonance Frequency');
xscale lin; xlabel 'reference k_t^2 [%]';
yscale lin; ylabel 'f [GHz]';
title 'f_r, f_a vs k_t^2';

%return;

%% Plot fit curve
fig cextr:fitcurve; clf;
ax = axes;
for k = 1:n
    plot(data(k).f, abs(data(k).zraw), '-'); skipcolor;
    plot(data(k).f, abs(data(k).zfit), '--');
    plot([fR(k) NaN fA(k)], [zR(k) NaN zA(k)], 'k+');
end
yscale log;
xscale log;

return;

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