function piezo_kt2Extract_sensitivity(parameter, values, parName, parUnit, parScale, axScale)
    opt axScale char lin;

    global debug_plot;
    debug_plot = false;

    
    %values = ([1.5 2.0 2.5]); % 1000um - 100um
%     values = (400:100:1000)*1e-9;
%     parameter = 'tPiezo';
%     parName = 'Piezolayer Thickness';
%     parUnit = 'nm';
%     parScale = 1e9;

    values = values*parScale;

    n = length(values);

    %values = fliplr(values);
    n = length(values);

    fig psens:kt2; clf;
    ax = axes;
    xscale log;
    yscale log;
    
    fig psens:kt2Y; clf; xscale log; yscale lin; title 'Admittance Y';
    fig psens:kt2phpase; clf; xscale log; yscale lin; title 'Phase of Z';
    
    %args = {'tSi', 1e-9, 'tSiO', 1e-9, 'tBotE', 1e-9, 'tTopE', 1e-9, 'fstep', .1e6};
    args = {'tSi', 0, 'tSiO', 0, 'tBotE', 0, 'tTopE', 0, 'fstep', .1e6};
    args = [args {'Qpiezo' 10}];

    for k = 1:n

        fprintf('Iteration %d:\n', k);
        
        [Zel, f, ~, k2] = rstor_model('E33', 2.3, parameter, values(k), args{:});
        %[f, ~, i] = resample(f, log(abs(Zel)), 50e6);
        %Zel = Zel(i);

        fig psens:kt2;
        hplot1 = plot_model(f, Zel);
        hplot1.DisplayName = sprintf('k_t^2 = %.2f%% (Piezo Only)', k2*100);
        %skiplegend; ax.ColorOrderIndex = ax.ColorOrderIndex -1;
        
        fig psens:kt2Y;
        hplot1 = plot(f, 1./imag(Zel));
        hplot1.DisplayName = sprintf('k_t^2 = %.2f%% (Piezo Only)', k2*100);
        
        fig psens:kt2phpase;
        hplot1 = plot(f, angle(Zel));
        hplot1.DisplayName = sprintf('k_t^2 = %.2f%% (Piezo Only)', k2*100);
         
%         ph = angle(Zel);
%         iph0 = crossing(ph);
%         fph0r(k) = f(iph0(1));
%         fph0a(k) = f(iph0(2));
        

    %     [Zel, f, ~, k2] = rstor_model('E33', values(k)); %('tPiezo', values(k), 'tSi', 725-6, 'tBotE', 100e-9);
    %     hplot = plot_model(f, Zel);
    %     hplot.DisplayName = sprintf('k_t^2 = %.1f%%', k2*100);

%         [ZelMax, ia] = max(Zel);
%         [ZelMin, ir] = min(Zel);
%         fr = f(ir);
%         fa = f(ia);
        [fr, fa] = resonance_find(f, abs(Zel), 1);
        fr = fr(1);
        fa = fa(1);

        data(k) = bvdfit(f, abs(Zel), 'MaxIter', 30);

        kt2_ref(k)   = k2;
        kt2_extr1(k) = pi^2/4 * (fa - fr)/fa; % paul's book (fa / fr swapped)
        kt2_extr2(k) = pi/2 * fr/fa / tan(pi/2 * fr/fa); % finnish thesis
        kt2_extr3(k) = pi/2 * fr/fa * tan(pi/2 * (fa - fr)/ fa);
        kt2_extr4(k) = data(k).kt2;

        %kt2_ph0extr(k) = 

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
    
    fig psens:kt2phpase; crosshair;

    %%
    fig psens:kt2extr; clf;
    plot(values/parScale, kt2_extr1*100, 'o-' , 'DisplayName', 'model 1');
    plot(values/parScale, kt2_extr2*100, 'o-' , 'DisplayName', 'model 2' );
    %plot(kt2_ref*100, kt2_extr3*100, 'o-' , 'DisplayName', 'model 3' );
    plot(values/parScale, kt2_extr4*100, 'o-' , 'DisplayName', 'BvD fit');

    xrefline = [min(values) max(values)]/parScale;
    yrefline = [kt2_ref(1) kt2_ref(1)]*100;
    plot(xrefline, yrefline, 'k--', 'DisplayName', 'Actual k_t^2'); skiplegend;
    %axis equal;
    %xlim([5 20]); ylim([5 20]); 
    fillmarkers;

    xscale(axScale); xlabel(sprintf('%s [%s]', parName, parUnit));
    yscale lin; ylabel 'extracted k_t^2 [%]';
    title(['Extracted k_t^2 vs ' parName]);

    legend show;
    legend location SouthEast;

    % fig psens:resfreq; clf;
    % plot(kt2_ref*100, fR/1e9, 'DisplayName', 'Resonance Frequency ');
    % plot(kt2_ref*100, fA/1e9, 'DisplayName', 'Antiresonance Frequency');
    % xscale lin; xlabel 'reference k_t^2 [%]';
    % yscale lin; ylabel 'f [GHz]';
    % title 'Resonance / Antiresonance Peaks vs k_t^2';

    %%
    idat = [data.initial];

    fig psens:fitcoeffs; clf;
    subplot(2, 3, 1);
    plot(values/parScale, [data.Q], 'DisplayName', 'Final Fit Value');
    plot(values/parScale, [idat.Q], '--', 'DisplayName', 'Initial Guess');
    title 'Q factor';
    xscale(axScale);
    legend show; legend location best;

    %fig psens:rmech; clf;
    subplot(2, 3, 2);
    plot(values/parScale, [data.R]);
    plot(values/parScale, [idat.R], '--');
    title 'Damping R_m';
    xscale(axScale);

    %fig psens:cmech; clf;
    subplot(2, 3, 3);
    plot(values/parScale, [data.Cm]*1e12, '--');

    title 'Compliance C_m';
    xscale(axScale);

    %fig psens:c0; clf;
    subplot(2, 3, 4);
    plot(values/parScale, [data.C0]*1e12);
    plot(values/parScale, [data.Cm]*1e12, '--');
    title 'C_0 [pF]';
    xscale(axScale); yscale log;

    subplot(2, 3, 5);
    %plot(values/parScale, [data.Rs]*1e12);
    %plot(values/parScale, [idat.Rs]*1e12, '--');
    %title 'R_s';
    plot(values/parScale, [data.Cm]./[data.C0], '-');
    title 'C_m / C_0'
    xscale(axScale); yscale log;

    subplot(2, 3, 6);
    plot(values/parScale, fR/1e9, 'DisplayName', 'f_r from fit');
    plot(values/parScale, fA/1e9, 'DisplayName', 'f_a from fit');
%     plot(values/parScale, fph0r/1e9, 'DisplayName', 'f_r from \phi=0');
%     plot(values/parScale, fph0a/1e9, 'DisplayName', 'f_a from \phi=0');
    xscale(axScale); xlabel 'reference k_t^2 [%]';
    yscale lin; ylabel 'f [GHz]';
    title 'f_r, f_a vs k_t^2';
    legend show; legend location best;

    fig psens:fitcurve; clf;
    ax = axes;
    for k = 1:n
        plot(data(k).f, abs(data(k).zraw), '-'); %skipcolor;
        %plot(data(k).f, abs(data(k).zfit), '--');
    end
    yscale log;
    xscale log;

    return;    
    insetax NorthEast Normal;
    a1 = gca;
    for k = 1:n
        df = (fR(k) - fA(k)) / 20;
        irange = data(k).f < fR(k) - df & data(k).f > fR(k) + df;

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
        irange = data(k).f < fA(k) - df & data(k).f > fA(k) + df;

        plot(data(k).f(irange), abs(data(k).zraw(irange)), ':', 'LineWidth', 3); %skipcolor;
        plot(data(k).f(irange), abs(data(k).zfit(irange)), 'k', 'LineWidth', 1);
        yscale log;
        axis tight; ylim(ylim - diff(ylim)/10);
        set(gca, 'XTick', [], 'YTick', []);
    end
    axes(a1);
end