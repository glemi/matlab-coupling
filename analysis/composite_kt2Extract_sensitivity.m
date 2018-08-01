function composite_kt2Extract_sensitivity(parameter, values, parName, parUnit, parScale, axScale)
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

    fig csens:kt2; clf;
    ax = axes;
    xscale log;
    yscale log;

    for k = 1:n

        fprintf('Iteration %d:\n', k);
        [Zel, f, ~, k2] = rstor_model('E33', 2.3, parameter, values(k));
        %[f, ~, i] = resample(f, log(abs(Zel)), 50e6);
        %Zel = Zel(i);

        fig csens:kt2;
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
    fig csens:kt2extr; clf;
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

    % fig csens:resfreq; clf;
    % plot(kt2_ref*100, fR/1e9, 'DisplayName', 'Resonance Frequency ');
    % plot(kt2_ref*100, fA/1e9, 'DisplayName', 'Antiresonance Frequency');
    % xscale lin; xlabel 'reference k_t^2 [%]';
    % yscale lin; ylabel 'f [GHz]';
    % title 'Resonance / Antiresonance Peaks vs k_t^2';

    %%
    idat = [data.initial];

    fig csens:fitcoeffs; clf;
    subplot(2, 3, 1);
    plot(values/parScale, [data.Q], 'DisplayName', 'Final Fit Value');
    plot(values/parScale, [idat.Q], '--', 'DisplayName', 'Initial Guess');
    title 'Q factor';
    xscale(axScale);
    legend show; legend location best;

    %fig csens:rmech; clf;
    subplot(2, 3, 2);
    plot(values/parScale, [data.R]);
    plot(values/parScale, [idat.R], '--');
    title 'Damping R_m';
    xscale(axScale);

    %fig csens:cmech; clf;
    subplot(2, 3, 3);
    plot(values/parScale, [data.Cm]*1e12, '-');
    plot(values/parScale, [data.C0]*1e12, '--');
    title 'Compliance C_m';
    xscale(axScale);

    %fig csens:c0; clf;
    subplot(2, 3, 4);
    plot(values/parScale, [data.C0]*1e12);
    plot(values/parScale, [idat.C0]*1e12, '--');
    title 'C_0 [pF]';
    xscale(axScale); yscale log;

    subplot(2, 3, 5);
    %plot(values/parScale, [data.Rs]*1e12);
    %plot(values/parScale, [idat.Rs]*1e12, '--');
    title 'R_s';
    xscale(axScale); yscale log;

    subplot(2, 3, 6);
    plot(values/parScale, fR/1e9, 'DisplayName', 'Resonance Frequency ');
    plot(values/parScale, fA/1e9, 'DisplayName', 'Antiresonance Frequency');
    xscale(axScale); xlabel 'reference k_t^2 [%]';
    yscale lin; ylabel 'f [GHz]';
    title 'f_r, f_a vs k_t^2';


    fig csens:fitcurve; clf;
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