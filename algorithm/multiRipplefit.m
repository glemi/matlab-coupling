%function fit = multiRipplefit(f, z, n, dfr, plotEnable)
function fit = multiRipplefit(f, z, n, dfr, plotEnable)
    opt n double 1;
    opt dfr double 5e6;
    opt plotEnable char off;
    hprev = gcf;
  
    df = mean(diff(f(1:10))); % sample spacing
    di = round(dfr/df);       % 
    [~, ia] = findpeaks(abs(z), 'MinPeakDistance', di);
    [~, ir] = findpeaks(1./abs(z),'MinPeakDistance', di);


    d = round(median(diff(ir))); % average #samples between two fr
    i = ir - round(d/3);         % starting sample of each window
    i = i(i>0);                  % can't be negative...
    N = length(i)-n;             % number of windows
    
    if N < 100
        disp('');
    end
    
    irange = i(1):i(n);
    bvdMultiFit(f(irange), z(irange), 100);

    warning off stats:nlinfit:IterationLimitExceeded;
    warning off stats:nlinfit:IllConditionedJacobian;
    warning off all;
    
    switch lower(plotEnable); case 'on'
        hf = fig('Ripple Fit Progress'); clf;
        hf.MenuBar = 'none'; hf.ToolBar = 'none';
        hf.Position(4) = 280;
    end
    
    tic;
    for k = 1:1:N
        irange = i(k):i(k+n);
        frange = f(irange);
        zrange = z(irange);
        
        tic;
        F = bvdMultiFit(frange, zrange, 51, 'recycle', 'MaxIter', 1);
        F.toc = toc;
        fit(k) = F;
        
        if fit(k).valid
            ffit = f(i(k)):20e3:f(i(k+1));
            zfit = fit(k).model(ffit);

            %% plot 
            switch lower(plotEnable); case 'on'
                cla;
                plot(frange, zrange, '.');
                ffit = frange(1):2e3:frange(end);
                zfit = F.model(ffit);

                ylim([min(zfit)*0.95 max(zfit)*1.05]);
                xlim([min(frange) max(frange)]);
                set(gca, 'XTick', [], 'YTick', []);
                progress = 100*(frange(1)-f(1))/(f(end)-f(1));
                title(sprintf('Fitting Ripples: %.0f%% (%.2fGHz)', progress, F.fr/1e9));

                plot(ffit, zfit);
                gap = diff(ylim)*0.1;
                plot(F.fmin, F.zmin-gap, '^', 'MarkerSize', 6);
                plot(F.fmax, F.zmax+gap, 'v', 'MarkerSize', 6); 
                plot([F.fr F.fa], [F.zr F.za], 'kx');
                fillmarkers;
                drawnow;
            end
            %frames(k) = getframe(hf);
        end
    end    
    toc;
    %movie2avi(frames(2:end), 'fitmovie1.avi');
    switch lower(plotEnable); case 'on'
        delete(hf);
    end
    
    warning on stats:nlinfit:IterationLimitExceeded;
    warning on stats:nlinfit:IllConditionedJacobian;
    
    return;
    fig ripplefit:coeffs; 
    
        subplot(3, 2, 1);
        fit     = fit([fit.valid]);
        fGHz    = [fit.fr]/1e9;
        initial = [fit.initial];

        %plot(fGHz, [initial.Q]);
        plot(fGHz, [fit.Q]);
        title Q

        subplot(3, 2, 2);
        %plot(fGHz, [initial.R]);
        plot(fGHz, [fit.R]);
        title R

        subplot(3, 2, 3);
        %plot(fGHz, [initial.C0]);
        plot(fGHz, [fit.C0]);
        title C0

        subplot(3, 2, 4);
        %plot(fGHz, [initial.dfr]);
        plot(fGHz, [fit.dfr]);
        title fr
        
        subplot(3, 2, 5);
        plot(fGHz, [fit.Cm]);
        title C_m
        
        subplot(3, 2, 6);
        plot(fGHz, [fit.Lm]);
        title L_m

    try; figure(hprev);end

%     fprintf('\n');
%     fprintf('Q = %.2e\n', fit(k).Q);
%     fprintf('Cm = %.2e\n', fit(k).Cm);
%     fprintf('Lm = %.2e\n', fit(k).Lm);
%     fprintf('C0 = %.2e\n', fit(k).C0);
%     fprintf('C0 = %.2e\n', fit(k).Rm);
end
