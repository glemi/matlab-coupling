%function fit = multiRipplefit(f, z, N, dfr, plotEnable)
function fit = HBAR_ripplefit(f, Z, N, dfr, plotEnable)
    opt N double inf;
    opt dfr double 5e6;
    opt plotEnable char off;
    hprev = gcf;
 
    warning off stats:nlinfit:IterationLimitExceeded;
    warning off stats:nlinfit:IllConditionedJacobian;
    warning off stats:nlinfit:ModelConstantWRTParam;
    
    df = mean(diff(f(1:10))); % sample spacing
    di = round(dfr/df);       % 
    [~, ia] = findpeaks(abs(Z), 'MinPeakDistance', di);
    [~, ir] = findpeaks(1./abs(Z),'MinPeakDistance', di);
    
    fr = f(ir);
    

    N = min([N length(ir)]);     % target number of ripples to fit
    nskip = round(length(ir)/(N*1.2)); % number of ripples to skip

    n = 3;                       % number of ripples in each window
    d = round(median(diff(ir))); % average #samples between two fr
    i = ir - round(d/3);         % starting sample of each window
    i = i(i>0);                  % can't be negative...
    N = length(i)-n;             % number of windows
    
    if N < 100
        disp('');
    end
    
    irange = i(1):i(n);
    %bvdMultiLsqFit(f(irange), Z(irange), 100, 'Algorithm', 'levenberg-marquardt');
    bvdMultiFit(f(irange), Z(irange), 100);

    switch lower(plotEnable); case 'on'
        hf = fig('Ripple Fit Progress'); clf;
        hf.MenuBar = 'none'; hf.ToolBar = 'none';
        hf.Position(4) = 280;
    end
    
    tic;
    for k = 1:nskip:N
        irange = i(k):i(k+n);
        frange = f(irange);
        zrange = Z(irange);
        
        tic;
        F = bvdMultiFit(frange, zrange, 51, 'recycle', 'MaxIter', 4);
        %F = bvdMultiLsqFit(frange, zrange, 51, 'recycle', 'MaxIterations', 3);
        F.toc = toc;
        fit(k) = F;
        
        if fit(k).valid
            ffit = f(i(k)):20e3:f(i(k+1));
            zfit = fit(k).model(ffit);

            %% plot 
            switch lower(plotEnable); case 'on'
                fig('Ripple Fit Progress'); cla;
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
    %toc;
    %movie2avi(frames(2:end), 'fitmovie1.avi');
    switch lower(plotEnable); case 'on'
        delete(hf);
    end
    
    warning on stats:nlinfit:IterationLimitExceeded;
    warning on stats:nlinfit:IllConditionedJacobian;
    
    try; figure(hprev);end %#ok

end
