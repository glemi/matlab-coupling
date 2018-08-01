%function [fr, fa, fit] = ripplefit(f, z, dfr)
function [fr, fa, fit] = ripplefit(f, z, n, dfr)
    opt n double 1;
    opt dfr double 5e6;
  
    df = mean(diff(f(1:10))); % sample spacing
    di = round(dfr/df);       % 
    [~, ia] = findpeaks(abs(z), 'MinPeakDistance', di);
    [~, ir] = findpeaks(1./abs(z),'MinPeakDistance', di);


    d = round(median(diff(ir)));
    m = round((n-1)/2);
    i = ir - m - round(d/3);
    
    N   = length(i);
    fr  = zeros(1,N);  
    dfr = zeros(1,N);
    cr  = zeros(1,N);
   
    irange = i(1):i(2);
    bvdMultiFit(f(irange), z(irange), n);

    warning off stats:nlinfit:IterationLimitExceeded;
    warning off stats:nlinfit:IllConditionedJacobian;
    
    for k = 1:1:N-1
        irange = i(k):i(k+1);
        frange = f(irange);
        zrange = z(irange);
            
        fr(k) = NaN;
        fa(k) = NaN;
        fit(k) = bvdfit(frange, zrange, 'MaxIter', 30);
        
        if fit(k).valid
            ffit = f(i(k)):20e3:f(i(k+1));
            zfit = fit(k).model(ffit);

            [~, jr] = min(zfit);
            [~, ja] = max(zfit);
            fr(k) = ffit(jr);
            fa(k) = ffit(ja);
        
            %% plot 
            %continue;
            cla;
            plot(frange, zrange, '.');
            ffit = frange(1):2e3:frange(end);
            zfit = fit(k).model(ffit);
            plot(ffit, zfit);
            plot(fr(k), fit(k).model(fr(k)), 'x');
            plot(fa(k), fit(k).model(fa(k)), 'x');
            ylim([min(zfit)*0.9 max(zfit)*1.1]);
            %[min(zfit)*0.9 max(zfit)*1.1]
            xlim([min(frange) max(frange)]);
            set(gca, 'XTick', [], 'YTick', []);
            title(sprintf('Fitting Ripples: %.0f%% (%.2fGHz)', 100*frange(1)./max(f), fr(k)/1e9));
            drawnow;
        end
    end    
%     fprintf('\n');
%     fprintf('Q = %.2e\n', fit(k).Q);
%     fprintf('Cm = %.2e\n', fit(k).Cm);
%     fprintf('Lm = %.2e\n', fit(k).Lm);
%     fprintf('C0 = %.2e\n', fit(k).C0);
%     fprintf('C0 = %.2e\n', fit(k).Rm);

    warning on stats:nlinfit:IterationLimitExceeded;
    warning on stats:nlinfit:IllConditionedJacobian;
end
