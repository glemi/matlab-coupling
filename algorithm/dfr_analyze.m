function data = dfr_analyze(fr, dfr, option)
    opt option char noplot;
    
    [fr, dfr] = rmoutliers(fr, dfr);

    [fr0, ~, ~, fwidth] = dfr_gaussfit(fr, dfr);
    [f0min, dfrmin, minfit] = dfrvalleyfit(fr, dfr, fr0);
    [f0max, dfrmax, maxfit] = dfrmaxfit(fr, dfr);

    data.fr     = fr;
    data.dfr    = dfr;
    data.f0min  = f0min;
    data.dfrmin = dfrmin;
    data.f0max  = f0max;
    data.dfrmax = dfrmax;
    data.maxfit = maxfit;
    data.minfit = minfit;
    data.fwidth = fwidth;

    switch option; case 'plot'
        %% plotting
        fminrange = f0min-200e6:1e6:f0min+200e6;
        fmaxrange = f0max-800e6:1e6:f0max+800e6;    
        fmindfrfit = feval(minfit, fminrange);
        fmaxdfrfit = feval(maxfit, fmaxrange);

        plot(fr, dfr, 'LineWidth', 2); 

        axis manual;
        plot(fminrange, fmindfrfit, 'k-', 'LineWidth', 1);
        plot(fmaxrange, fmaxdfrfit, 'k-', 'LineWidth', 1);
        plot([f0min f0max], [dfrmin dfrmax], 'o');
        fillmarkers;
    end
end

function [fr0, dfr0, dfrmin, fwidth] = dfr_gaussfit(fr, dfr)
    
    args = {'MinPeakProminence', 10e3, 'MinPeakWidth', 100e6};
    [~, fr0, w, p] = findpeaks(-dfr, fr, args{:});
    dfr0 = median(dfr, 'omitnan');
    
    function y = gauss(x, x0, y0, h, w)
        c = w/(2*sqrt(2*log(2)));
        y = y0 - h*exp( -(x-x0).^2 / (2*c^2) );
    end
    
    initial = [fr0, dfr0, p, w];
    model = @(c,x)gauss(x, c(1), c(2), c(3), c(4));
    final = nlinfit(fr, dfr, model, initial);

    dfrfit = model(final, fr);
%     plot(fr, dfrfit);
    
    fr0   = final(1);
    dfr0  = final(2);
    dfrmin = dfr0 - final(3);
    fwidth = final(4);
end


function [f0, dfr0, fitresult] = dfrvalleyfit(fr, dfr, fr0)
    width = 0.1e9;
    f1 = fr0 - width/2;
    f2 = fr0 + width/2;
    i = (fr > f1 & fr < f2);
    %  p1*x^2 + p2*x + p3
    fitresult = fit(fr(i)', dfr(i)', 'poly2');
    a = fitresult.p1;
    b = fitresult.p2;
    
    f0 = -b/(2*a);
    dfr0 = feval(fitresult, f0);
end

function [f0, dfr0, fitresult] = dfrmaxfit(fr, dfr)
    width = 0.2e9;
    [~, i] = min(dfr);
    [~, i] = max(dfr(1:i));
    f1 = fr(i) - width/2;
    f2 = fr(i) + width/2;
    i = (fr > f1 & fr < f2);
    %  p1*x^2 + p2*x + p3
    fitresult = fit(fr(i)', dfr(i)', 'poly2');
    a = fitresult.p1;
    b = fitresult.p2;
    
    f0 = -b/(2*a);
    dfr0 = feval(fitresult, f0);
end

function [fr, dfr] = rmoutliers(fr, dfr)

%     [Mfr1 Mfr2] = meshgrid(fr, fr);
%     Mdfr = abs(Mfr1 - Mfr2);
%     I = abs(Mdfr - median(diff(fr))) < 0.1e6;
%     i = sum(tril(I),1) ~= 1;
%     fr(i) = NaN;

    diff_fr = diff(fr);
    mdfr = median(diff_fr);
    i = abs(diff_fr - mdfr) < 0.1e6;
    i = [0 ~i];    

    dfr0 = median(dfr);
    dfrmin = dfr0 - 0.1e6;
    dfrmax = dfr0 + 0.05e6;
    j = dfr > dfrmax | dfr < dfrmin;
    dfr(j) = NaN;
    
    dfr = hampel(dfr,10);
    
    fr(i|j) = [];
    dfr(i|j) = [];
end