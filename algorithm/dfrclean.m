function dfr_clean = dfrclean(f, dfr)

    dfrdat = dfr_analyze(f, dfr);

%     fr  = dfrdat.fr;
%     dfr = dfrdat.dfr;
    fr = f;
    
    ipeak = fr > dfrdat.f0min - dfrdat.fwidth & fr < dfrdat.f0min + dfrdat.fwidth;
    
    
    n = length(fr);
    dfr_peak = smooth(dfr, 20, 'lowess');
    %dfr_rim = movmedian(dfr', n/5, 'omitnan');
    
    
    dfr_rim  = smooth(hampel(dfr, 50), n/5, 'lowess');
    ipeak = smooth(ipeak, n/5);
    
    
    
    
    dfr_clean = dfr_peak.*ipeak + dfr_rim.*(1-ipeak);
    
    dfr_clean = interp1(fr, dfr_clean, f, 'linear', 'extrap');
    
%     dfr_clean = movmedian(dfr_clean, 10, 'omitnan');
%     dfr_clean = smooth(dfr_clean, 10, 'lowess');
%     
%     dfr_clean = smooth(dfr, 20, 'sgolay');
    
    % dfr_clean = dfr;

%     fig dfrclean:eins;clf;
%     
%     plot(f, dfr);
%     plot(f, dfr_clean);
%     
%     ylim([5.74 5.88]*1e6);
    
end