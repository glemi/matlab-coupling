% function [fr, fa, ir, ia] = resonance_find(f, z, n, df)
%
%   f  : frequency vector
%   z  : impedance vector
%   n  : number of resonances to find
%   df : minimum spacing between resonance peaks
%
% usage:    
%           [fr, fa, ir, ia] = resonance_find(f, z, n, df);
%
function [fr, fa, ir, ia] = resonance_find(f, z, n, df)
    opt n double 1;
    opt df double 4e6;
    
    ir = NaN; ia = NaN;
    fr = NaN; fa = NaN;

    iSpacing = round(df/median(diff(f(1:10))));
    iSpacing = min(iSpacing, length(f)-2);
    %args = {'MinPeakProminence', .2, 'MinPeakWidth', 100e6/df};
    %args = {'MinPeakProminence', .09, 'NPeaks', 1, 'SortStr', 'descend'};
    args = {'SortStr', 'descend', 'MinPeakDistance', iSpacing};
    %args = {};
    
    logz = log(z);
    logy = log(1./z);
    
    logy = (logy - min(logy))/(max(logy) - min(logy));
    logz = (logz - min(logz))/(max(logz) - min(logz));
    
    [~, ia, ~, pa] = findpeaks(logz, args{:});
    [~, ir, ~, pr] = findpeaks(logy, args{:});
    
    % sort by peak prominence and keep only the n most prominent.
    [~, aOrder] = sort(pa, 'descend');
    [~, rOrder] = sort(pr, 'descend');
    ia = ia(aOrder);
    ir = ir(rOrder);
    
    % now sort peak indices by occurrence (i.e. ascending with frequency)
    [ir, ia] = align(ir, ia, n);
    
    fa = f(ia);
    fr = f(ir);
    
    %fprintf(' Resonance:      f_r = %.2f GHz\n', fr/1e9);
    %fprintf(' Anti-Resonance: f_a = %.2f GHz\n', fa/1e9);
    
    
%     fig resfind:eins;clf;
%     plot(f, logz);
%     plot(f, logy);
%     plot(fr, logy(ir), 'x');
%     plot(fa, logz(ia), 'x');
%     xscale log;
%     yscale lin;
end


function [x, y] = align(x, y, n)
    [X, Y] = meshgrid(x, y);
    D = Y-X;
    D(D<0) = inf;
    [d, i] = min(D, [], 1);
    y = y(i);
    
    x = x(d~=inf);
    y = y(d~=inf);
    
    n = min(n, length(x));
    x = x(1:n);
    y = y(1:n);
end
