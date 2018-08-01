function [zfit, data] = hbarfit_v2a(f, z, varargin)
    
    cnames =  {'cAlN' 'eAlN'  'epsAlN'   'tPiezo' 'tBotEl' 'tTopEl' };
    cunits =  {'P'    'C/m^2' '1'        'm'      'm'      'm'      }; 
    initial = [380e9   1.5     12       880e-9 ,  300e-9 , 200e-9   ];
    
%     cnames =  {'cAlN',   'eAlN',  'epsAlN'};
%     cunits =  {'P',      'C/m^2', '1',    }; 
%     initial = [380e9,     1.4,     10,   ];
    
    params = { 'tSubst', 725e-6, 'tOxide', 500e-9};
    %params = { 'Cd', 352e-15, 'tPiezo', 920e-9, 'tSubst', 725e-6, ...
    %     'tOxide', 500e-9, 'tBotEl', 306e-9, 'tTopEl', 175e-9};
    
    f = [f(:)];
    N = 200;
    [Mc, F] = z2M(f, z, N);  
    Mf = z2M(f, z);
    fig mplot; clf; Mplot(f, Mf, z);
    
    sinitial = scale(initial, initial);
    fine   = mkmodel(f, cnames, params, initial);
    coarse = mkmodel(f, cnames, params, initial, N);
    
    [sfinal, res, jac] = nlinmatrixfit(F, Mc, coarse, sinitial, 'MaxIter', 10);
    sfinal = abs(sfinal);
    sj = sum(jac,1)/length(jac);
    
    [M, hbardata] = fine(sfinal, f);
    
    zfit = M(:,1);
    
    data.names   = cnames;
    data.units   = cunits;
    data.initial = initial;
    data.final   = unscale(initial, sfinal);
    data.hbar    = hbardata;
    data.model   = @(f)coarse(final, f);
end

function model = mkmodel(f, cnames, params, initial, N)
    opt N double length(f);

    function [M, data] = hbarmodel(coeffs, ~)
        fprintf('%.4f ', coeffs); fprintf('\n');
        coeffs = unscale(initial, coeffs);
        values = num2cell(abs(coeffs));
        cpairs = [cnames; values];
        
        [Z, data] = HBAR_v2(f, cpairs{:}, params{:});
        [M, F] = z2M(f, Z, N);
        Mplot(F, M, [], 'k-');
        if any(isnan(M))
            disp(sum(sum(isnan(M))));
        end
    end
    
    model = @hbarmodel;
end

function [M, F] = z2M(f, Z, N)
    if nargin < 3
        N = [];
    end
    [a, au, al, F] = deripple(f, abs(Z), N);
    [p, pu, pl] = deripple(f, phase(Z), N);
    [r, ru, rl] = deripple(f, real(Z), N);
    [i, iu, il] = deripple(f, imag(Z), N);
    M = [a au al p pu pl r ru rl i iu il];
end

function scaled = scale(initial, coeffs)
    scalefact = 10.^(-floor(log10(initial)));
    scaled = scalefact.*coeffs;
end

function unscaled = unscale(initial, scaled)
    scalefact = 10.^(floor(log10(initial)));
    unscaled = scalefact.*scaled;
end

function Mplot(f, M, Z, ls)
    persistent p;
    if nargin < 4
        ls = '-';
    end
    
    fig mplot;
    set(gcf, 'DefaultAxesXScale', 'log');
    if nargin >= 3 && ~isempty(Z)
        subplot(2,2,1); title abs;   plot(f, [abs(Z)   M(:,1:3)]); yscale log;
        subplot(2,2,2); title phase; plot(f, [phase(Z) M(:,4:6)]);
        subplot(2,2,3); title real;  plot(f, [real(Z)  M(:,7:9)]);
        subplot(2,2,4); title imag;  plot(f, [imag(Z)  M(:,10:12)]);
    else 
        try delete(p); end;
        subplot(2,2,1); p(1,:) = plot(f, M(:,1:3), ls);
        subplot(2,2,2); p(2,:) = plot(f, M(:,4:6), ls);
        subplot(2,2,3); p(3,:) = plot(f, M(:,7:9), ls);
        subplot(2,2,4); p(4,:) = plot(f, M(:,10:12), ls);
    end
    drawnow;
end
