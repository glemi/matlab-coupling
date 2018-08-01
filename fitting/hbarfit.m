function [zfit, data] = hbarfit(f, z, varargin)
    
    cnames =  {'C33',   'E33', 'tPiezo', 'Qpiezo', 'tanDelta'};
    initial = [276.845e9 1.9    500e-9   10         0.002];
    
    cnames =  {'C33',   'E33', 'tPiezo', 'Qpiezo'};
    initial = [276.845e9 1.9    500e-9   1000    ];
    
    cnames =  {'C33',   'E33', 'tPiezo', 'tBotE', 'tTopE'};
    cunits =  {'P',   'C/m^2', 'm'     , 'm'    , 'm'}; 
    initial = [276.845e9 1.9    900e-9 , 100e-9 , 80e-9];
    
        
%     cnames =  {'E33', 'tPiezo', 'Qpiezo', 'tanDelta'};
%     initial = [1.9    500e-9   10         0.002];
    
    params = {'tSi', 725e-6};
    
    
    sinitial = scale(initial, initial);
    
    model = mkmodel(cnames, params, initial);

    [sfinal, res, jac] = nlinfit(f, z, model, sinitial);
    sfinal = abs(sfinal);
    sj = sum(jac,1)/length(jac);
    
    [zfit, hbardata] = model(sfinal, f);
    
    data.names   = cnames;
    data.units   = cunits;
    data.initial = initial;
    data.final   = unscale(initial, sfinal);
    data.hbar    = hbardata;
    data.model   = @(f)model(final, f);
end

function model = mkmodel(cnames, params, initial)
    
    function [z, data] = hbarmodel(coeffs, f)
        disp(coeffs);
        coeffs = unscale(initial, coeffs);
        values = num2cell(abs(coeffs));
        cpairs = [cnames; values];
        
        [z, ~, data] = HBAR(f, cpairs{:}, params{:});
        z = deripple(f, z);
        %z = abs(z);
        plot(f, z, 'k-');drawnow;
    end
    
    model = @hbarmodel;
end

function scaled = scale(initial, coeffs)
    scalefact = 10.^(-floor(log10(initial)));
    scaled = scalefact.*coeffs;
end

function unscaled = unscale(initial, scaled)
    scalefact = 10.^(floor(log10(initial)));
    unscaled = scalefact.*scaled;
end

