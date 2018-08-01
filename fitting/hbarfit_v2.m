function [zfit, data] = hbarfit_v2(f, z, varargin)
    
%     cnames =  {'C33',   'E33', 'tPiezo', 'Qpiezo', 'tanDelta'};
%     initial = [276.845e9 1.9    500e-9   10         0.002];
%     
%     cnames =  {'C33',   'E33', 'tPiezo', 'Qpiezo'};
%     initial = [276.845e9 1.9    500e-9   1000    ];
%      
%     cnames =  {'cAlN' 'eAlN'  'epsAlN'   'tPiezo' 'tBotEl' 'tTopEl' };
%     cunits =  {'P'    'C/m^2' '1'        'm'      'm'      'm'      }; 
%     initial = [380e9   1.5     12       880e-9 ,  300e-9 , 200e-9   ];
    
%     cnames =  { 'eAlN'  'epsAlN'   'tPiezo' 'tBotEl' 'tTopEl' };
%     cunits =  { 'C/m^2' '1'        'm'      'm'      'm'      }; 
%     initial = [ 1.5     12         880e-9   300e-9   200e-9  ];
    
%     cnames =  {'cAlN',   'eAlN',  'epsAlN', 'tPiezo'};
%     cunits =  {'P',      'C/m^2', '1',      'm'     }; 
%     initial = [233e9,     1.4     13,       880e-9 ];
    
    cnames =  {'cAlN',   'eAlN',  'epsAlN'};
    cunits =  {'P',      'C/m^2', '1',    }; 
    initial = [380e9,     1.4,     10,   ];
    
%     cnames =  {  'Cd',    'tBotEl' 'tTopEl'};
%     cunits =  {   'C/m^2'  'm'     'm'}; 
%     initial = [   1e-12    300e-9  100e-9];
        
%     cnames =  {'E33', 'tPiezo', 'Qpiezo', 'tanDelta'};
%     initial = [1.9    500e-9   10         0.002];
    
%    params = { 'tSubst', 725e-6, 'tOxide', 500e-9, 'tPiezo', 920e-9};
    params = { 'Cd', 352e-15, 'tPiezo', 920e-9, 'tSubst', 725e-6, ...
         'tOxide', 500e-9, 'tBotEl', 306e-9, 'tTopEl', 175e-9};
%     params = {'tSubst', 725e-6, 'tOxide', 500e-9, 'tBotEl', 300e-9, ...
%         'tTopEl', 200e-9};
    
    
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
        
        [z, data] = HBAR_v2(f, cpairs{:}, params{:});
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

