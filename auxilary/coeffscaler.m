function final = coeffscaler(algo, model, initial)
    
    function y = modelwrapper(coeffs, x)
        coeffs = unscale(initial, coeffs);
        y = model(coeffs, x);
    end        

    final = algo(@modelwrapper, scale(initial));
    final = unscale(initial, final);
end

function coeffs = scale(coeffs)
    scalefact = 10.^(-floor(log10(abs(coeffs))));
    scalefact(coeffs == 0) = 1;
    coeffs = scalefact.*coeffs;
end

function coeffs = unscale(initial, coeffs)
    scalefact = 10.^(floor(log10(abs(initial))));
    scalefact(initial == 0) = 1;
    coeffs = scalefact.*abs(coeffs).*sign(initial);
end