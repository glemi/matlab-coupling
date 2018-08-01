function lstruct = HBAR_layer(layid, matid, thickness, area, parameters)
    material = HBAR_materials(matid, parameters);
    override(layid, parameters);
    
    d = material.Density;
    v = material.Velocity;
    
    Zac     = area*d*v;
    phase   = @(f) 2*pi*f*thickness/v;
    TMatrix = @tmatrix;
    
    function T = tmatrix(f)
        ph = phase(f);
        ph = reshape(squeeze(ph), [1 1 numel(f)]);
        T  = [cos(ph), 1i*Zac*sin(ph); 1i/Zac*sin(ph), cos(ph)];
    end

    lstruct = var2struct(Zac, phase, TMatrix, thickness, area, material);
    if nargout == 0
        assignin('caller', layid, lstruct);
    end
end
