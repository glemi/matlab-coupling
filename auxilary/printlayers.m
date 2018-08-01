function printlayers(layers)
    
    %header  = ' # type name material  thick    c33D     e33         v         eps        kt2\n';
    %pattern = ' %d: %2s %6s  %5s   %5s  %5s  %8s  %9s  %10s  %6s\n';
    
    header  = ' # material  thick    c33D     e33         v         eps        kt2\n';
    pattern = ' %d: %5s   %5s  %5s  %8s  %9s  %10s  %6s\n';
    
    fprintf(header);
    
    n = length(layers);
    for k = 1:n
        layer = layers(k);
        material = layer.material;
        
        args{1} = k;
        %args{2} = layer.type;
        %args{3} = layer.name;
        args{4} = material.Symbol;
        args{5} = siPrefix(layer.thickness, 'm');
        args{6} = siPrefix(real(material.Stiffness), 'Pa');
        args{7} = siPrefix(real(material.PiezoCst), 'C/m²');
        args{8} = siPrefix(real(material.Velocity), 'm/s', '%.1f');
        args{9} = siPrefix(real(material.Permittivity), 'F/m', '%.1f');
        args{10} = sprintf('%.1f%%', real(material.Coupling*100));
        
        args(2:3) = [];
        
        if layer.material.PiezoCst ~= 0
            args{5} = siPrefix(real(material.PiezoCst), 'C/m²');
            args{8} = sprintf('%.1f%%', real(material.Coupling*100));
        else
            args{5} = '-    ';
            args{8} = '-   ';
        end
        
        fprintf(pattern, args{:});
        
    end
    
    
    
end