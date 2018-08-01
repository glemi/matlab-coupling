% function HBAR_Mplot(f, M, codes)
function h = HBAR_Mplot(f, M, codes, varargin)   
    prepareaxes(codes);
    
    n = length(codes);
    for k = 1:n
        code = codes{k};
        hax = getaxes(code);
        axes(hax); %#ok
        
        h(k) = HBAR_curveplot(f, M(:,k), code, varargin{:});
    end
end

function prepareaxes(codes)
    iabs = strncmp(codes, 'env:abs', 7);
    iph  = strncmp(codes, 'env:ph', 6);
    ire  = strncmp(codes, 'env:re', 6);
    iim  = strncmp(codes, 'env:im', 6);
    ifa  = strncmp(codes, 'fenv:abs', 8);
    ifi  = strncmp(codes, 'fenv:im', 7);
    ienv = [iabs; iph; ire; iim; ifa; ifi];
    nenv = sum(any(ienv,2));
    
    iripple = strncmp(codes, 'ripples:', 8);
    nripple = sum(iripple);
    
    irenv = strncmp(codes, 'renv:', 5);
    nrenv = sum(irenv);
    
    idirect = strncmp(codes, 'direct:', 7);
    ndirect = sum(idirect);
    
    ntotal = nenv+nripple+nrenv+ndirect;
    
    ncols = floor(sqrt(ntotal));
    nrows = ceil(ntotal/ncols);
    
    for k = 1:ntotal
        hax = subplot(nrows, ncols, k);
        
        if k > ntotal - ncols
            xlabel 'f [GHz]';
        end
    end
end

function haxes = getaxes(code)
    hfig = gcf;
    
    i = strfind([code ':'], ':')-1;
    tag = code(1:i(2));
    haxes = findobj(hfig, 'type',  'Axes', 'Tag', tag);
    
    if isempty(haxes)
        haxes = findobj(hfig, 'type',  'Axes', 'Tag', '');
        haxes = haxes(end);
        haxes.Tag = tag;
    end
end