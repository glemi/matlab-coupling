% function HBAR_Mplot(f, M, codes)
function h = HBAR_curveplot(f, Y, code, varargin)   
    fGHz = f/1e9;

    c = getplotconfig(code);
    title(c.title);
    xscale(c.xscale);
    yscale(c.yscale);

    y = Y/c.scale;
    h = plot(fGHz, y, varargin{:}); 
    xlim(round(2*[min(fGHz) max(fGHz)])/2);
end

function plotconfig = getplotconfig(code)
	switch code
        case 'direct:re',         title='Re\{Z\} [\Omega]';   scale=1;      yscale='lin'; xscale='lin';
        case 'direct:im',         title='Im\{Z\} [\Omega]';   scale=1;      yscale='lin'; xscale='lin';
        case 'direct:ph',         title='\phi\{Z\} [\circ]';  scale=pi/180; yscale='lin'; xscale='lin';
        case 'direct:abs',        title='|Z| [\Omega]';       scale=1;      yscale='log'; xscale='log';
            
		case 'env:re:avg',        title='Re\{Z\} [\Omega]';   scale=1;      yscale='lin'; xscale='lin';
		case 'env:re:upper',      title='Re\{Z\} [\Omega]';   scale=1;      yscale='lin'; xscale='lin';
		case 'env:re:lower',      title='Re\{Z\} [\Omega]';   scale=1;      yscale='lin'; xscale='lin';
		case 'env:im:avg',        title='Im\{Z\} [\Omega]';   scale=1;      yscale='lin'; xscale='lin';
		case 'env:im:upper',      title='Im\{Z\} [\Omega]';   scale=1;      yscale='lin'; xscale='lin';
		case 'env:im:lower',      title='Im\{Z\} [\Omega]';   scale=1;      yscale='lin'; xscale='lin';
		case 'env:abs:avg',       title='|Z| [\Omega]';       scale=1;      yscale='log'; xscale='log';
		case 'env:abs:upper',     title='|Z| [\Omega]';       scale=1;      yscale='log'; xscale='log';
		case 'env:abs:lower',     title='|Z| [\Omega]';       scale=1;      yscale='log'; xscale='log';
		case 'env:ph:avg',        title='\phi\{Z\} [\circ]';  scale=pi/180; yscale='lin'; xscale='lin';
		case 'env:ph:upper',      title='\phi\{Z\} [\circ]';  scale=pi/180; yscale='lin'; xscale='lin';
		case 'env:ph:lower',      title='\phi\{Z\} [\circ]';  scale=pi/180; yscale='lin'; xscale='lin';
            
        case 'fenv:im:avg',       title='Im\{Z\}/f [\Omega/GHz]';   scale=1e-9;  yscale='lin'; xscale='lin';
        case 'fenv:im:upper',     title='Im\{Z\}/f [\Omega/GHz]';   scale=1e-9;  yscale='lin'; xscale='lin';
        case 'fenv:im:lower',     title='Im\{Z\}/f [\Omega/GHz]';   scale=1e-9;  yscale='lin'; xscale='lin';
        case 'fenv:abs:avg',      title='|Z|\cdot{}f [\Omega GHz]'; scale=1e9;   yscale='lin'; xscale='log';
        case 'fenv:abs:upper',    title='|Z|\cdot{}f [\Omega GHz]'; scale=1e9;   yscale='lin'; xscale='log';
        case 'fenv:abs:lower',    title='|Z|\cdot{}f [\Omega GHz]'; scale=1e9;   yscale='lin'; xscale='log';    
            
        case 'renv:abs:u/l',       title='|Z_{up}|/|Z_{lo}|';  scale=1;      yscale='lin'; xscale='lin';
            
		case 'ripples:C0',        title='C_0 [pF]';           scale=1e-12;  yscale='lin'; xscale='lin';
		case 'ripples:Cm',        title='C_m [10^{-15}]';     scale=1e-15;  yscale='lin'; xscale='lin';
		case 'ripples:Lm',        title='L_m [10^{-3}]';      scale=1e-3;   yscale='lin'; xscale='lin';
		case 'ripples:Qm',        title='Q_m';                scale=1;      yscale='lin'; xscale='lin';
		case 'ripples:Rm',        title='R_m [10^3]';         scale=1e3;    yscale='lin'; xscale='lin';
		case 'ripples:keff',      title='k_{eff}^2 [%]';      scale=1e-2;   yscale='lin'; xscale='lin';
		case 'ripples:kt2',       title='k_t^2 [%]';          scale=1e-2;   yscale='lin'; xscale='lin';
		case 'ripples:fr',        title='f_r [GHz]';          scale=1e9;    yscale='lin'; xscale='lin';
		case 'ripples:fa',        title='f_a [GHz]';          scale=1e9;    yscale='lin'; xscale='lin';
		case 'ripples:dfr',       title='\Delta{f_r} [MHz]';  scale=1e6;    yscale='lin'; xscale='lin';
		case 'ripples:diff(fr)',  title='\Delta{f_r} [MHz]';  scale=1e6;    yscale='lin'; xscale='lin';
	end
    plotconfig = var2struct(title, scale, yscale, xscale);
end