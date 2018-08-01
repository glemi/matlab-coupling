function fit = bvdextfit_v2(f, Z, select, varargin)
    
	names  = {'fr' 'Rm' 'Qm' 'C0' 'td' 'R0' 'Rp' 'Rc' 'Lc' 'Cs' 'Cp' 'X1' 'X2' 'Y1' 'Y2'};
	always = {'fr' 'Rm' 'Qm' 'C0'};
	select = cellstr(select); select = select(:);
	
    select = union(select, always);
    select = templateorder(select, names);
	iselect = any(strmatcmp(names, select), 1);
	
	ivarble   = iselect;
    istatic   = ~ivarble;
	default   = getDefaultCoeffs(f, Z, names);
    initial   = getInitialCoeffs(f, Z, select);
    
    model     = mkmodel(default, ivarble);
    options   = statset('nlinfit');
    
    M = Z2M(Z); %progressPlot(f, M, true);
    nlmf = @(m,i)nlinmatrixfit(f, M, m, i, options, varargin{:});
    final = coeffscaler(nlmf, model, initial);
    
    final(ivarble) = final;
    final(istatic) = default(istatic);
    
    nvpairs = [names; num2cell(final)];
    fit = struct(nvpairs{:});
    fit.select = select;
    fit.Zfit   = M2Z(model(final(ivarble), f));
    fit.ZfitModel = @(f)M2Z(model(final(ivarble), f));
    
    %fit.Lc = -fit.Lc;
    %fit.X1 = -fit.X1;
    %fit.X2 = -fit.X2;
    fit.Cm = 1/(2*pi*fit.fr*fit.Qm*fit.Rm);
    fit.Lm = 1/(2*pi*fit.fr)*fit.Qm*fit.Rm;
    fit.fa = sqrt(1/(fit.Lm*fit.Cm)*(1 + fit.Cm/fit.C0))/(2*pi);
    
end

function Z = bvdmodel(f, coeffs)
    fr = coeffs(1);
    Rm = coeffs(2);
    Qm = coeffs(3);
    C0 = coeffs(4);
    td = coeffs(5);
    R0 = coeffs(6);
    Rp = coeffs(7);
    Rc = coeffs(8);
    Lc = coeffs(9);
    Cs = coeffs(10);
    Cp = coeffs(11);
    X1 = coeffs(12);
    X2 = coeffs(13);
    Y1 = coeffs(14);
    Y2 = coeffs(15);
    
    Cm = 1/(2*pi*fr*Qm*Rm);
    Lm = (Qm*Rm)/(2*pi*fr);
    
    ZC0 = 1./(2i*pi*f*C0*exp(-1i*td)); 
    ZCm = 1./(2i*pi*f*Cm);
    ZLm =     2i*pi*f*Lm;
    ZLc =     2i*pi*f*Lc;
    ZCp = 1./(2i*pi*f*Cp);
    ZCs = 1./(2i*pi*f*Cs);
    
    par = @(x,y)1./(1./x + 1./y);
    
    Zm = ZLm + ZCm + Rm;
    Z0 = ZC0 + R0;
    Z = par(Zm, Z0);
    Z = par(Rp, Z) + par(Rc + ZLc, ZCp) + ZCs;
    Z = Z +  1i*X1*2*pi*f  + 1i*(X2*2*pi*f).^2;
    Z = Z + 1./(Y1*2*pi*f) - 1./(Y2*2*pi*f).^2;
end

function default = getDefaultCoeffs(f, Z, names)
    fr = findMainResoance(f, Z);
    Rm = 20;
    Qm = 4;
    C0 = abs(mean(1./(2*pi*f(end).*Z(end))));
    
    td = 0;
    R0 = 0;
    Rp = inf; 
    Rc = 0;
    Lc = 0;
    Cp = 0;
    Cs = inf;
    X1 = 0;
    X2 = 0;
    Y1 = inf;
    Y2 = inf;

    default = nvassign(names);
end

function initial = getInitialCoeffs(f, Z, names)
    fr = findMainResoance(f, Z);
    Rm = 20;
    Qm = 4;
    C0 = abs(mean(1./(2*pi*f(end).*Z(end))));
    
    td =  0.04;
    R0 =  1;
    Rp =  1e6; 
    Rc =  2;
    Lc = -0.1e-12;
    Cp =  0.1e-12;
    Cs =  20e-12;
    X1 =  20e-12;
    X2 =  20e-12;
    Y1 = -20e-12;
    Y2 =  20e-12;

    initial = nvassign(names);
end

function model = mkmodel(default, i)
    model = @(c,f)wrapper(c,f);
    
    function M = wrapper(variable, f)
		coeffs    = default;
		coeffs(i) = variable;
        Z = bvdmodel(f, coeffs);
        M = Z2M(Z);
        
        %progressPlot(f, M, false);
    end
end

function fr = findMainResoance(f, Z)
    span = round(200e6/(diff(f(1:2))));
    P = smooth(angle(Z), span);
    args = {'MinPeakDistance', diff(f([1 end]))/2};
    [~, fr, ~, p] = findpeaks(P, f, args{:});
    [~, i] = max(p);
    fr = fr(i);
end

function out = nvassign(names, values)
    n = length(names);
    if nargin == 2 && nargout == 0
        for k = 1:n
            assignin('caller', names{k}, values(k));
        end
    elseif nargin == 1 && nargout == 1
        out = zeros(1,n);
        for k = 1:n
            out(k) = evalin('caller', names{k});
        end
    end
end

function M = Z2M(Z)
    Z = Z(:);
    M(:,1) = abs(Z);
    M(:,2) = angle(Z);
end

function Z = M2Z(M)
    [re, im] = pol2cart(M(:,2), M(:,1));
    Z = complex(re, im);
end


function progressPlot(f, M, keep)
    persistent h1 h2; try delete([h1 h2]); end; %#ok
    fig bvdextfit;
    if keep
        clf;
        subplot(2,1,1);      plot(f, M(:,1),  '-', 'LineWidth', 0.5);
        subplot(2,1,2);      plot(f, M(:,2),  '-', 'LineWidth', 0.5);
    else
        subplot(2,1,1); h1 = plot(f, M(:,1), 'k-', 'LineWidth', 0.5);
        subplot(2,1,2); h2 = plot(f, M(:,2), 'k-', 'LineWidth', 0.5);
    end
    drawnow;
end