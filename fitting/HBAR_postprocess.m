% function [M, F] = HBAR_postprocess(f, Z, codes, N, nskip)
%  
% Uses HBAR measured or model output data to obtain various derived
% property spectra. The data is arranged in the N x m matrix M, where m is
% the number of elements in the thrid argument 'codes'.  
%
% HBAR_postprocess resamples the data so that all resulting spectra are of
% equal length N. The Nx1 column vector F is a downsampled version of f.
% Usually the sample rate of Z needs to be 10 to 100 times higher to avoid
% aliasing of the HBAR ripples, but the derived spectra are vastly
% oversampled at this rate. A few hundred points should normally suffice. 
%
% The codes specify which spectra are computed and returned in M. 
%
% usage: 
%           codes = {'env:abs' 'env:ph' 'ripples:kt2'};
%           [M, F] = HBAR_postprocess(f, z, codes, 200); 
%
%           % Uses length(f) as the size for all column-data in M (not
%           % recommended!)
%           M = HBAR_postprocess(f, z, codes);
%
% input parameters: 
%           f      : frequency vector
%           z      : complex-valued HBAR input impedance
%           codes  : cell-array of strings of spectra to compute
%           N      : number of samples to compute for each spectrum 
%                    (optional, defaults to length(f))
%           noSkip : prevents skipping of ripples in the ripplefit
%                     algorithm
%  
% output parameters:
%           M      : Matrix of column vectors of computed spectra of size
%                    N x m, where m is the number of elements in codes
%           F      : Frequency vector of size N x 1 
%   
% list of valid codes:
% 	envelope spectra
%      env:re:avg       :  average of the envelope of Re{Z}
%      env:re:upper     :  upper part of the envelope of Re{Z}
%      env:re:lower     :  lower part of the envelope of Re{Z}
%      env:im:avg       :  average of the envelope of Im{Z}
%      env:im:upper     :  upper part of the envelope of Im{Z}
%      env:im:lower     :  lower part of the envelope of Im{Z}
%      env:abs:avg      :  average of the envelope of |Z|
%      env:abs:upper    :  upper part of the envelope of |Z|
%      env:abs:lower    :  lower part of the envelope of |Z|
%      env:ph:avg       :  average of the envelope of phase{Z}
%      env:ph:upper     :  upper part of the envelope of phase{Z}
%      env:ph:lower     :  lower part of the envelope of phase{Z}
%      
%	ripple spectra				  
%      ripples:C0       :  C0 ripple spectrum (static capacitance)
%      ripples:Cm       :  Cm ripple spectrum (mechanical capacitance/compliance)
%      ripples:Lm       :  Lm ripple spectrum (mechanical inductance/inertia)
%      ripples:Qm       :  Qm ripple spectrum (mechanical quality factor)
%      ripples:Rm       :  Rm ripple spectrum (mechanical resistance/damping)
%      ripples:keff     :  ripple effective coupling factor 
%      ripples:kt2      :  same as ripples:keff
%      ripples:fr       :  resonance frequencies
%      ripples:fa       :  antiresonance frequencies
%      ripples:dfr      :  resonance frequency spacing directly obtained from fit
%      ripples:diff(fr) :  resonance frequency spacing obtained from diff(fr)
%
% see also: 
%           HBAR_v3, bvdMultiFit, multiRippleFit, hbarmultifit
function [M, F, ripples] = HBAR_postprocess(f, Z, codes, N, varargin)
    opt N double length(f);

    % options
    args = varargin;
    noSkip = any(strcmp(args, 'noSkip'));
    noExtrap = any(strcmp(args, 'noExtrap'));
    
    F = linspace(f(1), f(end), N);
    f = f(:); F = F(:);
    
    rez = real(Z);
    imz = imag(Z);
    phz = angle(Z);
    abz = abs(Z);
    
    if any(strncmp(codes, 'direct:', 7))
        N = length(f);
        F = f;
    end 
    
    check = @(c, p)(any(cell2mat(strfind(c, p))));
    
    if check(codes, 'env:abs')
        [aabz, uabz, labz] = deripple(f, abz, N);
    end 
    if check(codes, 'env:ph')
        [aphz, uphz, lphz] = deripple(f, phz, N);
    end 
    if check(codes, 'env:re')
        [arez, urez, lrez] = deripple(f, rez, N);
    end 
    if check(codes, 'env:im')
        [aimz, uimz, limz] = deripple(f, imz, N);
    end 
    
    if any(strncmp(codes, 'ripples:', 8))
        if noSkip
            ripples = HBAR_ripplefit(f, abz);
        else
            ripples = HBAR_ripplefit(f, abz, N);
        end
        %fig toc; plot([ripples.toc]); 
        %title(sprintf('%f', sum([ripples.toc])));
    else
        ripples = struct([]);
    end
    
    rget = @(field)getRipples(ripples, field, F, noExtrap);
    
    m = length(codes);
    M = nan(N, m);
    for k = 1:m
        switch codes{k}
            case 'direct:re',           M(:,k) = rez;
            case 'direct:im',           M(:,k) = imz;
            case 'direct:abs',          M(:,k) = abz;
            case 'direct:ph',           M(:,k) = phz;
                                        
            case 'env:re:avg',          M(:,k) = arez;    
            case 'env:re:upper',        M(:,k) = urez;
            case 'env:re:lower',        M(:,k) = lrez;
            case 'env:im:avg',          M(:,k) = aimz;    
            case 'env:im:upper',        M(:,k) = uimz;
            case 'env:im:lower',        M(:,k) = limz;
            case 'env:abs:avg',         M(:,k) = aabz;    
            case 'env:abs:upper',       M(:,k) = uabz;
            case 'env:abs:lower',       M(:,k) = labz;
            case 'env:ph:avg',          M(:,k) = aphz;    
            case 'env:ph:upper',        M(:,k) = uphz;
            case 'env:ph:lower',        M(:,k) = lphz;      
                
            case 'fenv:im:avg',          M(:,k) = aimz./F;    
            case 'fenv:im:upper',        M(:,k) = uimz./F;
            case 'fenv:im:lower',        M(:,k) = limz./F;
            case 'fenv:abs:avg',         M(:,k) = aabz.*F;    
            case 'fenv:abs:upper',       M(:,k) = uabz.*F;
            case 'fenv:abs:lower',       M(:,k) = labz.*F;
                
            case 'renv:abs:u/l',        M(:,k) = uabz./labz;
                                     
			case 'ripples:C0',          M(:,k) = rget('C0');
            case 'ripples:Cm',          M(:,k) = rget('Cm');
			case 'ripples:Lm',          M(:,k) = rget('Lm');
			case 'ripples:Qm',          M(:,k) = rget('Qm');
			case 'ripples:Rm',          M(:,k) = rget('Rm');
			case 'ripples:keff',        M(:,k) = rget('kt2');
			case 'ripples:kt2',         M(:,k) = rget('kt2');
			case 'ripples:fr',          M(:,k) = rget('fr');
			case 'ripples:fa',          M(:,k) = rget('fa');
			case 'ripples:dfr',         M(:,k) = rget('dfr');
			case 'ripples:diff(fr)',    M(:,k) = rget('diff(fr)');
            case 'ripples:smooth(dfr)', M(:,k) = rget('smooth(dfr)');
                
			case 'ripples:1/Lm',        M(:,k) = 1./rget('Lm');  
			case 'ripples:1/Rm',        M(:,k) = 1./rget('Rm');
                
            otherwise
                error 'invalid hbar post-process code';
        end
    end
end


function Y = getRipples(ripples, field, F, noExtrap)
    f = [ripples.fr];
    v = [ripples.valid];
    
    if isfield(ripples, field)
        y = [ripples.(field)];
    elseif strncmp(field, 'diff(fr)', 8)
        y = diff(f)/round((f(2)-f(1))/ripples(1).dfr);
        f = f(2:end);
        v = v(1:end-1) & v(2:end);
    elseif strncmp(field, 'smooth(dfr)', 11)
        y = [ripples.dfr];
        y = dfrclean(f, y);
    else
        error 'unknown ripple data property field';
    end
    
    i = isnan(y) | isinf(y) | ~v;
    if noExtrap
        Y = interp1(f(~i), y(~i), F, 'linear');
    else
        Y = interp1(f(~i), y(~i), F, 'linear', 'extrap');
    end
    Y = hampel(Y,10); %eliminate outliers produced by extrapolation
end

