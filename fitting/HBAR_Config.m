classdef HBAR_Config < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        config_file = 'HBAR_config.txt';
        config;
    end
    
    methods
        function this = HBAR_Config(argument)
            if nargin == 0
                this.load(this.config_file);
            elseif ischar(argument)
                this.load(argument);
                this.config_file = argument;
            elseif isstruct(argument) 
                this.config = argument;
                this.config_file = '';
            end
        end
        function set_parameters(this, params)
            this.config = HBAR_parameters(this.config, params);
        end
        function values = get_parameters(this, names)
            values = HBAR_parameters(this.config, names);
        end
        function load(this, filename)
            this.config = HBAR_loadconfig(filename);
        end
        
        function process(this)
            n = length(this.config.layers);
            kt2 = this.config.override.kt2; % may be empty
            for k = 1:n 
                matdata = this.config.materials(k);
                matdata = procmat(matdata, kt2);
                this.config.materials(k) = matdata;
            end
        end
    end
end

function matdata = procmat(matdata, kt2)
    epsr = matdata.Permittivity;
    c    = matdata.Stiffness;
    e    = matdata.PiezoCst;
    nu   = matdata.Poisson;
    Q    = matdata.Qfactor;
    rho  = matdata.Density;

    % this has no effect if nu = 0;
    c33E = c*(1-nu)/((1+nu)*(1-2*nu)); 
    eps0 =  8.854187817e-12;
    epsf = epsr*eps0;
    
    if e > 0 && ~isempty(kt2)
        c33D = c33E/(1-kt2);
        epsS = epsf*(1-kt2);
        e    = sqrt(c33D*epsS*kt2);
    else
        epsS = epsf - e^2/c33E;
        c33D = c33E + e^2/epsS; 
        kt2  = e^2/(c33D*epsS);
    end
    
    c33D = c33D*exp(1i/(2*Q)); % losses 
    v    = sqrt(c33D/rho);
    
    matdata.Permittivity = epsS;
    matdata.Stiffness = c33D;
    matdata.PiezoCst = e;
    matdata.Velocity = v;
    matdata.Coupling = kt2;
end
