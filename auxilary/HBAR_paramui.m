classdef HBAR_paramui < handle
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        params = {'ePiezo' 'cPiezo' 'epsPiezo' 'tPiezo' 'tTopEl' 'tBotEl' 'tSubst' 'QSi' 'cSi' 'xC0'};
        values 
        config
        figure
        frequency
        Callback
    end
    
    properties (Access = private)
        scaledvalues
        scalefactors
        spinners = matlab.ui.container.internal.JavaWrapper.empty;
        block = true;
    end
    
    methods (Access = public)
        function this = HBAR_paramui(config)
            if nargin > 0 
                if isstruct(config)
                    this.config = config;
                elseif ischar(config) && logical(exist(config, 'file'))
                    this.config = HBAR_loadconfig(config);
                else
                    error 'Invalid config parameter!';
                end
            else
                this.config = HBAR_loadconfig('HBAR_config.txt');
            end
        end
        
        function start(this)
            this.values = HBAR_parameters(this.config, this.params);
            this.figure = fig('hbarui');
            this.uisetup();
            this.uiupdate();
            this.block = false;
        end
    end
    
    methods(Access = private)
        function spinnercallback(this, hs, event)
            n = length(this.params);
            for k = 1:n
                spinner = this.spinners(k);
                sv(k) = spinner.Value;
            end
            this.scaledvalues = sv;
            this.values = sv.*this.scalefactors;
            pairs = nvpairs(this.params, this.values);
            this.config = HBAR_parameters(this.config, pairs);
            if isa(this.Callback, 'function_handle') && ~isempty(this.Callback) &&  ~this.block
                this.Callback(this.config);
            end
        end
        function uiupdate(this)
            this.values = HBAR_parameters(this.config, this.params);
            [sv, sf] = scale(this.values);
            this.scaledvalues = sv;
            this.scalefactors = sf;
            
            n = length(this.params);
            for k = 1:n
                spinner = this.spinners(k);
                spinner.Value = round(sv(k), 3, 'significant');
            end
        end
        function uisetup(this)
            figure(this.figure);clf;   %#ok;     
            this.figure.Position = [1647 524 257 36*length(this.params)];
            this.figure.ToolBar = 'none';
            this.figure.MenuBar = 'none';
            
            hbox = uix.HBox('Parent', this.figure);
            vbox = uix.VBox('Parent', hbox);
            
            hbox.Padding = 10;
            vbox.Spacing = 10;
            
            n = length(this.params);
            for k = 1:n
                hbox = uix.HBox('Parent', vbox);
                
                hl = uicontrol(hbox, 'Style', 'text');
                hl.String = this.params{k};

                hs = uicomponent(hbox, 'Style', 'jspinner');             
                setupSpinner(hs, scale(this.values(k)));
                hs.StateChangedCallback = @(~,~)this.spinnercallback();
                this.spinners(k) = hs;
            end
            drawnow;
        end
    end
    
end

function setupSpinner(hs, value)
    value = round(value, 3, 'significant');
    vstep = round(value*0.02, 1, 'significant');
    vmin = round(value/10, 1, 'significant');
    vmax = round(value*1000, 1, 'significant');
    
    model = javax.swing.SpinnerNumberModel(value, vmin, vmax, vstep);
    hs.JavaComponent.setModel(model);
end
    
% getexponent = @(y)floor(log10(mean(y))) - mod(floor(log10(mean(y))), 3);
function [sv, sf] = scale(values)
    x = values;
    X = floor(log10(x));
    e = X - mod(X,3);
    sf = 10.^e;
    sv = x./sf;
end