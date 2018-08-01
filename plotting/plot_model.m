function h = plot_model(f, Zel, dname)
   
   h = plot(f/1e9, abs(Zel));
   
   ax = gca;
   xlabel 'frequency [GHz]';
   
   
   ax.YScale = 'log';
   ylabel 'Z_{el} [\Omega]';
   
   if nargin > 2
       h.DisplayName = dname;
   end
   
end