function data = cheeke_compute(params)
    
    % assign a local variable for each field of params;
    dissolve(params);
    
    % Mass fraction (density x thickness)
    mSubst = denSi*tSi;
    mTopEl = denPt*tTopE;
    mBotEl = denPt*tBotE;
    mPiezo = denPiezo*tPiezo;
    mOxide = denSiO*tSiO;
    mUpper = mPiezo + mTopEl + mBotEl;

    % Alternate names
    tSubst = tSi;
    tTopEl = tTopE;
    tBotEl = tBotE;
    denBotEl = denPt;
    denTopEl = denPt;
    denSubst = denSi;
    
    % Acoustic Time Constant (thickness / velocity) 
    TSubst = tSi/real(VSi);
    TTopEl = tTopE/real(VPt);
    TBotEl = tBotE/real(VPt);
    TPiezo = tPiezo/real(VPiezo);
    
    % Normalized Acoustic Impedances
    zSubst = ZSi/ZPiezo;
    zTopEl = ZPt/ZPiezo;
    zBotEl = ZPt/ZPiezo;
    
    delta_f0 = real(VSi)/(2*tSi);
    
    % delta_fN is dfr at the center of the first "normal region" (5)
    % gamma = pi/2: corresponds to the maximum of the SPRF
    delta_fN = delta_f0/(1 + mUpper/mSubst);
    
    % delta_fT is dfr at the center of the first "transition region" (6) 
    % where gamma = pi/2: corresponds to the minimum of the SPRF
    blah = TPiezo^2/mPiezo + TBotEl^2/mBotEl + mTopEl/(denPiezo^2*VPiezo^2);
    delta_fT = delta_f0/(1 + mSubst/TSubst*blah); 
    
    % "top resonator" and "composite substrate" resonant frequencies (A18, A19)
    f_upper = VPiezo / (tPiezo + denTopEl*tTopEl/denPiezo + .5*denBotEl*tBotEl/denPiezo );
    f_lower = VSi    / (tSi + .5*denBotEl*tBotEl/denSubst);
    
    % ratio of "half wavelength resonant frequencies (12)
    R = f_upper/f_lower;
    
    mN = round(R);
    mT = round(.5*(R-1));
    
    %Gamma = 1 - 2*denPiezo*VPiezo*( 1 + 2*pi*)
    
    % Frequency at mode order mN (A11) and mT (A28)
    fN = (mN*zSubst + 1)/2 / (zTopEl*TTopEl + zSubst*TSubst + zBotEl*TBotEl + TPiezo);
    fT = (mT + 1/2 + zSubst/2)/(2*zSubst) / (TSubst/zSubst + TBotEl/zBotEl + TPiezo + zTopEl*TTopEl);    
    
    % "normal" region multiplication factor to compute kt2 from keff2 (15)
    kt2Nfactor = mPiezo*(mPiezo + mTopEl + mBotEl + mSubst) / ...
                        (mPiezo + mTopEl + mBotEl/2)^2;
    
    % "tranisition" region multiplication factor to compute kt2 from keff2 (16)
    kt2Tfactor = mPiezo*(1 + mTopEl/mPiezo ...
                           + TSubst^2/mSubst * mPiezo/TPiezo^2 ...
                           + (denPiezo^2*VPiezo^2/(denSi^2*VSi^2) + 1) ...
                              * mBotEl/mPiezo/2) ...
                        / (mPiezo + mTopEl + mBotEl/2);
                    
    data = var2struct(delta_f0, delta_fN, delta_fT, f_upper, f_lower, ...
        R, mN, mT, fN, fT, kt2Nfactor, kt2Tfactor);
end

% params.C33         = C33;         % stifness constants AlScN
% params.Eps33       = Eps33;       % Permittivity constants
% params.E33         = E33;         % piezo constants AlScN
% params.CD33        = CD33;        % Stiffened stiffness 
% params.h33         = h33;   
% params.k2          = k2;          % coupling constant

% params.denPiezo    = denPiezo;  % AlScN mass density
% params.denSiO      = denSiO;    % Mass density of SiO2
% params.denSi       = denSi;     % Mass density of Si
% params.denPt       = denPt;     % Mass density of electrodes (Pt)
% 
% params.tPiezo      = tPiezo;    % AlScN layer thickness
% params.tTopE       = tTopE;     % Top electrode thickness
% params.tBotE       = tBotE;     % Bottom electrode thickness
% params.tSiO        = tSiO;      % SiO2 thickness
% params.tSi         = tSi;       % Si thickness% 

% params.phPiezo    = phPiezo;      %Phase delay AlScN
% params.phSiO      = phSiO;        %Phase delay of SiO2 longitudinal mode
% params.phSi       = phSi;         %Phase delay of Si longitudinal mode
% params.phTopEl    = phTopEl;      %Phase delay of top electrodes longitudinal mode
% params.phBotEl    = phBotEl;      %Phase delay of bottom electrodes longitudinal mode
% 
% params.ZPiezo     = ZPiezo;       %Longitudinal wave impedance of active layer 
% params.ZSiO       = ZSiO;         %Longitudinal wave impedance of SiO2 
% params.ZSi        = ZSi;          %Longitudinal wave impedance of Si
% params.ZPt        = ZPt;          %Longitudinal wave impedance of Pt
% 
% params.Ztop       = Ztop;
% params.Zbot       = Zbot;
% 
% params.VPiezo  	  = VPiezo;       %Longitudinal velocity in AlScN 
% params.VSi     	  = VSi;          %Longitudinal velocity in Si 
% params.VSiO    	  = VSiO;         %Longitudinal velocity in SiO2
% params.VPt     	  = VPt;          %Longitudinal velocity in Pt
% 
% params.gamma      = gamma;