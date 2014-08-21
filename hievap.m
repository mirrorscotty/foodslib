      function[hiarray] = hievap(Fme,Lme,evapdesc)

%  HIEVAP    Product side heat transfer coefficient models   
%
%            Synopsis: MATLAB function which determines inside tube heat
%            transfer coefficient. 
%
%            Call: hievap(Fme,Lme,evapdesc)
% 
%            Passed Arguments: 
%             Fme -  feed mass and energy information array  
%             Lme -  liquid concentrate mass and energy information array  
%             evapdesc - evaporator description array: 
%                  (1) - Din - inside tube diameter (m) 
%                  (2) - fthick - film thickness (m) 
%                  (3) - vmin - minimum avergae velocity (m/s)    
%                  (4) - length - tube length (m)
%                  (5) - Twall - tube wall temperature (K) 
%                  (6) - evaptype - evaporator model 
%                  (7) - mumodel - viscosity model 
%
%            Returned Arguments: 
%             hiarray(1) - heat transfer coefficient (kJ/h-m^2-K)  
%             hiarray(2) - Reynolds number 
%             hiarray(3) - Prandtl number
%             hiarray(4) - Nusselt number  
%
%            References:
%            [1] Geankoplis, C.J. 1983.  "Transport Processes and Unit
%                  Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.
%            [2] Billet, R. 1989. "Evaporation Technology: Principles,
%                  Applications, Economics," Verlagsgesellschaft, mbH, 
%                  Germany.
%            [3] Hartel, R.W. 1992. Evaporation and Freeze Concentration
%                 In "Handbook of Food Engineering," Heldman, D.R. and 
%                 Lund, D.B. (Eds.) Marcel Deckker, Inc. NY.  

%  **************************************************************************

%     Feed conditions: mass flow rate (kg/h); composition (w/w), 
%       component types, temperature (K)
      F = Fme(ulocate('mass'));
      Fcomp = Fme(ulocate('comp1st'):ulocate('complast')); 
      Fctype = Fme(ulocate('ctype1st'):ulocate('ctypelast'));  
      FtempK = Fme(ulocate('temp'));  

%     Liquid Concentrate conditions: mass flow rate (kg/h); 
%       composition (w/w), component types, temperature (K)
      L = Lme(ulocate('mass')); 
      Lcomp = Lme(ulocate('comp1st'):ulocate('complast')); 
      Lctype = Lme(ulocate('ctype1st'):ulocate('ctypelast'));
      LtempK = Lme(ulocate('temp')); 
 
%     Evaporator description variables 
      Din = evapdesc(1);         % inside tube diameter (m)
      fthick = evapdesc(2);      % film thickness (m)   
      vmin = evapdesc(3);        % minimum average velocity (m/s)   
      length = evapdesc(4);      % tube length (m)
      Twall = evapdesc(5);       % tube wall temperature (K) 
      evaptype = evapdesc(6);    % evaporator model 
      mumodel = evapdesc(7);     % viscosity model  

%     Gravity (m/s^2)
      g = gravity;

%     Bulk temperature (K) 
      Tbulk = (FtempK+LtempK)/2;  

%     Power law values for viscosity at bulk temperature  
      powerlaw = viscos(Fcomp,Fctype,Tbulk,mumodel);   
      Kbulk = powerlaw(1);       % Consistency Coefficient (Pa-s^n) 
      nbulk = powerlaw(2);       % Flow behavior index  

%     Bulk viscosity (Pa-s); taken at max velocity in tube  
      mubulk = Kbulk; 

%     Bulk density (kg/m^3)   
      rhobulk = rho(Fcomp,Fctype,Tbulk); 

%     Power law values for viscosity at wall temperature
      powerlaw = viscos(Fcomp,Fctype,Twall,mumodel);
      Kwall = powerlaw(1);       % Consistency Coefficient (Pa-s^n)
      nwall = powerlaw(2);       % Flow behavior index
 
%     Volumetric flow rate per tube (m^3/h)
      Vol = rhobulk*vmin*pi*((Din/2)^2-((Din-2*fthick)/2))*3600; 

%     Maximum shear rate (1/s) 
      if (fthick >= Din)  
        shearmax = shrwall(nwall,Vol,Din);
      else 
        shearmax = rhobulk*g*fthick/mubulk;  
      end  
 
%     Viscosity at tube wall (Pa-s)
      muwall = Kwall*shearmax^(nwall-1); 

%     Average viscosity (Pa-s)
      muavg = (mubulk+muwall)/2;

%     Reynolds number  
      NRe = reynolds(Din,fthick,vmin,rhobulk,Kbulk,nbulk,muavg);   

%     Bulk heat capacity (kJ/kg-K)
      cpbulk = cp(Fcomp,Fctype,Tbulk);  

%     Bulk thermal conductivity (kJ/h-m-K)
      kbulk = thermc(Fcomp,Fctype,Tbulk);  

%     Prandtl number  
      NPr = prandtl(cpbulk,mubulk,kbulk);

      if (evaptype == 1)
%       Natural circulation method: Boarts, Badger and Meisenburg [Ref. 2]
%       Nusselt Number
        NNu = 0.035*NRe^0.8*NPr^0.4;
%       Heat transfer coefficient (kJ/h-m^2-K)
        hi  = NNu*kbulk/Din;

      elseif (evaptype == 10) 
%       Forced circulation; conservative assumption no boiling in 
%        tubes [Ref 1] 
%       Nusselt Number 
        NNu = 0.028*NRe^0.8*NPr^(1/3)*(mubulk/muwall)^0.14;  
%       Heat transfer coefficient (kJ/h-m^2-K)  
        hi  = NNu*kbulk/Din; 

      elseif (evaptype == 11)  
%       Forced circulation method: Nusselt-Kraussold equation [Ref. 2]   
%       Good for Twall-Tin <= 6 K and NRe > 10000   
%       Nusselt Number   
        NNu = 0.032*(NRe^0.8)*(NPr^0.37)*(length/Din)^(-0.054); 
%       Heat transfer coefficient (kJ/h-m^2-K)   
        hi  = NNu*kbulk/Din;   
%       Check on conditions:  
%       if (Twall-Tbulk > 6) | (NRe < 1e4)  
%         disp('EVAPORATOR WARNING MESSAGE: Conditions of Twall-Tbulk < 6 K and ...
%         NRe > 10000 required for the inside heat transfer coefficient model for ...
%         a forced circulation evaporator are not being met.');
%         delT = Twall - Tbulk; 
%         NRe 
%       end   

      elseif (evaptype == 20) 
%       Falling film model; Chun and Seban [Ref. 3]
%       Laminar flow
        if NRe < 2100 
          C = 0.606; 
          n = -0.22;
          m = 0;
        else  % Turbulent flow 
          C = 0.0038; 
          n = 0.4; 
          m = 0.65;     
        end  
%       Dimensionless heat transfer coefficient (not really Nusselt number) 
        hdim = C * NRe^n * NPr^m;
        NNu = 0;
%       Heat transfer coefficient (W/m^2-K)
        hi = hdim*kbulk/((muavg^2/(rhobulk^2*g))^(1/3)); 

      elseif (evaptype == 21)

%       Falling film model 
%       Laminar flow; Angeletti and Moresi [Ref. 3] 
        if NRe < 2100
          C = 1.1; 
          n = -1/3; 
          m = 0;  
%         Dimensionless heat transfer coefficient (not really Nusselt number)
          hdim = C * NRe^n * NPr^m;
          NNu = hdim;
%         Heat transfer coefficient (W/m^2-K)
          hi = hdim*kbulk/((muavg^2/(rhobulk^2*g))^(1/3));  

        else % Turbulent flow; Narayana, Murphy & Sarma [Ref. 3]  
%         Dimensionless film thickness  
          fdim = exp(0.789 + 0.103*log(NRe) + 0.041*(log(NRe)^2));
%         Dimensionless heat transfer coefficient - numerator
          hdimnum = 0.89 * fdim^(1/3);
%         Dimensionless heat transfer coefficient - denominator
%         hdimden = 5 + [atan(2.73*sqrt(NPr)) - atan(0.455*sqrt(NPr))]/ ...
%                    [0.091*sqrt(NPr) + log(fdim/30)/(0.36*NPr)];
          hdimden = 5 + [atan(2.73*sqrt(NPr)) - atan(0.455*sqrt(NPr))]/ ...
                     [0.091*sqrt(NPr)];
%         Dimensionless heat transfer coefficient (not really Nusselt number)
%         hdim = hdimnum/hdimden;  
          hdim = hdimnum/hdimden + log(fdim/30)/(0.36*NPr);  
          NNu = hdim; 
%         Heat transfer coefficient (W/m^2-K)
          hi = hdim*kbulk/((muavg^2/(rhobulk^2*g))^(1/3));
        end 

      end  

      hiarray = [hi;NRe;NPr;NNu;];  
  
