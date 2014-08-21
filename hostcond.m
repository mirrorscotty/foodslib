      function[hoarray] = hostcond(Tsteam,Psteam,Twall,Do,length,vertical,Ntubes)

%  HOSTCOND  Steam side heat transfer coefficient for film-condensation 
%
%            Synopsis: MATLAB function which determines the steam side
%            heat transfer coefficient for film-condensation. All models  
%            assume saturated steam inlet conditions. Models include 
%            laminar flow and vertical tubes, turbulent flow and vertical
%            tubes, and laminar flow and vertical tier of horizontal tubes.  
%
%            Call: hostcond(Tsteam,Psteam,Twall,Do,length,vertical,Ntubes) 
%    
%            Passed Arguments: 
%             Tsteam - temperature (K)
%             Psteam - pressure (MPa)
%             Twall - tube wall temperature (K) 
%             Do - tube outer diameter (m) 
%             length - tube length (m) 
%             vertical - tube orientation: horizontal (0); veritcal (1)  
%             Ntubes - number of tubes in vertical tier of horizontal 
%                      tubes (value only required for horizontal config.) 
%
%            Returned Argument:
%             ho = steam side heat transfer coefficient (kJ/h-m^2-K)
% 
%            References:
%             Geankoplis, C.J. 1983.  "Transport Processes and Unit
%               Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.

%  *************************************************************************

%     Saturated Temperature (K) based on pressure (MPa) of steam  
      if (Psteam > 0) 
        Tsat = tpsteam(Psteam);  
      else 
        Tsat = Tsteam;
      end  

%     Average fluid temperature (K)
      Tf = (Tsat+Twall)/2;   
 
%     Density of saturated vapor (kg/m^3)   
      rhovap = rhosteam(Tf,0,1);  

%     Density of saturated liquid (kg/m^3)   
      rholiq = rhosteam(Tf,0,0);  

%     Enthalpy of saturated vapor (kJ/kg)
      Hvap = hsteam(Tsat,0,1);  

%     Enthalpy of saturated liquid (kJ/kg)  
      Hliq = hsteam(Tsat,0,0); 

%     Latent heat of vaporization (kJ/kg)   
      Hlat = Hvap-Hliq; 

%     Gravity (m/s^2)
      g = gravity;  

%     Thermal conductivity of saturated liquid (kJ/h-m-K) --> (kJ/s-m-K)   
      kliq = ksteam(Tf,0,0)/3600; 

%     Viscosity of saturated liquid (Pa-s)
      visliq = vissteam(Tf,0,0);

%     Temperature difference (K) between saturated vapor temperature (K) and 
%     tube wall temperature (K)      
      delT = Tsat - Twall;  

%     For vertical tubes
      if (vertical == 1) 

%       Assume laminar flow; Calculate Nusselt number and heat transfer
%       coefficient; Check laminar assumption  

%       Nusselt Number for laminar vertical film-condensation     
        NNu = 1.13*((rholiq*(rholiq-rhovap)*g*Hlat*length^3)/ ...
              (visliq*kliq*delT))^(0.25);

%       Heat transfer coefficient (kJ/s-m^2-K) --> (kJ/h-m^2-K)  
        ho = NNu*kliq/length*3600; 

%       Tube surface area (m^2)
        Area = pi*Do*length;

%       Total heat transfer per tube (kJ/h)
        Qtube = ho*Area*delT;  

%       Condensation mass flow rate per tube (kg/h) --> (kg/s)   
        C = Qtube/Hvap/3600;  

%       Reynolds Number 
        NRe = 4*C/(pi*Do*visliq); 

%       Check on laminar assumption 
        if (NRe > 1800)        % Turbulent flow 

%         Need simultaneous solution for Reynolds and Nusselt Numbers 
          iter = 0;  
          NRediff = NRe;
          while(NRediff > 1e1) 
            iter = iter + 1;   
            if (iter == 30) break; end  

%           Nusselt Number for turbulent vertical film-condensation     
            NNu = 0.0077*((g*rholiq^2*length^3)/visliq^2)^(1/3)*(NRe)^0.4;

%           Reynolds Number used in this iteration 
            NRelast = NRe;

%           Heat transfer coefficient (kJ/s-m^2-K) --> (kJ/h-m^2-K)
            ho = NNu*kliq/length*3600; 
 
%           Total heat transfer per tube (kJ/h) 
            Qtube = ho*Area*delT;
 
%           Condensation mass flow rate per tube (kg/h) --> (kg/s)    
            C = Qtube/Hvap/3600;   
 
%           Reynolds Number  
            NRe = 4*C/(pi*Do*visliq);  

%           While loop criteria 
            NRediff = NRe - NRelast; 

          end  % while loop for NRe 
        end  % turbulent flow model  

      end  % vertical tube models     

%     For vertical tier of horizontal tubes placed one below the other 
      if (vertical == 0) 

%       Most practical applications are laminar flow  
%       Reynolds Number  
        NRe = 0; 

%       Number of tubes in vertical tier; assumes arranged in a square
        Nvert = ceil(sqrt(Ntubes));

%       Nusselt Number for laminar horizontal film-condensation
        NNu = 0.725*((rholiq*(rholiq-rhovap)*g*Hvap*Do^3)/ ...
              (Nvert*visliq*kliq*delT))^(0.25);  
 
%       Heat transfer coefficient (kJ/s-m^2-K) --> (kJ/h-m^2-K)
        ho = NNu*kliq/Do*3600;   

      end  % end horizontal tube model   

      hoarray = [ho;NRe;NNu;];

 
 
