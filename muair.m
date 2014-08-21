      function [visair] = muair(Tair) 
    
%  MUAIR     Air Viscosity       
%
%            Synopsis: MATLAB function which returns viscosity of air 
%            based on air temperature at atmospheric pressure  
%            Linear interpolation of data in temperature range of 
%            -17.8 to 260 C is used    
%
%            Call: muair(Tair) 
%
%            Passed Argument: 
%             Tair - air temperature (K) 
%
%            Returned Argument:
%             visair - viscosity of air (Pa-s)   
%  
%            Reference: 
%             Geankoplis, C.J. 1983.  "Transport Processes and Unit
%              Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.
%

%  **************************************************************************
 
%     Temperature data (K)  

      T = [255.4; 273.2; 283.2; 311.0; 338.8; 366.5; 
           394.3; 422.1; 449.9; 477.6; 505.4; 533.2;];  

%     Viscosity data (Pa-s*10^5)  

      mu = [1.62; 1.72; 1.78; 1.90; 2.03; 2.15; 
            2.27; 2.37; 2.50; 2.60; 2.71; 2.80;];  

%     Viscosity (Pa-s)  
      mu = mu*1e-5;  

%     Use temperature (K) to find viscosity (Pa-s) in a 
%     1D-look-up table.
      if (Tair < min(T))
        visair = interp1(T,mu,min(T));
      elseif (Tair > max(T))
        visair = interp1(T,mu,max(T)); 
      else 
        visair = interp1(T,mu,Tair);   
      end   
 

