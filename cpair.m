      function [heatair] = cpair(Tair) 

%  CPAIR     Air heat capacity   
%
%            Synopsis: MATLAB function which returns heat capacity  
%            of air based on air temperature at atmospheric pressure
%            Linear interpolation of data in temperature range of
%            -17.8 to 260 C is used
%
%            Call: cpair(Tair)
%
%            Passed Argument:
%             Tair - air temperature (K)
%
%            Returned Argument:
%             heatair - heat capacity of air (kJ/kg-K)
%
%            Reference:
%             Geankoplis, C.J. 1983.  "Transport Processes and Unit
%              Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.
%

%  **************************************************************************
 
%     Temperature data (K)  

      T = [255.4; 273.2; 283.2; 311.0; 338.8; 366.5; 
           394.3; 422.1; 449.9; 477.6; 505.4; 533.2;];  

%     Heat capacity data (kJ/kg-K)  

      cp = [1.0048; 1.0048; 1.0048; 1.0048; 1.0090; 1.0090; 
            1.0132; 1.0174; 1.0216; 1.0258; 1.0300; 1.0341;];  

%     Use temperature (K) to find heat capacity (kJ/kg-K) in a 
%     1D-look-up table.
      if (Tair < min(T))
        heatair = interp1(T,cp,min(T));   
      elseif (Tair > max(T))
        heatair = interp1(T,cp,max(T));   
      else
        heatair = interp1(T,cp,Tair);   
      end 
 

