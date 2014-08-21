      function [tcair] = kair(Tair) 

%  KAIR      Air thermal conductivity  
%
%            Synopsis: MATLAB function which returns thermal conductivity
%            of air based on air temperature at atmospheric pressure
%            Linear interpolation of data in temperature range of
%            -17.8 to 260 C is used
%
%            Call: kair(Tair)
%
%            Passed Argument:
%             Tair - air temperature (K)
%
%            Returned Argument:
%             tcair - thermal conductivity of air (kJ/h-m-K)
%
%            Reference:
%             Geankoplis, C.J. 1983.  "Transport Processes and Unit
%              Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.
%

%  **************************************************************************

%     Temperature data (K)  

      T = [255.4; 273.2; 283.2; 311.0; 338.8; 366.5; 
           394.3; 422.1; 449.9; 477.6; 505.4; 533.2;];  

%     Thermal conductivity data (W/m-K)  

      k = [0.02250; 0.02423; 0.02492; 0.02700; 0.02925; 0.03115;
           0.03323; 0.03531; 0.03721; 0.03894; 0.04084; 0.04258;];  

%     Use temperature (K) to find thermal conductivity (W/m-K) in a 
%     1D-look-up table.
      if (Tair < min(T))
        thermair = interp1(T,k,min(T));   
      elseif (Tair > max(T))
        thermair = interp1(T,k,max(T));   
      else
        thermair = interp1(T,k,Tair);   
      end  

%     Thermal conductivity units conversion (W/m-K --> kJ/h-m-K)  
      tcair = thermair*3.6;   
 

