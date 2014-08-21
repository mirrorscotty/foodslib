      function [pinacl] = osmpnacl(conc) 
    
%  OSMPNACL  Osmotic pressure of sodium chloride solution        
%
%            Synopsis: MATLAB function which returns the osmotic 
%            pressure of a sodium chloride solution. Range is for 
%            solution concentrations of 0 to 2 g mol NaCl/kg water
%            at 25 C.
%
%            Call: osmpnacl(conc) 
%
%            Passed Argument: 
%             conc - NaCl concentration (g mol NaCl/kg water) 
%
%            Returned Argument:
%             osmpnacl - osmotic pressure (Pa)   
%  
%            Reference: 
%             Geankoplis, C.J. 1983.  "Transport Processes and Unit
%              Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.
%

%  **************************************************************************
 
%     NaCl concentration (g mol NaCl/kg water)  

      C = [0 0.01 0.1 0.5 1 2];  

%     Osmotic pressure (atm)    

      OP = [0 0.47 4.56 22.55 45.80 96.2]; 

%     Osmotic pressure (atm --> Pa)  
      OP = OP*1.01325e5;  

%     Use concentration to find osmotic pressure in a 
%     1D-look-up table.
      if (conc < min(C))
        visair = interp1(C,OP,min(C));
      elseif (conc > max(C))
        visair = interp1(C,OP,max(C)); 
      else 
        visair = interp1(C,OP,conc);   
      end   
 

