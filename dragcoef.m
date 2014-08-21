    function [coeff] = dragcoef(NRe)  

%  DRAGCOEF  Drag coefficient on a sphere 
%  
%            Synopsis: MATLAB function which determined the drag coefficient
%            on a sphere given Reynold Number
%
%            Call: dragcoef(NRe)
%
%            Passed Argument:
%             NRe - Reynolds Number 
%
%            Returned Argument:   
%             coeff - Drag coefficient   
%
%            References:
%             [1] Geankoplis, C.J. 1983.  "Transport Processes and Unit
%               Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.
%             [2] Perry's
%

%  ****************************************************************************  

%   Drag coefficient on droplet (Perry's)   
    if (NRe < 1e-1)  
      CD = 24/NRe;  
    elseif (NRe >= 1e-1) & (NRe < 1e3)  
      CD = 24/NRe*(1+0.4*NRe^0.7);  
    elseif (NRe >= 1e3) & (NRe < 3.5e5)  
      CD = 0.445;                                 
    elseif (NRe > 5e5) & (NRe <= 1e6)   % (Geankopolis)
      CD = 0.2;                            
    elseif (NRe > 1e6)     
      CD = 0.19 - 8e4/NRe;  
    else          % Linear interpolation in semi-turbulent region   
      CD = (NRe-3.5e5)*(0.2-0.445)/(5e5-3.5e5) + 0.445;       
    end  

    coeff = CD;  



