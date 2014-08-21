       function [NPr] = prandtl(cp,mu,k)

%  PRANDTL   Prandtl number  
%
%            Synopsis: MATLAB function to determine Prandtl number  
%
%            Call: prandtl(cp,mu,k) 
%
%            Passed Arguments: 
%             cp - heat capacity (kJ/kg-K)  
%             mu - viscosity (Pa-s)
%             k - thermal conductivity (kJ/h-m-K) 
%
%            Returned Arguments:
%             NPr - Prandtl number  
%
%            References: 
%            [1] Geankoplis, C.J. 1983.  "Transport Processes and Unit
%                  Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.

%  *************************************************************************

%     Prandtl number (Converting hr to sec in thermal conductivity)
      NPr = cp*mu*3600/k;


