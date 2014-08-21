       function [gwall] = shrwall(n,Q,Din)

%   SHRWALL  Shear rate at the pipe wall   
%
%            Synopsis: MATLAB function to determine the shear rate at 
%            the pipe wall. Accounts for Non-Newtonian fluids.    
%
%            Call: shrwall(n,Q,Din) 
%
%            Passed Arguments: 
%             n - flow behavior index  
%             Q - volumetric flow rate (m^3/h)   
%             Din - inner pipe diameter (m)   
%
%            Returned Arguments:
%             gwall - shear rate (1/s)     
%
%            References: 
%            [1] Steffe, J.F. 1992. Rheological Methods in Food Process
%                 Engineering. Freeman Press, East Lansing, Mich.  

%  *************************************************************************

%     Newtonian shear rate at the pipe wall (1/s)  
      Ngwall = (4*Q/3600)/(pi*(Din/2)^3);

%     Non-Newtonian sheat rate at the pipe wall (1/s) 
      gwall = (3*n+1)/(4*n)*Ngwall;  

 
