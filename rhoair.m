      function [denair] = rhoair(Tair,Pair) 

%  RHOAIR    Air density  
%
%            Synopsis: MATLAB function which determines density of air
%            based on deal gas law  
%
%            Call: rhoair(Tair,Pair)
%
%            Passed Argument:
%             Tair - air temperature (K)
%             Pair - air pressure (MPa) 
%   
%            Returned Argument:
%             rhoair - density of air (kg/m^3)
%

%  **************************************************************************

%     Gas law constant (kg-m^2/s^2-kgmol-K)  
      R = gaslawc;   

%     Molecular weight of air
      MWair = molwts(comptype('air'));   

%     Specific volume of air (m^3/kg)  
      vair = R*Tair/(MWair*Pair*1e6);    

%     Density of air (kg/m^3)   
      denair = 1.0/vair;    
  

