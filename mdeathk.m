       function [Nout] = mdeathk(Dvalm,DtempKm,Zvalm,Eam,Noin,Tbulk,time)

%  MDEATHK   Thermal death-rate kinetics of microorganisms   
%
%            Synopsis: MATLAB function to determine the viable 
%            microbial load following heat treatment for a specified
%            time at a specified temperature.  
%
%            Call: mdeathk(Dvalm,DtempKm,Zvalm,Eam,Noin,FtempK,time) 
%
%            Passed Arguments: 
%             Dvalm -  D value (min)
%             DtempKm - D value temperature basis (K) 
%             Zvalm - Z value (C)
%             Eam - activation energy (kcal/mol) 
%             Noin -  initial microbial load 
%             Tbulk - bulk processing temperature (K)  
%             time - minimum processing time (min)   
%
%            Returned Arguments:
%             Nout - final microbial load  
%
%            References: 
%            [1] Geankoplis, C.J. 1983.  "Transport Processes and Unit
%                  Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.

%  *************************************************************************

%     Microbial kinetics for case where the activation energy is given
      if (Eam > 0)   

%       Activation energy units conversion (kJ/mol --> kcal/mol)
        Eamn = Eam/4.1868;

%       Gas law constant (kcal/mol-K)
        R = 1.9872e-3;   

%       Reaction velocity constant (1/min)
        km = 2.303/Dvalm;

%       Arrhenius empirical constant (1/min)
        am = km/(exp(-Eamn/(R*DtempKm)));

%       Final microbial load   
        Nout = Noin/exp((am*(exp(-Eamn/(R*Tbulk))))*time);
%       Nout = exp(log(Noin)-(am*(exp(-Eamn/(R*Tbulk))))*time);

%     Microbial kinetics for case where D and z are specified 
      else 

%       Decimal reduction time at bulk temperature (min) 
        DTbulk = Dvalm * 10^((DtempKm - Tbulk)/Zvalm);  

%       Final microbial load  
        Nout = Noin / 10^(time/DTbulk);

      end   
 
