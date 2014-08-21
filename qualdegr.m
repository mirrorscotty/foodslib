       function [Qout] = qualdegr(Dvalq,DtempKq,Zvalq,Eaq,Qoin,Tbulk,time)

%  QUALDEG   Thermal degradation of quality    
%
%            Synopsis: MATLAB function to determine the final  
%            quality following heat treatment for a specified
%            time at a specified temperature.  
%
%            Call: qualdegr(Dvalq,DtempKq,Zvalq,Eaq,Qoin,Tbulk,time)  
%
%            Passed Arguments: 
%             Dvalq -  D value (min)
%             DtempKq - D value temperature basis (K) 
%             Zvalq - Z value (C)
%             Eaq - activation energy (kcal/mol) 
%             Qoin -  initial microbial load 
%             Tbulk - bulk processing temperature (K)  
%             time - minimum processing time (min)   
%
%            Returned Arguments:
%             Qout - final microbial load  
%
%            References: 
%            [1] Geankoplis, C.J. 1983.  "Transport Processes and Unit
%                  Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.

%  *************************************************************************

%     Microbial kinetics for case where the activation energy is given
      if (Eaq > 0)   

%       Activation energy units conversion (kJ/mol --> kcal/mol)
        Eaq = Eaq/4.1868;

%       Gas law constant (kcal/mol-K)
        R = 1.9872e-3;   

%       Reaction velocity constant (1/min)
        kq = 2.303/Dvalq;

%       Arrhenius empirical constant (1/min)
        aq = kq/(exp(-Eaq/(R*DtempKq)));

%       Final microbial load   
        Qout = exp(log(Qoin)-(aq*(exp(-Eaq/(R*Tbulk))))*time);

%     Microbial kinetics for case where D and z are specified 
      else 

%       Decimal reduction time at bulk temperature (min) 
        DTbulk = Dvalq * 10^((DtempKq - Tbulk)/Zvalq);  

%       Final microbial load  
        Qout = Qoin / 10^(time/DTbulk);

      end   
 
