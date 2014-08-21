
      function[Fk] = thcond(Fme) 

%  THCOND     Thermal conductivity of any stream type  
%
%             Synopsis: MATLAB function to determine the thermal  
%             conductivity of any stream. Specifically checks for pure 
%             water streams and treats them as steam streams. 
%    
%             Call: thcond(Fme)
%
%             Passed Arguement: 
%              Fme - mass and energy information array
%
%             Returned Argument
%              Fk - thermal conductivity (kJ/h-m-K)   

%  **********************************************************************  

%    Reference temperature (K)
     Tref = 273.15;

%    Re-assign passed information 
     F = Fme(1);                 % mass flow rate (kg/h)    
     Fcomp = Fme(2:11);          % component mass fractions (w/w)  
     Fctype = Fme(12:21);        % mass flow rate (kg/h)      
     FtempK = Fme(22);           % temperature (K)  
     Fpres = Fme(23);            % pressure (MPa) 
     Fx = Fme(24);               % quality  

%    Check for multi-component food stream 
     flagmc = purechk(Fcomp,Fctype);  
 
%    For pure water stream 
     if (flagmc == 1)
%      Use steam table algorithm hsteam.m to find enthalpy (kJ/kg) 
       Fk = ksteam(FtempK,Fpres,Fx);  
     end        

%    For multi-component streams, call heat capacity function cp.m
     if (flagmc == 0)
%      Call heat capacity function cp.m (kJ/kg-K)   
       Fk = thermc(Fcomp,Fctype,FtempK);     
     end    


