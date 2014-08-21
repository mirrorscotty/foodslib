
      function[Fh] = enthalpy(Fme) 

%  ENTHALPY   Enthalpy (kJ/kg) of any stream type  
%
%             Synopsis: MATLAB function to determine the enthalpy of 
%             any stream. Specifically checks for pure water streams 
%             and treats them as steam streams. 
%
%             Call: enthalpy(Fme)
%
%             Passed Arguement: 
%              Fme - mass and energy information array
%
%             Returned Argument:
%              Fh - enthalpy (kJ/kg)   

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
       Fh = hsteam(FtempK,Fpres,Fx);  
     end        

%    For air-water stream 
     if (flagmc == 2)  
%      Absolute humidity (kg h2o/kg dry air)  
       xh2o = findval1(Fctype,comptype('water'));
       Fw = Fcomp(xh2o)/(1-Fcomp(xh2o));  
%      Use airwater psychrometric airwater.m to find enthalpy (kJ/kg dry air) 
       Fhd = airwater(1,FtempK,Fpres,4,Fw); 
%      Air fraction 
       Fxair = 1-Fcomp(xh2o);  
%      Enthalpy (kJ/kg total) 
       Fh = Fhd*Fxair;       
     end   

%    For multi-component streams, call heat capacity function cp.m
     if (flagmc == 0)
%      Call heat capacity function cp.m (kJ/kg-K)   
       Fcp = cp(Fcomp,Fctype,FtempK);     
%      Enthalpy of stream  
       Fh = Fcp*(FtempK-Tref); 
     end    


