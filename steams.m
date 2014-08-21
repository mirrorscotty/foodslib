     function [sys,x0] = steams(t,x,u,flag)  

%  STEAMS    System input steam or water source variables u-array re-assignment  
%
%            For use with FOODS-LIB models  
%
%            Synopsis: S-Function builds standard steam/water stream array.
%            Steam property data is checked to determine if valid and 
%            consistent. Extends u_array from 5 inputs to "usizemax" 
%            variables. U_array elements are re-positioned exactly as they 
%            are for food streams. U_array elements 37-"usizemax" are 
%            available for passing new information.
%  
%            Call: steams(t,x,u,flag)
%
%            Passed Arguments:
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 information
%
%            Returned Arguments:
%             sys - SIMULINK S-function variable: for flag = 0,
%                   sys = sizes of parameters and initial conditions.
%                   For flag = 3, sys = system outputs which are the
%                   u-arrays for stream entering the model.

%  ********************* SIMULINK INTERFACE INFORMATION *********************
 
usize = usizemax;                   % Current u_array size for streams
noinputs = 5;                       % Current number of constant inputs
 
if flag == 0
   % sys(3) is number of outputs
   % sys(4) is number of inputs
   sys = [0, 0, usize, noinputs, 0, 1];
   x0 = ones(size(u));

elseif abs(flag) == 1
   % if flag=1, return state derivatives, xDot
     disp('Set Simulation-Parameters Stop Time to Zero.')
     break;

elseif flag == 3
   % if flag=3, return system outputs, y

%    Steam/water Stream Variables  
 
%  *********************** ENERGY BALANCE LEVEL ***********************
 
     StempC = u(1);            % Temperature (C)
     Spres  = u(2);            % Pressure (MPa)
     Sx     = u(3);            % Quality/State

     Tref = 273.15;            % Reference temperature (K)   
     StempK = StempC + Tref;   % Temperature (K); used for validating state  
 
%    Steam reference temperatures (K) and pressures (MPa) 
     Tlow = 273.16;  
     Thigh = 647.3;  
     Plow = 0.000611;  
     Phigh = 22.1;  

%    Saturated steam conditions check
     if (Sx >= 0) & (Sx <= 1) 
       if (Spres > Phigh)
         errortxt = [10 'STEAM/WATER SOURCE ERROR: Pressure exceeds the critical pressure' 10 ...
         'which is 22.1 MPa.' 10];  
         disp(errortxt);  
         return;
       end 
       if (StempK > Thigh)
         errortxt = [10 'STEAM/WATER SOURCE ERROR: Temperature exceeds the critical temperature ' 10 ... 
         'which is 374.15 C (647.3 K). ' 10]; 
         return;  
       end  
 
       if ((StempK <= Tlow) & (Spres > Plow)) StempK = tpsteam(Spres); end  
       if ((Spres <= 0) & (StempK > Tlow)) Spres = ptsteam(StempK); end  
       if (StempK > 0) & (Spres > 0) 
         Tsat = tpsteam(Spres); 
         if (Tsat-StempK > 1) 
           warntxt = [10 'STEAM/WATER SOUCE WARNING: Steam conditions are over specified' 10 ...
           'for steam in liquid-vapor region. Based on the given pressure, ' 10 ...  
           'a new saturated steam temperature has been determined.' 10];   
           disp(warntxt)
           StempK = Tsat; 
         end 
       end     
       if (Spres < Plow)
         errortxt = [10 'STEAM/WATER SOURCE ERROR: Pressure is below the triple point.' 10];
         disp(errortxt);
         return;  
       end   
       if (StempK < Tlow) 
         errortxt = [10 'STEAM/WATER SOURCE ERROR: Temperature is below the triple point.' 10]; 
         disp(errortxt);
         return; 
       end  
     end 

%    Superheated steam conditions check  
     if (Sx > 1)
       if (StempK == 0) | (Spres == 0)
         errortxt = [10 'STEAM/WATER SOURCE ERROR: Both temperature and pressure need to be' 10 ...
         'specified for superheated steam.' 10];
         disp(errortxt);
         return;
       end
       if (Spres > Phigh) 
         errortxt = [10 'STEAM/WATER SOURCE ERROR: Pressure exceeds the critical pressure.' 10];
         disp(errortxt);
         return; 
       end
       if (StempK > Thigh) 
         errortxt = [10 'STEAM/WATER SOURCE ERROR: Temperature exceeds the critical temperature.' 10]; 
         disp(errortxt);
         return;   
       end       
       if (Spres < Plow) 
         errortxt = [10 'STEAM/WATER SOURCE ERROR: Pressure is below the triple point.' 10];
         disp(errortxt);
         return;   
       end       
       if (StempK < Tlow)  
         errortxt = [10 'STEAM/WATER SOURCE ERROR: Temperature is below the triple point.' 10]; 
         disp(errortxt);
         return;  
       end      
       if (StempK > 0) & (Spres > 0) 
         Tsat = tpsteam(Spres); 
         if (Tsat > StempK) 
           errortxt = [10 'STEAM/WATER SOURCE ERROR: Temperature is too low for superheated' 10 ...
           'steam conditions.' 10];
           disp(errortxt);
           return;  
         end 
       end  
     end   

%    Compressed water conditions check
     if (Sx < 0)  
       errortxt = [10 'STEAM/WATER SOURCE ERROR: Compressed water properties are currently' 10 ...
       'unavailable.' 10];    
       disp(errortxt);
       return;
     end  
  
%  ************************* ECONOMICS LEVEL **************************
 
     Scost = u(4);            % Cost ($/kg)
     Sinfl = u(5);            % Annual cost inflation rate
 
%  ************** Re-ORGANIZATION &  EXTENSION OF U-ARRAY *************
 
%    Re-organization of u-array
     for I = ulocate('mass'):ulocate('ctypelast')  
       u(I) = 0;
     end
     u(ulocate('comp1st'))  = 1.0; 
                              % First component fraction is all water (w/w) 
     u(ulocate('ctype1st')) = 1.0;   
                              % First component fraction is of type water (=1) 
     u(ulocate('temp')) = StempK - Tref;          
     u(ulocate('pressure')) = Spres;         
     u(ulocate('quality')) = Sx;          
     u(ulocate('cost')) = Scost;        
     u(ulocate('inflation')) = Sinfl;  
 
%    Extension of u-array to "usizemax" elements
     for I = 37:usize  
       u(I) = 0;
     end 

     sys = [u;];

else
   % Otherwise, no need to return anything
   % since this is continuous system
   sys = [];  

end



 
