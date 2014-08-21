     function [sys,x0] = airs(t,x,u,flag)  

%  AIRS      System input air-water variables u-array re-assignment  
%
%            For use with FOODS-LIB models  
%
%            Synopsis: S-Function builds standard stream array using air
%            water inputs. Extends u_array from 3 inputs to "usizemax" 
%            variables. U_array elements are re-positioned exactly as they 
%            are for food streams. U_array elements 37-"usizemax" are 
%            available for passing new information.
%  
%            Call: airs(t,x,u,flag)
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
noinputs = 3;                       % Current number of constant inputs
 
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
 
     AtempC = u(1);            % Temperature (C)
     Apres  = u(2);            % Pressure (MPa)
     Arh    = u(3);            % Relative humidity   

     Tref = tref;              % Reference temperature (K)   
     AtempK = AtempC + Tref;   % Temperature (K)  

%    Absolute relative humidity (kg h2o/kg dry air)
     Aabsh = airwater(2,AtempK,Apres,1,Arh); 

%    Weight fraction of water 
     Axh2o = Aabsh/(1+Aabsh); 
%    Weight fraction of air  
     Axair = 1/(1+Aabsh);   

%  ************** Re-ORGANIZATION &  EXTENSION OF U-ARRAY *************
 
%    Re-organization of u-array
     for I = ulocate('mass'):usize    
       u(I) = 0;
     end
     u(ulocate('comp1st'))  = Axh2o; % Water (w/w) 
     u(ulocate('comp1st')+1)  = Axair; % Air (w/w) 
     u(ulocate('ctype1st')) = comptype('water');   
     u(ulocate('ctype1st')+1) = comptype('air');   
     u(ulocate('temp')) = AtempK - Tref;          
     u(ulocate('pressure')) = Apres;         
 
     sys = [u;];

else
   % Otherwise, no need to return anything
   % since this is continuous system
   sys = [];  

end



 
