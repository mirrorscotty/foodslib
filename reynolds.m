       function [NRe] = reynolds(D,fthick,v,rho,K,n,mu)

%  REYNOLDS  Reynolds number  
%
%            Synopsis: MATLAB function to determine Reynolds number  
%            using the generalized equation.   
%
%            Call: reynolds(D,fthick,v,rho,K,n,mu) 
%
%            Passed Arguments: 
%             D - characteristic diameter (m)  
%             fthick - film thickness (m)   
%             v - mass average velocity (m/s) 
%             rho - density (kg/m^3) 
%             K - consistency coefficient (Pa-s^n)
%             n - flow behaviour index  
%             mu - average viscosity (Pa-s) 
%
%            Returned Arguments:
%             NRe - Reynolds number 
%
%            References: 
%            [1] Geankoplis, C.J. 1983.  "Transport Processes and Unit
%                  Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.

%  *************************************************************************

%       For full tube 
        if (fthick >= D) 

%         Numerator for NRe determination 
          Num = D^n * v^(2-n) * rho;

%         Demoninator for NRe determination 
          Den = K * 8^(n-1) * ((3*n+1)/(4*n))^n;  

%         Reynolds number 
          NRe = Num/Den;  

%       For film in tube 
        else 

%         Gravity (m/s^2)
          g = gravity;  

%         Mass flow rate per wetted perimeter (kg/s-m) 
          Gamma = n/(2*n+1)*rho*(rho*g/K)^(1/n)*fthick^((2*n+1)/n);  

%         Reynold Number
          NRe = 4*Gamma/mu;  

        end   
 
