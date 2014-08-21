      function [cond] = ksteam(Tsteam, Psteam, xsteam)

%  KSTEAM    Thermal conductivity of steam
%
%            Synopsis: MATLAB function which determines the steam thermal 
%            conductivity (kJ/h-m-K) based on temperature, pressure, and 
%            quality.
%
%            Call: vissteam(Tsteam, Psteam, xsteam)
%
%            Passed Variables:
%             Tsteam - temperature (K)
%             Psteam - pressure (MPa)
%             xsteam - quality
%
%            Returned Variable:
%             cond = steam thermal conductivity (kJ/h-m-K)
%
%            Reference: Sengers, J.V. and Watson, J.T.R.  1985.  Improved
%            international formulation for the viscosity and thermal
%            conductivity of water substance.  J. Phys. Chem. Ref. Data
%            15(4): 1291-1314.
%            IAPS (International Association for the Properties of Steam)
%            Formulation 1985.

%  ***********************************************************************

%     For all steam thermal conductivity determination:
%     Critical temperature (K) 
      Tcrit = 647.27; 

%     Temperature (K) and pressure (MPa) combination 
%     Set variables for water and steam vapor 
      if (xsteam >= 0.0) & (xsteam <= 1.0) 
       if (Psteam > 0)
%         Saturation steam temperature (K); based on pressure
          Tsteam = tpsteam(Psteam);
        else
%         Saturation steam pressure (MPa); based on temperature
          Psteam = ptsteam(Tsteam);
        end
      end  

%     Dimensionless temperature 
      Tbar = Tsteam/Tcrit;  

%     Coefficients A(K)
      A(1) = 1.00;  
      A(2) = 6.978267;   
      A(3) = 2.599096; 
      A(4) = -0.998254;  

%     Intermediate calculation lambdao (dimensionless thermal conductivity
%     of steam in the ideal-gas limit)
      lambdao = 0.0; 
      for K = 1:4
        lambdao =  lambdao + A(K)/((Tbar)^(K-1));  
      end    
      lambdao = sqrt(Tbar)/lambdao;    

%     Coefficients B(I,J)
      B(1,1) = 1.3293046;  
      B(1,2) = -0.40452437;
      B(1,3) = 0.24409490;  
      B(1,4) = 0.018660751;  
      B(1,5) = -0.12961068;
      B(1,6) = 0.044809953;

      B(2,1) = 1.7018363;
      B(2,2) = -2.2156845;  
      B(2,3) = 1.6511057;
      B(2,4) = -0.76736002;
      B(2,5) = 0.37283344;
      B(2,6) = -0.11203160;

      B(3,1) = 5.2446158;  
      B(3,2) = -10.124111;  
      B(3,3) = 4.9874687;
      B(3,4) = -0.27297694;
      B(3,5) = -0.43083393;
      B(3,6) = 0.13333849;

      B(4,1) = 8.7127675;
      B(4,2) = -9.5000611;  
      B(4,3) = 4.3786606;  
      B(4,4) = -0.91783782;  
      B(4,5) = 0.0;  
      B(4,6) = 0.0;  

      B(5,1) = -1.8525999;
      B(5,2) = 0.93404690;
      B(5,3) = 0.0;  
      B(5,4) = 0.0;
      B(5,5) = 0.0; 
      B(5,6) = 0.0; 

%     Density of steam  (kg/m^3)  
      rhos = rhosteam(Tsteam, Psteam, xsteam); 

%     Critical steam density (kg/m^3)  
      rhocrit = 317.763;   

%     Dimensionless density
      rhobar = rhos/rhocrit;     

%     Intermediate lambda calculation
      lambda = 0.0;    
      for I = 1:5  
        for J = 1:6   
          lambda =  lambda + B(I,J)*((1.0D+00/Tbar - 1.0D+00)^(I-1))*                             ((rhobar - 1.0D+00)^(J-1));   
        end   
      end   
 
%     Intermediate calc
      lambda = exp(rhobar*lambda);    

%     Critical steam thermal conductivity (W/m-K)
      kcrit = 0.4945; 

      cond = lambdao*lambda*kcrit;   

%     Thermal conductivity (W/m-K) --> (kJ/h-m-K)  
      cond = cond*3.6;    

