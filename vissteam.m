      function [mu] = vissteam(Tsteam, Psteam, xsteam)

%  VISSTEAM  Dynamic viscosity of steam
%
%            Synopsis: MATLAB function which determines the steam dynamic
%            viscosity (Pa-s) based on temperature, pressure, and quality.
%
%            Call: vissteam(Tsteam, Psteam, xsteam)
%
%            Passed Variables:
%             Tsteam - temperature (K)
%             Psteam - pressure (MPa)
%             xsteam - quality
%
%            Returned Variable:
%             mu = steam dynamic viscosity (Pa-s)
%
%            Reference: Sengers, J.V. and Watson, J.T.R.  1985.  Improved
%            international formulation for the viscosity and thermal
%            conductivity of water substance.  J. Phys. Chem. Ref. Data
%            15(4): 1291-1314.
%            IAPS (International Association for the Properties of Steam)
%            Formulation 1985.

%  ***********************************************************************

%     For all steam dynamic viscosity determination:
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
      A(2) = 0.978197;
      A(3) = 0.579829;  
      A(4) = -0.202354;   

%     Intermediate calculation eta_o (dimensionless dynamic viscosity 
%     of steam in the ideal-gas limit)
      etao = 0.0;  
      for K = 1:4
        etao =  etao + A(K)/((Tbar)^(K-1));    
      end   
      etao = sqrt(Tbar)/etao;     

%     Coefficients B(I,J)
      B(1,1) = 0.5132047; 
      B(1,2) = 0.2151778;
      B(1,3) = -0.2818107;  
      B(1,4) = 0.1778064;
      B(1,5) = -0.04176610;
      B(1,6) = 0.0;
      B(1,7) = 0.0;

      B(2,1) = 0.3205656;
      B(2,2) = 0.7317883;  
      B(2,3) = -1.070786;
      B(2,4) = 0.460504;
      B(2,5) = 0.0;
      B(2,6) = -0.01578386;   
      B(2,7) = 0.0;

      B(3,1) = 0.0;  
      B(3,2) = 1.241044;  
      B(3,3) = -1.263184;
      B(3,4) = 0.234037;
      B(3,5) = 0.0;
      B(3,6) = 0.0;
      B(3,7) = 0.0;

      B(4,1) = 0.0;
      B(4,2) = 1.476783;  
      B(4,3) = 0.0;  
      B(4,4) = -0.4924179;  
      B(4,5) = 0.1600435;  
      B(4,6) = 0.0;  
      B(4,7) = -0.003629481;  

      B(5,1) = -0.7782567;
      B(5,2) = 0.0;
      B(5,3) = 0.0;  
      B(5,4) = 0.0;
      B(5,5) = 0.0; 
      B(5,6) = 0.0; 
      B(5,7) = 0.0; 

      B(6,1) = 0.1885447;  
      B(6,2) = 0.0;  
      B(6,3) = 0.0;  
      B(6,4) = 0.0;  
      B(6,5) = 0.0;
      B(6,6) = 0.0;
      B(6,7) = 0.0;

%     Density of steam  (kg/m^3)  
      rhos = rhosteam(Tsteam,Psteam,xsteam);   

%     Critical steam density (kg/m^3)  
      rhocrit = 317.763;   

%     Dimensionless density
      rhobar = rhos/rhocrit;    

%     Intermediate eta calculation
      eta = 0.0;    
      for I = 1:6
        for J = 1:7  
          eta =  eta + B(I,J)*((1.0/Tbar-1.0)^(I-1))*((rhobar-1.0)^(J-1)); 
        end    
      end     
  
%     Intermediate calc
      eta = exp(rhobar*eta);    

%     Critical steam viscosity (Pa-s)
      viscrit = 55.071e-6;    

      vissteam = etao*eta*viscrit;

      mu = vissteam;       
  


