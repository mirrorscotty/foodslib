      function [density] = rhosteam(Tsteam, Psteam, xsteam)

%  RHOSTEAM  Density of steam
%
%            Synopsis: MATLAB function which determines the steam
%            density (kg/m^3) based on temperature, pressure, and quality.
%
%            Call: rhosteam(Tsteam, Psteam, xsteam)
%
%            Passed Variables:
%             Tsteam - temperature (K)
%             Psteam - pressure (MPa)
%             xsteam - quality
%
%            Returned Variable:
%             density = steam density (kg/m^3)
%
%            Reference:  Irvine, T.F., and Liley, P.E.   1984.  "Steam
%             and Gas Tables with Computer Equations,"  Academic Press,
%             Inc., N.Y.

%  ***********************************************************************

%     For mixed water vapor steam specific volume determination:
      if (xsteam >= 0.0) & (xsteam <= 1.0)   
%       Critical temperature (K) 
        Tcrit = 647.3;    

%       If pressure of steam is given   
        if (Psteam > 0)  
%         Saturation steam temperature (K); based on pressure 
          Tsteam = tpsteam(Psteam);  
        else 
%         Saturation steam pressure (MPa); based on temperature 
          Psteam = ptsteam(Tsteam);    
        end      
 
%       Critical region effect
        Tc = (Tcrit - Tsteam)/Tcrit;   

%       Specific volume of saturated liquid (m^3/kg)

        if (Tsteam >= 273.16) & (Tsteam <= 647.3) 
          A = 1.0;
          B = -1.9153882;  
          C = 1.2015186e1;     
          D = -7.8464025;   
          E(1) = -3.888614;    
          E(2) = 2.0582238;       
          E(3) = -2.0829991;   
          E(4) = 8.2180004e-1;   
          E(5) = 4.7549742e-1;    
          E(6) = 0.0;    
          E(7) = 0.0;     
        end   
   
        Vfcr = 3.155e-3;       
        Esum = 0.0; 
        for N = 1:7  
           Esum = Esum+E(N)*(Tc^N);   
        end    

%       Exponents for Ys equation 
        F = 1.0/3.0;   
        G = 5.0/6.0;
        H = 7.0/8.0; 
        Ys = A + B*(Tc^(F)) + C*(Tc^(G)) + D*(Tc^(H)) + Esum;   

        Vsatliq = Ys*Vfcr;   

%       Specific volume of saturated steam (m^3/kg)
      
        A = 1.0;  
        B = 1.6351057;  
        C = 5.2584599e+1;    
        D = -4.4694653e+1;   
        E(1) = -8.9751114;  
        E(2) = -4.3845530e-1;  
        E(3) = -1.9179576e+1;    
        E(4) = 3.6765319e+1;    
        E(5) = -1.9462437e+1;    
        E(6) = 0.0;  
        E(7) = 0.0;   

        Pcr = 2.2089e+1;    
        Esum = 0.0;   
        for N = 1:7  
          Esum = Esum+E(N)*(Tc^N);    
        end     
        Ys = A + B*(Tc^(F)) + C*(Tc^(G)) + D*(Tc^(H)) + Esum;  

        Vsatgas = Ys*Vfcr*Pcr/Psteam;      
        
%       Specific volume of steam based on quality (m^3/kg)

        vsteam = Vsatliq + xsteam*(Vsatgas-Vsatliq);  

      end    

%     For superheated steam specific volume (m^3/kg)   
      if (xsteam > 1.0) 

%       Saturation temperature for the pressure specified
        Tsat = tpsteam(Psteam);    

        R = 4.61631e-4;    

        BSH(1) = 5.27993e-2;    
        BSH(2) = 3.75928e-3;    
        BSH(3) = 2.2e-2;    

        ASH(1) = -3.741378;    
        ASH(2) = -4.7838281e-3;    
        ASH(3) = 1.5923434e-5;    
      
        exsum = 0.0;  
        vsteam = 0.0;
        for N = 1:3
          exsum =  exsum + ASH(N)*(Tsat^(N-1));        
        end  

        M = 4.0e+1;    

        vsteam = R*Tsteam/Psteam;    
        vsteam = vsteam - BSH(1)*exp(-BSH(2)*Tsteam);    
        vsteam = vsteam+(1/(10*Psteam))*(BSH(3)-exp(exsum))*exp((Tsat-Tsteam)/M);     

      end    

%     Density of steam (kg/m^3)   
      density = 1.0/vsteam;    
  

