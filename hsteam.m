      function [enthalpy] = hsteam(Tsteam, Psteam, xsteam) 

%  HSTEAM    Enthalpy of steam   
%   
%            Synopsis: MATLAB function which determines the steam 
%            enthalpy (kJ/kg) based on temperature, pressure, and quality. 
%  
%            Call: hsteam(Tsteam, Psteam, xsteam)
%  
%            Passed Variables:
%             Tsteam - temperature (K) 
%             Psteam - pressure (MPa)  
%             xsteam - quality   
%   
%            Returned Variable:  
%             enthalpy = steam enthalpy (kJ/kg)    
%    
%            Reference:  Irvine, T.F., and Liley, P.E.   1984.  "Steam 
%             and Gas Tables with Computer Equations,"  Academic Press, 
%             Inc., N.Y.

%  ***********************************************************************

%     For saturated steam enthalpy determination:
      if xsteam >= 0 & xsteam <= 1  
%       Critical temperature (K) 
        Tcrit = 647.3;   

%       If pressure of steam is given   
        if (Psteam > 0) & (Tsteam == 0)  

%         Saturation steam temperature (K); based on pressure 
          Tsteam = tpsteam(Psteam);   
        end     
 
%       Critical region effect
        Tc = (Tcrit - Tsteam)/Tcrit;   

%       Enthalpy of saturated liquid at specified pressure (KJ/kg)

        if (Tsteam >= 273.16) & (Tsteam < 300) 
          A = 0.0;  
          B = 0.0;  
          C = 0.0;   
          D = 0.0;  
          E(1) = 6.24698837e2;   
          E(2) = -2.34385369e3;    
          E(3) = -9.50812101e3;  
          E(4) = 7.16287928e4;  
          E(5) = -1.63535221e5;   
          E(6) = 1.66531093e5;   
          E(7) = -6.47854585e4;   

        elseif (Tsteam >= 300) & (Tsteam < 600)  
          A = 8.839230108e-01;   
          B = 0.0;
          C = 0.0;
          D = 0.0;
          E(1) = -2.67172935;   
          E(2) = 6.22640035;     
          E(3) = -1.31789573e1;   
          E(4) = -1.91322436;   
          E(5) = 6.87937653e1;  
          E(6) = -1.24819906e2;   
          E(7) = 7.21435404e1;    

        elseif (Tsteam >= 600) & (Tsteam <= 647.3)  
          A = 1.0;    
          B = -4.41057805e-01;    
          C = -5.52255517;    
          D = 6.43994847;   
          E(1) = -1.64578795;    
          E(2) = -1.30574143;      
          E(3) = 0.0;   
          E(4) = 0.0;     
          E(5) = 0.0;   
          E(6) = 0.0;   
          E(7) = 0.0;  

        end     
   
        Hfcr = 2.0993e3;   
        Esum = 0.0; 
        for N = 1:7  
          Esum = Esum+E(N)*(Tc^N);   
        end   

%       Exponents for Ys equation 
        F = 1.0/3.0; 
        G = 5.0/6.0; 
        H = 7.0/8.0; 
        Ys = A + B*(Tc^(F)) + C*(Tc^(G)) + D*(Tc^(H)) + Esum;  

        Hsatliq = Ys*Hfcr;    

%       Enthalpy of saturated gas at specified pressure (KJ/kg)
      
        A = 1.0D+00;    
        B = 4.57874342e-1;   
        C = 5.08441288;    
        D = -1.48513244;     
        E(1) = -4.81351884;    
        E(2) = 2.69411792;    
        E(3) = -7.39064542;    
        E(4) = 1.04961689e1;  
        E(5) = -5.46840036;    
        E(6) = 0.0;     
        E(7) = 0.0;    
 
        Esum = 0.0;   
        for N = 1:7  
          Esum = Esum+E(N)*(Tc^N);   
        end
        Ys = A + B*(Tc^(F)) + C*(Tc^(G)) + D*(Tc^(H)) + Esum;   

        Hsatgas = Ys*Hfcr;  
        
%       Enthalpy of steam based on quality (kJ/kg)
        enth = Hsatliq + xsteam*(Hsatgas-Hsatliq);   
        enthalpy = enth; 
      end    

%     For superheated steam enthalpy
      if (xsteam > 1)  

%       Saturation temperature for the pressure specified
        Tsat = tpsteam(Psteam);  

        BSH(1) = 2.04121e3;    
        BSH(2) = -4.040021e1; 
        BSH(3) = -4.8095e-1;    

        ASH(1) = BSH(1) + BSH(2)*Psteam + BSH(3)*(Psteam^2);   
  
        BSH(4) = 1.610693;     
        BSH(5) = 5.472051e-2;   
        BSH(6) = 7.517537e-4;   

        ASH(2) = BSH(4) + BSH(5)*Psteam + BSH(6)*(Psteam^2);    
 
        BSH(7) = 3.383117e-4;   
        BSH(8) = -1.975736e-5;     
        BSH(9) = -2.87409e-7;   

        ASH(3) = BSH(7) + BSH(8)*Psteam + BSH(9)*(Psteam^2);      
 
        BSH(10) = 1.70782e3;   
        BSH(11) = -1.699419e1;   
        BSH(12) = 6.2746295e-2;      
        BSH(13) = -1.0284259e-4;       
        BSH(14) = 6.4561298e-8;      

        ASH(4) = BSH(10) + BSH(11)*Tsat + BSH(12)*(Tsat^2); 
        ASH(4) = ASH(4) + BSH(13)*(Tsat^3) + BSH(14)*(Tsat^4);  
       
        enth = 0.0;  
        for N = 1:3
          enth = enth + ASH(N)*(Tsteam^(N-1));    
        end 
 
        M = 4.5e1;    

        enth = enth - (ASH(4)*exp((Tsat-Tsteam)/M));   
        enthalpy = enth;  

      end  
  

