      function [Psat] = ptsteam(Tsteam)

%  TPSTEAM   Saturation temperature given steam pressure
%
%            Synopsis: MATLAB function which determines the steam
%            saturation pressure (MPa) based on steam temperature.
%
%            Call: ptsteam(Tsteam)
%
%            Passed Variables:
%             Tsteam - temperature (K)
%
%            Returned Variable:
%             Psat - pressure (MPa)
%
%            Reference:  Irvine, T.F., and Liley, P.E.   1984.  "Steam
%             and Gas Tables with Computer Equations,"  Academic Press,
%             Inc., N.Y.

%  ***********************************************************************

%     Requires check for 273.16 <= Tsteam <= 647.3 K 

      if (Tsteam >= 273.16) & (Tsteam <= 647.3)    
        A(1) = 0.104592e2;   
        A(2) = -0.404897e-2; 
        A(3) = -0.417520e-4;  
        A(4) = 0.368510e-6;  
        A(5) = -0.101520e-8;  
        A(6) = 0.865310e-12;   
        A(7) = 0.9036680e-15;  
        A(8) = -0.199690e-17;  
        A(9) = 0.779287e-21; 
        A(10) = 0.191482e-24;     
        A(11) = -0.396806e4;  
        A(12) = 0.395735e2;       

      end  

      intermed1 = 0.0;     
      for N = 1:10 
        intermed1 = intermed1 + (A(N)*(Tsteam^(N-1)));    
      end  

      intermed2 = intermed1 + (A(11)/(Tsteam-A(12)));   
      Psat = exp(intermed2);   
 

