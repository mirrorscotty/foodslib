      function [Tsat] = tpsteam(Psteam)

%  TPSTEAM   Saturation temperature given steam pressure
%
%            Synopsis: MATLAB function which determines the steam
%            saturation temperature (K) based on steam pressure.
%
%            Call: tpsteam(Psteam)
%
%            Passed Variables:
%             Psteam - pressure (MPa)
%
%            Returned Variable:
%             Tsat - Temperature (K)
%
%            Reference:  Irvine, T.F., and Liley, P.E.   1984.  "Steam
%             and Gas Tables with Computer Equations,"  Academic Press,
%             Inc., N.Y.

%  ***********************************************************************

%     Requires check at interface for 0.000611 <= Psteam <= 22.1 MPa

      if (Psteam >= 0.000611) & (Psteam < 12.33)  
        A = 0.426776e2;  
        B = -0.389270e4;   
        C = -0.948654e1;  

      elseif (Psteam >= 12.33) & (Psteam <= 22.1)  
        A = -0.387592e3;  
        B = -0.125875e5;  
        C = -0.152578e2;  
      end    

      Tsat = A + B/(log(Psteam) + C);     




