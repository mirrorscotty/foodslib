       function [mu] = viscos(Fcomp,Fctype,FtempK,vmodel)

%  VISCOS    Viscosity models
%
%            Synopsis: MATLAB function to determine viscosity of 
%            multi-component streams.  
%
%            Call: viscos(Fcomp,Fctype,FtempK,vmodel) 
%
%            Passed Arguments: 
%             Fcomp - component mass fractions (w/w)
%             Fctype - component types
%             FtempK - temperature (K)
%             vmodel - viscosity model  
%
%            Returned Arguments:
%             mu(1) - Consistency coefficient (Pa-s)
%             mu(2) - Flow behavior index    
%
%            References: 
%            [1] Steffe, J.F., Mohammed, I.O., and Ford, E.W. 1986. 
%                 Rheological properties of fluid foods: data compilation.
%                 In: "Physical and Chemical Properties of Foods," 
%                 M.R. Okos (Ed.), ASAE St Joseph, Mich. p 1-13.    
%            [2] Dale, M.C., Okos, M.R., and Nelson,P. 1982. Concentration of 
%                 Tomato Products: Analysis of Energy Savings Process
%                 Alternatives. J. Food Sci. 47(1853-1858).
%            [3] Angeletti, S. and Moresi, M. 1983. Modelling of Multiple-
%                 effect Falling-film Evaporators. J. Fd. Technol. 18(539-563)  
%            [4] Moresi, M. and Spinosi, M. 1984. Physical Properties of 
%                 Concentrated Apple Juice. In: "Engineering and Food. Vol. 1 - 
%                 Engineering Sciences in the Food Industry," B.M. McKenna
%                 (Ed.), Elsevier Applied Science Publishers, London.    


%  *************************************************************************

%       Reference temperature (K)
        Tref = tref;  

        if (vmodel == 1) 

%         Model based on total solids, protein, and fiber [Ref. 2]
%         Model actually determines viscosity not K and n 

%         Total solids (w/w) 
          Fts = tsolids(Fcomp,Fctype); 

%         Protein and fiber mass fractions (w/w)          
          Fprot = findcomp(Fcomp,Fctype,comptype('protein')); 
          Ffiber = findcomp(Fcomp,Fctype,comptype('fiber'));   

%         Consistency coefficient (Pa-s) 
          K = (1.44e-6*exp(1928/FtempK))*(1+Fprot*18)*(1+Fts*0.295)* ...
               (1+Ffiber*333.0);     
%         Flow behavior index
          n = 1;   


        elseif (vmodel == 2)  

%         Raw milk [Ref. 1]  
%         Temperature (C) --> (K)  
          TC = [0 5 10 20 25 30 35 40]; 
          TK = Tref + TC;

%         Power law index; Newtonian fluid; so mu = K 
          n = 1;
%         Consistency coefficient (Pa-s)  
          K = [0.00344 0.00305 0.00264 0.00199 0.0017 0.00149 0.00134 0.00123];

%         Linear interpolate in 1D-look-up table to find viscosity at specified
%         temperature 
          if (FtempK > max(TK))            % Specified temperature out of range
            K = interp1(TK,K,max(TK));  
          elseif (FtempK < min(TK));  
            K = interp1(TK,K,min(TK));  
          else 
            K = interp1(TK,K,FtempK);  
          end   

        elseif (vmodel == 3)   
 
%         Homogenized milk [Ref. 1]  
%         Temperature (C) --> (K)   
          TC = 20:10:80;  
          TK = Tref + TC;  
 
%         Power law index; Newtonian fluid; so mu = K 
          n = 1;  
%         Consistency coefficient (Pa-s) 
          K = [0.002 0.0015 0.0011 0.00095 0.000775 0.0007 0.0006]; 

%         Linear interpolate in 1D-look-up table to find viscosity 
          if (FtempK > max(TK))            % Specified temperature out of range 
            K = interp1(TK,K,max(TK)); 
          elseif (FtempK < min(TK)); 
            K = interp1(TK,K,min(TK));  
          else 
            K = interp1(TK,K,FtempK);      % At specified temperature  
          end  

        elseif (vmodel == 4) 

%         Apple Juice [Ref. 4] Good for 0<Fraction Sucrose Equivalent<0.748  
%         Total solids (w/w)  
          Fts = tsolids(Fcomp,Fctype); 
%         Mass solids to mass water ratio 
          Ftsw = Fts/(1-Fts);
%         Flow behavior index 
          n = 1 - 0.186*Ftsw + 0.047*Ftsw^2;  
%         Frequency Factor 
          C = -13.74 + 4.88*Ftsw - 2.38*Ftsw^2
%         Constant proportional to activation energy of flow 
          D = 1926.48 + 588*Ftsw - 873*Ftsw^2 + 700*Ftsw^3 - 106.7*Ftsw^4
%         Consistency coefficient (Pa-s^n)
          K = exp(C + D/FtempK);

        elseif (vmodel == 5)

%         Orange Juice [Ref. 3]  Good for 25 < FtempC < 50 
  
%         Total solids (w/w)
          Fts = tsolids(Fcomp,Fctype);
%         Regressions coefficient alpha (cp-K^Beta)  
          alpha = exp(34.668 - 20.24*Fts + 162*Fts^2) 
%         Regression coefficient Beta 
          Beta = -6.1055 + 3.96*Fts - 26.8*Fts^2 
%         Viscosity (cp --> Pa-s); taken as consistency coefficient 
          K = alpha * (FtempK^Beta) /1000
%         flow behavior index
          n = 1; 

        elseif (vmodel == 6)  

%         Water 
          K = vissteam(FtempK,0,0);
          n = 1;

        elseif (vmodel == 7) 
 
%         Air 
          K = muair(FtempK);
          n = 1;

        end   

        mu = [K;n];  
 
