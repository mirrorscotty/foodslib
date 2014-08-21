      function [Fk] =  thermc(Fcomp,Fctype,FtempK)

%  THERMC    Thermal conductivity of multi-component streams
%
%            Synopsis: MATLAB function to determine the thermal
%            conductivity of a multi-component stream (kJ/h-m-K).
%
%            Call: thermc(Fcomp,Fctype,FtempK)
%
%            Passed Variables:
%             Fcomp - component mass fractions (w/w)
%             Fctype - component types
%             FtempK - temperature (K)
%
%            Returned Variables:
%             Fk - thermal conductivity (kJ/h-m-K)    
%
%            Reference:
%             Choi, Y & Okos, M.R. 1986.  Effects of Temperature and
%               Composition on the thermal properties of foods. In, "Food
%               Engineering Applications. Vol. 1," Elsevier Applied Science
%               Publishers, N.Y.

%            Note:  Authorized changes by M.R. Okos to fat thermal conductivity
%            coefficient on TC as of 8/8/95.
%
%  ************************************************************************

%     Reference temperature (K)
      Tref = 273.15;

%     Convert T (K) to T (C)
      TC = FtempK - Tref;  

      for I = 1:length(Fcomp)

        if (Fctype(I) == comptype('water'))

%         Thermal conductivity of water component (W/m-K)  
          tconh2o = 0.57109+(1.7625e-3*TC)-(6.7036e-6)*(TC^2);  
%         Density of water component (kg/m^3)
          rhoh2o = 997.18 + 0.0031439*TC - (0.0037574)*(TC^2);

%         Thermal conductivity of ice (W/m-K)
%         tconice = 2.2196-(6.2489e-3*TC)-(1.0154e-4)*(TC^2);   
          tconice = 2.2196-(6.2489e-3*TC)+(1.0154e-4)*(TC^2);   
%         Density of ice fraction (kg/m^3)
          rhoice = 9.1689e2 - 1.3071e-1*TC;

%         Unfrozen and frozen water fractions (w/w)
          ice = unfroz(Fcomp,Fctype,FtempK);
          Fh2ouf = ice(1);  % unfrozen
          Fh2of = ice(2);   % frozen

%         Weighted thermal conductivity of water component 
          Fctc(I) = tconh2o*Fh2ouf/rhoh2o + tconice*Fh2of/rhoice;  
      
        elseif (Fctype(I) == comptype('protein'))

%         Thermal conductivity of protein component (W/m-K)  
          tconpro = 0.17881+(1.1958e-3*TC)-(2.7178e-6)*(TC^2);
%         Density of protein component (kg/m^3)
          rhopro = rho(1,Fctype(I),FtempK);
%         Weighted thermal conductivity of protein component
          Fctc(I) = tconpro*Fcomp(I)/rhopro; 

        elseif (Fctype(I) == comptype('carbohydrate'))

%         Thermal conductivity of carbohydrate component (W/m-K)  
          tconcarb = 0.20141+(1.3874e-3*TC)-(4.3312e-6)*(TC^2);  
%         Density of carbohydrate component (kg/m^3) 
          rhocarb = rho(1,Fctype(I),FtempK);
%         Weighted thermal conductivity of carbohydrate component 
          Fctc(I) = tconcarb*Fcomp(I)/rhocarb; 

        elseif (Fctype(I) == comptype('fiber')) 

%         Thermal conductivity of fiber component (W/m-K) 
          tconfib = 0.18331+(1.2497e-3*TC)-(3.1683e-6)*(TC^2);
%         Density of fiber component (kg/m^3)  
          rhofib = rho(1,Fctype(I),FtempK); 
%         Weighted thermal conductivity of fiber component 
          Fctc(I) = tconfib*Fcomp(I)/rhofib;  

        elseif (Fctype(I) == comptype('fat'))

%         Thermal conductivity of fat component (W/m-K)  
          tconfat = 0.18071-(2.7604e-4*TC)-(1.7749e-7)*(TC^2);   
%         Density of fat component (kg/m^3)   
          rhofat = rho(1,Fctype(I),FtempK); 
%         Weighted thermal conductivity of fiber component  
          Fctc(I) = tconfat*Fcomp(I)/rhofat;   

        elseif (Fctype(I) == comptype('ash'))

%         Thermal conductivity of ash component (W/m-K)  
          tconash = 0.32962+(1.4011e-3*TC)-(2.9069e-6)*(TC^2);  
%         Density of ash component (kg/m^3)    
          rhoash = rho(1,Fctype(I),FtempK);  
%         Weighted thermal conductivity of fiber component  
          Fctc(I) = tconash*Fcomp(I)/rhoash;    

        end 
      end   

%     Density of steam (kg/m^3)
      Frho = rho(Fcomp,Fctype,FtempK); 

%     Thermal conductivity (W/m-K)  
      tctot=sum(Fctc)*Frho;  

%     Thermal conductivity (W/m-K) --> (kJ/h-m-K)  
      Fk  = tctot*3.6;    


