      function [Falpha] =  thdif(Fcomp,Fctype,FtempK)

%  THDIF     Thermal diffusivity of multi-component streams
%
%            Synopsis: MATLAB function to determine the thermal
%            diffusivity of a multi-component stream (m^2/h).
%
%            Call: thdif(Fcomp,Fctype,FtempK)
%
%            Passed Variables:
%             Fcomp - component mass fractions (w/w)
%             Fctype - component types
%             FtempK - temperature (K)
%
%            Returned Variables:
%             Falpha - thermal diffusivity (m^2/h)    
%
%            Reference:
%             Choi, Y & Okos, M.R. 1986.  Effects of Temperature and
%               Composition on the thermal properties of foods. In, "Food
%               Engineering Applications. Vol. 1," Elsevier Applied Science
%               Publishers, N.Y.

%
%  ************************************************************************

%     Reference temperature (K)
      Tref = 273.15;

%     Convert T (K) to T (C)
      TC = FtempK - Tref;  

      for I = 1:length(Fcomp)

        if (Fctype(I) == comptype('water'))

%         Thermal diffusivity of water component (m^2/s)  
          tdifh2o = 0.13168+(6.2477e-4*TC)-(2.4022e-6)*(TC^2);  
%         Density of water component (kg/m^3)
          rhoh2o = 997.18 + 0.0031439*TC - (0.0037574)*(TC^2);

%         Thermal diffusivity of ice (m^2/s)
          tdifice = 1.1756-(6.0833e-3*TC)+(9.5037e-5)*(TC^2);   
%         Density of ice fraction (kg/m^3)
          rhoice = 9.1689e2 - 1.3071e-1*TC;

%         Unfrozen and frozen water fractions (w/w)
          ice = unfroz(Fcomp,Fctype,FtempK);
          Fh2ouf = ice(1);  % unfrozen
          Fh2of = ice(2);   % frozen

%         Weighted thermal diffusivity of water component 
          Fctd(I) = tdifh2o*Fh2ouf/rhoh2o + tdifice*Fh2of/rhoice;  
      
        elseif (Fctype(I) == comptype('protein'))

%         Thermal diffusivity of protein component (m^2/s)  
          tdifpro = 6.8714e-2+(4.7578e-4*TC)-(1.4646e-6)*(TC^2);
%         Density of protein component (kg/m^3)
          rhopro = rho(1,Fctype(I),FtempK);
%         Weighted thermal diffusivity of protein component
          Fctd(I) = tdifpro*Fcomp(I)/rhopro; 

        elseif (Fctype(I) == comptype('carbohydrate'))

%         Thermal diffusivity of carbohydrate component (m^2/s)  
          tdifcarb = 8.0842e-2+(5.3052e-4*TC)-(2.3218e-6)*(TC^2);  
%         Density of carbohydrate component (kg/m^3) 
          rhocarb = rho(1,Fctype(I),FtempK);
%         Weighted thermal diffusivity of carbohydrate component 
          Fctd(I) = tdifcarb*Fcomp(I)/rhocarb; 

        elseif (Fctype(I) == comptype('fiber')) 

%         Thermal diffusivity of fiber component (m^2/s) 
          tdiffib = 7.3976e-2+(5.1902e-4*TC)-(2.2202e-6)*(TC^2);
%         Density of fiber component (kg/m^3)  
          rhofib = rho(1,Fctype(I),FtempK); 
%         Weighted thermal diffusivity of fiber component 
          Fctd(I) = tdiffib*Fcomp(I)/rhofib;  

        elseif (Fctype(I) == comptype('fat'))

%         Thermal diffusivity of fat component (m^2/s)  
          tdiffat =  9.8777e-2-(1.2569e-4*TC)-(3.8286e-8)*(TC^2);   
%         Density of fat component (kg/m^3)   
          rhofat = rho(1,Fctype(I),FtempK); 
%         Weighted thermal diffusivity of fiber component  
          Fctd(I) = tdiffat*Fcomp(I)/rhofat;   

        elseif (Fctype(I) == comptype('ash'))

%         Thermal diffusivity of ash component (m^2/s)
          tdifash = 1.2461e-1+(3.7321e-4*TC)-(1.2244e-6)*(TC^2);  
%         Density of ash component (kg/m^3)    
          rhoash = rho(1,Fctype(I),FtempK);  
%         Weighted thermal diffusivity of fiber component  
          Fctd(I) = tdifash*Fcomp(I)/rhoash;    

        end 
      end   

%     Density of stream (kg/m^3)
      Frho = rho(Fcomp,Fctype,FtempK); 

%     Thermal diffusivity (m^2/s)  
      tdtot=sum(Fctd)*Frho;  

%     Thermal diffusivity (m^2/s) --> (m^2/h)  
      Falpha  = tdtot*3600*10e-7;    


