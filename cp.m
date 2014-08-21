      function [Fcp] =  cp(Fcomp,Fctype,FtempK)  

%  CP        Specific heat of multi-component streams   
%  
%            Synopsis: MATLAB function to determine the specific heat a 
%            multi-component stream (kJ/kg-K).
%  
%            Call: cp(Fcomp,Fctype,FtempK)
%        
%            Passed Variables: 
%             Fcomp - component mass fractions (w/w) 
%             Fctype - component types 
%             FtempK - temperature (K) 
%  
%            Returned Variables:  
%             Fcp - specific heat (kJ/kg-K)  
%   
%
%            Reference:  
%             Choi, Y & Okos, M.R. 1986.  Effects of Temperature and
%               Composition on the thermal properties of foods. In, "Food
%               Engineering Applications. Vol. 1," Elsevier Applied Science
%               Publishers, N.Y. 

%            Note:  Authorized changes by M.R. Okos to water heat capacity 
%            constants as of 8/8/95.  
%
%  ************************************************************************

%     Limits of use?

%     Current number of component models for cp
      cptypes = 6;   

%     Reference temperature (K)
      Tref = 273.15;  

%     Convert T (K) to T (C) 
      TC = FtempK - Tref; 

      for I = 1:length(Fcomp) 

        if (Fctype(I) == comptype('water')) 

%         Specific heat of water (kJ/kg-K) 
          if (TC <= 0)  
            cph2o = 4.1280 - (5.3062e-3)*TC + (9.9516e-4)*(TC^2);   
          else 
            cph2o = 4.1280 - (9.0864e-5)*TC + (5.4731e-6)*(TC^2);  
          end  

%         Heat capacity of ice (kJ/kg-K)
          cpice = 2.0623 + 6.0769e-3*TC;

%         Unfrozen and frozen water fractions (w/w)  
          ice = unfroz(Fcomp,Fctype,FtempK); 
          Fh2ouf = ice(1);  % unfrozen 
          Fh2of = ice(2);   % frozen   

%         Weighted specific heat of water component (KJ/kg-K) 
          Fccp(I) = Fh2ouf*cph2o + Fh2of*cpice; 

        elseif (Fctype(I) == comptype('protein')) 

%         Specific heat of protein component (kJ/kg-K) 
          cp = 2.0082 + 0.0012089*TC - (1.3129e-6)*(TC^2);   
%         Weighted specific heat of protein compenent (kJ/kg-K)  
          Fccp(I) = Fcomp(I)*cp;   

        elseif (Fctype(I) == comptype('carbohydrate')) 

%         Specific heat of carbohydrate component (kJ/kg-K)  
          cp = 1.5488 + 0.0019625*TC - (5.9399e-6)*(TC^2);  
%         Weighted specific heat of carbohydrate compenent (kJ/kg-K)   
          Fccp(I) = Fcomp(I)*cp;

        elseif (Fctype(I) == comptype('fiber')) 

%         Specific heat of fiber component (kJ/kg-K)  
          cp = 1.8459 + 0.0018306*TC - (4.6509e-6)*(TC^2);   
%         Weighted specific heat of fiber compenent (kJ/kg-K)   
          Fccp(I) = Fcomp(I)*cp;

        elseif (Fctype(I) == comptype('fat')) 

%         Specific heat of fat component (kJ/kg-K)  
          cp = 1.9842 + 0.0014733*TC - (4.8008e-6)*(TC^2);   
%         Weighted specific heat of fat compenent (kJ/kg-K)   
          Fccp(I) = Fcomp(I)*cp;

        elseif (Fctype(I) == comptype('ash')) 

%         Specific heat of ash component (kJ/kg-K)  
          cp = 1.0926 + 0.0018896*TC - (3.6817e-6)*(TC^2);   
%         Weighted specific heat of ash compenent (kJ/kg-K)   
          Fccp(I) = Fcomp(I)*cp;

        end   

%       For blank and erroneous component types   
        if ((Fctype(I) <= 0) | (Fctype > cptypes))   
          Fccp(I) = 0;  
        end  
      end   

%     Sum of weighted specific heats (kJ/kg-K)  
      Fcp = sum(Fccp);       

 




