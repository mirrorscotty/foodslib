      function [Fbpr] =  bpr(Fcomp,Fctype,WtempK,bprmodel)

%  BPR       Boiling point rise of multi-component stream   
%
%            Synopsis: MATLAB function to determine the boiling point rise 
%            of a multi-component stream (K).
%
%            Call: bpr(Fmass,model)
%
%            Passed Variables:
%             Fcomp - component mass fractions (w/w)
%             Fctype - component types
%             WtempK - boiling temperature of pure water (K) 
%             bprmodel - boiling point rise model   
%
%            Returned Variables:
%             Fbpr - boiling point rise (K)
%
%            Reference:
%            [1] Hui, Y (Ed.) Encyclopedia of Food Science & Technology. 
%            [2] Hartel, R.W. 1992. Evaporation and Freeze Concentration
%                  In "Handbook of Food Engineering"; D.R> Heldman & D.B.
%                  Lund (Eds.), Marcel Dekker, Inc. NY
%            [3] Geankoplis, C.J. 1983.  "Transport Processes and Unit
%                  Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.

%  ************************************************************************

%     Reference temperature (K)
      Tref = tref;

%     For no specified boiling point rise model 
      if (bprmodel == 0)     
        Fbpr = 0;

%     Model using molar concentrations [Refs. 1 & 2] 
      elseif (bprmodel == 1) 

%       Gas law constant (kJ/kgmol-K) 
        R = gaslawc/1000;

%       Enthalpy of saturated vapor (kJ/kg) 
        Hvap = hsteam(WtempK,0,1); 
%       Enthalpy of saturated liquid (kJ/kg) 
        hliq = hsteam(WtempK,0,0); 
%       Latent heat of vaporization (kJ/kg) 
        latent = Hvap - hliq;

        for I = 1:length(Fcomp)  

%         Molecular weight of component (kg/kgmol) 
          mw = molwts(Fctype(I));

%         Moles of component (kgmol) based on 1 kg total 
          moles(I) = Fcomp(I)/mw; 

        end  

%       Molar concentration of each component (kgmol/kgmol-total)
        molconc = moles./sum(moles);

        for I = 1:length(Fcomp);  

          if (Fctype(I) ~= comptype('water')) 

%           Boiling point rise for each non water component (K) 
            delT(I) = (R * (WtempK-Tref)^2 * molconc(I)) / latent;

          end 
        end  

%       Total boiling point rise (K) 
        Fbpr = sum(delT);

%     Model for sugar solution [Ref. 3] 
      elseif (bprmodel == 2)

%       Assumes all solids are sugar (w/w)    
        Fts = tsolids(Fcomp,Fctype);

%       Boiling point rise (K) 
        Fbpr = 1.78*Fts + 6.22*Fts^2;  
        
%     Model for milk [Ref. 1]
      elseif (bprmodel == 3)  

%       Milk solids data 
        Mts = [0.16; 0.275; 0.39; 0.49; 0.62; 0.69; 0.73;];

%       Boiling point Rise data (K)
        Mbpr = [0.5; 1; 1.5; 2; 3; 4; 5;];

%       Actual stream solids (w/w)   
        Fts = tsolids(Fcomp,Fctype);

%       Use solids and BPR data in 1-D look-up table
        if (Fts < min(Mts)) 
          Fbpr = interp1(Mts,Mbpr,min(Mts));  
        elseif (Fts > max(Mts))
          Fbpr = interp1(Mts,Mbpr,max(Mts));  
        else 
          Fbpr = interp1(Mts,Mbpr,Fts);  
        end   

      end   % BPR model selectio
        


