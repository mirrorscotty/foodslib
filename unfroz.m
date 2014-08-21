      function [froz] = unfroz(Fcomp,Fctype,FtempK)

%  UNFROZ    Frozen and unfrozen water fractions of stream
%  
%            Synopsis: MATLAB function to determine the unfrozen and 
%            frozen water fractions of stream based on composition (w/w) 
%            and temperature. 
%  
%            Passed Arguments:
%             Fcomp - component mass fractions (w/w)
%             Fctype - component types
%             FtempK - temperature (K)
%  
%            Returned Argument:
%             froz - [unfrozen water fraction; frozen water fraction]
%  
%            Reference:  
%             Choi, Y & Okos, M.R. 1986.  Effects of Temperature and
%              Composition on the thermal properties of foods. In, "Food
%              Engineering Applications. Vol. 1," Elsevier Applied Science
%              Publishers, N.Y.

%  *************************************************************************

%     Reference temperature (K)
      Tref = tref;  

%     Latent heat of fusion (J/gmole)
      lambda = 6013.4;
 
%     Gas law constant (J/gmole-K)
      R = gaslawc/1000;

%     Intermediate term in determination of frozen water fraction
      iterm = lambda/R*((1/Tref)-(1/FtempK));
 
%     Initialize sum of mole fractions of solids 
      Fmw = 0;  

%     Loop for summing mole fractions of solids: sum(Xi/Mi) 
      for I = 1:length(Fcomp)  

        if (Fctype(I) == 1)
%         Note water component mass fraction (w/w)  
          Fh2o = Fcomp(I);           

        elseif (Fctype(I) > 1) 

%         Molecular weight of component (kg/kgmol)     
          mw = molwts(Fctype(I));

%         Mole fraction of component (kg-mol of component/kg of component) 
          Fcmf = Fcomp(I)/mw; 

%         Sum of mole fractions of solids  
          Fmw = Fmw + Fcmf; 

        end    
      end  

%     Effective molecular weight of solids (kg/kgmol) 
      MWeff = 1/Fmw; 


      if (FtempK < Tref)

%       Molecular weight of water (kg/kgmol)  
        mwh2o = molwts(1);  

%       Unfrozen water fraction
        Fh2ouf = mwh2o*Fmw*exp(iterm)/(1-exp(iterm));   

%       Some high concentrations yield an unfrozen water fraction 
%        greater than the total water concentration; in this case
%        let the unfrozen water fraction equal the total water fraction.
        if (Fh2ouf > Fh2o)
          Fh2ouf = Fh2o;
        end 

%       Frozen water fraction (w/w)
        Fh2of = Fh2o - Fh2ouf;

      else                         % FtempK > Tref 
        Fh2ouf = Fh2o;             % Unfrozen fraction is all water  
        Fh2of = 0;
      end

      froz = [Fh2ouf;Fh2of];

 
