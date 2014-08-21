      function [Fmw] = molwts(Fctype)

%  MOLWTS    Molecular weights of stream components  
%  
%            Synopsis: MATLAB function which returns the molecular weight 
%            of stream components.  
%
%            Call: molwts(Fctype)   
%  
%            Passed Arguments:
%             Fctype - component types
%  
%            Returned Argument:
%             Fmw = molecular wight (kg/kgmol)
%  
%            Reference: (estimates)
%             Choi, Y & Okos, M.R. 1986.  Effects of Temperature and
%              Composition on the thermal properties of foods. In, "Food
%              Engineering Applications. Vol. 1," Elsevier Applied Science
%              Publishers, N.Y.

%  *************************************************************************

        if (Fctype == comptype('water'))
%         Molecular weight of water (kg/kgmol)  
          Fmw = 18;           

        elseif (Fctype == comptype('protein'))  
%         Molecular weight of protein (kg/kgmol)
          Fmw = 7.5e5;

        elseif (Fctype == comptype('carbohydrate'))   
%         Molecular weight of carbohydrate (kg/kgmol)
          Fmw = 342;

        elseif (Fctype == comptype('fiber'))   
%         Molecular weight of fiber (kg/kgmol)
          Fmw = 250000;

        elseif (Fctype == comptype('fat'))   
%         Molecular weight of fat (kg/kgmol)
          Fmw  = 450;

        elseif (Fctype == comptype('ash'))   
%         Molecular weight of ash (kg/kgmol)
          Fmw = 158;

        elseif (Fctype == comptype('air'))
%         Molecular weight of air (kg/kgmol)
          Fmw = 28.97;

        end 


 
