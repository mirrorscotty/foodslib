      function [properties]= airwater(prop,Tdb,press,method,state)

%  AIRWATER  Air-water (Psychrometric) Properties
%
%            Synopsis: MATLAB function to determine psychrometric properties 
%            for an air-water vapor mixtures (-100 to 200 C) where two 
%            psychrometric variables are specified.   
% 
%            Call: airwater(prop,Tdb,press,method,state) 
%   
%            Passed Variables:  
%             prop - variable determines which property value is returned:
%                    1. Enthalpy (kJ/kg dry air)
%                    2. Absolute humidity (kg water/kg dry air)
%                    3. Specific Volume (m^3/kg dry air)
%                    4. Relative Humidity (%)
%                    5. Wet bulb temperature (K)
%                    6. Dew Point temperature (K)
%             Tdb - dry bulb temeprature (K) 
%             press - pressure (MPa)  
%             method - variable determines which method is used to find prop:  
%                    1. Dry bulb temp and relative humidity ("state" in %)
%                    2. Dry bulb temp and wet bulb temperature (K)
%                    3. Dry bulb temp and dew point temperature (K)
%                    4. Dry bulb temp and absolute humidity (kg h2o/kg dry air)
%             state - value of second variable needed for method 
%  
%             CAUTION!!!: If algorithm is returning exceptionally large 
%             values or negative values, the function is being passed state 
%             values that indicate that the air is holding more water than 
%             is possible. 
%
%             Reference:  
%              ASHRAE Fundamentals SI Edition.
%

%  ***************************************************************************
 
%     Reference temperature (K)
      Tref = tref;  
%     Dry bulb temperature (C)
      TdbC = Tdb - Tref;   
%     Inititalize 
      Twbset = 0; 

%     Error message
      if (TdbC > 200) | (TdbC < -100)    
        errortxt = [10 'AIRWATER MATLAB FUNCTION ERROR: Temperature range for psychrometric properties' 10 ...  
        'is limited to -100 to 200 C.' 10];    
      end   

%     Total pressure (MPa) -> (kPa)
      p = press*1000;

%     Constants for determination of water vapor saturation pressure (Pa)
%     For temperature range of -100 to 0 C.
      C1 = -5674.5359;
      C2 = 6.3925247;
%     Value given by ASHRAE Funadamentals 1993 C(2)=-5.1523058e-01
      C3 = -0.9677843e-2;
      C4 = 0.62215701e-6;
      C5 = -0.20747825e-7;
      C6 = -0.9484024e-12;
      C7 = 4.1635019;

%     For temperature range of 0 to 200 C.
      C8 = -5.8002206e3;
      C9 = 1.3914993;
%     Value given by ASHRAE Fundamentals 1993 C(9)=-5.5162560
      C10 = -4.8640239e-2;
      C11 = 4.1764768e-5;
      C12 = -1.4452093e-8;
      C13 = 6.5459673;   

%     Intermediate saturation pressure calc
      if (TdbC <= 0.0) & (TdbC > -100)
        pwsint = C1/Tdb+C2+C3*Tdb+(C4*Tdb^2)+(C5*Tdb^3)+(C6*Tdb^4)+C7*log(Tdb);
      elseif (TdbC > 0.0) & (TdbC <= 200.0)
        pwsint = C8/Tdb+C9+C10*Tdb+(C11*Tdb^2)+(C12*Tdb^3)+C13*log(Tdb);  
      end  
%     Actual saturation pressure (kPa)
      if (TdbC >= -100) & (TdbC <= 200.0)
        pws = exp(pwsint)/1000;  
      else 
        pws = ptsteam(Tdb)*1000;   
      end  

%     Sat. humidity ratio at dry bulb temp (kg water/kg dry air)
      ws = 0.62198*pws/(p-pws);
 
%     Gas law constant divided by Mol. Wt. of air
      Ra = 8.31434/28.9645;
 
%     For properties based on dry bulb temperature and relative humidity
      if method==1 
%       Percent relative humidity
        RH = state; 
%       Partial pressure of water vapor (Pa)
        pw = RH*pws/100; 
%       Humidity ratio (kg water/kg dry air)
        w = 0.62198*pw/(p-pw); 
      end 
%
%     For properties based on dry bulb temperature and wet bulb temperature
      if method==2  
        Twbset = 1;   
%       Wet bulb temperature (K)
        Twb = state;
%       Wet bulb temperature (C)
        TwbC = Twb - Tref;
 
%       Water saturation pressure at wet bulb temperature
%       Intermediate calc.
        pwswbint = C8/Twb+C9+C10*Twb+(C11*Twb^2)+(C12*Twb^3)+C13*log(Twb);
%       Actual saturation pressure (Pa)
        pwswb = exp(pwswbint)/1000;
 
%       Humidity ratio of saturated air (kg water/kg dry air)
        wswb = 0.62198*(pwswb/(p-pwswb));
%       Humidity ratio (kg water/kg dry air)
        w = ((2501.0-2.381*TwbC)*wswb-(TdbC-TwbC))/(2501.0+1.805*TdbC-4.186*TwbC);
%       Partial pressure of water vapor (Pa)  
        pw = p*w/(0.61298 + w);  
%       Degree of saturation
        mu = w/ws; 
%       Relative humidity (%)
        RH = mu*100/(1.0-(1.0-mu)*(pws/p)); 

     end 
 
%     For properties based on dry bulb temperature and dew point temperature
      if method==3 
        Tdp = state;
        Tdpset = 1;   
%       Water saturation pressure at dew point temperature
%       Intermediate calc.
        pwsdpint = C8/Tdp+C9+C10*Tdp+(C11*Tdp^2)+(C12*Tdp^3)+C13*log(Tdp);
        pwsdp = exp(pwsdpint)/1000;
%       Partial pressure of water vapor (Pa)  
        pw = pwsdp;  
%       Humidity ratio (kg water/kg dry air)
        w = 0.62198*pwsdp/(p-pwsdp);
%       Degree of saturation
        mu = w/ws;
%       Relative humidity (%)
        RH = mu*100/(1.0-(1.0-mu)*(pws/p));
      end

%     For properties based on dry bulb temperature and absolute humidity
      if method==4 
        w = state;

%       Partial pressure of water vapor (Pa)
        pw = p*w/(0.62198+w);
%       Relative humidity (%)
        RH = (pw/pws)*100;
      end

%     Enthalpy (kJ/kg dry air)
      H = (1.006*TdbC)+w*(2500.45+1.805*TdbC);  

%     Return specified property value
      if prop==1  
        awprop = H;
      elseif prop==2 
        awprop = w;
      elseif prop==3
%       Specific volume (m^3/kg dry air)
        v = Ra*Tdb/(p-pw);
        awprop = v;
      elseif prop==4 
        awprop = RH;
      elseif prop==5 
%       Wet bulb temperature (K)
        if prop==5
          if (Twbset == 1)
            awprop = Twb;
          else 
            TwblowC = -99.9;
%           The highest wet bulb temp possible is at sat. temp for pressure
            TwbhighC = tpsteam(press)-tref;
            diffH = 1; 
            while (diffH > 1e-2)    
              TwbC = (TwblowC+TwbhighC)/2; 
              Twb = TwbC + Tref;   
%             Intermediate calc.
              if (TwbC <= 0) & (TwbC > -100) 
                pwswbint = C1/Twb+C2+C3*Twb+(C4*Twb^2)+(C5*Twb^3)+(C6*Twb^4)+C7*log(Twb);
              elseif (TwbC > 0) & (TwbC <= 200)    
                pwswbint = C8/Twb+C9+C10*Twb+(C11*Twb^2)+(C12*Twb^3)+C13*log(Twb);
              end  
%             Actual saturation pressure (Pa)
              if (TwbC >= -100) & (TwbC <= 200.0)
                pwswb = exp(pwswbint)/1000;
              else
                pwswb = ptsteam(Twb)*1000;
              end
       
%             Humidity ratio of saturated air (kg water/kg dry air)
              wswb = 0.62198*(pwswb/(p-pwswb));
%             Enthalpy (kJ/kg dry air) 
              Hwb = (1.006*TwbC)+wswb*(2500.45+1.805*TwbC);

              diffH = Hwb-H;
              if (diffH > 0)  
                TwbhighC = TwbC;  
              else 
                TwblowC = TwbC;  
              end
              diffH = abs(diffH);   
            end
            awprop = Twb;   
          end    
        end
      elseif (prop == 6) 
%       Dew point temperature (K)
        if (Tdpset == 1) 
          awprop = Tdp;  
        else   
%         Constants for determination  
          C14 = 6.54;  
          C15 = 14.526;  
          C16 = 0.7389;  
          C17 = 0.09486;  
          C18 = 0.4569;  
          alpha = log(pw); 
%         For dew point temperatures between 0 and 93 C
          TdpC = C14 + C15*alpha + C16*alpha^2 + C17*alpha^3 + C18*pw^0.1984; 
%         For dew point temperatures below 0 C
          if (TdpC < 0)
            TdpC = 6.09 + 12.608*alpha + 0.4959*alpha^2;  
          end  
          Tdp = TdpC + Tref;  
          awprop = Tdp; 
        end   
      end 
 
    properties=awprop; 
      
 


