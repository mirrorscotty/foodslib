      function [sys,x0] = hxs(t,x,u,flag,Uname,HIname,CIname,HOname,COname)

%  HXS       Heat exchanger model    
%
%            For use with FOODS-LIB models   
%
%            Synopsis: S-function for heat exchanger steady-state 
%            model. May determine mass and energy balances, perform design at 
%            basic transport level. 
%   
%            Call: unitops(t,x,u,flag,Uname,HIname,CIname,HOname,COname)
%   
%            Passed Arguments:  
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 and unit operation operating information  
%             Uname - unit operation name
%             HIname - hot stream inlet name 
%             CIname - cold stream inlet name 
%             HOname - hot stream outlet name 
%             COname - cold stream outlet name 
%  
%            Returned Arguments: 
%             sys - SIMULINK S-function variable: for flag = 0, 
%                   sys = sizes of parameters and initial conditions. 
%                   For flag = 3, sys = system outputs which are the 
%                   u-arrays for each stream exiting the unit operation 
%
%            References: 
%             Geankoplis, C.J. 1983.  "Transport Processes and Unit
%              Operations," 2nd ed.  Allyn & Bacon, Inc., Boston.
% 

%  ********************* SIMULINK INTERFACE INFORMATION ********************* 

SI = 2;                            % Number of system input streams
SO = 2;                            % Number of system output streams  
usize = usizemax;                  % Current u_array size for streams
UOinputs = 26;                     % Number of heat exchanger model inputs 

if flag == 0
   % sys(3) is number of outputs
   % sys(4) is number of inputs
   sys = [0, 0, (SO*usize), ((SI*usize)+UOinputs), 0, 1];
   x0 = ones(size(u));

elseif flag == 3 | flag == 9
   % if flag=3, return system outputs
  
%     Breakup of SIMULINK constructed u-array into stream arrays and 
%      operating information array  
 
      HIuarray = u(1:usize);          % full hot stream inlet array 
      CIuarray = u(usize+1:2*usize);  % full cold stream inlet array 
      UOuarray = u(SI*usize+1:SI*usize+UOinputs); 
                                     % full operating inputs array

%  *********************** USER-DEFINED DESIGN LEVEL ************************

      dlevel = UOuarray(1);           % design level
      mandq = UOuarray(2);            % micro and quality assessment
      econ = UOuarray(3);             % economic evaluation  

%  ******************* MASS BALANCE LEVEL VARIABLE ASSIGNMENTS **************

%     Hot stream inlet information  
      HI = HIuarray(ulocate('mass')); % mass flow rate (kg/h)
      HIinit = HI;
      HIcomp = HIuarray(ulocate('comp1st'):ulocate('complast'));
                                    % component mass fractions (w/w)
      HIctype = HIuarray(ulocate('ctype1st'):ulocate('ctypelast')); 
                                    % component types  
      HImass = [HI;HIcomp;HIctype];     % mass information array 

%     Cold stream inlet information 
      CI = CIuarray(ulocate('mass')); % mass flow rate (kg/h)
      CIinit = CI;
      CIcomp = CIuarray(ulocate('comp1st'):ulocate('complast'));
                                    % component mass fractions (w/w)
      CIctype = CIuarray(ulocate('ctype1st'):ulocate('ctypelast'));
                                    % component types  
      CImass = [CI;CIcomp;CIctype];     % mass information array


%  ***************************  MASS BALANCE LEVEL ************************** 

%     Overall Mass Balances (kg/h) 
      HO = HI;
      CO = CI; 

%     Component mass fractions (w/w) and types 
      HOcomp = HIcomp;
      COcomp = CIcomp;
      HOctype = HIctype;
      COctype = CIctype; 

%     Exit stream mass arrays; used for print results table  
      HOmass = [HO;HOcomp;HOctype];
      COmass = [CO;COcomp;COctype];

%     System output arrays 
      HOuarray = HOmass;  
      COuarray = COmass;   

      endtxt = [10 'END OF HEAT EXCHANGER (',Uname,') MASS BALANCE' 10];  
      disp(endtxt); 

%  **************** ENERGY BALANCE LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 2)  

%     Reference temperature (K)
      Tref = tref;

%     Hot stream inlet information  
      HItempC = HIuarray(ulocate('temp'));% temperature (C) 
      HItempK = HItempC + Tref;            % temperature (K) 
      HIpres = HIuarray(ulocate('pressure')); % pressure (MPa) 
      HIx = HIuarray(ulocate('quality')); % quality  
      HIme = [HImass;HItempK;HIpres;HIx]; 
                                         % mass and energy information array  

%     Cold stream inlet information
      CItempC = CIuarray(ulocate('temp')); % temperature (C) 
      CItempK = CItempC + Tref;            % temperature (K)  
      CIpres = CIuarray(ulocate('pressure'));    % pressure (MPa)   
      CIx = CIuarray(ulocate('quality'));   % quality
      CIme = [CImass;CItempK;CIpres;CIx];
                                         % mass and energy information array 

%     Unit operation information for energy balances  
      fconfig = UOuarray(4); % flow configuration: co-current (0); counter (1)
      delTmin = UOuarray(5); % minimum temperature differential (C)
      THOminC = UOuarray(6); % minimum hot stream outlet temperature (C) 
      Thoutmin = THOminC + Tref;   % (K) 
      HOxmin = UOuarray(7);  % minimum hot stream exit quality (if steam) 
      TCOmaxC = UOuarray(8); % maximum cold stream outlet temperature (C) 
      Tcoutmax = TCOmaxC + Tref;   % (K) 
      COxmax = UOuarray(9);  % maximum cold stream exit quality (if water)  
      Pcold = UOuarray(10);  % cold side maximum operating pressure (MPa) 

%  *********************** ENERGY BALANCE LEVEL *****************************

%     Feed stream checks for recycle problems
      if ((HI <= 0) & (CI <= 0)) | (HItempK <= Tref) | (CItempK <= Tref)
 
%  ((HI <= 0) & (HItempK <= Tref)) | ((CI <= 0) & (CItempK <= Tref))
        sys =[HIuarray;CIuarray];
        return;
      end
 
%     Check that minimum hot stream outlet temp is less than or equal to hot
%      stream inlet temp and greater than the cold stream inlet plus the
%      temperature differential
      if (Thoutmin > 0)
        if (Thoutmin > HItempK) | ((Thoutmin < (CItempK + delTmin)) & (Thoutmin > Tref))
          errortxt = [10 'HEAT EXCHANGER (',Uname,') ERROR: Minimum hot stream outlet temperature' 10 ...
          'is not less than or equal to hot stream inlet temperature and greater than' 10 ...
          'or equal to the cold stream inlet temperature plus the temperature differential.' 10];
          disp(errortxt);  
          HItempC
          TcplusdT = CItempC + delTmin
%         return; 
        end  

%     Check that maximum cold stream outlet temp is greater than or equal to
%      cold stream inlet temp and less than the hot stream inlet less
%      the temperature differential
      elseif (Tcoutmax > 0) & (Thoutmin <= 0)
        if (Tcoutmax < CItempK) | (Tcoutmax > (HItempK - delTmin))
          errortxt = [10 'HEAT EXCHANGER (',Uname,') ERROR: Maximum cold stream outlet temperature' 10 ... 
          'is not greater than or equal to cold stream inlet temperature and less than' 10 ... 
          'the hot stream inlet temperature less the temperature differential.' 10];
          disp(errortxt);
          CItempC 
          ThlessdT = HItempC - delTmin
%         return; 
        end  
      end 

%     Enthalpy of hot stream inlet (kJ/kg)  
      HIh = enthalpy(HIme);
%     Enthalpy of cold stream inlet (KJ/kg) 
      CIh = enthalpy(CIme);  
 
%     Hot stream type: type = 0 (food) or 1 (water/steam)
      hottype = 0;
      if (tsolids(HIcomp,HIctype) == 0 ) hottype = 1; end
%     Hot food stream water content in equals water content out (w/w)
      if (tsolids(HIcomp,HIctype) >= 0) HOx = HIx; end
%     Cold stream type: type = 0 (food) or 1 (water/steam)
      coldtype = 0;
      if (tsolids(CIcomp,CIctype) == 0) coldtype = 1; end
%     Cold food stream water content in equals water content out (w/w)
      if (tsolids(CIcomp,CIctype) >= 0) COx = CIx; end

%     Initializing iteration controllers  
      stop = 0;
      it = 1;

%     CASE 1 & 2: Hot and cold stream mass flowrates are both known
      if (HI > 0) & (CI > 0)

%       CASE 1: Cocurrent; hot and cold stream mass flows known
        if (fconfig == 0)
%         CASE 1A: Minimum hot stream outlet temperature is given
          if (Thoutmin > 0)
            while (stop == 0)  % Criteria for stopping search for HOtempK
%             Hot stream oulet enthalpy (kJ/kg)
              if (hottype == 1)
                if (it == 1)
                  if (Thoutmin == HItempK)
                    HOx = HOxmin;
                  else HOx = 0; end
                  it = 2;  
                end
                HOpres = ptsteam(Thoutmin);
                HOh = hsteam(Thoutmin,0,HOx);
              else
                HOx = HIx;
                HOcp = cp(HOcomp,HOctype,Thoutmin);
                HOh = HOcp*(Thoutmin-Tref);
                HOpres = HIpres;
              end

%             Cold stream outlet enthalpy (kJ/kg)
              COh = (HI*HIh+CI*CIh-HO*HOh)/CO;

%             Determination of cold stream outlet temperature (K)
%             Initial conditions for Bisection Method; want to find a cold
%             stream temperature that yields calculated cold stream enthalpy.
              hdiff = 1;
              Tlow = CItempK - 0.01;
              Thigh = Thoutmin - delTmin + 0.01;
              COtempK = CItempK;
              xlow = CIx;
              xhigh = COxmax;
              COx = CIx;
              if (coldtype == 1)      % For water/steam cold stream
                while (abs(hdiff) > 0.1e-8)
                  Pvap = ptsteam(COtempK);
                  if (COx == 0) & (Pvap < Pcold)
                    COtempK = (Tlow+Thigh)/2;
                    COpres = ptsteam(COtempK);
%                   Enthalpy based on last iteration temperature (K)
                    hcldout2 = hsteam(0,COpres,COx);
%                   Test condition assignment
                    hdiff = (hcldout2-COh)/COh;
%                   Reassignment of interval
                    if (hdiff > 0)
                      Thigh = COtempK;
                    elseif (hdiff < 0)
                      Tlow = COtempK;
                    end
                    hdiff = abs(hdiff);
                    if (COtempK > Thoutmin - delTmin)
                      if (Thoutmin < HItempK)
                        Thoutmin = Thoutmin + 1;
                      elseif (Thoutmin >= HItempK) & (hottype == 1)
                        Thoutmin = HItempK;
                        HOx = HOx + 0.01;
                      end
                      break;
                    end
                  elseif ((COx >= 0) & (Pvap >= Pcold)) | (CIx > 0)
                    COx = (xlow+xhigh)/2;
                    if (Pvap >= Pcold)
                      COpres = Pcold;
                      COtempK = tpsteam(COpres);
                    else
                      COpres = Pvap;
                      COtempK = CItempK;
                    end
%                   Enthalpy based on last iteration quality
                    hcldout2 = hsteam(0,COpres,COx);
%                   Test condition assignment
                    hdiff = (hcldout2-COh)/COh;
%                   Reassignment of interval
                    if (hdiff > 0)
                      xhigh = COx;
                    elseif (hdiff < 0)
                      xlow = COx;
                    end
                    hdiff = abs(hdiff);
                    if ((1-COx) < 1e-4) & (hdiff > 0.1e-8)
                      if (Thoutmin < HItempK)
                        Thoutmin = Thoutmin + 1;
                      elseif (Thoutmin >= HItempK) & (hottype == 1)
                        Thoutmin = HItempK;
                        HOx = HOx + 0.01;
                      end
                      break;
                    end
                  end
                end
 
              else                    % For cold food stream
                while (abs(hdiff) > 0.1e-8)  
                  COtempK = (Tlow+Thigh)/2;
%                 Specific heat of exit heat transfer stream (kJ/kg(stream)-K)
                  COcp = cp(COcomp,COctype,COtempK);
%                 Enthalpy based on last iteration temperature (K)
                  hcldout2= COcp*(COtempK-Tref);
%                 Test condition assignment
                  hdiff = (hcldout2-COh)/COh;
%                 Reassignment of interval
                  if (hdiff > 0)
                    Thigh = COtempK;
                  elseif (hdiff < 0)
                    Tlow = COtempK;
                  end
                  hdiff = abs(hdiff);  
                  if (COtempK > Thoutmin - delTmin)
                    if (Thoutmin < HItempK)
                      Thoutmin = Thoutmin + 1;
                    elseif (Thoutmin >= HItempK) & (hottype == 1)
                      Thoutmin = HItempK;
                      HOx = HOx + 0.01;
                    end
                    break;
                  end
                end % while hdiff
                COpres = CIpres;
              end % if coldtype
              if (Thoutmin > HItempK) | (HOx > HIx)
                disp('HEAT EXCHANGER ERROR: Hot stream flow rate is too high');           
                return;
              end
              if (hdiff <= 0.1e-8)
                HOtempK = Thoutmin;
                stop = 1;
              end
            end
 
%         CASE 1B: Maximum cold stream outlet temperature given
          elseif (Tcoutmax > 0) & (Thoutmin <= 0)
            while (stop == 0)  % Criteria to stop searching for COtempK
 
%             Cold stream oulet enthalpy (kJ/kg)
              if (coldtype == 1)
                Pvap = ptsteam(Tcoutmax);  
                if (Pvap > Pcold)
                  Tcoutmax = tpsteam(Pcold);
                  COx = COxmax;
                end
                if (it == 1)
                 if (Tcoutmax == CItempK)
                    COx = CIx;
                  else COx = CIx; end
                  it = 2;  
                end
                COpres = ptsteam(Tcoutmax);
                COh = hsteam(0,COpres,COx);
              else
                COx = CIx;
                COcp = cp(COcomp,COctype,Tcoutmax);
                COh = COcp*(Tcoutmax-Tref);
                COpres = CIpres;
              end
 
%             Cold stream outlet enthalpy (kJ/kg)
              HOh = (HI*HIh+CI*CIh-CO*COh)/H0;
 
%             Determination of hot stream outlet temperature (K)
%             Initial conditions for Bisection Method; want to find hot
%             stream temperature that yields calculated hot stream enthalpy.
              hdiff = 1;
              Tlow = Tcoutmax + delTmin - 0.01;
              Thigh = HItempK;
              HOtempK = HItempK;
              xlow = 0;
              xhigh = HIx;
              HOx = HIx;
              if (hottype == 1)      % For water/steam cold stream
                while (hdiff > 0.1e-8)
                  if (HOx == 0) & (HOtempK <= HItempK)
                    HOtempK = (Tlow+Thigh)/2;
                    HOpres = ptsteam(HOtempK);
%                   Enthalpy based on last iteration temperature (K)
                    HOh2 = hsteam(HOtempK,0,HOx);
%                   Test condition assignment
                    hdiff = (HOh2-HOh)/HOh;
%                   Reassignment of interval
                    if (hdiff > 0)
                      Thigh = HOtempK;
                    elseif (hdiff < 0)
                      Tlow = HOtempK;
                    end
                    hdiff = abs(hdiff);
                    if (HOtempK < Tcoutmax + delTmin)
                      if (Tcoutmax > CItempK)
                        Tcoutmax = Tcoutmax - 1;
                      elseif (Tcoutmax <= CItempK) & (coldtype == 1)
                        Pvap = ptsteam(Tcoutmax);  
                        if (Pvap > Pcold)
                          Tcoutmax = tpsteam(Pcold);
                          COx = COx + 0.01;
                        end
                        if (Tcoutmax <= CItempK)
                          Tcoutmax = CItempK;
                          COx = COx + 0.01;
                        end
                      end
                      break;
                    end
                  elseif (HOx >= 0) & (HOtempK == HItempK)
                    HOx = (xlow+xhigh)/2;
                    if (HOx < 1e-4) HOx = 0; end
                    HOpres = ptsteam(HOtempK);
%                   Enthalpy based on last iteration quality
                    HOh2 = hsteam(HOtempK,0,HOx);
%                   Test condition assignment
                    hdiff = (HOh2-HOh)/HOh;
 
%                   Reassignment of interval
                    if (hdiff > 0)
                      xhigh = HOx;
                    elseif (hdiff < 0)
                      xlow = HOx;
                    end
                    hdiff = abs(hdiff);
                    if ((1-HOx) < 1e-4) & (hdiff > 0.1e-8)
                      if (Tcoutmax > CItempK)
                        Tcoutmax = Tcoutmax - 1;
                      elseif (Tcoutmax <= CItempK) & (coldtype == 1)
                        Pvap = ptsteam(Tcoutmax);
                        if (Pvap > Pcold)
                          Tcoutmax = tpsteam(Pcold)
                          COx = COx + 0.01;
                        end
                        if (Tcoutmax <= CItempK)
                          Tcoutmax = CItempK;
                          COx = COx + 0.01;
                        end
                      end
                      break;
                    end
                  end
                end % while hdiff
                COpres = ptsteam(Tcoutmax);
              else                    % For hot food stream
                while (hdiff > 0.1e-8)
                  HOtempK = (Tlow+Thigh)/2;
%                 Specific heat of hot outlet stream (kJ/kg(stream)-K)
                  HOcp = cp(HOcomp,HOctype,HOtempK);
%                 Enthalpy based on last iteration temperature (K)
                  HOh2= HOcp*(HOtempK-Tref);
%                 Test condition assignment
                  hdiff = (HOh2-HOh)/HOh;
%                 Reassignment of interval
                  if (hdiff > 0)
                    Thigh = HOtempK;
                  elseif (hdiff < 0)
                    Tlow = HOtempK;
                  end
                  hdiff = abs(hdiff);
                  if (HOtempK < Tcoutmax + delTmin)

                    if (Tcoutmax > CItempK)
                      Tcoutmax = Tcoutmax - 1;
                    end
                    if (coldtype == 1)
                      Pvap = ptsteam(Tcoutmax);
                      if (Pvap > Pcold)
                        Tcoutmax = tpsteam(Pcold)
                        COx = COx + 0.01;
                      end
                      if (Tcoutmax <= CItempK)
                        Tcoutmax = CItempK;
                        COx = COx + 0.01;
                      end
                    end
                    COpres = ptsteam(Tcoutmax);
                    break;
                  end
                end % while hdiff
                HOpres = HIpres;
              end
              if (Tcoutmax < CItempK) | (COx > 1) stop = 1; end
              if (hdiff <= 0.1e-8)
                COtempK = Tcoutmax;
                stop = 1;
              end
            end % while stop == 0
          end % CASE 1B

%       CASE 2: Counter-current; hot and cold stream mass flows known
        elseif (fconfig == 1)
%         CASE 2A: Minimum hot stream outlet temperature given
          if (Thoutmin > 0)
            while (stop == 0)  % Criteria to stop search for HOtempK
%             Hot stream oulet enthalpy (kJ/kg)
              if (hottype == 1)
                if (it == 1)
                  if (Thoutmin == HItempK)
                    HOx = HOxmin;
                  else HOx = 0; end
                  it = 2;  
                end
                HOpres = ptsteam(Thoutmin);
                HOh = hsteam(Thoutmin,0,HOx);
              else
                HOx = HIx;
                HOcp = cp(HOcomp,HOctype,Thoutmin);
                HOh = HOcp*(Thoutmin-Tref);  
                HOpres = HIpres;
              end

%             Cold stream outlet enthalpy (kJ/kg)
              COh = (HI*HIh+CI*CIh-HO*HOh)/CO

%             Determination of cold stream outlet temperature (K)
%             Initial conditions for Bisection Method; want to find a cold
%             stream temperature that yields calculated cold stream enthalpy.
              hdiff = 1;
              Tlow = CItempK - 0.01;
              Thigh = HItempK - delTmin + 0.01;
              COtempK = CItempK;
              xlow = CIx;
              xhigh = COxmax;
              COx = CIx;
              if (coldtype == 1)      % For water/steam cold stream
                while (hdiff > 0.1e-8)
                  Pvap = ptsteam(COtempK);  
                  if (COx == 0) & (Pvap < Pcold)
                    COtempK = (Tlow+Thigh)/2
                    COpres = ptsteam(COtempK);
%                   Enthalpy based on last iteration temperature (K)
                    hcldout2 = hsteam(0,COpres,COx)
%                   Test condition assignment
                    hdiff = (hcldout2-COh)/COh;
%                   Reassignment of interval
                    if (hdiff > 0)
                      Thigh = COtempK;
                    elseif (hdiff < 0)
                      Tlow = COtempK;
                    end
                    hdiff = abs(hdiff);
                    if (COtempK > HItempK - delTmin)
                      if (Thoutmin < HItempK)
                        Thoutmin = Thoutmin + 1;
                      elseif (Thoutmin >= HItempK) & (hottype == 1)
                        Thoutmin = HItempK;
                        HOx = HOx + 0.01;
                      end
                      break;
                    end
                  elseif ((COx >= 0) & (Pvap >= Pcold)) | (CIx > 0)
                    COx = (xlow+xhigh)/2;
                    if (Pvap >= Pcold)
                      COpres = Pcold;
                      COtempK = tpsteam(COpres);
                    else
                      COpres = Pvap;
                      COtempK = CItempK;
                    end
%                   Enthalpy based on last iteration quality
                    hcldout2 = hsteam(COtempK,0,COx);
%                   Test condition assignment
                    hdiff = (hcldout2-COh)/COh;
%                   Reassignment of interval
                    if (hdiff > 0)
                      xhigh = COx;
                    elseif (hdiff < 0)
                      xlow = COx;
                    end
                    hdiff = abs(hdiff);
                    if ((1-COx) < 1e-4) & (hdiff > 0.1e-8)
                      if (Thoutmin < HItempK)
                        Thoutmin = Thoutmin + 1;
                      elseif (Thoutmin >= HItempK) & (hottype == 1)
                        Thoutmin = HItempK;
                        HOx = HOx + 0.01;
                      end
                      break;
                    end
                  end
                end
 
              else                    % For food stream
                while (hdiff > 0.1e-8)
                  COtempK = (Tlow+Thigh)/2;
%                 Specific heat of exit heat transfer stream (kJ/kg(stream)-K)
                  COcp = cp(COcomp,COctype,COtempK);
%                 Enthalpy based on last iteration temperature (K)
                  hcldout2= COcp*(COtempK-Tref);
%                 Test condition assignment
                  hdiff = (hcldout2-COh)/COh;
%                 Reassignment of interval
                  if (hdiff > 0)
                    Thigh = COtempK;
                  elseif (hdiff < 0)
                    Tlow = COtempK;
                  end
                  hdiff = abs(hdiff);
                  if (COtempK > HItempK - delTmin)
                    if (Thoutmin < HItempK)
                      Thoutmin = Thoutmin + 1;
                    elseif (Thoutmin >= HItempK) & (hottype == 1)
                      Thoutmin = HItempK;
                      HOx = HOx + 0.01;
                    end
                    break;
                  end
                end % while hdiff
                COpres = CIpres;
              end % if coldtype
              if (Thoutmin > HItempK) | (HOx > HIx)
                disp('HEAT EXCHANGER ERROR: Hot stream flow rate is too high');    
                return;
              end
              if (hdiff <= 0.1e-8)
                HOtempK = Thoutmin;
                stop = 1;
              end
            end
 
%         CASE 2B: Maximum cold stream outlet temperature given
          elseif (Tcoutmax > 0) & (Thoutmin <= 0)
 
            while (stop == 0)  % Criteria to stop searching for COtempK
 
%             Cold stream oulet enthalpy (kJ/kg)
              if (coldtype == 1)
                Pvap = ptsteam(Tcoutmax);  
                if (Pvap > Pcold)
                  Tcoutmax = tpsteam(Pcold);
                  COx = COxmax;
                end
                if (it == 1)
                 if (Tcoutmax == CItempK)
                    COx = CIx;
                  else COx = CIx; end
                  it = 2;  
                end
                COpres = ptsteam(Tcoutmax);
                COh = hsteam(0,COpres,COx);
              else
                COx = CIx;
                COcp = cp(COcomp,COctype,Tcoutmax);
                COh = COcp*(Tcoutmax-Tref);
                COpres = CIpres;
              end
 
%             Cold stream outlet enthalpy (kJ/kg)
              HOh = (HI*HIh+CI*CIh-CO*COh)/HO;  
 
 
              if (HOh < 0)
                disp('HEAT EXCHANGER ERROR: Hot stream flow rate is too low.');
                return;
              end
 
%             Determination of hot stream outlet temperature (K)
%             Initial conditions for Bisection Method; want to find hot
%             stream temperature that yields calculated hot stream enthalpy.
              hdiff = 1;
              Tlow = CItempK + delTmin - 0.01;
              Thigh = HItempK;
              HOtempK = HItempK;
              xlow = 0;
              xhigh = HIx;
              HOx = HIx;
              if (hottype == 1)      % For water/steam cold stream
                while (hdiff > 0.1e-8)
                  if (HOx == 0) & (HOtempK <= HItempK)
                    HOtempK = (Tlow+Thigh)/2;  
                    HOpres = ptsteam(HOtempK);
%                   Enthalpy based on last iteration temperature (K)
                    HOh2 = hsteam(HOtempK,0,HOx)
%                   Test condition assignment
                    hdiff = (HOh2-HOh)/HOh;
%                   Reassignment of interval
                    if (hdiff > 0)
                      Thigh = HOtempK;
                    elseif (hdiff < 0)
                      Tlow = HOtempK;
                    end
                    hdiff = abs(hdiff);
                    if (HOtempK < CItempK + delTmin)
                      if (Tcoutmax > CItempK)
                        Tcoutmax = Tcoutmax - 1;
                      elseif (Tcoutmax <= CItempK) & (coldtype == 1)
                        Pvap = ptsteam(Tcoutmax);  
                        if (Pvap > Pcold)
                          Tcoutmax = tpsteam(Pcold);
                          COx = COx + 0.01;
                        end
                        if (Tcoutmax <= CItempK)
                          Tcoutmax = CItempK;
                          COx = COx + 0.01;
                        end
                      end
                      break;
                    end
                  elseif (HOx >= 0) & (HOtempK == HItempK)
                    HOx = (xlow+xhigh)/2;  
                    if (HOx < 1e-4) HOx = 0; end
                    HOpres = ptsteam(HOtempK);
%                   Enthalpy based on last iteration quality
                    HOh2 = hsteam(HOtempK,0,HOx)
%                   Test condition assignment
                    hdiff = (HOh2-HOh)/HOh;
%                   Reassignment of interval
                    if (hdiff > 0)
                      xhigh = HOx;
                    elseif (hdiff < 0)
                      xlow = HOx;
                    end
                    hdiff = abs(hdiff);
 
                    if ((1-HOx) < 1e-4) & (hdiff > 0.1e-8)
                      if (Tcoutmax > CItempK)
                        Tcoutmax = Tcoutmax - 1;
                      elseif (Tcoutmax <= CItempK) & (coldtype == 1)
                        Pvap = ptsteam(Tcoutmax);
                        if (Pvap > Pcold)
                          Tcoutmax = tpsteam(Pcold)
                          COx = COx + 0.01;
                        end
                        if (Tcoutmax <= CItempK)
                          Tcoutmax = CItempK;
                          COx = COx + 0.01;
                        end
                      end
                      break;
                    end
                  end
                end % while hdiff
                COpres = ptsteam(Tcoutmax);
              else                    % For food stream
                while (hdiff > 0.1e-8)
                  HOtempK = (Tlow+Thigh)/2;  
%                 Specific heat of hot outlet stream (kJ/kg(stream)-K)
                  HOcp = cp(HOcomp,HOctype,HOtempK);
%                 Enthalpy based on last iteration temperature (K)
                  HOh2= HOcp*(HOtempK-Tref);  
%                 Test condition assignment
                  hdiff = (HOh2-HOh)/HOh;  
%                 Reassignment of interval
                  if (hdiff > 0)
                    Thigh = HOtempK;
                  elseif (hdiff < 0)
                    Tlow = HOtempK;
                  end
                  hdiff = abs(hdiff);
                  if (HOtempK < CItempK + delTmin)
                    if (Tcoutmax > CItempK)
                      Tcoutmax = Tcoutmax - 1;
                    elseif (Tcoutmax <= CItempK) & (coldtype == 1)
                      Pvap = ptsteam(Tcoutmax);
                      if (Pvap > Pcold)
                        Tcoutmax = tpsteam(Pcold)
                        COx = COx + 0.01;
                      end
                      if (Tcoutmax <= CItempK)
                        Tcoutmax = CItempK;
                        COx = COx + 0.01;
                      end
                    end
                    COpres = ptsteam(Tcoutmax);
                    break;
                  end
                end
                HOpres = HIpres;
              end % hottype
              if (Tcoutmax < CItempK) | (COx > 1) stop = 1; end
              if (hdiff <= 0.1e-8)
                COtempK = Tcoutmax;  
                stop = 1;
              end
            end % while stop = 0
          end % CASE 2B
        end % CASE 2

      end % CASE 1 & 2
 
%     CASE 3 & 4: Cold stream mass flowrate unknown
      if (CI == 0) & (HI > 0)
 
%       CASE 3A and 4A: Hot stream outlet temperature is given
        if (Thoutmin > 0)
          
%         Hot stream oulet enthalpy (kJ/kg)
          if (hottype == 1)
            if (Thoutmin == HItempK);
%             Minimum hot stream outlet temperature equals the hot stream 
%             inlet temperature. So use min outlet quality. 
              HOx = HOxmin;  
            else 
              HOx = 0; 
            end
            HOpres = ptsteam(Thoutmin);
            HOh = hsteam(Thoutmin,0,HOx);
          else
            HOx = HIx;
            HOcp = cp(HOcomp,HOctype,Thoutmin);
            HOh = HOcp*(Thoutmin-Tref);  
            HOpres = HIpres;
          end
 
          HOtempK = Thoutmin;
 
          if (coldtype == 0) | ((coldtype == 1) & (CIx == 0))
%           CASE 3A: Cocurrent
            if (fconfig == 0)
%             Cold stream outlet temperature (K)
              if (Tcoutmax > Tref) & (Tcoutmax < Thoutmin - delTmin)
                COtempK = Tcoutmax;
              else
                COtempK = Thoutmin - delTmin;
              end
%           CASE 4A: Counter-current
            elseif (fconfig == 1)
%             Cold stream outlet temperature (K)
              if (Tcoutmax > Tref) & (Tcoutmax < HItempK - delTmin)
                COtempK = Tcoutmax;
              else
                COtempK = HItempK - delTmin;
              end
            end
            if (coldtype == 1)
              Pvap = ptsteam(COtempK);  
              if (Pvap > Pcold)
                COtempK = tpsteam(Pcold);
%               Cold stream outlet temperature exceeds the boiling point. 
%               So use cold stream max outlet quality 
                COx = COxmax;
              else
                COx = CIx;  
              end
            end
          elseif (coldtype == 1) & (CIx > 0)
            COtempK = CItempK; 
%           Cold stream outlet quality exceeds zero. So use cold stream
%           max outlet quality  
            COx = COxmax;
          end
 
%         Cold stream outlet enthalpy (kJ/kg)
          if (coldtype == 1)
            COpres = ptsteam(COtempK);
            COh = hsteam(0,COpres,COx);
          else
            COx = CIx;
            COcp = cp(COcomp,COctype,COtempK);
            COh = COcp*(COtempK-Tref);  
            COpres = CIpres;
          end
 
%         Cold stream mass flow rate (kg/h)
          CI = HO*(HOh-HIh)/(CIh-COh);
          CO = CI;  

        end % CASE 3A & 4A  
          
%       CASE 3B and 4B: Cold stream outlet temperature is given
        if (Tcoutmax > 0) & (Thoutmin <= 0)
          
%         Cold stream oulet temperature (K)
          COtempK = Tcoutmax;

%         Cold stream oulet enthalpy (kJ/kg)
          if (coldtype == 1)
            Pvap = ptsteam(Tcoutmax);  
            if (Tcoutmax == CItempK) | (Pvap > Pcold);
%             Maximum cold stream outlet temperature equals the cold stream 
%             inlet temperature. So use cold stream max outlet quality
              COx = COxmax 
              if (Pvap > Pcold) COtempK = tpsteam(Pcold); end 
            else 
              if (CIx > 0) 
                COtempK = CItempK;
                COx = COxmax;
              else 
                COx = CIx;  
              end
            end
            COpres = ptsteam(Tcoutmax);
            COh = hsteam(0,COpres,COx);
          else
            COx = CIx;
            COcp = cp(COcomp,COctype,Tcoutmax);
            COh = COcp*(Tcoutmax-Tref)
            COpres = CIpres;
          end
 
          if (hottype == 0) | ((hottype == 1) & (HIx == 0))
            HOx = HIx
%           CASE 3B: Cocurrent
            if (fconfig == 0)
              HOtempK = COtempK + delTmin;
%           CASE 4B: Counter-current
            elseif (fconfig == 1)
%             Hot stream outlet temperature (K)
              HOtempK = CItempK + delTmin;
            end
          end
          if (hottype == 1) & (HIx > 0)
            HOtempK = HItempK; 
%           Hot stream exit temperature equals hot stream inlet temperature. 
%           So use hot stream mn outlet quality
            HOx = HOxmin;
          end
 
%         Hot stream outlet enthalpy (kJ/kg)
          if (hottype == 1)
            HOpres = ptsteam(HOtempK);
            HOh = hsteam(HOtempK,0,HOx);
          else
            HOcp = cp(HOcomp,HOctype,HOtempK);
            HOh = HOcp*(HOtempK-Tref);  
            HOpres = HIpres;
          end
 
%         Cold stream mass flow rate (kg/h)
          CI = HO*(HOh-HIh)/(CIh-COh);
          CO = CI;  
        end % CASE 3B & 4B  
          
      end % CASE 3 & 4  

%     CASE 5 & 6: Hot stream mass flowrate unknown
      if (HI == 0) & (CI > 0)

%       CASE 5A and 6A: Hot stream outlet temperature is given
        if (Thoutmin > 0)

%         Hot stream oulet enthalpy (kJ/kg)
          if (hottype == 1)
            if (Thoutmin == HItempK);
%             Minimum hot stream outlet temperature equals the hot stream
%             inlet temperature. So, use min hot stream outlet quality
              HOx = HOxmin;   
            else
              if (HIx > 0)
                HOtempK = HItempK;  
                HOx = HOxmin;  
              else
                HOx = HIx;  
              end
            end  
            HOpres = ptsteam(Thoutmin);
            HOh = hsteam(Thoutmin,0,HOx);  
          else
            HOx = HIx;
            HOcp = cp(HOcomp,HOctype,Thoutmin);
            HOh = HOcp*(Thoutmin-Tref);  
            HOpres = HIpres;
          end

          HOtempK = Thoutmin;
 
          if (coldtype == 0) | ((coldtype == 1) & (CIx == 0))
            COx = CIx;
%           CASE 5A: Cocurrent
            if (fconfig == 0)
%             Cold stream outlet temperature (K)
              if (Tcoutmax > Tref) & (Tcoutmax < Thoutmin - delTmin)
                COtempK = Tcoutmax;
              else
                COtempK = Thoutmin - delTmin;
              end
%           CASE 6A: Counter-current
            elseif (fconfig == 1)
%             Cold stream outlet temperature (K)
              if (Tcoutmax > Tref) & (Tcoutmax < HItempK - delTmin)
                COtempK = Tcoutmax;
              else
                COtempK = HItempK - delTmin;
              end
            end
            if (coldtype == 1)
              Pvap = ptsteam(COtempK);  
              if (Pvap > Pcold)
                COtempK = tpsteam(Pvap);
%               Cold stream outlet temperature exceeds the boiling point.
%               So use max cold stream outlet quality
                COx = COxmax;   
              else
                COx - CIx;
              end
            end
          elseif (coldtype == 1) & (CIx > 0)
            COtempK = CItempK;  
%           Cold stream outlet quality exceeds zero. So use max cold stream
%           outlet quality
            COx = COxmax; 
          end
 
%         Cold stream outlet enthalpy (kJ/kg)
          if (coldtype == 1)
            COpres = ptsteam(COtempK);
            COh = hsteam(0,COpres,COx);
          else
            COx = CIx;
            COcp = cp(COcomp,COctype,COtempK);
            COh = COcp*(COtempK-Tref);  
            COpres = CIpres;
          end
 
%         Hot stream mass flow rate (kg/h)
          HI = CI*(CIh-COh)/(HOh-HIh);
          HO = HI;    
 
        end % CASE 5A & 6A
 
%       CASE 5B and 6B: Cold stream outlet temperature is given
        if (Tcoutmax > 0) & (Thoutmin <= 0)
 
%         Cold stream oulet enthalpy (kJ/kg)
          if (coldtype == 1)
            Pvap = ptsteam(Tcoutmax);  
            if (Tcoutmax == CItempK) | (Pvap > Pcold);
%             Maximum cold stream outlet temperature equals the cold stream
%             inlet temperature. So use max cold outlet quality
              COx = COxmax
              if (Pvap > Pcold) Tcoutmax = tpsteam(Pcold); end
            else
              if (CIx > 0)
                COtempK = CItempK;  
                COx = COxmax;  
              else
                COx = CIx;  
              end
            end
            COpres = ptsteam(Tcoutmax);
            COh = hsteam(0,COpres,COx);
          else
            COx = CIx;
            COcp = cp(COcomp,COctype,Tcoutmax);
            COh = COcp*(Tcoutmax-Tref);  
            COpres = CIpres;
          end
 
          COtempK = Tcoutmax;
 
          if (hottype == 1) & (HIx > 0)
%           Hot stream exit temperature equals hot stream inlet temperature.
%           So use min hot stream outlet quality
            HOx = HOxmin;
            HOtempK = HItempK;
          else
%           CASE 5B: Cocurrent
            if (fconfig == 0)
%             Hot stream outlet temperature (K)
              HOtempK = COtempK + delTmin;
%           CASE 6B: Counter-current
            elseif (fconfig == 1)
%             Hot stream outlet temperature (K)
              HOtempK = CItempK + delTmin;
            end
          end
 
%         Hot stream outlet enthalpy (kJ/kg)
          if (hottype == 1)
            if (HOtempK ~= HItempK) 
              HOpres = ptsteam(HOtempK);
              HOh = hsteam(HOtempK,0,HOx);
            else 
              HOpres = HIpres;  
              HOh = hsteam(HOtempK,0,HOx);
            end
          else
            HOx = HIx;
            HOcp = cp(HOcomp,HOctype,HOtempK);
            HOh = HOcp*(HOtempK-Tref);  
            HOpres = HIpres;
          end
 
%         Hot stream mass flow rate (kg/h)
          HI = CI*(CIh-COh)/(HOh-HIh);
          HO = HI;  
 
        end % CASE 5B & 6B
 
      end % CASE 5 & 6
 
%     Update mass information arrays
      HImass = [HI;HIcomp;HIctype];   
      HOmass = [HO;HOcomp;HOctype];   
      CImass = [CI;CIcomp;CIctype];
      COmass = [CO;COcomp;COctype];

%     Update mass and energy info arrays 
      HIme = [HImass;HItempK;HIpres;HIx]; 
      HOme = [HOmass;HOtempK;HOpres;HOx]; 
      CIme = [CImass;CItempK;CIpres;CIx]; 
      COme = [COmass;COtempK;COpres;COx]; 

%     System output arrays
      HOuarray = HOme;
      COuarray = COme;

      endtxt = [10 'END OF HEAT EXCHANGER (',Uname,') ENERGY BALANCE' 10];
      disp(endtxt);

      end % Energy Balance Level  

%  *************** BASIC TRANSPORT LEVEL VARIABLE ASSIGNMENTS ***************

      if (dlevel >= 3)

%     Unit operation information for basic transport (sizing) 
      U = UOuarray(11);        % overall heat transfer coefficient  
      plwidth = UOuarray(12);  % plate width (m)     
      pllength = UOuarray(13);   % plate length (m)    
      plgap = UOuarray(14);    % plate gap (m)    

      end % Basic Transport Level Assignments  

%  ************************** BASIC TRANSPORT LEVEL *************************

      if (dlevel >= 3) 

%     Heat transfer on hot side (kJ/h)
      Qhot = HI*(HIh - HOh);  

%     Log mean temperature differential for Cocurrent
      if (fconfig == 0)
        Tlogmean = ((HItempK-CItempK)-(HOtempK-COtempK))/ ... 
                   log((HItempK-CItempK)/(HOtempK-COtempK));  
%     Log mean temperature differential for Counter-current
      else
        Tlogmean = ((HItempK-COtempK)-(HOtempK-CItempK))/ ... 
                   log((HItempK-COtempK)/(HOtempK-CItempK));  
      end
      Tdiff = Tlogmean;

%     In case of Tlogmean is NaN; use algebraic mean difference  
      if (isnan(Tlogmean)) 
        if (fconfig == 0)
          Tdiff = ((HItempK - CItempK) + (HOtempK - COtempK))/2; 
        elseif (fconfig == 1) i
          Tdiff = ((HItempK - COtempK) + (HOtempK - CItempK))/2;  
        end 
      end

%     Heat Transfer Area (m^2)
      areahx = Qhot/(U*Tdiff);  

%     Average density (kg/m^3)
      Crho = (density(CIme)+density(COme))/2; 
      Hrho = (density(HIme)+density(HOme))/2; 
 
%     Volumetric flow rates (m^3/h)
      CVq = CI/Crho;   
      HVq = HI/Hrho;  

%     Path length through heat exchanger (m)
      length = areahx/(2*plwidth);  

%     Number of plates 
      Nplates = length/pllength; 

      endtxt = [10 'END OF UNIT (',Uname,') BASIC TRANSPORT' 10];
      disp(endtxt);

      end  % Basic Transport Level 

%  ************* ADVANCED TRANSPORT LEVEL VARIABLE ASSIGNMENTS **************

      if (dlevel >= 4) 

%     Unit operation inforamtion for advanced transport  
      Hmumodel = UOuarray(15);  % Hot stream viscosity model
      Cmumodel = UOuarray(16);  % Cold stream viscosity model

%  ************************ ADVANCED TRANSPORT LEVEL ************************

%     Hydraulic diameter (m)  
      Dh = 4*(plwidth*plgap)/(2*plwidth+2*plgap);

%     Volume (m^3) 
      Vol = plwidth*plgap*length;

%     Hot stream inlet & outlet density (kg/m^3) 
      HIrho = density(HIme);
      HOrho = density(HOme);

%     Hot stream volumetric flow rate (m^3/h)
      HIVq = HI/HIrho;  
      HOVq = HO/HOrho;  

%     Hot stream average velocity (m/s) 
      HIvel = HIVq/(plwidth*plgap*3600); 
      HOvel = HOVq/(plwidth*plgap*3600); 

%     Hot stream inlet & outlet viscosity (Pa-s) 
      if (hottype == 1) 
        HIK = vissteam(HItempK,0,HIx);
        HIn = 1;
        HOK = vissteam(HOtempK,0,HOx); 
        HOn = 1; 
      else
        HIplaw = viscos(HIcomp,HIctype,HItempK,Hmumodel);  
        HIK = HIplaw(1);
        HIn = HIplaw(2); 
        HOplaw = viscos(HOcomp,HOctype,HOtempK,Hmumodel);
        HOK = HOplaw(1); 
        HOn = HOplaw(2);    
      end 

%     Hot stream inlet & outlet Reynolds Number
      HINRe = reynolds(Dh,Dh,HIvel,HIrho,HIK,HIn,0);
      HONRe = reynolds(Dh,Dh,HOvel,HOrho,HOK,HOn,0);

%     Cold stream inlet & outlet density (kg/m^3)
      CIrho = density(CIme);
      COrho = density(COme);
 
%     Cold stream volumetric flow rate (m^3/h)
      CIVq = CI/CIrho;
      COVq = CO/COrho;
 
%     Cold stream average velocity (m/s)  
      CIvel = CIVq/(plwidth*plgap*3600);
      COvel = COVq/(plwidth*plgap*3600);  
 
%     Cold stream inlet & outlet viscosity (Pa-s)
      if (hottype == 1) 
        CIK = vissteam(CItempK,0,CIx);
        CIn = 1;   
        COK = vissteam(COtempK,0,COx); 
        COn = 1;   
      else
        CIplaw = viscos(CIcomp,CIctype,CItempK,Cmumodel);  
        CIK = CIplaw(1);
        CIn = CIplaw(2); 
        COplaw = viscos(COcomp,COctype,COtempK,Cmumodel);  
        COK = COplaw(1); 
        COn = COplaw(2);    
      end
 
%     Cold stream inlet & outlet Reynolds Number
      CINRe = reynolds(Dh,Dh,CIvel,CIrho,CIK,CIn,0);  
      CONRe = reynolds(Dh,Dh,COvel,COrho,COK,COn,0);  

      endtxt = [10 'END OF UNIT (',Uname,') ADVANCED TRANSPORT' 10]; 
      disp(endtxt); 

      end  % Advanced Transport Level 

%  *************** MICROBIAL & QUALITY VARIABLE ASSIGNMENTS *****************

%     Hot Feed Microbial Death Kinetics
      HDvalm = HIuarray(ulocate('microDval'));   % D value (min)
      HDtempCm = HIuarray(ulocate('microDtemp'));% D value temperature basis (C)
      HDtempKm = HDtempCm+Tref; % D value temperature basis (K)
      HZvalm =  HIuarray(ulocate('microZval'));   % Z value (C)
      HEam =  HIuarray(ulocate('microEa'));       % activation energy (kcal/mol)
      HNoin =  HIuarray(ulocate('microload'));     % microbial load
      HImdk = [HDvalm;HDtempCm;HZvalm;HEam;HNoin];  % microbial array

%     Hot Feed Quality Kinetics
      HDvalq = HIuarray(ulocate('qualDval'));    % D value (min)
      HDtempCq = HIuarray(ulocate('qualDtemp')); % D value temperature basis (C)
      HDtempKq = HDtempCq+Tref; % D value temperature basis (K)
      HZvalq = HIuarray(ulocate('qualZval'));     % Z value (C)
      HEaq = HIuarray(ulocate('qualEa'));       % activation energy (kcal/mole)  
      HQoin = HIuarray(ulocate('qualload'));     % quality load
      HIqual = [HDvalq;HDtempCq;HZvalq;HEaq;HQoin];  % quality array

%     Cold Feed Microbial Death Kinetics
      CDvalm = CIuarray(ulocate('microDval'));   % D value (min)
      CDtempCm = CIuarray(ulocate('microDtemp'));% D value temperature basis (C)
      CDtempKm = CDtempCm+Tref; % D value temperature basis (K)
      CZvalm =  CIuarray(ulocate('microZval'));   % Z value (C)
      CEam =  CIuarray(ulocate('microEa'));       % activation energy (kcal/mol)
      CNoin =  CIuarray(ulocate('microload'));     % microbial load
      CImdk = [CDvalm;CDtempCm;CZvalm;CEam;CNoin];  % microbial array
 
%     Cold Feed Quality Kinetics
      CDvalq = CIuarray(ulocate('qualDval'));    % D value (min)
      CDtempCq = CIuarray(ulocate('qualDtemp')); % D value temperature basis (C)
      CDtempKq = CDtempCq+Tref; % D value temperature basis (K)
      CZvalq = CIuarray(ulocate('qualZval'));     % Z value (C)
      CEaq = CIuarray(ulocate('qualEa'));       % activation energy (kcal/mole) 
      CQoin = CIuarray(ulocate('qualload'));     % quality load
      CIqual = [CDvalq;CDtempCq;CZvalq;CEaq;CQoin];  % quality array
 
%  ********************* MICROBIAL & QUALITY ASSESSMENT *********************

%     Pass microbial and quality assessment variables to exit streams
      HOmdk = HImdk;
      COmdk = CImdk;
      HOqual = HIqual;
      COqual = CIqual;
 
      if (mandq == 1) & (dlevel >= 3)

%       Bulk temperature for microbial assessment (K)
        HTbulk = (HOtempK+HItempK)/2;
 
%       Minimum Residence Time (min)
        Htmin = (length*plwidth*plgap)/HVq*60;

%       Temperature for quality assessment (K)
        HTqual = HItempK;

%       Maximum processing time = Minimum processing time (min)
        Htqmax = Htmin;

        if (HNoin > 0) 
 
%         Final microbial load
          HNout = mdeathk(HDvalm,HDtempKm,HZvalm,HEam,HNoin,HTbulk,Htmin);  
 
%         Log reduction in microbial load
          HDlogNred = log10(HNoin/HNout);

%         Exit stream microbial array
          HOmdk = HImdk;  
          HOmdk(size(HImdk,1)) = HNout;

        end 

        if (HQoin > 0)  
 
%         Final quality load
          HQout = qualdegr(HDvalq,HDtempKq,HZvalq,HEaq,HQoin,HTqual,Htqmax);
 
%         Log reduction in quality
          HDlogQred = log10(HQoin/HQout);
 
%         Exit stream quality array
          HOqual = HOqual;
          HOqual(size(HOqual,1)) = HQout;

        end  

%       Bulk temperature for microbial assessment (K)
        CTbulk = (CItempK+COtempK)/2;
 
%       Minimum Residence Time (min)
        Ctmin = (length*plwidth*plgap)/CVq*60;

%       Temperature for quality assessment (K)
        CTqual = COtempK;
 
%       Maximum processing time = Minimum processing time (min)
        Ctqmax = Ctmin;

        if (CNoin > 0)  
 
%         Final microbial load
          CNout = mdeathk(CDvalm,CDtempKm,CZvalm,CEam,CNoin,CTbulk,Ctmin); 
 
%         Log reduction in microbial load
          CDlogNred = log10(CNoin/CNout);

%         Exit stream microbial array
          COmdk = CImdk;  
          COmdk(size(CImdk,1)) = CNout;

        end 

        if (CQoin > 0)  
 
%         Final quality load
          CQout = qualdegr(CDvalq,CDtempKq,CZvalq,CEaq,CQoin,CTqual,Ctqmax);
 
%         Log reduction in quality
          CDlogQred = log10(CQoin/CQout);
 
%         Exit stream quality array
          COqual = COqual;
          COqual(size(COqual,1)) = CQout;

        end  

        endtxt = [10 'END OF UNIT (',Uname,') MICROBIAL & QUALITY ASSESSMENT' 10];
        disp(endtxt); 

      end  % microbial and quality assessment 

%     Update system output array
      HOuarray = [HOme;HOmdk;HOqual];
      COuarray = [COme;COmdk;COqual];

%  ****************** ECONOMIC EVALUATION VARIABLE ASSIGMENTS ***************

      if (econ >= 1)
        ECON = UOuarray(17:26); % economic analysis array 
      end  

%  ********************* SYSTEM OUTPUT ARRAY ASSIGNMENT *********************

%     Size stream info arrays to full system output u-arrays  
      HOuarray = ustd(HOuarray);
      COuarray = ustd(COuarray);

      sys =[HOuarray;COuarray];     
 
%  ****************************** PRINTED OUTPUT ****************************

  if flag == 9    % print only after a solution is found 

%   Print results to specified file name; 
%    call filename manager function rfile.m 
    fid = rfile;

%   Name of unit operation; system input streams; system output streams  
    UNAME = Uname; 
    SINAMES = str2mat(HIname,CIname);
    SONAMES = str2mat(HOname,COname);

%   System inputs array 
    SIARRAY = [HImass CImass]; 
%   System outputs array 
    SOARRAY = [HOmass COmass];

%   No intermediate streams
    SNNAMES = [];
    SNARRAY = []; 

%   Print mass balance results table
    msprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    if (dlevel >= 2) 
%     System inputs array  
      SIARRAY = [HIme CIme]; 
%     System outputs array 
      SOARRAY = [HOme COme];

%     Print energy balance results table
      eprint(UNAME,SINAMES,SONAMES,SNNAMES,SIARRAY,SOARRAY,SNARRAY); 

    end  
    
    % THIS IS WHERE THE CHANGE WAS MADE 8/21/01
    fid = fopen(fid,'a+');
    % END CHANGE 8/21/01
    
    if (dlevel >= 2) 
    fprintf(fid,'UNIT MODEL INPUTS\n'); 
    end     
    fprintf(fid,'     \n');

    if (dlevel >= 2)
      if (fconfig == 0) flowc = 'Co-current'; end
      if (fconfig == 1) flowc = 'Counter-current'; end
      fprintf(fid,'  Flow configuration                     %s \n',flowc); 
      fprintf(fid,'  Minimum Temperature Differential           %5.2f C\n',delTmin); 
      if (Thoutmin > 0) 
      fprintf(fid,'  Minimum Hot Stream Exit Temp.             %6.2f C\n',THOminC); 
      end 
      if (hottype == 1) 
      fprintf(fid,'  Minimum Hot Stream Exit Quality            %5.3f C\n',HOxmin);
      end  
      if (Tcoutmax > 0) & (Thoutmin <= 0)   
      fprintf(fid,'  Maximum Cold Stream Exit Temp.            %6.2f C\n',TCOmaxC); 
      end
      if (coldtype == 1)  
      fprintf(fid,'  Maximum Cold Stream Exit Quality           %5.3f \n',COxmax);
      end  
      fprintf(fid,'  Maximum Cold Side Operating Pressure   %5.3e MPa\n',Pcold);
 
    end

    if (dlevel >= 3) 
      fprintf(fid,'  Overall Heat Transfer Coeff.           %5.3e kJ/h-m^2-K\n',U); 
      fprintf(fid,'  Plate width                                %5.3f m\n',plwidth); 
      fprintf(fid,'  Plate length                               %5.3f m\n',pllength); 
      fprintf(fid,'  Plate gap                                 %5.4f m\n',plgap); 
    end  

    if (dlevel >= 4)
      fprintf(fid,'  Hot Stream Viscosity Model                 %5.1f m\n',Hmumodel); 
      fprintf(fid,'  Cold Stream Viscosity Model                %5.1f m\n',Cmumodel); 
    end 

    if (dlevel >= 2)  
    fprintf(fid,'     \n');
    fprintf(fid,'UNIT MODEL OUTPUTS\n');     
    fprintf(fid,'     \n');
    end  

    fprintf(fid,'  Heat Transfer Rate                     %5.3e kJ/h\n',(HI*HIh-HO*HOh));

    if (dlevel >= 3)  
      fprintf(fid,'  Heat Transfer Area                     %5.3e m^2\n',areahx);
      fprintf(fid,'  Number of Plates                           %5.1f \n',Nplates);
      fprintf(fid,'  Hot Stream Avg. Density                  %7.2f kg/m^3\n',Hrho);
      fprintf(fid,'  Hot Stream Avg. Volumetric Flow Rate       %5.2f m^3/h\n',HVq);
      fprintf(fid,'  Cold Stream Avg. Density                 %7.2f kg/m^3\n',Crho);
      fprintf(fid,'  Cold Stream Avg. Volumetric Flow Rate      %5.2f m^3/h\n',CVq);
    end 

    if (dlevel >=4) 
      fprintf(fid,'     \n');
      fprintf(fid,'  Hydraulic Diameter                     %5.3e m\n',Dh);
      fprintf(fid,'  Hot Stream Inlet:    \n');
      fprintf(fid,'    Density                              %5.3e kg/m^3\n',HIrho);
      fprintf(fid,'    Consistency Coefficent               %5.3e Pa-s^n\n',HIK);
      fprintf(fid,'    Flow Behavior Index                  %5.3e \n',HIn);
      fprintf(fid,'    Velocity                             %5.3e m/s\n',HIvel);
      fprintf(fid,'    Reynold Number                       %5.3e \n',HINRe);
      fprintf(fid,'  Hot Stream Outlet:    \n');
      fprintf(fid,'    Density                              %5.3e kg/m^3\n',HOrho);
      fprintf(fid,'    Consistency Coefficent               %5.3e Pa-s^n\n',HOK);
      fprintf(fid,'    Flow Behavior Index                  %5.3e \n',HOn);
      fprintf(fid,'    Velocity                             %5.3e m/s\n',HOvel);
      fprintf(fid,'    Reynold Number                       %5.3e \n',HONRe);
      fprintf(fid,'  Cold Stream Inlet:    \n');
      fprintf(fid,'    Density                              %5.3e kg/m^3\n',CIrho);
      fprintf(fid,'    Consistency Coefficent               %5.3e Pa-s^n\n',CIK);
      fprintf(fid,'    Flow Behavior Index                  %5.3e \n',CIn);
      fprintf(fid,'    Velocity                             %5.3e m/s\n',CIvel);
      fprintf(fid,'    Reynold Number                       %5.3e \n',CINRe);
      fprintf(fid,'  Cold Stream Outlet:    \n');
      fprintf(fid,'    Density                              %5.3e kg/m^3\n',COrho);
      fprintf(fid,'    Consistency Coefficent               %5.3e Pa-s^n\n',COK);
      fprintf(fid,'    Flow Behavior Index                  %5.3e \n',COn);
      fprintf(fid,'    Velocity                             %5.3e m/s\n',COvel);
      fprintf(fid,'    Reynold Number                       %5.3e \n',CONRe);
    
    endtxt = [10 'END OF (',Uname,') TRANSPORT RESULTS PRINT TO FILE' 10];  
    disp(endtxt);  

    end  

%   Microbial and Quality Assessment results printout 
    if (mandq == 1) & (dlevel >= 3)


    fprintf(fid,'     \n');

    if (HNoin > 0)

    fprintf(fid,' Hot Stream    \n');
    fprintf(fid,'  Microbial (Worst Case) Assessment:   \n');
    fprintf(fid,'   Minimum Residence Time                   %5.4f min \n',Htmin);
    fprintf(fid,'   Minimum Processing Temperature            %5.2f C \n',HTbulk-Tref);
    fprintf(fid,'   Log Reduction in Microbial Load        %5.2e \n',HDlogNred);

    end 

    if (CNoin > 0)

    fprintf(fid,'     \n');
    fprintf(fid,' Cold Stream    \n');
    fprintf(fid,'  Microbial (Worst Case) Assessment:   \n');
    fprintf(fid,'   Minimum Residence Time                   %5.4f min \n',Ctmin);
    fprintf(fid,'   Minimum Processing Temperature            %5.2f C \n',CTbulk-Tref);
    fprintf(fid,'   Log Reduction in Microbial Load        %5.2e \n',CDlogNred);

    end  

    if (HQoin > 0) 
    fprintf(fid,'     \n');
    fprintf(fid,' Hot Stream   \n');
    fprintf(fid,'  Quality (Worst Case) Assessment:   \n');
    fprintf(fid,'   Maximum Residence Time                   %5.4f min \n',Htqmax);
    fprintf(fid,'   Maximum Processing Temperature            %5.2f C \n',HTqual-Tref);
    fprintf(fid,'   Log Reduction in Quality Load          %5.2e \n',HDlogQred);

    end  

    if (CQoin > 0)
 
    fprintf(fid,'     \n');
    fprintf(fid,' Cold Stream   \n');
    fprintf(fid,'  Quality (Worst Case) Assessment:   \n');
    fprintf(fid,'   Maximum Residence Time                   %5.4f min \n',Ctqmax);
    fprintf(fid,'   Maximum Processing Temperature            %5.2f C \n',CTqual-Tref);
    fprintf(fid,'   Log Reduction in Quality Load          %5.2e \n',CDlogQred);
    fprintf(fid,'     \n');

    end 
 
%     System input and exit stream names for microbial and quality assessment 
      SIMQNAME = str2mat(HIname,CIname);
      SOMQNAME = str2mat(HOname,COname); 

%     System input microbial array
      SIMARRAY = [HOmdk CImdk];
%     System outputs microbial array
      HOmdkpr = [HOmdk;Htmin;HTbulk-Tref;];
      COmdkpr = [COmdk;Ctmin;CTbulk-Tref;];
      SOMARRAY = [HOmdkpr COmdkpr];

%     No intermediate streams
      SNMQNAME = [];
      SNMARRAY = [];
      SNQARRAY = [];

      mdkprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIMARRAY,SOMARRAY,SNMARRAY);

%     System input quality array
      SIQARRAY = [HIqual CIqual];
%     System outputs quality array
      HOqualpr = [HOqual;Htqmax;HTqual-Tref;];
      COqualpr = [COqual;Ctqmax;CTqual-Tref;];
      SOQARRAY = [HOqualpr COqualpr];

      qprint(Uname,SIMQNAME,SOMQNAME,SNMQNAME,SIQARRAY,SOQARRAY,SNQARRAY);
    end 

    if (econ >= 1) 
      if (hottype == 1) & (HIinit <=0)
        if (HIx > 0)
          S = HI;
          W = 0;
        else
          S = 0;
          W = HI;
        end
      elseif (coldtype == 1) & (CIinit <=0)
        if (CIx > 0)
          S = CI;
          W = 0;
        else
          S = 0;
          W = CI;
        end  
      else 
        S = 0;
        W = 0; 
      end  
      Uutil = [S W 0 0 0 0];  
      Usize = areahx;
      Unum = 1; 
      Uecon = ECON;  
      uofile(Uutil,Usize,Unum,Uecon,Uname);

    end  


    fprintf(fid,'     \n');
    fprintf(fid,'**********************************************************************************************\n');
    fprintf(fid,'**********************************************************************************************\n');
    fprintf(fid,'     \n');

    % THIS IS WHERE THE CHANGE WAS MADE 8/21/01
    fclose('all');
    % END CHANGE 8/21/01

    sys = []; 
  end % flag 9 print 

else
   % Otherwise, no need to return anything since this is continuous system
   sys = [];

end
 
