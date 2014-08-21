   function [sys,x0] = econs(t,x,u,flag,filesc,fileuo)

%  ECONS     Economic analysis model 
%
%            For use with FOODS-LIB models 
% 
%            Synopsis: S-function for economic analysis of entire food
%            process models. Performs stream costing, utility costing,
%            equipment costing, and time-value of money analysis.   
%
%            Call: econs(t,x,u,flag,filesc,fileuo)  
%           
%            Passed Arguments:
%             t,x,flag - SIMULINK S-function time, state vector, flag
%             u - SIMULINK S-function input vector; u-array carries stream
%                 and unit operation operating information
%             filesc - stream cost data file
%             fileuo - unit operation cost data file 
%
%            Returned Arguments:
%             sys - SIMULINK S-function variable: for flag = 0,
%                   sys = sizes of parameters and initial conditions.
%                   For flag = 3, sys = is an empty vector  
%
%            Reference: 
%             Peters, M.S. & Timmerhaus, K.D. 1991. "Plant Design and
%              Economics for Chemical Engineers," 4th Ed. McGraw-Hill, 
%              Inc., NY 

%  ************** SIMULINK INTERFACE & DATA FILES INFORMATION ***************

ECONins = 52;                     % Number of economic analysis inputs
 
if flag == 0
   % sys(3) is number of outputs
   % sys(4) is number of inputs
   sys = [0, 0, 0, ECONins, 0, 0];
   x0 = ones(size(u));
      
elseif flag == 3 | flag == 9
   % if flag=3, return no system outputs
   % if flag=9, print results 

%     Breakup of SIMULINK constructed u-array into economic analysis 
%     levels arrays    

      GENarray = u(1:4);             % General information array  
      OParray = u(5:12);             % Operating info array
      CAParray = u(13:25);           % Capital investment array 
      MANarray = u(26:38);           % Manufacturing-Cost ratios array  
      TVarray = u(39:52);            % Time-value of money array 

%   ******** GENERAL ANALYSIS INFORMATION VARIABLE ASSIGNMENTS **************

      elevel = GENarray(1);          % economic analysis level 
      hpery = GENarray(2);           % hours of opertion per year (h/y)
      curMSI = GENarray(3);          % current MSI  
      plife = GENarray(4);           % project life (y)  

%   *************** STREAM COST & VALUE VARIABLE ASSIGNMENTS ****************  

%     Open specified stream cost/value data file
      fid = fopen(filesc,'r');
      fid
%     Read in the data in specified format into a vector
      SCvector = fscanf(fid,'%g %g %g %g %s %g');
%     Close data file
      fclose(fid);

%     Find index for end of each stream cost/value record
      J = 0;
      for I = 1:length(SCvector)
        if (SCvector(I) < 0)
          J = J + 1;
          endstrm(J) = I;
        end
      end

      for I = 1:length(endstrm)    % Number of streams in econ file
        if I == 1
          SCtype(I) = SCvector(1);   % Stream type
          SC(I) = SCvector(2);       % Stream mass flow rate (kg/h)
          SCcost(I) = SCvector(3);   % Stream cost/value ($/kg)
          SCinfl(I) = SCvector(4);   % Stream inflation (dec percent)
          SCnamel = endstrm(I)-1-4;  % Stream name length and name
          SCname(I,1:SCnamel) = setstr(SCvector(5:(endstrm(I)-1)))';
        else
          SCtype(I) = SCvector(endstrm(I-1)+1);
          SC(I) = SCvector(endstrm(I-1)+2);
          SCcost(I) = SCvector(endstrm(I-1)+3);
          SCinfl(I) = SCvector(endstrm(I-1)+4);
          SCnamel = endstrm(I) -1 - endstrm(I-1) - 4;
          SCname(I,1:SCnamel) = setstr(SCvector((endstrm(I-1)+5):(endstrm(I)-1)))';
        end
      end  

      R = 0;
      rflows = 0;
      rcosts = 0;
      rinfl = 0;  
      rname = []; 
      P = 0; 
      pflows = 0; 
      pcosts = 0;  
      pinfl = 0;
      pname = [];   
      W = 0;
      wflows = 0; 
      wcosts = 0;  
      winfl = 0; 
      wname = []; 

%     Assignment of raw material, product, and waste stream cost/value info
      for I = 1:length(endstrm)
%       Raw materials information
        if (SCtype(I) == 1)
          R = R + 1;
          rflows(R) = SC(I);         % Flow rate of raw material (kg/h)
          rcosts(R) = SCcost(I);     % Cost of raw material ($/kg)
          rinfl(R) = SCinfl(I);      % Inflation rate for raw materials (dec percent))
          rname(R,:) = SCname(I,:);  % Name of raw material stream
%       Products generating income information
        elseif (SCtype(I) == 2)
          P = P + 1;
          pflows(P) = SC(I);         % Flow rate of product (kg/h)
          pcosts(P) = SCcost(I);     % Cost of product ($/kg)
          pinfl(P) = SCinfl(I);      % Inflation rate for product (dec percent))
          pname(P,:) = SCname(I,:);  % Name of product stream
%       Wastes information
        elseif (SCtype(I) == 3)
          W = W + 1;
          wflows(W) = SC(I);         % Flow rate of waste (kg/h)
          wcosts(W) = SCcost(I);     % Cost of waste ($/kg)
          winfl(W) = SCinfl(I);      % Inflation rate for waste (dec percent)
          wname(W,:) = SCname(I,:);  % Name of waste stream
        end
      end
 
      raws = R;                    % Number of raw material streams
      prods = P;                   % Number of product streams
      wastes = W;                  % Number of waste streams
 
%   ********** LEVEL 1 ECONOMIC ANALYSIS: STREAM COST & VALUE ***************

%     Simple raw material costs 
%     100% capacity flow rate (kg/y)
      rflowcap = hpery .* rflows;    
%     100% capacity annual cost ($/y)
      rcostcap = rcosts .* rflowcap;   
%     Total raw material costs ($/y)   
      rcosttot = sum(rcostcap);    

%     Simple income 
%     100% capacity flow rate (kg/y) 
      pflowcap = hpery .* pflows;  
%     100% capacity annual cost ($/y) 
      pcostcap = pcosts .* pflowcap; 
%     Total income   
      pcosttot = sum(pcostcap);    

%     Simple waste costs  
%     100% capacity flow rate (kg/y) 
      wflowcap = hpery .* wflows;  
%     100% capacity annual cost ($/y) 
      wcostcap = wcosts .* wflowcap; 
%     Total waste costs   
      wcosttot = sum(wcostcap);  

%   ************* LEVEL 1 PRINT RESULTS: STREAM COST & VALUE ****************

    if flag == 9
%   Print to results to specified file name
    global fileea;  
    filehere = isstr(fileea);  
    if (filehere == 0)
      % CHANGE MADE 4/21/01; MADE STRING INPUT POSSIBLE
      fileea = input('Enter filename for economic analysis results:   ','s');
      % END CHANGE 4/21/01
      fid = fopen(fileea,'a+');
    else
      fid = fopen(fileea,'a+');
    end

%   Limiting or sizing all stream names to 15 characters
    rname = namesize(rname,15);
    pname = namesize(pname,15);
    wname = namesize(wname,15);

    fprintf(fid,'     \n');
    fprintf(fid,'ECONOMIC REPORT\n');
    fprintf(fid,'     \n');
    fprintf(fid,'SYSTEM STREAM COST ANALYSIS - FULL CAPACITY\n');
    fprintf(fid,'     \n');
    fprintf(fid,'BASIS: Operating Time - %5.0f h/y \n',hpery);  
    fprintf(fid,'     \n');
    fprintf(fid,'Item                            Rate          Value/         Annual          Annual\n');  
    fprintf(fid,'                                               Cost        Value/Cost      Value/Cost\n');
    fprintf(fid,'                              (units/y)      ($/unit)         ($/y)           ($/y)\n');  
    fprintf(fid,'--------------------------------------------------------------------------------------\n');   
    fprintf(fid,'PROCESS PRODUCTS (kg/y)  \n');  
    for P = 1:prods  
    fprintf(fid,' %s           %12.1f        %2.4f   %13.0f    \n', ... 
            pname(P,:),pflowcap(P),pcosts(P),pcostcap(P));    
    end 
    fprintf(fid,'TOTAL INCOME                                                         %14.0f \n',pcosttot); 
    fprintf(fid,'--------------------------------------------------------------------------------------\n');   
    fprintf(fid,'RAW MATERIALS (kg/y)  \n');           
    for R = 1:raws   
    fprintf(fid,' %s           %12.1f        %2.4f   %13.0f    \n', ... 
            rname(R,:),rflowcap(R),rcosts(R),rcostcap(R));
    end   
    fprintf(fid,'TOTAL RAW MATERIALS                                                  %14.0f \n',rcosttot);
    fprintf(fid,'--------------------------------------------------------------------------------------\n');
    fprintf(fid,'WASTE STREAMS (kg/y)  \n');           
    for W = 1:wastes     
    fprintf(fid,' %s           %12.1f        %2.4f   %13.0f    \n', ... 
            wname(W,:),wflowcap(W),wcosts(W),wcostcap(W));
    end   
    fprintf(fid,'TOTAL WASTE CHARGES                                                  %14.0f \n',wcosttot);
    fprintf(fid,'--------------------------------------------------------------------------------------\n');
    if elevel == 1
    fprintf(fid,'GROSS EARNINGS                                                       %14.0f \n',pcosttot-rcosttot-wcosttot);  
    fprintf(fid,'--------------------------------------------------------------------------------------\n');
    end 
    end  

%   ****************** UNIT OPERATION VARIABLE ASSIGNMENT *******************
 
      if (elevel >= 2) 

%     Open specified unit operation cost data file
      fid = fopen(fileuo);
%     Read in the data in specified format into a vector
      UOvector = fscanf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %s %g');
%     Close data file
      fclose(fid);

%     Find index for end of each unit operation cost record
      U = 0;
      for I = 1:length(UOvector)
        if (UOvector(I) < 0)
          U = U + 1;
          enduo(U) = I;
        end
      end  

      for I = 1:length(enduo)    % Number of unit operations in file
        if I == 1
%         Utilities demand:  steam demand (kg/h), water demand (kg/h),
%         electricity demand (kW), natural gas demand (m^3/h),
%         fuel oil demand (L/h), compressed air demand (m^3/h)
          uutil(I,1:6) = UOvector(1:6)';
          usize(I) = UOvector(7);        % unit op calculated size
          uopnum(I) = UOvector(8);       % quantity to be purchased
          ubaseMSI(I) = UOvector(9);   % base MSI
          ubassize(I) = UOvector(10);  % base size
          ubascost(I) = UOvector(11);  % base cost ($)
          uscaling(I) = UOvector(12);  % scaling index for size
          uminsize(I) = UOvector(13);  % minimum size for scaling range
          umaxsize(I) = UOvector(14);  % maximum size for scaling range
          usalvage(I) = UOvector(15);  % salvage value (decimal percent
                                              %  of initial cost)
          ueqlife(I) = UOvector(16);   % equipment life (y)
          uoplabor(I) = UOvector(17);  % operational labor
                                              %  (man-h/h operation)
          umlabor(I) = UOvector(18);   % maintenance labor
                                              %  (man-h/h operation)
          UOnamel = enduo(I)-1-18;     % Unit Op name length and name
          uname(I,1:UOnamel) = setstr(UOvector(19:(enduo(I)-1)))';
        else
          uutil(I,1:6) = UOvector((enduo(I-1)+1):(enduo(I-1)+6))';
          usize(I) = UOvector(enduo(I-1)+7);
          uopnum(I) = UOvector(enduo(I-1)+8);
          ubaseMSI(I) = UOvector(enduo(I-1)+9);
          ubassize(I) = UOvector(enduo(I-1)+10);
          ubascost(I) = UOvector(enduo(I-1)+11);
          uscaling(I) = UOvector(enduo(I-1)+12);
          uminsize(I) = UOvector(enduo(I-1)+13);
          umaxsize(I) = UOvector(enduo(I-1)+14);
          usalvage(I) = UOvector(enduo(I-1)+15);
          ueqlife(I) = UOvector(enduo(I-1)+16);
          uoplabor(I) = UOvector(enduo(I-1)+17);
          umlabor(I) = UOvector(enduo(I-1)+18);
          UOnamel = enduo(I)-1-enduo(I-1)-18;
          uname(I,1:UOnamel) = setstr(UOvector((enduo(I-1)+19):(enduo(I)-1)))';
        end
      end

      utils = 6;                       % Number of utilities 
      unitops = U;                     % Number of unit operations

      cutil(1) = OParray(1);           % steam cost ($/kg)
      cutil(2) = OParray(2);           % water cost ($/kg)
      cutil(3) = OParray(3);           % electricity cost ($/kW-h)
      cutil(4) = OParray(4);           % natural gas cost ($/m^3)
      cutil(5) = OParray(5);           % fuel oil cost ($/L)
      cutil(6) = OParray(6);           % compressed air cost ($/m^3)
      coplabor = OParray(7);           % operating labor cost ($/man-h)
      cmlabor =  OParray(8);           % maintenance labor cost ($/man-h)

%  *************** LEVEL 2 ECONOMIC ANALYSIS: STREAM & UTILITY **************

%   Simple utility costs
%   100% capacity utility demands (units/y); each utility demand summed 
%    over all unit ops    
    if (U == 1) 
      utild = uutil .* hpery;  
    else 
      utild = sum(uutil) .* hpery;    % NEED TO CHECK WITH U > 1 
    end    
 
%   100% capacity utility costs ($/y) for each utility    
    utilc = utild .* cutil; 
%   100% capacity total utility cost ($/y)  
    utilctot = sum(utilc); 

%   Simple labor demand (man-h/y) summed over all unit operations 
    oplabtot = sum(uoplabor)*hpery;
    mlabtot = sum(umlabor)*hpery;
%   100% capacity labor cost ($/y)
    coplabt = oplabtot*coplabor;
    cmlabt = mlabtot*cmlabor;


%  ***************** LEVEL 2 PRINT RESULTS: STREAM & UTILITY ****************

    if flag == 9
%   Continue full capacity table 
    fid = fileea;   
    fprintf(fid,'UTILITIES  \n'); 
    fprintf(fid,' steam (kg/y)              %12.1f        %2.4f   %13.0f    \n',utild(1),cutil(1),utilc(1));
    fprintf(fid,' water (kg/y)              %12.1f        %2.4f   %13.0f    \n',utild(2),cutil(2),utilc(2));
    fprintf(fid,' electric (kW-h/y)         %12.1f        %2.4f   %13.0f    \n',utild(3),cutil(3),utilc(3));
    fprintf(fid,' natural gas               %12.1f        %2.4f   %13.0f    \n',utild(4),cutil(4),utilc(4));
    fprintf(fid,' fuel oil                  %12.1f        %2.4f   %13.0f    \n',utild(5),cutil(5),utilc(5));
    fprintf(fid,' compressed air            %12.1f        %2.4f   %13.0f    \n',utild(6),cutil(6),utilc(6));
    fprintf(fid,'TOTAL UTILITIES                                                      %14.0f \n',utilctot);                                                            
    fprintf(fid,'--------------------------------------------------------------------------------------\n');   
    fprintf(fid,'LABOR  \n'); 
    fprintf(fid,' operating (man-h/y)       %12.1f       %2.4f   %13.0f    \n',oplabtot,coplabor,coplabt);                                   
    fprintf(fid,' maintenance (man-h/y)     %12.1f       %2.4f   %13.0f    \n',mlabtot,cmlabor,cmlabt);           
    fprintf(fid,'TOTAL LABOR                                                          %14.0f \n',coplabt+cmlabt);

    fprintf(fid,'--------------------------------------------------------------------------------------\n');   
    fprintf(fid,'TOTAL GROSS EARNINGS                                                 %14.0f \n',pcosttot-rcosttot-wcosttot-utilctot-coplabt-cmlabt);  
    fprintf(fid,'--------------------------------------------------------------------------------------\n');   
    fprintf(fid,' \n');   

%   Limiting or sizing all stream names to 15 characters
    uname = namesize(uname,15); 

    fprintf(fid,' \n');   
    fprintf(fid,'UTILITIES DEMAND BY UNIT OPERATION: \n'); 
    fprintf(fid,' \n');   
    fprintf(fid,'BASIS: Operating Time - %5.0f h/y \n',hpery); 
    fprintf(fid,' \n');   
    fprintf(fid,'Unit Operation             Steam           Water           Elec.         Natural         Fuel         Compressed   \n');
    fprintf(fid,'                                                                           Gas            Oil             Air     \n');
    fprintf(fid,'                          (kg/y)           (kg/y)        (kW-h/y)         \n');
    fprintf(fid,'---------------------------------------------------------------------------------------------------------------------  \n');
    for U = 1:unitops
    fprintf(fid,'%s    %12.1f    %12.1f    %12.1f    %12.1f    %12.1f    %12.1f  \n', ...
    uname(U,:),uutil(U,1)*hpery,uutil(U,2)*hpery,uutil(U,3)*hpery,uutil(U,4)*hpery,uutil(U,5)*hpery,uutil(U,6)*hpery);
    end  
    fprintf(fid,'---------------------------------------------------------------------------------------------------------------------  \n');       
    fprintf(fid,'TOTAL DEMANDS:   %14.1f  %14.1f  %14.1f  %14.1f  %14.1f  %14.1f  \n', ...
     utild(1),utild(2),utild(3),utild(4),utild(5),utild(6));  
    fprintf(fid,'---------------------------------------------------------------------------------------------------------------------  \n');
    fprintf(fid,' \n');  

    fprintf(fid,' \n');   
    fprintf(fid,'LABOR DEMAND BY UNIT OPERATION: \n'); 
    fprintf(fid,' \n');   
    fprintf(fid,'BASIS: Operating Time - %5.0f h/y \n',hpery); 
    fprintf(fid,' \n');   
    fprintf(fid,'Unit Operation          Operating       Maintenance  \n');   
    fprintf(fid,'                          Labor            Labor  \n');          
    fprintf(fid,'                        (man-h/y)        (man-h/y)  \n');   
    fprintf(fid,'-------------------------------------------------------  \n');
    for U = 1:unitops
    fprintf(fid,'%s          %7.1f         %7.1f  \n', ... 
            uname(U,:),uoplabor(U)*hpery,umlabor(U)*hpery);
    end  
    fprintf(fid,'-------------------------------------------------------  \n');     
    fprintf(fid,'TOTAL DEMANDS:         %9.1f       %9.1f  \n',oplabtot,mlabtot);
    fprintf(fid,'-------------------------------------------------------  \n');
    fprintf(fid,' \n');

    end % flag = 9 

    end % LEVEL 2

%   ******* LEVEL 3 ECONOMIC ANALYSIS: STREAM, UTILITY, EQUIPMENT **********

    if (elevel >= 3)

%   Cost of each unit operation ($)  
    uopcost = uopnum .* (curMSI ./ ubaseMSI) .* ubascost .* ... 
               (usize ./ ubassize) .^ uscaling;       

%   Number of five year blocks  
    T = ceil(plife/5);  
%   Extension of arrays to include full 5 year blocks 
    extend = T*5 - plife;  
%   Index for year zero formation in arrays  
    yzero = plife + extend + 1;    
%   Number of capital investment factors
    cinvinfo = length(CAParray);  
%   Initialize investment statement array  
    cinvest = zeros(cinvinfo+1,yzero);     
%   Total equipment cost in year 0 ($)  
    cinvest(1,yzero) = sum(uopcost);

%   ********** LEVEL 3 PRINT RESULTS: STREAM, UTILITY, EQUIPMENT ************

    if flag == 9 
    fprintf(fid,' \n');   
    fprintf(fid,'PURCHASED EQUIPMENT COSTS: \n');  
    fprintf(fid,' \n');   
    fprintf(fid,'BASIS: Current MSI - %5.0f \n',curMSI);  
    fprintf(fid,' \n');   
    fprintf(fid,'Unit Operation     Size     Quant.     Base     Base      Base      Low Range   High Range   Scaling      Total  \n');  
    fprintf(fid,'                                       MSI      Size    Cost ($)      Size         Size       Index      Cost ($)  \n');   
    fprintf(fid,'---------------------------------------------------------------------------------------------------------------------  \n');  
    for U = 1:unitops
    fprintf(fid,'%s  %8.2f     %4.1f    %5.0f    %8.3f %9.0f   %8.2f     %8.2f      %4.3f    %9.0f  \n', ...
    uname(U,:),usize(U),uopnum(U),ubaseMSI(U),ubassize(U),ubascost(U),uminsize(U),umaxsize(U),uscaling(U),uopcost(U));  
    end  
    fprintf(fid,'---------------------------------------------------------------------------------------------------------------------  \n');  
    fprintf(fid,'TOTAL COSTS:                                                                                         %12.0f \n',cinvest(1,yzero)); 
    fprintf(fid,'---------------------------------------------------------------------------------------------------------------------  \n');  
    fprintf(fid,' \n');   
    end % flag = 9  

    end % LEVEL 3

%   *************** TIME-VALUE OF MONEY VARIABLE ASSIGNMENTS ***************

    if (elevel >= 4) 
 
%   Capital Investment Factors
    cinvrat(1) = 1;                  % Purchased equipment ratio factor
    cinvrat(2:(cinvinfo+1)) = CAParray(1:cinvinfo);

%   Maufacturing-Cost Ratio Factors
    cmaninfo = length(MANarray);     % Number of manufacturing-cost factors
    cmanurat(1:cmaninfo) = MANarray(1:cmaninfo);

%   Time-Value of Money Information
    capacity(1) = TVarray(1)/100;     % first year capacity (d%)
    capacity(plife) = TVarray(2)/100; % last year capacity (d%)
    interest = TVarray(3)/100;        % interest rate (d%)
    inflat = TVarray(4)/100;          % inflation rate (d%)
    utilinfl(1) = TVarray(5)/100;     % steam inflation rate (d%)
    utilinfl(2) = TVarray(6)/100;     % water inflation rate (d%)
    utilinfl(3) = TVarray(7)/100;     % electricity inflation rate (d%)
    utilinfl(4) = TVarray(8)/100;     % natural gas inflation rate (d%)
    utilinfl(5) = TVarray(9)/100;    % fuel oil inflation rate (d%)
    utilinfl(6) = TVarray(10)/100;    % compressed air inflation rate (d%)
    oplinfl = TVarray(11)/100;     % Operating labor inflation (d%)
    mnlinfl = TVarray(12)/100;    % Maintenance labor inflation (d%)
    dorcomp = TVarray(13);         % Discrete (0) or Compounded (1) Models for
                                  %  inflation and interest
    corptax = TVarray(14)/100;         % corporate tax rate (d%)



%   ********** LEVEL 4 ECONOMIC ANALYSIS: TIME-VALUE OF MONEY ***************

    cmanu = zeros(20,plife);  
    csubtots = zeros(20,yzero); 

%   Yearly capacity calculation - linear distribution over Y years   
    capacity = linspace(capacity(1),capacity(plife),plife);    

%   Yearly income and direct production costs  
%   Raw material costs for each year of project  
%   Actual capacity in year Y (kg/y)  
    rawcap = rflowcap' * capacity;
%   Raw material cost in year Y ($/kg) 
%   rcosty = zeros(R,plife+extend); 
    rcosty = 0;
    for Y = 1:plife  
      for R = 1:raws  
        rcosty(R,Y) = inflate(rcosts(R),rinfl(R),Y,dorcomp); 
      end  
    end  
%   Actual capacity annual cost in year Y ($/y) 
    Rcosts = rcosty .* rawcap; 
%   Total raw material costs in year Y ($/y) 
    if (raws == 1) 
      cmanu(1,:) = Rcosts;      
    else 
      cmanu(1,:) = sum(Rcosts);        
    end   

%   Waste discharge costs for each year of project  
%   Actual capacity in year Y (kg/y)  
    wastecap = wflowcap' * capacity;   
%   Waste discharge cost in year Y ($/kg)   
    wcosty = 0;
    for Y = 1:plife  
      for W = 1:wastes    
        wcosty(W,Y) = inflate(wcosts(W),winfl(W),Y,dorcomp);
      end   
    end   
%   Actual capacity annual cost in year Y ($/y) 
    Wcosts = wcosty .* wastecap;
%   Total waste discharge costs in year Y ($/y)       
    if (wastes == 1)  
      cmanu(15,:) = Wcosts;   
    elseif (W > 1)  
      cmanu(15,:) = sum(Wcosts);    
    end    

%   Income for each year of project
%   Actual capacity in year Y (kg/y)
    prodcap = pflowcap' * capacity;   
%   Product income in year Y ($/kg)
    pcosty = 0;
    for Y = 1:plife
      for P = 1:prods  
        pcosty(P,Y) = inflate(pcosts(P),pinfl(P),Y,dorcomp);
      end
    end  
%   Actual capacity product income in year Y ($/y)
    Pcosts = pcosty .* prodcap;   
%   Total product income in year Y ($/y)
    if (prods == 1)
      csubtots(11,1:plife) = Pcosts;   
    elseif (P > 1)  
      csubtots(11,1:plife) = sum(Pcosts);  
    end

%   Utility costs for each year of project
%   Actual demand in year Y (units/y)
    utildem = utild' * capacity;     
%   Utility costs in year Y ($/kg)
    for Y = 1:plife
      for UT = 1:utils  
        utcosty(UT,Y) = inflate(cutil(UT),utilinfl(UT),Y,dorcomp);% NEED TO CHECK U > 1
      end
    end 
%   Actual utility cost in year Y ($/y)
    UTcosts = utcosty .* utildem;   
%   Total cost of each utility in year Y ($/y)
    cmanu(10,:) = sum(UTcosts);  

%   Depreciation - Straight-Line Method - for equipment and
%    associated installation, intrumentation, piping, electrical
%   Initialize arrays   
    deprtot = zeros(1,plife);
    salvtot = zeros(1,plife);  
    for U = 1:unitops  
%     Number of times a piece of equipment must be replaced
%      during project life
      noequip = ceil(plife/ueqlife(U));    
      for N = 1:noequip  
%       Depreciable investment assoc. with equipment  
        deprinv = 0; 
%       Depreciation each year for 1st equipment life  
        if (N == 1)  
%         New equipment cost 
          newequip = uopcost(U);           
%         Salvage value of new equipment   
          newsalv = usalvage(U) * uopcost(U);  
          for I = 2:5  
            deprinv = deprinv + newequip*cinvrat(I);  
          end   
          for Y = 1:ueqlife(U)  
             cmanu(12,Y) = cmanu(12,Y)+(newequip-newsalv+deprinv)/ueqlife(U); 
          end   
%         First equipment life salvage value 
          csubtots(5,ueqlife(U)) = csubtots(5,ueqlife(U)) + newsalv; 
%       For intermediate equipment life             
        elseif (N > 1) & (N < noequip)
%         New equipment cost  
          newequip = inflate(uopcost(U),inflat,(N-1)*ueqlife(U),dorcomp); 
%         Sum equipment purchased in year of replacement 
          cinvest(1,(N-1)*ueqlife(U)) = cinvest(1,(N-1)*ueqlife(U)) + newequip; 
%         New salvage value  
          newsalv = usalvage(U) * newequip; 
          for I = 2:5 
            deprinv = deprinv + newequip*cinvrat(I);  
            cinvest(I,(N-1)*ueqlife(U)) = cinvest(I,(N-1)*ueqlife(U))+newequip*cinvrat(I);   
          end  
          for Y = ((N-1)*ueqlife(U)+1):N*ueqlife(U) 
             cmanu(12,Y) = cmanu(12,Y) + (newequip-newsalv+deprinv)/ueqlife(U);  
          end  
%         Second through second-to last equipment life salvage value
          csubtots(5,N*ueqlife) = csubtots(5,N*ueqlife) + newsalv;  
%       For last equipment life
        elseif (N == noequip) 
%         Cost of new equipment in year Y
          newequip = inflate(uopcost(U),inflat,(N-1)*ueqlife(U),dorcomp);  
%         Sum equipment purchased in year of replacement
          cinvest(1,(N-1)*ueqlife(U)) = cinvest(1,(N-1)*ueqlife(U))+newequip;  
%         New salvage value  
          newsalv = usalvage(U) * newequip;  
          for I = 2:5
            deprinv = deprinv + newequip*cinvrat(I);  
            cinvest(I,(N-1)*ueqlife(U)) = cinvest(I,(N-1)*ueqlife(U))+newequip*cinvrat(I);  
          end
          sumdepr = 0; 
          for Y = ((N-1)*ueqlife(U)+1),plife
            cmanu(12,Y) = cmanu(12,Y) + (newequip-newsalv+deprinv)/ueqlife(U); 
            sumdepr = sumdepr + (newequip-newsalv)/ueqlife(U); 
          end   
%         Last equipment salvage value
          csubtots(5,plife) = csubtots(5,plife) + newequip - sumdepr;  
        end   
      end
    end   

%   Determine Direct & Indirect Costs (year 0)
    for I = 2:cinvinfo+1  
      cinvest(I,yzero) = cinvest(1,yzero) .* cinvrat(I);       
    end    

%   Depreciation of Yard, Service Facilities, Buildings
    deprinv = sum(cinvest(6:8,yzero));  
    deprtot = deprtot + deprinv/plife; 

%   Direct Costs Total (DCI); Investment made in year zero
    csubtots(1,:) = sum(cinvest(1:9,:));  
 
%   Indirect Costs Total (ICI); Investment made in year zero
    csubtots(2,:) = sum(cinvest(10:13,:));  
 
%   Fixed-capital Investment (FCI = DCI + ICI);
%   Investment made in year zero
    csubtots(3,:) = csubtots(1,:) + csubtots(2,:);      
 
%   Total Capital Investment  (FCI + working capital)
%   Investment made in year zero
    csubtots(4,:) = csubtots(3,:) + cinvest(14,:);  

%  ************* LEVEL 4 PRINTED RESULTS: TIME-VALUE OF MONEY ***************

    if flag == 9

%   Printed tables for estimated capital-investment  
    T = ceil(plife/5); 
    for Y = 1:T 
      if (Y == 1)       
    fprintf(fid,' \n');      
    fprintf(fid,'ESTIMATED CAPITAL-INVESTMENT STATEMENT (YEAR 0 - %1.0f) \n',Y+4);  
    fprintf(fid,' \n');  
    fprintf(fid,'Item                                     Ratio      Year %1.0f      Year %1.0f      Year %1.0f      Year %1.0f      Year %1.0f      Year %1.0f \n',Y-1,Y,Y+1,Y+2,Y+3,Y+4); 
    fprintf(fid,'                                         Factor       ($)         ($)         ($)         ($)         ($)         ($) \n');   
    fprintf(fid,'------------------------------------------------------------------------------------------------------------------------ \n');   
    fprintf(fid,'  Purchased Equipment (delivered), PE    %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(1),cinvest(1,yzero),cinvest(1,Y:Y+4));    
    fprintf(fid,'  Purchased Equipment Installation       %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(2),cinvest(2,yzero),cinvest(2,Y:Y+4));                      
    fprintf(fid,'  Instrumentation and Controls           %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(3),cinvest(3,yzero),cinvest(3,Y:Y+4));                      
    fprintf(fid,'  Piping (installed)                     %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(4),cinvest(4,yzero),cinvest(4,Y:Y+4));                      
    fprintf(fid,'  Electrical (installed)                 %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(5),cinvest(5,yzero),cinvest(5,Y:Y+4));                      
    fprintf(fid,'  Buildings (including services)         %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(6),cinvest(6,yzero),cinvest(6,Y:Y+4));                      
    fprintf(fid,'  Yard Improvements                      %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(7),cinvest(7,yzero),cinvest(7,Y:Y+4));                      
    fprintf(fid,'  Service Facilities                     %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(8),cinvest(8,yzero),cinvest(8,Y:Y+4));                      
    fprintf(fid,'  Land                                   %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(9),cinvest(9,yzero),cinvest(9,Y:Y+4));                      
    fprintf(fid,' DIRECT CAPITAL INVESTMENT                       %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(1,yzero),csubtots(1,Y:Y+4)); 
    fprintf(fid,'  Engineering and Supervision            %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(10),cinvest(10,yzero),cinvest(10,Y:Y+4));                      
    fprintf(fid,'  Construction Expenses                  %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(11),cinvest(11,yzero),cinvest(11,Y:Y+4));                      
    fprintf(fid,'  Contractor Fee                         %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(12),cinvest(12,yzero),cinvest(12,Y:Y+4));                      
    fprintf(fid,'  Contingency                            %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(13),cinvest(13,yzero),cinvest(13,Y:Y+4));                      
    fprintf(fid,' INDIRECT CAPITAL INVESTMENT                     %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(2,yzero),csubtots(2,Y:Y+4)); 
    fprintf(fid,' FIXED-CAPITAL INVESTMENT,FCI                    %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(3,yzero),csubtots(3,Y:Y+4)); 
    fprintf(fid,' WORKING CAPITAL                         %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(14),cinvest(14,yzero),cinvest(14,Y:Y+4));                      
    fprintf(fid,'TOTAL CAPITAL INVESTMENT                         %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(4,yzero),csubtots(4,Y:Y+4)); 
    fprintf(fid,'------------------------------------------------------------------------------------------------------------------------ \n');   
    fprintf(fid,' \n');   

      else  
    fprintf(fid,' \n');
    fprintf(fid,'ESTIMATED CAPITAL-INVESTMENT STATEMENT (YEAR %2.0f - %2.0f) \n',5*Y-4,5*Y);
    fprintf(fid,' \n');
    fprintf(fid,'Item                                     Ratio      Year %1.0f      Year %1.0f      Year %1.0f      Year %1.0f      Year %1.0f   \n',5*Y-4,5*Y-3,5*Y-2,5*Y-1,5*Y);
    fprintf(fid,'                                         Factor       ($)          ($)         ($)         ($)         ($)  \n');
    fprintf(fid,'------------------------------------------------------------------------------------------------------------ \n');
    fprintf(fid,'  Purchased Equipment (delivered), PE    %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(1),cinvest(1,5*Y-4:5*Y));
    fprintf(fid,'  Purchased Equipment Installation       %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(2),cinvest(2,5*Y-4:5*Y));
    fprintf(fid,'  Instrumentation and Controls           %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(3),cinvest(3,5*Y-4:5*Y));
    fprintf(fid,'  Piping (installed)                     %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(4),cinvest(4,5*Y-4:5*Y));
    fprintf(fid,'  Electrical (installed)                 %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(5),cinvest(5,5*Y-4:5*Y));
    fprintf(fid,'  Buildings (including services)         %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(6),cinvest(6,5*Y-4:5*Y));
    fprintf(fid,'  Yard Improvements                      %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(7),cinvest(7,5*Y-4:5*Y));
    fprintf(fid,'  Service Facilities                     %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(8),cinvest(8,5*Y-4:5*Y));
    fprintf(fid,'  Land                                   %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(9),cinvest(9,5*Y-4:5*Y));
    fprintf(fid,' DIRECT CAPITAL INVESTMENT                       %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(1,5*Y-4:5*Y));
    fprintf(fid,'  Engineering and Supervision            %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(10),cinvest(10,5*Y-4:5*Y));
    fprintf(fid,'  Construction Expenses                  %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(11),cinvest(11,5*Y-4:5*Y));
    fprintf(fid,'  Contractor Fee                         %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(12),cinvest(12,5*Y-4:5*Y));
    fprintf(fid,'  Contingency                            %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(13),cinvest(13,5*Y-4:5*Y));
    fprintf(fid,' INDIRECT CAPITAL INVESTMENT                     %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(2,5*Y-4:5*Y));
    fprintf(fid,' FIXED-CAPITAL INVESTMENT,FCI                    %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(3,5*Y-4:5*Y));
    fprintf(fid,' WORKING CAPITAL                         %4.3f   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvrat(14),cinvest(14,5*Y-4:5*Y));
    fprintf(fid,'TOTAL CAPITAL INVESTMENT                         %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(4,5*Y-4:5*Y));
    fprintf(fid,'------------------------------------------------------------------------------------------------------------ \n');
    fprintf(fid,' \n');
      end  
    end  
    end % flag = 9
       
%   Operating and Maintenance Labor Cost
%   Actual labor demand in year Y 
    oplabcap = oplabtot .* capacity; 
    mlabcap = mlabtot .* capacity; 
%   Actual cost in year Y ($/man-h)
    for Y = 1:plife 
      olcosty(Y) = inflate(coplabor,oplinfl,Y,dorcomp);  
      mlcosty(Y) = inflate(cmlabor,mnlinfl,Y,dorcomp);  
    end  
%   Actual annual labor costs in year Y ($/y)
    cmanu(2,:) = oplabcap .* olcosty; 
    cmanu(5,:) = mlabcap .* mlcosty; 

%   Operating Supervision Cost ($/y) - ratio of oper. labor
    cmanu(3,:) = cmanurat(1) .* cmanu(2,:);

%   Maintenance Materials Cost ($/y) - ratio of purchased equipment
    for Y = 1:plife 
      if (Y == 1) 
        cmanu(6,Y) = cmanurat(2)*cinvest(1,yzero);  
      else  
        cmanu(6,Y) = inflate(cmanu(6,1),inflat,Y,dorcomp);  
      end 
    end  

%   Maintenance Supervision Cost ($/y) - ratio of maint. labor
    cmanu(7,:) = cmanurat(3) .* cmanu(5,:);  
   
%   Total Maintenance & Repair Cost
    cmanu(4,:) = cmanu(5,:)+cmanu(6,:)+cmanu(7,:);   

%   Operating Supplies Cost ($/y) - ratio of maint. items
    cmanu(8,:) = cmanurat(4) .* cmanu(4,:);  

%   Lab Charges Cost ($/y) - ratio of oper. labor
    cmanu(9,:) = cmanurat(9) .* cmanu(2,:);  

    for Y = 1:plife   
%     Local Taxes ($/y) - ratio of FCI
%     Insurance ($/y) - ratio of FCI
%     Financing ($/y) - ratio of FCI
      if (Y == 1) 
        cmanu(13,Y) = cmanurat(7) .* csubtots(3,yzero);  
        cmanu(14,Y) = cmanurat(8) .* csubtots(3,yzero);  
        cmanu(20,Y) = cmanurat(13) .* csubtots(3,yzero);  
      else   
        cmanu(13,Y) = inflate(cmanu(13,1),inflat,Y,dorcomp);  
        cmanu(14,Y) = inflate(cmanu(14,1),inflat,Y,dorcomp);  
        cmanu(20,Y) = inflate(cmanu(20,1),inflat,Y,dorcomp);  
      end 
    end   
 
%   Plant Overhead Costs ($/y) - ratio of oper. & main.
    cmanu(16,:) = cmanurat(9) .* (cmanu(2,:)+cmanu(3,:)+cmanu(4,:));   

%   Total Product Cost (TPC) ($/y)  
    subtotpn = sum(cmanu(1:10,:)) + sum(cmanu(13:16,:)); 
    subtotpd = 1 - cmanurat(6) - cmanurat(10) - cmanurat(11) - cmanurat(12);
    csubtots(10,1:plife) = subtotpn/subtotpd;   

%   Patent & Royalties Cost ($/y) - ratio of TPC
    cmanu(11,:) = cmanurat(6) .* csubtots(10,1:plife);  

%   Administrative Cost($/y) - ratio of TPC
    cmanu(17,:) = cmanurat(10) .* csubtots(10,1:plife);  

%   Distribution & Marketing Cost ($/y) - ratio of TPC
    cmanu(18,:) = cmanurat(11) .* csubtots(10,1:plife);  
 
%   Research & Development Cost ($/y) - ratio of TPC
    cmanu(19,:) = cmanurat(12) .* csubtots(10,1:plife);  

%   Direct Production Costs ($/y)
    csubtots(6,1:plife) = sum(cmanu(1:11,:));   
    csubtots(6,1:plife) = csubtots(6,1:plife)-cmanu(4,:)+cmanu(15,:);  
 
%   Fixed Charges ($/y)
    csubtots(7,1:plife) = sum(cmanu(13:14,:));  
 
%   Manufacturing Costs ($/y)
    csubtots(8,1:plife) = csubtots(6,1:plife)+csubtots(7,1:plife)+cmanu(16,1:plife);   
 
%   General Expenses ($/y)
    csubtots(9,1:plife) = sum(cmanu(17:20,:));  
 
%   Annual operating income ($/y)
    csubtots(12,1:plife) = csubtots(11,1:plife) - csubtots(10,1:plife);  

%   Income before taxes ($/y)
    csubtots(13,1:plife) = csubtots(12,1:plife) - cmanu(12,:);  
    untaxed = 0; 
    for Y = 1:plife  
      if (csubtots(13,Y) < 0) 
        untaxed = untaxed + csubtots(13,Y);  
      end   
 
%     Income taxes ($/y)
      if (csubtots(13,Y) > 0) 
        taxable = csubtots(13,Y) + untaxed;  
        if (taxable > 0) 
          csubtots(14,Y) = corptax*taxable;  
          untaxed = 0;
        else   
          csubtots(14,Y) = 0;
          untaxed = taxable; 
        end   
      else   
        csubtots(14,Y) = 0;   
      end 
    end  
 
%   Income after taxes ($/y)
    csubtots(15,1:plife) = csubtots(13,1:plife)-csubtots(14,1:plife); 
 
%   Annual cash income ($/y)
    csubtots(16,1:plife) = cmanu(12,:) + csubtots(15,1:plife);   
 
%   Annual cash flow ($/y)
%   Year zero is only capital investment
    csubtots(17,yzero) = csubtots(4,yzero);  
%   Years 1 through project life is sum of capital investment
%    and annual cash income
    csubtots(17,1:plife) = csubtots(16,1:plife);   
 
%   Discount factors in year Y to give present worth occuring uniformly 
%    over one-year periods after the reference point
%   Combined interest rate - uses interest rate and inflation rate
    icomb = inflat + interest + inflat*interest;  
    for Y = 1:plife  
      csubtots(18,Y) = (exp(icomb)-1)/icomb*exp(-icomb*Y);  
%     Discount factors in year Y to give present worth occuring in an 
%      instant at a point in time after the reference point
      csubtots(19,Y) = exp(-icomb*Y);   
 
%     Annual present value ($/y)
      if (Y < plife) 
%       Sum of Discounted annual cash flow and discounted salvage value 
%        & total acpital investment
        csubtots(20,Y) = csubtots(18,Y)*csubtots(17,Y)+csubtots(19,Y)*(csubtots(5,Y)-csubtots(4,Y));   
      else   
%       End of project add back discounted working capital
        csubtots(20,Y) = csubtots(18,Y)*csubtots(17,Y)+csubtots(19,Y)*(csubtots(5,Y)-csubtots(4,Y)+cinvest(14,yzero));   
       end  
     end 

%    Place working capital in final profect year  
     cinvest(14,plife) = cinvest(14,yzero); 
 
%    Total present value of annual cash flows ($)
     csubtots(20,yzero) = csubtots(20,yzero)+sum(csubtots(20,1:plife));  

%    Net present worth ($)
     present = csubtots(20,yzero)-csubtots(17,yzero);  

    if flag == 9

%   Printed tables for estimated cash flow  

%   Extend arrays for print
    if (extend > 0) 
      capacity(1,plife+extend) = 0; 
      cmanu(1,plife+extend) = 0;   
    end 

    for Y = 1:T
      if (Y == 1)
    fprintf(fid,' \n');
    fprintf(fid,'ESTIMATED CASH FLOW STATEMENT (YEAR 0 - %1.0f) \n',Y+4);
    fprintf(fid,' \n');
    fprintf(fid,'Item                                     Ratio      Year %1.0f      Year %1.0f      Year %1.0f      Year %1.0f      Year %1.0f      Year %1.0f \n',Y-1,Y,Y+1,Y+2,Y+3,Y+4);
    fprintf(fid,'                                         Factor       ($)         ($)         ($)         ($)         ($)         ($) \n');
    fprintf(fid,'------------------------------------------------------------------------------------------------------------------------ \n');
    fprintf(fid,'PERCENT OPERATING CAPACITY                                   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',100 .*capacity(Y:Y+4));
    fprintf(fid,'------------------------------------------------------------------------------------------------------------------------ \n');
    fprintf(fid,'FIXED CAPITAL INVESTMENT                         %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(3,yzero),csubtots(1,Y:Y+4));
    fprintf(fid,'WORKING CAPITAL                                  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvest(14,yzero),cinvest(14,Y:Y+4));
    fprintf(fid,'TOTAL CAPITAL INVESTMENT                         %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(4,yzero),csubtots(4,Y:Y+4));
    fprintf(fid,'TOTAL INCOME (SALES)                                         %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(11,Y:Y+4));
    fprintf(fid,'  Raw Materials                                              %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanu(1,Y:Y+4));
    fprintf(fid,'  Operating Labor, OL                                        %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanu(2,Y:Y+4));
    fprintf(fid,'  Oper. Supervision & Cleric, OS        %3.2f OL              %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(1),cmanu(3,Y:Y+4));
    fprintf(fid,'  Opererating Supplies                  %3.2f M               %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(4),cmanu(8,Y:Y+4));
    fprintf(fid,'  Maintenance Labor, ML                                      %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanu(5,Y:Y+4));
    fprintf(fid,'  Maint. Supervision & Cleric           %3.2f ML              %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(3),cmanu(7,Y:Y+4));
    fprintf(fid,'  Maintenance Supplies                  %3.2f PE              %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(2),cmanu(6,Y:Y+4));
    fprintf(fid,'  Power & Utilities                                          %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanu(10,Y:Y+4));
    fprintf(fid,'  Waste Discharge Charges                                    %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanu(15,Y:Y+4));
    fprintf(fid,'  Laboratory Charges                    %3.2f OL              %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(5),cmanu(9,Y:Y+4));
    fprintf(fid,'  Patents & Royalties                   %3.2f TPC             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(6),cmanu(11,Y:Y+4));
    fprintf(fid,' DIRECT PRODUCTION COST                                      %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(6,Y:Y+4));
    fprintf(fid,'  Local (Property) Taxes                %3.2f FCI             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(7),cmanu(13,Y:Y+4));
    fprintf(fid,'  Insurance                             %3.2f FCI             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(8),cmanu(14,Y:Y+4));
    fprintf(fid,' FIXED CHARGES                                               %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(7,Y:Y+4));
    fprintf(fid,' PLANT-OVERHEAD COST                  %3.2f OL+OS+M           %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(9),cmanu(16,Y:Y+4));
    fprintf(fid,'TOTAL MANUFACTURING COST                                     %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(8,Y:Y+4));
    fprintf(fid,'  Administrative Costs                  %3.2f TPC             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(10),cmanu(17,Y:Y+4));
    fprintf(fid,'  Distribution and Selling Costs        %3.2f TPC             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(11),cmanu(18,Y:Y+4));
    fprintf(fid,'  Research and Development              %3.2f TPC             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(12),cmanu(19,Y:Y+4));
    fprintf(fid,'  Financing                             %3.2f FCI             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(13),cmanu(20,Y:Y+4));
    fprintf(fid,'GENERAL EXPENSES                                             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(9,Y:Y+4));
    fprintf(fid,'TOTAL PRODUCT COST,TPC                                       %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(10,Y:Y+4));
    fprintf(fid,'ANNUAL OPERATING INCOME                                      %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(8,Y:Y+4));
    fprintf(fid,'DEPRECIATION                                                 %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanu(12,Y:Y+4));
    fprintf(fid,'INCOME BEFORE TAXES                                          %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(13,Y:Y+4));
    fprintf(fid,'INCOME TAXES                                                 %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(14,Y:Y+4));
    fprintf(fid,'INCOME AFTER %4.3f TAXES                                     %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',corptax,csubtots(15,Y:Y+4));
    fprintf(fid,'ANNUAL CASH INCOME                                           %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(16,Y:Y+4));
    fprintf(fid,'ANNUAL CASH FLOW                                 %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(17,yzero),csubtots(17,Y:Y+4));
    fprintf(fid,' Discount Factors Continuous Uniform %4.3f INTEREST   %4.3f       %4.3f       %4.3f       %4.3f       %4.3f       %4.3f  \n',interest,1.0,csubtots(18,Y:Y+4));
    fprintf(fid,'SALVAGE                                                      %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(5,Y:Y+4));
    fprintf(fid,' Discount Factors Finite Year End    %4.3f INTEREST   %4.3f       %4.3f       %4.3f       %4.3f       %4.3f       %4.3f  \n',interest,1.0,csubtots(19,Y:Y+4));
    fprintf(fid,'ANNUAL PRESENT VALUE                             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(17,yzero),csubtots(20,Y:Y+4));
    fprintf(fid,'------------------------------------------------------------------------------------------------------------------------ \n');
    fprintf(fid,'TOTAL PRESENT VALUE (YEARS 1 TO LIFE)            %10.0f\n',csubtots(20,yzero));
    fprintf(fid,'NET PRESENT WORTH                                %10.0f\n',present);
    fprintf(fid,'------------------------------------------------------------------------------------------------------------------------ \n');
    fprintf(fid,'M = sum of all maintenance costs\n');  
    fprintf(fid, ' \n');   

      else 
    fprintf(fid,' \n');
    fprintf(fid,'ESTIMATED CASH FLOW STATEMENT (YEAR %2.0f - %2.0f) \n',5*Y-4,5*Y);
    fprintf(fid,' \n');
    fprintf(fid,'Item                                     Ratio      Year %1.0f      Year %1.0f      Year %1.0f      Year %1.0f      Year %1.0f \n',5*Y-4,5*Y-3,5*Y-2,5*Y-1,5*Y);
    fprintf(fid,'                                         Factor       ($)         ($)         ($)         ($)         ($)  \n');
    fprintf(fid,'------------------------------------------------------------------------------------------------------------ \n');
    fprintf(fid,'PERCENT OPERATING CAPACITY                       %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',100 .*capacity(5*Y-4:5*Y));
    fprintf(fid,'------------------------------------------------------------------------------------------------------------ \n');
    fprintf(fid,'FIXED CAPITAL INVESTMENT                         %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(1,5*Y-4:5*Y));
    fprintf(fid,'WORKING CAPITAL                                  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cinvest(14,5*Y-4:5*Y));
    fprintf(fid,'TOTAL CAPITAL INVESTMENT                         %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(4,5*Y-4:5*Y));
    fprintf(fid,'TOTAL INCOME (SALES)                             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f    \n',csubtots(11,5*Y-4:5*Y));
    fprintf(fid,'  Raw Materials                                  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f    \n',cmanu(1,5*Y-4:5*Y));
    fprintf(fid,'  Operating Labor, OL                            %10.0f  %10.0f  %10.0f  %10.0f  %10.0f     \n',cmanu(2,5*Y-4:5*Y));
    fprintf(fid,'  Oper. Supervision & Cleric, OS        %3.2f OL  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(1),cmanu(3,5*Y-4:5*Y));
    fprintf(fid,'  Opererating Supplies                  %3.2f M   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(4),cmanu(8,5*Y-4:5*Y));
    fprintf(fid,'  Maintenance Labor, ML                          %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanu(5,5*Y-4:5*Y));
    fprintf(fid,'  Maint. Supervision & Cleric           %3.2f ML  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(3),cmanu(7,5*Y-4:5*Y));
    fprintf(fid,'  Maintenance Supplies                  %3.2f PE  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(2),cmanu(6,5*Y-4:5*Y));
    fprintf(fid,'  Power & Utilities                              %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanu(10,5*Y-4:5*Y));
    fprintf(fid,'  Waste Discharge Charges                        %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanu(15,5*Y-4:5*Y));  
    fprintf(fid,'  Laboratory Charges                    %3.2f OL  %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(5),cmanu(9,5*Y-4:5*Y));
    fprintf(fid,'  Patents & Royalties                   %3.2f TPC %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(6),cmanu(11,5*Y-4:5*Y));
    fprintf(fid,' DIRECT PRODUCTION COST                          %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(6,5*Y-4:5*Y));
    fprintf(fid,'  Local (Property) Taxes                %3.2f FCI %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(7),cmanu(13,5*Y-4:5*Y));
    fprintf(fid,'  Insurance                             %3.2f FCI %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(8),cmanu(14,5*Y-4:5*Y));
    fprintf(fid,' FIXED CHARGES                                   %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(7,5*Y-4:5*Y));
    fprintf(fid,' PLANT-OVERHEAD COST                 %3.2f OL+OS+M%10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(9),cmanu(16,5*Y-4:5*Y));
    fprintf(fid,'TOTAL MANUFACTURING COST                         %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(8,5*Y-4:5*Y));
    fprintf(fid,'  Administrative Costs                  %3.2f TPC %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(10),cmanu(17,5*Y-4:5*Y));
    fprintf(fid,'  Distribution and Selling Costs        %3.2f TPC %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(11),cmanu(18,5*Y-4:5*Y));
    fprintf(fid,'  Research and Development              %3.2f TPC %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(12),cmanu(19,5*Y-4:5*Y));
    fprintf(fid,'  Financing                             %3.2f FCI %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanurat(13),cmanu(20,5*Y-4:5*Y));
    fprintf(fid,'GENERAL EXPENSES                                 %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(9,5*Y-4:5*Y));
    fprintf(fid,'TOTAL PRODUCT COST,TPC                           %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(10,5*Y-4:5*Y));
    fprintf(fid,'ANNUAL OPERATING INCOME                          %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(8,5*Y-4:5*Y));
    fprintf(fid,'DEPRECIATION                                     %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',cmanu(12,5*Y-4:5*Y));
    fprintf(fid,'INCOME BEFORE TAXES                              %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(13,5*Y-4:5*Y));
    fprintf(fid,'INCOME TAXES                                     %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(14,5*Y-4:5*Y));
    fprintf(fid,'INCOME AFTER %4.3f TAXES                         %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',corptax,csubtots(15,5*Y-4:5*Y));
    fprintf(fid,'ANNUAL CASH INCOME                               %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(16,5*Y-4:5*Y));
    fprintf(fid,'ANNUAL CASH FLOW                                 %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(17,5*Y-4:5*Y));
    fprintf(fid,' Discount Factors Continuous Uniform %4.3f INTEREST   %4.3f       %4.3f       %4.3f       %4.3f       %4.3f \n',interest,csubtots(18,5*Y-4:5*Y));
    fprintf(fid,'SALVAGE                                          %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(5,5*Y-4:5*Y));
    fprintf(fid,' Discount Factors Finite Year End    %4.3f INTEREST   %4.3f       %4.3f       %4.3f       %4.3f       %4.3f \n',interest,csubtots(19,5*Y-4:5*Y));
    fprintf(fid,'ANNUAL PRESENT VALUE                             %10.0f  %10.0f  %10.0f  %10.0f  %10.0f  \n',csubtots(20,5*Y-4:5*Y));
    fprintf(fid,'------------------------------------------------------------------------------------------------------------ \n');

 
      end 
    end % flag = 9 

    end % LEVEL 4

  end % flag = 3 or 9

  sys = [];

else
   % Otherwise, no need to return anything since this is continuous system
   sys = [];

end 
 
