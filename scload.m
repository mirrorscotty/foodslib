     function[Fname] = scload(SCfile,strmtype)  

%  SCLOAD    Loads user-specified stream costs data file for economic
%            analysis  
%
%            Synopsis: MATLAB function which sets the filename for the 
%            FOODS-LIB model design results.
% 
%            Call: scload(SCfile)  
%
%            Returned Argument:
%             filesc - user-specified stream cost data file name 

%  *************************************************************************  

%   Open specified stream cost/value data file 
    fid = fopen(SCfile);
%   Read in the data in specified format into a vector 
    SCvector = fscanf(fid,'%g %g %g %g %s %g');   
%   Close data file 
    fclose(fid);

%   Find index for end of each stream cost/value record
    J = 0;   
    for I = 1:length(SCvector)
      if (SCvector(I) < 0)
        J = J + 1;  
        endstrm(J) = I;  
      end
    end


    for I = 1:length(endstrm)    % Number of streams in econ file
      if I == 1
        SCtype(I) = SCvector(1)    % Stream type  
        SC(I) = SCvector(2)        % Stream mass flow rate (kg/h) 
        SCcost(I) = SCvector(3)    % Stream cost/value ($/kg)  
        SCinfl(I) = SCvector(4)    % Stream inflation (decimal %)  
        SCnamel = endstrm(I)-1-4   % Stream name length and name  
        SCname(I,1:SCnamel) = setstr(SCvector(5:(endstrm(I)-1)))'
      else
        SCtype(I) = SCvector(endstrm(I-1)+1)
        SC(I) = SCvector(endstrm(I-1)+2)
        SCcost(I) = SCvector(endstrm(I-1)+3)
        SCinfl(I) = SCvector(endstrm(I-1)+4)
        SCnamel = endstrm(I) -1 - endstrm(I-1) - 4
        SCname(I,1:SCnamel) = setstr(SCvector((endstrm(I-1)+5):(endstrm(I)-1)))'
      end
    end
 
    J = 0;  

    for I = 1:length(endstrm)
      if (SCtype(I) == strmtype) 
        J = J + 1; 
        F(J,1) = SC(I);
        Fcost(J) = SCcost(I);
        Finfl(J) = SCinfl(I); 
        Fname(J,:) = SCname(I,:); 
      end 
    end 

    if J == 0
      F = [];
      Fcost = [];
      Finfl = [];
      Fname = [];
    end 

%   F = F(:,1:size(Fname,2))
 
%   FC = [F; Fname]

