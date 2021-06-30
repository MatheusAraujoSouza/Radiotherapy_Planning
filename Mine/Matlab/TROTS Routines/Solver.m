    function Solver(Path, Analysis, Method, p, dataID_TGoal, dataID_CUB, dataID_RUB)
%
% This routine inlcude the CPLEX solver and the analysis window created by
% Allen Holder.
%
% OBSERVATION: "Regular" stands for "Good tissue"
%
[TGoal, RUB, CUB, ATumor, ACritical, ARegular] = Load_data(Path, dataID_TGoal, dataID_CUB, dataID_RUB);



%
% tol - the amount that we are going to increase the tumor lower
%       bound, TLB, and then compute \omega to guarantee the
%       tumor uniformity
% 
tol = 10 ^ -4;


%
% Grab some sizes.
%
% tr = tumor row, tc = tumor colum
% gr = critical row, gc = critical colum
% rr = regular row, rc = regular colum

[tr,tc] = size(ATumor);
[gr,gc] = size(ACritical);
[rr,rc] = size(ARegular);

% Some matrices are not in double precision. We need to put it into double
% precision to use sparse.

ATumor    = sparse(double(ATumor));
ACritical = sparse(double(ACritical));
ARegular  = sparse(double(ARegular));

%
% Extract the upper and lower bounds from the prescription matrix
%

TUB = TGoal * (1 + p);
TLB = (TGoal * (1 - p)) + tol;


%
% Calculate the weight omega that guarantees the tumor uniformity.
%
% The bound presented in the papaer leads to numerical instability
% with the interior point method always thinking that it is infeasible.
% This value seems to work well.
%

kappa = (norm(TLB, inf));
omega = kappa / tol;

%
% Build the data for the solver
%
%
%                Ux = u
%      
% ( A_T    0    0    0 ) (  x   )    ( TUB )
% (-A_T   -L    0    0 ) (\alpha)    (-TUB )
% ( A_C    0  -U_c   0 ) (\beta ) <= ( CUB )
% ( A_R    0    0  -U_r) (\gamma)    ( RUB )
%
%
%     Lower bounds             Upper bounds
%
% (  x   )    ( inf )        (  x   )    (  0  )
% (\alpha)    ( TLB )        (\alpha)    (  0  )
% (\beta ) <= ( inf )        (\beta ) >= (-CUB )
% (\gamma)    ( inf))        (\gamma)    (  0  )
% 
%

if Analysis == "Average"

  U  = sparse([ ATumor,    zeros(tr, tr), zeros(tr, gr), zeros(tr, rr);
               -ATumor,      -eye(tr),    zeros(tr, gr), zeros(tr, rr);
                ACritical, zeros(gr, tr),      -eye(gr), zeros(gr, rr);
                ARegular,  zeros(rr, tr), zeros(rr, gr),    -eye(rr)]);
  u  = sparse([TUB; -TLB; CUB; RUB]);
  lb = sparse([zeros(tc, 1); zeros(tr, 1); -CUB; zeros(rr,1)]);
  ub = sparse([inf * ones(tc, 1); TLB; inf * ones(gr, 1); inf * ones(rr, 1)]);
  c  = sparse([zeros(tc,1); omega * ones(tr, 1); ones(gr, 1); ones(rr, 1)]);

elseif Analysis == "Absolute"

  U  = sparse([ ATumor,    zeros(tr, 1), zeros(tr, 1), zeros(tr, 1);
               -ATumor,     -ones(tr,1), zeros(tr, 1), zeros(tr, 1);
                ACritical, zeros(gr, 1),  -ones(gr,1), zeros(gr, 1);
                ARegular,  zeros(rr, 1), zeros(rr, 1), -ones(rr,1)]);
  u  = sparse([TUB; -TLB; CUB; RUB]);
  lb = sparse([zeros(tc, 1); 0; -min(CUB); 0]);
  ub = sparse([inf * ones(tc, 1); min(TLB); inf; inf]);
  c  = sparse([zeros(tc,1); omega; 1; 1]);

else

  disp('Analysis Type is Not Selected')
  return

end


% Call CPLEX solver and set some parameters:
% Time limit of 900s
% Optimality tolerance = 1e^-04

StartTime = clock;

if Method == "Simplex"
    opt = cplexoptimset('cplex');
    opt.timelimit = 900;
    opt.lpmethod = 1;
    opt.simplex.tolerances.optimality = 1e-04;
    [x,~,~,output] = cplexlp(c, U, u, [], [],lb,ub,[],opt);
elseif Method == "Interior point"
    opt = cplexoptimset('cplex');
    opt.timelimit = 900;
    opt.lpmethod = 4;
    opt.simplex.tolerances.optimality = 1e-04;
    [x,~,~,output] = cplexlp(c, U, u, [], [],lb,ub,[],opt);
end

TotalTime = etime(clock,StartTime);

%
% Partition the solution up into the different types of
% decision variables.
%
           
if Analysis == "Average"

  alpha = x(tc+1:tc+tr);
  beta  = x(tc+tr+1:tc+tr+gr);
  gamma = x(tc+tr+gr+1:tc+tr+gr+rr);
  x     = x(1:tc);

elseif Analysis == "Absolute"

  alpha = x(tc+1);
  beta  = x(tc+2);
  gamma = x(tc+3);
  x     = x(1:tc);

end

A = [ATumor;ACritical;ARegular];
RadLevel = A*x;

% 
% We now supply some usefull information from the solution.
% Start by opening the window for the analysis to appear.
% 

AnalysisWindow

%
% Display the type of Analysis being conducted.
%

if Analysis == "Average"
  set(findobj('Tag','AnalType'),'String','An Average Analysis was used.')
else
  set(findobj('Tag','AnalType'),'String','An Absolute Analysis was used.')
end

%
% Display Information about whether the tumor is recieveing an
% appropriate dose.
%

if Analysis == "Average"

  if norm(alpha,1) > tol
	set(findobj('Tag','AnalTumor'),'String', ...
	   sprintf('On average, the Tumors Prescription is not attainable.'))
  else
	AverageTumorDose = norm(ATumor*x,1)/tr;
	set(findobj('Tag','AnalTumor'),'String', ...
	   sprintf('The average tumor dose is within the prescribed uniformity setting. The average tumor dose is %2.2f.',AverageTumorDose))
  end

else

  if alpha > tol
        set(findobj('Tag','AnalTumor'),'String', ...
           sprintf('The Tumors Prescription is not attainable.'))
  else
        MinDose = min(ATumor*x);
	MaxDose = max(ATumor*x);
        set(findobj('Tag','AnalTumor'),'String', ...
           sprintf('The tumor dose is within the prescribed uniformity setting. The min. (max.) tumor dose is %2.2f (%2.2f).',MinDose,MaxDose))
  end

end

%
% Display How the Critical Structures fair under our plan.
%

if Analysis == "Average"

  AverageCritDesc = mean(beta);
  if AverageCritDesc > tol
	set(findobj('Tag','AnalCrit'),'String', ...
	   sprintf('The Average Critical Structure Dose is Over its Prescribed Dose by %2.2f',AverageCritDesc))
  else
 	set(findobj('Tag','AnalCrit'),'String', ...
	   sprintf('The Average Critical Structure Dose is Under its Prescribed Dose by %2.2f',abs(AverageCritDesc)))
   end

else

  if beta > tol
	set(findobj('Tag','AnalCrit'),'String', ...
	   sprintf('The Critical Structures Must Recieve an Excessive Dose of at least %2.2f',beta))
  else
	set(findobj('Tag','AnalCrit'),'String', ...
           sprintf('The Critical Structures Recieved at least %2.2f Under Their Prescribed Dose.',abs(beta)))
  end

end

%
% State How We Did on The Regular Tissue
%

if (Analysis == "Average") && (rr ~= 0)

  AverageRegDesc = mean(gamma);

  if AverageRegDesc > tol
	set(findobj('Tag','AnalHot'),'String', ...
	   sprintf('On Average, the Restricted Regular Tissue is Above its Prescribed Dose by %2.2f.',AverageRegDesc))
  elseif AverageRegDesc < tol
        set(findobj('Tag','AnalHot'),'String', ...
           sprintf('On Average, the Restricted Regular Tissue is Below its Prescribed Dose.'))

  end

elseif rr ~= 0

  if gamma > tol
	set(findobj('Tag','AnalHot'),'String', ...
	   sprintf('The Restricted Regular Tissue is Guaranteed to be Above its Prescribed Dose by %2.2f.',gamma))
  elseif gamma < tol
        set(findobj('Tag','AnalHot'),'String', ...
           sprintf('The Restricted Regular Tissue is Guaranteed to be Below its Prescribed Dose.'))

  end

else

  set(findobj('Tag','AnalHot'),'String', ...
	sprintf('There were no restrictions on the Regular Tissue.'))

end

%
% State the Solver Type
%

if Method == "Interior Point"
   set(findobj('Tag','AnalSolver'),'String', ...
	sprintf('An Interior Point Method Was Used.'))
else
   set(findobj('Tag','AnalSolver'),'String', ...
	sprintf('A Simplex Method Was Used.'))
end

%
% State the Time required by the solver.
%

set(findobj('Tag','AnalTime'),'String', ...
	sprintf('The solution time was %4.2f seconds.',TotalTime))

%
% State any notes.
%


if max(RadLevel) <= 1.1*max(TUB)
  set(findobj('Tag','AnalNotes'),'String', ...
	sprintf('There are no hot spots.'))

else

  set(findobj('Tag','AnalNotes'),'String', ...
  sprintf('There are hot spots. You may want to make further restrictions.'))

end
