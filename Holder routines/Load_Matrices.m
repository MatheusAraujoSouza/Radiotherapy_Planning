%
% DesignPlan version 0.1 - 4/5/00
%
%
% This function recieves the prescription, (TLB, TUB, CUB, RUB),
% from the modeling GUI.  The following model is solved:
%
%  min \omega l^T \alpha + u_C^T \beta + u_G^T \gamma
%
%  such that
%
%     TLB - L \alpha <=  A_T x   <= TUB
%                        A_C x   <= CUB + U_C \beta  
%                        A_G x   <= RUB + U_G \gamma (RUB is GUB in paper)
%                  0 <= L \alpha <= TLB
%               -RUB <= U_C \beta
%                  0 <= U_G \gamma
%                  0 <= x
%
%  Complete details of this model are found in "Designing Radiotherapy
%  Plans with Elastic Constraints and Interior Point Methods" by Allen
%  Holder, Trinity University Mathematics, Technical Report #59, 2000.
%
%  Allen Holder - aholder@trinity.edu
%                 http://www.trinity.edu/aholder
%

function Load_Matrices

% Define some global variables that are used in GUI.
% The global variables are used mostly for speed.
%
%
% prescription - an N x N x 2 matrix;
%       (:,:,1) contains the prescribed limits,
%       (:,:,2) is 0,1,2
%               0 - no prescribed limit
%               1 - indicates a tumor pixel
%               2 - indicates a regular tissue
%
% DisplayPrescription - same as the first level of prescription
% 	but that all regular tissue levels are set at 5 for viewing
%	purposes.
%
% UniLevel - Uniformity Level for the for the tumor
%

%global DisplayPrescription prescription UniLevel
load pre

%
% Load in the dosematrix A.  This is a link to one of
% many files, each containing a different geometry for the
% problem.  Change the link to either add or remove complexity.
%
%Tirei os "flops" porque não funcionavam mais

%dosematrix_1_64_32: o primeiro digito significa a discretização dos
%ângulos, o segundo é o grid da imagem, e o terceiro é o número de pencils.

%if get(findobj('Tag','15_10_MatrixBtn'),'Value') == 1
       load dosematrix_1_64_32 %(alteração Vinicius Jameli)
%elseif get(findobj('Tag','5_32_MatrixBtn'),'Value') == 1
%        load /export/home0/aurelio/rad/matrices/dosematrix_5_64_32.mat
%elseif get(findobj('Tag','1_32_MatrixBtn'),'Value') == 1
%        load /export/home0/aurelio/rad/matrices/dosematrix_1_64_32.mat
%else
%        disp('No matrix to load.')
%        return;
%end

%
% NumARows - number of rows in A, N x N
% NumACols - number of columns in A, \Theta * \eta
%

[NumARows,NumACols] = size(A);

%
% This is used later to set up a display matrix
%

n = length(prescription);

%
% p is the precentage of deviation allowed over the tumor.
%

UniLevel= 2;
p = UniLevel / 100;

%
% tol - the amount that we are going to increase the tumor lower
%       bound, TLB, and then compute \omega to guarantee the
%       tumor uniformity
%

tol = 10 ^ -4;


%
% These index vectors show where the three types of tissue are.
%

TumorRowIndex = find(prescription(:, : , 2) == 1);
CriticalRowIndex = find(prescription(:, : , 2) == 2);
%RegularRowIndex = find(prescription(:, : , 2) == 3);
RegularRowIndex= find(prescription(:,:,2) == 3 | prescription(:,:,2) == 0);

%
% Partition the dose deposition matrix A into those rows correpsonding
% to tumor pixels, Critical structure pixels, and regular tissue pixels.
%

ATumor    = sparse(A(TumorRowIndex,:));
ACritical = sparse(A(CriticalRowIndex,:));
ARegular  = sparse(A(RegularRowIndex,:));

%
% Grab some sizes.
%

[tr,tc] = size(ATumor);
[gr,gc] = size(ACritical);
[rr,rc] = size(ARegular);

%
% Remove the columns corresponding to the sub-beams that do not
% strike the tumor.
%

CheckColSum = ones(1,tr)*ATumor;
TumorColIndex = find(CheckColSum >= 0.1);
ATumor    = sparse(ATumor(:,TumorColIndex));
ACritical = sparse(ACritical(:,TumorColIndex));
ARegular  = sparse(ARegular(:,TumorColIndex));

%
% Extract the upper and lower bounds from the prescription matrix
%

TGoal = prescription(TumorRowIndex);
TUB = TGoal * (1 + p);
TLB = (TGoal * (1 - p)) + tol;
CUB = prescription(CriticalRowIndex);
RUB = prescription(RegularRowIndex);

%
% Grab the size of the reduced ATumor
%

[tr,tc] = size(ATumor);

%
% Set the Analysis Type.
%	0 - Average Analysis
%	1 - Absolute Analysis
%

%AnalysisType = get(findobj('Tag','AbsoluteBtn'),'Value');
AnalysisType = 0;

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

if AnalysisType == 0

  U  = sparse([ ATumor,    zeros(tr, tr), zeros(tr, gr), zeros(tr, rr);
               -ATumor,      -eye(tr),    zeros(tr, gr), zeros(tr, rr);
                ACritical, zeros(gr, tr),      -eye(gr), zeros(gr, rr);
                ARegular,  zeros(rr, tr), zeros(rr, gr),    -eye(rr)]);
  b  = sparse([TUB; -TLB; CUB; RUB]);
  lbounds = sparse([zeros(tc, 1); zeros(tr, 1); -CUB; zeros(rr,1)]);
  ubounds = sparse([inf*ones(tc, 1); TLB; inf*ones(gr, 1); inf*ones(rr, 1)]);
  c  = sparse([zeros(tc,1); omega * ones(tr, 1); ones(gr, 1); ones(rr, 1)]);

elseif AnalysisType == 1

  U  = sparse([ ATumor,    zeros(tr, 1), zeros(tr, 1), zeros(tr, 1);
               -ATumor,     -ones(tr,1), zeros(tr, 1), zeros(tr, 1);
                ACritical, zeros(gr, 1),  -ones(gr,1), zeros(gr, 1);
                ARegular,  zeros(rr, 1), zeros(rr, 1), -ones(rr,1)]);
  b  = sparse([TUB; -TLB; CUB; RUB]);
  lbounds = sparse([zeros(tc, 1); 0; -min(CUB); 0]);
  ubounds = sparse([inf * ones(tc, 1); min(TLB); inf; inf]);
  c  = sparse([zeros(tc,1); omega; 1; 1]);

else

  disp('Analysis Type is Not Selected')
  return

end;
  A= [U speye(length(b))];
  lbounds = [lbounds; zeros(length(b),1)];
  ubounds = [ubounds; inf*ones(length(b),1)];
  c = [c; zeros(length(b),1)];
BIG= 1e10; %1e32
NAME='dosematrix';

save radio A b c lbounds ubounds BIG NAME
return

[m,n] = size(A); y0 = zeros(m,1);
pmin = max(bnrm/100, 100);
dmin = cnrm*.425; dmin = max(dmin, pmin/40);
pmin = min(pmin, 1.e+4); dmin = min(dmin, 1.e+3);

rho = min(100,bnrm);
x0 = A'*((A'*A)\(b-rho*A*e)) + rho*e;
pmin = max(pmin, -min(x0)); x0 = max(pmin,x0);
z0 = full((c+dmin).*(c > 0) + dmin*(-dmin < c & c <= 0) - c.*(c <= -dmin));
save init x0 y0 z0
