% -- Minha versão com comentários

function meu_gera_resolve_plota_cplex
% Define some global variables that are used in GUI.
% The global variables are used mostly for speed.
%
%
% prescription - an N x N x 2 matrix;
%       (:,:,1) contains the prescribed limits,
%       (:,:,2) is 1,2,3
%               1 - no prescribed limit
%               2 - indicates a tumor pixel
%               3 - indicates a regular tissue
%
% -- Ou seja, prescription (n,n,1), pra cada par de pixels, contém uma
% upperbound, e (n,n,2) para cada part de píxes, ele pode não ter nenhuma
% prescrição, ter uma prescrição do tumor, ou ser um tecido normal.
%
% DisplayPrescription - same as the first level of prescription
% 	but that all regular tissue levels are set at 5 for viewing
%	purposes.
%
% UniLevel - Uniformity Level for the for the tumor
%

%global DisplayPrescription prescription UniLevel
load ('C:\Principal\IC\Estruturas_e_rotinas\Holder_cecilia\Dados Holder\pre.mat');

%
% Load in the dosematrix A.  This is a link to one of
% many files, each containing a different geometry for the
% problem.  Change the link to either add or remove complexity.
%

%if get(findobj('Tag','15_10_MatrixBtn'),'Value') == 1
        load ('C:\Principal\IC\Estruturas_e_rotinas\Holder_cecilia\Dados Holder\dosematrix_1_64_32.mat'); %(alteração Vinicius Jameli)
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
% NumACols - number of columns in A, \Theta * \eta -- Número de ângulos
% vezes o número de discretização dos beamlets.
%

[NumARows,NumACols] = size(A);

%
% This is used later to set up a display matrix -- Lenght é a maior
% dimensão, que é 64 nesse caso.
%

n = length(prescription);

%
% p is the precentage of deviation allowed over the tumor.
% -- é aquele TUB = (1 + tol)TG, só que esse "tol" é o "p". O próximo
%"tol" é pra outra coisa

UniLevel= 2;
p = UniLevel / 100;

%
% tol - the amount that we are going to increase the tumor lower
%       bound, TLB, and then compute \omega to guarantee the
%       tumor uniformity
% -- Detalhe triloco da teoria


tol = 10 ^ -4;


%
% These index vectors show where the three types of tissue are.
%

% -- O "Find" percorre de cima pra baixo. Primeiro percorre toda a primeira
% coluna, depois passa pra segunda. Fazer o find(prescription(:,:,2) vai
% percorrer a imagem 64x64 de cima pra baixo e da esquerda pra direita
% listando os píxels. No exemplo do paper, ele percorre da esquerda pra
% direita.

TumorRowIndex = find(prescription(:, : , 2) == 1);
CriticalRowIndex = find(prescription(:, : , 2) == 2);
RegularRowIndex = find(prescription(:, : , 2) == 3);

%
% Partition the dose deposition matrix A into those rows correpsonding
% to tumor pixels, Critical structure pixels, and regular tissue pixels.
%

%-- A matriz A, que já foi criada também desse jeito de percorrer de cima
%pra baixo, tem as linhas sendo os píxels e as colunas sendo a dose do
%x-ésimo beamlet

ATumor    = sparse(A(TumorRowIndex,:));
ACritical = sparse(A(CriticalRowIndex,:));
ARegular  = sparse(A(RegularRowIndex,:));

%
% Grab some sizes.
%

% -- 
%tr = tumor row, tc = tumor colum
%gr = critical row, gc = critical colum
%rr = regular row, rc = regular colum

[tr,tc] = size(ATumor);
[gr,gc] = size(ACritical);
[rr,rc] = size(ARegular);

%
% Remove the columns corresponding to the sub-beams that do not
% strike the tumor.
%

%-- Isso aqui vai criar um vetor 1 x (numero de colunas da matriz Atumor)
%sendo cada componente do vetor a soma da coluna de Atumor.
CheckColSum = ones(1,tr)*ATumor;

%--
%Agora, veja que, se a soma de uma coluna é zero, signifca que aquele
%beamlet não bateu em nenhum píxel do tumor. Então podemos tirar esta
%coluna. Então o comando find vai encontrar as  linhas em que são não nulas
%e guardar o índice delas.
TumorColIndex = find(CheckColSum >= 0.1);

%Com os índices, agora ele vai redefinir as matrizes Atumor,ACritial e
%ARegular só para os beamlets que intersectam o tumor.
ATumor    = sparse(ATumor(:,TumorColIndex));
ACritical = sparse(ACritical(:,TumorColIndex));
ARegular  = sparse(ARegular(:,TumorColIndex));

%
% Extract the upper and lower bounds from the prescription matrix
%

%-- Veja que aqui é "Row", não é o mesmo "TumorColIndex" de antes
TGoal = prescription(TumorRowIndex);
TUB = TGoal * (1 + p);
TLB = (TGoal * (1 - p)) + tol;
CUB = prescription(CriticalRowIndex);
RUB = prescription(RegularRowIndex);

%
% Grab the size of the reduced ATumor
%

% -- 
%tr = tumor row, tc = tumor colum
[tr,tc] = size(ATumor);

%
% Set the Analysis Type.
%	0 - Average Analysis
%	1 - Absolute Analysis
%

%AnalysisType = get(findobj('Tag','AbsoluteBtn'),'Value');
AnalysisType = 1    ;

%
% Calculate the weight omega that guarantees the tumor uniformity.
%
% The bound presented in the papaer leads to numerical instability
% with the interior point method always thinking that it is infeasible.
% This value seems to work well.
%

%--Aquelas doidera de convergência pra qualquer omega
kappa = (norm(TLB, inf));
omega = kappa / tol;

%
% Build the data for the solver
%

% ( A_T    0    0    0 ) (  x   )    ( TUB )
% (-A_T   -L    0    0 ) (\alpha)    (-TUB )
% ( A_C    0  -U_c   0 ) (\beta ) == ( CUB )
% ( A_R    0    0  -U_r) (\gamma)    ( RUB )

% (  x   )    ( inf )
% (\alpha)    ( TLB )
% (\beta ) <= ( inf )
% (\gamma)    ( inf))
% 
% (  x   )    (  0  )
% (\alpha)    (  0  )
% (\beta ) >= (-CUB )
% (\gamma)    (  0  )


if AnalysisType == 0

  U  = sparse([ ATumor,    zeros(tr, tr), zeros(tr, gr), zeros(tr, rr);
               -ATumor,      -eye(tr),    zeros(tr, gr), zeros(tr, rr);
                ACritical, zeros(gr, tr),      -eye(gr), zeros(gr, rr);
                ARegular,  zeros(rr, tr), zeros(rr, gr),    -eye(rr)]);
  u  = sparse([TUB; -TLB; CUB; RUB]);
  lb = sparse([zeros(tc, 1); zeros(tr, 1); -CUB; zeros(rr,1)]);
  ub = sparse([inf * ones(tc, 1); TLB; inf * ones(gr, 1); inf * ones(rr, 1)]);
  c  = sparse([zeros(tc,1); omega * ones(tr, 1); ones(gr, 1); ones(rr, 1)]);

elseif AnalysisType == 1

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

%
% Set the Solver Type
%	0 - Interior Point Method
%	1 - Simplex Method
%

%if get(findobj('Tag','SimplexBtn'),'Value') == 1
%  SolverType = 1;
%elseif get(findobj('Tag','InteriorBtn'),'Value') == 1
  SolverType = 0;
%else
%  disp('No Solver Choosen')
%  return
%end

%
% Grab the the clock
%

StartTime = clock;


%
% Off to the Solver
%

% if SolverType == 1
%   [x, fval] = linprog(c, U, u, [], [], lb, ub, [], ...
% 	optimset('LargeScale', 'off','MaxIter',10^10) );
% elseif SolverType == 0
%   [x, fval,flag,output] = linprog(c, U, u, [], [], lb, ub, [], ...
% 	optimset('LargeScale', 'on','MaxIter',10^10,'TolFun',1e-4) );
% end


if SolverType == 1
    opt = cplexoptimset('cplex');
    opt.timelimit = 300;
    opt.lpmethod = 1;
    [x,~,~,output] = cplexlp(c, U, u, [], [],lb,ub,[],opt);
elseif SolverType == 0
    opt = cplexoptimset('cplex');
    opt.timelimit = 300;
    opt.lpmethod = 4;
    [x,~,~,output] = cplexlp(c, U, u, [], [],lb,ub,[],opt);
end

%
% Calculuate the time required by the solver.
%


output

TotalTime = etime(clock,StartTime);

%
% Partition the solution up into the different types of
% decision variables.
%
           
if AnalysisType == 0

  alpha = x(tc+1:tc+tr);
  beta  = x(tc+tr+1:tc+tr+gr);
  gamma = x(tc+tr+gr+1:tc+tr+gr+rr);
  x     = x(1:tc);

elseif AnalysisType == 1

  alpha = x(tc+1);
  beta  = x(tc+2);
  gamma = x(tc+3);
  x     = x(1:tc);

end;

%
% Embedd the smaller plan, x, into a larger one and
% calculate the radiation levels over the image.
%

ThePlan = zeros(NumACols,1);
%Daquelas 11520 beamlets (colunas) iniciais, só tem solução diferente de
%zero os que interseccionam o tumor.
ThePlan(TumorColIndex) = x;
%Assim, finalmente, RadLevel é a quantidade de dose em cada um dos 4096
%pixels.
RadLevel = A*ThePlan;
DisplayMatrix = zeros(n,n);

%
% This lines up the Radiation Levels into the corresponding
% matrix for displaying the image.
%

%Só percorre o vetor RadLevel 4096x1 e distribui na displaymatrix 64x64 de
%cima pra baixo, da esquerda pra direita.
for j=1:n
   for i=1:n
	DisplayMatrix(i,j)=RadLevel(64*(j-1)+i);
   end
end

%
% Pop-up some nifty pictures. ViewPrescription is used so that the
% regular restricted regions do not show up on the contour.
%


%--Começa a montar a figura
figure;
%-- Aqui, "ViewPrescription" se torna a matriz 64x64 das prescrições de
%dose.
ViewPrescription = prescription(:,:,1);
%-- Tiramos o tecido normal da jogada
ViewPrescription(RegularRowIndex) = 0;
%--Fazemos o "contour" disso.
contour(ViewPrescription)
hold on;
contour(DisplayMatrix)
hold off;

figure;
mesh(DisplayMatrix)

%
% We now supply some usefull information from the solution.
%
% Start by opening the window for the analysis to appear
%

%-- Análises da solução baseada nuns negócio do paper. Vamo ve o que vai
%mudar depois.

AnalysisWindow

%
% Display the type of Analysis being conducted.
%

if AnalysisType == 0
  set(findobj('Tag','AnalType'),'String','An Average Analysis was used.')
else
  set(findobj('Tag','AnalType'),'String','An Absolute Analysis was used.')
end

%
% Display Information about whether the tumor is recieveing an
% appropriate dose.
%

if AnalysisType == 0

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

if AnalysisType == 0

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

if (AnalysisType == 0) & (rr ~= 0)

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

if SolverType == 0
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
