function Table = Solver(Difficulty, Table, Counter, Case, Path, Analysis, Method, p, dataID_TGoal, dataID_CUB, dataID_RUB)
%
% This routine inlcudes the:
%   - CPLEX solver of the LP
%   - Analysis of the solution
%   - Table to store the results
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

if Analysis == "Average"
    
    
  U  = sparse([ ATumor,    zeros(tr, tr), zeros(tr, gr), zeros(tr, rr);
               -ATumor,         -eye(tr), zeros(tr, gr), zeros(tr, rr);
                ACritical, zeros(gr, tr),      -eye(gr), zeros(gr, rr);
                ARegular,  zeros(rr, tr), zeros(rr, gr),    -eye(rr)]);
  u  = sparse([TUB; -TLB; CUB; RUB]);
  lb = sparse([zeros(tc, 1); zeros(tr, 1); -CUB; zeros(rr,1)]);
  ub = sparse([inf * ones(tc, 1); TLB; inf * ones(gr, 1); inf * ones(rr, 1)]);
  c  = sparse([zeros(tc,1); (1/tc) * omega * ones(tr, 1); (1/gc) * ones(gr, 1); (1/rc) * ones(rr, 1)]);

elseif Analysis == "Absolute"

  U  = sparse([ ATumor,    zeros(tr, 1), zeros(tr, 1), zeros(tr, 1);
               -ATumor,     -ones(tr,1), zeros(tr, 1), zeros(tr, 1);
                ACritical, zeros(gr, 1),  -ones(gr,1), zeros(gr, 1);
                ARegular,  zeros(rr, 1), zeros(rr, 1), -ones(rr,1)]);
  u  = sparse([TUB; -TLB; CUB; RUB]);
  lb = sparse([zeros(tc, 1); 0; -min(CUB); 0]);
  ub = sparse([inf * ones(tc, 1); min(TLB); inf; inf]);
  c  = sparse([zeros(tc,1); omega; 1; 1]);
  
end

%
% Getting the size of the problem
%

[Number_of_Constraints, Number_of_Variables] = size(U);

% Call CPLEX solver and set some parameters:
% Time limit of 18000s - 5 hours

if Method == "Automatic"
    
    opt = cplexoptimset('cplex');
    opt.timelimit = 18000;
    [x,~,~,output] = cplexlp(c, U, u, [], [],lb,ub,[],opt);
elseif Method == "Interior point"
    
    opt = cplexoptimset('cplex');
    opt.timelimit = 18000;
    opt.lpmethod = 4;
    [x,~,~,output] = cplexlp(c, U, u, [], [],lb,ub,[],opt);
end

save("x_" + num2str(Counter),'x')

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

% 
% We now supply some usefull information from the solution.
%


A = [ATumor;ACritical;ARegular];
RadLevel = A*x;

%
% Information about whether the tumor is recieveing an
% appropriate dose.
%

if Analysis == "Average"                                                            %Evite repetição de código, da linha 152 a 187 está fazendo a mesma coisa, crie 3 funções uma para testar o que é e outras duas para cada um dos testes com suas respostas
    if norm(alpha,1) > tol
        Tumor_dose = "On average, the tumors prescription is not attainable.";      
    else
        AverageTumorDose = norm(ATumor*x,1)/tr;
        Tumor_dose = "The average tumor dose is within the prescribed uniformity setting. The average tumor dose is " + num2str(AverageTumorDose);
    end
else
    if alpha > tol
        Tumor_dose = "The tumors prescription is not attainable.";
    else
        MinDose = min(ATumor*x);
        MaxDose = max(ATumor*x);
        Tumor_dose = "The tumor dose is within the prescribed uniformity setting. The min and max tumor dose are " + num2str(MinDose) + " and " + num2str(MaxDose);
    end
end

%
% How the Critical Structures fair under our plan.
%

if Analysis == "Average"                                                                                     %Dado que você já fez o teste do que é você só precisa criar mais duas funções, não tenha medo de fazer funções pequenas
    AverageCritDesc = mean(beta);                                                                            % se você fizer duas funções cada uma fica com 3 linhas de código
    if AverageCritDesc > tol
        Critical_dose = "The average citical sructure dose is over its prescribed dose by " + num2str(AverageCritDesc);
    else
        Critical_dose = "The average critical structure dose is under its prescribed dose by " + num2str(AverageCritDesc);
    end

else
    if beta > tol
        Critical_dose = "The critical structures must recieve an excessive dose of at least " + num2str(beta);
    else
        Critical_dose = "The critical structures recieved at least "  + num2str(beta)+ " under their prescribed dose.";
    end
end


%
% State How We Did on The Regular Tissue
%
%do 193 ao 209 jogue dentro de uma função, coloque um nome que explica todos esses tests depois disso reduza esses blocos internos, não tenha medo de usar returns esqueça isso da programação estruturada que te fala que só pode ter um return.
if (Analysis == "Average") && (rr ~= 0)  %siga um padrão nos nomes das suas variaveis, nunca comece com letra maiuscula com nome de variavel o correto: averageRegDesc e nunca use "_" o correto é regularDose                                                        
    AverageRegDesc = mean(gamma);
    if AverageRegDesc > tol %assumindo que nunca vai ser igual, você só tem dois resultados não precisa usar if e else ou é um ou é outro, você tem certeza que nunca vai ser igual?, pode usar um operadores ternario e reduzir isso tudo para 3 linhas 
	    Regular_dose = "On average, the restricted regular tissue is above its prescribed dose by " + num2str(AverageRegDesc);
    elseif AverageRegDesc < tol
        Regular_dose = "On average, the restricted regular tissue is below its prescribed dose.";
    end
elseif rr ~= 0
    if gamma > tol     %se nunca vai ser igual faça igual construa uma função e use um operador ternario para reduzir tudo isso para 3 linhas dentro de uma função.
        Regular_dose = "The restricted regular tissue is guaranteed to be above its prescribed dose by " + num2str(gamma);
    elseif gamma < tol
        Regular_dose = "The restricted regular tissue is guaranteed to be below its prescribed dose.";
    end
else
    Regular_dose = "There were no restrictions on the regular tissue.";
end

%
% About hot spots
%

if max(RadLevel) <= 1.1*max(TUB)       %mude para radLevel
    Hot_spots = "There are no hot spots.";    % não use "_", mude para hotSpots
else
    Hot_spots = "There are hot spots. You may want to make further restrictions.";
end
%Mude para assim se não quiser usar uma função para isso, mas com uma função fica com só 2 linhas também.
%hotSplots = if (max(RadLevel) <= 1.1*max(TUB)); casetrue("There are no hot spots."); else; casefalse("There are hot spots. You may want to make further restrictions."); end



%
% Add rows to the table
%

Table(Counter,:) = {Counter, Case, Difficulty, Analysis, Method, Number_of_Constraints, Number_of_Variables, output.cplexstatus, output.cplexstatusstring, output.message, output.iterations, output.time, Tumor_dose, Critical_dose, Regular_dose, Hot_spots}; 

 end
%Comece o nome de todas as variaveis com minusculo, boas praticas!!!
