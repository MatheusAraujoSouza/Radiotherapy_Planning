using SparseArrays
using LinearAlgebra
using JuMP
using CPLEX
#Uma função não deve fazer mais de uma coisa, Solver não traz muita informação do que realmente ela faz, não tenha medo de colocar nomes grandes SolverForEinstenEquations(...)
function Solver(Path, Analysis, Method, p, dataID_TGoal, dataID_CUB, dataID_RUB)  #Você poderia criar uma struct ou uma classe e jogar todos esses parametros nela, assim sua função iria receber apenas 
#                                                                                  um unico parametro e quando for usar dentro da função, você faz assim dataParameters.Path etc
# This routine inlcude the CPLEX solver & the analysis window created by           imagina se você precisa aumentar a quantidade desses parametros o tamanho que isso vai chegar.
# Allen Holder.
#
# OBSERVATION: "Regular" stands for "Good tissue"
#

(TGoal, RUB, CUB, ATumor, ACritical, ARegular) = Load_data(Path, dataID_TGoal, dataID_CUB, dataID_RUB)



#
# tol - the amount that we are going to increase the tumor lower()
#       bound; TLB; & then compute \omega to guarantee the
#       tumor uniformity
#
tol = 10 ^ -4 #certo

#Mude o nome das variáveis para tumorRow, tumorColum, critalRow etc, são boas praticas seguidas por todos os programadores, ninguém quer ficar voltando para ler o que era no seu comentario
#
# Grab some sizes.
#
# tr = tumor row; tc = tumor colum
# gr = critical row; gc = critical colum
# rr = regular row; rc = regular colum

(tr,tc) = size(ATumor) #certo
(gr,gc) = size(ACritical)
(rr,rc) = size(ARegular)

# Some matrices are not in double precision. We need to put it into double()
# precision to use sparse.

ATumor    = spzeros(double(ATumor)) #certo
ACritical = spzeros(double(ACritical))
ARegular  = spzeros(double(ARegular))

#
# Extract the upper & lower bounds from the prescription matrix
#

TUB = TGoal * (1 + p)
TLB = (TGoal * (1 - p)) + tol


#
# Calculate the weight omega that guarantees the tumor uniformity.
#
# The bound presented in the papaer leads to numerical instability
# with the interior point method always thinking that it is infeasible.
# This value seems to work well.
#

kappa = (norm(TLB, Inf))
omega = kappa / tol

#
# Build the data for the solver
#
#
#                Ux = u
#
# ( A_T    0    0    0 ) (  x   )    ( TUB )
# (-A_T   -L    0    0 ) (\alpha)    (-TUB )
# ( A_C    0  -U_c   0 ) (\beta ) .== ( CUB )
# ( A_R    0    0  -U_r) (\gamma)    ( RUB )
#
#
#     Lower bounds             Upper bounds
#
# (  x   )    ( inf )        (  x   )    (  0  )
# (\alpha)    ( TLB )        (\alpha)    (  0  )
# (\beta ) <= ( inf )        (\beta ) >= (-CUB )
# (\gamma)    ( inf))        (\gamma)    (  0  )
#
#

  
  #se quando não for um é o outro, da para reduzir isso para um if, logo se é Avarage entra se não já é Absolute coloque tudo isso dentro de oura função. Boas práticas!! 
if Analysis == "Average" #certo

  U  = spzeros([ ATumor    zeros(tr, tr) zeros(tr, gr) zeros(tr, rr);
                -ATumor      -I          zeros(tr, gr) zeros(tr, rr);
                 ACritical zeros(gr, tr)    -I         zeros(gr, rr);
                 ARegular  zeros(rr, tr) zeros(rr, gr)     -I       ])
  u  = spzeros([TUB; -TLB; CUB; RUB])
  lb = spzeros([zeros(tc, 1); zeros(tr, 1); -CUB; zeros(rr,1)])
  ub = spzeros([inf * ones(tc, 1); TLB; inf * ones(gr, 1); inf * ones(rr, 1)])
  c  = spzeros([zeros(tc,1); omega * ones(tr, 1); ones(gr, 1); ones(rr, 1)])

elseif Analysis == "Absolute"

  U  = spzeros([ ATumor    zeros(tr, 1) zeros(tr, 1) zeros(tr, 1);
                -ATumor     -ones(tr,1) zeros(tr, 1) zeros(tr, 1);
                 ACritical zeros(gr, 1)  -ones(gr,1) zeros(gr, 1);
                 ARegular  zeros(rr, 1) zeros(rr, 1) -ones(rr,1)])
  u  = spzeros([TUB; -TLB; CUB; RUB;])
  lb = spzeros([zeros(tc, 1); 0; -minimum(CUB); 0])
  ub = spzeros([Inf * ones(tc, 1); minimum(TLB); Inf; Inf])
  c  = spzeros([zeros(tc,1); omega; 1; 1]')

end


# Call CPLEX solver & set some parameters:
# Time limit of 900s
# Optimality tolerance = 1e^-04

model = Model(CPLEX.Optimizer)

if Method == "Simplex"                                  #O mesmo para issa parte, coloque dentro de outra função e reduza para um if 
#    opt = cplexoptimset["cplex"]
#    opt.timelimit = 900
#    opt.lpmethod = 1
#    opt.simplex.tolerances.optimality = 1e-04
#    (x,~,~,output) = cplexlp(c, U, u, [], [],lb,ub,[],opt)
elseif Method == "Interior point"
    opt = cplexoptimset["cplex"]
    opt.timelimit = 900
    opt.lpmethod = 4
    opt.simplex.tolerances.optimality = 1e-04
    (x,~,~,output) = cplexlp(c, U, u, [], [],lb,ub,[],opt)
end

TotalTime = etime(clock,StartTime)

#
# Partition the solution up into the different types of
# decision variables.
#

if Analysis == "Average"                                                    #Da para reduzir tudo isso em um if também, faça o mesmo, crie uma função use o primeiro if e retorne caso verdade se não ele vai sair e já executar o proximo                                                              #se for retorne as variaveis se não aplica o que tem no elseif, outra coisa use if else não tem necessidade                                                                    #de usar else if
  alpha = x[tc+1:tc+tr]                                                     #se acha que a informação de bater no nome do elseif relevante coloque no titulo dessa função, tudo isso vai deixar seu código muito melhor!
  beta  = x[tc+tr+1:tc+tr+gr]
  gamma = x[tc+tr+gr+1:tc+tr+gr+rr]
  x     = x[1:tc]

elseif Analysis == "Absolute"

  alpha = x[tc+1]
  beta  = x[tc+2]
  gamma = x[tc+3]
  x     = x[1:tc]

end

A = [ATumor ACritical ARegular]'
RadLevel = A*x

#Fazer retornar um logfile

end


#
# We now supply some usefull information from the solution.
# Start by opening the window for the analysis to appear.
#

#AnalysisWindow

#
# Display the type of Analysis being conducted.
#

#if Analysis .== "Average"
#  set(findobj("Tag','AnalType'),'String','An Average Analysis was used.")
#else()
#  set(findobj("Tag','AnalType'),'String','An Absolute Analysis was used.")
#end

#
# Display Information about whether the tumor is recieveing an
# appropriate dose.
#

#if Analysis .== "Average"

#  if norm(alpha,1) .> tol
#	set(findobj("Tag','AnalTumor'),'String", ...
#	   sprintf("On average, the Tumors Prescription is not attainable."))
#  else()
#	AverageTumorDose = norm(ATumor*x,1)/tr
#	set(findobj("Tag','AnalTumor'),'String", ...
#	   sprintf("The average tumor dose is within the prescribed uniformity setting. The average tumor dose is #2.2f.",AverageTumorDose))
#  end

#else()

#  if alpha .> tol
#        set(findobj("Tag','AnalTumor'),'String", ...
#           sprintf("The Tumors Prescription is not attainable."))
#  else()
#        MinDose = min(ATumor*x)
#	MaxDose = max(ATumor*x)
#        set(findobj("Tag','AnalTumor'),'String", ...
#           sprintf("The tumor dose is within the prescribed uniformity setting. The min. (max.) tumor dose is #2.2f [#2.2f].",MinDose,MaxDose))
# end

#end

#
# Display How the Critical Structures fair under our plan.
#

#if Analysis .== "Average"

#  AverageCritDesc = mean(beta)
# if AverageCritDesc .> tol
#	set(findobj("Tag','AnalCrit'),'String", ...
#	   sprintf("The Average Critical Structure Dose is Over its Prescribed Dose by #2.2f",AverageCritDesc))
#  else()
# 	set(findobj("Tag','AnalCrit'),'String", ...
#	   sprintf("The Average Critical Structure Dose is Under its Prescribed Dose by #2.2f",abs(AverageCritDesc)))
#   end

#else()

#  if beta .> tol
#	set(findobj("Tag','AnalCrit'),'String", ...
#	   sprintf("The Critical Structures Must Recieve an Excessive Dose of at least #2.2f",beta))
# else()
#	set(findobj("Tag','AnalCrit'),'String", ...
#           sprintf("The Critical Structures Recieved at least #2.2f Under Their Prescribed Dose.",abs(beta)))
 # end

#end

#
# State How We Did on The Regular Tissue
#

#if (Analysis .== "Average") && (rr ~= 0)

#  AverageRegDesc = mean(gamma)

#  if AverageRegDesc .> tol
#	set(findobj("Tag','AnalHot'),'String", ...
#	   sprintf("On Average, the Restricted Regular Tissue is Above its Prescribed Dose by #2.2f.",AverageRegDesc))
#  elseif AverageRegDesc .< tol
#        set(findobj("Tag','AnalHot'),'String", ...
#           sprintf("On Average, the Restricted Regular Tissue is Below its Prescribed Dose."))

#  end

#elseif rr ~= 0

#  if gamma .> tol
#	set(findobj("Tag','AnalHot'),'String", ...
#	   sprintf("The Restricted Regular Tissue is Guaranteed to be Above its Prescribed Dose by #2.2f.",gamma))
#  elseif gamma .< tol
#        set(findobj("Tag','AnalHot'),'String", ...
#           sprintf("The Restricted Regular Tissue is Guaranteed to be Below its Prescribed Dose."))

#  end

#else()

#  set(findobj("Tag','AnalHot'),'String", ...
#	sprintf("There were no restrictions on the Regular Tissue."))

#end

#
# State the Solver Type
#

#if Method .== "Interior Point"
#   set(findobj("Tag','AnalSolver'),'String", ...
#	sprintf("An Interior Point Method Was Used."))
#else()
#   set(findobj("Tag','AnalSolver'),'String", ...
#	sprintf("A Simplex Method Was Used."))
#end

#
# State the Time required by the solver.
#

#set(findobj("Tag','AnalTime'),'String", ...
#	sprintf("The solution time was #4.2f seconds.",TotalTime))

#
# State any notes.
#


#if maximum(RadLevel) <= 1.1*maximum(TUB)
#  set(findobj("Tag','AnalNotes'),'String", ...
#	sprintf("There are no hot spots."))

#else()

#  set(findobj("Tag','AnalNotes'),'String", ...
#  sprintf("There are hot spots. You may want to make further restrictions."))

#end

#end
