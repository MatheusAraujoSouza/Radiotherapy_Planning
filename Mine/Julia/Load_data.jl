using MAT
function Load_data(Path, dataID_TGoal, dataID_CUB, dataID_RUB)
# This function load the data from the database.
#
# Each structure of the patient on the database has an pencil-beam dose
# matrix A & its respective prescription dose. Three categories are
# created for those structures: tumor; critical; & regular tissue.
#
# **OBSERVATION** - "Regular" stands for "Good tissue"
#
# So; supposing there are six structures in total: two "tumors"; three
# "criticals"; & one "regular" structures. Those matrices are
# respectively:
#
# Tumor    -> AT_1 & AT_2
# Critical -> AC_1; AC_2 & AC_3
# Regular  -> AR_1
#
# And its prescriptions vectors [where P stands for "prescription"]:
#
# Tumor    -> PT_1 & PT_2
# Critical -> PC_1; PC_2 & PC_3
# Regular  -> PR_1
#
# In order to compile this data compatibly with the math model; the lines of
# each prescription vector & pencil-beam matrix are concatenaded:
#
# ATumor = [AT_1] , ACritical = [AC_1] , ARegular = [AR_1]
#          [AT_2]               [AC_2]
#                               [AC_3]
#
# TGoal = [PT_1] , CUB = [PC_1] , RUB = [PR_1]
#         [PT_2}         [PC_2]
#                        [PC_3
#
###########################################################################

#DADOS PARA TESTAR APENAS

direct = matread(Path)

#
# Load the TGoal prescription vector & the ATumor matrix
#

println("passou o matread")

(~ , size_p) = size(direct["problem"]["dataID"]) #certo
ATumor = Array{Float64}(undef, 0, Int(direct["data"]["misc"]["size"]))
TGoal = Array{Float64}(undef, 0, 1)
for i = dataID_TGoal
    if i != 0 #certo
        (n_rows_A ,~) = size(direct["data"]["matrix"]["A"][i]) #certo
        ATumor = [ATumor; direct["data"]["matrix"]["A"][i]] #PARECE QUE FAZER ATumor = [] NAO FUNFA
        for j = 1:size_p
            if direct["problem"]["dataID"][j] == i && direct["problem"]["IsConstraint"][j] == 1
                TGoal = [TGoal; ones(n_rows_A,1)*direct["problem"]["Objective"][j]]
            end
        end
    end
end

#
# Load the CUB prescription vector & the ACritical matrix
#

ACritical = Array{Float64}(undef, 0, Int(direct["data"]["misc"]["size"]))
CUB = Array{Float64}(undef, 0, 1)
for i = dataID_CUB
    if i != 0 #certo
        (n_rows_A ,~) = size(direct["data"]["matrix"]["A"][i]) #certo
        ACritical = [ACritical; direct["data"]["matrix"]["A"][i]] #PARECE QUE FAZER ATumor = [] NAO FUNFA
        for j = 1:size_p
            if direct["problem"]["dataID"][j] == i && direct["problem"]["IsConstraint"][j] == 1
                CUB = [CUB; ones(n_rows_A,1)*direct["problem"]["Objective"][j]]
            end
        end
    end
end

# Load the RUB prescription vector & the ARegular matrix

ARegular = Array{Float64}(undef, 0, Int(direct["data"]["misc"]["size"]))
RUB = Array{Float64}(undef, 0, 1)
for i = dataID_RUB
    if i != 0 #certo
        (n_rows_A ,~) = size(direct["data"]["matrix"]["A"][i]) #certo
        ARegular = [ARegular; direct["data"]["matrix"]["A"][i]] #PARECE QUE FAZER ATumor = [] NAO FUNFA
        for j = 1:size_p
            if direct["problem"]["dataID"][j] == i && direct["problem"]["IsConstraint"][j] == 1
                RUB = [RUB; ones(n_rows_A,1)*direct["problem"]["Objective"][j]]
            end
        end
    end
end

return (TGoal, RUB, CUB, ATumor, ACritical, ARegular)

end
