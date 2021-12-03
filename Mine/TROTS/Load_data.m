function [TGoal, RUB, CUB, ATumor, ACritical, ARegular] = Load_data(Path, dataID_TGoal, dataID_CUB, dataID_RUB)
% This function load the data from the database.
%
% Each structure of the patient on the database has a pencil-beam dose
% matrix A and its respective prescription dose. Three categories are
% created for those structures: tumor, critical, and regular tissue.
%
% OBSERVATION: "Regular" stands for "Good tissue"
%
% So, supposing there are six structures in total: two "tumors", three
% "criticals", and one "regular" structure. Those matrices are,
% respectively:
% 
% Tumor    -> AT_1 and AT_2
% Critical -> AC_1, AC_2 and AC_3 
% Regular  -> AR_1
%
% And its prescriptions vectors (where P stands for "prescription"):
%
% Tumor    -> PT_1 and PT_2
% Critical -> PC_1, PC_2 and PC_3 
% Regular  -> PR_1
%
% In order to compile this data compatibly with the math model, each
% prescription vector's rows and pencil-beam matrices are concatenated:
%
% ATumor = [AT_1] , ACritical = [AC_1] , ARegular = [AR_1]
%          [AT_2]               [AC_2]
%                               [AC_3]
%
% TGoal = [PT_1] , CUB = [PC_1] , RUB = [PR_1]
%         [PT_2}         [PC_2]
%                        [PC_3]
%
%
% ***Some matrices are not in double precision. We need to put it into double
% precision to use sparse.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% The "easy cases" has no regular tissue considered. In light of this, we
% initially make RUB = [].
%

RUB=[];

%
% Load the path of the case (prostate, liver or h_a_n)
%

load(Path);

%
% Load the TGoal prescription vector and the ATumor matrix
%

sum_r_A = 0; 
[~ , size_p] = size(problem);
ATumor = [];
for i = dataID_TGoal
    if i ~= 0
        [n_rows_A ,~] = size(data.matrix(i).A);
        ATumor = [ATumor; sparse(double(data.matrix(i).A))];
        for j = 1:size_p
            if problem(j).dataID == i && problem(j).IsConstraint == 1      
                TGoal(sum_r_A + 1:sum_r_A + n_rows_A ,1) = ones(n_rows_A,1)*problem(j).Objective;
                sum_r_A = sum_r_A + n_rows_A;
            end
        end
    end
end

%
% Load the CUB prescription vector and the ACritical matrix
%

sum_r_A = 0; 
ACritical = [];
for i = dataID_CUB 
    if i ~= 0
        [n_rows_A ,~] = size(data.matrix(i).A);
        ACritical = [ACritical; sparse(double(data.matrix(i).A))];
        for j = 1:size_p
            if problem(j).dataID == i && problem(j).IsConstraint == 1      
                CUB(sum_r_A + 1:sum_r_A + n_rows_A ,1) = ones(n_rows_A,1)*problem(j).Objective;
                sum_r_A = sum_r_A + n_rows_A;
            end
        end
    end
end

% Load the RUB prescription vector and the ARegular matrix

sum_r_A = 0;
ARegular = [];
for i = dataID_RUB
    if i ~= 0
        [n_rows_A ,~] = size(data.matrix(i).A);
        ARegular = [ARegular; sparse(double(data.matrix(i).A))];
        for j = 1:size_p
            if problem(j).dataID == i && problem(j).IsConstraint == 1      
                RUB(sum_r_A + 1:sum_r_A + n_rows_A ,1) = ones(n_rows_A,1)*problem(j).Objective;
                sum_r_A = sum_r_A + n_rows_A;
            end
        end
    end
end
end

