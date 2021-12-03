% This is the main routine to run the model.
% All the parameters of the Solver function are changed here automatically
% with the purpose of doesn't waste time doing it manually.
% However, it still needed a "logfile" in the Solver function to save each
% result.
%  
% The parameters of the solver are:
%  
% - Path - Corresponds to the case in question. The parameter carries the
% path where the data of the respective case is saved. The cases are:
% Prostate, liver, head and neck, and head and neck altered.
% 
% - Analysis - Corresponds to the solution analysis, which can be the
% absolute or the average one.
% 
% - Method - Corresponds to the resolution method of the LP, which can be
% Simplex or Interior Points.
% 
% - p - Corresponds to the upper and lower bound variation of the prescription
% dose in tumor. Values usually range between 2% and 15%
% 
% - Total_number_of_tests - Corresponds to the number of tests to be done with
% each case. See that the database has 30 prostate, 10 liver, and
% 30 head and neck tests available for testing.
% 
% - Structures of each case - This one can be a little tricky.
% For each case, it is necessary to choose the patient structures that will be
% considered in the method. These are: PTV, OAR and regular tissues (RT).
% In the database used, each of these structures for each case
% has a different ID. In light of this, some vectors carrying those ID's were
% created. Finally, if the ID == 0, there's no structure. Example:
% 
%     dataID_TGoal = [1 9];        <-- PTV
%     dataID_CUB = [2 14 16 22];   <-- OAR
%     dataID_RUB = 8;              <-- RT
%
%
% OBSERVATION: The number of each case was changed from "01" to "1" in 
% order to facilitate the following loops, which find the 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Counter of tests
%

Counter = 0;

%
% Preallocate the table which saves all the results:
%

vartype = {'int8', 'string', 'string', 'string', 'string', 'int16', 'int16', 'int8', 'string', 'string', 'int16', 'double', 'string', 'string', 'string', 'string'};
Table = table('Size', [36,16], 'VariableTypes', vartype);
Table.Properties.VariableNames = {'Counter', 'Case', 'Dificulty', 'Analysis', 'Method', 'Number_of_constraints', 'Number_of_variables', 'Solver_status', 'Solver_status_string', 'Solver_message', 'Iterations', 'Time', 'Tumor_dose', 'Critical_tissue_dose', 'Regular_tissue_dose', 'Hot_spots'};


%
% The total number of tests regarding each case.
%

total_number_of_test = 1;

for current_test = 1:total_number_of_test
    prostate = "C:\Principal\IC\Dados_grandes_IC\Dados_TROTS\Prostate_CK\Prostate_CK_" + num2str(current_test) + ".mat";
    liver = "C:\Principal\IC\Dados_grandes_IC\Dados_TROTS\Liver\Liver_" + num2str(current_test) + ".mat";
    h_a_n = "C:\Principal\IC\Dados_grandes_IC\Dados_TROTS\Head-and-Neck\Head-and-Neck_" + num2str(current_test) + ".mat";
%   h_a_n_alt = "C:\Principal\IC\Estruturas_e_rotinas\TROTS\Dados\Head-and-Neck_Alt\Head-and-Neck-Alt_" + num2str(current_test) + ".mat";

    folder = [prostate liver h_a_n];
        for Path = folder
            for Difficulty = "hard"
                if Path == prostate
                    Case = "Prostate";
%                     if Difficulty == "easy"
%                         dataID_TGoal = 1;
%                         dataID_CUB = 3;
%                         dataID_RUB = 0;
%                         
%                     elseif Difficulty == "moderate"
%                         dataID_TGoal = [1 18];
%                         dataID_CUB = [3 4];
%                         dataID_RUB = 16;  

                    if Difficulty == "hard"
                        dataID_TGoal = [1 13 18];
                        dataID_CUB = [2 3 4 10];
                        dataID_RUB = 16;
                    end
                        
                elseif Path == liver
                    Case = "Liver";
%                     if Difficulty == "easy"
%                         dataID_TGoal = 2;
%                         dataID_CUB = 3;
%                         dataID_RUB = 0;
%                         
%                     elseif Difficulty == "moderate"
%                         dataID_TGoal = 2;
%                         dataID_CUB = [3 4 5];
%                         dataID_RUB = 9;

                    if Difficulty == "hard"
                        dataID_TGoal = 2;
                        dataID_CUB = [3 4 5 6 8];
                        dataID_RUB = [9 24];
                    end
                    
%                 elseif Path == h_a_n
%                     Case = "Head and Neck";
% %                     if Difficulty == "easy"
% %                         dataID_TGoal = 1;
% %                         dataID_CUB = 22;
% %                         dataID_RUB = 0;
% %                         
% %                     elseif Difficulty == "moderate"
% %                         dataID_TGoal = [1 9];
% %                         dataID_CUB = [2 14 16 22];
% %                         dataID_RUB = 8;
% 
%                     if Difficulty == "hard"
%                         dataID_TGoal = [1 9];
%                         dataID_CUB = [2 3 12 13 14 16 18 20 22 25 35]; 
%                         dataID_RUB = 8;
                end
            end
                    for Analysis = ["Absolute","Average"]
                        for Method = ["Interior point","Automatic"]
                            for p = 0.02
                                Counter = Counter +1
                                Table = Solver_consertado(Difficulty, Table, Counter, Case, Path, Analysis, Method, p, dataID_TGoal, dataID_CUB, dataID_RUB);
                            end
                        end
                    end
        end
end
%
% Export my MATLAB's table to a csv file       
%

writetable(Table, 'Resultados_hard.csv')