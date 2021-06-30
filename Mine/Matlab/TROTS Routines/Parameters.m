% This is the main routine to run the model.
%
% All the parameters of the Solver function are changed here automatically
% with the purpose of doesn't waste time doing it manually.
% However, it still needed a "logfile" in the Solver function to save each
% result.
% 
% The parameters of the solver are:
% 
% - Path - Corresponds to the case in question. The parameter carry the
% path where the data of the respective case is saved. The cases are:
% Prostate, liver, head and neck, and head and neck alterated.
%
% - Analysis - Corresponds to the solution analysis, which can be the
% absolute or the average one.
%
% - Method - Corresponds to the resolution method of the LP, which can be
% Simplex or Interior Points.
%
% - p - Corresponds to the upper and lower bound variation of the prescription
% dose in tumor. Values ??usually range between 2% and 15%
%
% - Total_number_of_tests - Corresponds to the number of tests that will be done with
% each case. See that the
% database has 30 prostate, 10 liver, and 30 head and neck tests
% available for testing.
%
% - structures of each case - This one can be a little tricky.
% For each case, it is necessary to choose the patient structures that will be
% considered in the method. These are: PTV, OAR and regular tissues (RT).
% In the database used, each of these structures for each case
% has a different ID. In light of this, a matrix called Volume is created,
% where the Id's (different from zero, that is, if an ID is zero, the algorithm ignores it)
% of PTV are in the first line, those of OAR are in the
% second and the RT are in the third. An example of two different volumes is
% shown in the prostate case.
%
% PROSTATE:
%   Volume_1: [1 13  0 0 ;   <-- PTV
%              3  4  5 0 ;   <-- OAR
%             16  0  0 0 ;]  <-- RT
%
%   Volume_2: [1 13 18  0 ;  <-- PTV
%              3  4  5 10 ;  <-- OAR
%             16  0  0  0 ;] <-- RT
%
% This second volume show that it was added a critical structure with the
% ID = 10, and the tumor structure which the ID is 18 was removed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    total_number_of_test = 1;

    %
    % The number of each case was changed from "01" to "1" in order to
    % facilitate the following loops 
    %
    for current_test = 1:total_number_of_test
        prostate = convertCharsToStrings(['C:\Principal\IC\Estruturas_e_rotinas\TROTS\Dados\Prostate_CK\Prostate_CK_' num2str(current_test) '.mat']);
        liver = convertCharsToStrings(['C:\Principal\IC\Estruturas_e_rotinas\TROTS\Dados\Liver\Liver_' num2str(current_test) '.mat']);
        h_a_n = convertCharsToStrings(['C:\Principal\IC\Estruturas_e_rotinas\TROTS\Dados\Head-and-Neck\Head-and-Neck_' num2str(current_test) '.mat']);
        h_a_n_alt = convertCharsToStrings(['C:\Principal\IC\Estruturas_e_rotinas\TROTS\Dados\Head-and-Neck_Alt\Head-and-Neck-Alt_' num2str(current_test) '.mat']);

        folder = [prostate liver h_a_n h_a_n_alt];
        for Path = folder
            if Path == prostate
                Matrix_ID = [1 0 0;
                             3 0 0;
                             16 0 0];

            elseif Path == liver
                Matrix_ID = [2 0 0 0 0 0;
                             3 4 5 6 7 8;
                             9 0 0 0 0 0];

            elseif Path == h_a_n
                Matrix_ID = [1 9 18 0 0 0 0 0 0 0 0 0 0 0;
                             2 3 12 13 14 16 18 20 22 25 27 29 31 33;
                             8 0 0 0 0 0 0 0 0 0 0 0 0 0];

            elseif Path == h_a_n_alt
                Matrix_ID = [1 9 18 0 0 0 0 0 0 0 0 0 0 0;
                             2 3 12 13 14 16 18 20 22 25 27 29 31 33;
                             8 0 0 0 0 0 0 0 0 0 0 0 0 0];         

            end
                    for Analysis = ["Absolute","Average"]
                        for Method = ["Interior point","Simplex"]
                            for p = [0.02 0.07 0.10]
                                Solver(Path, Analysis, Method, p, Matrix_ID(1,:), Matrix_ID(2,:), Matrix_ID(3,:))
                            end
                        end
                    end
        end
    end