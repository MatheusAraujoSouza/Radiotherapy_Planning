# Radiotherapy Planning
This is the biggest repository in this profile. It contains a bunch of things regarding radiotherapy planning. In the next paragraphs, you will find a brief description of each folder above.

# CORT
CORT is a database created by the Massachusetts General Hospital for radiotherapy planning. The folder contains two implementations showed in the paper "Data_set_harvard". You can download the database at https://www.dropbox.com/sh/gny44brgsfdejr7/AAAJi0B-ZGJs9YNTiS2Ju4iya?dl=0

# Holder Routines
These are the implementations made by Allen Holder in his publication "Designing Radiotherapy Plans with Elastic Constraints and Interior Point Methods", where the data used are stored at the following link: https://www.ime.unicamp.br/~aurelio/dose/

# TG119
The purpose of the task group 119 were to define standard IMRT planning “problems” that physicists can use to test the accuracy of their radiotherapy treatments. The folder contains some papers and a MATLAB routine to import this dada created by Guilherme Giacomini.

# TROTS
As CORT, TROTS is also a database for the same purpose. Although, TROTS is way bigger than CORT, and has some interesting routines associated. The folder contains the TROTS paper and routine. More informations at the following link https://sebastiaanbreedveld.nl/trots/

# Mine
Here is where the fun begins. This folder contains my implementations of the Holder model using the TROTS database AND the holder data itself. The code in MATLAB is already done, but the code made in Julia doesn't. See "report_1.pdf" for more information about it and without any math. But if you want some, see "report_2.pdf". Both reports are in Portuguese.
