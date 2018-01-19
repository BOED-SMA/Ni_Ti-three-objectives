clear
clc



% 1) DEFINE INPUT SPACE 
VR1_SPACE=0:0.005:0.1;  %Variable 1
VR2_SPACE=50.2:0.1:51.2; % Variable 2
VR3_SPACE=50.2:0.1:51.2; % Variable 3


nV=2; % No OF VARIABLES

% 2) DEFINE OBJECTIVE FUNCTIONS
objectfun = @(y) [abs(y(:, 2)-303),  abs(y(:, 2)-y(:, 3)-40), y(:, 4)]; % Define Objective functions
constrfun = @(y) [y(:, 3)]; % Define Objective functions
boundpoint = [83, 136, 360];
nO=3; % NUMBER OF OBJECTIVES
nC=0; % NUMBER OF CONSTRAINTS


% 3) DEFINE nE and nI
nE=10; % Preselected number of evaluations of he BOED process. Each evaluation correspods to an experiment
nI=1; % Number of randomly selected experiments
% if want to parallelize the micromechanical simulations you have to either
% initiate a local parpool HERE (BEFORE THE CODE BEGINNING) by using the command : parpool(workers_number) or you
% have to initiate a parpool in the supercomputer using the apropriate
% commands. For the TAMU ada you have to use:
% matlabsubmit -t 150:00 -w 4 -p 1 -g -s 20 speed3b.m (or alternativly matlabsubmit -t 150:00 -w 4 -p 1 -g -x "-x", -x means utilize all the availible cores of 1 node)in order to initiate the analysis. 
% THe given examples will initiate 4 workers which are loaded in 1 node each,using 20 cores and also use gpu nodes
% THE PARPOOL IS USED ONLY TO EXECUTE THE 4-5 NEEDED MICROMECHANICAL
% SIMULATIONS UNDER DIFFERENT THERMO-MECHANICAL LOADS AT THE SAME TIME

% CODE BEGINNING 
XSPACE=DXSPACE(VR1_SPACE,VR2_SPACE,VR3_SPACE);
[X, Y, XSPACE, qnum, dnum] = DEVELOP_INITIAL_IO_DATABASE(objectfun,constrfun,XSPACE,nI,nV,nO,nC); % DEVELOP INITIAL INPUT-OUTPUT DATABASE
[X, Y, XSPACE, PARETO_FRONT, qnum, dnum]= MAIN_BOED_ITERATIVE_LOOP(objectfun,constrfun,XSPACE,X,Y,nI,nE,qnum, dnum,nV,nO,nC, boundpoint); % MAIN BOED ITERATIVE LOOP
% PLOT_FIGURES(X,Y,XSPACE,PARETO_FRONT,qnum, dnum)
%##############################################################################################################################################################################################################   
%##############################################################################################################################################################################################################      
%######################################################################################  MAIN BOED FUNCTIONS    ###############################################################################################
%##############################################################################################################################################################################################################      
%##############################################################################################################################################################################################################   
%##############################################################################################################################################################################################################  


function [X, Y, XSPACE, qnum, dnum] = DEVELOP_INITIAL_IO_DATABASE(objectfun,constrfun,XSPACE,nI,nV,nO,nC)
%qnum the number of quried rows of YSPACE 
%dnum the number of discarded data 888 999
sprintf(['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$','\n','  DEVELOP INITIAL INPUT-OUTPUT DATABASE  ','\n','$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'])    

%INITIALIZE VARIABLES AND MATRICES
qnum = 0;
dnum = 0;
Y = nan(1, nO+nC); %INITIALIZE Y MATRIX
X = nan(1, nV); %INITIALIZE X MATRIX

m=1;     
while m<=nI 
            
    UNEXPLORED_XSPACE=XSPACE(qnum+1:size(XSPACE,1)-dnum,2:nV+1);
    n=size(UNEXPLORED_XSPACE,1); % NUMBER OF REMAINING POTENTIAL EXPERIMENTS
    
    if (n==0) % if m reaches the size of XSPACE -  the number of discarded data +1 break the while loop
    sprintf(['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$','\n','  NO MORE INITIAL EXPERIMENTS CAN BE PERFORMED ','\n','$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'])    
    break
    end   
    clear data
    
    % SELECTION OF INPUT VARIABLES WITHOUT REPLACEMENT FROM THE PREDIFINED INPUT SPACE (XSPACE)
    INPUT_VALUES=datasample(UNEXPLORED_XSPACE,1);

      
    % CALCULATION OF OPERATIONAL OBJECTIVES
    [OUTPUT_VALUES, ERROR] = SYSTEMS_MODEL2(INPUT_VALUES,XSPACE,qnum, dnum,nV,nO,nC);
    [OBJ_CONSTR_VALUES]=CALCULATE_OBJECTIVES_CONSTRAINTS(objectfun,constrfun,OUTPUT_VALUES,ERROR,nV,nO,nC);
    
    % UPDATE INPUT-OUTPUT DATABASE
    [X,Y]=UPDATE_IO_DATABASE(OBJ_CONSTR_VALUES,X,Y,INPUT_VALUES,ERROR,nV,nO,nC);
    [XSPACE, dnum, qnum ] =SORT_XSPACE(XSPACE,dnum,qnum,INPUT_VALUES,ERROR,nV,nO,nC);                                            
    
    
    [idx, CASE_NUMBER]=IDENTIFY_IDX_AND_CASE_NUMBER(INPUT_VALUES,XSPACE,nV,nO,nC); % ONLY FOR PRINTING PURPOSES
    if ERROR==0
    sprintf(['RANDOMLY SELECTED EXPERIMENT: ',num2str(m),' FOR CASE NUMBER: ',num2str(CASE_NUMBER),' SUCCESFULLY COMPLETED'])
    m=m+1;
    else
    sprintf(['RANDOMLY SELECTED EXPERIMENT: ',num2str(m),' FOR CASE NUMBER: ',num2str(CASE_NUMBER),' NOT COMPLETED - PROBLEMATIC INPUT VALUE. SELECT NEXT MATERIAL TO TEST FROM THE REDUCED INPUT SPACE'])    
    end
        
end   

sprintf(['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$','\n','  INITIAL_DATABASE COMPLETED ','\n','$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'])

end

function [OBJ_CONSTR_VALUES]=CALCULATE_OBJECTIVES_CONSTRAINTS(objectfun,constrfun,OUTPUT_VALUES,ERROR,nV,nO,nC)


IMPORTANT_OUTPUT_VARIABLES=OUTPUT_VALUES(:,13:22); % THE MODEL OUTPUTS 22 VARIABLES, HOWEVER ONLY THE 13 TO 22 ARE USEFULL FOR THIS PROCESS. HENCE I AM SAVING THOSE VARIABLES IN THE IMPORTANT VARIABLES MATRIX

  if ERROR~=0
        OBJ_CONSTR_VALUES(1, 1:nO+nC) = 666; % DEFINE 666 IN ORDER TO RECOGNISE WHEN THERE WAS AN ERROR
  else
        OBJ_CONSTR_VALUES(1, 1:nO) = objectfun(IMPORTANT_OUTPUT_VARIABLES);
        if nC~=0
        OBJ_CONSTR_VALUES(1, nO+1:nO+nC) = constrfun(IMPORTANT_OUTPUT_VARIABLES);
        end
  end

end

function [X, Y] =UPDATE_IO_DATABASE(OBJ_CONSTR_VALUES,X,Y,INPUT_VALUES,ERROR,nV,nO,nC)

m=size(X,1);

if ERROR==0 % UPDATE ONLY WHEN THERE IS NO ERROR
    
    if m==1 & isnan(X(1,1))==1 % For the first time that X,Y are updated
        X(1,1:nV)=INPUT_VALUES(1:nV);
        Y(1, 1:nO+nC) = OBJ_CONSTR_VALUES;
    elseif m>=1 & isnan(X(1,1))==0
        X(m+1,1:nV)=INPUT_VALUES(1:nV);
        Y(m+1, 1:nO+nC) = OBJ_CONSTR_VALUES;
    end

end
   

end

function [X, Y,XSPACE, PARETO_FRONT, qnum, dnum]=MAIN_BOED_ITERATIVE_LOOP(objectfun,constrfun,XSPACE,X,Y,nI,nE,qnum, dnum,nV,nO,nC, referencepoint)
sprintf(['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$','\n','  BOED MAIN ITERATIVE LOOP  ','\n','$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'])

nT=size(XSPACE,1);
nRT=nT-nI-dnum; %  Remaining INPUT SPACE that can be covered during the evaluations of the BOED method
i=1;
 
while i<=nE  && i <=nRT    % IT WILL STOP WHEN THE STOPNUM IS REACHED OR WHEN THE ENTIRE SPACE IS COVERED. THE ENTIRE SPACE IS DEFINED AS THE INITIAL_SPACE-THE INITIAL RUNS-NOT ACCEPTABLE POINTS
   
    %MACHINE LEARNING AND SELECTOR STEPS - DEFINITION OF THE NEXT POINT TO TEST    
    boundpoint = referencepoint;
    [Y_MEAN, Y_STD]= MACHINE_LEARNING(X,Y,XSPACE, qnum, dnum,nV,nO,nC); % Construct Surrogate model based on the input-output data (X,Y)
    [INPUT_VALUES p_integral]= SELECTOR(X,Y, Y_MEAN, Y_STD,XSPACE, qnum, dnum, boundpoint,nV,nO,nC); % Select the next material to test (INPUT)
     
    % CALCULATION OF OPERATIONAL OBJECTIVES
    [OUTPUT_VALUES, ERROR] = SYSTEMS_MODEL2(INPUT_VALUES,XSPACE,qnum, dnum,nV,nO,nC);
    [OBJ_CONSTR_VALUES]=CALCULATE_OBJECTIVES_CONSTRAINTS(objectfun,constrfun,OUTPUT_VALUES,ERROR,nV,nO,nC);
    
    % UPDATE INPUT-OUTPUT DATABASE
    [X,Y]=UPDATE_IO_DATABASE(OBJ_CONSTR_VALUES,X,Y,INPUT_VALUES,ERROR,nV,nO,nC);
    [XSPACE, dnum, qnum ] =SORT_XSPACE(XSPACE,dnum,qnum,INPUT_VALUES,ERROR,nV,nO,nC);                                                
    
    % UPDATE AUXILARY VARIABLES
     nRT=nT-nI-dnum;
    
    [idx, CASE_NUMBER]=IDENTIFY_IDX_AND_CASE_NUMBER(INPUT_VALUES,XSPACE,nV,nO,nC); % ONLY FOR PRINTING PURPOSES
    if ERROR==0
    sprintf(['EVALUATION: ',num2str(i),' FOR CASE NUMBER: ',num2str(CASE_NUMBER),' SUCCESFULLY COMPLETED'])
    i=i+1;
    else
    sprintf(['EVALUATION ',num2str(i),' FOR CASE NUMBER: ',num2str(CASE_NUMBER),' PROBLEMATIC INPUT VALUE, PROCEED TO NEXT CASE_NUMBER'])    
    end
    

clear PARETO_FRONT
[PARETO_FRONT, ~]=Find_pareto_front_multi2(Y);
[~, I] = sort(PARETO_FRONT(:, 2));
PARETO_FRONT = PARETO_FRONT(I, :);


% WRITE PARETO FRONT IN TXT FILE --- ONLY FOR DEBUGGING PURPOSES
root_directory=pwd;
OUTNAME1=[root_directory,'/OUTPUT_PARETO_FRONT.txt'];
if isempty(PARETO_FRONT)==1;
PARETO_FRONT=zeros(1,2);
CASE_NO_OF_PARETO_FRONT_POINTS=0;
else
    
INDEX2=size(PARETO_FRONT);
for l=1:INDEX2(1)
[row column]=find(Y(:,:)==PARETO_FRONT(l,:));
indices_we_want(l)=row(1);
end
    
    
CASE_NO_OF_PARETO_FRONT_POINTS=XSPACE([indices_we_want],1);
end
if ERROR==0
    K=i-1;
else
    K=i;
end
INDEX2=size(PARETO_FRONT);
fid = fopen(OUTNAME1,'a');
Res=fprintf(fid, '%i \n',K); % IF YOU SEE IN THE TXT MULTIPLE CASES OF SAME I AND DIFFERENT CASE NUMBER. THAT MEANS THAT,THE INPUT OF THAT CASE NUMBER WASNTS VALID
for j=1:INDEX2(1);
Res=fprintf(fid, '%f %f %f \n',CASE_NO_OF_PARETO_FRONT_POINTS(j),PARETO_FRONT(j,:));
end 
fclose(fid);
    
end
sprintf(['$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$','\n','  BOED PROCESS COMPLETED  ','\n','$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'])




end

function [Y_MEAN, Y_STD]= MACHINE_LEARNING(X,Y,XSPACE, qnum, dnum,nV,nO,nC)
XSPACER=XSPACE(:,2:nV+1);

MEAN=mean(XSPACER);
STD=std(XSPACER);

X_NORM = (X-MEAN)./(ones(size(X, 1),1)*STD);
left_idx = qnum+1:size(XSPACER, 1)-dnum;% from current number to the end of xspace
X_TEST=XSPACER(left_idx, :);
X_TEST_NORM=(X_TEST-MEAN)./(ones(size(X_TEST, 1),1)*STD);
[Y_MEAN, Y_STD] = GPR_SURROGATE_MODEL(X_NORM,Y,X_TEST_NORM);
end

function [Y_MEAN, Y_STD] = GPR_SURROGATE_MODEL(X, Y, X_TEST)
% ymean,  size R*C 
% ysd,  Standard deviation, sigma. ysd^2 = variance 
C = size(Y, 2);
R = size(X_TEST, 1);
Y_MEAN = zeros(R, C);
Y_STD = zeros(R, C);
for n = 1:C
    gprMd = fitrgp(X, Y(:,n));%, 'BasisFunction', 'linear');%
    for m = 1:R
        [Y_MEAN(m, n), Y_STD(m, n), ~] = predict(gprMd, X_TEST(m, :));% 
    end
end
end

function [INPUT p_integral]= SELECTOR(X,Y,Y_MEAN, Y_STD,XSPACE, qnum, dnum, boundpoint,nV,nO,nC)
XSPACER=XSPACE(:,2:nV+1);
left_idx = qnum+1:size(XSPACER, 1)-dnum;% from current number to the end of xspace
X_TEST=XSPACER(left_idx, :);
m=size(left_idx,2);
[PARETO_FRONT, ~] = Find_pareto_front_multi2(Y);% value of y object, sorted
max_integral = 0;

for i = 1:m
    
            
%     p_integral = EHVI(Y_MEAN(i,:), Y_STD(i,:), PARETO_FRONT, boundpoint);
    p_integral = EHVI_multi_objective_prob(Y_MEAN(i,:), Y_STD(i,:), PARETO_FRONT, boundpoint);
    
    if p_integral >= max_integral
        max_integral = p_integral;
        INPUT_VALUES=X_TEST(i,:);
%         NI_COMPO=X_TEST(i,1);
%         VF=X_TEST(i,2);
        
        [idx, CASE_NUMBER]=IDENTIFY_IDX_AND_CASE_NUMBER(INPUT_VALUES,XSPACE,nV,nO,nC);;
        next_idx = idx;
    end
    
end

INPUT=XSPACER(next_idx,:);


end

function p_integral = EHVI(mu, sigma, PARETO_FRONT, extreme_point)

Phi = @(s) 1/2*(1+erf(s/sqrt(2)));
bigpsi = @(a, b, mu, sigma) sigma*normpdf((b-mu)/sigma)+(a-mu)*Phi((b-mu)/sigma);
n = size(PARETO_FRONT, 1);

part1 = zeros(1, n);
part2 = zeros(1, n);

y = zeros(n+2, 2);
y(1, :) = [extreme_point(1), -inf];
y(2:end-1, :) = PARETO_FRONT;
y(end, :) = [-inf, extreme_point(2)];

% y is the objective front array
for i = 2:size(y, 1)
    part1(i-1) = Phi((y(i,1)-mu(1))/sigma(1))*(y(i-1, 1)-y(i, 1))*...
        bigpsi(y(i, 2), y(i,2), mu(2), sigma(2));
    part2(i-1) = ( bigpsi(y(i-1, 1), y(i-1, 1), mu(1), sigma(1)) - ...
        bigpsi(y(i-1, 1), y(i, 1), mu(1), sigma(1)) )*...
        bigpsi(y(i, 2), y(i, 2), mu(2), sigma(2));
end
part1(end) = 0;
p_integral = sum(part1)+sum(part2);

end


%##############################################################################################################################################################################################################   
%##############################################################################################################################################################################################################      
%#######################################################################  MODEL (IN THIS CASE SMA MICROMECHANICAL MODEL)   ####################################################################################
%##############################################################################################################################################################################################################      
%##############################################################################################################################################################################################################   
%##############################################################################################################################################################################################################  


function [OUTPUT_VALUES, ERROR] = SYSTEMS_MODEL2(INPUT_VALUES,XSPACE,qnum, dnum,nV,nO,nC)

% DEFINE YOUR VARIABLES NAMES
NI_COMPO=INPUT_VALUES(1);
VF=INPUT_VALUES(2);

% !!!!!!!!!!! RUN YOUR MODEL !!!!!!!!!!!!!!!!!
MODEL = xlsread('SMA_MICROMECHANICAL_MODEL_SIMULATOR.xlsx');  % THIS SIMULATES THE MICROMECHANICAL MODEL
idx = 0;
for m = 1:size(MODEL, 1)
    if MODEL(m, 4) == NI_COMPO && MODEL(m, 9) == VF
        idx = m;
    end
end
if idx == 0
    error('something wrong')
end

OUTPUT_VALUES=MODEL(idx, :);  
if MODEL(idx, 13)~=888 & MODEL(idx,13)~=999
ERROR=0;
else
ERROR=1; % any value than zero will do
end
% !!!!!!!!!!! END OF MODEL !!!!!!!!!!!!!!!!!

end


%##############################################################################################################################################################################################################   
%##############################################################################################################################################################################################################      
%#################################################################################     AUXILARRY FUNCTIONS         ############################################################################################
%##############################################################################################################################################################################################################      
%##############################################################################################################################################################################################################   
%##############################################################################################################################################################################################################   

function [OUTPUT_PROPERTIES OUTPUT_PARAMETERS ERROR SHEET1 SHEET2 SHEET3 ANALYSIS_ALREADY_EXECUTED]= CHECK_FOR_ALREADY_EXECUTED_ANALYSIS(NI_COMPO,HT_TIME,HT_TEMP,VF,CASE_NUMBER,ISOB_NO)

% HERE I AM USING THE CSV FILES THAT I CREATE DURING THE EXECUTION OF THE BOED CODE IN ORDER TO SEE IF I HAVE ALREADY EXECUTED A SIMULATION WITH THE THE SAME INPUTV ARIABLES VF AND NI_COMPO.
% IF IN THIS CHECK I FAIL TO DETECT AN ALREADY EXECUTED SIMULATION THEN I AM ALSO CHECKING A SECONDARY CSV FILE WHERE I HAVE A DATBASE FROM OLDER RUNS FROM THE BOED CODE. THIS LIBRARY BASICALLY IS USEFULL FOR RESTARTING THE PROCES...
idx1=0;
idx2=0;
idx3= exist ('OUTPUT1.csv'); 
idx4= exist ('OUTPUT1_OLD.csv');
idx5=0;

if (idx3)~=0 % CHECK IF THERE IS ANALYSIS ALREADY EXECUTED,IN THE CURRENT OUTPUT1,2 AND 3 DATABASES, WITH THE SAME INPUT VARIABLES VF AND NI_COMPO
sprintf(['CHECK IF THERE IS ANALYSIS ALREADY EXECUTED WITH THE SAME INPUT VARIABLES VF AND NI_COMPO - CHECK IN THE CURRENT DATABASES'])
A1=csvread('OUTPUT1.csv'); % READ THE OUTPUT FILE THAT I SAVE AFTER EACH ANALYSIS THAT IS PERFORMED. I AM USING THIS BECAUSE ONLY THERE I CAN FIND THE VF AND NI_COMPO. NOTE THAT AT THIS POINT WE ARE ALREADY IN THE ROOT DIRECTORY.
A2=csvread('OUTPUT2.csv'); % READ THE OUTPUT FILE THAT I SAVE AFTER EACH ANALYSIS THAT IS PERFORMED. I AM USING THIS BECAUSE ONLY THERE I CAN FIND THE VF AND NI_COMPO. NOTE THAT AT THIS POINT WE ARE ALREADY IN THE ROOT DIRECTORY.
A3=csvread('OUTPUT3.csv'); % READ THE OUTPUT FILE THAT I SAVE AFTER EACH ANALYSIS THAT IS PERFORMED. I AM USING THIS BECAUSE ONLY THERE I CAN FIND THE VF AND NI_COMPO. NOTE THAT AT THIS POINT WE ARE ALREADY IN THE ROOT DIRECTORY.
n1=size(A1,1);
idx1 = 0;
for m = 1:n1
if A1(m, 4) == NI_COMPO && A1(m, 9) == VF
        idx1 = m;
end
end
end

if idx4~=0  % CHECK IF THERE IS ANALYSIS ALREADY EXECUTED WITH THE SAME INPUT VARIABLES VF AND NI_COMPO
sprintf(['CHECK IF THERE IS ANALYSIS ALREADY EXECUTED WITH THE SAME INPUT VARIABLES VF AND NI_COMPO - CHECK IN THE OLD DATABASES'])
B1=csvread('OUTPUT1_OLD.csv'); % READ THE OUTPUT FILE THAT I SAVE AFTER EACH ANALYSIS THAT IS PERFORMED. I AM USING THIS BECAUSE ONLY THERE I CAN FIND THE VF AND NI_COMPO. NOTE THAT AT THIS POINT WE ARE ALREADY IN THE ROOT DIRECTORY.
B2=csvread('OUTPUT2_OLD.csv'); % READ THE OUTPUT FILE THAT I SAVE AFTER EACH ANALYSIS THAT IS PERFORMED. I AM USING THIS BECAUSE ONLY THERE I CAN FIND THE VF AND NI_COMPO. NOTE THAT AT THIS POINT WE ARE ALREADY IN THE ROOT DIRECTORY.
B3=csvread('OUTPUT3_OLD.csv'); % READ THE OUTPUT FILE THAT I SAVE AFTER EACH ANALYSIS THAT IS PERFORMED. I AM USING THIS BECAUSE ONLY THERE I CAN FIND THE VF AND NI_COMPO. NOTE THAT AT THIS POINT WE ARE ALREADY IN THE ROOT DIRECTORY.
n1=size(B1,1);
idx2 = 0;
for m = 1:n1
if B1(m, 4) == NI_COMPO && B1(m, 9) == VF
        idx2 = m;
end
end
end


if(idx1~=0 && idx2~=0)
sprintf(['THERE IS MICROMECHANICAL ANALYSIS ALREADY EXECUTED, IN THE OLD AND CURRENT DATABASE, WITH THE SAME INPUT VARIABLES VF AND NI_COMPO']) 
ANALYSIS_ALREADY_EXECUTED=3;
C1=A1; % you can use also C1=B1;
C2=A2;
C3=A3;
idx5=idx1;
elseif(idx1~=0 && idx2==0)
sprintf(['THERE IS MICROMECHANICAL ANALYSIS ALREADY EXECUTED, IN THE CURRENT DATABASE, WITH THE SAME INPUT VARIABLES VF AND NI_COMPO']) 
ANALYSIS_ALREADY_EXECUTED=1;
C1=A1;
C2=A2;
C3=A3;
idx5=idx1;
elseif(idx1==0 && idx2~=0)
sprintf(['THERE IS MICROMECHANICAL ANALYSIS ALREADY EXECUTED, IN THE OLD DATABASE, WITH THE SAME INPUT VARIABLES VF AND NI_COMPO']) 
ANALYSIS_ALREADY_EXECUTED=2;
C1=B1;
C2=B2;
C3=B3;
idx5=idx2;
elseif(idx1==0 && idx2==0)
sprintf(['THERE IS NO MICROMECHANICAL ANALYSIS ALREADY EXECUTED, IN THE OLD AND CURRENT DATABASE, WITH THE SAME INPUT VARIABLES VF AND NI_COMPO']) 
ANALYSIS_ALREADY_EXECUTED=0;
idx5=0;
C1=0;
C2=0;
C3=0;
end


if ANALYSIS_ALREADY_EXECUTED~=0
sprintf([' LOAD THE RESULTS FROM THE ALREADY EXECUTED ANALYSIS AND PROCEED'])   
idx6=(idx5-1)*ISOB_NO+1;
idx7=(idx5-1)*ISOB_NO*3+1;

SHEET1=C1(idx5,:);
SHEET1(1,1)=CASE_NUMBER; % SAVE THE CURRENT CASE NUMBER
SHEET1(1,2)=HT_TEMP; % SAVE THE CURRENT HEAT TREATMENT TEMP
SHEET1(1,3)=HT_TIME; % SAVE THE CURRENT HEAT TREATMENT TIME

SHEET2=C2(idx6:(idx6-1)+ISOB_NO,:);
SHEET2(1:ISOB_NO,1)=CASE_NUMBER; % SAVE THE CURRENT CASE NUMBER
SHEET2(1:ISOB_NO,2)=HT_TEMP; % SAVE THE CURRENT HEAT TREATMENT TEMP
SHEET2(1:ISOB_NO,3)=HT_TIME; % SAVE THE CURRENT HEAT TREATMENT TIME

SHEET3=C3(:,idx7:(idx7-1)+3*ISOB_NO);
for i=1:ISOB_NO
SHEET3(1,i*3-1)=CASE_NUMBER;
end

root_directory=pwd;
OUTPUT_WRITE_UPDATE_CURRENT_DATABASE(root_directory,SHEET1,SHEET2,SHEET3)
OUTPUT_PROPERTIES=SHEET1;
OUTPUT_PARAMETERS=SHEET2;
% OUTPUT_VALUES=OUTPUT_PROPERTIES;
if OUTPUT_PROPERTIES(13)~=888 & OUTPUT_PROPERTIES(13)~=999
ERROR=0;
else
ERROR=1; % any value than zero will do
end
else
OUTPUT_PROPERTIES=0;
OUTPUT_PARAMETERS=0;
ERROR=0;
SHEET1=0;
SHEET2=0;
SHEET3=0;
end

end


function OUTPUT=OUTPUT_WRITE_UPDATE_CURRENT_DATABASE(root_directory,SHEET1,SHEET2,SHEET3)

CSV_NAME1=[root_directory,'/OUTPUT1.csv'];
CSV_NAME2=[root_directory,'/OUTPUT2.csv'];
CSV_NAME3=[root_directory,'/OUTPUT3.csv'];
OUTNAME1=[root_directory,'/OUTPUT_MATERIAL_PARAMETERS_1.txt'];
OUTNAME2=[root_directory,'/OUTPUT_MATERIAL_PARAMETERS_2.txt'];

if exist (CSV_NAME1)==0
INDEX1(1:3)=0;
INDEX2(1:3)=0;
INDEX3(1:3)=0; 
INDEX4(1:3)=0;
else
A1=csvread('OUTPUT1.csv');
A2=csvread('OUTPUT2.csv');
A3=csvread('OUTPUT3.csv');
INDEX1=size(A1);
INDEX2=size(A2);
INDEX3=size(A3);
end


INDEX4=size(SHEET1);
INDEX5=size(SHEET2);
INDEX6=size(SHEET3);
A1(INDEX1(1)+1,:)=SHEET1;
A2(INDEX2(1)+1:INDEX2(1)+INDEX5(1),:)=SHEET2;
A3(:,INDEX3(2)+1:INDEX3(2)+INDEX6(2))=SHEET3;
  

dlmwrite(CSV_NAME1, A1,'precision', 19)
dlmwrite(CSV_NAME2, A2,'precision', 19)
dlmwrite(CSV_NAME3, A3,'precision', 19)

fid = fopen(OUTNAME1,'a');
Res=fprintf(fid, '%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %e %f\n',SHEET1);
fclose(fid);     

for i=1:INDEX5
fid = fopen(OUTNAME2,'a');
Res=fprintf(fid, '%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n',SHEET2(i,:));
fclose(fid);       
end 


OUTPUT=1;


end

function OUTPUT=OUTPUT_WRITE_UPDATE_OLD_DATABASE(root_directory,SHEET1,SHEET2,SHEET3)

CSV_NAME1=[root_directory,'/OUTPUT1_OLD.csv'];
CSV_NAME2=[root_directory,'/OUTPUT2_OLD.csv'];
CSV_NAME3=[root_directory,'/OUTPUT3_OLD.csv'];

if exist (CSV_NAME1)==0
INDEX1(1:3)=0;
INDEX2(1:3)=0;
INDEX3(1:3)=0; 
INDEX4(1:3)=0;
else
A1=csvread('OUTPUT1_OLD.csv');
A2=csvread('OUTPUT2_OLD.csv');
A3=csvread('OUTPUT3_OLD.csv');
INDEX1=size(A1);
INDEX2=size(A2);
INDEX3=size(A3);
end


INDEX4=size(SHEET1);
INDEX5=size(SHEET2);
INDEX6=size(SHEET3);
A1(INDEX1(1)+1,:)=SHEET1;
A2(INDEX2(1)+1:INDEX2(1)+INDEX5(1),:)=SHEET2;
A3(:,INDEX3(2)+1:INDEX3(2)+INDEX6(2))=SHEET3;
  
dlmwrite(CSV_NAME1, A1,'precision', 19)
dlmwrite(CSV_NAME2, A2,'precision', 19)
dlmwrite(CSV_NAME3, A3,'precision', 19)

OUTPUT=1;


end


function [idx, CASE_NUMBER]=IDENTIFY_IDX_AND_CASE_NUMBER(INPUT_VALUES,XSPACE,nV,nO,nC)
% IDENTIFIES WHICH ROW OF XSPACE REPRESENTS THE INPUT VARIABLES VR1,VR2
% IT RETURNS THE CASE_NUMBER THAT CORRESPONDS TO THE IDENTIFIED ROW
% DONT CONFUSE IDX WITH THE CASE_NUMBER.
XSPACER=XSPACE(:,2:nV+1);

for i=1:nV
VR(i)=INPUT_VALUES(i);
end


idx = 0;
for m = 1:size(XSPACE, 1)
    
   
if XSPACER(m, :) == VR(:)' 
        idx = m;
        CASE_NUMBER=XSPACE(m, 1);
end
end
  

if idx == 0
error('something wrong')
end


end

function XSPACE=DXSPACE(VR1_SPACE,VR2_SPACE,VR3_SPACE)
delete('CASES.csv');

XSPACE1=combvec(VR2_SPACE,VR1_SPACE)'; % DEFINE XSPACE MATRIX BASED ON THE INPUT SPACE
XSPACE(:,1)=1:1:size(XSPACE1,1); % THIS COLUMN NOTES THE CASE NUMBERS.
XSPACE(:,2)=XSPACE1(:,1);
XSPACE(:,3)=XSPACE1(:,2);
XSPACE(:,2)=round(XSPACE(:,2),3);
XSPACE(:,3)=round(XSPACE(:,3),3);

dlmwrite('CASES.csv', XSPACE,'precision', 19)


end

function A=COMB(A1,A2)
% TAKES AS INPUT TWO MATRIXES A1 AND A2 AND RETURNS THE COMBINATIONS OF
% THEM

if size(A1,1)==1 % WHEN THE A1 IT IS VECTOR AND IT IS DEFIND AS A ROW I HAVE TO SWITCH IT TO COLUMN
A1=A1';
end


index1=max(size(A1));
index2=min(size(A1));
index3=max(size(A2));
index4=min(size(A2));


k=0;
for i =1:index1;
    for j=1:index3;
    start=(j-1)*index1+1;
    finish=(j)*index1;    
    A(start:finish,1:index2)=A1(:,1:index2);
    A(start:finish,index2+1)=A2(j);
    k=k+1;
    end
end


end

function [ XSPACE, dnum, qnum ] =SORT_XSPACE(XSPACE,dnum,qnum,INPUT_VALUES,ERROR,nV,nO,nC)

[idx, CASE_NUMBER]=IDENTIFY_IDX_AND_CASE_NUMBER(INPUT_VALUES,XSPACE,nV,nO,nC);

    if ERROR~=0
        XSPACE([end-dnum, idx], :) = XSPACE([idx, end-dnum], :);
        dnum = dnum+1;
    else
        qnum = qnum+1;
        XSPACE([qnum, idx], :) = XSPACE([idx, qnum], :);
    end
   
end


function PLOT_FIGURES(X,Y,XSPACE,PARETO_FRONT,qnum, dnum)

OBJECTIVE1_max=max(Y(:,1));
OBJECTIVE1_min=min(Y(:,1));
OBJECTIVE2_max=max(Y(:,2));
OBJECTIVE2_min=min(Y(:,2));
scalex=1.24;
scaley=1.15;


figure(1)
% PLOT FORMATING
x0=0;
y0=0;
width=500;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
set(gca,'fontsize',28)
set(gca,'fontWeight','bold')
set(gca,'fontName','Arial')
set(gca,'linewidth',2)
set(gca, 'box', 'on')
points=500;% THE NUMBER OF THE POINTS THAT I WILL DISCRITIZE THE CMAP
xlabel('|A_f-30| (^oC)','FontSize', 36)
ylabel('|A_f-M_s-40| (^oC)','FontSize', 36)
axis([OBJECTIVE1_min OBJECTIVE1_max*scalex OBJECTIVE2_min OBJECTIVE2_max*scaley])
daspect([OBJECTIVE1_max*scalex OBJECTIVE2_max*scaley 1])

% DEFINITIO OF THE VARIABLES
hold on
NI_COMPO_max=max(X(:,1));
NI_COMPO_min=min(X(:,1));
VF_max=max(X(:,2));
VF_min=min(X(:,2));


% STORE EVERYTHING IN THE OUT VARIABLES FOR DEBUGGING
clear out
for i=1:qnum
NI_COMPO=X(i,1);
VF=X(i,2);
color=lin(NI_COMPO,NI_COMPO_max,NI_COMPO_min,1,0);
rad=lin(VF,VF_max,VF_min,80,20);    
out(i,1)=Y(i,1);
out(i,2)=Y(i,2);
out(i,3)=VF;
out(i,4)=NI_COMPO;
out(i,5)=fin(color,points);
end

% DEFINITION OF THE COLOR MAP
clear cmap
cmap = jet(points);
%cmap=flip(cmap)
c=colormap(cmap);
caxis([NI_COMPO_min NI_COMPO_max]);
c=colorbar;
c.Label.String = 'Ni Concentration (% at.)';
c.Label.FontSize = 36;
c.Label.VerticalAlignment='top';

% PLOT ENTIRE OBJECTIVE SPACE
for i=1:qnum
NI_COMPO=X(i,1);
VF=X(i,2);
color=lin(NI_COMPO,NI_COMPO_max,NI_COMPO_min,1,0);
rad=lin(VF,VF_max,VF_min,90,30);    
hold on
plot(Y(i,1), Y(i,2),'MarkerSize', rad,'marker','.','Color',cmap(fin(color,points),:))% Show the objectives values for all the executed analysis
hold on
plot(Y(i,1), Y(i,2),'MarkerSize', rad/3.3,'marker','o','Color','black')% Show the objectives values for all the executed analysis
end


%PLOT PARETO FRONT
for i=1:size(PARETO_FRONT,1)
    
for j=1:qnum
    if PARETO_FRONT(i,1)==Y(j,1) && PARETO_FRONT(i,2)==Y(j,2)
    index=j;
    end
end
NI_COMPO=X(index,1);
VF=X(index,2);
color=lin(NI_COMPO,NI_COMPO_max,NI_COMPO_min,1,0);
rad=lin(VF,VF_max,VF_min,90,30);
hold on
plot(PARETO_FRONT(i,1), PARETO_FRONT(i,2)); %  Show Front Array
end
plot(PARETO_FRONT(:,1), PARETO_FRONT(:,2),'-.black','LineWidth', 2); %  Show Front Array

end

function value=lin(X_EVAL,X_MAX,X_MIN,Y_MAX,Y_MIN)
% DEFINE A LINEAR FUNCTION WHICH AT X_MAX -->Y_MAX AND X_MIN-->Y_MIN
% CALCULATE THE VALUE OF THE FUNCTION AT X_EVAL
   alpha=(Y_MAX-Y_MIN)/(X_MAX-X_MIN);
   bhta=Y_MIN-alpha*X_MIN;
   value=alpha*X_EVAL+bhta;

end

function output=fin(INPUT,L)
%DEFINES AN XX SPACE FROM 0 TO 1 WITH POINTS EQUAL TO THE ONES SPECIFIED BY L.
%THEN IT CHECKS THE DEFINED INPUT BETWEEN WHICH POINTS BELONGS IN THE XX
%SPACE AND IF FOR EXAMPLE IT BELONGS BETWEEN 4 AND 5 IT RETURNS THE 4 (THE FIRST OF THE TWO).
% THIS IS USEFULL FOR THE DEFINITION OF THE COLORS. THE CMAP IS DEFINED AS
% A COLOR MAP WITH L POINTS. THE 1ST POINT IS THE LOW THE LAST IS THE HIGH.
% HENCE FIRST I MAP THE COMPOSITION FROM 1 TO 0 AND THEN I DEFINE USING THE
% FIN THE COLOR THAT CORRESPONDS TO THAT COMPOSITION. FOR EXAMPLE IF
% COMPOSITION 50.3 CORRESPONDS TO COLOR 0.64. THEN I SEARCH IN THE XX TO
% FIND BETWEEN WHICH POINTS THE 0.64 LIES AND THEN I PICK THE INDEX OF THE
% FIRST OF THE TWO. THEN I DEFINE AS COLOR THE COLOR FROM THE CMAP THAT
% CORRESPONDS TO THAT INDEX.


xx=0:1/(L-1):1;

for i=1:L

if (INPUT>=xx(i))  
  output=i;  
end
    
end

end


