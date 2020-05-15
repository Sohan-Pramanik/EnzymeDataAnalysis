%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sohan Pramanik, Pramanis@purdue.edu
%
% Program Description 
% This is the main function that runs function M4_Algorithm_LC1_03 and
% M4_Regression_LC1_03.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ____________________
%% INITIALIZATION

% importing data from csv file
data = readmatrix('Data_PGOX50_enzyme.csv');
data1 = readmatrix('Data_nextGen_KEtesting_allresults.csv');

% Get valid data from two csv files
PGOX50 = data(5:end, 2:end);
enzyme_data = data1(3:end, 2:end);

% Call function for PGO-X50
[Vmax_PGOX50, Km_PGOX50, V0i_PGOX50] = M4_Algorithm_LC1_03(PGOX50);

%Generate enzyme data  for enzyme A by combining test data and duplicate
%test data
enzymeA = zeros (length(enzyme_data),10);
for i = 11:20    %Nested loop to go through substrate concentrations for test and duplicate test
    enzymeA(:,(i-10)) = (enzyme_data(:,i) + enzyme_data(:,(i-10)))/2;
end
% Call function for enzyme A
[Vmax_A, Km_A, V0i_A] = M4_Algorithm_LC1_03(enzymeA);
% edit data to fit next enzyme
enzyme_data = enzyme_data(:,21:end);

%Generate enzyme data  for enzyme B by combining test data and duplicate
%test data
enzymeB = zeros (length(enzyme_data),10);
for i = 11:20    %Nested loop to go through substrate concentrations for test and duplicate test
    enzymeB(:,(i-10)) = (enzyme_data(:,i) + enzyme_data(:,(i-10)))/2;
end
% Call function for enzyme B
[Vmax_B, Km_B, V0i_B] = M4_Algorithm_LC1_03(enzymeB);
% edit data to fit next enzyme
enzyme_data = enzyme_data(:,21:end);

%Generate enzyme data  for enzyme C by combining test data and duplicate
%test data
enzymeC = zeros (length(enzyme_data),10);
for i = 11:20    %Nested loop to go through substrate concentrations for test and duplicate test
    enzymeC(:,(i-10)) = (enzyme_data(:,i) + enzyme_data(:,(i-10)))/2;
end
% Call function for enzyme C
[Vmax_C, Km_C, V0i_C] = M4_Algorithm_LC1_03(enzymeC);
% edit data to fit next enzyme
enzyme_data = enzyme_data(:,21:end);

%Generate enzyme data  for enzyme D by combining test data and duplicate
%test data
enzymeD = zeros (length(enzyme_data),10);
for i = 11:20    %Nested loop to go through substrate concentrations for test and duplicate test
    enzymeD(:,(i-10)) = (enzyme_data(:,i) + enzyme_data(:,(i-10)))/2;
end
% Call function for enzyme D
[Vmax_D, Km_D, V0i_D] = M4_Algorithm_LC1_03(enzymeD);
% edit data to fit next enzyme
enzyme_data = enzyme_data(:,21:end);

%Generate enzyme data  for enzyme E by combining test data and duplicate
%test data
enzymeE = zeros (length(enzyme_data),10);
for i = 11:20    %Nested loop to go through substrate concentrations for test and duplicate test
    enzymeE(:,(i-10)) = (enzyme_data(:,i) + enzyme_data(:,(i-10)))/2;
end
% Call function for enzyme E
[Vmax_E, Km_E, V0i_E] = M4_Algorithm_LC1_03(enzymeE);

% Generate array of Km containing Km for all 
Km = [Km_A Km_B Km_C Km_D Km_E];
% Call function to calculate Price for enzyme A-E using Km
Price = M4_Regression_LC1_03(Km);

%% ____________________
%% COMMAND WINDOW OUTPUT
%{
M4_Main_function

The reference value for PGOX50: Vmax is 1.610, Km is 214.280, and V0i is
 
v_ref =

    0.0280    0.0550    0.1100    0.1900    0.3380    0.6130    0.9170    1.2010    1.2820    1.5700

The SSE for reference V0i is 0.025.

The calculated value for PGOX50: Vmax is 1.614, Km is 218.771, and V0i is
 
V0i =

    0.0272    0.0533    0.1088    0.1837    0.3376    0.6132    0.9174    1.1757    1.2820    1.5704

The SSE for V0i is 0.026.

For PGOX50 the theoretical initial velocity is
 
v =

    0.0272    0.0535    0.1036    0.1947    0.3698    0.5870    0.8609    1.1229    1.3245    1.4551

The SSE for this prediction is 0.026.

The calculated value for enzymeA: Vmax is 0.997, Km is 155.781, and V0i is
 
V0i =

    0.0234    0.0458    0.0875    0.1685    0.2851    0.4413    0.5979    0.7325    0.8512    0.9107

The SSE for this prediction is 0.001482.

The calculated value for enzymeB: Vmax is 0.908, Km is 354.162, and V0i is
 
V0i =

    0.0095    0.0192    0.0371    0.0702    0.1385    0.2245    0.3639    0.5148    0.6885    0.7464

The SSE for this prediction is 0.001519.

The calculated value for enzymeC: Vmax is 1.225, Km is 188.260, and V0i is
 
V0i =

    0.0240    0.0469    0.0880    0.1730    0.3149    0.4992    0.7230    0.8808    1.0399    1.1265

The SSE for this prediction is 0.000913.

The calculated value for enzymeD: Vmax is 1.619, Km is 294.458, and V0i is
 
V0i =

    0.0204    0.0398    0.0768    0.1557    0.2967    0.4849    0.7551    1.0365    1.2773    1.4294

The SSE for this prediction is 0.001557.

The calculated value for enzymeE: Vmax is 1.661, Km is 169.272, and V0i is
 
V0i =

    0.0360    0.0700    0.1411    0.2335    0.4615    0.7125    1.0020    1.2464    1.4402    1.5663

The SSE for this prediction is 0.002147.

Regression Data:
The calculated linearized model is: y = -2.4042 * x + 7.8954
The calculated SSE for linearized data is: 0.1397
The calculated SST for linearized data is: 4.4451
The calculated r^2 for linearized data is: 0.9686
The calculated general model for these data is: y = 78594899.5504 * x ^ (-2.4042)
The price prediction for enzyme:
Enzyme A = $420.7833
Enzyme B = $58.4115
Enzyme C = $266.8872
Enzyme D = $91.0470
Enzyme E = $344.6161
%}
%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.



