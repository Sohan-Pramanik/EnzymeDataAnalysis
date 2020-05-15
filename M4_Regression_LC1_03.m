function [Price] = M4_Regression_LC1_03(Km)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sohan Pramanik, Pramanis@purdue.edu
%
% Program Description 
% This UDF generates the general model of the price for enzyme, and returns
% the price of enzyme giving the Km of the enzyme.
%
% Function Call
% function [Price] = M4_Regression_LC1_03(Km)
%
% Input Arguments
% Km: the Km value of the enzymes that need to calculate the price
%
% Output Arguments
% Price: the price of each enzyme given the Km value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ____________________
%% INITIALIZATION

% Assign data to each variable as needed
data = readmatrix('Data_NaturalCatalysts_priceCatalog.csv');
Michaelis_Constant = data(:,1);
Price_data = data(:,2);

%% ____________________
%% CALCULATIONS
% We have determined the data is exponential using several methods
logPrice = log10(Price_data);   % linearize the price data
logConstant = log10(Michaelis_Constant); % linearize Michaelis Constant
coeffs = polyfit(logConstant,logPrice,1);   % use polyfit function to find the parameter for linearized model
yVal = polyval(coeffs, logConstant); % the prices
M = coeffs(1);     
B = coeffs(2);
% determine the parameter for general model
m = M;          
b = 10 ^ (B);

logP_predict = M * Michaelis_Constant + B;      % calculated the log Price value using the linearized model
SSE_lin = sum((yVal - logPrice) .^2);   % calculate the SSE for linearized model
SST_lin = sum((logPrice - mean(logPrice)) .^2); % calculate the SST for linearized model
r_square = 1 - SSE_lin/SST_lin;                 % calculate the r^2 for linearized model
Price = b .* (Km .^ m);                         % calculate the price using the general model and the given Km for each enzyme
M_Constant = linspace(150, 350, length(Michaelis_Constant)); % spacing values for Km
Price_general_model = b .* (M_Constant .^ m); % calculate the price using the general model and all Km in the data set

%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS

% plot of the linearized data 
figure(10)
plot(Michaelis_Constant,Price_data,'b.')
hold on
plot(M_Constant,Price_general_model,'r-')
xlabel('Michaelis Constant (uM)')
ylabel('Price (USD($)/lb)')
title('Plot of raw data with general model overlaid ')
grid on 
legend('general model','raw data','location','Northeast')
equation = sprintf('y = %.4f * x^{%.4f}', b, m);
xt = (260);
yt = (150);
text(xt, yt, equation, 'Color', 'r', 'FontWeight', 'bold');

%% ____________________
%% COMMAND WINDOW OUTPUT

% displays the calculated values
fprintf('\nRegression Data:\n');
fprintf('The calculated linearized model is: y = %.4f * x + %.4f\n', M, B);
fprintf('The calculated SSE for linearized data is: %.4f\n', SSE_lin);
fprintf('The calculated SST for linearized data is: %.4f\n', SST_lin);
fprintf('The calculated r^2 for linearized data is: %.4f\n', r_square);
fprintf('The calculated general model for these data is: y = %.4f * x ^ (%.4f)\n', b, m);
fprintf('The price prediction for enzyme:\n');
fprintf('Enzyme A = $%.4f\n', Price(1));
fprintf('Enzyme B = $%.4f\n', Price(2));
fprintf('Enzyme C = $%.4f\n', Price(3));
fprintf('Enzyme D = $%.4f\n', Price(4));
fprintf('Enzyme E = $%.4f\n', Price(5));

%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.



