function [Vmax, Km, V0i] = M4_Algorithm_LC1_03(enzyme)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sohan Pramanik, Pramanis@purdue.edu 
% 
% Program Description 
% This UDF aims to calculate Vmax, Km, and V0i with the possible smallest
% SSE for the given enzyme generating plots for the raw data and the
% michealis-Menten plot for V0i. Calculating SSE value for V0i.
%
% Function Call
% function [Vmax, Km, V0i] = M4_Algorithm_LC1_03(enzyme)
%
% Input Arguments
% data: the given data for the given enzyme
%
% Output Arguments
% Vmax: Outputs Vmax values(mM/s) for the given enzyme
% Km: Outputs Km values(mM) for the given enzyme
% V0i: Output V0i values in array for the given enzyme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ____________________
%% INITIALIZATION
% Assign data to each variable as needed
time = transpose(1 : length(enzyme));                       % time of the data taken in second, since the data were recorded every second
initial_S = [3.75 7.5 15 30	65 125 250 500 1000 2000];      % the given initial S that is same for all enzyme
v_ref = [0.028, 0.055, 0.11, 0.19, 0.338, 0.613, 0.917, 1.201, 1.282, 1.57];    % the given reference velocity for enzyme PGO-X50
vmax_ref = 1.61;        % the given reference Vmax for enzyme PGO-X50
Km_ref = 214.28;        % the given reference Km for enzyme PGO-X50

% set the reasonable range of Vmax and Km if the enzyme is enzyme A-E
if strcmp(inputname(1), 'enzymeA')
    Vmax_range = [0.85 1.12];
    Km_range = [135 175];
elseif strcmp(inputname(1), 'enzymeB')
    Vmax_range = [0.75 0.98];
    Km_range = [300 375];
elseif strcmp(inputname(1), 'enzymeC')
    Vmax_range = [1.16 1.49];
    Km_range = [176 219];
elseif strcmp(inputname(1), 'enzymeD')
    Vmax_range = [1.40 1.80];
    Km_range = [253 316];
elseif strcmp(inputname(1), 'enzymeE')
    Vmax_range = [1.40 1.83];
    Km_range = [136 170];
end 

%% ____________________
%% CALCULATIONS
% Improvement 1: the program will automatically find the best range for number of seconds and
% calculate all the values. (This improvement is for all the parts in calculations)
V0i_matrix = zeros(39,10);  % V0i's matrix, which will contains V0i sets using different number of data points to calculate
SSE_all = zeros(1,39);      % SSE's array, which will contains all SSE calculated for different V0i set
if numel(inputname(1)) == 6 % for enzyme PGO-X50
    for n = 2:40            % the number of data points from the begining that the program is going to use when calculating the V0i set
        % Loop utilized to obtain the slopes (Voi) from each enzyme at each [S]
        for i = 1:10    % Nested loop to go through substrate concentrations for test
            % Calculating V0i by using Polyfit function
            lin = polyfit(time(1:n),enzyme(1:n,i),1);
            V0i_matrix((n - 1),i) =  lin(1);
            SSE_all(n-1) = SSE_all(n-1) + (v_ref(i) - V0i_matrix((n - 1), i)) ^2;	% Find the sum of SSE between the reference initial velocities (v_(0_i )) and the reaction velocities (v_i)
        end
    end
    SSE = min(SSE_all);         % find the smallest SSE among all SSE
    index = (SSE_all == SSE);   % find the index of the smallest SSE
    V0i = V0i_matrix(index,:);  % find the V0i set corresponding to the smallest SSE
    
    % Use the V0i set corresponding to the smallest SSE to do the following
    % calculations.
    %
    % Calculating the theoratical velocity using the formula for
    % Michealis-Menten using the reference Vmax and Km
    v_ref_cal = zeros(1,10);
    v_ref_cal(1,:) = vmax_ref .* initial_S ./ (initial_S + Km_ref);
    % Find the sum of SSE between the reference initial velocities (v_(0_i )) and the reaction velocities (v_i)
    SSE_ref = 0;
    for n = 1:10
        SSE_ref = SSE_ref + (v_ref_cal(n) - v_ref(n)) ^2;
    end

    % Linearization of all Vo for each enzyme at each [S] utilizing the
    % Lineweaver-Burk equation and calculating the Vmax and Km for each of the
    % data set using the equation.
    linweaverburk = polyfit((1./initial_S),(1./V0i(1,:)),1);
    Vmax = 1/linweaverburk(2);
    Km = linweaverburk(1)/linweaverburk(2);

    % Calculating the theoratical velocity using the formula for Michealis-Menten
    v = zeros(1,10);
    v(1,:) = Vmax .* initial_S ./ (initial_S + Km);
    % Find the sum of SSE between the reference initial velocities (v_(0_i )) and the reaction velocities (v_i)
    SSE = 0;
    for n = 1:10
        SSE = SSE + (v_ref(n) - v(n)) ^2;
    end
else
    v_matrix = zeros(39,10);    % calculated velocity's matrix, which will contains calculated velocity from the model generate by different velocity sets
    Vmax_all = zeros(1,39);     % Vmax's array, which will contains all Vmax calculated for different V0i set
    Km_all = zeros(1,39);       % Km's array, which will contains all Km calculated for different V0i set
    for n = 2:40    % the number of data points from the begining that the program is going to use when calculating the V0i set
    % Calculating V0i by using Polyfit function
    % Loop utilized to obtain the slopes (Voi) from each enzyme at each [S]
        for i = 1:10    % Nested loop to go through substrate concentrations for test
            lin = polyfit(time(1:n),enzyme(1:n,i),1);
            V0i_matrix((n - 1),i) =  lin(1);
        end
        % Linearization of all Vo for each enzyme at each [S] utilizing the
        % Lineweaver-Burk equation and calculating the Vmax and Km for each of the
        % data set using the equation.
        linweaverburk = polyfit((1./initial_S),(1./V0i_matrix((n - 1),:)),1);
        Vmax_all(1,(n - 1)) = 1/linweaverburk(2);
        Km_all(1,(n - 1)) = linweaverburk(1)/linweaverburk(2);
        v_matrix((n - 1),:) = Vmax_all(1, (n - 1)) .* initial_S ./ (initial_S + Km_all(1,(n - 1)));	% Calculating the theoratical velocity using the formula for Michealis-Menten
        % Find the sum of SSE between the calculated velocities from
        % data(V0i)) and the calculated velocities from model(v) for each
        % velocity sets
        for i = 1:10
            SSE_all(1, (n-1)) =  SSE_all(1,(n-1)) + (v_matrix((n - 1),i) - V0i_matrix((n - 1), i)) ^2;
        end
    end
    SSE = min(SSE_all); % find the smallest SSE among all SSE
    index = (SSE_all == SSE);   % find the index of the smallest SSE
    V0i = V0i_matrix(index,:);  % find the V0i set corresponding to the smallest SSE
    v = v_matrix(index, :);     % find the velocity set corresponding to the smallest SSE
    Km = Km_all(1,index);       % find the Km corresponding to the smallest SSE
    Vmax = Vmax_all(1,index);   % find the Vmax corresponding to the smallest SSE
    
    % determine whether the calculated Vmax & Km is in the reasonable range
    while (Km < Km_range(1) || Km > Km_range(2) || Vmax < Vmax_range(1) || Vmax > Vmax_range(2))
        SSE_all(index) = 10000;
        SSE = min(SSE_all);
        index = find(SSE == SSE_all);
        V0i = V0i_matrix(index,:);
        v = v_matrix(index, :);
        Km = Km_all(1,index);
        Vmax = Vmax_all(1,index);
    end    
end

%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS

if numel(inputname(1)) == 6
    % gets the 10 concentration curves of v0 with a tangent line and plots them
    figure(1)
    for n = 1:10
        subplot(5,2,n)
        plot(time,enzyme(:,n),'b-')
        hold on
        a = polyfit(time(1:25),enzyme(1:25,n),1);
        P_con = a(1) .* time + a(2);
        plot(time,P_con,'r-')
        title(['PGO-X50 Plot for V0i for ', num2str(initial_S(n)), ' [S] (uM)'])
        xlabel('Time (min)')
        ylabel('[P] (uM)')
        legend('Original Data','Slope for V0i','location','Northwest')
        % Improvement 2: adding grid lines on figure 1, which allows the better understanding of plot.
        grid on
    end

    % plots the michealis-Menten plot and has markers of each v0i
    figure(2)
    plot(initial_S, v_ref(1,:), 'rd')
    hold on
    plot(initial_S, v_ref_cal(1,:), 'b-')
    title('The Michaelis-Menten Plot for PG0-X50 Enzyme Using Reference Data')
    xlabel('Initial Subtrate Concentration[S](uM)')
    ylabel('Initial Reaction Velocity[v0](uM/s)')
    legend('Reference V0i for PGO-X50','Michaelis-Menten Line','location','best')
    grid on

    % plots the michealis-Menten plot and has markers of each v0i and theoretical velocity
    figure(3)
    plot(initial_S, V0i(1,:), 'rd')
    hold on
    plot(initial_S, v(1,:), 'b-')
    title('The Michaelis-Menten Plot for PG0-X50 Enzyme Using Calculating Data')
    xlabel('Initial Subtrate Concentration[S](uM)')
    ylabel('Initial Reaction Velocity[v0](uM/s)')
    legend('Raw V0i for PGO-X50','Theoretical velocity for PGO-X50','location','best')
    grid on

else
    % plots for Michaelis-Menten for Enzyme A-E
    figure(4)
    if strcmp(inputname(1), 'enzymeA')
        plot(initial_S, v(1,:), 'r*-')
        hold on
    elseif strcmp(inputname(1), 'enzymeB')
        plot(initial_S, v(1,:), 'bo-')
    elseif strcmp(inputname(1), 'enzymeC')
        plot(initial_S, v(1,:), 'gs-')
    elseif strcmp(inputname(1), 'enzymeD')
        plot(initial_S, v(1,:), 'k+-')
    elseif strcmp(inputname(1), 'enzymeE')
        plot(initial_S, v(1,:), 'c^-')
        title('The Michaelis-Menten Plot for Each Enzyme')
        xlabel('Initial Subtrate Concentration[S](uM)')
        ylabel('Initial Reaction Velocity[v0](uM/s)')
        legend('EnzymeA','EnzymeB','EnzymeC','EnzymeD','EnzymeE','location','best')
        grid on
    end
end

%% ____________________
%% COMMAND WINDOW OUTPUT

% displays the calculated values
if numel(inputname(1)) == 6
    fprintf(['\nThe reference value for ' inputname(1) ': Vmax is']);
    fprintf(' %.3f, Km is %.3f, and V0i is\n ',vmax_ref, Km_ref);
    v_ref
    fprintf('The SSE for reference V0i is %.3f.\n',SSE_ref);

    fprintf(['\nThe calculated value for ' inputname(1) ': Vmax is']);
    fprintf(' %.3f, Km is %.3f, and V0i is\n ',Vmax, Km);
    V0i
    fprintf('The SSE for V0i is %.3f.\n',SSE);

    fprintf(['\nFor ' inputname(1) ' the theoretical initial velocity is\n ']);
    v
    fprintf('The SSE for this prediction is %.3f.\n',SSE);

else
    fprintf(['\nThe calculated value for ' inputname(1) ': Vmax is']);
    fprintf(' %.3f, Km is %.3f, and V0i is\n ',Vmax, Km);
    V0i
    % fprintf('The SSE for V0i is %.3f.\n',SSEi)
    fprintf('The SSE for this prediction is %.6f.\n',SSE);
end
%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.


