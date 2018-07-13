% Dose Response Example
% This code imports a single dose response curve (dose versus viability)
% and fits it to the sigmoidal dose response curve given by:
% V(dose) = Vmax/(1 + exp(m(dose-LD50)))
% where m is the slope parameter (1/m is proportional to the standard
% deviation of the normal distribution of viability as a function of dose)
% and LD50 is the dose at which half of the cells die (center of the normal
% distribution)
% Note: We extract the maximum cell viability at the 0 dose to normalize for
% naturally occuring cell death.

% This code calls the scripts fitsinglepop (which puts the parameters into
% the model equation, and outputs an error vector of the difference between
% the value at the model with the current parameters and the data at that
% dose. This is used by lsqnonlin for the optimization

% The find errorBS code is called to find the bootstrapped parameter
% estimate intervals. In this function, we input the model predicted values
% at each dose, and then we generate new simulated data sets by randomly
% selecting from the residuals at each dose and adding those to the model
% values at each dose. Each simulated data set is fit to the model,
% generating 500 parameter sets. We take the 95th percentile of this set to
% find the 95 % confidence intervals.


close all; clear all; clc

% load in dataset (this one is just dose in 1st column, viability in 2nd
% column
data = xlsread('../data/dose_response_stiff.xls');
% you will have to change this part depending on how your data is set up in
% excel
dose = data(:,1);
viability = data(:,2);
stiffness = data(:,3);
uniqstiff = unique(stiffness);
Vmaxall = [];
% make a vector the length of your data of the mean viability at dose =0
for i = 1:length(uniqstiff)
    ind = find(stiffness==uniqstiff(i));
    i0 = find(dose(ind) == 0);
    V= viability(ind);
    Vmax = V(i0); % finds the viability only where dose =0
    n = length(V);
    Vmaxmean(i)= mean(Vmax); % find mean of the viabilities
% note that when doing multiple conditions, we will want to make sure that
% the Vmax changes for each condition, which is why we make the vector
% below
    Vmaxvec = repmat(Vmaxmean(i), n,1);  % makes vector of Vmaxmean at each dose
% for different conditions, the value of Vmaxall will change to be the one
% corresponding to that condition.
Vmaxall = vertcat(Vmaxall, Vmaxvec);
end


%% Plot raw data
figure;
for i = 1:length(uniqstiff)
    ind = stiffness== uniqstiff(i);
    plot(dose(ind), viability(ind), 'o', 'LineWidth',2)
    hold on
    %legend(['Stiffness= ', num2str(uniqstiff(i))])
end
xlabel('dose (\muM)')
ylabel('Viability')
title('Dose Response Raw Data')
legend(['Stiffness= ', num2str(uniqstiff(1)),' Pa'], ['Stiffness= ', num2str(uniqstiff(2)),' Pa'], ['Stiffness= ', num2str(uniqstiff(3)),' Pa'])
legend boxoff



%% Use lsqnonlin to perform fitting
% set conditions of lsqnonlin 
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
% two parameters to fit-- parameter 1 is slope, parameter 2 is LD50
% take initial guess and set bounds for parameter estimates
params0 = [0.01 180];
LB = [0 0];
UB = [ 1 Inf];
uniqstiff = unique(stiffness);

for i = 1:length(uniqstiff)
    ind = stiffness== uniqstiff(i);
[P_fit, resnorm, residuals] = lsqnonlin(@fitsinglepop,...
    params0,...
    LB,...
    UB,...
    options,...
    dose(ind),...
    viability(ind),...
    Vmaxall(ind));
% LSQNONLIN takes the function that outputs the error vector and gives it
% the initial guess, bounds, options, and other inputs of fit singlepop(in
% this case it is dose, viability (from data) and the Vmaxall vector, but
% their could be others)

% best fit parameters
m_fit(i)= P_fit(1); 
LD50_fit(i) = P_fit(2);

% Plug parameters back into model
Vmodel = Vmaxall(ind)./(1+exp(m_fit(i).*(dose(ind)-LD50_fit(i))));
% for plotting (don't make doses repeat)
dmod = 0:20:max(dose);
V = Vmaxmean(i)./(1+exp(m_fit(i).*(dmod-LD50_fit(i))));
end
colors = {'b', 'g', 'r'};
figure;
for i = 1:length(uniqstiff)
ind = stiffness== uniqstiff(i);
% Plug parameters back into model
Vmodel = Vmaxall(ind)./(1+exp(m_fit(i).*(dose(ind)-LD50_fit(i))));
% for plotting (don't make doses repeat)
dmod = 0:20:max(dose);
V = Vmaxmean(i)./(1+exp(m_fit(i).*(dmod-LD50_fit(i))));
plot(dose(ind), viability(ind), 'o','color', colors{i}, 'LineWidth',2)
hold on
plot(dmod, V, '-','color', colors{i}, 'LineWidth',2)
xlabel('dose (\muM)')
ylabel('viability')
title('Dose response data and sigmoid fit')
end
legend ('data 5 mM', 'fit 5 mM , LD50 = 38.7 uM', 'data 20 mM', 'fit 20 mM LD50 = 45.4 uM', 'data 25 mM', 'fit 25 mM LD50 = 33.1 uM ')
legend boxoff
%% How to find error bars using boot strapping method of replacement

% give the function the residuals from lsqnonlin, the model fit to the
% doses we have data for, the number of reps and the number of doses for
% making Vmaxall vector for each simulated dataset
[bootCI, betaboot] = finderrorBS(residuals,Vmodel, dose, nreps,ndose);
% this gives you the 95th percentile ranges of the parameter estimates for
% the slope and LD50
% ( 1st column is lower and upper bound of the slope, second column is
% % lower and upper bound of the LD50) 


