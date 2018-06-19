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
data = xlsread('../data/dose_response_example.xls');
% you will have to change this part depending on how your data is set up in
% excel
dose = data(:,1);
viability = data(:,2);

% make a vector the length of your data of the mean viability at dose =0
ind = find(dose == 0);
Vmax = viability(ind); % finds the viability only where dose =0
nreps = length(Vmax); % finds the number of replicates of dose response by 
                       % (by counting how many times dose =0)
ndose=12; % tell it the number of doses you used 
Vmaxmean= mean(Vmax); % find mean of the viabilities
% note that when doing multiple conditions, we will want to make sure that
% the Vmax changes for each condition, which is why we make the vector
% below
Vmaxall = repmat(Vmaxmean, ndose.*nreps,1);  % makes vector of Vmaxmean at each dose
% for different conditions, the value of Vmaxall will change to be the one
% corresponding to that condition.


%% Plot raw data
figure;
plot(dose, viability, 'o', 'LineWidth',2)
xlabel('dose (\muM)')
ylabel('Viability')
title('Dose Response Raw Data')


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


[P_fit, resnorm, residuals] = lsqnonlin(@fitsinglepop,...
    params0,...
    LB,...
    UB,...
    options,...
    dose,...
    viability,...
    Vmaxall);
% LSQNONLIN takes the function that outputs the error vector and gives it
% the initial guess, bounds, options, and other inputs of fit singlepop(in
% this case it is dose, viability (from data) and the Vmaxall vector, but
% their could be others)

% best fit parameters
m_fit= P_fit(1); 
LD50_fit = P_fit(2);

% Plug parameters back into model
Vmodel = Vmaxall./(1+exp(m_fit.*(dose-LD50_fit)));
% for plotting (don't make doses repeat)
dmod = 0:20:max(dose);
V = Vmaxmean./(1+exp(m_fit.*(dmod-LD50_fit)));

figure;
plot(dose, viability, 'ro', 'LineWidth',2)
hold on
plot(dmod, V, 'b-', 'LineWidth',2)
xlabel('dose (\muM)')
ylabel('Viability')
title('Dose Response Data and Model Fit')

%% How to find error bars using boot strapping method of replacement

% give the function the residuals from lsqnonlin, the model fit to the
% doses we have data for, the number of reps and the number of doses for
% making Vmaxall vector for each simulated dataset
[bootCI, betaboot] = finderrorBS(residuals,Vmodel, dose, nreps,ndose);
% this gives you the 95th percentile ranges of the parameter estimates for
% the slope and LD50
% ( 1st column is lower and upper bound of the slope, second column is
% % lower and upper bound of the LD50) 


