function[bootCI betaboot] = finderrorBS(residuals,model, dose, nreps, ndose)



nboot = 500; % number of bootstrapping with replacement data sets
[~, bootIndices] = bootstrp(nboot, [], residuals); % randomly generates indices
bootResiduals = residuals(bootIndices); % uses indices to sample from residuals with replacement
viaBoot = repmat(model,1,nboot) + bootResiduals; % creates simulated data sets with randomly added residuals
% viaBoot stands for bootstrap simulated viability
% Code below puts the bounds on the viability to make it between 0 and 1
for j = 1: size(viaBoot,1)
    for k = 1:size(viaBoot,2)
        if viaBoot(j,k) > 1 
            viaBoot(j,k) = 1;
        end
         if viaBoot(j,k) < 0
             viaBoot(j,k) = 0;
         end
             
    end
end
% build up the bootstrap data set
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
%psingle0 = horzcat( [180 0.1 35 0.1], 0.5.*ones(1, nsamp));%  LD50res, sloperes, LD50sens, slopesens, fres
paramslb = [ 0 0];
paramsub = [ 1 Inf];
psingle0 = [ 0.1 70]; % initial guess 
 for i = 1:nboot
     viasim = viaBoot(:,i);
     % Find the average viability for the simulated data set
     ind = find(dose == 0);
     Vmax = viasim(ind);
     Vmaxallboot =[];
     Vmaxmean= mean(Vmax);
     Vmaxallboot = repmat(Vmaxmean, ndose.*nreps,1); 
        betaboot(:,i) = lsqnonlin(@fitsinglepop,...
        psingle0,...
        paramslb,...
        paramsub,...
        options,...
        dose,...
        viasim,...
        Vmaxallboot);
end

 bootCI = prctile(betaboot', [2.5 97.5]);
end