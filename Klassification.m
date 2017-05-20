clear;
load('bbgdata.mat');

%% Variables
%
% P Q <- Switching Probability
% uu1 uu2 <- World Excess Return
% th1 th2 <- World Volatility
%

%% Assumptions 
%
% uuz = uRf
uuz = mean(Returns(:,1))/1200;

wRet = Returns(:,5)/100; % Annualised world return per month
wExRet = wRet - Returns(:,1)/1200;

%% Output Here

addpath('MS_Regress-Matlab')
addpath('MS_Regress-Matlab/m_Files'); % add 'm_Files' folder to the search path

constVec=ones(length(wExRet),1); % A constant vector in mean equation (just an example of how to do it)
indep = [constVec];
k=2;                                % Number of States
S=[1 1];                        % Defining which parts of the equation will switch states (column 1 and variance only)
advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details

[Spec_Out]=MS_Regress_Fit(wExRet,indep,k,S,advOpt) % Estimating the model

%% Plot Return & Probabilities
val = [1];
subplot(2,1,2);
plot(Spec_Out.filtProb(:,1));
subplot(2,1,1);
hold on
for i=1:length(wExRet)
    val = [val val(end)*exp(wExRet(i)/12)];
    if(Spec_Out.filtProb(i,1)<.5)
        plot(i,val(i),'r+');
    else
        plot(i,val(i),'g+');
    end
    
end

P = Spec_Out.Coeff.p(1,1);
Q = Spec_Out.Coeff.p(2,2);

%% Mu and Sigma for Each Regimes
uus = Spec_Out.Coeff.S_Param{1};
sis = [sqrt(Spec_Out.Coeff.covMat{1}) sqrt(Spec_Out.Coeff.covMat{2})];

%% Betas of Assets
[betas, ses] = getBetas;

%% Asset Expected Returns
% Regime 1
ER = uuz + (uus'-uuz)*betas;

