function[] = mvOpt()
% This function conducts a mean-variance optimization and plots 
% the efficient frontier
%
%
% Required Inputs are:
% - (T x N) matrix of returns, where T = number of observations, 
%    N = number of asset classes
%
%
% February 2014, Ph. Rohner


%% input parameters

method = 'constr'; % 'constr' for short sales constraints; 'unconstr' for unconstrained optimization
numPort = 100; % Number of portfolios on efficient frontier
retFreq = 12; % Return frequency 12=monthly


%% load data file

close all

load o.mat
numAC = 10;



%% calculate return and CovMat

ExpRet = o.ER1; 
CovMat = o.Sig1;


%% optimization

[frontWts, frontRet, frontVol] = MeanVarianceOptimization(ExpRet, CovMat, numPort, method);

ExpRet = o.ER2; 
CovMat = o.Sig2;

[frontWts2, frontRet2, frontVol2] = MeanVarianceOptimization(ExpRet, CovMat, numPort, method);

load ../new.mat
Returns = Returns(:,2:11);
ExpRet = mean(Returns)';
CovMat = cov(Returns);
[frontWts3, frontRet3, frontVol3] = MeanVarianceOptimization(ExpRet, CovMat, numPort, method);
%% outputs

    % plot efficient frontier

        portRet = ((1 + frontRet).^(retFreq) - 1);
        portVol = frontVol*sqrt(retFreq);
        
        portRet2 = ((1 + frontRet2).^(retFreq) - 1);
        portVol2 = frontVol2*sqrt(retFreq);
        
        
        portRet3 = ((1 + frontRet3).^(retFreq) - 1);
        portVol3 = frontVol3*sqrt(retFreq);
        
        ACvols = sqrt(diag(CovMat)*retFreq);

        figure(1)
        hold on
        plot(portVol, portRet,'LineWidth',2)
        plot(portVol2, portRet2,'LineWidth',2)
        plot(portVol3, portRet3,'LineWidth',2)
        xlabel('Expected Risk (Volatility)','FontSize',12)
        ylabel('Expected Return','FontSize',12)
        title('Mean-Variance Efficient Frontier' ,'FontSize',12)
        grid on
        box on
        xmax = max(ACvols) + 0.02;
        ymax = max(portRet) + 0.02;
        xmin = 0;
        ymin = 0.0;
        axis([xmin,xmax, ymin, ymax]);



function[frontWts, frontRet, frontVol] = MeanVarianceOptimization(ExpRet, CovMat, numPort, method)

    numAssets = length(ExpRet);
    V0 = zeros(1, numAssets);
    V1 = ones(1, numAssets);

    % Set constraints
    
    UB = ones(numAssets,1);
    
    switch method
        case 'unconstr'
            LB = ones(numAssets,1)*(-10);
            
        case 'constr'
            LB = zeros(numAssets,1);
    end
    
    % Define inequality constraints
    
    A1 = zeros(1,numAssets); 
    A2 = zeros(1,numAssets); 
    A = [A1; A2];
    b = [1; 1];
    
    
    % Find minimum variance return
    minVarWts = quadprog(CovMat, V0, A, b, V1, 1, LB);
    minVarRet = minVarWts'*ExpRet;
    minVarVol = sqrt(minVarWts'*CovMat*minVarWts);
    
    % Find maximum return
    maxRetWts = quadprog([], -ExpRet, A, b, V1, 1, LB);
    maxRet = maxRetWts'*ExpRet;
    
    
    % Calculate frontier portfolios
    targetRet = linspace(minVarRet, maxRet, numPort);
    frontRet = zeros(numPort,1);
    frontVol = zeros(numPort,1);
    frontWts = zeros(numAssets, numPort);
    
    frontRet(1,1) = minVarRet;
    frontVol(1,1) = minVarVol;
    frontWts(:,1) = minVarWts;
    
    Aeq = [V1; ExpRet'];

    
    for i = 2:numPort
        beq = [1; targetRet(i)];
        weights = quadprog(CovMat, V0, A, b, Aeq, beq, LB, UB);
        frontRet(i,1) = weights'*ExpRet;
        frontVol(i,1) = sqrt(weights'*CovMat*weights);
        frontWts(:,i) = weights;
    end
    
        
        





    

