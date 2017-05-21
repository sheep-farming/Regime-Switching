function[] = mvOpt()

% This function conducts a mean-variance optimization and plots 
% the mean-variance efficient frontier and the frontier allocations.
%
%
% Required Inputs are:
% - (T x N) matrix of returns, where T = number of observations, 
%    N = number of asset classes
% - (T x 1) vector with dates
% - (N x 1) vector of asset class names
%
%
% October 2014, Ph. Rohner


%% Input Parameters

close all

% set dates for estimation period
dateStart = datenum('31.01.1994','dd.mm.yyyy');
dateEnd = datenum('28.02.2017','dd.mm.yyyy');

method = 'constr'; % 'constr' for short sales constraints; 'unconstr' for unconstrained optimization
retFreq = 12; % DO NOT CHANGE! must be 12 for monthly return data, 52 for weekly return data
showTS = 'off'; % Display chart with asset class time series: 'on' or 'off'
showFrontier = 'on'; % Display efficient frontier: 'on' or 'off'
showAssetClasses = 'on'; % Display Asset Classes: 'on' or 'off'
showFrontierAllocations = 'on'; % Display frontier allocations: 'on' or 'off'
numPort = 30; % Number of portfolios on efficient frontier
tickSpace = 5; % number of dates on x-axis


%% load data file
load '../new.mat'
Returns = Returns(:,2:11);
[dummy numAC] = size(Returns);
datesVal = datenum(Dates);



%% Select Sample Period and calculate Prices
for j = 1:numAC
    dummy = Returns(:,j);
    Ret(:,j) = dummy((datesVal >= dateStart) & (datesVal <= dateEnd));
end
Dat = datesVal((datesVal >= dateStart) & (datesVal <= dateEnd));
[numVal numAC] = size(Ret);


Ind100(1,1:numAC) = 100;
for j = 2:numVal
    Ind100(j,:) = Ind100(j-1,:).*(1 + Ret(j,:));
end

switch showTS
    case 'on'
        plot(Ind100,'LineWidth',1.5);
        legend(Names,'Location','Northwest')
        ticks = round(1:((numVal - 1)/tickSpace):numVal);
        xTickNames = Dat(ticks);
        set(gca,'XTick',ticks)
        set(gca,'XTickLabel',datestr(xTickNames,12))
        axis([1,length(Ind100), 0, max(max(Ind100))])
    case 'off'
    otherwise
end

%% Calculations

ExpRet = mean(Ret)'; % use median to cope with outliers (robust statistic)
%ExpRet = expReturns;
CovMat = cov(Ret);

[frontWts, frontRet, frontVol] = MeanVarianceOptimization(ExpRet, CovMat, numPort, method);


%% Plot Efficient Frontier

switch showFrontier
    case 'on'
        portRet = ((1 + frontRet).^(retFreq) - 1);
        portVol = frontVol*sqrt(retFreq);
        ACvols = sqrt(diag(CovMat)*retFreq);

        figure(2)
        hold on
        plot(portVol, portRet,'LineWidth',2)
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
    case 'off'
end


%% Plot Asset Classes    
AssetReturns = ((1 + ExpRet).^(12) - 1);
AssetVols = diag(sqrt(CovMat))*sqrt(12);

switch showAssetClasses
    case 'on'
        for h = 1:length(ExpRet)
                plot(AssetVols(h), AssetReturns(h),'d','markersize',2,'MarkerEdgeColor','black','MarkerFaceColor','black','linewidth',3,'linestyle','none');
                if xmin+(xmax-xmin)/2 > AssetVols(h)*0.01
                    halign = 'left';
                    text(AssetVols(h) + 0.005, AssetReturns(h), ['   ' Names(h) '   '], 'horizontalalignment',halign, 'FontSize', 12);
                else
                    halign = 'right';
                    text(AssetVols(h) + 0.005, AssetReturns(h), ['   ' Names(h) '   '], 'horizontalalignment',halign, 'FontSize', 12);
                end
        end
    case 'off'
end

% Plot GMV Portfolio
figure(2)
hold on
plot(portVol(1,1),portRet(1,1),'o','markersize',2,'MarkerEdgeColor','red','MarkerFaceColor','black','linewidth',3,'linestyle','none')
text(portVol(1,1) + 0.005,portRet(1,1),'GMV','horizontalalignment','left', 'FontSize', 12)



%% plot Frontier Allocations

switch showFrontierAllocations
    case 'on'
        figure(3)
        area(frontWts')
        legend(Names, 'Location', 'NorthEastOutside') 
        title('Frontier Allocations','FontSize',12)
        ylabel('Portfolio fraction','FontSize',12)
        xlabel('Frontier portfolios','FontSize',12)
        axis([1, numPort, 0, 1]);
        box on
    case 'off'
end





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
    A1(1,1) = 1; % exclude asset class 1 (cash)
    A = [A1; A2];
    b = [0;  0];
    
    
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
    frontWts = zeros(numAssets, numPort)
    
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
    
        
        





    

