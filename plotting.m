%% Sharpe & Max Drawdown for Non-Regime-Dependent Strategy
prets = [];

for i=50:length(Returns)-1
    Returnz = Returns(1:50,2:11);
    er = mean(Returnz);
    cv = cov(Returnz);
    w = inv(cv)*er';
    w=w+abs(w); % Short-sell constraint
    w=w/sum(w);
    
    pret = Returns(i,2:11)*w;
    prets = [prets pret];
    
end
%input('Press Enter to Calculate Non-RS Strategy: Monthly Rebalancing with Short Sell Constraint')

%plot(ret2price(prets),'rx');

NonRegimeDependentMean = mean(prets)*12
NonRegimeDependentVolatility = std(prets)*sqrt(12)
NRDSharpe = (NonRegimeDependentMean-Returns(end,1)*12) / NonRegimeDependentVolatility
NRDMaxDrawDown = maxdrawdown(ret2price(prets))

%% PLOTTING: Mean Variance Portfolio & RS Strategy

clf 
plot(ret2price(rets),'k-','LineWidth',3)
hold on

plot(ret2price(Returns(51:end,12)),'g-','LineWidth',2);
plot(ret2price(prets),'k:','LineWidth',2);
plot(ret2price(Returns(51:end,6)),'r-','LineWidth',2);
plot(ret2price(Returns(51:end,4)),'b-','LineWidth',2);



Legends = {'RS Strategy',  'World','Non-Regime-Dependent Strategy','Swiss SPI','US S&P500'};
lgd = legend(Legends,'Location','Northwest','FontSize',20);
set(gca, 'Xtick', [0:12:240])
set(gca, 'XtickLabel', [1998:1:2018],'FontSize',12)


input('Press Enter')

%% PLOTTING: Probabilities


clf
hold on
plot(out.Spec_Out{end}.smoothProb(50:end,1),'k-');
plot(out.Spec_Out{end}.filtProb(50:end,1),'b--');
set(gca, 'Xtick', [0:12:240])
set(gca, 'XtickLabel', [1998:1:2018],'FontSize',12)


legend({'Smooth Probability','Filter Probability'},'Location','Southeast','FontSize',20)


input('Press Enter')


%% PLOTTING: Strategy & WORLD
clf
pld = [pld(:,2:end)];

plot(pld,'-','LineWidth',1);
hold on
%plot(ret2price(rets),'k-','LineWidth',3)

hold on

set(gca, 'Xtick', [0:12:240])
set(gca, 'XtickLabel', [1998:1:2018],'FontSize',12)
hold on
%Set Y ticks
a=[cellstr(num2str(get(gca,'ytick')'*100))];
pct = char(ones(size(a,1),1)*'%');
new_yticks = [char(a),pct];
%set(gca,'yticklabel',new_yticks,'FontSize',14)
ylabel('Cumulative Return','FontSize',24)

hold on
plot(ret2price(rets),'k-','LineWidth',3);
Legends = {'China','Japan','US','EMBI', 'Switzerland','North America', 'EU', 'UK', 'Pacific', 'EM','World','RS Strategy'};
lgd = legend(Legends,'Location','Northwest','FontSize',14);





