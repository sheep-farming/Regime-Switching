load('neu.mat')
load('../new.mat')
out.uuz = Returns(:,1);
%% Eq 2

P = out.P;
Q = out.Q;

%% Eq 3, 4
uu1 = out.uu1;
uu2 = out.uu2;
uuz = out.uuz;

ew1 = P.*uu1+(1-P).*uu2;
ew2 = (1-Q).*uu1+Q.*uu2;

%% Eq 5, 6


%% Eq 7
betas = out.betas;

e1t = {};
e2t = {};

for i=50:length(betas)
    e1t{i} = (ones(length(betas{i}),1) - betas{i}) * uuz(i) + betas{i} * ew1(i);
    e2t{i} = (ones(length(betas{i}),1) - betas{i}) * uuz(i) + betas{i} * ew2(i);
end

%% Eq 8
sig1 = out.sig1;
sig2 = out.sig2;
vols = out.vols;

Pi1 = {};
Pi2 = {};
Sigma1 = {};
Sigma2 = {};
w1s = [];
w2s = [];

for i=50:length(betas)
    de = e1t{i} - e2t{i};
    V = diag((vols{i}(2:end-1)).^2);
    Pi1{i} = betas{i}*betas{i}'*sig1(i)^2 + V;  %% Eq 8
    Pi2{i} = betas{i}*betas{i}'*sig2(i)^2 + V;  %% Eq 8
    Sigma1{i} = P(i)*Pi1{i} + (1-P(i))*Pi2{i} + P(i)*(1-P(i))*de*de'; %% Eq 9
    Sigma2{i} = (1-Q(i))*Pi1{i} + Q(i)*Pi2{i} + Q(i)*(1-Q(i))*de*de'; %% Eq 9
    
    %% Eq 10
    w1 = inv(Sigma1{i})*e1t{i};
    w2 = inv(Sigma2{i})*e2t{i};
    w1 = w1 / sum(w1);
    w2 = w2 / sum(w2);
    w1s = [w1s w1];
    w2s = [w2s w2];
    
end

%% Backtesting i=50:end

ret1 = sum((w1s'.*Returns(50:end,2:end-1))');
ret2 = sum((w2s'.*Returns(50:end,2:end-1))');
prob = out.prob(50:end);

rets = [];

for i=1:length(prob)-20
    if(prob(i)>0.5)
        rets = [rets ret1(i+1)];
    else
        rets = [rets ret2(i+1)];
    end
end



pld = [];
for i=1:12
    pld = [pld ret2price(Returns(51:end,i))]
end

%% Mean Variance Portfolio & RS Strategy

clf 

plot(ret2price(rets),'k-','LineWidth',3)
hold on
plot(ret2price(Returns(52:end,2:11)*mvw'),'b:','LineWidth',2)
input('Enter')

clf
pld = [pld(:,2:end)]

plot(pld,'-','LineWidth',1);
hold on
%plot(ret2price(rets),'k-','LineWidth',3)


Legends = {'China','Japan','US','EMBI', 'Switzerland','North America', 'EU', 'UK', 'Pacific', 'EM','World','RS Strategy'};
lgd = legend(Legends,'Location','Northwest','FontSize',14);
set(lgd,'FontSize',14);
set(gca, 'Xtick', [0:12:240])
set(gca, 'XtickLabel', [1998:1:2018],'FontSize',12)

%Set Y ticks
a=[cellstr(num2str(get(gca,'ytick')'*100))];
pct = char(ones(size(a,1),1)*'%');
new_yticks = [char(a),pct];
%set(gca,'yticklabel',new_yticks,'FontSize',14)
ylabel('Cumulative Return','FontSize',24)


plot(ret2price(rets),'k-','LineWidth',3)
Legends = {'China','Japan','US','EMBI', 'Switzerland','North America', 'EU', 'UK', 'Pacific', 'EM','World','RS Strategy'};
lgd = legend(Legends,'Location','Northwest','FontSize',14);


input('Enter')


%% Figure B

clf
hold on
plot(out.Spec_Out{end}.smoothProb(50:end,1),'k-');
plot(out.Spec_Out{end}.filtProb(50:end,1),'b--');


legend({'Smooth Probability','Filter Probability'},'Location','Southeast','FontSize',20)

input('Enter')
%% Figure C

