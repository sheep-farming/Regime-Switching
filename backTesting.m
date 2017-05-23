load('new.mat') % RAW RETURN DATA FILE
load('est.mat') % ESTIMATED PARAMETERS DATA FILE FROM mainEstimation.m
%load('unconOptimal.mat')

out.uuz = Returns(:,1); % Assume Risk Free Rate = JPM 3 Mo Cash Total Return Index

%% See Paper https://www0.gsb.columbia.edu/faculty/aang/papers/inquire.pdf

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
    %% Eq 8
    Pi1{i} = betas{i}*betas{i}'*sig1(i)^2 + V;  %% Eq 8
    Pi2{i} = betas{i}*betas{i}'*sig2(i)^2 + V;  %% Eq 8
    %% Eq 9
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

%% Backtesting from Month 50 to end


prob = out.prob(50:end);
ret1 = sum(w1s(:,1:end-1).*Returns(51:end,2:11)');
ret2 = sum(w2s(:,1:end-1).*Returns(51:end,2:11)');

rets = [];
ws = [];
for i=1:length(w1s')-1
    if(prob(i)>0.5)
        
        rets = [rets ret1(i)];
        ws =[ws, w1s(i)];    
    else
        rets = [rets ret2(i)];
        ws =[ws, w2s(i)];
    end
end



pld = [];
for i=1:12
    pld = [pld ret2price(Returns(51:end,i))];
end




%% Analyse Turnover

tTurnover = 0;
for i=1:length(ws)-1

    moTurnover = sum(abs(ws(:,i)-ws(:,i+1)))/2;
    tTurnover = tTurnover + moTurnover;

end
yearlyTurnover = tTurnover*12/(length(ws)) %% Annualised Turnover of RS Strategy


%% Analyse Sharpe
meanReturn = mean(rets)*12
stdReturn = std(rets)*sqrt(12)
RSSharpe = (meanReturn- Returns(end,1))/ stdReturn

meanReturn = mean(Returns(50:end-1,12))*12
stdReturn = std(Returns(50:end-1,12))*sqrt(12)
worldSharpe = (meanReturn - Returns(end,1))/ stdReturn





%% Plotting
input('Press Enter to Plot Figures');


plotting


