clear;clc;clf;

load('new.mat');

%Fuck All Bonds

%Returns(:,2)=Returns(:,5);
%Returns(:,3)=Returns(:,6);


addpath('MS_Regress-Matlab')
addpath('MS_Regress-Matlab/m_Files'); % add 'm_Files' folder to the search path

%dateStart = datenum(Dates(1),'dd/mm/yyyy');
%dateEnd = datenum(Dates(end),'dd/mm/yyyy');
%% Variables
%
% P Q <- Switching Probability
% uu1 uu2 <- World Excess Return
%

%% Assumptions 
%
% uuz = uRf
Returnz = Returns;
% figure(3);
% hold on;
% for i=1:12
%     plot(1:length(Returns)+1,ret2price(Returns(:,i),100),'--');
% end

%% Regression Parameters
    

k=2;                                % Number of States
S=[1 1];                        % Defining which parts of the equation will switch states (column 1 and variance only)
advOpt.distrib='Normal';            % The Distribution assumption ('Normal', 't' or 'GED')
advOpt.std_method=1;                % Defining the method for calculation of standard errors. See pdf file for more details
%advOpt.printIter='1';
%advOpt.printOut='1';

%% Variables in Iteration


mVal = ones(1,50);
sVal = ones(1,50);
prob = ones(1,50);
%input('Enter');

%% Iteration to Oldie

out = {};
out.P=[];
out.Q=[];
out.p1=[];
out.p2=[];
out.uu1=[];
out.uu2=[];
out.sig1=[];
out.sig2=[];
out.w1={};
out.w2={};
out.ER={};
out.betas = {};
out.betases = {};
out.covAMs={};
out.vols={};
out.uuz = Returns()

out.condER=[];
out.condER=[];

out.idioVols = {};
out.prob=[];
out.Spec_Out={};


i=278;
ws = [];
while 0
end
%for i=50:length(Returnz)
    
    Returns = Returnz(1:i,:);

    uuz = Returns(:,1); % Risk Free Rate, assumed to be Zero Beta Premium

    wRet = Returns(:,12); 
    wExRet = wRet - Returns(:,1);
    constVec=ones(length(wExRet),1); % A constant vector in mean equation (just an example of how to do it)
    indep = [constVec];
    % Regression

    [Spec_Out]=MS_Regress_Fit(wExRet,indep,k,S,advOpt); % Estimating the model
    prob = [prob Spec_Out.filtProb(end,1)];
    %% Plot Return & Probabilities
    mVal=[mVal mVal(end)*(1+Returns(end,end))]
    

%     if(Spec_Out.filtProb(end,1)<.5)
%         hold on
%         plot(i,mVal(end),'r+');
%     else
%         hold on
%         plot(i,mVal(end),'g+');
%     end

    P = Spec_Out.Coeff.p(1,1);
    Q = Spec_Out.Coeff.p(2,2);
    
    out.prob(i)=Spec_Out.filtProb(end,1);
    out.P(i)=P;
    out.Q(i)=Q;
    out.Spec_Out{i}=Spec_Out;
    
%     figure(2);
%     subplot(2,1,1);
%     plot(1:i+1,mVal);
%     subplot(2,1,2);
%     plot(1:i+1,prob);
    

    


    %% Mu and Sigma for Each Regimes

    uus = Spec_Out.Coeff.S_Param{1};
    out.uu1(i)=uus(1);
    out.uu2(i)=uus(2);
    
    sis = [sqrt(Spec_Out.Coeff.covMat{1}) sqrt(Spec_Out.Coeff.covMat{2})];
    out.sig1(i)=sqrt(Spec_Out.Coeff.covMat{1});
    out.sig2(i)=sqrt(Spec_Out.Coeff.covMat{2});

    %% Betas of Assets

    [betas, ses, vols,covAMs] = getBetas(Returns);
    out.betas{i}=betas(2:end-1)';
    out.betases{i}=ses(2:end-1)';
    out.vols{i}=vols';
    out.covAMs{i} = covAMs;
    
    
    

    %% Resulted Values

    condER = [P*uus(1)+(1-P)*uus(2) (1-Q)*uus(1)+Q*uus(2)];
    condVar = [P*sis(1)^2+(1-P)*sis(2)^2+P*(1-P)*(uus(1)-uus(2))^2 (1-Q)*sis(1)^2+Q*sis(2)^2+Q*(1-Q)*(uus(1)-uus(2))^2];
    condStd = sqrt(condVar);
    out.condER(i,:)=condER;
    out.condStd(i,:)=condStd;

    %% Asset Expected Returns for Next Period given Current State of Regime 

    ER = uuz(end) + (condER'-uuz(end))*betas;

    %% Asset Idiosyncratic Risk (StD) by Deducting Market Vol * Indi Beta from Indi Vol
    mVol = vols(12);

    idioVols = sqrt(var(Returns - wRet * betas));
    
    
    
    % Sigs = condStd'*betas+ones(2,1)*idioVols

    %% Estimate Covariance Matrix
    % Kill 1 2

    idioVols = idioVols(2:end-1);
    out.idioVols{i}=idioVols;
    vbetas = betas(2:end-1)';
    ER = ER(:,2:end-1);

    V = diag(idioVols.^2);

    Om1 = (vbetas*vbetas')*(sis(1))^2+V;
    Om2 = (vbetas*vbetas')*(sis(2))^2+V;

    %% Above OK


    Sigma1 = P*Om1+(1-P)*Om2+P*(1-P)*(ER(1,:)-ER(2,:))'*(ER(1,:)-ER(2,:))
    Sigma2 = (1-Q)*Om1+Q*Om2+Q*(1-Q)*(ER(1,:)-ER(2,:))'*(ER(1,:)-ER(2,:))

    w1 = inv(Sigma1)*ER(1,:)';
    w2 = inv(Sigma2)*ER(2,:)';
    w1=(ones(1,length(w1))*w1)^-1*w1; % w1/sum(w1)
    w2=(ones(1,length(w2))*w2)^-1*w2; % w1/sum(w1)
    out.w1{i}=w1;
    out.w2{i}=w2;
    out.ER{i}=ER;
    
    ws = [ws,w1,w2];

    
    %figure(4)


%     if(prob(end)<.5)
%         sret = Returnz(i+1,2:end-1)*w1;
%     else
%         sret = Returnz(i+1,2:end-1)*w2;
%     end

    %sVal=[sVal sVal(end)*(1+sret(end))];

%     plot(sVal,'g-');
%     hold on;
%     plot(mVal,'k--');

    i=i+1
%         if(mod(i,48)==0)
%             %input('5 Times Break');
%         end
%end


%% ER: Expected Returns, Sigs: Expected Volatilities for Assets (1,2,..., 12) under Regime 1 / 2


%% Plot Strategy on One



