%% RAmen!
%                            _ooOoo_
%                           o8888888o
%                           88" . "88
%                           (| -_- |)
%                            O\ = /O
%                        ____/`---'\____
%                      .   ' \\| |// `.
%                       / \\||| : |||// \
%                     / _||||| -:- |||||- \
%                       | | \\\ - / | |
%                     | \_| ''\---/'' | |
%                      \ .-\__ `-` ___/-. /
%                   ___`. .' /--.--\ `. . __
%                ."" '< `.___\_<|>_/___.' >'"".
%               | | : `- \`.;`\ _ /`;.`/ - ` : | |
%                 \ \ `-. \_ __\ /__ _/ .-` / /
%         ======`-.____`-.___\_____/___.-`____.-'======
%                            `=---='
%
%         .............................................
%                Dear Mr. Buddha No Bug Please
%
%  Created by Du Pupu on 19/5/17.
%

load('returns.mat')

meanReturns = mean(Returns)*12

%% Regression Test for MXEM v MXWO
% r: reg value; m: slope; b: offset

rAst = Returns(:,1)';
rMkt = Returns(:,12)';

[corr,astBeta,astSlope]=regression(rMkt,rAst)

%plot(rMkt,rAst,'r.')
%hold on
%plot([max(rMkt),min(rMkt)],[max(rMkt)*astBeta+astSlope,min(rMkt)*astBeta+astSlope])

betas = [];

for i=1:11
    rAst = Returns(:,i)';
    [corr,astBeta,astSlope]=regression(rMkt,rAst);
    betas = [betas astBeta];
end


betas = [betas 1];
expReturns = betas * meanReturns(12);
plot(betas,meanReturns, 'rx')
hold on
[corr,astBeta,astSlope]=regression(betas,meanReturns)
plot([max(betas),min(betas)],[max(betas)*astBeta+astSlope,min(betas)*astBeta+astSlope])