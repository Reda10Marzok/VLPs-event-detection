%% Using the CircHist Class
%% Plot Distribution Data
% Generate a noisy sample (von Mises distribution with |theta| == 90°).
rng default
sDist = rad2deg(circ_vmrnd(pi/2,2,100)); % convert to degree
nBins = 36; % number of bins, makes bin size of 10°
%%
% Plot the circular histogram:
obj1 = CircHist(sDist,nBins);
%%
% Adjust appearance:
obj1.colorBar = 'k'; 
obj1.avgAngH.LineStyle = '--';
obj1.avgAngH.LineWidth = 1;
obj1.colorAvgAng = [.5 .5 .5];
% remove offset between bars and plot-center
rl = rlim;
obj1.setRLim([0,rl(2)]);
% draw circle at r == 0.5 (where r == 1 would be the outer plot edge)
rl = rlim;
obj1.drawCirc((rl(2)-rl(1))/2,'--b','LineWidth',2)
obj1.scaleBarSide = 'right';
obj1.polarAxs.ThetaZeroLocation = 'right';
obj1.setThetaLabel('Direction','bottomleft');
% draw resultant vector r as arrow
delete(obj1.rH)
obj1.drawArrow(obj1.avgAng,obj1.r * range(rl),'HeadWidth',10,'LineWidth',2,'Color','r')
% Change theta- and rho-axis ticks
obj1.polarAxs.ThetaAxis.MinorTickValues = [];
thetaticks(0:90:360);
rticks(0:4:20);
obj1.drawScale; % update scale bar
%% Plot Multi-Sample Distribution
% Generate another noisy sample with a different distribution-width |kappa|.
rng default
s2Dist = rad2deg(circ_vmrnd(pi/2,1.5,100));
sCell = {sDist,s2Dist}; % pack both samples into a cell-array
figure
CircHist(sCell,nBins);
%% Combine Multiple Histograms in One Figure
% Create subplot, note that the created |axes| must be a |polaraxes|.
fH = figure;
subplot(1,2,1,polaraxes);
CircHist(sDist,nBins);
subplot(1,2,2,polaraxes);
CircHist(s2Dist,nBins);
% Adjust figure-window size
drawnow
fH.Position(3) = 1100; % width
fH.Position(4) = 500; % height
%% Plot Already-Binned Data
% Bin the generated multi-sample distribution before plotting.
edges = 0:10:360;
histData = histcounts(mod([sDist;s2Dist],360),edges);
figure
CircHist(histData,edges,'dataType','histogram');
%% Axial Data
% Copy the von Mises data with an offset of 180° to generate an axial, bimodal
% distribution.
sAxial = [sDist;sDist+180];
%%
% Call |CircHist| with |'areAxialData'| specified as |true|.
figure
CircHist(sAxial,nBins,'areAxialData',true);
%%
% Note that now the average angle is indicated by an axis that halves the diagram at this
% angle.
%% Draw Arrows
rng default
arrowLen = randn(numel(sDist),1); % random arrow lengths
arrowLen = arrowLen / max(arrowLen);
arrowLen = arrowLen + abs(min(arrowLen));
figure
obj2 = CircHist([1,2],36); % dummy data
delete([obj2.avgAngH;obj2.avgAngCiH(:);obj2.barH(:);obj2.rH]); % remove dummy data
title('');
obj2.scaleBar.Label.String = 'Vector length';
obj2.polarAxs.ThetaAxis.MinorTickValues = [];
thetaticks(0:90:360);
arrowH = obj2.drawArrow(sDist,arrowLen);
%%
% Change visual properties and add another arrow.
set(arrowH,'HeadStyle','plain','HeadWidth',3)
% Draw a single arrow ending at the outer plot edge
obj2.drawArrow(100,[],'Color','r','LineWidth',3)
%%
drawnow % (else, the last figure is not shown in the published version for some reason)