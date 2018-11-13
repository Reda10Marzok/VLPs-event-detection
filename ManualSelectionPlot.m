% MANUAL SELECTION PLOT  :  This script plots all a manuallly selected
%                           events obtained from ManualSelection.
%                           A stack of selected events is given to
%                           visualize an example.

clear;close all;clc

f='./20080515-000000-ETNA-ECPN-N.sac'; %.sac file of the specific station and component.

    K=rsac(f);
    time=K(:,1);
    x=K(:,2);
    header=K(:,3);

fs=round(1/(time(2)-time(1)));   % Sampling frequency

% Downsampling factor 10.

fsp=10;  %Hz

r=fs/fsp;
xsp=downsample(x,r);
tsp=0:1/fsp:time(end);


%% Butterworth 4 poles Filtering
[B,A]=butter(4,[1/30 1/6]/(fsp/2));
[H,f]=freqz(B,A,2048,fsp);


% Filtering
ysp=filter(B,A,xsp);
% Mean subtract
ysp=ysp-mean(ysp);

% Load here the selected events in .mat variable. This variable must
% contain the start and end times for all selected events. 
data=load('VLPsTempECPNStack2.mat')
variables=fields(data)
VLPsTemp=data.(variables{1})



VLP=[]

for i=1:length(VLPsTemp)

ini=round(VLPsTemp(i,1)*fsp);
fin=round(VLPsTemp(i,2)*fsp);
[Ymax,tmax]=max(ysp(ini:fin));
Maxvalue=tmax+ini;
ini=(Maxvalue-30*fsp); 
fin=(Maxvalue+30*fsp);

VLPAux=[ini fin];
 VLP=[VLP ; VLPAux];
 
 t1=tsp(ini:fin);
 t1=t1-t1(1);
 
hold on
plot(t1,ysp(ini:fin))
ylim([-7000 7000]);
xlim([t1(1) t1(end)]);
ylabel('Velocity (m/s)');
set(0,'defaultfigurecolor',[1 1 1])
xlabel('Time(s)')
ylabel('Amplitude')
pause
i=i+1;

end
title('Superimposed events manually selected')

% University of Granada - Final project of the Telecommunication engineering 
% degree - Signal Theory, Telematics and Communications Department (TSTC).
% Student : Reda Marzok Ben Omar.


