% VLP DETECTION : This  script  implements  the  VLP  events detection
%                 given a three component continuous signal associated 
%                 to a specific seismic station.



clear;close all;clc
%% Pre-processing. Filtering and Decimation.
%% E component 
clear;close all;clc

f='./20080515-000000-ETNA-ECPN-E.sac';  %.sac file of the specific station and component.

    K=rsac(f);
    time=K(:,1);
    x_E=K(:,2);
    header=K(:,3);

fs=round(1/(time(2)-time(1)));          % Sampling frequency

% Downsampling factor 10.

fsp=10;                                 % Decimation factor can be modified here.

    r=fs/fsp;
    xsp_E=downsample(x_E,r);
    tsp=0:1/fsp:time(end);

%% Butterworth 4 poles Filtering. (E)

[B,A]=butter(4,[1/30 1/6]/(fsp/2));     % The specific bandpass of the filter is configured here
[H,f]=freqz(B,A,2048,fsp);              % to analyse VLP events. In order to examine other type 
                                        % of events, modify here the desired frenquency band.


% Filtering
ysp_E=filter(B,A,xsp_E);
% Mean subtract                         % Removing linear trends of the signal
ysp_E=ysp_E-mean(ysp_E);
Componente_E=ysp_E;

%% N component

f='./20080515-000000-ETNA-ECPN-N.sac';

    K=rsac(f);
    time=K(:,1);
    x_N=K(:,2);
    header=K(:,3);

fs=round(1/(time(2)-time(1)));          % Sampling frequency

% Downsampling 

    xsp_N=downsample(x_N,r);
    tsp=0:1/fsp:time(end);


%% Butterworth 4 poles Filtering. (N)
[B,A]=butter(4,[1/30 1/6]/(fsp/2));
[H,f]=freqz(B,A,2048,fsp);

% Filtering
ysp_N=filter(B,A,xsp_N);
% Mean subtract to remove linear trend 
ysp_N=ysp_N-mean(ysp_N);
Componente_N=ysp_N;


%%  Z component

f='./20080515-000000-ETNA-ECPN-Z.sac';

    K=rsac(f);
    time=K(:,1);
    x_Z=K(:,2);
    header=K(:,3);

fs=round(1/(time(2)-time(1)));          % Sampling frequency


% Downsampling

    xsp_Z=downsample(x_Z,r);
    tsp=0:1/fsp:time(end);

%% Butterworth 4 poles Filtering. (Z)
[B,A]=butter(4,[1/30 1/6]/(fsp/2));
[H,f]=freqz(B,A,2048,fsp);

% Filtering
ysp_Z=filter(B,A,xsp_Z);
% Mean subtract to remove linear trend
ysp_Z=ysp_Z-mean(ysp_Z);
Componente_Z=ysp_Z;

%% ENVELOPE DISCRIMINATION FILTER : VLP CANDIDATES SELECTION

 e=abs(hilbert(Componente_Z));                      % Creation of the envelope. We specify here in which
%e=abs(hilbert(Componente_E));                      % component we decide to do the detection. Z Component 
%e=abs(hilbert(Componente_N));                      % is set as default.
                                        
                                
% Creation of the impulse response : Defined by a mexican hat Wavelet.

L=20;                                               % Central lobe length of the mexican hat wavelet.  
                                                    % It must match with  the maximum length of
                                                    % the events that we we want to detect.                     
                                                    
th=-2*L:1/fsp:2*L;
h=(1-(th/(L/2)).^2).*exp(-0.5*(th/(L/2)).^2);
h(h<0)=h(h<0)/-min(h(h<0));
beta=1.4;                                           % Penalization factor. It can be modified in order
                                                    % to obtain more or less event candidates.                                                                                                      
h(h<0)=beta*h(h<0);   
figure(1);
plot(th,h,'k');
title(['Impulse response beta = ' num2str(beta)]);

CF=filter(h,1,e);                                   % Characteristic function result from the convolution.
M=(length(h)-1)/2;                                  % group delay
CF=[CF(M+1:end) ; zeros(M,1)];                      % delay compensation 
CF(1:M)=0;
CF(CF<0)=0;

% It is necessary to find now the maximum results obtained from the characteristic
% function, and the correponding times. (VLP Event candidates)
[pks,locs] = findpeaks(CF,'MinPeakDistance',round(40*fsp));
figure(2);
plot(tsp,CF,'k',tsp(locs),pks,'r*');
title('Main Peaks of the output filter signal')
ylabel('Amplitude');
xlabel('Time(s)');
Tp=tsp(locs);                                        % VLP event candidates time vector
n=length(Tp);                                        % Number of VLP event candidates   

fprintf('Number of event candidates : %d\n',n)


Vtemps=[];

% Selection of 60 seconds signal segment for every candidate obtained:
for n=1:length(Tp)    
    ini=max(1,round((Tp(n)-30)*fsp));
    fin=min(length(Componente_Z),round((Tp(n)+30)*fsp));
    t1=tsp(ini:fin);
    Vtemps=[Vtemps; t1];
    t1=t1-t1(1);
    % Uncomment the following paragraph to visualize every event candidate
    % and its spectrogram :
%           figure(3)
%           subplot(311)
%           plot(t1,Componente_Z(ini:fin),'color','k');
%           xlim([t1(1) t1(end)])
%           xlabel('Time(s)')
%           ylabel('Amplitude')
%           set(gca, 'fontsize', 12)
%           title(['Evento: ' num2str(n) '/' num2str(length(Tp))])
%           subplot(312);
%           spectrogram(Componente_Z(ini:fin),128,126,1024,fsp,'yaxis');
%           ylim([0 1]);
%           subplot(313);
%           psd=abs(fft(Componente_Z(ini:fin),1024));
%           psd=psd(1:512);
%           psd=psd/max(psd);
%           freq=0:5/512:5-1/512;
%           plot(freq,psd,'k','Linewidth',1.2)
%           xlim([0 1])
%           pause
end

%% POLARIZATION FILTER : Discrimination of event candidates 
v=length(Vtemps);
W=20;                               % Window size 
for j=1:n

    ini=round(Vtemps(j,1))*fsp;
    fin=round(Vtemps(j,v))*fsp;

    t1=tsp(ini:fin);
    t1=t1-t1(1);
 
    Length=(fin-ini)/fsp ;

    % Window selection for every candidate 
    [Ymax,tmax]=max(Componente_Z(ini:fin));
    Maxvalue=tmax+ini;
    iniWindow=round(Maxvalue-(W/2)*fsp); 
    finWindow=round(Maxvalue+(W/2)*fsp);
    t2=tsp(iniWindow:finWindow);
    t2=t2-t2(1);
    
    % Uncomment the following paragraph to visualize every selected window
    % for every candidate:
    
%     figure(4);
%     plot(t2,Componente_Z(iniWindow:finWindow),'k')
%     xlim([t2(1) t2(end)]);
%     title(['Selected Window, Candidate : '  num2str(j) '/' num2str(n)]);
%     xlabel('Time(s)');
%     set(gca, 'fontsize', 12)
%     pause


%% Looking for the polarization parameters for the selected Window

ini_bucle=iniWindow;
fin_bucle=finWindow;


[Cov_Matrix V D]=CovarianceMatrix(Componente_E(ini_bucle:fin_bucle),Componente_N(ini_bucle:fin_bucle),Componente_Z(ini_bucle:fin_bucle));


Lambda1=D(3,3);
Lambda2=D(2,2);
Lambda3=D(1,1);

% Rectilinearity RL
RL(j)=1-((Lambda2+Lambda3)/(2*Lambda1));



% P-azimuth
Largest_Eigenvector=V(:,3);
U11=Largest_Eigenvector(1,1);
U21=Largest_Eigenvector(2,1);
U31=Largest_Eigenvector(3,1);

Modulo=norm(Largest_Eigenvector);
U11=U11/Modulo;
U21=U21/Modulo;
U31=U31/Modulo;

%Azimuth
P_Azimuth(j)=atand(U11/U21);
if (U11*U21)<0
    P_Azimuth(j)=P_Azimuth(j)+180;
end

%Back_Azimuth
Back_Azimuth(j)=P_Azimuth(j)+180;

% Incidence angle

P_incidence(j)=acosd(U31);


 end
 
% Plotting the azimuth histogram
alpha = P_Azimuth;
    alpha=deg2rad(alpha);
    figure(5)
    set(0,'defaultfigurecolor',[1 1 1])
    circ_plot(alpha,'hist',[],50,false,false,'linewidth',1,'color','r');
    title('Azimuth histogram','Fontsize',8)
    view([-90 90])
    
    % Plotting the Backazimuth histogram
alpha = Back_Azimuth;
    alpha=deg2rad(alpha);
    figure(6)
    set(0,'defaultfigurecolor',[1 1 1])
    circ_plot(alpha,'hist',[],50,false,false,'linewidth',1,'color','r');
    title('BackAzimuth histogram','Fontsize',8)
    view([-90 90])
    
    % Plotting the incidence histogram
    
alpha = P_incidence;
    alpha=deg2rad(alpha);
    figure(7)
    set(0,'defaultfigurecolor',[1 1 1])
    circ_plot(alpha,'hist',[],50,false,false,'linewidth',1,'color','r');
    title('Incidence histogram','Fontsize',8)
    view([-90 90])
    
%% Creation of a string with datetime format to plot the results 
VLP_Time={};
for u=1:n
VLP_Time(u,1)={datestr(seconds(Tp(1,u)),'2008-05-15 HH:MM:SS')}; %Change here the year if it is necessary
end

t = datetime(VLP_Time,'InputFormat','yyyy-MM-dd HH:mm:ss');



% Rectilinearity of all event candidates, showing the threeshold used in the
% automatic detection
figure(8)
scatter(t,RL,20,'k','filled')
x = [t(1) t(end)];
y = [0.8 0.8];
line(x,y,'Color','red','LineStyle','--','Linewidth',1.4)
 ylim([0 1])
 ylabel('Rectilinearity');
 set(gca, 'fontsize', 12)
 
 
 
% Azimuth of all event candidates
figure(9)
scatter(t,P_Azimuth,20,'m','filled')
ylim([0 360])
ylabel('Azimuth');
set(gca, 'fontsize', 12)

% Incidence of all event candidates
figure(10)
scatter(t,P_incidence,20,'g','filled')
ylim([0 360])
ylabel('Incidence');
set(0,'defaultfigurecolor',[1 1 1])

% Backazimuth of all event candidates, showing the threshold used in the
% automatic detection
figure(11)
scatter(t,Back_Azimuth,20,'filled','MarkerFaceColor',[118/255 6/255 69/255])
ylim([0 360])
ylabel('Backazimuth');
x = [t(1) t(end)];
y = [300 300];
line(x,y,'Color','red','LineStyle','--','Linewidth',1.4)                    % Threshold ECPN 
x1 = [t(1) t(end)];
y1 = [350 350];
line(x1,y1,'Color','red','LineStyle','--','Linewidth',1.4)                  % Threshold ECPN
set(0,'defaultfigurecolor',[1 1 1])
set(gca, 'fontsize', 12)
 
 
 
RLcounter=0;
for i=1:n
    if RL(i)>0.8
        RLcounter=RLcounter+1;
    end
end
RLcounter;
%% Polarization threshold 
BA=Back_Azimuth;
VLPnumbers=0;
VLPindex=[];
RLdetected=[];
BAdetected=[];
for i=1:n
    if ((BA(i)<= 350)&&(BA(i) >= 300)&&(RL(i)>=0.8))                    % BackAzimuth threshold is configured here to detect events in ECPN seismic station.
        VLPnumbers=VLPnumbers+1;                                        % If the detection is going to be applied in other station this threshold must be changed.
        VLPindex=[VLPindex ; i ];                                       % EBEL seismic station Threshold: if((P_Azimuth(i)<= 60)&&(P_Azimuth(i) >= 10)&&(RL(i)>=0.8))
        RLdetected=[ RLdetected ; RL(i) ];                              % EPLC seismic station Threshold: if ((BA(i)<= 245)&&(BA(i) >= 195)&&(RL(i)>=0.8))
        BAdetected=[ BAdetected ; BA(i) ];                              % EPDN seismic station Threshold: if((P_Azimuth(i)<= 165)&&(P_Azimuth(i) >= 115)&&(RL(i)>=0.8))    
    end
end
    

fprintf('Number of VLPs detected : %d\n',VLPnumbers)

VLPnumbers;                                                             % Numbers of VLP detected
VLPindex;                                                               % Index of the VLPs detected (With respect to all the candidates obatined)


% University of Granada - Final project of the Telecommunication engineering 
% degree - Signal Theory, Telematics and Communications Department (TSTC).
% Student : Reda Marzok Ben Omar.
