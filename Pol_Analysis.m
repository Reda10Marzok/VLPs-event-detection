% POLARIZATION ANALYSIS : This script performs a complete polarization 
%                         analysis from a three component continuous signal
%                         associated to a specific seismic station.



%% Pre-processing. Filtering and Decimation.
%% E component
clear;close all;clc

f='./20080515-000000-ETNA-ECPN-E.sac';  %.sac file of the specific station and component.

    K=rsac(f);
    time=K(:,1);
    x_E=K(:,2);
    header=K(:,3);

fs=round(1/(time(2)-time(1)));      % Sampling frequency

% Downsampling factor 10.

fsp=10;  % Decimation factor can be modified here.

    r=fs/fsp;
    xsp_E=downsample(x_E,r);
    tsp=0:1/fsp:time(end);


%% Butterworth 4 poles Filtering. (E)

[B,A]=butter(4,[1/30 1/6]/(fsp/2)); % The specific bandpass of the filter is configured here
[H,f]=freqz(B,A,2048,fsp);          % to analyse VLP events. In order to examine other type 
                                    % of events, modify here the desired frenquency band.
                                                                
% Filtering
ysp_E=filter(B,A,xsp_E);
% Mean subtract                      % Removing linear trends of the signal
ysp_E=ysp_E-mean(ysp_E);



% Segment of analysis 
    ini=1*fsp;                      % The segment of analysis is set here to examine a whole day
    fin=86400*fsp;                  % (86400 seconds). The segment of analysis can be modified as
                                    % desired, just changing the numbers.
                                    
      
Length=(fin-ini)/fsp;               % Lenght of the downsampled signal

YM=max(abs(ysp_E));
Componente_E=ysp_E(ini:fin);


%% Same process will be aplied to the rest of components 
%% N component

f='./20080515-000000-ETNA-ECPN-N.sac';

    K=rsac(f);
    time=K(:,1);
    x_N=K(:,2);
    header=K(:,3);

fs=round(1/(time(2)-time(1)));      % Sampling frequency

% Downsampling 
    
    xsp_N=downsample(x_N,r);
    tsp=0:1/fsp:time(end);


%% Butterworth 4 poles Filtering. (N)

% Filtering
ysp_N=filter(B,A,xsp_N);
% Mean subtract to remove linear trend 
ysp_N=ysp_N-mean(ysp_N);


YM=max(abs(ysp_N));
Componente_N=ysp_N(ini:fin);


%% Z component

f='./20080515-000000-ETNA-ECPN-Z.sac';

    K=rsac(f);
    time=K(:,1);
    x_Z=K(:,2);
    header=K(:,3);

fs=round(1/(time(2)-time(1)));       % Sampling frequency

% Downsampling

    xsp_Z=downsample(x_Z,r);
    tsp=0:1/fsp:time(end);
    
%% Butterworth 4 poles Filtering. (Z)

% Filtering
ysp_Z=filter(B,A,xsp_Z);
% Mean subtract to remove linear trend
ysp_Z=ysp_Z-mean(ysp_Z);


YM=max(abs(ysp_Z));
Componente_Z=ysp_Z(ini:fin);




%% Covariance matrix and polarization parameters

W=20;                                % Window size. Ajustable parameter.                   
    ini_bucle=1;                     
    fin_bucle=ini_bucle+W*fsp;
    i=1;


while fin_bucle<Length*fsp
    
    

[Cov_Matrix V D]=CovarianceMatrix(Componente_E(ini_bucle:fin_bucle),Componente_N(ini_bucle:fin_bucle),Componente_Z(ini_bucle:fin_bucle));

    % Eigenvalues. 
    Lambda1=D(3,3);
    Lambda2=D(2,2);
    Lambda3=D(1,1);

    % It must be verified that Lambda1 >Lambda2 >Lambda3. Otherwise
    % execution will be interrupted :
    
    if (Lambda2 > Lambda1)
        fprintf('Error in eigenvalues\n')
        pause
    end

    if (Lambda3 > Lambda1)
        fprintf('Error in eigenvalues\n')
        pause
    end

    if (Lambda3 > Lambda2)
        fprintf('Error in eigenvalues\n')
        pause
    end


    % Rectilinearity. Jurkevics 1988.
    RL(i)=1-((Lambda2+Lambda3)/(2*Lambda1));    



    % P-azimuth. Largest eigenvector corresponding to largest Eigenvalue
    % must be used.
    
    Largest_Eigenvector=V(:,3);
    U11=Largest_Eigenvector(1,1);
    U21=Largest_Eigenvector(2,1);
    U31=Largest_Eigenvector(3,1);

    Modulo=norm(Largest_Eigenvector);
    U11=U11/Modulo;
    U21=U21/Modulo;
    U31=U31/Modulo;

    %Azimuth
    P_Azimuth(i)=atand(U11/U21);
    
        if (U11*U21)<0
            P_Azimuth(i)=P_Azimuth(i)+180;
        end

    %Back_Azimuth
    Back_Azimuth(i)=P_Azimuth(i)+180;

    %Incidence angle

    P_incidence(i)=acosd(U31);




ini_bucle=ini_bucle+((W/2)*fsp);     % Increment is set here to be half of the window.       
fin_bucle=ini_bucle+W*fsp;
i=i+1;

end

%% Plotting polarization parameters

    t=ini/fsp:W/2:fin/fsp;                  
    t=t-t(1);
    t(1)=[];
    t(end)=[];

    % X axis will be shown in 'YYYY-MM-DD HH:MM:SS' format.
    Time={};
    n=length(t);
        for u=1:n
        Time(u,1)={datestr(seconds(t(u)),'2008-05-15 HH:MM:SS')}; %Modify here the day of analysis if it is needed.
        end

t = datetime(Time,'InputFormat','yyyy-MM-dd HH:mm:ss');

figure(1)
subplot(311)
plot(t,RL,'k.');
title('Rectilinearity');
set(gca, 'fontsize', 12,'XTick', [])
set(gca, 'YTick', 0:0.2:1)
xlim([t(1) t(end)])

subplot(312)
plot(t,P_Azimuth,'.','Color', [88/255 7/255  7/255]);
ylim([0 200])
xlim([t(1) t(end)])
title('Azimuth');
set(gca, 'YTick', 0:40:200)
set(gca, 'fontsize', 12,'XTick', [])

subplot(313)
plot(t,P_incidence,'.','Color', [0/255 2/255  53/255]);
title('Incidence');
xlim([t(1) t(end)])
set(gca, 'YTick', 0:40:200)
set(0,'defaultfigurecolor',[1 1 1]);
set(gca, 'fontsize', 12);



%% Plotting the azimuth, Back azimuth and incidence histogram

% Circular Statistics Toolbox for Matlab is used here.
% Circ_plot function By Philipp Berens & Marc J. Velasco, 2009.
% Fourth argument of circ_plot specify the number of bins desired.
    figure (2)

    alpha = P_Azimuth;
    alpha=deg2rad(alpha);
    circ_plot(alpha,'hist',[],120,false);
    title('Azimuth Histogram')
    view([-90 90])
    set(gca, 'fontsize', 12)
    
    
    figure (3)

    alpha = Back_Azimuth;
    alpha=deg2rad(alpha);
    circ_plot(alpha,'hist',[],120,false);
    title('Backazimuth Histogram')
    view([-90 90])
    set(gca, 'fontsize', 12)
    
    figure (4)

    alpha = P_incidence;
    alpha=deg2rad(alpha);
    circ_plot(alpha,'hist',[],120,false);
    title('Incidence histogram')
    view([-90 90])
    set(gca, 'fontsize', 12)


% University of Granada - Final project of the Telecommunication engineering 
% degree - Signal Theory, Telematics and Communications Department (TSTC).
% Student : Reda Marzok Ben Omar.

