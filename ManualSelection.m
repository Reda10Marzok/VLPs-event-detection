% MANUAL SELECTION : This script performs a manual selection of events.
%                    Firstly, a specified component from a seismic station 
%                    will be pre-processed (Downsampling and Filtering). 
%                    After this, a window of size W will be shown to the 
%                    user in order to select the VLP events. VLPsTemp is
%                    the variable that will contain the times chosen by the
%                    user. This is not an automatic detection.
%                    
%                   

clear;close all;clc

f='./20080515-000000-ETNA-ECPN-Z.sac'; %.sac file of the specific station and component

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


%% Butterworth 4 poles Filtering. 
    [B,A]=butter(4,[1/30 1/6]/(fsp/2));
    [H,f]=freqz(B,A,2048,fsp);


%Filtering
ysp=filter(B,A,xsp);
%Mean subtract 
ysp=ysp-mean(ysp);

%Windowing
    W=60;       % Window size in seconds
    ini=1*fsp;
    fin=ini+W*fsp;
 
    fig=figure(1);
    YM=max(abs(ysp));
    VLPsTemp=[];
    
while fin<length(ysp)
    
    plot(tsp(ini:fin),ysp(ini:fin),'b');
    ylim([-YM YM]);
    xlim([tsp(ini) tsp(fin)]);
    choice = menu('Is there any VLP event ?','Yes','No','Exit') ;
    
    if choice==1
        [x,y]=ginput(2);
        VLPsTempAux=[x(1) x(2)];
        VLPsTemp=[VLPsTemp ; VLPsTempAux];
        ini=ini+ (W*fsp);
        fin=ini+W*fsp;
        
    elseif choice==2
            ini=ini+ (W*fsp);
            fin=ini+W*fsp;
    else
        close(fig)
         break 
    end
    
end


% University of Granada - Final project of the Telecommunication engineering 
% degree - Signal Theory, Telematics and Communications Department (TSTC).
% Student : Reda Marzok Ben Omar.