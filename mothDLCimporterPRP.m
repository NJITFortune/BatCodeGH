% c = readmatrix('C:\Users\DeepLabCut\Desktop\fem12-2022-test\fem12_bt_25DLC_resnet50_Fem12restestFeb8shuffle1_110000_filtered.csv');
% c = readmatrix('C:\Users\DeepLabCut\Desktop\PamelaMoth2022\UltraMoth-PRP-2022-02-15\videos\fem7_tp_20-008DLC_resnet50_UltraMothFeb15shuffle1_210000.csv');

[filer,pather] = uigetfile('*.csv');
c = readmatrix(fullfile(pather,filer));


%% Put the data into a human-readable structure
Fs = 1200; % The framerate of the high-speed video capture

% Make a temporary time series - this will be replaced by the TTL markers
% in the Spike2 file that captured the audio recording.
tim = 1/Fs:1/Fs:length(c(:,1))/Fs;

% Read the data into a human-readable structure
% Columns are X, Y, Probability

% Reye 2,3,4
% Leye 5,6,7
% Rantenna 8,9,10
% Lantenna 11,12,13
% Rlegtip 14,15,16
% Llegtip 17,18,19
% Rwingbase 20,21,22
% Rwingtip 23,24,25
% Rwingsplit 26,27,28
% Lwingbase 29,30,31
% Lwingtip 32,33,34
% Lwingsplit 35,36,37
% Tail 38,39,40

%idx.reye = 2; idx.leye = 5;
    m.reye(:,1) = c(:,2); m.reye(:,2) = c(:,3); mreye(:,3) = c(:,4);
    m.leye(:,1) = c(:,5); m.leye(:,2) = c(:,6); mleye(:,3) = c(:,7);

%idx.rAntenna = 8; idx.lAntenna = 11;    
    m.rAntenna(:,1) = c(:,8); m.rAntenna(:,2) = c(:,9); m.rAntenna(:,3) = c(:,10);
    m.lAntenna(:,1) = c(:,11); m.lAntenna(:,2) = c(:,12); m.lAntenna(:,3) = c(:,13);

%idx.rlegTip = 14; idx.llegTip = 17;
    m.rlegTip(:,1) = c(:,14); m.rlegTip(:,2) = c(:,15); m.rlegTip(:,3) = c(:,16);
    m.llegTip(:,1) = c(:,17); m.llegTip(:,2) = c(:,18); m.llegTip(:,3) = c(:,19);
    
%idx.rWingbase = 20; idx.rWingtip = 23; idx.rWingsplit = 26;   
    m.rWingbase(:,1) = c(:,20); m.rWingbase(:,2) = c(:,21); m.rWingbase(:,3) = c(:,22);
    m.rWingtip(:,1) = c(:,23); m.rWingtip(:,2) = c(:,24); m.rWingtip(:,3) = c(:,25);
    m.rWingsplit(:,1) = c(:,26); m.rWingsplit(:,2) = c(:,27); m.rWingsplit(:,3) = c(:,28);

%idx.lWingtip = 29; idx.lWingbase = 32; idx.lWingsplit = 35; 
    m.lWingbase(:,1) = c(:,29); m.lWingbase(:,2) = c(:,30); m.lWingbase(:,3) = c(:,31);
    m.lWingtip(:,1) = c(:,32); m.lWingtip(:,2) = c(:,33); m.lWingtip(:,3) = c(:,34);
    m.lWingsplit(:,1) = c(:,35); m.lWingsplit(:,2) = c(:,36); m.lWingsplit(:,3) = c(:,37);

%idx.Tail = 38;
    m.Tail(:,1) = c(:,38); m.Tail(:,2) = c(:,39); m.Tail(:,3) = c(:,40);


figure(2); clf; hold on;
    plot(m.reye(:,1),-m.reye(:,2), '.')
    plot(m.leye(:,1),-m.leye(:,2), '.')
    plot(m.rWingtip(:,1),-m.rWingtip(:,2), '.')
    plot(m.lWingtip(:,1),-m.lWingtip(:,2), '.')
    
    tailidx = find(m.Tail(:,3) > 0.9);
    plot(m.Tail(tailidx,1),-m.Tail(tailidx,2), '.')

    plot(m.rlegTip(:,1), -m.rlegTip(:,2), '.');
    plot(m.llegTip(:,1), -m.llegTip(:,2), '.');
    

    
figure(3); clf;
    ax(1) = subplot(311); hold on; ylabel('Wingtip Y with tail')
        plot(tim, m.rWingtip(:,2)); 
        plot(tim, m.lWingtip(:,2));
        plot(tim(tailidx), m.Tail(tailidx,2), 'k.')        
        
    ax(2) = subplot(312); 
        specgram(m.lWingtip(:,2) - mean(m.lWingtip(:,2)), 1024, 1200, [], 1000);
        ylim([0 30]); colormap(flipud(gray)); caxis([75 90]);
        ylabel('Wingbeat freq, Hz');
    ax(3) = subplot(313); hold on;
        R = m.rlegTip(:, 1:2);
        L = m.llegTip(:, 1:2);
        for j = length(R):-1:1
            dist(j) = pdist2(R(j,:),L(j,:));
        end
        R = m.rlegTip(:, 1:2);
        L = m.llegTip(:, 1:2);
        ff = 0;
        for j = length(R):-1:1
            if m.rlegTip(j,3) > 0.5 && m.llegTip(j,3) > 0.5
                ff = ff+1;
                dist(j) = pdist2(R(j,:),L(j,:));
            end
        end
        plot(tim(find(dist)), dist, '.');
        ylabel('Feet distance, px')    
    linkaxes(ax, 'x');

    
%% Initial plot
figure(1); clf; hold on;

for j = 2:3:length(c(1,:))
    
    plot(c(:,j), -c(:,j+1), '.');

end

