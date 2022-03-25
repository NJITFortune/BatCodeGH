%c = readmatrix('c:\Users\DeepLabCut\Downloads\Fem12-PRP-2021-10-22\videos\fem12_lf_10goodDLC_resnet152_Fem12Oct22shuffle1_1030000_filtered.csv');
%c = readmatrix('c:\Users\DeepLabCut\Downloads\Fem12-PRP-2021-10-22\videos\fem12_tp_09goodDLC_resnet50_Fem12Oct22shuffle1_850000_filtered.csv');
c = readmatrix('c:\Users\DeepLabCut\Desktop\fem12\fem12_tp_05DLC_resnet50_Fem12Oct22shuffle1_150000_filtered.csv');

%% Put the data into a human-readable structure
Fs = 1200; % The framerate of the high-speed video capture

% Make a temporary time series - this will be replaced by the TTL markers
% in the Spike2 file that captured the audio recording.
tim = 1/Fs:1/Fs:length(c(:,1))/Fs;

% Read the data into a human-readable structure
% Columns are X, Y, Probability

idx.reye = 2; idx.leye = 5;
    m.reye(:,1) = c(:,2); m.reye(:,2) = c(:,3); mreye(:,3) = c(:,4);
    m.leye(:,1) = c(:,5); m.leye(:,2) = c(:,6); mleye(:,3) = c(:,7);

idx.rAntenna = 8; idx.lAntenna = 11;    
    m.rAntenna(:,1) = c(:,8); m.rAntenna(:,2) = c(:,9); m.rAntenna(:,3) = c(:,10);
    m.lAntenna(:,1) = c(:,11); m.lAntenna(:,2) = c(:,12); m.lAntenna(:,3) = c(:,13);

idx.rlegTip = 14; idx.llegTip = 20;
    m.rlegTip(:,1) = c(:,14); m.rlegTip(:,2) = c(:,15); m.rlegTip(:,3) = c(:,16);
    m.llegTip(:,1) = c(:,20); m.llegTip(:,2) = c(:,21); m.llegTip(:,3) = c(:,22);
    
idx.rlegKnee = 17; idx.llegKnee = 23;
    m.rlegKnee(:,1) = c(:,17); m.rlegKnee(:,2) = c(:,18); m.rlegKnee(:,3) = c(:,19);
    m.llegKnee(:,1) = c(:,23); m.llegKnee(:,2) = c(:,24); m.llegKnee(:,3) = c(:,25);

idx.rWingtip = 26; idx.rWingbase = 29; idx.rWingcrease = 32; idx.rWingCM = 35;   
    m.rWingtip(:,1) = c(:,26); m.rWingtip(:,2) = c(:,27); m.rWingtip(:,3) = c(:,28);
    m.rWingbase(:,1) = c(:,29); m.rWingbase(:,2) = c(:,30); m.rWingbase(:,3) = c(:,31);
    m.rWingcrease(:,1) = c(:,32); m.rWingcrease(:,2) = c(:,33); m.rWingcrease(:,3) = c(:,34);
    m.rWingCM(:,1) = c(:,35); m.rWingCM(:,2) = c(:,36); m.rWingCM(:,3) = c(:,37);

idx.lWingtip = 38; idx.lWingbase = 41; idx.lWingcrease = 44; idx.lWingCM = 47;
    m.lWingtip(:,1) = c(:,38); m.lWingtip(:,2) = c(:,39); m.lWingtip(:,3) = c(:,40);
    m.lWingbase(:,1) = c(:,41); m.lWingbase(:,2) = c(:,42); m.lWingbase(:,3) = c(:,43);
    m.lWingcrease(:,1) = c(:,44); m.lWingcrease(:,2) = c(:,45); m.lWingcrease(:,3) = c(:,46);
    m.lWingCM(:,1) = c(:,47); m.lWingCM(:,2) = c(:,48); m.lWingCM(:,3) = c(:,49);

idx.Tail = 50;
    m.Tail(:,1) = c(:,50); m.Tail(:,2) = c(:,51); m.Tail(:,3) = c(:,52);
    
idx.rFootie = 53; idx.lFootie = 56;   
    m.rFootie(:,1) = c(:,53); m.rFootie(:,2) = c(:,54); m.rFootie(:,3) = c(:,55);
    m.lFootie(:,1) = c(:,56); m.lFootie(:,2) = c(:,57); m.lFootie(:,3) = c(:,58);

figure(2); clf; hold on;
    plot(m.reye(:,1),-m.reye(:,2), '.')
    plot(m.leye(:,1),-m.leye(:,2), '.')
    plot(m.rWingtip(:,1),-m.rWingtip(:,2), '.')
    plot(m.lWingtip(:,1),-m.lWingtip(:,2), '.')
    
    tailidx = find(m.Tail(:,3) > 0.5);
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

