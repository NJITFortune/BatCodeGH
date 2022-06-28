function EigenPairAnime2022(im, dlc, snd, idx, lim, mname)
% mothwinganime(im, data, snd, ttl, lim, mname)
% 'im' is the Video Structure 
% 'data' includes the measurements from DeepLabCut
% 'snd' is an electric channel structure from Spike2 
% 'lim' is the frame limits - for example, [800 1800]
% 'mname'  specifies an filename for movie output, 'filename.avi'

    trailframeno = 6;
    timwin = 1;

    numFrames = im.NumFrames;
    vFs = 10; % 10 frames per second (NSB2022)
    vtim = 1/vFs:1/vFs:numFrames/vFs;

    writerObj = VideoWriter(mname);
    writerObj.FrameRate = 10;
    open(writerObj);

% Some constants    
    Fs = snd(idx).Fs;
    tim = 1/Fs:1/Fs:length(snd(idx).data)/Fs;
    %env = lpf(abs(real(hilbert(snd.values))), Fs, [500 5]);
    velims = [min(dlc.relVel), max(dlc.relVel)];
    distlims = [min(dlc.reldistance(lim(1):lim(2))), max(dlc.reldistance(lim(1):lim(2)))];
    fishvelims = [-10, max([max(dlc.vel1(lim(1):lim(2))), max(dlc.vel2(lim(1):lim(2)))])];

% Use figure 1
figure(1);

%% The main loop

for i = lim(1)+8:lim(2)
    
    clf;

    curtim = i/vFs;

    % Our image
    axes('position', [0.1875 0.1500 0.7750 0.8150]); cla;
    tmpFrame = read(im, i, "native");
    imshow(tmpFrame); %freezeColors;

    hold on;
    % Plot nose contrails
            plot(dlc.fish1(idx).nose(i-trailframeno:i,1),dlc.fish1(idx).nose(i-trailframeno:i,2), 'm.-', 'LineWidth', 1, 'MarkerSize', 4);
            plot(dlc.fish2(idx).nose(i-trailframeno:i,1),dlc.fish2(idx).nose(i-trailframeno:i,2), 'b.-', 'LineWidth', 1, 'MarkerSize', 4);
    hold off;

% Velocity plot    
    axes('position', [0.05 0.7 0.25 0.15]); cla;

    tvt = find(vtim > curtim - timwin & vtim < curtim + timwin);
        tvt(end+1) = tvt(end)+1;
    hold on;
    plot(vtim(tvt) - vtim(tvt(1)), dlc.vel1(tvt), 'm-');
    plot(vtim(tvt) - vtim(tvt(1)), dlc.vel2(tvt), 'b-');
    xlim([0 timwin*2])
    ylim(fishvelims)
    hold on; plot([timwin timwin], distlims, 'k-', 'LineWidth', 1);
    title('Fish Velocity')
    set(gca, 'xtick', [], 'ytick', []);
    
    % Vertical line
    
    plot([8 8],[600 -800], 'k-');
hold off;

% Relative distance plot
    axes('position', [0.05 0.5 0.25 0.15]); cla;
    plot(vtim(tvt) - vtim(tvt(1)), dlc.reldistance(tvt), 'LineWidth', 4);
    hold on; plot([timwin timwin], distlims, 'k-', 'LineWidth', 1);
    ylim(distlims); xlim([0 timwin*2])
    title('Distance between fish')
    set(gca, 'xtick', [], 'ytick', []);
    
% Specgram plot
    tt = find(tim > curtim - timwin & tim < curtim + timwin);
    axes('position', [0.05 0.3 0.25 0.15]); cla;
    specgram(snd(idx).data(tt), 1024*8, Fs, [], floor(1024*8*0.9));
    ylim([535 675]);
    colormap('HOT'); 
    caxis([0 50]); 
    hold on; plot([timwin timwin],[535 675], 'w-', 'LineWidth', 2);
    title('EODs')
    set(gca, 'xlabel', [], 'ylabel', [],'xtick', [], 'ytick', []);

% capture move for movie, if requested 


    out= getframe(gcf);
    writeVideo(writerObj,out);    


    pause(0.0001);
    %close(1);
    
end


  close(writerObj);
