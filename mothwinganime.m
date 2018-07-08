function mothwinganime(im, data, snd, ttl, lim, offset, mname)
% mothwinganime(im, data, snd, ttl, lim, mname)
% 'im' is the image structure from frontwingtrack
% 'data' includes the measurements from frontwingtrack
% 'snd' is the sound channel structure from Spike2 
% 'lim' is the frame limits - for example, [800 1800]
% 'mname' is optional. This specifies an filename for movie output, 'filename.avi'

%% Preparations
% Only make a movie if the user asks for it
if nargin > 6;
    writerObj = VideoWriter(mname);
    writerObj.FrameRate = 30;
    open(writerObj);
end;

% Some constants    
Fs = 1/snd.interval;
tim = 1/Fs:1/Fs:snd.length/Fs;
env = lpf(abs(real(hilbert(snd.values))), Fs, [500 5]);

% Use figure 1
figure(1); 

%% The main loop
for i = lim(1)+8:lim(2);

    ii = i - offset;
    
    clf;
    
    % Our image
    axes('position', [0.1875 0.1500 0.7750 0.8150]); cla;
    imshow(im(ii).orig); freezeColors;

    hold on;
    % Loop for the 7 previous data points for the wing contrails
    for j = 0:7;
        %if ~isempty(find(data.right.vel(i) < 0, 1)); % Downstrokes
            plot(data.right.X2(ii-(j+1):ii-j),data.right.Y2(ii-(j+1):ii-j), 'r*-', 'LineWidth', 2);
        %else % Upstrokes
            plot(data.right.X2(ii-(j+1):ii-j),data.right.Y2(ii-(j+1):ii-j), 'r*-', 'LineWidth', 1);
        %end;
        
        %if ~isempty(find(data.left.vel(i) < 0, 1)); % Downstrokes
            plot(data.left.X2(ii-(j+1):ii-j),data.left.Y2(ii-(j+1):ii-j), 'g*-', 'LineWidth', 2);
        %else % Upstrokes
            plot(data.left.X2(ii-(j+1):ii-j),data.left.Y2(ii-(j+1):ii-j), 'g*-', 'LineWidth', 1);
        %end;    
    end;
    hold off;

% Position / Velocity plot    
    axes('position', [0.05 0.7 0.25 0.15]); cla;

    % position
    plot(1:36, -data.right.Y(ii-7:ii+28) -200, 'r', 1:36, -data.left.Y(ii-7:ii+28) -200, 'g'); 
    xlim([1 36]);    % was 22, 14
    ylim([-800 600]);
hold on;
    % velocity
    plot(1:36, -data.right.vel(ii-7:ii+28) +300, 'r', 1:36, -data.left.vel(ii-7:ii+28) +300, 'g');

    text(20,-500, 'Position');
    text(20,100, 'Velocity');
    set(gca, 'xtick', [], 'ytick', []);
    
    % Vertical line
    
    plot([8 8],[600 -800], 'k-');
hold off;

% Envelope plot
    axes('position', [0.05 0.5 0.25 0.15]); cla;
    tt = find(tim > ttl.times(i-7) & tim < ttl.times(i+28));
    plot(tim(tt), env(tt));
    xlim([tim(tt(1)) tim(tt(end))]);
    ylim([-0.015 0.1]);
    text(tim(tt(1))+0.017, 0.01, 'AM');
    timlin = tim(tt(1)) + (tim(tt(end)) - tim(tt(1)))/5;
    hold on; plot([timlin timlin],[-0.015 0.1], 'k-');    
    set(gca, 'xtick', [], 'ytick', []);
    
    
% Specgram plot
    axes('position', [0.05 0.3 0.25 0.15]); cla;
    specgram(snd.values(tt), 256, Fs, [], 250);
    ylim([10000 60000]);
    colormap('HOT'); 
    caxis([-40 10]); 
    set(gca, 'xlabel', [], 'ylabel', [],'xtick', [], 'ytick', []);
    hold on; plot([timlin-tim(tt(1)) timlin-tim(tt(1))],[11000 59000], 'w-', 'LineWidth', 2);

% capture move for movie, if requested 

if nargin > 6;
    out= getframe(gcf);
    writeVideo(writerObj,out);    
end;

    pause(0.0001);
    %close(1);
    
end;

if nargin > 6;
    close(writerObj);
end;
