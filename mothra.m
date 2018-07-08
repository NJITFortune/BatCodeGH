function out = mothra(snd, ttl)
% out = mothra(snd, ttl)
% Where snd is the structure for the channel with the audio recording
% and ttl is the structure for the frame timings

Fs = 1/snd.interval; % Get the sample frequency
tim = 1/Fs:1/Fs:snd.length/Fs;

fsnd = wrenfilter(snd.values,Fs,[10000 120000 3]);

startim = ttl.times(1);

figure(1); 

ax(1) = subplot(312);
    plot(tim, fsnd);   

ax(2) = subplot(311);
    specgram(fsnd, 512, Fs, [], 320);
    caxis([-30 10]);
    colormap('HOT');
    ylim([10000 80000]);
    

xlim([startim+2 startim+5]);

% Plot little red lines where every frame was shot.
subplot(312);
    hold on;
    tops = max(abs(fsnd))*0.1;
    for j=1:ttl.length; plot([ttl.times(j) ttl.times(j)], [tops -tops], 'r-'); end;
 
hlsnd = lpf(abs(real(hilbert(fsnd))), Fs, [500, 5]);
hhsnd = hpf(abs(real(hilbert(fsnd))), Fs, [40000, 5]);
hhsnd = lpf(hhsnd, Fs, [65000, 5]);

ax(3) = subplot(313);
    plot(tim, hlsnd + 0.06, 'g-'); 
    hold on;
    plot(tim, hhsnd, 'm-');

linkaxes(ax,'x');

out.hsnd = hlsnd;
out.hhsnd = hhsnd;
out.Fs = Fs;

numclicks = input('How many clicks do you want? ');

if numclicks > 0 ;
    [tt, ~] = ginput(numclicks);
    
    for k = 1:numclicks
        
        tmp = find(ttl.times < tt(k)); 
        frameno(k) = tmp(end);
    end;
end;

frameno

%[x, ~] = ginput(2);
%    xlim([x(1) x(2)]);
    

