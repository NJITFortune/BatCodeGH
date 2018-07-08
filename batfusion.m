function out = batfusion(vidata, pulsetims, snd, Fs, ident, nums, lims)
% Usage: out = batfusion(vidata, pulsetims, snd, Fs, ident, num, lims)
% num is [freq amplitude position]
% freq is in kHz; vol 3 or 5; pos 1 top, 2 bottom, 3 left; 


%% Preparations

    if nargin < 7; lims = [0 0]; end; % If the user sets the lims use those, otherwise we use ginput later.

    startim = pulsetims(vidata.left.imnum(1)); % The video pulses are a subset of the audio recording, so we use those as the start and end times.
    endtim = pulsetims(vidata.left.imnum(end));

    pt = pulsetims(vidata.left.imnum(2:end-1)); % I don't know why I use this - skips first and last pulses.

    tim = 1/Fs:1/Fs:length(snd)/Fs; % The time base for the sound
    tt = find(tim > startim & tim < endtim); % The indices for the sound in relation to the pulses.

    vFs = 1/mean((pulsetims(2:end) - pulsetims(1:end-1))); % Average duration between video pulses
    % vtim = 1/vFs:1/vFs:length(vidata.left.Y)/vFs; % Time base for the video

    out.screamtims = 0;
    
%% Initial Calculations

    % High-pass filter the sound
    [b,a] = butter(3, 20000/Fs, 'high');
    hpsnd = filtfilt(b,a,snd);
    % Low-pass filter the wing movements
    [d,c] = butter(3, 100/vFs, 'low');
    leftfilt = filtfilt(d,c,(vidata.left.Y - (mean(vidata.left.Y))));
    rightfilt = filtfilt(d,c,(vidata.right.Y - (mean(vidata.right.Y))));

%% Get the zero crossings

left = zeros(1, length(leftfilt)); right = zeros(1, length(rightfilt));

% Usual trick with diff of a square wave from original signal
left(find(leftfilt < 0)) = 0; left(find(leftfilt > 0)) = 1; left = diff(left); leftidx = find(left == 1);
if leftidx(1) == 1; leftidx(1) = 2; end;
    %figure(3); subplot(211); plot(vtim, leftfilt,'k'); hold on; plot(vtim(leftidx), zeros(1,length(leftidx)), 'g*');
right(find(rightfilt < 0)) = 0; right(find(rightfilt > 0)) = 1;  right = diff(right); rightidx = find(right == 1);
    %figure(3); subplot(212); plot(vtim, rightfilt,'k'); hold on; plot(vtim(rightidx), zeros(1,length(rightidx)), 'g*');
if rightidx(1) == 1; rightidx(1) = 2; end;

% Get a more precise estimate of the zero crossing time by adding fractions of samples

left_amp_step = leftfilt(leftidx) - leftfilt(leftidx-1); % amplitude step between samples just below and just above zero
left_amp_frac = leftfilt(leftidx) ./ left_amp_step;    % time fraction of the amplitude step that is below zero
left_frac = leftidx - left_amp_frac;                % add fraction to time stamp just below zero crossing to get the "real" zero crossing
left_times = cumsum(diff(left_frac)/vFs) + startim; % Make a new time base
    %figure(3); subplot(211); plot(left_times - pulsetims(vidata.left.imnum(1)), zeros(1,length(left_times)), 'm*');
left_durs = left_times(2:end) - left_times(1:end-1); % Get the durations...
left_freqs = 1 ./ (left_times(2:end) - left_times(1:end-1)); % ...which gives us the instantaneous frequency

right_amp_step = rightfilt(rightidx) - rightfilt(rightidx-1);   
right_amp_frac = rightfilt(rightidx) ./ right_amp_step;   
right_frac = rightidx - right_amp_frac;                
right_times = cumsum(diff(right_frac)/vFs) + startim;
    %figure(3); subplot(212); plot(right_times - pulsetims(vidata.right.imnum(1)), zeros(1,length(right_times)), 'm*');
right_durs = right_times(2:end) - right_times(1:end-1);
right_freqs = 1 ./ (right_times(2:end) - right_times(1:end-1));



%% Get the amplitudes for each cycle


for kk = length(left_times):-1:2;
    tmptimidx = find(pulsetims(vidata.left.imnum) > left_times(kk-1) & pulsetims(vidata.left.imnum) < left_times(kk));
    leftamplitude(kk-1) = max(leftfilt(tmptimidx)) - min(leftfilt(tmptimidx));
end;
for kk = length(right_times):-1:2;
    tmptimidx = find(pulsetims(vidata.right.imnum) > right_times(kk-1) & pulsetims(vidata.right.imnum) < right_times(kk));
    rightamplitude(kk-1) = max(rightfilt(tmptimidx)) - min(rightfilt(tmptimidx));
end;

% deend = min([length(leftidx) length(rightidx)]); % These should be the same length... but it failed here for some reason
% for kk = deend:-1:2;
%     leftamplitude(kk-1) = max(leftfilt(leftidx(kk-1):leftidx(kk))) - min(leftfilt(leftidx(kk-1):leftidx(kk)))
%     rightamplitude(kk-1) = max(rightfilt(rightidx(kk-1):rightidx(kk))) - min(rightfilt(rightidx(kk-1):rightidx(kk)));
% end;

%% Get the mean wingbeat frequency for the whole trial

    wingfreq = mean([right_freqs left_freqs]);    

%% Plotting

ident

figure(1); clf;

figure(1); ax(1) = subplot(511); 
        plot(pt-pt(1), medfilt1(vidata.left.vel,3), 'g', pt-pt(1), medfilt1(vidata.right.vel,3), 'r'); xlim([0 pt(end)-pt(1)]); ylabel('Velocity'); % Velocity
figure(1); ax(2) = subplot(512); 
        plot(left_times(2:end) - pt(1), left_freqs, 'g', right_times(2:end) - pt(1), right_freqs, 'r'); xlim([0 pt(end)-pt(1)]); ylim([5 35]); ylabel('Inst. Wingbeat Freq.'); % Wingbeat Freq
figure(1); ax(3) = subplot(513); 
        plot(pt-pt(1), vidata.left.Y(1:end-1), 'g', pt-pt(1), vidata.right.Y(1:end-1), 'r'); ylabel('Position'); % Wing Y
figure(1); ax(4) = subplot(514); 
    plot(left_times(2:end)-pt(1), leftamplitude, 'g', right_times(2:end) - pt(1), rightamplitude, 'r'); ylabel('Amplitude'); % Wingbeat Amp
figure(1); ax(5) = subplot(515); 
        specgram(hpsnd(tt),256,Fs); ylim([0000 150000]); colormap('HOT'); caxis([-40 -20]); hold on; % Specgram of microphone
        if lims(2) == 0; 
            fprintf('Click start and end of stimulus. \n');
            [lims(1), ~] = ginput(1); plot([lims(1) lims(1)], [10000 140000], 'g', 'LineWidth', 2); 
            [lims(2), ~] = ginput(1); plot([lims(2) lims(2)], [10000 140000], 'm', 'LineWidth', 2);
        end;
        if lims(2) ~= 0; 
            plot([lims(1) lims(1)], [10000 140000], 'g', 'LineWidth', 2); 
            plot([lims(2) lims(2)], [10000 140000], 'm', 'LineWidth', 2);
        end;
        
linkaxes(ax,'x');

                %hold on;
                %plot([lims(1) lims(1)], [0 max(smEnv)], 'g', 'LineWidth', 2); 
                %plot([lims(2) lims(2)], [0 max(smEnv)], 'm', 'LineWidth', 2);

% figure(2); subplot(212); plot(pow.fftfreq(ff), pow.fftdata(ff)); xlim([10 50]);


% Leftover figure crapola 
% figure(1); ax(2) = subplot(512); plot(left_times(2:end), left_freqs, 'g', right_times(2:end), right_freqs, 'r'); % Wingbeat Freq
% figure(1); ax(3) = subplot(313); plot(tim(tt), hpsnd(tt)); % Oscillogram
% figure(1); ax(1) = subplot(311); plot(left_times(2:end), left_durs, 'g*-', right_times(2:end), right_durs, 'r*-');
% linkaxes(axa,'x');
% figure(2); subplot(211); plot(tim(tt), env, 'k'); xlim([tim(tt(1)) tim(tt(end))]); 

%% Output analysis

% Left amp and freq data

pretim = find((left_times-pt(1)) < lims(1));
stimtim = find((left_times-pt(1)) > lims(1) & (left_times-pt(1)) < lims(2));
postim = find((left_times-pt(1)) > lims(2));

out.l.pre.amp = nanmean(leftamplitude(pretim));
out.l.pre.freq = nanmean(left_freqs(pretim));
out.l.stim.amp = nanmean(leftamplitude(stimtim));
out.l.stim.freq = nanmean(left_freqs(stimtim));
out.l.post.amp = nanmean(leftamplitude(postim(1):end));
out.l.post.freq = nanmean(left_freqs(postim(1):end));

pretim = find((right_times-pt(1)) < lims(1));
stimtim = find((right_times-pt(1)) > lims(1) & (right_times-pt(1)) < lims(2));
postim = find((right_times-pt(1)) > lims(2));

% Right amp and freq data
out.r.pre.amp = nanmean(rightamplitude(pretim));
out.r.pre.freq = nanmean(right_freqs(pretim));
out.r.stim.amp = nanmean(rightamplitude(stimtim));
out.r.stim.freq = nanmean(right_freqs(stimtim));
out.r.post.amp = nanmean(rightamplitude(postim(1):end));
out.r.post.freq = nanmean(right_freqs(postim(1):end));

%% Wing echo data

% Extract the amplitude envelope of the microphone recording


    [q,w] = butter(3, 200/Fs, 'low');

    reductionfactor = 10;
    fullEnv = filtfilt(q,w,abs(hilbert(hpsnd(tt))));
    fullsmEnv = fullEnv(1:reductionfactor:end); smFs = Fs/reductionfactor; 
    fullsmTim = 1/smFs:1/smFs:length(fullsmEnv)/smFs; % Downsample the envelope
    
    protectiontime = 0.1; % This is extra time before and after the stimulus to make sure that we got the whole thing.
    ttt = find(tim > lims(1)+pt(1)-protectiontime & tim < lims(2)+pt(1)+protectiontime);
    env = filtfilt(q,w,abs(hilbert(hpsnd(ttt)))); % Take the low-pass filtered Hilbert of the original sound chunk (the envelope)
    smEnv = env(1:reductionfactor:end); smFs = Fs/reductionfactor; 
    smTim = 1/smFs:1/smFs:length(smEnv)/smFs; % Downsample the envelope

    % figure(10); plot(smTim, smEnv);
    
    pow = fftbatmachine(smEnv, smFs); % Perform an FFT on the wingbeat focussed on when the sound played
    wop = fftbatmachine(fullsmEnv, smFs); % FFT using the entire dataset.  This is because the above is sensitive to the edges.
    
    figure(2); clf; 
        subplot(311); plot(fullsmTim, fullsmEnv, 'k'); xlim([fullsmTim(1) fullsmTim(end)]); xlabel('Time(sec)'); ylabel('Amplitude');
        hold on;
        plot([lims(1) lims(1)], [0 max(fullsmEnv)], 'g', 'LineWidth', 2); 
        plot([lims(2) lims(2)], [0 max(fullsmEnv)], 'm', 'LineWidth', 2);

    % Calculate the 1/f distribution
    uu = find(pow.fftfreq >= 1);
        uulen = length(uu);
        aFreq = pow.fftfreq(uu(1)):max(pow.fftfreq)/(uulen):max(pow.fftfreq);
        aMag = 1 ./ aFreq;
        aMag = aMag * max(pow.fftdata(uu));

        %uulen
        length(aMag);
        length(pow.fftfreq(uu));

        fftlen = min([length(aFreq) uulen]);
        
    figure(2); subplot(312); 
    % Plot the power spectrum with 1/f power law
    %        semilogy(pow.fftfreq(uu), pow.fftdata(uu), 'k'); hold on; semilogy(aFreq, aMag, 'm'); xlim([0 100]); xlabel('AM Frequency Hz'); ylabel('Power');

    % Plot the raw microphone trace
            plot(tim(tt), hpsnd(tt)); xlim([tim(tt(1)) tim(tt(end))]);   

    % Pick the peak power at the wingbeat frequency from the video
    noFnoiseDat = pow.fftdata(1:fftlen) - aMag(1:fftlen)'; 
    wingfreqidx = find(aFreq >= wingfreq-2 & aFreq <= wingfreq+2); % Look around the wingfreq
    [peakwingamp, peakwingidx] = max(noFnoiseDat(wingfreqidx));

    % Pick the peak power in the range of possible frequencies
    newff = find(aFreq >= 15 & aFreq <= 45);
    [peakamp, peakidx] = max(noFnoiseDat(newff));
    
    figure(2); subplot(313); 
        plot(aFreq, noFnoiseDat); xlim([15 45]);
    
    %[peakamp, peakidx] = max(noFnoiseDat(newff)); 
    figure(2); subplot(313);
        hold on; 
        plot(aFreq(newff(peakidx)), peakamp, 'r*'); 
        plot(aFreq(wingfreqidx(peakwingidx)), peakwingamp, 'go', 'Markersize', 10); xlabel('AM Frequency Hz'); ylabel('Power');

    out.vidwingfreq = wingfreq;

    out.wingecho.peakamp = peakamp;
    out.wingecho.peakfreq = aFreq(newff(peakidx));
    out.wingecho.wingamp = peakwingamp;
    out.wingecho.wingfreq = aFreq(wingfreqidx(peakwingidx));
    
    out.wingecho.peakdB = real(20 * log(peakamp/mean(noFnoiseDat(newff))));
    out.wingecho.wingdB = real(20 * log(peakwingamp/mean(noFnoiseDat(newff))));
    
    out.wingecho.fft = pow;
    out.wingecho.hilbert.env = fullsmEnv;
    out.wingecho.hilbert.tim = fullsmTim;

%% Get both chirp times
yn = 0;
yn = input('Are there any moth calls? 1 for yes, enter or 0 for no. ');

if yn == 1;
    figure(3); specgram(hpsnd(tt),1024,Fs); ylim([0000 100000]); colormap('HOT'); % Specgram of microphone
    caxis([-40 -20]); hold on; plot([lims(1) lims(1)], [10000 90000], 'c-'); plot([lims(2) lims(2)], [10000 90000], 'g-');
    [chirptims, ~] = ginput;
    ys = ones(1,length(chirptims)*2);
    for i=1:length(chirptims); 
        plot([chirptims(i) chirptims(i)], [10000 90000], 'g-', 'LineWidth', 2);
    end;   
    out.screamtims = chirptims - lims(1);
end;

out.ident = ident;
out.sndinfo = nums;
out.lims = lims;

 
if yn == 1; pause(0.5); close(3); end;



