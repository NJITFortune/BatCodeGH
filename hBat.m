function out = hBat(sig, Fs)

% Make a time series
tim = 1/Fs:1/Fs:length(sig)/Fs;
z = zeros(1, length(sig)); 

% Make the filters
    [b,a] = butter(3, 200/(Fs/2), 'low'); % Low Pass filter for AMs
    [d,c] = butter(3, 15000/(Fs/2), 'high'); % High Pass filter for raw signal
    %[f,e] = butter(3, 30000/(Fs/2), 'high'); % High Pass filter for raw signal

% Filter the signal
hpSIG = filtfilt(d,c, sig); % High-pass filtered signal
hHP = 2*filtfilt(b,a,real(hilbert(abs(hpSIG)))); % Envelope of the signal

figure(1); clf;
    ax(1) = subplot(211); specgram(hpSIG, 256, Fs, [], 250);
        colormap('HOT'); caxis([-10 20]);
    ax(2) = subplot(212); plot(tim, [hpSIG, hHP]);
    linkaxes(ax, 'x');
    
fprintf('Click the level in top subplot. \n');    
    
[~, ylevel] = ginput(1);

z(hHP > ylevel) = 1;

    startIDXs = find(diff(z) == 1);
    endIDXs = find(diff(z) == -1);
    
    if startIDXs(1) > endIDXs(1); endIDXs = endIDXs(2:end); end
    if startIDXs(end) < endIDXs(end); startIDXs = startIDXs(1:end-1); end
    
fprintf('Found %i bat events. \n', length(startIDXs));    
    
    
% Construct the output structure: Data for each vocalization
for j = length(startIDXs):-1:1
   
    out(j).tims = [tim(startIDXs(j)), tim(endIDXs(j))];
    out(j).idxs = [startIDXs(j), endIDXs(j)];
    
    out(j).fft = fftcalculator(hpSIG(startIDXs(j):endIDXs(j)), Fs);
   
    % Make the frequency trace
    stepsize =   0.0005; 
    windowidth = 0.005;
    
    startims = out(j).tims(1) - windowidth:stepsize:out(j).tims(2);
    
    for k = length(startims):-1:1
        ttt = find(tim > startims(k) & tim < startims(k) + windowidth);
        tmp = fftcalculator(hpSIG(ttt), Fs);
        out(j).traceTIM(k) = mean(tim(ttt)) - out(j).tims(1);
            [out(j).traceAMP, idx] = max(tmp.fftdata);
        out(j).traceFREQ(k) = tmp.fftfreqs(idx);    
    end
end

% Plot random selection of up to 9 events

if length(startIDXs) > 9
    listofevents = randn(9, 


function fftout = fftcalculator(data, Fs)
% Compute the FFT (Fast Fourier Transform)
% out = fftmachine(data, Fs, smoothwindow);
% Where out is a strucutre with fftfreq and fftdata
% The smoothwindow is for a medfilt1 low-pass filtering
% of the fft data itself.  This should generally be low and
% odd, 9 or less.

data = data(isnan(data) == 0);
L = length(data);
smoothwindow = 20;

NFFT = 2^nextpow2(L); % Next power of 2 from length of the data

fftdata = fft(data, NFFT) / L;

%f = Fs/2*linspace(0,1,NFFT/2+1);

% fftdata = fft(data);

% We use only half of the data, hence fftdata(1:round(end/2));
% And we take the absolute value of the real component and filter
% that so that it is smooth

fftout.fftdata = 2*abs(fftdata(1:NFFT/2+1));
fftout.fftdata = medfilt1( fftout.fftdata, smoothwindow);

% Now we need to generate the X values - which are the frequencies

fftout.fftfreqs = Fs/2*linspace(0,1,NFFT/2+1);

% Sometimes the rounding makes it so that the lengths of the
% data and the frequency values are off by one.  Let us correct that.

minlen = min([length(fftout.fftfreqs) length(fftout.fftdata)]);
fftout.fftfreqs = fftout.fftfreqs(1:minlen);
fftout.fftdata = fftout.fftdata(1:minlen);

end

end

