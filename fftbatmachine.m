function out = fftbatmachine(data, Fs, smoothwindow)
% Compute the FFT (Fast Fourier Transform)
% out = fftmachine(data, Fs, smoothwindow);
% Where out is a strucutre with fftfreq and fftdata
% The smoothwindow is for a medfilt1 low-pass filtering
% of the fft data itself.  This should generally be low and
% odd, 9 or less.

if nargin < 2
	fprintf('Usage: fftmachine(data, SampleRate, [window])');
	exit 0
end

data = data(isnan(data) == 0);

L = length(data);

NFFT = 2^nextpow2(L); % Next power of 2 from length of the data

fftdata = fft(data,NFFT)/L;

%f = Fs/2*linspace(0,1,NFFT/2+1);

% fftdata = fft(data);

% We use only half of the data, hence fftdata(1:round(end/2));
% And we take the absolute value of the real component and filter
% that so that it is smooth

out.fftdata = 2*abs(fftdata(1:NFFT/2+1));
% out.fftdata = abs(real(fftdata(1:round(end/2))));

if nargin == 3;
	out.fftdata = medfilt1( out.fftdata, smoothwindow);
end;

% Now we need to generate the X values - which are the frequencies

%stepsize = Fs/round(length(data));
%out.fftfreq = stepsize:stepsize:Fs/2;
out.fftfreq = Fs/2*linspace(0,1,NFFT/2+1);

% Sometimes the rounding makes it so that the lengths of the
% data and the frequency values are off by one.  Let us correct that.

minlen = min([length(out.fftfreq) length(out.fftdata)]);
out.fftfreq = out.fftfreq(1:minlen);
out.fftdata = out.fftdata(1:minlen);

% Finally, we can plot

%plot(out.fftfreq, out.fftdata);

