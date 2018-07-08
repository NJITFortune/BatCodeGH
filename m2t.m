fileID = fopen('~/foobar.txt', 'w');

formatSpec = '%s %d %d %d %d %8.8f %8.8f %8.8f %8.8f %8.8f\n';

j = 1:length(bat);
kk = find(j ~= 241 & j ~= 242 & j ~= 243 & j ~= 244 & j ~= 245 & j ~= 246 & j ~= 247 & j ~= 248 & j ~= 249 & j ~= 250 & j ~= 251 & j ~= 252 & j ~= 253 & j ~= 254);
j = j(kk);

for i=1:length(j); 
    j(i)
lfreqdiff = (bat(j(i)).l.stim.freq - bat(j(i)).l.pre.freq);
rfreqdiff = (bat(j(i)).r.stim.freq - bat(j(i)).r.pre.freq);
prefreqmean = mean([bat(j(i)).l.pre.freq bat(j(i)).r.pre.freq]);

lampdiff = (bat(j(i)).l.stim.amp - bat(j(i)).l.pre.amp);
rampdiff = (bat(j(i)).r.stim.amp - bat(j(i)).r.pre.amp);
preampmean = mean([bat(j(i)).l.pre.amp bat(j(i)).r.pre.amp]);

    changewingfreq = mean([lfreqdiff rfreqdiff]);
    changewingfreqpct = changewingfreq / prefreqmean;
    
    changewingamp = mean([lampdiff rampdiff]);
    changewingamppct = changewingamp / preampmean;

    
    
 fprintf(fileID, formatSpec, bat(j(i)).ident, bat(j(i)).sndinfo(1), bat(j(i)).sndinfo(2), bat(j(i)).sndinfo(3), ...
     length(find(bat(j(i)).screamtims ~= 0)), changewingfreq, changewingfreqpct, changewingamp, changewingamppct, bat(j(i)).wingecho.wingamp*1000); 

end;

fclose(fileID);
