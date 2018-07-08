function [spss, names] = analanal(fdat)
% Usage: out = analanal(fdat)

% Prepare the figures
figure(1); clf; subplot(211); hold on; subplot(212); hold on;
figure(2); clf; subplot(211); hold on; subplot(212); hold on;
figure(3); clf; subplot(121); hold on; subplot(122); hold on;

% Set some colors
colrl(1,:)='r-'; colrl(2,:)='m-'; colrl(3,:)='b-'; colrl(4,:)='g-'; 
colrl(5,:)='c-'; colrl(6,:)='k-'; colrl(7,:)='y-'; 
colrs(1,:)='r*'; colrs(2,:)='m*'; colrs(3,:)='b*'; colrs(4,:)='g*'; 
colrs(5,:)='c*'; colrs(6,:)='k*'; colrs(7,:)='y*'; 

% Extract the sound and microphone properties
for i = 1:length(fdat)-1;
    if isempty(fdat(i).sndinfo) ~= 1;
        sndfreq(i) = fdat(i).sndinfo(1);
        sndamp(i) = fdat(i).sndinfo(2);
        sndpos(i) = fdat(i).sndinfo(3);
    end;
end;

idx.pos1idx = []; idx.pos2idx = []; idx.pos3idx = [];
idx.amp3idx = []; idx.amp5idx = [];
z = 0;

%% Giant loop with the two amps and the three positions - calculate for all freqs
for amp = 3:2:5;

    for pos = 1:3;

    freqs = 10:5:100;

    for jj = 1:length(freqs);
    
        kk = find(sndfreq == freqs(jj) & sndamp == amp & sndpos == pos);
        foo(jj).freq = freqs(jj);

    if isempty(kk) ~= 1;
        for i = 1:length(kk);
            foo(jj).lprefreq(i) = fdat(kk(i)).l.pre.freq;
            foo(jj).lstimfreq(i) = fdat(kk(i)).l.stim.freq;
            foo(jj).lpostfreq(i) = fdat(kk(i)).l.post.freq;
    
            foo(jj).lpreamp(i) = fdat(kk(i)).l.pre.amp;
            foo(jj).lstimamp(i) = fdat(kk(i)).l.stim.amp;
            foo(jj).lpostamp(i) = fdat(kk(i)).l.post.amp;
    
    
            foo(jj).wingdB(i) = fdat(kk(i)).wingecho.wingdB;
            foo(jj).wingamp(i) = fdat(kk(i)).wingecho.wingamp;
    
            foo(jj).screamnum(i) = length(find(fdat(kk(i)).screamtims > 0));
    
            z = z+1;
            
            if amp == 3; idx.amp3idx(end+1) = z; end;
            if amp == 5; idx.amp5idx(end+1) = z; end;
            if pos == 1; idx.pos1idx(end+1) = z; end;
            if pos == 2; idx.pos2idx(end+1) = z; end;
            if pos == 3; idx.pos3idx(end+1) = z; end;

            %out(z).wingamp = foo(jj).wingamp(i);
            out(z).deltamp = foo(jj).lstimamp(i) - foo(jj).lpreamp(i);
            out(z).deltapercentamp = (foo(jj).lstimamp(i) - foo(jj).lpreamp(i)) / foo(jj).lpreamp(i);
            
            out(z).deltapercentfreq = (foo(jj).lstimfreq(i) - foo(jj).lprefreq(i)) / foo(jj).lprefreq(i);
            out(z).deltafreq = (foo(jj).lstimfreq(i) - foo(jj).lprefreq(i));
            
            out(z).screamnum = foo(jj).screamnum(i);
            
            spss(z,:) = [fdat(kk(i)).sndinfo, foo(jj).screamnum(i), out(z).deltafreq, out(z).deltapercentfreq, out(z).deltamp, out(z).deltapercentamp, foo(jj).wingamp(i)];
            names{z} = fdat(kk(i)).ident;
            
            
            figure(1); subplot(211); plot(foo(jj).freq, foo(jj).lstimfreq(i) - foo(jj).lprefreq(i), colrs(amp-2,:));
            figure(1); subplot(212); plot(foo(jj).freq, foo(jj).lstimamp(i) - foo(jj).lpreamp(i), colrs(amp-2,:));
        
            figure(2); subplot(211); plot(foo(jj).freq, foo(jj).wingdB(i), colrs(pos,:));
            figure(2); subplot(212); plot(foo(jj).freq, foo(jj).wingamp(i), colrs(pos,:));

            figure(3); subplot(121); plot(foo(jj).freq, foo(jj).screamnum(i), colrs(amp-2,:));        
            figure(3); subplot(122); plot(foo(jj).freq, foo(jj).screamnum(i), colrs(pos,:));

        end;
        
%%%% ASCII Code
% save -ascii maledata.txt out
% fileID = fopen('malenames.dat','w');
% formatSpec = '%s\n';
% for i=1:234; fprintf(fileID, formatSpec, idx{i}); end;
% fclose(fileID)
        
 fileID = fopen('tmp.dat','w');
 formatSpec = '%s %d %d %d %d %2.1f %2.1f %2.5f %2.5f %2.5f \n';
 for i=1:max(z); fprintf(fileID, formatSpec, fdat(i).ident, spss(i,1), spss(i,2),spss(i,3),spss(i,4),spss(i,5),spss(i,6),spss(i,7),spss(i,8)*10,spss(i,9)*1000); end;
 fclose(fileID);

%     chngfreq = foo(jj).lstimfreq - foo(jj).lprefreq;
%     chngamp = foo(jj).lstimamp - foo(jj).lpreamp;
% 
%     hAmp(jj) = ttest(chngamp);
%     hFreq(jj) = ttest(chngfreq);
%     hEcho(jj) = ttest(foo(jj).wingamp);

    end;
    
    %figure(1); subplot(211); plot(foo(jj).freq-1, mean(chngfreq), 'c*', 'MarkerSize', 8);
    %figure(1); subplot(212); plot(foo(jj).freq-1, mean(chngamp), 'c*', 'MarkerSize', 8);

end;

end;
end;


% for i = 1:length(mdat);
%     Mnum(i) = mdat(i).num;
% end;
% 
% freqs = 10:5:100;
% 
% figure(4); clf; subplot(211); hold on; subplot(212); hold on;
% figure(5); clf; subplot(211); hold on; subplot(212); hold on;
% figure(6); clf; hold on;
% 
% for jj = 1:length(freqs);
%     
%     kk = find(Mnum == freqs(jj));
%     oof(jj).freq = freqs(jj);
% 
% for i = 1:length(kk);
%     oof(jj).lprefreq(i) = mdat(kk(i)).l.pre.freq;
%     oof(jj).lstimfreq(i) = mdat(kk(i)).l.stim.freq;
%     oof(jj).lpostfreq(i) = mdat(kk(i)).l.post.freq;
%     
%     oof(jj).lpreamp(i) = mdat(kk(i)).l.pre.amp;
%     oof(jj).lstimamp(i) = mdat(kk(i)).l.stim.amp;
%     oof(jj).lpostamp(i) = mdat(kk(i)).l.post.amp;
%     
%     oof(jj).wingdB(i) = mdat(kk(i)).wingecho.wingdB;
%     oof(jj).wingamp(i) = mdat(kk(i)).wingecho.wingamp;
%     
%     oof(jj).screamnum(i) = length(find(mdat(kk(i)).screamtims > 0));
% 
%     
%         figure(4); subplot(211); plot(oof(jj).freq, oof(jj).lstimfreq(i) - oof(jj).lprefreq(i), 'bo');
%         figure(4); subplot(212); plot(oof(jj).freq, oof(jj).lstimamp(i) - oof(jj).lpreamp(i), 'bo');
%         
%         figure(5); subplot(211); plot(oof(jj).freq, oof(jj).wingdB(i), 'bo');
%         figure(5); subplot(212); plot(oof(jj).freq, oof(jj).wingamp(i), 'bo');
% 
%         figure(6); plot(oof(jj).freq, oof(jj).screamnum(i), 'bo');
%                 
% end;
% 
%     chngfreq = oof(jj).lstimfreq - oof(jj).lprefreq;
%     chngamp = oof(jj).lstimamp - oof(jj).lpreamp;
% 
%     mhAmp(jj) = ttest(chngamp);
%     mhFreq(jj) = ttest(chngfreq);
%     mhEcho(jj) = ttest(oof(jj).wingamp);
% 
%     figure(4); subplot(211); plot(oof(jj).freq, mean(chngfreq), 'r*', 'MarkerSize', 8);
%     figure(4); subplot(212); plot(oof(jj).freq, mean(chngamp), 'r*', 'MarkerSize', 8);
% 
% end;

