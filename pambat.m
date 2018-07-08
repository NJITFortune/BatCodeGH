%% Read the files from the images
flist = dir('*.jpg');

%% Convert the image to black and white / yellow and pink
% for ff = 2:length(flist);
strt = 1;
for ff = strt:250;

thresh = 0.02; % Threshold level

     im(ff).orig = imread(flist(ff).name); 
     
     bw = im2bw(im(ff).orig, thresh); % Make the BW image
 
       for y=1:3; imyg(:,:,y) = double(bw); end;
%       imyg(:,:,3) = ones(600,800);                   
%       imyg(:,:,2) = zeros(600,800);                  
       imyg(:,:,3) = ones(480,640);                   
       imyg(:,:,2) = zeros(480,640);                  

% figure(1); clf;
% 
% for i = 1:length(im);
%     subplot(5,5,i);
%     imshow(im(i).bw);
% end;        

%% Get clicks

% Start by clicking just 17 frames
% rclks = [];
% lclks = [];
% 
% for j = 1:17;
%     figure(1);    
%     figprop = get(gcf,'Position'); 
%     set(gcf,'Position',[figprop(1) figprop(2) 1000 400]);    
% 
%     subplot(121);
%     imshow(im(j).orig);
% 
%     if isempty(rclks) == 0;
%         hold on;
%         plot(rclks(:,1), rclks(:,2), 'r-*');
%         plot(lclks(:,1), lclks(:,2), 'g-*');
%         subplot(122);
%         imshow(im(j).yg);
%         tmp = ginput(1);
%         rclks = [rclks; tmp];        
%         tmp = ginput(1);
%         lclks = [lclks; tmp];
%     else
%         subplot(122); 
%         imshow(im(j).yg);
%         rclks = ginput(1);
%         lclks = ginput(1);
%     end;
% end;    

%% Auto wing track
if ff == strt;

figure(3);

imshow(imyg);

% Step one, get medial limits

[rlim, rhinge] = ginput(1);
    hold on; 
    plot([rlim, rlim], [1 599], 'w-', 'LineWidth', 2);
    plot([1 rlim], [rhinge rhinge], 'y-', 'LineWidth', 2);
[llim, lhinge] = ginput(1);
    hold on; 
    plot([llim, llim], [1 599], 'w-', 'LineWidth', 2);
    plot([llim, 799], [lhinge, lhinge], 'y-', 'LineWidth', 2);
    
    vlim = 565; % Ventral limit is to prevent problems with numeric labels

% Step two, get the leg/antenna limits

    [rAnt, ~] = ginput(1);
    plot([rAnt rAnt], [1 599], 'y-', 'LineWidth', 4);
    [lAnt, ~] = ginput(1);
    plot([lAnt lAnt], [1 599], 'y-', 'LineWidth', 4);


% Step three, get the current wing position    

    [Xrc, Yrc] = ginput(1);
        plot(Xrc, Yrc, 'r*');    
    [Xlc, Ylc] = ginput(1);
        plot(Xlc,Ylc, 'g*');
    
    boxlim = [90 120];
    cnt = 5;
    
    plot([Xrc-boxlim(1) Xrc-boxlim(1)], [Yrc-boxlim(2) Yrc+boxlim(2)], 'w-', 'LineWidth', 2);
    plot([Xlc+boxlim(1) Xlc+boxlim(1)], [Ylc-boxlim(2) Ylc+boxlim(2)], 'w-', 'LineWidth', 2);
    
    plot([Xrc-boxlim(1) rlim], [Yrc-boxlim(2) Yrc-boxlim(2)], 'w-', 'LineWidth', 2);
    plot([Xrc-boxlim(1) rlim], [Yrc+boxlim(2) Yrc+boxlim(2)], 'w-', 'LineWidth', 2);
    
    plot([llim Xlc+boxlim(1)], [Ylc-boxlim(2) Ylc-boxlim(2)], 'w-', 'LineWidth', 2);
    plot([llim Xlc+boxlim(1)], [Ylc+boxlim(2) Ylc+boxlim(2)], 'w-', 'LineWidth', 2);

    rbox = [Xrc-boxlim(1), Yrc-boxlim(2);...
        rlim, Yrc-boxlim(2);...
        rlim, Yrc+boxlim(2);...
        Xrc-boxlim(1), Yrc+boxlim(2)];
    lbox = [llim, Ylc-boxlim(2);...
        Xlc+boxlim(1), Ylc-boxlim(2);...
        Xlc+boxlim(1), Ylc+boxlim(2);...
        llim, Ylc+boxlim(2)];
if rbox(1,2) <= 0; rbox(1,2) = 1; rbox(2,2) = 1; end;
if lbox(1,2) <= 0; lbox(1,2) = 1; lbox(2,2) = 1; end;
if lbox(2,1) > 800; lbox(2,1) = 800; lbox(3,1) = 800; end;
if rbox(3,2) > vlim; rbox(3,2) = vlim; rbox(4,2) = vlim; end;
if lbox(3,2) > vlim; lbox(3,2) = vlim; lbox(4,2) = vlim; end;
rbox = int16(rbox);
lbox = int16(lbox);

    Xrc(2:ff) = Xrc(1);
    Xlc(2:ff) = Xlc(1);
    Yrc(2:ff) = Yrc(1);
    Ylc(2:ff) = Ylc(1);

end;



% Right side X
X=[];
     for m = rbox(1,1):rbox(2,1);
         X(m) = sum(bw(rbox(1,2):rbox(3,2),m));
     end;
     rxdist = find(X > cnt);
     Xrc(ff) = rxdist(1);
% Left side X   
X=[];
     for m = lbox(1,1):lbox(2,1);
         X(m) = sum(bw(lbox(1,2):lbox(3,2),m));
     end;
     lxdist = find(X > cnt); 
     Xlc(ff) = lxdist(end);

% Right side Y
Y=[];     
    if rbox(3,2) < rhinge; % Top quadrant
     for m = rbox(1,2):rbox(3,2);
         Y(m) = sum(bw(m,rbox(1,1):rbox(2,1)));
     end;
     rydist = find(Y > cnt);
        Yrc(ff) = rydist(1); 
    end;
    if rbox(1,2) > rhinge; % Bottom quadrant
     for m = rbox(1,2):rbox(3,2);
         Y(m) = sum(bw(m,rbox(1,1):rAnt));
     end;
     rydist = find(Y > cnt);
        Yrc(ff) = rydist(end);
    end;
    if length(Yrc) < ff; % Middle quadrant
     for m = rbox(1,2):rbox(3,2);
         Y(m) = sum(bw(m,rbox(1,1):rAnt));
     end;
     rydist = find(Y > cnt);

        if Yrc(ff-2) < Yrc(ff-1); % Down stroke
            Yrc(ff) = rydist(end);
        else
            Yrc(ff) = rydist(1);
        end;
    end;

% Left side Y
Y=[];
    if lbox(3,2) < lhinge; % Top quadrant
     for m = lbox(1,2):lbox(3,2);
         Y(m) = sum(bw(m,lbox(1,1):lbox(2,1)));
     end;
     lydist = find(Y > cnt);
        Ylc(ff) = lydist(1); 
    end;
    if lbox(1,2) > lhinge; % Bottom quadrant
     for m = lbox(1,2):lbox(3,2);
         Y(m) = sum(bw(m,lAnt:lbox(2,1)));
     end;
     lydist = find(Y > cnt);
        Ylc(ff) = lydist(end);
    end;
    if length(Ylc) < ff; % Middle quadrant
     for m = lbox(1,2):lbox(3,2);
         Y(m) = sum(bw(m,lAnt:lbox(2,1)));
     end;
     lydist = find(Y > cnt);
     if Ylc(ff-2) < Ylc(ff-1); % Down stroke
            Ylc(ff) = lydist(end);
        else
            Ylc(ff) = lydist(1);
        end;
    end;

%     plot(Xrc(k),Yrc(k),'r*', Xlc(k), Ylc(k),'g*');
%     plot([Xrc(k)-boxlim(1) Xrc(k)-boxlim(1)], [Yrc(k)-boxlim(2) Yrc(k)+boxlim(2)], 'w-', 'LineWidth', 2);
%     plot([Xlc(k)+boxlim(1) Xlc(k)+boxlim(1)], [Ylc(k)-boxlim(2) Ylc(k)+boxlim(2)], 'w-', 'LineWidth', 2);
%     
%     plot([Xrc(k)-boxlim(1) rlim], [Yrc(k)-boxlim(2) Yrc(k)-boxlim(2)], 'w-', 'LineWidth', 2);
%     plot([Xrc(k)-boxlim(1) rlim], [Yrc(k)+boxlim(2) Yrc(k)+boxlim(2)], 'w-', 'LineWidth', 2);
%     
%     plot([llim Xlc(k)+boxlim(1)], [Ylc(k)-boxlim(2) Ylc(k)-boxlim(2)], 'w-', 'LineWidth', 2);
%     plot([llim Xlc(k)+boxlim(1)], [Ylc(k)+boxlim(2) Ylc(k)+boxlim(2)], 'w-', 'LineWidth', 2);
%     pause(0.1);

    % Update the boxes
        rbox = [Xrc(ff)-boxlim(1), Yrc(ff)-boxlim(2);...
        rlim, Yrc(ff)-boxlim(2);...
        rlim, Yrc(ff)+boxlim(2);...
        Xrc(ff)-boxlim(1), Yrc(ff)+boxlim(2)];
    
        lbox = [llim, Ylc(ff)-boxlim(2);...
        Xlc(ff)+boxlim(1), Ylc(ff)-boxlim(2);...
        Xlc(ff)+boxlim(1), Ylc(ff)+boxlim(2);...
        llim, Ylc(ff)+boxlim(2)];
if rbox(1,2) <= 0; rbox(1,2) = 1; rbox(2,2) = 1; end;
if lbox(1,2) <= 0; lbox(1,2) = 1; lbox(2,2) = 1; end;
if lbox(2,1) > 800; lbox(2,1) = 800; lbox(3,1) = 800; end;
if rbox(3,2) > vlim; rbox(3,2) = vlim; rbox(4,2) = vlim; end;
if lbox(3,2) > vlim; lbox(3,2) = vlim; lbox(4,2) = vlim; end;
rbox = int16(rbox);
lbox = int16(lbox);
    
    
end;    
    
%     rbox = [Xrc-boxlim(1), Yrc+boxlim(2);...
%         rlim, Yrc+boxlim(2);...
%         rlim, Yrc-boxlim(2);...
%         Xrc-boxlim(1), Yrc-boxlim(2)];
%     lbox = [llim, Ylc+boxlim(2);...
%         Xlc+boxlim(1), Ylc+boxlim(2);...
%         Xlc+boxlim(1), Ylc-boxlim(2);...
%         llim, Ylc-boxlim(2)];
%     
% end;    
%     
%     figure(2); 
%     
%     for zz = 2:length(im);
% 
%     hold off; imshow(im(zz).yg); hold on;
%     plot([rlim, rlim], [1 599], 'w-', 'LineWidth', 2);
%     plot([1 rlim], [rhinge rhinge], 'y-', 'LineWidth', 2);
%     plot([llim, llim], [1 599], 'w-', 'LineWidth', 2);
%     plot([llim, 799], [lhinge, lhinge], 'y-', 'LineWidth', 2);
% 
%     % Get the xlimits of the wings (easy)
%      for m = 1:800;
%          X(m) = sum(im(zz).bw(1:vlim,m));
%      end;
%      xdist = find(X > 10); 
%      Xr = xdist(1);
%      Xl = xdist(end);
%      plot([Xr Xr], [1 599], 'c-');
%      plot([Xl Xl], [1 599], 'c-');    
%     pause(0.5);
% 
%     end;
%     
%     plot([Xr-boxlim(1) Xr-boxlim(1)], [Yr-boxlim(2) Yr+boxlim(2)], 'w-', 'LineWidth', 2);
%     plot([Xl+boxlim(1) Xl+boxlim(1)], [Yl-boxlim(2) Yl+boxlim(2)], 'w-', 'LineWidth', 2);
%     
%     plot([Xr-boxlim(1) rlim], [Yr-boxlim(2) Yr-boxlim(2)], 'w-', 'LineWidth', 2);
%     plot([Xr-boxlim(1) rlim], [Yr+boxlim(2) Yr+boxlim(2)], 'w-', 'LineWidth', 2);
%     
%     plot([llim Xl+boxlim(1)], [Yl-boxlim(2) Yl-boxlim(2)], 'w-', 'LineWidth', 2);
%     plot([llim Xl+boxlim(1)], [Yl+boxlim(2) Yl+boxlim(2)], 'w-', 'LineWidth', 2);

% for q = length(im):-1:1;
% 
%     % Get the xlimits of the wings (easy)
%     for m = 1:800;
%         X(m) = sum(im(q).bw(m,:);
%     end;
% 
%     xdist = find(X > 10); 
%     Xr = xdist(1);
%     Xl = xdist(end);
% 
%     figure(1); subplot(122); hold on; 
%     plot([Xr Xr], [1 599], 'c-');
%     plot([Xl Xl], [1 599], 'c-');

    % Get the ylimits of the wings (not as easy)
    

% for q = 1:17;
% 
%     
%     
%     
% end;
%     