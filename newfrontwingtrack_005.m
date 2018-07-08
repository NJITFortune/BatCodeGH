function [im, dat] = newfrontwingtrack_005(lims)
% [im dat] = newfrontwingtrack(lims) (now just "dat")
% cd into the directory with the jpeg files
% Run this code with the first frame and last frame. Leave empty if you
% want all frames.
% THE KEY IS TO start with wings high or low in the first frame.
% Change some defaults if necessary:
%    
% 

%% Defaults and setup

flist = dir('fem3/*.jpg');

     xbuff = 40; % This is the X-window width for looking for Y position of wingtips.
     ybuff = 20; % This is the X-window width for looking for Y position of wingtips.

    % If the user does not specify the limits, do the whole list of jpegs.
    if nargin == 0; % zero 'arguments'
        lims(1) = 2; % Always skip the first frame...
        lims(2) = length(flist); % Go to the end.
    end;
        
    % Always avoid the first frame of the series (catch user errors)
    if lims(1) == 1; lims(1) = 2; lims(2) = lims(2)-1; end;

%% First click

    cd fem3; firstimage = rgb2gray(imread(flist(1).name)); cd ..;

    figure(1); clf; imshow(firstimage);
    fprintf('Click between the eyes.\n');
    [xcenter, ycenter] = ginput(1);
    hold on; plot(xcenter, ycenter, 'g*');
    clear firstimage;
    
%% Read the image files 
    
fprintf('Reading images.\n');

    % How many frames will we process?
    numframes = (lims(2) - lims(1)) + 1;
    
    % Read each frame
    cd fem3;
    for k = lims(2):-1:lims(1); 
        pp = k-lims(1)+1; %This is used to renumber the frames to start at 1
        g(:, :, pp) = rgb2gray(imread(flist(k).name)); 
    end; %It converts the images to gray scale
    cd ..;
    
    % get a diff for each pair of frames
    qNumFrames = lims(2)-lims(1);
    for k = qNumFrames:-1:1; 
        d(:, :, k) = imabsdiff(g(:, :, k), g(:, :, k+1)); 
    end; %It finds the differences between each frame, 
    %it saves the movement in d    

    fprintf('Thesholding differences.\n');

    % Threshold the difference images
    thresh = 0.03;
    bw = (d >= thresh * 255); % THRESHOLD IS CRITICAL
    % Get rid of little speckles and stuff
    bw2 = bwareaopen(bw, 20, 8); % Gets rid of all white areas that are less than 20 pixels large

    fprintf('Tracking the wings.\n');
    
%% Find X and Y of outer edges of wings

for k = qNumFrames:-1:1;
    
    % NOTE it is really sum(bw2(Y,X,framenumber));
    for j=30:610; hsum(k,j) = sum(bw2(ybuff:end-ybuff,j,k)); end; %This is summing the vertical lines to find the 1s
    hindx{k} = find(hsum(k,:));
    
    % Get the X boundaries of the wings in index numbers
     minX(k) = min(hindx{k}); % Finding the min of hsum, created in hindx, Left Side, Right Wing
     maxX(k) = max(hindx{k}); % Right Side, Left Wing

     % Make a vertical strip buff wide inside of both maxX and minX and sum
     % each horizontal row.
     for j=ybuff:480-ybuff; tmpLeftVsum(j) = sum(bw2(j,minX(k):minX(k)+xbuff,k)); end; % Left Side, Right Wing
     for j=ybuff:480-ybuff; tmpRightVsum(j) = sum(bw2(j,maxX(k)-xbuff:maxX(k),k)); end; % Right Side, Left Wing

     % Get the tops and bottoms for each side
     LeftTopIdx(k) = find(tmpLeftVsum, 1, 'last');
     LeftBottomIdx(k) = find(tmpLeftVsum, 1, 'first');
     RightTopIdx(k) = find(tmpRightVsum, 1, 'last' );
     RightBottomIdx(k) = find(tmpRightVsum, 1, 'first');

     
     % Loop to find most distant point in our wingbox
     
%% Moth's Left wing (our Right side) [maxX(k)-xbuff maxX(k)], [RightBottomIdx(k) RightTopIdx(k)]
     
     % Detect the vertical location of the box
        tmpDistRight = [];
        
        if RightTopIdx(k) < ycenter; % Above the middle point (most frames)
            clear tmpDistRight; 
            for rr = RightBottomIdx(k):RightBottomIdx(k)+xbuff; % Cycle through each row of our box
                if sum(bw2(rr,maxX(k)-xbuff:maxX(k),k)) > 1; % Some rows will be all black, so we avoid those
                    furthestX(rr) = find(bw2(rr, maxX(k)-xbuff:maxX(k), k), 1, 'last') + maxX(k)-xbuff; % Take the last white pixel
                    tmpDistRight(rr) = pdist([xcenter, ycenter; furthestX(rr), rr]); % Calculate the distance to the center
                    matchIdxR(rr) = rr;
                end;
            end;
            if exist('tmpDistRight','var') == 1;
                [~, tmpidx] = max(tmpDistRight);
                LeftX(k) = furthestX(tmpidx);
                LeftY(k) = matchIdxR(tmpidx);
            end;
            
        elseif RightBottomIdx(k) > ycenter; % Below the middle point (some frames)
            clear tmpDistRight;
            for rr = RightTopIdx(k)-xbuff:RightTopIdx(k); % Cycle through each row of our box
                if sum(bw2(rr,maxX(k)-xbuff:maxX(k),k)) > 1; % Some rows will be all black, so we avoid those
                    furthestX(rr) = find(bw2(rr, maxX(k)-xbuff:maxX(k), k), 1, 'last') + maxX(k)-xbuff; % Take the last white pixel
                    tmpDistRight(rr) = pdist([xcenter, ycenter; furthestX(rr), rr]); % Calculate the distance to the center
                    matchIdxR(rr) = rr;
                end;
            end;
            if exist('tmpDistRight','var') == 1;
                [~, tmpidx] = max(tmpDistRight);
                LeftX(k) = furthestX(tmpidx);
                LeftY(k) = matchIdxR(tmpidx);            
            end;
                
        elseif RightBottomIdx(k) < ycenter && RightTopIdx(k) > ycenter;
            clear tmpDistRight;
            for rr = RightBottomIdx(k):RightTopIdx(k); % Cycle through each row of our box
                if sum(bw2(rr,maxX(k)-xbuff:maxX(k),k)) > 1; % Some rows will be all black, so we avoid those
                    furthestX(rr) = find(bw2(rr, maxX(k)-xbuff:maxX(k), k), 1, 'last') + maxX(k)-xbuff; % Take the last white pixel
                    tmpDistRight(rr) = pdist([xcenter, ycenter; furthestX(rr), rr]); % Calculate the distance to the center
                    matchIdxR(rr) = rr;
                end;
            end;
            if exist('tmpDistRight','var') == 1;
                [~, tmpidx] = max(tmpDistRight);
                LeftX(k) = furthestX(tmpidx);
                LeftY(k) = matchIdxR(tmpidx);                        
            end;
            
        end;
     
%% Moth's Right wing (our Left side) [minX(k) maxX(k)+xbuff], [LefttBottomIdx(k) LefttTopIdx(k)]
     
     % Detect the vertical location of the box
        tmpDistLeft = [];
        
        if LeftTopIdx(k) < ycenter; % Above the middle point (most frames)
            clear tmpDistLeft; 
            for rr = LeftBottomIdx(k):LeftBottomIdx(k)+xbuff; % Cycle through each row of our box
                if sum(bw2(rr,minX(k):minX(k)+xbuff,k)) > 1; % Some rows will be all black, so we avoid those
                    furthestX(rr) = find(bw2(rr, minX(k):minX(k)+xbuff, k), 1, 'first') + minX(k); % Take the first white pixel
                    tmpDistLeft(rr) = pdist([xcenter, ycenter; furthestX(rr), rr]); % Calculate the distance to the center
                    matchIdxL(rr) = rr;
                end;
            end;
            
            if exist('tmpDistLeft','var') == 1;
                [~, tmpidx] = max(tmpDistLeft);
                RightX(k) = furthestX(tmpidx);
                RightY(k) = matchIdxL(tmpidx);
            end;
                
        elseif LeftBottomIdx(k) > ycenter; % Below the middle point (some frames)
            clear tmpDistLeft;
            for rr = LeftTopIdx(k)-xbuff:LeftTopIdx(k); % Cycle through each row of our box
                if sum(bw2(rr,minX(k):minX(k)+xbuff,k)) > 1; % Some rows will be all black, so we avoid those
                    furthestX(rr) = find(bw2(rr, minX(k):minX(k)+xbuff, k), 1, 'first') + minX(k); % Take the first white pixel
                    tmpDistLeft(rr) = pdist([xcenter, ycenter; furthestX(rr), rr]); % Calculate the distance to the center
                    matchIdxL(rr) = rr;
                end;
            end;
            if exist('tmpDistLeft','var') == 1;
                [~, tmpidx] = max(tmpDistLeft);
                RightX(k) = furthestX(tmpidx);
                RightY(k) = matchIdxL(tmpidx);            
            end;
                
        elseif LeftBottomIdx(k) < ycenter && LeftTopIdx(k) > ycenter;
            clear tmpDistLeft;
            for rr = LeftBottomIdx(k):LeftTopIdx(k); % Cycle through each row of our box
                if sum(bw2(rr,minX(k):minX(k)+xbuff,k)) > 1; % Some rows will be all black, so we avoid those
                    furthestX(rr) = find(bw2(rr, minX(k):minX(k)+xbuff, k), 1, 'first') + minX(k); % Take the first white pixel
                    tmpDistLeft(rr) = pdist([xcenter, ycenter; furthestX(rr), rr]); % Calculate the distance to the center
                    matchIdxL(rr) = rr;
                end;
            end;
            
            if exist('tmpDistLeft','var') == 1;
                [~, tmpidx] = max(tmpDistLeft);
                RightX(k) = furthestX(tmpidx);
                RightY(k) = matchIdxL(tmpidx);                     
            end;
            
        end;
        
% OLD: YY These finds pull out where the wing is within our verticle slices...
%      yMinidx(k) = mean(find(tmpLeftVsum)); % Left Side, Right Wing
%      yMaxidx(k) = mean(find(tmpRightVsum)); % Right Side, Left Wing     
     
     % Calculate the angle from the center of the moth to the putative
     % wingtips


plotme = 919;     
        if plotme == 99;
             figure(2); 
             clf; 
             subplot(131);     imshow(bw2(:,:,k));
             hold on; 
             text(100, 100, num2str(k), 'Color', 'w');
             plot([maxX(k)-xbuff maxX(k)], [RightTopIdx(k) RightTopIdx(k)],'c');
             plot([maxX(k)-xbuff maxX(k)], [RightBottomIdx(k) RightBottomIdx(k)],'b');
             plot([maxX(k)-xbuff maxX(k)-xbuff], [RightTopIdx(k) RightBottomIdx(k)],'g');
             plot([maxX(k) maxX(k)], [RightTopIdx(k) RightBottomIdx(k)],'g');

             plot([minX(k) minX(k)+xbuff], [LeftTopIdx(k) LeftTopIdx(k)],'c');
             plot([minX(k) minX(k)+xbuff], [LeftBottomIdx(k) LeftBottomIdx(k)],'b');
             plot([minX(k)+xbuff minX(k)+xbuff], [LeftTopIdx(k) LeftBottomIdx(k)],'r');
             plot([minX(k) minX(k)], [LeftTopIdx(k) LeftBottomIdx(k)],'r');

             plot([50 610], [ycenter ycenter], 'm'); plot(xcenter,ycenter,'g*');

             plot(LeftX(k), LeftY(k), 'g*');
             plot(RightX(k), RightY(k), 'r*');

             subplot(132);
             imshow(bw2(LeftBottomIdx(k):LeftTopIdx(k), minX(k):minX(k)+xbuff, k));
             subplot(133);
             imshow(bw2(RightBottomIdx(k):RightTopIdx(k), maxX(k)-xbuff:maxX(k), k));

             pause(0.01);
        end;     

end;


%% Output to structure

for k=qNumFrames:-1:1;         
     
        if RightY(k) == 0; RightY(k) = mean([RightY(k-1) RightY(k+1)]); end;
        if LeftY(k) == 0; LeftY(k) = mean([LeftY(k-1) LeftY(k+1)]); end;
        
         dat.right.Y(k) = RightY(k);
         dat.left.Y(k)  = LeftY(k);
         
         dat.right.X(k) = minX(k);
         dat.left.X(k)  = maxX(k);
         
     Rangle(k) = atan2d(RightX(k)-ycenter, minX(k)-xcenter);
     Langle(k) = atan2d(LeftX(k)-ycenter, maxX(k)-xcenter);
         
         dat.right.angle(k) = Rangle(k);
         dat.left.angle(k) = Langle(k);
                  
end

    dat.xcenter = xcenter; dat.ycenter = ycenter;

%%
    fprintf('Plotting.\n');

     figure(3); clf;

        imshow(bw2(:,:,k)); hold on;
        plot(dat.right.X, dat.right.Y, '*-r');
        plot(dat.left.X, dat.left.Y, '*-g');

        %plot(LeftX, LeftY, '*w');

        dat.LeftX = LeftX; dat.LeftY = LeftY;
        
%%
    fprintf('Final calculations and generation of output movie.\n');

dat.left.imnum = lims(1):lims(2);
dat.right.imnum = lims(1):lims(2);
 
% Determine whether the wing is going up or down on both sides
 
leftup = find(dat.left.Y(1:end-1) - dat.left.Y(2:end) >= 0); %+ 1; what matt said
rightup = find(dat.right.Y(1:end-1) - dat.right.Y(2:end) >= 0); %+ 1;
 
% Calculate wing velocities
for p = 1:length(dat.left.Y)-1;    
     dat.right.vel(p) = -1 * pdist([dat.right.X(p),dat.right.Y(p); dat.right.X(p+1), dat.right.Y(p+1)]);
     dat.left.vel(p) = -1 * pdist([dat.left.X(p),dat.left.Y(p); dat.left.X(p+1), dat.left.Y(p+1)]);
end;
% 
dat.right.vel(rightup) = -dat.right.vel(rightup);
dat.left.vel(leftup) = -dat.left.vel(leftup);
% 
% hold on;

    % Pre-allocate the main variables
   xlen = floor(length(g(:,1,1))/2);
   ylen = floor(length(g(1,:,1))/2);
   
   im(numframes).orig = zeros(xlen,ylen); % This is the biggest memory hog - the stack of images.
   im(qNumFrames).bw = zeros(xlen,ylen); % This is the biggest memory hog - the stack of images.

   for k=1:numframes; 
       im(k).orig = g(1:2:end,1:2:end,k); 
   end;
   for k=1:qNumFrames;
       im(k).bw = bw2(1:2:end,1:2:end,k);
   end;

    dat.right.X2 = dat.right.X / 2;
    dat.left.X2 = dat.left.X / 2;
    dat.right.Y2 = dat.right.Y / 2;
    dat.left.Y2 = dat.left.Y / 2;
end
   

%% Award-winning code

% lastGoodRightY = [];
% lastGoodLeftY = [];
% lastGoodRightX = [];
% lastGoodLeftX = []; 

     % making sure the Y value isnt too low (meaning at the top of the
     % graph) and if it is something weird, uses that last good (x,y) value
%      if ceil(mean(tmpMinidx)) <= 40  %ceil rounds up to the nearest whole number
%         tmpMinidx = lastGoodRightY;
%         minX = lastGoodRightX;
%      end
% 
%      if ceil(mean(tmpMaxidx)) <= 40
%         tmpMaxidx = lastGoodLeftY;
%         maxX = lastGoodLeftX;
%      end
%      
%      % And we take the mean of that to pick a single point (in the middle)
%      lastGoodRightY = mean(tmpMinidx);
%         dat.right.Y(k) = mean(tmpMinidx);
%         if isnan(dat.right.Y(k)) == 1; 
%             figure;
%             imshow(bw2(:,:,k)); pause(0.5);
%             hold on; plot([maxX maxX], [1 440], 'g'); plot([minX minX], [1 440], 'r'); hold off;
%         end;
%      lastGoodLeftY  = mean(tmpMaxidx);
%       
%      lastGoodRightX = minX;
%      dat.right.X(k) = minX;
%      lastGoodLeftX  = maxX;
%      dat.left.X(k)  = maxX;     
     

     %sanity-checks to make sure the min and max values are on the correct
     %side of the moth
%      if minX >= 200
%         dat.right.Y(k) = lastGoodRightY;
%         dat.left.Y(k) = lastGoodLeftY;
%         dat.right.X(k) = lastGoodRightX;
%         dat.left.X(k) = lastGoodLeftX; 
%         continue;
%      end
% 
%      if maxX <= 300
%         dat.right.Y(k) = lastGoodRightY;
%         dat.left.Y(k) = lastGoodLeftY;
%         dat.right.X(k) = lastGoodRightX;
%         dat.left.X(k) = lastGoodLeftX; 
%         continue;
%      end
     
     %making sure there are actually values found in min and max
%      if isnan(minX) == 1;         
%          hold on;
%          for i=1:5;
%          plot(dat.right.X(k-i), dat.right.Y(k-i), 'r*-');
%          plot(dat.left.X(k-i), dat.left.Y(k-i), 'g*-');
%          end;
%         hold off;
%         tmp = input('PROBLEM!');
%         dat.right.Y(k) = lastGoodRightY;
%         dat.left.Y(k) = lastGoodLeftY;
%         dat.right.X(k) = lastGoodRightX;
%         dat.left.X(k) = lastGoodLeftX; 
%         continue;
%      end

   
