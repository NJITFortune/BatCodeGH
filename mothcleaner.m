function out = mothcleaner(in, im)
% Usage: out = mothcleaner(dat, im)
% Employs medfilt1 to remove outliers.
% Allows the user to manually fix the mojority of remaining outliers
% Recalculates velocity and angle
% Version: 3 Nov 2015

%% Transfer data to our output structure
out.left.imnum = in.left.imnum;
out.right.imnum = in.right.imnum;
out.xcenter = in.xcenter;
out.ycenter = in.ycenter;

%% Use a median filter to clean up the data
medlim = 3;

out.left.X = medfilt1(in.left.X,medlim);
out.left.Y = medfilt1(in.left.Y,medlim);
out.right.X = medfilt1(in.right.X,medlim);
out.right.Y = medfilt1(in.right.Y,medlim);

out.left.X2 = medfilt1(in.left.X2,medlim);
out.left.Y2 = medfilt1(in.left.Y2,medlim);
out.right.X2 = medfilt1(in.right.X2,medlim);
out.right.Y2 = medfilt1(in.right.Y2,medlim);

%% Plot the data on the first image
figure(1); clf;
imshow(im(1).orig);
    hold on;
    plot(out.left.X2, out.left.Y2, 'g*-');
    plot(out.right.X2, out.right.Y2, 'r*-');
    plot(in.xcenter/2, in.ycenter/2, 'c*');

%% The main loop
% Initialize our variable for controlling the fixing loop    
loopr = 17;

% Loop until the user gives us the stop signal, which is 99
while loopr ~= 99;

    loopr = input('1:leftDown 2:medialLeft 3:rightDown 4:medialRight 5:leftUp 6:rightUp 99: Done - ');

    % If we aren't quitting, get a click
    if loopr ~= 99; [xout, yout] = ginput(1); plot(xout,yout, 'y*'); end;

    
% Lateral Left %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loopr == 1;
    
    fixpts = find(out.left.X2 > xout & out.left.Y2 > yout);
    goodpts = setdiff(1:length(out.left.X2), fixpts); 
    plot(out.left.X2(fixpts), out.left.Y2(fixpts), 'k*'); % Mark the points to be fixed with black

    for j = 1:length(fixpts);
       % Get the closest valid points on either side of the bad point 
            pg = goodpts(find(goodpts < fixpts(j), 1, 'last' ));
            ng = goodpts(find(goodpts > fixpts(j), 1 ));
%           pg = fixpts(j) - find(out.left.X2(max([1, fixpts(j)-10]):max([1, fixpts(j)-1])) < xout, 1, 'first'); % Positive direction 
%           ng = fixpts(j) + find(out.left.X2(min([fixpts(j)+1, length(out.left.X2)]):min([fixpts(j)+10, length(out.left.X2)])) < xout, 1, 'first'); % Negative direciton
       % Put the fixed data into the output variable             
            out.left.X2(fixpts(j)) = mean([out.left.X2(pg) out.left.X2(ng)]);
            out.left.Y2(fixpts(j)) = mean([out.left.Y2(pg) out.left.Y2(ng)]);
            out.left.X(fixpts(j)) = mean([out.left.X(pg) out.left.X(ng)]);
            out.left.Y(fixpts(j)) = mean([out.left.Y(pg) out.left.Y(ng)]);
       % Plot the fixed data points
            plot(out.left.X2(fixpts(j)), out.left.Y2(fixpts(j)), 'm*'); 
            plot([out.left.X2(pg) out.left.X2(ng)] , [out.left.Y2(pg) out.left.Y2(ng)], 'w*-');        
    end;
end;

% Medial Left %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loopr == 2;
    
    fixpts = find(out.left.X2 < xout & out.left.Y2 > yout); 
    goodpts = setdiff(1:length(out.left.X2), fixpts); 
    plot(out.left.X2(fixpts), out.left.Y2(fixpts), 'k*');
    
    for j = 1:length(fixpts);
       % Get the closest valid points on either side of the bad point
            pg = goodpts(find(goodpts < fixpts(j), 1, 'last' ));
            ng = goodpts(find(goodpts > fixpts(j), 1 ));
            %pg = fixpts(j) - find(out.left.X2(max([1, fixpts(j)-10]):max([1, fixpts(j)-1])) < yout, 1, 'first'); 
            %ng = fixpts(j) + find(out.left.X2(min([fixpts(j)+1, length(out.left.X2)]):min([fixpts(j)+10, length(out.left.X2)])) < yout, 1, 'first'); 

            % Put the fixed data into the output variable             
            out.left.X2(fixpts(j)) = mean([out.left.X2(pg) out.left.X2(ng)]);
            out.left.Y2(fixpts(j)) = mean([out.left.Y2(pg) out.left.Y2(ng)]);
            out.left.X(fixpts(j)) = mean([out.left.X(pg) out.left.X(ng)]);
            out.left.Y(fixpts(j)) = mean([out.left.Y(pg) out.left.Y(ng)]);
        % Plot the fixed data points
            plot(out.left.X2(fixpts(j)), out.left.Y2(fixpts(j)), 'm*');  
            plot([out.left.X2(pg) out.left.X2(ng)] , [out.left.Y2(pg) out.left.Y2(ng)], 'w-*');
    end;
end;

% Lateral Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loopr == 3;

        fixpts = find(out.right.X2 < xout & out.right.Y2 > yout);
        goodpts = setdiff(1:length(out.right.X2), fixpts); 
        plot(out.right.X2(fixpts), out.right.Y2(fixpts), 'k*');
   
    for j = 1:length(fixpts); 
        
       % Get the closest valid points on either side of the bad point 
            pg = goodpts(find(goodpts < fixpts(j), 1, 'last' ));
            ng = goodpts(find(goodpts > fixpts(j), 1 ));
%           pg = fixpts(j) - find(out.right.X2(max([1, fixpts(j)-10]):max([1, fixpts(j)-1])) < xout, 1, 'first') ;
%           ng = fixpts(j) + find(out.right.X2(min([fixpts(j)+1, length(out.right.X2)]):min([fixpts(j)+10, length(out.right.X2)])) < xout, 1, 'first') ;
       % Put the fixed data into the output variable             
            out.right.X2(fixpts(j)) = mean([out.right.X2(pg) out.right.X2(ng)]);
            out.right.Y2(fixpts(j)) = mean([out.right.Y2(pg) out.right.Y2(ng)]);
            out.right.X(fixpts(j)) = mean([out.right.X(pg) out.right.X(ng)]);
            out.right.Y(fixpts(j)) = mean([out.right.Y(pg) out.right.Y(ng)]);
       % Plot the fixed data points
            plot(out.right.X2(fixpts(j)), out.right.Y2(fixpts(j)), 'm*');   
            plot([out.right.X2(pg) out.right.X2(ng)] , [out.right.Y2(pg) out.right.Y2(ng)], 'w*-');
    end;
end;

% Medial Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loopr == 4;
    
    fixpts = find(out.right.X2 > xout & out.right.Y2 > yout);
    goodpts = setdiff(1:length(out.right.X2), fixpts); 
    plot(out.right.X2(fixpts), out.right.Y2(fixpts), 'k*');
    
    for j = 1:length(fixpts);
        
       % Get the closest valid points on either side of the bad point 
             pg = goodpts(find(goodpts < fixpts(j), 1, 'last' ));
             ng = goodpts(find(goodpts > fixpts(j), 1 ));
%            pg = fixpts(j) - find(out.right.Y2(max([1, fixpts(j)-10]):max([1, fixpts(j)-1])) < yout, 1, 'first'); 
%            ng = fixpts(j) + find(out.right.Y2(min([fixpts(j)+1, length(out.right.X2)]):min([fixpts(j)+10, length(out.right.X2)])) < yout, 1, 'first');
       % Put the fixed data into the output variable             
            out.right.X2(fixpts(j)) = mean([out.right.X2(pg) out.right.X2(ng)]);
            out.right.Y2(fixpts(j)) = mean([out.right.Y2(pg) out.right.Y2(ng)]);
            out.right.X(fixpts(j)) = mean([out.right.X(pg) out.right.X(ng)]);
            out.right.Y(fixpts(j)) = mean([out.right.Y(pg) out.right.Y(ng)]);
       % Plot the fixed data points
            plot(out.right.X2(fixpts(j)), out.right.Y2(fixpts(j)), 'm*');     
            plot([out.right.X2(pg) out.right.X2(ng)] , [out.right.Y2(pg) out.right.Y2(ng)], 'w*-');
    end;
end;

% Upper Left %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loopr == 5;
    
    fixpts = find(out.left.X2 > xout & out.left.Y2 < yout);
    goodpts = setdiff(1:length(out.left.X2), fixpts); 
    plot(out.left.X2(fixpts), out.left.Y2(fixpts), 'k*');
    
    for j = 1:length(fixpts);
        
       % Get the closest valid points on either side of the bad point 
             pg = goodpts(find(goodpts < fixpts(j), 1, 'last' ));
             ng = goodpts(find(goodpts > fixpts(j), 1 ));
%            pg = fixpts(j) - find(out.right.Y2(max([1, fixpts(j)-10]):max([1, fixpts(j)-1])) < yout, 1, 'first'); 
%            ng = fixpts(j) + find(out.right.Y2(min([fixpts(j)+1, length(out.right.X2)]):min([fixpts(j)+10, length(out.right.X2)])) < yout, 1, 'first');
       % Put the fixed data into the output variable             
            out.left.X2(fixpts(j)) = mean([out.left.X2(pg) out.left.X2(ng)]);
            out.left.Y2(fixpts(j)) = mean([out.left.Y2(pg) out.left.Y2(ng)]);
            out.left.X(fixpts(j)) = mean([out.left.X(pg) out.left.X(ng)]);
            out.left.Y(fixpts(j)) = mean([out.left.Y(pg) out.left.Y(ng)]);
       % Plot the fixed data points
            plot(out.left.X2(fixpts(j)), out.left.Y2(fixpts(j)), 'm*');     
            plot([out.left.X2(pg) out.left.X2(ng)] , [out.left.Y2(pg) out.left.Y2(ng)], 'w*-');
    end;
end;

% Upper Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loopr == 6;
    
    fixpts = find(out.right.X2 < xout & out.right.Y2 < yout);
    goodpts = setdiff(1:length(out.right.X2), fixpts); 
    plot(out.right.X2(fixpts), out.right.Y2(fixpts), 'k*');
    
    for j = 1:length(fixpts);
        
       % Get the closest valid points on either side of the bad point 
             pg = goodpts(find(goodpts < fixpts(j), 1, 'last' ));
             ng = goodpts(find(goodpts > fixpts(j), 1 ));
%            pg = fixpts(j) - find(out.right.Y2(max([1, fixpts(j)-10]):max([1, fixpts(j)-1])) < yout, 1, 'first'); 
%            ng = fixpts(j) + find(out.right.Y2(min([fixpts(j)+1, length(out.right.X2)]):min([fixpts(j)+10, length(out.right.X2)])) < yout, 1, 'first');
       % Put the fixed data into the output variable             
            out.right.X2(fixpts(j)) = mean([out.right.X2(pg) out.right.X2(ng)]);
            out.right.Y2(fixpts(j)) = mean([out.right.Y2(pg) out.right.Y2(ng)]);
            out.right.X(fixpts(j)) = mean([out.right.X(pg) out.right.X(ng)]);
            out.right.Y(fixpts(j)) = mean([out.right.Y(pg) out.right.Y(ng)]);
       % Plot the fixed data points
            plot(out.right.X2(fixpts(j)), out.right.Y2(fixpts(j)), 'm*');     
            plot([out.right.X2(pg) out.right.X2(ng)] , [out.right.Y2(pg) out.right.Y2(ng)], 'w*-');
    end;
end;



end;

%% Recalculate velocity

% Determine whether the wing is going up or down on both sides
 
leftup = find(out.left.Y(1:end-1) - out.left.Y(2:end) >= 0); %+ 1; what matt said
rightup = find(out.right.Y(1:end-1) - out.right.Y(2:end) >= 0); %+ 1;
 
% Calculate wing velocities
for p = 1:length(out.left.Y)-1;    
     out.right.vel(p) = -1 * pdist([out.right.X(p),out.right.Y(p); out.right.X(p+1), out.right.Y(p+1)]);
     out.left.vel(p) = -1 * pdist([out.left.X(p),out.left.Y(p); out.left.X(p+1), out.left.Y(p+1)]);
end;
% 
    out.right.vel(rightup) = -out.right.vel(rightup);
    out.left.vel(leftup) = -out.left.vel(leftup);

%% Recalculate angle

for k = 1:length(out.right.Y);
     out.right.angle(k) = atan2d(out.right.Y(k)-out.ycenter, out.right.X(k)-out.xcenter);
     out.left.angle(k) = atan2d(out.left.Y(k)-out.ycenter, out.left.X(k)-out.xcenter);
end;

kk = find(out.right.angle > 160); out.right.angle(kk) = -180 + (out.right.angle(kk) - 180);
