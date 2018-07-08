function [ output_args ] = mothfunplot( im, dat, savmov )
%mothfunplot Summary of this function goes here
%   Detailed explanation goes here

trailen = 20;


figure(1); clf;
subplot(3,2,1); 
subplot(3,2,2); axis([0, 640, -480, 0]);
subplot(3,1,2); 

Ymin = min([min(dat.left.Y) min(dat.right.Y)]);
Ymax = max([max(dat.left.Y) max(dat.right.Y)]);
subplot(3,1,2); axis([0, 1+trailen*4, Ymin, Ymax]);

lvel = medfilt1(dat.left.vel, 4);
rvel = medfilt1(dat.right.vel, 4);
Vmin = min([min(lvel) min(rvel)]);
Vmax = max([max(lvel) max(rvel)]);
subplot(3,1,3); axis([0, 1+trailen*4, Vmin, Vmax]);

%% Main loop for each frame

for k = trailen+1:length(im)-((3*trailen)+2);

subplot(3,2,1); % The image with wing trails

    imshow(im(k).orig); 
    hold on;
    plot(dat.left.X2((k-trailen):k), dat.left.Y2((k-trailen):k), 'g*-');
    plot(dat.right.X2((k-trailen):k), dat.right.Y2((k-trailen):k), 'r*-');
    hold off;
    
subplot(3,2,2); % Stick figure
    plot([dat.xcenter dat.left.X(k)], [-dat.ycenter -dat.left.Y(k)], 'g*-');
    axis([0, 640, -480, 0]);
    hold on;
    plot([dat.xcenter dat.right.X(k)], [-dat.ycenter -dat.right.Y(k)], 'r*-');
    plot(dat.xcenter, -dat.ycenter, 'k*');    
    hold off;
    
 subplot(3,1,2); % Y Position
 
    plot(dat.left.Y(k-trailen:k+(trailen*3)),'g'); 
    axis([0, 1+trailen*4, Ymin, Ymax]);
    hold on;
    plot(dat.right.Y(k-trailen:k+(trailen*3)),'r'); 
    plot([trailen trailen], [Ymin Ymax], 'k');
    text(trailen+1, Ymax/2, 'Y-Position');
    hold off;
    
 subplot(3,1,3); % Velocity
 
    plot(lvel(k-trailen:k+(trailen*3)),'g'); 
    axis([0, 1+trailen*4, Vmin, Vmax]);
    hold on;
    plot(rvel(k-trailen:k+(trailen*3)),'r'); 
    plot([trailen trailen], [Vmin Vmax], 'k');
    text(trailen+1, Vmax/2, 'Wingtip Velocity');
    hold off;
    
    pause(0.01);
    
end

