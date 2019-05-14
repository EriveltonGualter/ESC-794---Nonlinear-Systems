clear all
close all
load('matlab.mat')
close all

%% Animation
h=figure(1); box on;

O = [0, 0]; % Origin
axis(gca, 'equal');
if(x(1,end) >l+r) 
    axis([-(l+r) (x(1,end)*2) 0 1.5*(l+r)])
else
    axis([(x(1,end)*2) (l+r) 0 1.5*(l+r)])

end
x = [ -states(1,:)*r; -states(1,:)*r - l*sin(states(2,:)) ];
 
y = [r*ones(size(states(1,:))) ; r*ones(size(states(1,:))) + l*cos(states(2,:)) ];
 

idx = 1:length(time); 
idxq = linspace(min(idx), max(idx), 100);    % Interpolation Vector

timei = interp1(idx, time, idxq, 'linear');       % Downsampled Vector
xi(1,:) = interp1(idx, x(1,:), idxq, 'linear');       % Downsampled Vector
yi(1,:) = interp1(idx, y(1,:), idxq, 'linear');       % Downsampled Vector
xi(2,:) = interp1(idx, x(2,:), idxq, 'linear');       % Downsampled Vector
yi(2,:) = interp1(idx, y(2,:), idxq, 'linear');       % Downsampled Vector

f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,length(xi(1,:))) = 0;
first_frame = true;

 for t=1:length(xi(1,:))
     cla

     % Calculating joint coordinates for animation purpose
    
    pendulum = line([xi(1,t) xi(2,t)], [yi(1,t) yi(2,t)] );
    
    wheel = viscircles([xi(1,t) yi(1,t)],r,'LineStyle','-');
    title(sprintf('Time: %0.2f sec', timei(t)));
    pause(0.0001)
       
    drawnow;    
    
    % gif utilities
    set(gcf,'color','w');   % set figure background to white
    frame = getframe(1);    % get desired frame
    im = frame2im(frame);   % Transform frame to image
    [imind,cm] = rgb2ind(im,256);  % Convert RGB image to indexed image
    outfile = 'simulation2.gif';    % GIF is the BEST. However, you can modify the extensions.

    % On the first loop, create the file. In subsequent loops, append.
    if first_frame
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf); 
        first_frame = false;
    else
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
    end
 end

%imwrite(im,map,'WheeledPendulum.gif','DelayTime',0,'LoopCount',inf) %g443800

% unlive - every quarter 
% casual leave - 2.5 
% sick leave

% figure(3)
% subplot(2,1,1)
% plot(time, states(1:2,:), 'LineWidth', 2);
% hh1(1) = line(time(1), states(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
% hh1(2) = line(time(1), states(2,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
% xlabel('time (sec)'); ylabel('angle (deg)');
% 
% subplot(2,1,2)
% hh2 = plot([O(1), x(1,1);x(1,1), x(2,1)], [O(2), y(1,1);y(1,1), y(2,1)], ...
%       '.-', 'MarkerSize', 20, 'LineWidth', 2);
% axis equal
% axis([-2*(l+r) 2*(l+r) -2*(l+r) 2*(l+r)]);
% ht = title(sprintf('Time: %0.2f sec', time(1)));
% 
% tic;     % start timing
% for id = 1:length(time)
%    % Update XData and YData
%    set(hh1(1), 'XData', time(id)          , 'YData', states(1, id));
%    set(hh1(2), 'XData', time(id)          , 'YData', states(2, id));
%    set(hh2(1), 'XData', [0, x(1, id)]  , 'YData', [0, y(1, id)]);
%    set(hh2(2), 'XData', x(:, id)       , 'YData', y(:, id));
%    set(ht, 'String', sprintf('Time: %0.2f sec', time(id)));
% 
%    drawnow;
% end
% fprintf('Animation (Smart update): %0.2f sec\n', toc);
