% wyswietl wykres dla pierwszego skanu
function [theta_end, rho_end, marker_points] =  get_lines_from_scan(laser_range_data, print_figures, rotate)
theta_end = [0,0];
rho_end = [0,0];
% theta from -135 to 135 deg
theta=[-2.35619449615:degtorad(0.5):2.35619449615];
[points(1,:),points(2,:)] = pol2cart(theta, laser_range_data(1,:));
%max(x)
%%// Remove frames
%plot_gca = gca;
if rotate == -1
    transf_points = [0 0];
    my_rot_z = rotz(-45);
    [size_1, size2] = size(points(1,:));
    for i = 1 : size2
        transf_points(1,i) = my_rot_z(1,1) * points(1,i);
        transf_points(1,i) = transf_points(1,i) + my_rot_z(1,2) * points(2,i);

        transf_points(2,i) = my_rot_z(2,1) * points(1,i);
        transf_points(2,i) = transf_points(2,i) + my_rot_z(2,2) * points(2,i);
    end
    points = transf_points;
end

if rotate == 1
    transf_points = [0 0];
    my_rot_z = rotz(45);
    [size_1, size2] = size(points(1,:));
    for i = 1 : size2
        transf_points(1,i) = my_rot_z(1,1) * points(1,i);
        transf_points(1,i) = transf_points(1,i) + my_rot_z(1,2) * points(2,i);

        transf_points(2,i) = my_rot_z(2,1) * points(1,i);
        transf_points(2,i) = transf_points(2,i) + my_rot_z(2,2) * points(2,i);
    end
    points = transf_points;
end
%%// Get the figure as a uint8 variable
%F = getframe(gcf);
%im = frame2im(F);
%imshow(im)
my_figure = figure;
plot(points(1,:),points(2,:),'.')
%axis([-5 5 -5 5])

ransac_thDist = 0.02;

input_points = ginput;
points_size = size(points(1,:));
my_points = [0;0];
for i= 1:points_size(2)
    if points(1,i) > input_points(1,1)
        if points(2,i) < input_points(1,2)
            my_points = horzcat(points(:,i),my_points);
        end
    end
end
my_points_2 = [0;0];
%plot(my_points(1,:),my_points(2,:),'.');
%axis([-5 5 -5 5])

points_size = size(my_points(1,:));
for i= 1:points_size(2)
    if my_points(1,i) < input_points(2,1)
        if my_points(2,i) > input_points(2,2)
            my_points_2 = horzcat(my_points(:,i),my_points_2);
        end
    end
end

 [theta_end(1),rho_end(1)]=ransac(my_points_2,10,ransac_thDist,0.001);

 if print_figures == 1
	plot(my_points_2(1,:),my_points_2(2,:),'.');
    p1 = (rho_end(1) - input_points(1,1)* sin(theta_end(1)) )/ cos(theta_end(1)); 
    p2 = (rho_end(1) - input_points(2,1)* sin(theta_end(1)) )/ cos(theta_end(1)); 
    hold on;
    line([input_points(1,1) input_points(2,1)],[p1 p2], 'Color','green');
 end

my_points_3 = [0;0];
[my_points_2_size_1,my_points_2_size_2] = size(my_points_2(1,:));
for i = 1:my_points_2_size_2
    dist = abs(sin(theta_end(1))*my_points_2(1,i)+cos(theta_end(1))*my_points_2(2,i)...
        - rho_end(1) )/sqrt(sin(theta_end(1))^2 + cos(theta_end(1))^2);
    if dist > ransac_thDist
        my_points_3 = horzcat(my_points_2(:,i),my_points_3);
    end
    
end
	 if (print_figures == 1)
        [theta_end(2),rho_end(2)]=ransac(my_points_3,10,ransac_thDist,0.001);
        p1 = (rho_end(2) - input_points(1,1)* sin(theta_end(2)) )/ cos(theta_end(2));
        p2 = (rho_end(2) - input_points(2,1)* sin(theta_end(2)) )/ cos(theta_end(2));
        hold on;
        line([input_points(1,1) input_points(2,1)],[p1 p2],'Color', 'red');

        %  axis([0 1.6 0 1.4])
     end
    
    marker_points = my_points_2;
    
    
end


