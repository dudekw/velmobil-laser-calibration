close all;
file_front = 'laser_front.csv';
file_rear = 'laser_rear.csv';
laser_front_data = csvread(file_front,1,11,[1 11 302 551]);
laser_rear_data = csvread(file_rear,1,11,[1 11 276 551]);
theta_lines = [0 0;
               0 0];
rho_lines = [0 0;
             0 0];
init_dist = [0 0 0];
init_rot  = rotz(180);

%%
%  IMPORTANT  !!!!!
%
% use average value of computed rotations
average = 1;

% IF don't use the average, specify REAR LASER line ID to align 
% FROT LASER line

% align to line 1 or line 2
%align = 1;
%align = 2;

% File format
% file_format = 'yaml';
file_format = 'xacro';

%%
% Calculate lines for front laser

[theta_lines(1,:),rho_lines(1,:), marker_points_1] = get_lines_from_scan(laser_front_data,1,0);
%%
% get coss point for front laser

R1 = (rho_lines(1,1));
R2 = (rho_lines(1,2));
sin_1 = sin(theta_lines(1,1));
cos_1 = cos(theta_lines(1,1));
sin_2 = sin(theta_lines(1,2));
cot_1 = cot(theta_lines(1,1));
cos_2 = cos(theta_lines(1,2));

y_cross_front = (R2 - (R1 * sin_2/sin_1 ))/(-cot_1*sin_2+cos_2);
x_cross_front = -(y_cross_front*cos_2 - R2)/sin_2;
%%
% get two additional points (one for each line) for x = x_cross_front - 0.5
sec_point_1_1 = [0 0];
sec_point_1_1(2) = -2;
sec_point_1_1(1) = (R1 - sec_point_1_1(2) * cos_1)/sin_1;

sec_point_1_2 = sec_point_1_1;
sec_point_1_2(2) = -2;
sec_point_1_2(1) = (R2 - sec_point_1_2(2) * cos_2)/sin_2;
%%
% print cross point (o) and additional points (x)
hold on;
plot(x_cross_front,y_cross_front, 'o');
hold on;
plot(sec_point_1_1(1),sec_point_1_1(2), 'x');
hold on;
plot(sec_point_1_2(1),sec_point_1_2(2), 'x');

%%
% Calculate lines for rear laser

[theta_lines(2,:),rho_lines(2,:),marker_points_2] = get_lines_from_scan(laser_rear_data,1,0);
%%
% get coss point for rear laser
R1 = (rho_lines(2,1));
R2 = (rho_lines(2,2));
sin_1 = sin(theta_lines(2,1));
cos_1 = cos(theta_lines(2,1));
sin_2 = sin(theta_lines(2,2));
cot_1 = cot(theta_lines(2,1));
cos_2 = cos(theta_lines(2,2));

y_cross_rear = (R2 - (R1 * sin_2/sin_1 ))/(-cot_1*sin_2+cos_2);
x_cross_rear = -(y_cross_rear*cos_2 - R2)/sin_2;

%%
% get two additional points (one for each line) for x = x_cross_front - 0.5
sec_point_2_1 = [0 0];
sec_point_2_1(2) = 2;
sec_point_2_1(1) = (R1 - sec_point_2_1(2) * cos_1)/sin_1;

sec_point_2_2 = sec_point_2_1;
sec_point_2_2(2) = 2;
sec_point_2_2(1) = (R2 - sec_point_2_2(2) * cos_2)/sin_2;
%%
% print cross point (o) and additional points (x)
hold on;
plot(x_cross_rear,y_cross_rear, 'o');
hold on;
plot(sec_point_2_1(1),sec_point_2_1(2), 'x');
hold on;
plot(sec_point_2_2(1),sec_point_2_2(2), 'x');

%%
% compose matchedPoints 2  -> poits from LASER FRONT (cross, line1, line2)
matchedPoints1 = [x_cross_rear y_cross_rear ;
                  sec_point_2_1;
                  sec_point_2_2];
matchedPoints2 = [0 0 0;
                  0 0 0;
                  0 0 0];  



matchedPoints2(1,1:3) = ([x_cross_front y_cross_front 0 ]);
matchedPoints2(2,1:3) = [sec_point_1_1(1) sec_point_1_1(2) 0];
matchedPoints2(3,1:3) = [sec_point_1_2(1) sec_point_1_2(2)  0];

%%
% Rotate points by the initial rotation (180 deg)
matchedPoints2(1,1:3) = (init_rot * matchedPoints2(1,1:3)')';
matchedPoints2(2,1:3) = (init_rot * matchedPoints2(2,1:3)')';
matchedPoints2(3,1:3) = (init_rot * matchedPoints2(3,1:3)')';

plot(matchedPoints2(1:3,1),matchedPoints2(1:3,2),'*');

%%
% GET inclination factor of lines

% lines from REAR LASER
a_1_1 = ( matchedPoints1(1,2) - matchedPoints1(2,2) ) / ( matchedPoints1(1,1) - matchedPoints1(2,1));
a_1_2 = ( matchedPoints1(1,2) - matchedPoints1(3,2) ) / ( matchedPoints1(1,1) - matchedPoints1(3,1));
% lines from FROT LASER
a_2_1 = ( matchedPoints2(1,2) - matchedPoints2(2,2) ) / ( matchedPoints2(1,1) - matchedPoints2(2,1));
a_2_2 = ( matchedPoints2(1,2) - matchedPoints2(3,2) ) / ( matchedPoints2(1,1) - matchedPoints2(3,1));

%%
% find corresponding lines

if abs(a_1_1 - a_2_1) < abs( a_1_1 - a_2_2)
    disp('line 1 of REAR LASER corresponds to the line 1 of FRONT LASER')
    rotation1 = atan(a_1_1) - atan(a_2_1);
    rotation2 = atan(a_1_2) - atan(a_2_2);
    
else
    disp('line 1 of REAR LASER corresponds to the line 2 of FRONT LASER')
    rotation1 = atan(a_1_1) - atan(a_2_2);
    rotation2 = atan(a_1_2) - atan(a_2_1);
end


print_rear_moved = [0 0;
                    0 0;
                    0 0];

b_1_1_moved = matchedPoints2(1,2) - a_1_1*matchedPoints2(1,1);
b_1_2_moved = matchedPoints2(1,2) - a_1_2*matchedPoints2(1,1);
y_1_moved = 2;
print_rear_moved(1,1:2) = [matchedPoints2(1,1) matchedPoints2(1,2)];
print_rear_moved(2,1:2) = [(y_1_moved-b_1_1_moved)/a_1_1 y_1_moved];
print_rear_moved(3,1:2) = [(y_1_moved-b_1_2_moved)/a_1_2 y_1_moved];

line([matchedPoints2(1,1) matchedPoints2(2,1)],[ matchedPoints2(1,2) matchedPoints2(2,2)], 'Color', 'blue', 'LineStyle' , '--');
line([matchedPoints2(1,1) matchedPoints2(3,1)],[ matchedPoints2(1,2) matchedPoints2(3,2)], 'Color', 'black', 'LineStyle', '--');
                             
line([print_rear_moved(1,1) print_rear_moved(2,1)],[ print_rear_moved(1,2) print_rear_moved(2,2)], 'Color', 'green', 'LineStyle' , '--');
line([print_rear_moved(1,1) print_rear_moved(3,1)],[ print_rear_moved(1,2) print_rear_moved(3,2)], 'Color', 'red', 'LineStyle', '--');
%%
% Get additional rotation value [rad]( rotation between lasers = init + rotation)
if average == 1
    rotation = (rotation1+rotation2)/2;
elseif average == 0 
    if align == 1 
          rotation = rotation1;
    elseif align == 2
          rotation = rotation2;
    end
end
%%
% Get distance between lasers

% rotate points from the FRONT LASER by additional angle
additional_rot = rotz(rad2deg(rotation));
matchedPoints2(1,1:3) = (additional_rot * matchedPoints2(1,1:3)')';
matchedPoints2(2,1:3) = (additional_rot * matchedPoints2(2,1:3)')';
matchedPoints2(3,1:3) = (additional_rot * matchedPoints2(3,1:3)')';

% find distances
dist_x = matchedPoints2(1,1) - x_cross_rear;
dist_y = matchedPoints2(1,2) - y_cross_rear;

% apply distances
matchedPoints2 = matchedPoints2 - [dist_x dist_y 0;
                                   dist_x dist_y 0;
                                   dist_x dist_y 0];

fprintf('\n \n')
disp('COMPUTED:')
fprintf('   X: %g [m]',dist_x )
fprintf('   Y: %g [m]',dist_y )
fprintf('   theta: %g [rad]',rotation )
fprintf('\n')
fprintf('\n')
disp('LASER front:')
fprintf('   X: %g [m]',0 )
fprintf('   Y: %g [m]',0 )
fprintf('   theta: %g [rad]',0 )
fprintf('\n')
disp('LASER rear:')
fprintf('   X: %g [m]',-dist_x )
fprintf('   Y: %g [m]',-dist_y )
fprintf('   theta: %g [rad]',deg2rad(180)+rotation )
fprintf('\n')

if strcmp(file_format,'yaml')
  %%
  % write yaml file

  fileID = fopen('velmobil_laser_poses.yaml','w');
  fprintf(fileID,'front_laser_pose:\n');
  fprintf(fileID,'   x: %g\n',0);
  fprintf(fileID,'   y: %g\n',0);
  fprintf(fileID,'   theta: %g\n',0);

  % quat = SpinCalc('EA321toQ',[rotation/2 0 0 ], 0.001,1);
  % fprintf(fileID,'   orientation:\n');
  % fprintf(fileID,'      x: %g\n', quat(1));
  % fprintf(fileID,'      y: %g\n', quat(2));
  % fprintf(fileID,'      z: %g\n', quat(3));
  % fprintf(fileID,'      w: %g\n', quat(4));

  fprintf(fileID,'rear_laser_pose:\n');
  fprintf(fileID,'   x: %g\n',-dist_x);
  fprintf(fileID,'   y: %g\n',-dist_y);
  fprintf(fileID,'   theta: %g\n', deg2rad(180)+rotation);

  % quat = SpinCalc('EA321toQ',[deg2rad(180)+rotation/2 0 0 ], 0.001,1);
  % fprintf(fileID,'   orientation:\n');
  % fprintf(fileID,'      x: %g\n', quat(1));
  % fprintf(fileID,'      y: %g\n', quat(2));
  % fprintf(fileID,'      z: %g\n', quat(3));
  % fprintf(fileID,'      w: %g\n', quat(4));
  fclose(fileID);

elseif strcmp(file_format,'xacro')
  %%
  % write xacro file

  fileID = fopen('velmobil_laser_poses.xacro','w');
  fprintf(fileID,['<xacro:property name=\"front_laser_origin_from_calibration\"> \n' ...
                  '   <origin xyz=\"0 0 0\" rpy=\"0 0 0\" /> \n' ...
                  '</xacro:property> \n']);
  fprintf(fileID,['<xacro:property name=\"rear_laser_origin_from_calibration\"> \n' ...
                  '   <origin xyz=\"%g %g 0\" rpy=\"0 0 %g\" /> \n' ...
                  '</xacro:property> \n',-dist_x, -dist_y, deg2rad(180)-rotation]);
  fclose(fileID);
end
%%
%  TEST

% show FRONT LASER maching points after transformations
plot(matchedPoints2(1:3,1),matchedPoints2(1:3,2),'*');
% draw lines defined by above points
line([matchedPoints2(1,1) matchedPoints2(2,1)],[ matchedPoints2(1,2) matchedPoints2(2,2)], 'Color', 'blue');
line([matchedPoints2(1,1) matchedPoints2(3,1)],[ matchedPoints2(1,2) matchedPoints2(3,2)], 'Color', 'black');
