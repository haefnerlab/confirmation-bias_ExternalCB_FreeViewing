function [hex_X,hex_Y,x_mid_in,y_mid_in,x_mid_out,y_mid_out] = Hex_grid_generator(GaborData)

% , dd, GaborData.last_X,GaborData.last_Y

%% Function to Create Hexagonal Grid
%% Author: Tauseef Gulrez, Footscray, Melbourne, Australia, 09-06-2018

% "Measuring Behaviour as Perceptual Control Using Human Simulations of
% Predator and Prey Pursuit." SS, TG, MP, WM.
% Measuring Behavior 2018, Manchester, UK.

% Inputs:
% last_x = last coordinate point on x-axis for last grid in x-direction,
% last_y = last coordinate point on y-axis for last grid in y-direction,

% Outputs:
% hex_X,hex_Y - Corner points of hexagons
% x_mid_in,y_mid_in - Mid (x,y) points inside of Hexagon
% x_mid_out,y_mid_out - Mid (x,y) points outside of Hexagon

% Example:
% [hex_X,hex_Y,x_mid_in,y_mid_in,x_mid_out,y_mid_out] = Hex_grid_generator(20,30);

% Copyright (c) 2018, Tauseef Gulrez



last_X = GaborData.last_X;
last_Y = GaborData.last_Y;

point_X = GaborData.point_X;
point_Y = GaborData.point_Y;

dd = GaborData.grid_length;
hex_X = [];
hex_Y = [];
x_mid_in = [];
y_mid_in = [];
x_mid_out = [];
y_mid_out = [];
hex_px = [];
hex_py = [];
k = 0;

for j = 1:1:length(point_X)
    for q = 1:1:length(point_Y)
        
        hex_X_temp = [];
        hex_Y_temp = [];
        
        kk = 0;
        
        for i = [60,120,180,240,300,360,60]
            
            k = k+1;
            kk = kk+1;
            
            hex_01_x = dd * cosd(i);
            hex_01_y = dd * sind(i);
            
            hex_X(k,:) = point_X(j) + hex_01_x;         %kth row of matrix hex_X
            hex_Y(k,:) = point_Y(q) + hex_01_y;
            
            hex_X_temp(kk,:) = point_X(j) + hex_01_x;
            hex_Y_temp(kk,:) = point_Y(q) + hex_01_y;
            
            
        end
        % Middle Point
        % Let us Make two lines:
        % Line 1 (first and fourth point of hexagon - in counterclockwise)
        x1 = [hex_X_temp(1) hex_X_temp(4)];
        y1 = [hex_Y_temp(1) hex_Y_temp(4)];
        
        % Line2 (second and fifth point of hexagon - in counterclockwise)
        x2 = [hex_X_temp(2) hex_X_temp(5)];
        y2 = [hex_Y_temp(2) hex_Y_temp(5)];
        
        %fit linear polynomial
        p1 = polyfit(x1,y1,1);
        p2 = polyfit(x2,y2,1);
        
        %calculate intersection
        x_intersect = fzero(@(x) polyval(p1-p2,x),3);
        y_intersect = polyval(p1,x_intersect);
        
        x_mid_in = [x_mid_in ; x_intersect];
        y_mid_in = [y_mid_in ; y_intersect];

        % Let us Make two more lines:
        % Line 3 (first and second point of hexagon - in counterclockwise)
        x3 = [hex_X_temp(1) hex_X_temp(2)];
        y3 = [hex_Y_temp(1) hex_Y_temp(2)];
        
        % Line 4 (fifth and sixth point of hexagon - in counterclockwise)
        x4 = [hex_X_temp(6) hex_X_temp(5)];
        y4 = [hex_Y_temp(6) hex_Y_temp(5)];
        
        %fit linear polynomial
        p3 = polyfit(x3,y3,1);
        p4 = polyfit(x4,y4,1);
        
        %calculate intersection
        x_mid = fzero(@(x) polyval(p3-p4,x),3);
        y_mid = polyval(p3,x_mid);
        
        x_mid_out = [x_mid_out ; x_mid];
        y_mid_out = [y_mid_out ; y_mid];
        
        % Save all 360 Points
        hex_px = [hex_X_temp(6) ; hex_px];
        hex_py = [hex_Y_temp(6) ; hex_py];
    end
end

end