function [current_position_pixel] = deg2pixel(current_position_deg,center,varargin)
%Screen_width_pixel
swp = 1920;
%Screen_width_inch
swi = 21;
%Screen_distance_inch
sdi = 43;

% Interpret the user parameters
k = 1;
while k <= length(varargin) && ischar(varargin{k})
    switch (lower(varargin{k}))
        case 'screen_width_pixel'
            swp = varargin{k+1};
        case 'screen_width_inch'
            sw = varargin{k+1};
        case 'screen_distance_inch'
            sdi = varargin{k+1};
        otherwise
            error( 'saccades = GetSacs( eyeTrace, varargin )' );
    end
    k = k + 2;
end
if size(current_position_deg,2)==2
%Pixel per inch
ppi = swp/swi;
center_rep=repmat(center,size(current_position_deg,1),1);
%Position distance between the center and the current position in pixel
%Computing the current position in degrees given pixel per inch,
%position distance and screen distance.
pdp=tand(current_position_deg)*sdi*ppi;
current_position_pixel=[pdp(:,1)+center_rep(:,1) pdp(:,2)+center_rep(:,2)];
end
end