function genCircles(is_template,contrast,wPtr,w,xc,yc,category)
% diff = 100;
% shft = w/3;
% gray = [127 127 127];
white = [255 255 255];
black = [0 0 0];
% diff_color = 255-127;
number_of_images = 10;
dark = 127-127*contrast;
light = 127+127*contrast;
for i=1:number_of_images
    ax=w;
    bx=2*xc-w;
    t1 = (bx-ax).*rand(1) + ax;
    
    locations(i,1) = t1;
    
    ay=w;
    by=2*yc-w;
    t2 = (by-ay).*rand(1) + ay;
    
    locations(i,2) = t2;
    
    %%prevent overlap
end
if is_template==1
    stimulus_bbox_dark = ptbCenteredRect([xc+w yc], [w w]);
    stimulus_bbox_light = ptbCenteredRect([xc-w yc], [w w]);
    Screen('FillOval', wPtr, white ,stimulus_bbox_dark )
    Screen('FillOval',  wPtr, black, stimulus_bbox_light);
    
    
else
    if category==-1
        for i=1:number_of_images
            stimulus_bbox = ptbCenteredRect([locations(i,1) locations(i,2)], [w w]);
            Screen('FillOval',wPtr,light,stimulus_bbox);
        end
    else
        for i=1:number_of_images
            stimulus_bbox = ptbCenteredRect([locations(i,1) locations(i,2)], [w w]);
            Screen('FillOval',wPtr,dark,stimulus_bbox);
        end
    end
    
    
    
end

end