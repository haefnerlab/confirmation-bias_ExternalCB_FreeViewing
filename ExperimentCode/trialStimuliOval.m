function [image_properties,locations,tracker_info, quit] = trialStimuliOval(GaborData, wPtr, tracker_info,settings)
KbName('UnifyKeyNames');
exitKey = KbName(settings.keyExit);
leftKey = KbName(settings.keyLeft);
rightKey = KbName(settings.keyRight);
quit = false;
tracker_info = Eyelink.startTrialPrepare(tracker_info);
eye_tracker_points = [];
% broke_fixation = false;
black = [0 0 0];
gray = [127 127 127];
white = [255 255 255];
center=[settings.screenSize(3)/2 settings.screenSize(4)/2]; % Get the side of the vertical axis


xcenter = settings.screenSize(3)/2; % Get the middle of the horizontal axis
ycenter = settings.screenSize(4)/2; % Get the middle of the vertical axis 
stimulus_base_horiz = ptbCenteredRect([xcenter - 350 ycenter], [120 40]); % show choice 
stimulus_base_vert = ptbCenteredRect([xcenter + 350 ycenter], [40 120]);
stimulus = GaborData.contrast(GaborData.current_trial) * black;

[~,~,x_mid_in,y_mid_in,x_mid_out,y_mid_out] = Hex_grid_generator(GaborData); % determine locations 
A = [x_mid_in; x_mid_out]; B = [y_mid_in; y_mid_out];
locations = [A, B];

Screen('FillRect', wPtr, gray);
% drawTrialNo();
% Screen('DrawLine', wPtr ,white, xcenter, ycenter-15, xcenter, ycenter+15 ,3);
% Screen('DrawLine', wPtr ,white, xcenter-15, ycenter, xcenter+15, ycenter ,3);
tracker_info = Eyelink.startTrial(tracker_info);
Eyelink.drawFixationSymbol(tracker_info,xcenter,ycenter, wPtr);
Eyelink.clearBuffer(tracker_info.drained);
[~, onset_fixation] = Screen('Flip', wPtr);
WaitSecs(GaborData.prefixation_duration);
tstart = GetSecs();
count_queue_drawn = 0;
while GetSecs() - tstart < GaborData.fixation_duration
    Eyelink.drawFixationSymbol(tracker_info,xcenter,ycenter, wPtr);
    [~, onset_fixation] = Screen('Flip', wPtr);
    tracker_info = Eyelink.getQueue(tracker_info);
    [eyex, eyey , ~] = Eyelink.getGazePoint(tracker_info);
    [Fixation]=isFixation(tracker_info.pre_fixationRadius, [eyex eyey],[xcenter,ycenter], 0);
    if Fixation ~= 1
    image_properties.choice = nan;
    return;
    end
    
end
% Prep for first stimulus frame by clearing the drawStimulusFrame.

% Plot a hexagonal grid

for i=1:GaborData.number_of_images   % we need number of images
    stimulus_bbox = ptbCenteredRect([locations(i,1) locations(i,2)], [GaborData.xaxis(GaborData.current_trial,i) GaborData.yaxis(GaborData.current_trial,i)]);
    Screen('FillOval',wPtr,stimulus,stimulus_bbox);
end



disp('Stim Starts');
Eyelink.clearBuffer(tracker_info.drained);
[~, start_time] = Screen('Flip', wPtr, onset_fixation + GaborData.fixation_duration);
tstart = GetSecs();
count_queue_drawn = 0;
while GetSecs() - tstart < GaborData.stim_duration
    tracker_info = Eyelink.getQueue(tracker_info);
    count_queue_drawn = count_queue_drawn + 1;
    [gx, gy, time_stamp(count_queue_drawn)] = Eyelink.getGazePoint(tracker_info);
    hist_check(:,count_queue_drawn) = [gx gy]';
    for i=1:GaborData.number_of_images   % we need number of images
        stimulus_bbox = ptbCenteredRect([locations(i,1) locations(i,2)], [GaborData.xaxis(GaborData.current_trial,i) GaborData.yaxis(GaborData.current_trial,i)]);
        Screen('FillOval',wPtr,stimulus,stimulus_bbox);
    end
    %Eyelink.drawFixationSymbol(tracker_info,gx,gy, wPtr);
    Screen('Flip', wPtr);
end

tracker_info.count_queue_drawn = count_queue_drawn;
tracker_info.hist_check = hist_check;
tracker_info.time_stamp = time_stamp;

% Computing a saccade information.
hist_check_deg=pixel2deg(hist_check',center);
[~, fixations] = simpleSaccadeDetect(time_stamp, hist_check_deg(:,1), hist_check_deg(:,2));
%storing the relevent values
fixation_info=struct2cell(fixations);
onset_offset=squeeze(cell2mat(fixation_info(1:2,:,:)))';
duration=squeeze(cell2mat(fixation_info(3,:,:)))';
if size(cell2mat(fixation_info(4,:,:)),3)>1
xy_locations=deg2pixel(squeeze([cell2mat(fixation_info(4,:,:)) cell2mat(fixation_info(5,:,:))])',center);
else
xy_locations=deg2pixel(squeeze([cell2mat(fixation_info(4,:,:)) cell2mat(fixation_info(5,:,:))]),center);
end
tracker_info.onset_offset=onset_offset;
tracker_info.duration=duration;
tracker_info.xylocations=xy_locations;
tracker_info.locationtodraw=locations;
Screen('FillRect', wPtr, gray);
[~, start_blank] = Screen('Flip', wPtr, start_time + GaborData.stim_duration);
Screen('FillOval', wPtr, black, stimulus_base_horiz); %horizontal
Screen('FillOval', wPtr, black, stimulus_base_vert); %vertical
Screen('Flip', wPtr, start_blank + GaborData.go_cue_time);
[key, rt, timeout] = ptbWaitKey([leftKey, rightKey, exitKey], 1.8);

% Close textures to avoid memory problems.
% for i = 1:total_frames
%     Screen('Close', image_texture(i));
% end
% Screen('Close', show_left_patch);
% Screen('Close', show_right_patch);

if key == exitKey
    quit = true;
    image_properties.choice = nan;
end

if timeout
    image_properties.choice = nan;
else
    image_properties.reaction = rt * 1000;
    if key == leftKey
        image_properties.choice = 0;
    elseif key == rightKey
        image_properties.choice = 1;
    end
end