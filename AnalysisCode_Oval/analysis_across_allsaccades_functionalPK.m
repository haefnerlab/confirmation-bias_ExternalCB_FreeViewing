function [params_boot,best_hprs,data,accuracy,num_image,oval_sig,...
    num_trial,num_trial_sanity,acc_evidence,choice_in_fav,random_landing,...
    acc_evidence_sanity,choice_in_fav_sanity,random_landing_sanity,...
    fixations_per_trial,eyelink,saccade_code,fixation_dist,saccades_used,saccades_used_sanity,saccade_dist_mid,...
    log_bernoulli, spatial_pk, temporal_fixation_pk, temporal_duration_pk, temporal_pk_across_trials]...
    = analysis_across_allsaccades_functionalPK(subjectID, expt_type, boot_n,...
    fix_cond,saccade_types,saccade_types_sanity,fix_cluster_dist,peripheryPKbound,sanity_bnd,hprs_ridge)

% initialize
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');
% loading data
[data,~] = LoadAllSubjectData(subjectID,expt_type,datadir);
disp('Data for the subject loaded!');

%note that number images is 13, some images did not get displayed, taken care of by shown_pos
choice = data.choice(1:data.current_trial);
accuracy = data.accuracy(1:data.current_trial);
num_trial= data.current_trial;
shown_pos = [1 2 3 4 5 6 7 8 9 10 11 13 14];
num_image = length(shown_pos);%data.number_of_images(1);
fixation_dist_all = [];
bins_pk = [1920 1080];
trial_total_length = data.stim_duration*1000; % in ms.
for trial=1:num_trial
    %trying to get fixations from eyelink!!!!
    fix_x = [];
    fix_y = [];
    fix_dur = [];
    fix_pupil_size = [];
    marker = 0;
    count = 0;
    temp_x = 0;
    temp_y = 0;
    temp_dur = 0;
    temp_pupil_size = 0;
    flag = 0;
    for ev=1:data.eye_tracker_points{trial}.eventNum
        if ev==1 && data.eye_tracker_points{trial}.events(2,ev)==9
            flag = 1;
        end
        if data.eye_tracker_points{trial}.events(2,ev)==7
            % disp(['X cordinate at start of fixation: ' num2str(data.eye_tracker_points{trial}.events(19,ev))]);
            flag = 1;
        end
        %         if data.eye_tracker_points{trial}.events(2,ev)==8
        %             flag = 0;
        %         end
        if flag==1 && data.eye_tracker_points{trial}.events(2,ev)~=7
            count = count + 1;
            marker = 1;
            temp_x = temp_x + data.eye_tracker_points{trial}.events(19,ev);
            temp_y = temp_y + data.eye_tracker_points{trial}.events(20,ev);
            temp_dur = temp_dur + (data.eye_tracker_points{trial}.events(6,ev) - data.eye_tracker_points{trial}.events(5,ev));
            if count==1
                temp_onset  = data.eye_tracker_points{trial}.events(5,ev);
            end
            temp_pupil_size = temp_pupil_size + data.eye_tracker_points{trial}.events(21,ev);
        end
        if flag==0 || ev==data.eye_tracker_points{trial}.eventNum
            if marker==1
                fix_x(end+1) = temp_x/count;
                fix_y(end+1) = temp_y/count;
                fix_dur(end+1) = temp_dur/count;
                fix_onset = temp_onset;
                fix_pupil_size(end+1) = temp_pupil_size/count;
                marker = 0;
                count = 0;
                temp_x = 0;
                temp_y = 0;
                temp_dur = 0;
                temp_onset = 0;
                temp_pupil_size = 0;
            end
        end
        if data.eye_tracker_points{trial}.events(2,ev)==8
            flag = 0;
        end
    end
    %   disp(['Number of fixations: ' num2str(length(fix_x))]);
    
    eyelink.fixations{trial} = [fix_x' fix_y'];
    eyelink.fixations_duration{trial} = fix_dur;
    eyelink.fixations_onset_time_stamps{trial} = fix_onset/trial_total_length;
    eyelink.fixations_pupil_size{trial} = fix_pupil_size;
    
    % trying to get fixations from saccade detection code!!!!
    %getting duration of fixations
    duration = data.eye_tracker_points{trial}.duration;
    % subtract initial timestamp to ensure time stamps within 1500 ms
    onset_time_stamps = data.eye_tracker_points{trial}.onset_offset(:,1)- data.eye_tracker_points{trial}.time_stamp(1);
    %getting locations of saccades
    location_x = data.eye_tracker_points{trial}.xylocations(:,1)';
    location_y = data.eye_tracker_points{trial}.xylocations(:,2)';
    %getting stimulus location
    oval_location = [data.eye_tracker_points{trial}.locationtodraw(shown_pos,1)'; data.eye_tracker_points{trial}.locationtodraw(shown_pos,2)'];
    %getting signal
    signal = squeeze(data.ideal_frame_signals(trial,shown_pos));
    %eliminating null saccade data
    null_index = [];
    for i=1:length(location_x)
        if duration(i)<0
            null_index(end+1) = i;
        end
        if location_x(i)<0 || location_y(i)<0
            null_index(end+1) = i;
        end
    end
    location_x_pre = location_x;
    location_y_pre = location_y;
    duration([null_index])=[];
    location_x([null_index])=[];
    location_y([null_index])=[];
    onset_time_stamps([null_index])=[];
    %concatenating the saccade information (locations, onset time and so on)
    ind_remove = [];
    for x=1:length(location_x)-1
        dist_val = sqrt(((location_x(x+1)-location_x(x))^2 + (location_y(x+1)-location_y(x))^2));
        if dist_val< fix_cluster_dist % based on pixels distance between center of two ellipses
            ind_remove(end+1) = x+1;
        end
    end
    duration([ind_remove])=[];
    location_x([ind_remove])=[];
    location_y([ind_remove])=[];
    onset_time_stamps([ind_remove])=[];
    saccade_code.fixations{trial} = [location_x' location_y'];
    saccade_code.fixations_duration{trial} = duration;
    saccade_code.fixations_onset_time_stamps{trial} = onset_time_stamps/trial_total_length;
    
    if strcmp(fix_cond,'saccade_code')
        fixations_per_trial(trial) = length(location_x);
        for ln=1:fixations_per_trial(trial)-1
            fixation_dist{trial}(ln) = sqrt(((saccade_code.fixations{trial}(ln+1,1)-saccade_code.fixations{trial}(ln,1))^2 + (saccade_code.fixations{trial}(ln+1,2)-saccade_code.fixations{trial}(ln,1))^2));
            fixation_dist_all(end+1) = fixation_dist{trial}(ln);
        end
    else
        fixations_per_trial(trial) = length(fix_x);
        for ln=1:fixations_per_trial(trial)-1
            fixation_dist{trial}(ln) = sqrt(((eyelink.fixations{trial}(ln+1,1)-eyelink.fixations{trial}(ln,1))^2 + (eyelink.fixations{trial}(ln+1,2)-eyelink.fixations{trial}(ln,1))^2));
            fixation_dist_all(end+1) = fixation_dist{trial}(ln);
        end
    end
    %trying to visualize saccades selected
    % %         figure();
    %         subplot(1,4,4)
    %         plot(data.eye_tracker_points{trial}.events(2,:),'-ok','linewidth',2);
    %         axis image;
    %         title('Event points of eyelink')
    %
    %         subplot(1,4,1)
    %         plot(oval_location(1,:),oval_location(2,:),'rx','linewidth',4);
    %         hold on;
    %         plot(location_x_pre,location_y_pre,'-ob','linewidth',2);
    %         hold on;
    %         axis image;
    %         xlim([0 1920])
    %         ylim([0 1080])
    %         title(['Saccade from code before prune(' num2str(length(location_x_pre)) ' fixations)'])
    %
    %         subplot(1,4,2)
    %         plot(oval_location(1,:),oval_location(2,:),'rx','linewidth',4);
    %         hold on;
    %         plot(location_x,location_y,'-ob','linewidth',2);
    %         hold on;
    %         axis image;
    %         xlim([0 1920])
    %         ylim([0 1080])
    %         title(['Saccade from code(' num2str(length(location_x)) ' fixations)'])
    %
    %         subplot(1,4,3)
    %         plot(oval_location(1,:),oval_location(2,:),'rx','linewidth',4);
    %         hold on;
    %         plot(fix_x,fix_y,'-ok','linewidth',2);
    %         xlim([0 1920])
    %         ylim([0 1080])
    %         title(['Eyelink points(' num2str(length(fix_x)) ' fixations)'])
    %         pause;
    %         close all;
end

for ss=1:length(saccade_types)
    saccade_dist_edges{ss} = linspace(min(fixation_dist_all),max(fixation_dist_all),saccade_types(ss)+1);
    for i=1:length(saccade_dist_edges{ss})-1
        saccade_dist_mid{ss}(i) = (saccade_dist_edges{ss}(i) + saccade_dist_edges{ss}(i+1))/2;
    end
end
saccade_dist_edges_sanity = linspace(min(fixation_dist_all),max(fixation_dist_all),saccade_types_sanity+1);
disp('Saccades obtained and other preprocessing done!');

disp('Computing saccade distances and learning functional weights ...');
max_fix = max(fixations_per_trial);
sorted_oval_all = zeros(num_trial,num_image,max_fix);
sorted_oval_location = zeros(num_trial,2,num_image,max_fix);
fixation_onset_time_stamps = zeros(num_trial,max_fix);
%sorting by distance
for trial = 1:num_trial
    signal = squeeze(data.ideal_frame_signals(trial,shown_pos));
    for i=1:fixations_per_trial(trial)
        if strcmp(fix_cond,'saccade_code')
            [peri_stim_distance{trial}(i,:),order] = sort(compute_distance(oval_location,saccade_code.fixations{trial}(i,:)));
        else
            [peri_stim_distance{trial}(i,:),order] = sort(compute_distance(oval_location,eyelink.fixations{trial}(i,:)));
        end
        sorted_oval_all(trial,1:peripheryPKbound,i) = signal(order(1:peripheryPKbound));
        if strcmp(fix_cond,'saccade_code')
            sorted_oval_distance(trial,:,1:peripheryPKbound,i) = squeeze(oval_location(:,order(1:peripheryPKbound)))-squeeze(saccade_code.fixations{trial}(i,:))';
            fixation_onset_time_stamps(trial,i) = saccade_code.fixations_onset_time_stamps{trial}(i);
        else
            sorted_oval_distance(trial,:,1:peripheryPKbound,i) = squeeze(oval_location(:,order(1:peripheryPKbound)))-squeeze(eyelink.fixations{trial}(i,:))';  
            fixation_onset_time_stamps(trial,i) = eyelink.fixations_onset_time_stamps{trial}(i);
        end
    end
end
sorted_oval_distance(:,1,:,:) = sorted_oval_distance(:,1,:,:)/1920;
sorted_oval_distance(:,2,:,:) = sorted_oval_distance(:,1,:,:)/1080;
oval_sig = reshape(sorted_oval_all,num_trial,max_fix*num_image);
ordered_oval_distance = reshape(sorted_oval_distance,num_trial,2,max_fix*num_image);

disp('Searching ridge based hyperparameters ...');
[best_hprs, ~] = CustomRegression.xValidatePK_with_lapse_functionalPK(oval_sig, choice, ordered_oval_distance, fixation_onset_time_stamps, max_fix, hprs_ridge, 0, 10 );
% best_hprs = 0;
disp(['Best hyperparameters found: ' num2str(best_hprs)]);
%Analysis of PK after cross validation on all trials with non empty bins
disp('Computing the kernels...');
for j = 1:boot_n
    [signal_result(j,:,:), choice_result(j,:),ordered_oval_distance_result(j,:,:,:),fixation_onset_time_stamps_result(j,:,:)] = bootstrap(oval_sig, choice, ordered_oval_distance, fixation_onset_time_stamps, num_trial);
end
for j = 1:boot_n
    if mod(j,10)==0 || j==1
        disp(['Running bootstrap step ' num2str(j) '/' num2str(boot_n) ' ...']);
    end
    [params_boot(j,:)] = CustomRegression.PsychophysicalKernelwithlapse_functionalPK(squeeze(signal_result(j,:,:)), squeeze(choice_result(j,:)), squeeze(ordered_oval_distance_result(j,:,:,:)), squeeze(fixation_onset_time_stamps_result(j,:,:)), max_fix, best_hprs, 0);
end

disp('Computing bias in saccade...');
%Obtaining bias in saccade
for ss=1:length(saccade_types)
    [acc_evidence{ss},choice_in_fav{ss},random_landing{ss},saccades_used{ss}] = ...
        compute_extCB_all_saccades_functionalPK(data,fixations_per_trial,sorted_oval_all,...
        saccade_dist_edges{ss},ordered_oval_distance,fixation_onset_time_stamps,...
        params_boot,fixation_dist,peripheryPKbound);
end
disp('Computing bias in saccade for sanity check ...');
%Obtaining bias in saccade
[acc_evidence_sanity,choice_in_fav_sanity,random_landing_sanity,saccades_used_sanity,num_trial_sanity] = ...
    compute_extCB_all_saccades_sanity_functionalPK(data,fixations_per_trial,sorted_oval_all,...
    saccade_dist_edges_sanity,ordered_oval_distance, fixation_onset_time_stamps, ...
    params_boot,fixation_dist,peripheryPKbound,sanity_bnd);

disp('Computing log odds, spatial kernel and temporal kernel ...')
[log_bernoulli, spatial_pk, temporal_fixation_pk, temporal_duration_pk, temporal_pk_across_trials] = compute_log_odds_and_pks(oval_sig, ordered_oval_distance, fixation_onset_time_stamps, max_fix, params_boot, bins_pk, data.stim_duration);

    
    
    function [signals, choices, ordered_oval_distances, onset_time_stamps] = bootstrap(signals_raw, choices_raw, ordered_oval_distance_raw, fixation_onset_time_stamps_raw, trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signals = [];
        choices = [];
        for z = 1:trials
            trial_num = sample_nums(z);
            signals = [signals; signals_raw(trial_num,:)];
            choices = [choices choices_raw(trial_num)];
            ordered_oval_distances(z,:,:) = squeeze(ordered_oval_distance_raw(trial_num,:,:));
            onset_time_stamps(z,:) = fixation_onset_time_stamps_raw(z,:);
        end
    end

    function  [distance] = compute_distance (xy_point1,xy_point2)
        distance=[];
        for k=1:length(xy_point1)
            distance(end+1) = sqrt((xy_point1(1,k)-xy_point2(1))^2 + (xy_point1(2,k)-xy_point2(2))^2);
        end
    end

end