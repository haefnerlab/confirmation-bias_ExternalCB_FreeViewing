function [params_boot,sobl,abbl,best_hprs,data,accuracy,num_image,oval_sig,...
    num_trial,mean_bin_index,prob,mean_bin_index_random,prob_random,acc_evidence,choice_in_fav,random_landing,fixations_per_trial,fixations_cases,eyelink,saccade_code,saccade_length,fixation_dist]...
    = analysis_across_allsaccades(subjectID, expt_type, boot_n, hpr1, hpr2, bin,fix_cond,saccade_types,fix_cluster_dist)

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
            %             disp(['X cordinate at start of fixation: ' num2str(data.eye_tracker_points{trial}.events(19,ev))]);
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
            temp_pupil_size = temp_pupil_size + data.eye_tracker_points{trial}.events(21,ev);
        end
        if flag==0 || ev==data.eye_tracker_points{trial}.eventNum
            if marker==1
                fix_x(end+1) = temp_x/count;
                fix_y(end+1) = temp_y/count;
                fix_dur(end+1) = temp_dur/count;
                fix_pupil_size(end+1) = temp_pupil_size/count;
                marker = 0;
                count = 0;
                temp_x = 0;
                temp_y = 0;
                temp_dur = 0;
                temp_pupil_size = 0;
            end
        end
        if data.eye_tracker_points{trial}.events(2,ev)==8
            flag = 0;
        end
    end
    %     disp(['Number of fixations: ' num2str(length(fix_x))]);
    
    eyelink.fixations{trial} = [fix_x' fix_y'];
    eyelink.fixations_duration{trial} = fix_dur;
    eyelink.fixations_pupil_size{trial} = fix_pupil_size;
    
    % trying to get fixations from saccade detection code!!!!
    %getting duration of saccades
    duration = data.eye_tracker_points{1,trial}.duration;
    %getting locations of saccades
    location_x = data.eye_tracker_points{1,trial}.xylocations(:,1)';
    location_y = data.eye_tracker_points{1,trial}.xylocations(:,2)';
    %getting stimulus location
    oval_location = [data.eye_tracker_points{1,trial}.locationtodraw(shown_pos,1)'; data.eye_tracker_points{1,trial}.locationtodraw(shown_pos,2)'];
    %getting signal
    signal = squeeze(data.ideal_frame_signals(trial,shown_pos));
    %eliminating null saccade data
    null_index = [];
    trial_total_length = data.stim_duration*1000; % in ms.
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
    saccade_code.fixations{trial} = [location_x' location_y'];
    saccade_code.fixations_duration{trial} = duration;
    
    
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
for ll=1:length(saccade_types)
saccade_length(ll) = saccade_types(ll)*max(fixation_dist_all);%median(fixation_dist_all);%0.25*max(fixation_dist_all);%
fixations_cases(ll) = sum(fixation_dist_all<=saccade_length(ll));
end
max_fix = max(fixations_per_trial);

sorted_oval_all = zeros(num_trial,max_fix,num_image);
sorted_oval_closest = zeros(num_trial,max_fix);
sorted_oval_meantwoclosest = zeros(num_trial,max_fix);
sorted_oval_twoclosest = zeros(num_trial,max_fix,2);

%sorting by distance
for trial = 1:num_trial
    signal = squeeze(data.ideal_frame_signals(trial,shown_pos));
    for i=1:fixations_per_trial(trial)
        if strcmp(fix_cond,'saccade_code')
            [~,order] = sort(compute_distance(oval_location,saccade_code.fixations{trial}(i,:)));
        else
            [~,order] = sort(compute_distance(oval_location,eyelink.fixations{trial}(i,:)));
            
        end
        sorted_oval_all(trial,i,:) = signal(order);
        sorted_oval_closest(trial,i) = signal(order(1));
        sorted_oval_twoclosest(trial,i,:) = signal(order(1:2));
%         sorted_oval_twoclosest(trial,i,1) = signal(order(1));
        sorted_oval_meantwoclosest(trial,i) = mean(signal(order(1:2)));
    end
end

oval_sig = reshape(sorted_oval_twoclosest,num_trial,max_fix*2);


%Analysis of PK after cross validation on all trials with non empty bins
[best_hprs, ~] = CustomRegression.xValidatePK_with_lapse(reshape(sorted_oval_twoclosest,num_trial,max_fix*2), choice, max_fix, hpr1, 0, hpr2, 0, 10);
disp('Best hyperparameters found!');
for j = 1:boot_n
    [signal_result, choice_result] = bootstrap(reshape(sorted_oval_twoclosest,num_trial,max_fix*2), choice, num_trial);
    [sobl(j,:), ~] = CustomRegression.LinearPK_with_lapse(signal_result, choice_result, 0);
    %     sobl(j,:) = zeros(1,4);
    if mod(j,100)==0 || j==1
        disp(['Bootstrap step ' num2str(j) ' ...']);
    end
    [params_boot(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapse(signal_result, choice_result, max_fix, best_hprs(1), 0, best_hprs(3), 0);
    [abbl(j,:), ~, ~] = CustomRegression.ExponentialPK_with_lapse(signal_result, choice_result, 0);
    %     abbl(j,:) = zeros(1,4);
end

disp('Computing bias in saccade...');
%Obtaining bias in saccade
for bn=1:length(bin)
    [mean_bin_index{bn},prob{bn},mean_bin_index_random{bn},prob_random{bn},acc_evidence{bn},choice_in_fav{bn},random_landing{bn}] = compute_extCB_all_saccades(data,fixations_per_trial,sorted_oval_twoclosest,bin(bn),boot_n,params_boot,fixation_dist,saccade_length);
end

    function [signals, choices] = bootstrap(signals_raw, choices_raw, trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signals = [];
        choices = [];
        for z = 1:trials
            trial_num = sample_nums(z);
            signals = [signals; signals_raw(trial_num,:)];
            choices = [choices choices_raw(trial_num)];
        end
    end

    function  [distance] = compute_distance (xy_point1,xy_point2)
        distance=[];
        for k=1:length(xy_point1)
            distance(end+1) = sqrt((xy_point1(1,k)-xy_point2(1))^2 + (xy_point1(2,k)-xy_point2(2))^2);
        end
    end

end