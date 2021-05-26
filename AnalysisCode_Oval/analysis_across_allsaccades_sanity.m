function [params_boot,data,accuracy,num_image,oval_sig,...
    num_trial,mean_bin_index,prob,mean_bin_index_random,prob_random,acc_evidence,choice_in_fav,random_landing,fixations_per_trial]...
    = analysis_across_allsaccades_sanity(subjectID, expt_type, boot_n, params_boot, bin,fix_cond,fix_cluster_dist,peripheryPKbound,bnd)
% initialize
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');
% loading data
[data,~] = LoadAllSubjectData(subjectID,expt_type,datadir);

%note that number images is 13, some images did not get displayed, taken care of by shown_pos
choice = data.choice(1:data.current_trial);
accuracy = data.accuracy(1:data.current_trial);
trial_index = zeros(1,data.current_trial);%data.ratio>=0.5 & data.ratio<=0.65;
shown_pos = [1 2 3 4 5 6 7 8 9 10 11 13 14];
for trial=1:data.current_trial
    if abs(sum(data.frame_categories(trial,shown_pos)==1)/length(shown_pos)-min(bnd(1),bnd(2)))<=abs(bnd(1)-bnd(2))
        trial_index(trial) = 1;
    end
end
num_trial= sum(trial_index);
num_image = length(shown_pos);%data.number_of_images(1);
trial = 1;
for tr=1:length(trial_index)
    if trial_index(tr)==1
        %trying to get fixations from eyelink
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
        duration([null_index])=[];
        location_x([null_index])=[];
        location_y([null_index])=[];
        %concatenating the saccade information (locations, onset time and so on)
        ind_remove = [];
        for x=1:length(location_x)-1
            dist_val = sqrt(((location_x(x+1)-location_x(x))^2 + (location_y(x+1)-location_y(x))^2));
            if dist_val<fix_cluster_dist % pixel distance between center of two ellipses
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
        else
            fixations_per_trial(trial) = length(fix_x);
        end
        
        %trying to visualize saccades selected
        %     figure();
        %     subplot(1,2,1)
        %     plot(data.eye_tracker_points{trial}.events(2,:),'-ok');
        %     subplot(1,2,2)
        %     plot(oval_location(1,:),oval_location(2,:),'rx');
        %     hold on;
        %     plot(location_x,location_y,'-o');
        %     hold on;
        %     plot(fix_x,fix_y,'--og','linewidth',2);
        %     xlim([0 1920])
        %     ylim([0 1080])
        %     pause;
        %     close all;
        trial = trial + 1;
    end
end
max_fix = max(0,max(fixations_per_trial));
sorted_oval_all = zeros(num_trial,num_image,max_fix);

%sorting by distance
trial = 1;
for tr = 1:length(trial_index)
    if trial_index(tr)==1
        signal = squeeze(data.ideal_frame_signals(trial,shown_pos));
        for i=1:fixations_per_trial(trial)
            
            if strcmp(fix_cond,'saccade_code')
                [~,order] = sort(compute_distance(oval_location,saccade_code.fixations{trial}(i,:)));
            else
                [~,order] = sort(compute_distance(oval_location,eyelink.fixations{trial}(i,:)));
            end
            sorted_oval_all(trial,1:peripheryPKbound,i) = signal(order(1:peripheryPKbound));
        end
        trial = trial + 1;
    end
end

oval_sig = reshape(sorted_oval_all,num_trial,max_fix*num_image);

%Obtaining bias in saccade
[mean_bin_index,prob,mean_bin_index_random,prob_random,acc_evidence,choice_in_fav,random_landing] = compute_extCB_all_saccades_sanity(data,fixations_per_trial,sorted_oval_all,bin,boot_n,params_boot,trial_index,peripheryPKbound);


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