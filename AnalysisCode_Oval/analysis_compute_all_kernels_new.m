function [params_boot,sobl,abbl,best_hprs,data,accuracy,num_image,mean_oval_signal1,acc_evi_all,pro_choice_all,...
    acc_evi_ideal1,pro_choice_ideal1,acc_evi_ideal2,pro_choice_ideal2,acc_evi_ideal_random2,pro_choice_ideal_random2,acc_evi_ideal1avg,pro_choice_ideal1avg,num_trial_used,num_trial_done,tr_ratio1,choice1,fixations_per_trial]...
    = analysis_compute_all_kernels_new(subjectID, expt_type, boot_n, hpr1, hpr2, split_parameter, bin, choose_abs_PK, learn_num_PK)
% initialize
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');
% loading data
[data,~] = LoadAllSubjectData(subjectID,expt_type,datadir);
disp('Data loaded for the subject!');
%note that number images is 13, some images did not get displayed, taken
%care of by shown_pos
choice = data.choice(1:data.current_trial);
accuracy = data.accuracy(1:data.current_trial);
num_trial= data.current_trial;
num_trial_done = data.current_trial;
ind_tr = [];
ind_tr1 = [];
shown_pos = [1 2 3 4 5 6 7 8 9 10 11 13 14];
mean_oval_signal = zeros(data.current_trial,2*length(shown_pos));
num_image = length(shown_pos);%data.number_of_images(1);


for trial=1:num_trial
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
% %     disp(['Number of fixations: ' num2str(length(fix_x))]);
    
     
    %trying to look at trials within a range of ratios
    %     if abs(sum(data.frame_categories(trial,shown_pos)==1)/length(shown_pos)-0.5)<=0.1
    %getting onset_time of saccades
    onset_realtime = data.eye_tracker_points{1, trial}.onset_offset(:,1)';
    %getting duration of saccades
    duration = data.eye_tracker_points{1,trial}.duration;
    %getting locations of saccades
    location_x = data.eye_tracker_points{1,trial}.xylocations(:,1)';
    location_y = data.eye_tracker_points{1,trial}.xylocations(:,2)';
    %getting stimulus location
    oval_location = [data.eye_tracker_points{1,trial}.locationtodraw(shown_pos,1)'; data.eye_tracker_points{1,trial}.locationtodraw(shown_pos,2)'];
    %getting signal
    signal = data.ideal_frame_signals(trial,shown_pos);
    %eliminating null saccade data
    null_index = [];
    trial_total_length = data.stim_duration*1000; % in ms.
    for i=length(onset_realtime):-1:1
        %converting the time to real time
        onset_realtime(i) = onset_realtime(i)-onset_realtime(1);
        if duration(i)<0
            null_index(end+1) = i;
        end
        if location_x(i)<0 || location_y(i)<0
            null_index(end+1) = i;
        end
    end
    duration([null_index])=[];
    onset_realtime([null_index])=[];
    location_x([null_index])=[];
    location_y([null_index])=[];
    
    %concatenating the saccade informatoin (locations, onset time and so on)
    
    location_x1 = location_x;
    location_y1 = location_y;
    ind_remove = [];
            for x=1:length(location_x1)-1
                dist_val = sqrt(((location_x1(x+1)-location_x1(x))^2 + (location_y1(x+1)-location_y1(x))^2));
                if dist_val<135 % 270 pixels is the distance between center of two ellipses
                    location_x1(x+1) = location_x1(x);
                    location_y1(x+1) = location_y1(x);
                    ind_remove(end+1) = x+1;
                end
            end
    duration([ind_remove])=[];
    onset_realtime([ind_remove])=[];
    location_x([ind_remove])=[];
    location_y([ind_remove])=[];
    fixations_per_trial(trial) = length(location_x);


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
%     xlim([0 1920]);
%     ylim([0 1080]);
%     pause;
%     close all;
%     
    saccade_info = vertcat(location_x,location_y,onset_realtime,duration);
    %making time reference for each bin
    split_time_reference = [trial_total_length/split_parameter:trial_total_length/split_parameter:trial_total_length];
    %allocating the saccades to appropriate bin
    for i=1:length(split_time_reference)
        if i==1
            location_by_split{i} = [location_x(find(onset_realtime<=split_time_reference(i)));location_y(find(onset_realtime<=split_time_reference(i)))];
        else
            location_by_split{i} = [location_x(find(onset_realtime<=split_time_reference(i) & onset_realtime>split_time_reference(i-1)));location_y(find(onset_realtime<=split_time_reference(i) & onset_realtime>split_time_reference(i-1)))];
        end
    end
    %sorting by distance
    for i=1:size(location_by_split,2)
        for j=1:size(location_by_split{i},2)
            [result,order] = sort(compute_distance(oval_location,location_by_split{i}(:,j)));
            sorted_oval_information{i}{j} = [result; signal(order)];
        end
    end
    
    for i=1:size(location_by_split,2)
        sum_sig = 0;
        for j=1:size(location_by_split{i},2)
            %summing the signal
            sum_sig = sum_sig + sorted_oval_information{i}{j}(2,:);
        end
        %identifying trials containing empty bin
        if isempty(location_by_split{i})
            ind_tr(end+1) = trial;
            ind_tr1(end+1) = trial;
            break;
            %             mean_oval_signal(trial,data.number_of_images*(i-1)+1:data.number_of_images*(i)) = mean_oval_signal(trial,data.number_of_images*(i-2)+1:data.number_of_images*(i-1));
        else
            %averaing the sorted signal in the bin
            mean_oval_signal(trial,length(sorted_oval_information{i}{j}(2,:))*(i-1)+1:length(sorted_oval_information{i}{j}(2,:))*(i)) = sum_sig/size(location_by_split{i},2);
        end
    end
    %     else
    % %         identifying trials outside preferred ratio range
    %         ind_tr(end+1) = trial;
    %     end
    
end

%keeping a copy fr PK on all trials without empty bins of saccade
choice1 = choice;
mean_oval_signal1 = mean_oval_signal;
choice1(ind_tr1) = [];
mean_oval_signal1(ind_tr1,:) = [];
tr_ratio1 = data.true_ratio;
tr_ratio1(ind_tr1) = [];

% eliminating the trials with emtpy bin of saccade or outside preferred range of ratio
tr_ratio = data.true_ratio;
tr_ratio(ind_tr) = [];
choice(ind_tr) = [];
mean_oval_signal(ind_tr,:) = [];
data.current_trial = size(choice,2);
num_trial_used = data.current_trial;



%Analysis of PK after cross validation on all trials with non empty bins
[best_hprs, ~] = CustomRegression.xValidatePK_with_lapse_duration(mean_oval_signal1, choice1, split_parameter, num_image, hpr1, 0, hpr2, 0, 10,learn_num_PK);
% [best_hprs, ~] = xValidatePK(mean_oval_signal, choice, split_parameter, num_image, hpr1, 0, hpr2, 0, 10,learn_num_PK);
disp('Best hyperparameters found!');
trials = size(mean_oval_signal1, 1);
for j = 1:boot_n
    [signal_result, choice_result] = bootstrap(mean_oval_signal1, choice1, trials);
     [sobl(j,:), ~] = CustomRegression.LinearPK_with_lapse_duration(signal_result, choice_result, 0);
    sobl(j,:) = zeros(1,4);
    if j==1 || mod(j,100)==0
        disp(['Bootstraping step ' num2str(j) ' ...']);
    end
    [params_boot(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapse_duration(signal_result, choice_result, split_parameter, num_image, learn_num_PK, best_hprs(1), 0, best_hprs(3), 0);%,hprs, 0, hprs, 1);
    abbl(j,:) = zeros(1,4);
    [abbl(j,:), ~, ~] = CustomRegression.ExponentialPK_with_lapse_duration(signal_result, choice_result, 0);
end
disp('Computing the bias from data...');
%Obtaining bias in saccade
[acc_evi_all,pro_choice_all,acc1] = compute_extCB_all_ovals(data,mean_oval_signal,split_parameter,bin,boot_n,params_boot,choose_abs_PK);
[acc_evi_ideal1,pro_choice_ideal1] = compute_extCB_ideal_closest(data,mean_oval_signal,split_parameter,bin,boot_n,params_boot,choose_abs_PK);
[acc_evi_ideal2,pro_choice_ideal2,acc_evi_ideal_random2,pro_choice_ideal_random2,acc2] = compute_extCB_ideal_two_closest(data,mean_oval_signal,split_parameter,bin,boot_n,params_boot,choose_abs_PK);
[acc_evi_ideal1avg,pro_choice_ideal1avg] = compute_extCB_ideal_avg_closest(data,mean_oval_signal,split_parameter,bin,boot_n,params_boot,choose_abs_PK);



    function [signals, choices] = bootstrap(signals_raw, choices_raw, trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signals = [];
        choices = [];
        for z = 1:trials
            trial_num = sample_nums(z);
            signals = [signals; signals_raw(trial_num, :)];
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