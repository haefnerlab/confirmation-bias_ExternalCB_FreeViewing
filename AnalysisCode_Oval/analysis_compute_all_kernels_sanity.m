function [params_boot,sobl,abbl,best_hprs,data,accuracy,num_image,mean_oval_signal,acc_evi_all,pro_choice_all,...
    acc_evi_ideal1,pro_choice_ideal1,acc_evi_ideal2,pro_choice_ideal2,acc_evi_ideal1avg,pro_choice_ideal1avg,num_trial]...
    = analysis_compute_all_kernels_sanity(data, choice, num_image,num_trial,boot_n, hpr1, hpr2, split_parameter, bin, choose_abs_PK, learn_num_PK)


% choice = data.choice(1:data.current_trial);
% accuracy = data.accuracy(1:data.current_trial);
% num_image = data.number_of_images(1);
mean_oval_signal = data;
accuracy = 0;
% num_trial= data.current_trial;
[best_hprs, ~] = xValidatePK_with_lapse(mean_oval_signal, choice, split_parameter, num_image, hpr1, 0, hpr2, 0, 10,learn_num_PK);
trials = size(mean_oval_signal, 1);


for j = 1:boot_n
    [signal_result, choice_result] = bootstrap(mean_oval_signal, choice, trials);
%     [sobl(j,:), ~] = LinearPK_with_lapse(signal_result, choice_result, 0);
    sobl(j,:) = zeros(1,4);
    disp(j)
    [params_boot(j,:), ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(signal_result, choice_result, split_parameter, num_image, learn_num_PK, best_hprs(1), 0, best_hprs(3), 0);%,hprs, 0, hprs, 1);
    abbl(j,:) = zeros(1,4);
    %     [abbl(j,:), ~, ~] = ExponentialPK_with_lapse(signal_result, choice_result, 0);
end

[acc_evi_all,pro_choice_all] = compute_extCB_all_ovals(data,mean_oval_signal,split_parameter,bin,boot_n,params_boot,choose_abs_PK);
[acc_evi_ideal1,pro_choice_ideal1] = compute_extCB_ideal_closest(data,mean_oval_signal,split_parameter,bin,boot_n,params_boot,choose_abs_PK);
[acc_evi_ideal2,pro_choice_ideal2] = compute_extCB_ideal_two_closest(data,mean_oval_signal,split_parameter,bin,boot_n,params_boot,choose_abs_PK);
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