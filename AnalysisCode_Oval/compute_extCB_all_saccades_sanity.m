function [mean_bin_index,prob,mean_bin_index_random,prob_random,acc_evidence,choice_in_fav,random_landing] = compute_extCB_all_saccades_sanity(data,fixations_per_trial,oval_signal,bin,boot_n,pk, trials_chosen,peripheryPKbound)
% initialize

temporal_kernel = prctile(pk(:,1:end-2), 50);
bias = prctile(pk(:,end-1), 50);
num_image = 13;%data.number_of_images(1);
% setting accumulated evidence matrix
shown_pos = [1 2 3 4 5 6 7 8 9 10 11 13 14];
max_fix = max(fixations_per_trial);
acc_evidence_all = [];
choice_all = [];
random_landing_all = [];
for trial=1:size(oval_signal,1)
    if trials_chosen(trial)==1
        for i=1:fixations_per_trial(trial)-1
            sg = squeeze(oval_signal(trial,:,1:i));
            sg = sg(:);
            acc_evidence{trial}(i) = bias + sum(sg .* temporal_kernel(1:i*peripheryPKbound)');
            acc_evidence_all(end+1) = acc_evidence{trial}(i);
            choice_in_fav{trial}(i) = sign_comparison(acc_evidence{trial}(i),oval_signal(trial,i+1,1));
            %         choice_in_fav{trial}(i) = sign_comparison(acc_evidence{trial}(i),mean(oval_signal(trial,i+1,1),oval_signal(trial,i+1,2)));
            choice_all(end+1) = choice_in_fav{trial}(i);
            random_landing{trial}(i) = mean(data.frame_categories(trial,shown_pos)==sign(acc_evidence{trial}(i)));
            random_landing_all(end+1) = random_landing{trial}(i);
        end
    end
end

%reshaping into a column vetor
choice_in_fav = reshape(choice_all,[],1);
acc_evidence = reshape(acc_evidence_all,[],1);
%absolute value of acc_evidence
acc_evidence = abs(acc_evidence);
trials = size(acc_evidence, 1);

for j = 1:boot_n
    
    [signal_result, choice_result] = bootstrap(acc_evidence_all,choice_all, length(acc_evidence_all));
    [signal_random, choice_random] = bootstrap(acc_evidence_all,random_landing_all, length(acc_evidence_all));
    
    %bootstrap from here
    %sorting the signal, correspondingly sorting choice in fav
    [sorted_acc_evidence,order] = sort(signal_result);
    b_choice = choice_result(order);
    %binning the signal
    max_value = max(sorted_acc_evidence);
    min_value = min(sorted_acc_evidence);
    bin_index = linspace(min_value,max_value,bin+1);
    allocation = discretize(sorted_acc_evidence,bin_index);
    mean_bin_index(j,:) = linspace((bin_index(1)+bin_index(2))/2,(bin_index(length(bin_index)-1)+bin_index(length(bin_index)))/2,bin);
    %computing choice in favor
    for i=1:size(mean_bin_index,2)
        prob(j,i) = sum(b_choice(allocation==i))/length(b_choice(allocation==i));
    end
    
    
    [sorted_acc_evidence_random,order] = sort(signal_random);
    b_choice_random = choice_random(order);
    %binning the signal
    max_value_random = max(sorted_acc_evidence_random);
    min_value_random = min(sorted_acc_evidence_random);
    bin_index_random = linspace(min_value_random,max_value_random,bin+1);
    allocation_random = discretize(sorted_acc_evidence_random,bin_index_random);
    mean_bin_index_random(j,:) = linspace((bin_index_random(1)+bin_index_random(2))/2,(bin_index_random(length(bin_index_random)-1)+bin_index_random(length(bin_index_random)))/2,bin);
    %computing choice in favor
    for i=1:size(mean_bin_index,2)
        prob_random(j,i) = sum(b_choice_random(allocation==i))/length(b_choice_random(allocation==i));
    end
    
end
end


function [result]=sign_comparison(a,b)
if b>0
    if a+b>b
        result=1;
    else
        result=0;
    end
elseif b<0
    if a+b<b
        result=1;
    else
        result=0;
    end
else
    if a==0
        result=1;
    else
        result=0;
    end
end

end

function [signals, choices] = bootstrap(signals_raw, choices_raw, trials)
sample_nums = randsample(trials, trials, true); % random sample with replacement
signals = [];
choices = [];
for z = 1:trials
    trial_num = sample_nums(z);
    signals = [signals; signals_raw(trial_num)];
    choices = [choices choices_raw(trial_num)];
end
end



