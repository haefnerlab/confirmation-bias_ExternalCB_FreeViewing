function [mean_bin_index,prob,mean_bin_index_random,prob_random,acc_evidence] = compute_extCB_ideal_two_closest(data,mean_oval_signal,split_parameter,bin,boot_n,pk,choose_abs_PK)
% initialize
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');
temporal_kernel = prctile(pk(:,1:end-2), 50);
bias = prctile(pk(:,end-1), 50);
num_image = 13;%data.number_of_images(1);
shown_pos = [1 2 3 4 5 6 7 8 9 10 11 13 14];
% setting accumulated evidence matrix
acc_evidence = zeros(size(mean_oval_signal,1),split_parameter-1);
for trial=1:size(mean_oval_signal,1)
    temp = 0;
    for i=1:split_parameter-1
        weight = temporal_kernel(((i-1)*num_image+1):((i-1)*num_image+2));
        if choose_abs_PK==1
            weight(weight<0) = 0;
        end
        if i==1
            temp = temp + bias + sum(mean_oval_signal(trial,((i-1)*num_image+1):((i-1)*num_image+2)).*weight);
        else
            temp = temp + sum(mean_oval_signal(trial,((i-1)*num_image+1):((i-1)*num_image+2)).*weight);
        end
        acc_evidence(trial,i) = temp;
    end
end
%setting choice in favor logical matrix
for trial=1:size(acc_evidence,1)
    for i=1:split_parameter-1
        random_landing(trial,i) = mean(data.frame_categories(trial,shown_pos)==sign(acc_evidence(trial,i)));
        choice_in_fav(trial,i) = sign_comparison(mean_oval_signal(trial,i*num_image+1),acc_evidence(trial,i));
    end
end

%reshaping into a column vetor
choice_in_fav = reshape(choice_in_fav,[],1);
acc_evidence = reshape(acc_evidence,[],1);
%absolute value of acc_evidence
acc_evidence=abs(acc_evidence);


trials = size(acc_evidence, 1);


for j = 1:boot_n
    
    [signal_result, choice_result] = bootstrap(acc_evidence,choice_in_fav, trials);
    [signal_random, choice_random] = bootstrap(acc_evidence,random_landing, trials);
    
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

    function [result]=sign_comparison(a,b)
        if b>0
            if a+b>b
                result=1;
            else
                result=0;
            end
        else
            if a+b<b
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
            signals = [signals; signals_raw(trial_num, :)];
            choices = [choices choices_raw(trial_num)];
        end
    end

end


