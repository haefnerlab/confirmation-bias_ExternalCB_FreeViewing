function [mean_bin_index,prob] = compute_extCB_ideal_avg_closest(data,mean_oval_signal,split_parameter,bin,boot_n,pk,choose_abs_PK)
% initialize
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');
temporal_kernel = prctile(pk(:,1:end-2), 50);
bias = prctile(pk(:,end-1), 50);
num_image = 13;%data.number_of_images(1);
acc_evidence = zeros(size(mean_oval_signal,1),split_parameter-1);
% setting accumulated evidence matrix
for trial=1:size(mean_oval_signal,1)
    temp = 0;
    for i=1:split_parameter-1
        weight = mean(temporal_kernel(([1:split_parameter]-1)*num_image+1));
        if choose_abs_PK==1
            weight(weight<0) = 0;
        end
        if i==1
            temp = temp + bias + mean_oval_signal(trial,(i-1)*num_image+1)*weight;
        else
            temp = temp + mean_oval_signal(trial,(i-1)*num_image+1)*weight;
        end
        acc_evidence(trial,i) = temp;
    end
end
%setting choice in favor logical matrix
for trial=1:size(acc_evidence,1)
    for i=1:split_parameter-1
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
            signals = [signals; signals_raw(trial_num, :)];
            choices = [choices choices_raw(trial_num)];
        end
    end

end


