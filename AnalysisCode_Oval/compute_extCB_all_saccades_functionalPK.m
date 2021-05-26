function [acc_evidence_cases,choice_cases,weight_cases,sp_weight_cases,random_landing_cases,saccades_used] = ...
    compute_extCB_all_saccades_functionalPK(data,fixations_per_trial,...
    oval_signal,saccade_dist_edges,oval_coordinates,time_stamps,...
    params,saccade_dist,peripheryPKbound)
% initialize
eps = 0.01;
scaling = prctile(params(:,1).^2,50);
dist_param = prctile(exp(params(:,2)),50);
mu = prctile(-1*(params(:,3).^2),50);
sigma = prctile(params(:,4).^2 + eps,50);
temporal_fixation = prctile(params(:,5),50);% * eps,50);
bias = prctile(params(:,6),50);
% lapse = 1e-4+(1-1e-4)*sigmoid(params(7));%
% lapse = prctile(params(:,7).^2,50);
temporal_duration = prctile(params(:,8),50);% * eps,50);

num_frames = size(oval_signal,3);
num_im = size(oval_signal,2);
% setting accumulated evidence matrix
shown_pos = [1 2 3 4 5 6 7 8 9 10 11 13 14];
% max_fix = max(fixations_per_trial);

for k=1:length(saccade_dist_edges)-1
    random_landing_cases{k} = [];
    acc_evidence_cases{k} = [];
    choice_cases{k} = [];
    weight_cases{k} = [];
    sp_weight_cases{k} = [];
end
for tr=1:size(oval_signal,1)
    sg = squeeze(oval_signal(tr,:,:));
    sg = sg(:);
    for fx=1:num_frames
        temporal_weights(tr,((fx-1) * num_im + 1):fx * num_im) = (ones(1,num_im) * exp(temporal_fixation * fx)) * exp(temporal_duration * time_stamps(tr,fx));
    end
    elliptical_pixel_dist = squeeze(sqrt(oval_coordinates(tr,1,:).^2 + dist_param * oval_coordinates(tr,2,:).^2));
    spatial_weights(tr,:) = exp((-1*(elliptical_pixel_dist-mu).^2)./(sigma^2)) * scaling;
    logits(tr) = dot(sg,spatial_weights(tr,:) .* temporal_weights(tr,:)) + bias;
end
for trial=1:size(oval_signal,1)
    for i=1:fixations_per_trial(trial)-1
        sg = squeeze(oval_signal(trial,:,1:i));
        sg = sg(:);
        acc_evidence{trial}(i) = dot(sg,spatial_weights(trial,1:i*peripheryPKbound) .* temporal_weights(trial,1:i*peripheryPKbound)) + bias;
        choice_in_fav{trial}(i) = sign_comparison(acc_evidence{trial}(i),squeeze(oval_signal(trial,1,i+1))); %closest element(1) for the next saccade(i+1)
        random_landing{trial}(i) = mean(data.frame_categories(trial,shown_pos)==sign(acc_evidence{trial}(i)));
        weight{trial}(i) = squeeze(spatial_weights(trial,i*peripheryPKbound+1) * temporal_weights(trial,i*peripheryPKbound+1));
        sp_weight{trial}(i) = squeeze(spatial_weights(trial,i*peripheryPKbound+1));
        for k=1:length(saccade_dist_edges)-1
            if (saccade_dist{trial}(i)>=saccade_dist_edges(k) && saccade_dist{trial}(i)<saccade_dist_edges(k+1))
                acc_evidence_cases{k}(end+1) = acc_evidence{trial}(i);
                choice_cases{k}(end+1) = choice_in_fav{trial}(i);
                weight_cases{k}(end+1) = weight{trial}(i);
                sp_weight_cases{k}(end+1) = sp_weight{trial}(i);
                random_landing_cases{k}(end+1) = random_landing{trial}(i);
            end
        end
    end
end
for k=1:length(saccade_dist_edges)-1
    saccades_used(k) = length(acc_evidence_cases{k});
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

end



