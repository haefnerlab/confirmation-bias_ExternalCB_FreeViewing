close all; clear all; clc;
warning('off');
fix_cond = 'saccade_code';%'eyelink';
boots = 50;
bin_pk = [1081 1921];
fix_cluster_dist = 240; % value used to cluster fixations closeby as one point!! prevents bias from artifacts!!
boots_thresh = 1; %bootstraps to get threshold, takes too much time so set to 1
threshold_performance = [0.5 0.75];
hprs_ridge = logspace(-3,3,7);
folds = 10;

subjects = ...
    {'OvalEyeMovement-subject04'; 'OvalEyeMovement-subject03';...
    'OvalEyeMovement-subject07'; 'OvalEyeMovement-subject06';...
    'OvalEyeMovement-subject08'; 'OvalEyeMovement-subject09';...
    'OvalEyeMovement-subject10'; 'OvalEyeMovement-subject11';...
    'OvalEyeMovement-subject12'; 'OvalEyeMovement-subject13'};
% subjects = {'OvalEyeMovement-subject12'; 'OvalEyeMovement-subject13'};

[num_sub,~] = size(subjects);
sanity_bnd = [0.5 0.55];
saccade_types = [1 2 3 5 7]; % number of saccade bins
evidence_types = [2 3 5 7];
saccade_types_sanity = 1;
saccade_types_thresh = 1;
expt_type = 1; % ratio change
peripheryPKbound = 13; % number of peripheral images to compute PK on based on distance from fixation
format short g;
disp_fig = 0;
screen_mid = [1920/2 1080/2];
standardize = 0;

disp('Starting search of best hyper parameters across all subjects...');
[best_hpr_ridge]= CustomRegression.combined_cross_validated_hprs(subjects, expt_type,...
    fix_cond,fix_cluster_dist,peripheryPKbound,hprs_ridge,standardize,folds);
disp(['Best ridge parameters across subjects is: ' num2str(best_hpr_ridge)]);
for sub=1:num_sub
    tic;
    disp(['Starting analysis for subject ' num2str(sub) ' ...']);
    [params, data, accuracy, num_image, subj_performance_bnd, mean_oval_signal,...
        num_trial, num_trial_sanity, num_trial_thresh,...
        acc_evidence, choice_in_favor, weights_on_saccades, sp_weights_on_saccades, random_landing,...
        acc_evidence_sanity, choice_in_favor_sanity, random_landing_sanity,...
        acc_evidence_thresh, choice_in_favor_thresh, weights_on_saccades_thresh, sp_weights_on_saccades_thresh, random_landing_thresh,...
        fixations_num, eyelink, saccade_code, fixation_dist,...
        saccades_used, saccades_used_sanity, saccades_used_thresh, saccade_dist_mid,...
        logits, spatial_pk, temporal_fixation_pk, temporal_duration_pk, temporal_pk_across_trials]...
        = analysis_across_allsaccades_functionalPK_SingleHyprs(subjects{sub}, expt_type, boots, fix_cond,...
        saccade_types, saccade_types_sanity, saccade_types_thresh, fix_cluster_dist, peripheryPKbound,...
        sanity_bnd, boots_thresh, threshold_performance, best_hpr_ridge, standardize);
    
    each_trial_duration = data.stim_duration;
    floor(sub) = subj_performance_bnd(1);
    thresh(sub) = subj_performance_bnd(2);
    
    eps = 0.1;
    scaling(sub,1) = prctile(params(:,1).^2,50);
    scaling(sub,2) = prctile(params(:,1).^2,50)-prctile(params(:,1).^2,16);
    scaling(sub,3) = prctile(params(:,1).^2,84)-prctile(params(:,1).^2,50);
    dist_param(sub,1) = prctile(exp(params(:,2)),50);
    dist_param(sub,2) = prctile(exp(params(:,2)),50)-prctile(exp(params(:,2)),16);
    dist_param(sub,3) = prctile(exp(params(:,2)),84)-prctile(exp(params(:,2)),50);
    mu(sub,1) = prctile(-1*(params(:,3).^2),50);
    mu(sub,2) = prctile(-1*(params(:,3).^2),50)-prctile(-1*(params(:,3).^2),16);
    mu(sub,3) = prctile(-1*(params(:,3).^2),84)-prctile(-1*(params(:,3).^2),50);
    sigma(sub,1) = prctile(params(:,4).^2 + eps,50);
    sigma(sub,2) = prctile(params(:,4).^2 +  eps,50)-prctile(params(:,4).^2 + eps,16);
    sigma(sub,3) = prctile(params(:,4).^2 + eps,84)-prctile(params(:,4).^2 + eps,50);
    temporal_fixation(sub,1) = prctile(params(:,5),50);
    temporal_fixation(sub,2) = prctile(params(:,5),50)-prctile(params(:,5),16);
    temporal_fixation(sub,3) = prctile(params(:,5),84)-prctile(params(:,5),50);
    bias(sub,1) = prctile(params(:,6),50);
    bias(sub,2) = prctile(params(:,6),50)-prctile(params(:,6),16);
    bias(sub,3) = prctile(params(:,6),84)-prctile(params(:,6),50);
    lapse(sub,1) = prctile(params(:,7).^2,50);
    lapse(sub,2) = prctile(params(:,7).^2,50)-prctile(params(:,7).^2,16);
    lapse(sub,3) = prctile(params(:,7).^2,84)-prctile(params(:,7).^2,50);
    %     lapse(sub,1) = prctile((1e-4+(1-1e-4)*sigmoid(params(:,7))),50);
    %     lapse(sub,2) = prctile((1e-4+(1-1e-4)*sigmoid(params(:,7))),50)-prctile((1e-4+(1-1e-4)*sigmoid(params(:,7))),16);
    %     lapse(sub,3) = prctile((1e-4+(1-1e-4)*sigmoid(params(:,7))),84)-prctile((1e-4+(1-1e-4)*sigmoid(params(:,7))),50);
    temporal_duration(sub,1) = prctile(params(:,8),50);
    temporal_duration(sub,2) = prctile(params(:,8),50)-prctile(params(:,8),16);
    temporal_duration(sub,3) = prctile(params(:,8),84)-prctile(params(:,8),50);
    
    trial_num(sub) = num_trial;
    trial_num_sanity(sub) = num_trial_sanity;
    trial_num_thresh(sub) = num_trial_thresh;
    fixations{sub} = fixations_num;
    max_fixations(sub) = max(fixations_num);
    eyelink_per_subj{sub} = eyelink;
    saccade_code_per_subj{sub} = saccade_code;
    fixation_dist_per_subj{sub} = fixation_dist;
    oval_signal{sub} = mean_oval_signal;
    
    spatial_kernel{sub} = spatial_pk;
    temporal_fixation_kernel{sub} = temporal_fixation_pk;
    temporal_duration_kernel{sub} = temporal_duration_pk;
    temporal_kernel_across_trials{sub} = temporal_pk_across_trials;
    log_bernoulli{sub} = logits;
    
    saccades_per_subj_sanity(sub) = saccades_used_sanity;
    saccades_per_subj_thresh(sub) = saccades_used_thresh;
    
    for bn=1:length(saccade_types)
        saccades_per_subj{bn}(sub,:) = saccades_used{bn};
        saccade_dist_per_subj{bn}(sub,:) = saccade_dist_mid{bn};
        for ss=1:saccade_types(bn)
            acc_evi_mid{bn}(sub,ss) = mean(acc_evidence{bn}{ss});
            acc_evi_mid_err{bn}(sub,ss) = sqrt(var(acc_evidence{bn}{ss})/length(acc_evidence{bn}{ss}));
            
            prob_choice_favor{bn}(sub,ss) = mean(choice_in_favor{bn}{ss});
            prob_choice_favor_err{bn}(sub,ss) = sqrt((prob_choice_favor{bn}(sub,ss) * (1 - prob_choice_favor{bn}(sub,ss)))/length(choice_in_favor{bn}{ss}));
            
            prob_choice_favor_random{bn}(sub,ss) = mean(random_landing{bn}{ss});
            prob_choice_favor_random_err{bn}(sub,ss) = sqrt(var(random_landing{bn}{ss})/length(random_landing{bn}{ss}));
        end
    end
    acc_evi_sanity(sub,1) = mean(acc_evidence_sanity{1});
    acc_evi_sanity(sub,2) = sqrt(var(acc_evidence_sanity{1})/length(acc_evidence_sanity{1}));
    
    prob_choice_favor_sanity(sub,1) = mean(choice_in_favor_sanity{1});
    prob_choice_favor_sanity(sub,2) = sqrt((prob_choice_favor_sanity(sub,1) * (1 - prob_choice_favor_sanity(sub,1)))/length(choice_in_favor_sanity{1}));
    
    prob_choice_favor_sanity_random(sub,1) = mean(random_landing_sanity{1});
    prob_choice_favor_sanity_random(sub,2) = sqrt(var(random_landing_sanity{1})/length(random_landing_sanity{1}));
    
    acc_evi_thresh(sub,1) = mean(acc_evidence_thresh{1});
    acc_evi_thresh(sub,2) = sqrt(var(acc_evidence_thresh{1})/length(acc_evidence_thresh{1}));
    
    prob_choice_favor_thresh(sub,1) = mean(choice_in_favor_thresh{1});
    prob_choice_favor_thresh(sub,2) = sqrt((prob_choice_favor_thresh(sub,1) * (1 - prob_choice_favor_thresh(sub,1)))/length(choice_in_favor_sanity{1}));
    
    prob_choice_favor_thresh_random(sub,1) = mean(random_landing_thresh{1});
    prob_choice_favor_thresh_random(sub,2) = sqrt(var(random_landing_thresh{1})/length(random_landing_thresh{1}));
    
    thresh_saccade_weights{sub} = weights_on_saccades_thresh{1};
    thresh_choice_saccades{sub} = choice_in_favor_thresh{1};
    thresh_random_landing_weights{sub} = random_landing_thresh{1};
    saccade_weights{sub} = weights_on_saccades{1}{1};
    choice_saccades{sub} = choice_in_favor{1}{1};
    random_landing_weights{sub} = random_landing{1}{1};
    
    thresh_saccade_sp_weights{sub} = sp_weights_on_saccades_thresh{1};
    saccade_sp_weights{sub} = sp_weights_on_saccades{1}{1};
    
    tr_ratio{sub} = data.true_ratio;
    ratio_used = data.ratio;
    choice{sub} = data.choice;
    
    disp(['Psychometric analysis for subject ' num2str(sub) ' ...']);
    ratios = linspace(0, 1, 101);%used in plots later
    [pm_fit, uniq_vals, yvals, stderrs] = GaborPsychometric(data, 1);
    subject_pm_curve_psig(sub,:) = (1-pm_fit.Fit(3)-pm_fit.Fit(4))*arrayfun(@(x) pm_fit.options.sigmoidHandle(x,pm_fit.Fit(1),pm_fit.Fit(2)), ratios)+pm_fit.Fit(4);
    
    disp('Analysis of bias w.r.t magnitude of evidence accumulated ...');
    [all_evidence_acc, evi_order] = sort(abs(acc_evidence{1}{1}));
    all_choice_subj = choice_in_favor{1}{1}(evi_order);
    all_prob_baseline = random_landing{1}{1}(evi_order);
    for ev_bn=1:length(evidence_types)
        evi_bin_edges = linspace(min(all_evidence_acc), max(all_evidence_acc),evidence_types(ev_bn)+1);
        for bb=1:evidence_types(ev_bn)
            binned_acc_evidence{ev_bn}(sub,bb) = mean(all_evidence_acc(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1)));
            binned_prob_bias{ev_bn}(sub,bb) = mean(all_choice_subj(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1)));
            binned_prob_bias_err{ev_bn}(sub,bb) = sqrt((binned_prob_bias{ev_bn}(sub,bb)*(1-binned_prob_bias{ev_bn}(sub,bb)))/length(all_choice_subj(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1))));
            binned_baseline{ev_bn}(sub,bb) = mean(all_prob_baseline(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1)));
            binned_baseline_err{ev_bn}(sub,bb) = sqrt(var(all_prob_baseline(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1)))/length(all_prob_baseline(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1))));
        end
    end
    
    disp(['All Analysis for subject ' num2str(sub) ' complete!']);
    disp('-----------------------------------------------------------------------------------------------------');
    toc;
    
    if disp_fig==1
        figure();
        subplot(2,3,1)
        pk = spatial_kernel{sub};
        imagesc(pk);set(gca,'YDir','normal'); h = colorbar; xlabel(h,'Spatial Weights','fontsize',8);h.LineWidth=1.5; h.FontSize=15;
        hold on;
        scatter(screen_mid(1), screen_mid(2), 50, '+', 'r' );
        ylabel('Screen y-axis','fontsize',20);
        xlabel('Screen x-axis','fontsize',20);
        xticks([1 fix(1920/2) 1920]);
        yticks([1 fix(1080/2) 1080]);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(2,3,2)
        temporal_fixation_pk = temporal_fixation_kernel{sub};
        plot(1:length(temporal_fixation_pk),temporal_fixation_pk,'-ob','LineWidth',2)
        hold on;
        ylabel('Fixation weights','fontsize',20);
        xlabel('Fixation Number','fontsize',20);
        xticks([1 fix(max_fixations(sub)/2) max_fixations(sub)]);
        xline(0.0,'k','LineWidth',0.75);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(2,3,3)
        temporal_duration_pk = temporal_duration_kernel{sub};
        plot(linspace(0,each_trial_duration,length(temporal_duration_pk)),temporal_duration_pk,'-ob','LineWidth',2)
        hold on;
        ylabel('Temporal weights','fontsize',20);
        xlabel('Time stamp of fixation','fontsize',20);
        xticks([0 each_trial_duration/2 each_trial_duration]);
        xline(0.0,'k','LineWidth',0.75);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(2,3,5)
        llo_mn = [];
        choice_mn = [];
        err_ch = [];
        err_ch_mn = [];
        llo_mean = [];
        choice_mean = [];
        bin_num = 50;
        [sorted_llo,order_llo] = sort(log_bernoulli{sub});
        choice_used = choice{sub}(order_llo);
        bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
        for bn=1:length(bin_edges)-1
            llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
            choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
            err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
        end
        [llo_mn, ord] = sort(llo_mean);
        choice_mn = choice_mean(ord);
        err_ch_mn = err_ch(ord);
        errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
        %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
        ylabel('Probability chose vertical','Fontsize',20);
        xlabel('Log likelihood odds','Fontsize',20);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([0.0 1]);
        xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
        yline(0.5,'-k','linewidth',2);
        hold on;
        xline(0.0,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(2,3,4)
        ratios = linspace(0, 1, 101);
        uniq_vals=linspace(-0.05,1.05,12);
        for tt=1:(length(uniq_vals)-1)
            subj_resp(sub,tt)=mean(choice{sub}(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)));
            ntrial_subj(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
        end
        uniq_vals=linspace(0,1,11);
        errorbar(uniq_vals,subj_resp(sub,:),sqrt((subj_resp(sub,:)).*(1-subj_resp(sub,:))./ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
        uniq_vals=linspace(-0.05,1.05,12);
        for tt=1:(length(uniq_vals)-1)
            pred_resp = lapse(1)*0.5+(1-lapse(1))*sigmoid(log_bernoulli{sub}(((tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)))));
            subj_resp_pred(sub,tt)=mean(pred_resp);
            ntrial_subj_pred(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
            subj_resp_pred_err(sub,tt) = sqrt(subj_resp_pred(sub,tt)*(1-subj_resp_pred(sub,tt))/ntrial_subj_pred(sub,tt));
        end
        uniq_vals=linspace(0,1,11);
        hold on;
        errorbar(uniq_vals,subj_resp_pred(sub,:),subj_resp_pred_err(sub,:),'--or','Linewidth',2);
        ylabel('Probability chose vertical','Fontsize',20);
        xlabel({'Ratio of vertical ellipses'},'Fontsize',20);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([0.0 1]);
        yline(0.5,'-k','linewidth',2);
        hold on;
        xline(0.5,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(2,3,6)
        bar(1, prob_choice_favor{1}(sub),'FaceColor','b','EdgeColor','k','LineWidth',0.75);
        hold on;
        bar(2, prob_choice_favor_random{1}(sub),'FaceColor','r','EdgeColor','k','LineWidth',0.75);
        hold on;
        errorbar([1 2],[prob_choice_favor{1}(sub) prob_choice_favor_random{1}(sub)],[squeeze(prob_choice_favor_err{1}(sub)) squeeze(prob_choice_favor_random_err{1}(sub))],'ok','LineWidth',2,'linestyle','none');
        xlim([0.25 2.75]);
        text(1,0.9,['Trial Num: ' num2str(trial_num(sub))],'Fontsize',12);
        hold on;
        text(1,0.85,['Sacc Num: ' num2str(sum(fixations{sub}))],'Fontsize',12);
        xticks([1 2]);
        xticklabels({'Subject'; 'Baseline'});
        ylabel('Probability of chose in favor','Fontsize',20);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([0.4 1]);
        yline(0.5,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        sgtitle(['Subject ' num2str(sub)],'fontsize',30);
        
        %         pause;
        close all;
    end
end
%%
disp('Plotting stuff may take a while ...............................');
for sub=1:num_sub
    fig = figure(sub);
    set(fig,'defaultLegendAutoUpdate','off');
    
    subplot(2,4,5)
    %     [floor(sub), thresh(sub), ~] = getThresholdWindow(data,1, 0.5, 0.7);
    plot(ratio_used,'-ob');
    xlabel('Trial number','Fontsize',20);
    ylabel({'Ratio of vertical ellipses'},'Fontsize',20);
    yline(0.5,'k','LineWidth',1.5);
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.5 1.0]);
    xlim([1 length(ratio_used)]);
    %     yline(thresh(sub),'-k','linewidth',2)
    %     hold on;
    %     yline(floor(sub),'-k','linewidth',2)
    %     yline(0.5,'-k','linewidth',2)
    %     ylim([0 1]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Staircase','fontsize',20);
    
    subplot(2,4,1)
    pk = spatial_kernel{sub};
    imagesc(pk);set(gca,'YDir','normal'); h = colorbar; xlabel(h,'Spatial Weights','fontsize',8);h.LineWidth=1.5; h.FontSize=15;
    hold on;
    scatter(screen_mid(1), screen_mid(2), 50, '+', 'r' );
    ylabel('Screen y-axis','fontsize',20);
    xlabel('Screen x-axis','fontsize',20);
    xticks([1 fix(1920/2) 1920]);
    yticks([1 fix(1080/2) 1080]);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Spatial Weights','fontsize',20);
    
    
    subplot(2,4,2)
    temporal_fixation_pk = temporal_fixation_kernel{sub};
    plot(1:length(temporal_fixation_pk),temporal_fixation_pk,'-ob','LineWidth',2)
    hold on;
    ylabel('Fixation weights','fontsize',20);
    xlabel('Fixation Number','fontsize',20);
    xticks([1 fix(max_fixations(sub)/2) max_fixations(sub)]);
    xlim([1 max_fixations(sub)]);
    xline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Fixation weights','Fontsize',20);
    
    subplot(2,4,3)
    temporal_duration_pk = temporal_duration_kernel{sub};
    plot(linspace(0,each_trial_duration,length(temporal_duration_pk)),temporal_duration_pk,'-ob','LineWidth',2)
    hold on;
    ylabel('Temporal weights','fontsize',20);
    xlabel('Timestamps','fontsize',20);
    xticks([0 each_trial_duration/2 each_trial_duration]);
    xline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Timestamp weights','Fontsize',20)
    
    subplot(2,4,4)
    temporal_pk = temporal_kernel_across_trials{sub};
    plot(1:length(temporal_pk),temporal_pk,'-ob','LineWidth',2)
    hold on;
    ylabel('Combined temporal weights','fontsize',20);
    xlabel('Fixation Number','fontsize',20);
    xticks([1 fix(max_fixations(sub)/2) max_fixations(sub)]);
    xlim([1 max_fixations(sub)]);
    xline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Combined temporal','Fontsize',20);
    
    subplot(2,4,6)
    llo_mn = [];
    choice_mn = [];
    err_ch = [];
    err_ch_mn = [];
    llo_mean = [];
    choice_mean = [];
    bin_num = 20;
    [sorted_llo,order_llo] = sort(log_bernoulli{sub});
    choice_used = choice{sub}(order_llo);
    bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
    for bn=1:length(bin_edges)-1
        llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
    end
    [llo_mn, ord] = sort(llo_mean);
    choice_mn = choice_mean(ord);
    err_ch_mn = err_ch(ord);
    errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
    %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
    ylabel('Probability chose vertical','Fontsize',20);
    xlabel('Log likelihood odds','Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('LLO Sanity','fontsize',20);
    
    subplot(2,4,7)
    ratios = linspace(0, 1, 101);
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp(sub,tt)=mean(choice{sub}(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)));
        ntrial_subj(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    errorbar(uniq_vals,subj_resp(sub,:),sqrt((subj_resp(sub,:)).*(1-subj_resp(sub,:))./ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        pred_resp = lapse(1)*0.5+(1-lapse(1))*sigmoid(log_bernoulli{sub}(((tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)))));
        subj_resp_pred(sub,tt) = mean(pred_resp);
        ntrial_subj_pred(sub,tt) = sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
        subj_resp_pred_err(sub,tt) = sqrt(subj_resp_pred(sub,tt)*(1-subj_resp_pred(sub,tt))/ntrial_subj_pred(sub,tt));
    end
    uniq_vals=linspace(0,1,11);
    hold on;
    errorbar(uniq_vals,subj_resp_pred(sub,:),subj_resp_pred_err(sub,:),'--or','Linewidth',2);
    ylabel('Probability chose vertical','Fontsize',20);
    xlabel({'Ratio of vertical ellipses'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Match to data','fontsize',20);
    
    subplot(2,4,8)
    hold on;
    LH3(1)=bar(1, prob_choice_favor{1}(sub),'FaceColor','b','EdgeColor','k','LineWidth',0.75);
    L3{1} = 'Subject';
    hold on;
    LH3(2)=bar(2, prob_choice_favor_random{1}(sub),'FaceColor','r','EdgeColor','k','LineWidth',0.75);
    L3{2} = 'Baseline';
    hold on;
    errorbar([1 2],[prob_choice_favor{1}(sub) prob_choice_favor_random{1}(sub)],[squeeze(prob_choice_favor_err{1}(sub)) squeeze(prob_choice_favor_random_err{1}(sub))],'ok','LineWidth',2,'linestyle','none');
    xlim([0.25 2.75]);
    text(1,0.9,['Trial Num: ' num2str(trial_num(sub))],'Fontsize',12);
    hold on;
    text(1,0.85,['Sacc Num: ' num2str(sum(fixations{sub}))],'Fontsize',12);
    xticks([1 2]);
    xticklabels({'Subject'; 'Baseline'});
    ylabel('Probability of chose in favor','Fontsize',20);
    hold on;
    legend(LH3,L3, 'Fontsize',15, 'Box','off');
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('External CB','fontsize',20);
    
    sgtitle(['All analysis for Subject ' num2str(sub)],'fontsize',30);
end

%%
figure();
for sub=1:num_sub
    subplot(2,5,sub)
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp(sub,tt)=mean(choice{sub}(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)));
        ntrial_subj(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    LH1(1) = errorbar(uniq_vals,subj_resp(sub,:),sqrt((subj_resp(sub,:)).*(1-subj_resp(sub,:))./ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
    L1{1} = 'Data';
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        pred_resp = lapse(1)*0.5+(1-lapse(1))*sigmoid(log_bernoulli{sub}(((tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)))));
        subj_resp_pred(sub,tt)=mean(pred_resp);
        subj_resp_pred_err(sub,tt) = sqrt(subj_resp_pred(sub,tt)*(1-subj_resp_pred(sub,tt))/ntrial_subj_pred(sub,tt));
        ntrial_subj_pred(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    hold on;
    LH1(2) = errorbar(uniq_vals,subj_resp_pred(sub,:),subj_resp_pred_err(sub,:),'-or','Linewidth',2);
    L1{2} = 'Fit with PK';
    if sub==1 || sub==6
        ylabel('Probability chose vertical','Fontsize',20);
    end
    if sub==3 || sub==8
        xlabel({'Ratio of vertical ellipses'},'Fontsize',20);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    if sub==1
        legend(LH1,L1, 'Fontsize',12, 'Box','off','Location','northwest');
    end
    title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('PS fit w.r.t learned Weights/PK compared to real data','fontsize',30)
%%
bin_num=20;
figure();
llo_mn = [];
choice_mn = [];
err_ch = [];
err_ch_mn = [];
llo_mean = [];
choice_mean = [];
for sub=1:num_sub
    subplot(2,5,sub)
    %     llo = llo( ~any( isnan( llo ) | isinf( llo ), 2 ),: );
    [sorted_llo,order_llo] = sort(log_bernoulli{sub});
    choice_used = choice{sub}(order_llo);
    bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
    for bn=1:length(bin_edges)-1
        llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
    end
    [llo_mn, ord] = sort(llo_mean);
    choice_mn = choice_mean(ord);
    err_ch_mn = err_ch(ord);
    errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
    %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
    if sub==1 || sub==6
        ylabel('Probability chose vertical','Fontsize',20);
    end
    if sub==3 || sub==8
        xlabel('Log likelihood odds','Fontsize',20);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Check how good weights predict behavior','fontsize',30);
%%
figure();
ratios = linspace(0, 1, 101);
for sub=1:num_sub
    plot(ratios, subject_pm_curve_psig(sub,:), 'color','k','LineWidth', 0.5);
    hold on;
end
yline(0.5,'-k','linewidth',2);
hold on;
xline(0.5,'-k','linewidth',2);
hold on;
ylabel('Probability chose vertical','Fontsize',20);
xlabel({'Ratio of vertical ellipses'},'Fontsize',20);
hold on;
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ylim([0.0 1]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
plot(ratios, mean(subject_pm_curve_psig,1), 'k','LineWidth', 2);
title('Psychometric Curves for all subjects','fontsize',30);
%%
figure();
max_pk = 0;
min_pk = 100000;
max_fx_pk = 0;
min_fx_pk = 100000;
max_dur_pk = 0;
min_dur_pk = 100000;
max_comb_pk = 0;
min_comb_pk = 100000;
for sub=1:num_sub
    max_pk = max(max_pk,max(spatial_kernel{sub}(:)));
    min_pk = min(min_pk,min(spatial_kernel{sub}(:)));
    max_fx_pk = max(max_fx_pk,max(temporal_fixation_kernel{sub}(:) * scaling(sub,1)));
    min_fx_pk = min(min_fx_pk,min(temporal_fixation_kernel{sub}(:) * scaling(sub,1))) ;
    max_dur_pk = max(max_dur_pk,max(temporal_duration_kernel{sub}(:)));% * scaling(sub);
    min_dur_pk = min(min_dur_pk,min(temporal_duration_kernel{sub}(:)));% * scaling(sub);
    max_comb_pk = max(max_comb_pk,max(temporal_kernel_across_trials{sub}(:)));
    min_comb_pk = min(min_comb_pk,min(temporal_kernel_across_trials{sub}(:)));
end
for sub=1:num_sub
    subplot(3,4,sub)
    pk = spatial_kernel{sub};%./scaling(sub,1);
    imagesc(pk);set(gca,'YDir','normal'); h = colorbar; xlabel(h,'Spatial Weights','fontsize',8);h.LineWidth=1.5; h.FontSize=15;
    axis image;
    hold on;
    scatter(screen_mid(1), screen_mid(2), 50, '+', 'r' );
    if sub==1 || sub==5 || sub==9
        ylabel('Screen y-axis','fontsize',20);
    end
    if sub==2 || sub==6 || sub==10
        xlabel('Screen x-axis','fontsize',20);
    end
    xticks([1 fix(1920/2) 1920]);
    yticks([1 fix(1080/2) 1080]);
    hold on;
    caxis([min_pk, max_pk]);
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Spatial weights w.r.t fixation in center of screen','fontsize',30);

figure();
for sub=1:num_sub
    %     subplot(3,4,sub)
    temporal_fixation_pk = temporal_fixation_kernel{sub} * scaling(sub,1);
    plot(1:length(temporal_fixation_pk),temporal_fixation_pk,'-ob','LineWidth',2)
    hold on;
    if sub==1 || sub==5 || sub==9
        ylabel('Temporal Weights','fontsize',20);
    end
    if sub==2 || sub==6 || sub==10
        xlabel('Fixation Number','fontsize',20);
    end
    %     if sub==7
    %         ylim([0 150]);
    %     else
    %         ylim([0 5]);
    %     end
    %     xticks([1 fix(max_fixations(sub)/2) max_fixations(sub)]);
    xlim([1 max([max_fixations(1) max_fixations(3:end)])]);
    xticks(1:max([max_fixations(1) max_fixations(3:end)]));
    ylim([0.0 max_fx_pk]);
    %     ylim([min_fx_pk max_fx_pk]);
    xline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    %     title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Temporal fixation based weights','fontsize',30);

figure();
for sub=1:num_sub
    %     subplot(3,4,sub)
    temporal_duration_pk = temporal_duration_kernel{sub};% * scaling(sub);
    plot(linspace(0,each_trial_duration,length(temporal_duration_pk)),temporal_duration_pk,'-ob','LineWidth',2)
    hold on;
    if sub==1 || sub==5 || sub==9
        ylabel('Temporal weights','fontsize',20);
    end
    if sub==2 || sub==6 || sub==10
        xlabel('Time stamp of fixation','fontsize',20);
    end
    xticks([0 each_trial_duration/2 each_trial_duration]);
    xlim([0 each_trial_duration]);
    ylim([0.0 max_dur_pk]);
    %     ylim([min_dur_pk max_dur_pk]);
    xline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    %     title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Time stamp based weights','fontsize',30)

figure();
for sub=1:num_sub
    subplot(3,4,sub)
    temporal_pk = temporal_kernel_across_trials{sub};
    plot(1:length(temporal_pk),temporal_pk,'-ob','LineWidth',2)
    hold on;
    if sub==1 || sub==5 || sub==9
        ylabel('Combined temporal weights','fontsize',20);
    end
    if sub==2 || sub==6 || sub==10
        xlabel('Fixation Number','fontsize',20);
    end
    xticks([1 fix(max_fixations(sub)/2) max_fixations(sub)]);
    xline(0.0,'k','LineWidth',0.75);
    xlim([1 max_fixations(sub)]);
    ylim([min_comb_pk max_comb_pk]);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end
sgtitle('Combined fixation and time stamp based weights','fontsize',30);

%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, prob_choice_favor{1}(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
        L6{1} = 'Subject';
    else
        bar(k1, prob_choice_favor{1}(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, prob_choice_favor_random{1}(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
        L6{2} = 'Baseline';
    else
        bar(k2, prob_choice_favor_random{1}(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    errorbar([k1 k2],[prob_choice_favor{1}(sub) prob_choice_favor_random{1}(sub)],[prob_choice_favor_err{1}(sub) prob_choice_favor_random_err{1}(sub)],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.25,0.9,[num2str(trial_num(sub)) ' trials'],'Fontsize',15);
    text(k1-0.25,0.8,[num2str(saccades_per_subj{1}(sub)) ' sac.'],'Fontsize',15);
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of chose';' in favor of accumulated evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 2.5;
    k2 = k2 + 2.5;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k+0.5);
legend(LH6,L6, 'Fontsize',20, 'Box','off');
sgtitle('All saccades for all subjects','fontsize',30);

%%
%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    bar(k1, prob_choice_favor{1}(sub)-prob_choice_favor_random{1}(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    hold on;
    errorbar([k1],[prob_choice_favor{1}(sub)-prob_choice_favor_random{1}(sub)],[sqrt(prob_choice_favor_err{1}(sub)^2+prob_choice_favor_random_err{1}(sub)^2)],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.5,0.175,[num2str(trial_num(sub)) ' trials'],'Fontsize',15);
    text(k1-0.5,0.15,[num2str(saccades_per_subj{1}(sub)) ' sac.'],'Fontsize',15);
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of chose';' in favor of accumulated evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([-0.05 0.2]);
    yline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 1;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k);
% legend(LH6,L6, 'Fontsize',20, 'Box','off');
sgtitle('All saccades for all subjects','fontsize',30);


%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, prob_choice_favor_sanity(sub,1),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
        L6{1} = 'Subject';
    else
        bar(k1, prob_choice_favor_sanity(sub,1),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, prob_choice_favor_sanity_random(sub,1),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
        L6{2} = 'Baseline';
    else
        bar(k2, prob_choice_favor_sanity_random(sub,1),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    errorbar([k1 k2],[prob_choice_favor_sanity(sub,1) prob_choice_favor_sanity_random(sub,1)],[prob_choice_favor_sanity(sub,2) prob_choice_favor_sanity_random(sub,2)],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.25,0.9,[num2str(trial_num_sanity(sub)) ' trials'],'Fontsize',15)
    text(k1-0.25,0.8,[num2str(saccades_per_subj_sanity(sub)) ' sac.'],'Fontsize',15);
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of chose';' in favor of accumulated evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 2.5;
    k2 = k2 + 2.5;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k+0.5);
legend(LH6,L6, 'Fontsize',20, 'Box','off');
sgtitle(['Sanity check figure with ratio > ' num2str(min(sanity_bnd(1),sanity_bnd(2))) ' and < ' num2str(max(sanity_bnd(1),sanity_bnd(2)))],'fontsize',30);

%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, prob_choice_favor_thresh(sub,1),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
        L6{1} = 'Subject';
    else
        bar(k1, prob_choice_favor_thresh(sub,1),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, prob_choice_favor_thresh_random(sub,1),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
        L6{2} = 'Baseline';
    else
        bar(k2, prob_choice_favor_thresh_random(sub,1),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    errorbar([k1 k2],[prob_choice_favor_thresh(sub,1) prob_choice_favor_thresh_random(sub,1)],[prob_choice_favor_thresh(sub,2) prob_choice_favor_thresh_random(sub,2)],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.25,0.9,[num2str(trial_num_thresh(sub)) ' trials'],'Fontsize',15)
    text(k1-0.25,0.75,[num2str(saccades_per_subj_thresh(sub)) ' sac.'],'Fontsize',15);
    text(k1-0.25,0.85,'Thresh ratio: ','Fontsize',15);
    text(k1-0.25,0.82, [num2str(floor(sub),'%1.2f') '     ' num2str(thresh(sub),'%1.2f')],'Fontsize',15);
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of chose';' in favor of accumulated evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 2.5;
    k2 = k2 + 2.5;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k+0.5);
legend(LH6,L6, 'Fontsize',20, 'Box','off');
sgtitle('Saccade bias of all subjects for threshold trials','fontsize',30);

%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    bar(k1, prob_choice_favor_thresh(sub,1)-prob_choice_favor_thresh_random(sub,1),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    hold on;
    errorbar(k1,prob_choice_favor_thresh(sub,1)-prob_choice_favor_thresh_random(sub,1),[sqrt(prob_choice_favor_thresh(sub,2)^2 + prob_choice_favor_thresh_random(sub,2)^2)],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.5,0.175,[num2str(trial_num(sub)) ' trials'],'Fontsize',15);
    text(k1-0.5,0.15,[num2str(saccades_per_subj{1}(sub)) ' sac.'],'Fontsize',15);
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of chose';' in favor of accumulated evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([-0.1 0.2]);
    yline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 1;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k);
% legend(LH6,L6, 'Fontsize',20, 'Box','off');
sgtitle('All saccades from threshold trials for all subjects','fontsize',30);

%%
for bn=2:length(saccade_types)
    fig=figure();
    set(fig,'defaultLegendAutoUpdate','off');
    for sub=1:num_sub
        subplot(2,5,sub)
        if sub==1
            %             LH9(1) = errorbar(1:saccade_types(bn), prob_choice_favor{bn}(sub,:),prob_choice_favor_err{bn}(sub,:),'-ob','LineWidth',2);
            LH9(1) = errorbar(saccade_dist_per_subj{bn}(sub,:), prob_choice_favor{bn}(sub,:),prob_choice_favor_err{bn}(sub,:),'-ob','LineWidth',2);
            L9{1} = 'Subject';
        else
            %             errorbar(1:saccade_types(bn), prob_choice_favor{bn}(sub,:),prob_choice_favor_err{bn}(sub,:),'-ob','LineWidth',2);
            errorbar(saccade_dist_per_subj{bn}(sub,:), prob_choice_favor{bn}(sub,:),prob_choice_favor_err{bn}(sub,:),'-ob','LineWidth',2);
        end
        hold on;
        if sub==1
            %             LH9(2) = errorbar(1:saccade_types(bn), prob_choice_favor_random{bn}(sub,:),prob_choice_favor_random_err{bn}(sub,:),'-or','LineWidth',2);
            LH9(2) = errorbar(saccade_dist_per_subj{bn}(sub,:), prob_choice_favor_random{bn}(sub,:),prob_choice_favor_random_err{bn}(sub,:),'-or','LineWidth',2);
            L9{2} = 'Baseline';
            legend(LH9,L9, 'Fontsize',20, 'Box','off');
        else
            errorbar(saccade_dist_per_subj{bn}(sub,:), prob_choice_favor_random{bn}(sub,:),prob_choice_favor_random_err{bn}(sub,:),'-or','LineWidth',2);        end
        hold on;
        if sub==3 || sub==8
            xlabel('Saccade length','Fontsize',30);
        end
        if sub==1 || sub==6
            ylabel({'Probability of chose in favor';' of accumulated evidence'},'Fontsize',30);
        end
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([0.4 1]);
        %         xlim([1 saccade_types(bn)]);
        yline(0.5,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        title(['Subject ' num2str(sub)],'Fontsize',20);
    end
    sgtitle(['Bias w.r.t saccade lengths of ' num2str(saccade_types(bn)) ' bins for all trials'],'fontsize',30);
end

%%
for ev_bn=1:length(evidence_types)
    fig10=figure();
    set(fig10,'defaultLegendAutoUpdate','off');
    for sub=1:num_sub
        subplot(2,5,sub)
        if sub==1
            LH9(1) = errorbar(binned_acc_evidence{ev_bn}(sub,:), binned_prob_bias{ev_bn}(sub,:),binned_prob_bias_err{ev_bn}(sub,:),'-ob','LineWidth',2);
            L9{1} = 'Subject';
        else
            errorbar(binned_acc_evidence{ev_bn}(sub,:), binned_prob_bias{ev_bn}(sub,:),binned_prob_bias_err{ev_bn}(sub,:),'-ob','LineWidth',2);
        end
        hold on;
        if sub==1
            LH9(2) = errorbar(binned_acc_evidence{ev_bn}(sub,:), binned_baseline{ev_bn}(sub,:),binned_baseline_err{ev_bn}(sub,:),'-or','LineWidth',2);
            L9{2} = 'Baseline';
            legend(LH9,L9, 'Fontsize',20, 'Box','off');
        else
            errorbar(binned_acc_evidence{ev_bn}(sub,:), binned_baseline{ev_bn}(sub,:),binned_baseline_err{ev_bn}(sub,:),'-or','LineWidth',2);
        end
        hold on;
        if sub==3 || sub==8
            xlabel('Evidence Accumulated','Fontsize',30);
        end
        if sub==1 || sub==6
            ylabel({'Probability of chose in favor';' of accumulated evidence'},'Fontsize',30);
        end
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([0.4 1]);
        yline(0.5,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        title(['Subject ' num2str(sub)],'Fontsize',20);
    end
    sgtitle(['Bias w.r.t evidence accumulated of ' num2str(evidence_types(ev_bn)) ' bins for all trials'],'fontsize',30);
end

%%
weight_bins = 3;
mn_weights_thresh = [];
prob_chose_in_favor_thresh_weights = [];
prob_chose_in_favor_thresh_weights_err = [];
prob_random_thresh_weights = [];
prob_random_thresh_weights_err = [];

mn_weights = [];
prob_chose_in_favor_weights = [];
prob_chose_in_favor_weights_err = [];
prob_random_weights = [];
prob_random_weights_err = [];

for sub=1:num_sub
    [sorted_thresh_saccade_weights, order_thresh_sort] = sort(thresh_saccade_weights{sub}(:));
    sorted_thresh_choice_saccades = thresh_choice_saccades{sub}(order_thresh_sort);
    sorted_thresh_random_landing_weights = thresh_random_landing_weights{sub}(order_thresh_sort);
    [sorted_saccade_weights, order_sort] = sort(saccade_weights{sub}(:));
    sorted_choice_saccades = choice_saccades{sub}(order_sort);
    sorted_random_landing_weights = random_landing_weights{sub}(order_sort);
    bin_edges_weights_thresh = linspace(sorted_thresh_saccade_weights(1),sorted_thresh_saccade_weights(end),weight_bins+1);
    bin_edges_weights = linspace(sorted_saccade_weights(1),sorted_saccade_weights(end),weight_bins+1);
    
    for bn=1:weight_bins
        indx_thresh = (sorted_thresh_saccade_weights>=bin_edges_weights_thresh(bn) & sorted_thresh_saccade_weights<bin_edges_weights_thresh(bn+1));
        mn_weights_thresh(sub,bn) = (bin_edges_weights_thresh(bn) + bin_edges_weights_thresh(bn+1))/2;
        prob_chose_in_favor_thresh_weights(sub,bn) = mean(sorted_thresh_choice_saccades(indx_thresh));
        prob_chose_in_favor_thresh_weights_err(sub,bn) = sqrt((prob_chose_in_favor_thresh_weights(bn)*(1-prob_chose_in_favor_thresh_weights(bn)))/length(sorted_thresh_choice_saccades(indx_thresh)));
        prob_random_thresh_weights(sub,bn) = mean(sorted_thresh_random_landing_weights(indx_thresh));
        prob_random_thresh_weights_err(sub,bn) = sqrt((prob_random_thresh_weights(bn)*(1-prob_random_thresh_weights(bn)))/length(sorted_thresh_random_landing_weights(indx_thresh)));
        
        mn_weights(sub,bn) = (bin_edges_weights(bn) + bin_edges_weights(bn+1))/2;
        indx = (sorted_saccade_weights>=bin_edges_weights(bn) & sorted_saccade_weights<bin_edges_weights(bn+1));
        prob_chose_in_favor_weights(sub,bn) = mean(sorted_choice_saccades(indx));
        prob_chose_in_favor_weights_err(sub,bn) = sqrt((prob_chose_in_favor_weights(bn)*(1-prob_chose_in_favor_weights(bn)))/length(sorted_choice_saccades(indx_thresh)));
        prob_random_weights(sub,bn) = mean(sorted_random_landing_weights(indx));
        prob_random_weights_err(sub,bn) = sqrt((prob_random_weights(bn)*(1-prob_random_weights(bn)))/length(sorted_random_landing_weights(indx)));
    end
end

mn_sp_weights_thresh = [];
prob_chose_in_favor_thresh_sp_weights = [];
prob_chose_in_favor_thresh_sp_weights_err = [];
prob_random_thresh_sp_weights = [];
prob_random_thresh_sp_weights_err = [];

mn_sp_weights = [];
prob_chose_in_favor_sp_weights = [];
prob_chose_in_favor_sp_weights_err = [];
prob_random_sp_weights = [];
prob_random_sp_weights_err = [];

for sub=1:num_sub
    [sorted_thresh_saccade_sp_weights, sp_order_thresh_sort] = sort(thresh_saccade_weights{sub}(:));
    sp_sorted_thresh_choice_saccades = thresh_choice_saccades{sub}(sp_order_thresh_sort);
    sp_sorted_thresh_random_landing_weights = thresh_random_landing_weights{sub}(sp_order_thresh_sort);
    [sp_sorted_saccade_weights, sp_order_sort] = sort(saccade_sp_weights{sub}(:));
    sp_sorted_choice_saccades = choice_saccades{sub}(sp_order_sort);
    sp_sorted_random_landing_weights = random_landing_weights{sub}(sp_order_sort);
    bin_edges_sp_weights_thresh = linspace(sorted_thresh_saccade_sp_weights(1),sorted_thresh_saccade_sp_weights(end),weight_bins+1);
    bin_edges_sp_weights = linspace(sp_sorted_saccade_weights(1),sp_sorted_saccade_weights(end),weight_bins+1);
    
    for bn=1:weight_bins
        sp_indx_thresh = (sorted_thresh_saccade_sp_weights>=bin_edges_sp_weights_thresh(bn) & sorted_thresh_saccade_sp_weights<bin_edges_sp_weights_thresh(bn+1));
        mn_sp_weights_thresh(sub,bn) = (bin_edges_sp_weights_thresh(bn) + bin_edges_sp_weights_thresh(bn+1))/2;
        prob_chose_in_favor_thresh_sp_weights(sub,bn) = mean(sp_sorted_thresh_choice_saccades(sp_indx_thresh));
        prob_chose_in_favor_thresh_sp_weights_err(sub,bn) = sqrt((prob_chose_in_favor_thresh_sp_weights(bn)*(1-prob_chose_in_favor_thresh_sp_weights(bn)))/length(sp_sorted_thresh_choice_saccades(sp_indx_thresh)));
        prob_random_thresh_sp_weights(sub,bn) = mean(sp_sorted_thresh_random_landing_weights(sp_indx_thresh));
        prob_random_thresh_sp_weights_err(sub,bn) = sqrt((prob_random_thresh_sp_weights(bn)*(1-prob_random_thresh_sp_weights(bn)))/length(sp_sorted_thresh_random_landing_weights(sp_indx_thresh)));
        
        mn_sp_weights(sub,bn) = (bin_edges_sp_weights(bn) + bin_edges_sp_weights(bn+1))/2;
        sp_indx = (sp_sorted_saccade_weights>=bin_edges_sp_weights(bn) & sp_sorted_saccade_weights<bin_edges_sp_weights(bn+1));
        prob_chose_in_favor_sp_weights(sub,bn) = mean(sp_sorted_choice_saccades(sp_indx));
        prob_chose_in_favor_sp_weights_err(sub,bn) = sqrt((prob_chose_in_favor_sp_weights(bn)*(1-prob_chose_in_favor_sp_weights(bn)))/length(sp_sorted_choice_saccades(sp_indx_thresh)));
        prob_random_sp_weights(sub,bn) = mean(sp_sorted_random_landing_weights(sp_indx));
        prob_random_sp_weights_err(sub,bn) = sqrt((prob_random_sp_weights(bn)*(1-prob_random_sp_weights(bn)))/length(sp_sorted_random_landing_weights(sp_indx)));
    end
end


%%
fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    if sub==1
        LH9(1) = errorbar(mn_weights_thresh(sub,:), prob_chose_in_favor_thresh_weights(sub,:),prob_chose_in_favor_thresh_weights_err(sub,:),'-ob','LineWidth',2);
        L9{1} = 'Subject';
    else
        errorbar(mn_weights_thresh(sub,:), prob_chose_in_favor_thresh_weights(sub,:),prob_chose_in_favor_thresh_weights_err(sub,:),'-ob','LineWidth',2);
    end
    hold on;
    if sub==1
        LH9(2) = errorbar(mn_weights_thresh(sub,:), prob_random_thresh_weights(sub,:),prob_random_thresh_weights_err(sub,:),'-or','LineWidth',2);
        L9{2} = 'Baseline';
        legend(LH9,L9, 'Fontsize',20, 'Box','off');
    else
        errorbar(mn_weights_thresh(sub,:), prob_random_thresh_weights(sub,:),prob_random_thresh_weights_err(sub,:),'-or','LineWidth',2);
    end
    hold on;
    if sub==3 || sub==8
        xlabel('Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor';' of accumulated evidence'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t weights for threshold trials'],'fontsize',30);


fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    if sub==1
        LH9(1) = errorbar(mn_weights(sub,:), prob_chose_in_favor_weights(sub,:),prob_chose_in_favor_weights_err(sub,:),'-ob','LineWidth',2);
        L9{1} = 'Subject';
    else
        errorbar(mn_weights(sub,:), prob_chose_in_favor_weights(sub,:),prob_chose_in_favor_weights_err(sub,:),'-ob','LineWidth',2);
    end
    hold on;
    if sub==1
        LH9(2) = errorbar(mn_weights(sub,:), prob_random_weights(sub,:),prob_random_weights_err(sub,:),'-or','LineWidth',2);
        L9{2} = 'Baseline';
        legend(LH9,L9, 'Fontsize',20, 'Box','off');
    else
        errorbar(mn_weights(sub,:), prob_random_weights(sub,:),prob_random_weights_err(sub,:),'-or','LineWidth',2);
    end
    hold on;
    if sub==3 || sub==8
        xlabel('Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor';' of accumulated evidence'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t weights for all trials'],'fontsize',30);

%%
fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    if sub==1
        LH9(1) = errorbar(mn_sp_weights(sub,:), prob_chose_in_favor_sp_weights(sub,:),prob_chose_in_favor_sp_weights_err(sub,:),'-ob','LineWidth',2);
        L9{1} = 'Subject';
    else
        errorbar(mn_sp_weights(sub,:), prob_chose_in_favor_sp_weights(sub,:),prob_chose_in_favor_sp_weights_err(sub,:),'-ob','LineWidth',2);
    end
    hold on;
    if sub==1
        LH9(2) = errorbar(mn_sp_weights(sub,:), prob_random_sp_weights(sub,:),prob_random_sp_weights_err(sub,:),'-or','LineWidth',2);
        L9{2} = 'Baseline';
        legend(LH9,L9, 'Fontsize',20, 'Box','off');
    else
        errorbar(mn_sp_weights(sub,:), prob_random_sp_weights(sub,:),prob_random_sp_weights_err(sub,:),'-or','LineWidth',2);
    end
    hold on;
    if sub==3 || sub==8
        xlabel('Spatial Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor';' of accumulated evidence'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t spatial weights for all trials'],'fontsize',30);


close all; clear all; clc;
warning('off');
fix_cond = 'saccade_code';%'eyelink';
boots = 50;
bin_pk = [1081 1921];
fix_cluster_dist = 240; % value used to cluster fixations closeby as one point!! prevents bias from artifacts!!
boots_thresh = 1; %bootstraps to get threshold, takes too much time so set to 1
threshold_performance = [0.5 0.75];
hprs_ridge = logspace(-3,3,7);
folds = 10;

subjects = ...
    {'OvalEyeMovement-subject04'; 'OvalEyeMovement-subject03';...
    'OvalEyeMovement-subject07'; 'OvalEyeMovement-subject06';...
    'OvalEyeMovement-subject08'; 'OvalEyeMovement-subject09';...
    'OvalEyeMovement-subject10'; 'OvalEyeMovement-subject11';...
    'OvalEyeMovement-subject12'; 'OvalEyeMovement-subject13'};
% subjects = {'OvalEyeMovement-subject12'; 'OvalEyeMovement-subject13'};

[num_sub,~] = size(subjects);
sanity_bnd = [0.5 0.55];
saccade_types = [1 2 3 5 7]; % number of saccade bins
evidence_types = [2 3 5 7];
saccade_types_sanity = 1;
saccade_types_thresh = 1;
expt_type = 1; % ratio change
peripheryPKbound = 13; % number of peripheral images to compute PK on based on distance from fixation
format short g;
disp_fig = 0;
screen_mid = [1920/2 1080/2];
standardize = 0;

disp('Starting search of best hyper parameters across all subjects...');
[best_hpr_ridge]= CustomRegression.combined_cross_validated_hprs(subjects, expt_type,...
    fix_cond,fix_cluster_dist,peripheryPKbound,hprs_ridge,standardize,folds);
disp(['Best ridge parameters across subjects is: ' num2str(best_hpr_ridge)]);
for sub=1:num_sub
    tic;
    disp(['Starting analysis for subject ' num2str(sub) ' ...']);
    [params, data, accuracy, num_image, subj_performance_bnd, mean_oval_signal,...
        num_trial, num_trial_sanity, num_trial_thresh,...
        acc_evidence, choice_in_favor, weights_on_saccades, sp_weights_on_saccades, random_landing,...
        acc_evidence_sanity, choice_in_favor_sanity, random_landing_sanity,...
        acc_evidence_thresh, choice_in_favor_thresh, weights_on_saccades_thresh, sp_weights_on_saccades_thresh, random_landing_thresh,...
        fixations_num, eyelink, saccade_code, fixation_dist,...
        saccades_used, saccades_used_sanity, saccades_used_thresh, saccade_dist_mid,...
        logits, spatial_pk, temporal_fixation_pk, temporal_duration_pk, temporal_pk_across_trials]...
        = analysis_across_allsaccades_functionalPK_SingleHyprs(subjects{sub}, expt_type, boots, fix_cond,...
        saccade_types, saccade_types_sanity, saccade_types_thresh, fix_cluster_dist, peripheryPKbound,...
        sanity_bnd, boots_thresh, threshold_performance, best_hpr_ridge, standardize);
    
    each_trial_duration = data.stim_duration;
    floor(sub) = subj_performance_bnd(1);
    thresh(sub) = subj_performance_bnd(2);
    
    eps = 0.1;
    scaling(sub,1) = prctile(params(:,1).^2,50);
    scaling(sub,2) = prctile(params(:,1).^2,50)-prctile(params(:,1).^2,16);
    scaling(sub,3) = prctile(params(:,1).^2,84)-prctile(params(:,1).^2,50);
    dist_param(sub,1) = prctile(exp(params(:,2)),50);
    dist_param(sub,2) = prctile(exp(params(:,2)),50)-prctile(exp(params(:,2)),16);
    dist_param(sub,3) = prctile(exp(params(:,2)),84)-prctile(exp(params(:,2)),50);
    mu(sub,1) = prctile(-1*(params(:,3).^2),50);
    mu(sub,2) = prctile(-1*(params(:,3).^2),50)-prctile(-1*(params(:,3).^2),16);
    mu(sub,3) = prctile(-1*(params(:,3).^2),84)-prctile(-1*(params(:,3).^2),50);
    sigma(sub,1) = prctile(params(:,4).^2 + eps,50);
    sigma(sub,2) = prctile(params(:,4).^2 +  eps,50)-prctile(params(:,4).^2 + eps,16);
    sigma(sub,3) = prctile(params(:,4).^2 + eps,84)-prctile(params(:,4).^2 + eps,50);
    temporal_fixation(sub,1) = prctile(params(:,5),50);
    temporal_fixation(sub,2) = prctile(params(:,5),50)-prctile(params(:,5),16);
    temporal_fixation(sub,3) = prctile(params(:,5),84)-prctile(params(:,5),50);
    bias(sub,1) = prctile(params(:,6),50);
    bias(sub,2) = prctile(params(:,6),50)-prctile(params(:,6),16);
    bias(sub,3) = prctile(params(:,6),84)-prctile(params(:,6),50);
    lapse(sub,1) = prctile(params(:,7).^2,50);
    lapse(sub,2) = prctile(params(:,7).^2,50)-prctile(params(:,7).^2,16);
    lapse(sub,3) = prctile(params(:,7).^2,84)-prctile(params(:,7).^2,50);
    %     lapse(sub,1) = prctile((1e-4+(1-1e-4)*sigmoid(params(:,7))),50);
    %     lapse(sub,2) = prctile((1e-4+(1-1e-4)*sigmoid(params(:,7))),50)-prctile((1e-4+(1-1e-4)*sigmoid(params(:,7))),16);
    %     lapse(sub,3) = prctile((1e-4+(1-1e-4)*sigmoid(params(:,7))),84)-prctile((1e-4+(1-1e-4)*sigmoid(params(:,7))),50);
    temporal_duration(sub,1) = prctile(params(:,8),50);
    temporal_duration(sub,2) = prctile(params(:,8),50)-prctile(params(:,8),16);
    temporal_duration(sub,3) = prctile(params(:,8),84)-prctile(params(:,8),50);
    
    trial_num(sub) = num_trial;
    trial_num_sanity(sub) = num_trial_sanity;
    trial_num_thresh(sub) = num_trial_thresh;
    fixations{sub} = fixations_num;
    max_fixations(sub) = max(fixations_num);
    eyelink_per_subj{sub} = eyelink;
    saccade_code_per_subj{sub} = saccade_code;
    fixation_dist_per_subj{sub} = fixation_dist;
    oval_signal{sub} = mean_oval_signal;
    
    spatial_kernel{sub} = spatial_pk;
    temporal_fixation_kernel{sub} = temporal_fixation_pk;
    temporal_duration_kernel{sub} = temporal_duration_pk;
    temporal_kernel_across_trials{sub} = temporal_pk_across_trials;
    log_bernoulli{sub} = logits;
    
    saccades_per_subj_sanity(sub) = saccades_used_sanity;
    saccades_per_subj_thresh(sub) = saccades_used_thresh;
    
    for bn=1:length(saccade_types)
        saccades_per_subj{bn}(sub,:) = saccades_used{bn};
        saccade_dist_per_subj{bn}(sub,:) = saccade_dist_mid{bn};
        for ss=1:saccade_types(bn)
            acc_evi_mid{bn}(sub,ss) = mean(acc_evidence{bn}{ss});
            acc_evi_mid_err{bn}(sub,ss) = sqrt(var(acc_evidence{bn}{ss})/length(acc_evidence{bn}{ss}));
            
            prob_choice_favor{bn}(sub,ss) = mean(choice_in_favor{bn}{ss});
            prob_choice_favor_err{bn}(sub,ss) = sqrt((prob_choice_favor{bn}(sub,ss) * (1 - prob_choice_favor{bn}(sub,ss)))/length(choice_in_favor{bn}{ss}));
            
            prob_choice_favor_random{bn}(sub,ss) = mean(random_landing{bn}{ss});
            prob_choice_favor_random_err{bn}(sub,ss) = sqrt(var(random_landing{bn}{ss})/length(random_landing{bn}{ss}));
        end
    end
    acc_evi_sanity(sub,1) = mean(acc_evidence_sanity{1});
    acc_evi_sanity(sub,2) = sqrt(var(acc_evidence_sanity{1})/length(acc_evidence_sanity{1}));
    
    prob_choice_favor_sanity(sub,1) = mean(choice_in_favor_sanity{1});
    prob_choice_favor_sanity(sub,2) = sqrt((prob_choice_favor_sanity(sub,1) * (1 - prob_choice_favor_sanity(sub,1)))/length(choice_in_favor_sanity{1}));
    
    prob_choice_favor_sanity_random(sub,1) = mean(random_landing_sanity{1});
    prob_choice_favor_sanity_random(sub,2) = sqrt(var(random_landing_sanity{1})/length(random_landing_sanity{1}));
    
    acc_evi_thresh(sub,1) = mean(acc_evidence_thresh{1});
    acc_evi_thresh(sub,2) = sqrt(var(acc_evidence_thresh{1})/length(acc_evidence_thresh{1}));
    
    prob_choice_favor_thresh(sub,1) = mean(choice_in_favor_thresh{1});
    prob_choice_favor_thresh(sub,2) = sqrt((prob_choice_favor_thresh(sub,1) * (1 - prob_choice_favor_thresh(sub,1)))/length(choice_in_favor_sanity{1}));
    
    prob_choice_favor_thresh_random(sub,1) = mean(random_landing_thresh{1});
    prob_choice_favor_thresh_random(sub,2) = sqrt(var(random_landing_thresh{1})/length(random_landing_thresh{1}));
    
    thresh_saccade_weights{sub} = weights_on_saccades_thresh{1};
    thresh_choice_saccades{sub} = choice_in_favor_thresh{1};
    thresh_random_landing_weights{sub} = random_landing_thresh{1};
    saccade_weights{sub} = weights_on_saccades{1}{1};
    choice_saccades{sub} = choice_in_favor{1}{1};
    random_landing_weights{sub} = random_landing{1}{1};
    
    thresh_saccade_sp_weights{sub} = sp_weights_on_saccades_thresh{1};
    saccade_sp_weights{sub} = sp_weights_on_saccades{1}{1};
    
    tr_ratio{sub} = data.true_ratio;
    ratio_used = data.ratio;
    choice{sub} = data.choice;
    
    disp(['Psychometric analysis for subject ' num2str(sub) ' ...']);
    ratios = linspace(0, 1, 101);%used in plots later
    [pm_fit, uniq_vals, yvals, stderrs] = GaborPsychometric(data, 1);
    subject_pm_curve_psig(sub,:) = (1-pm_fit.Fit(3)-pm_fit.Fit(4))*arrayfun(@(x) pm_fit.options.sigmoidHandle(x,pm_fit.Fit(1),pm_fit.Fit(2)), ratios)+pm_fit.Fit(4);
    
    disp('Analysis of bias w.r.t magnitude of evidence accumulated ...');
    [all_evidence_acc, evi_order] = sort(abs(acc_evidence{1}{1}));
    all_choice_subj = choice_in_favor{1}{1}(evi_order);
    all_prob_baseline = random_landing{1}{1}(evi_order);
    for ev_bn=1:length(evidence_types)
        evi_bin_edges = linspace(min(all_evidence_acc), max(all_evidence_acc),evidence_types(ev_bn)+1);
        for bb=1:evidence_types(ev_bn)
            binned_acc_evidence{ev_bn}(sub,bb) = mean(all_evidence_acc(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1)));
            binned_prob_bias{ev_bn}(sub,bb) = mean(all_choice_subj(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1)));
            binned_prob_bias_err{ev_bn}(sub,bb) = sqrt((binned_prob_bias{ev_bn}(sub,bb)*(1-binned_prob_bias{ev_bn}(sub,bb)))/length(all_choice_subj(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1))));
            binned_baseline{ev_bn}(sub,bb) = mean(all_prob_baseline(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1)));
            binned_baseline_err{ev_bn}(sub,bb) = sqrt(var(all_prob_baseline(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1)))/length(all_prob_baseline(all_evidence_acc>=evi_bin_edges(bb) & all_evidence_acc<=evi_bin_edges(bb+1))));
        end
    end
    
    disp(['All Analysis for subject ' num2str(sub) ' complete!']);
    disp('-----------------------------------------------------------------------------------------------------');
    toc;
    
    if disp_fig==1
        figure();
        subplot(2,3,1)
        pk = spatial_kernel{sub};
        imagesc(pk);set(gca,'YDir','normal'); h = colorbar; xlabel(h,'Spatial Weights','fontsize',8);h.LineWidth=1.5; h.FontSize=15;
        hold on;
        scatter(screen_mid(1), screen_mid(2), 50, '+', 'r' );
        ylabel('Screen y-axis','fontsize',20);
        xlabel('Screen x-axis','fontsize',20);
        xticks([1 fix(1920/2) 1920]);
        yticks([1 fix(1080/2) 1080]);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(2,3,2)
        temporal_fixation_pk = temporal_fixation_kernel{sub};
        plot(1:length(temporal_fixation_pk),temporal_fixation_pk,'-ob','LineWidth',2)
        hold on;
        ylabel('Fixation weights','fontsize',20);
        xlabel('Fixation Number','fontsize',20);
        xticks([1 fix(max_fixations(sub)/2) max_fixations(sub)]);
        xline(0.0,'k','LineWidth',0.75);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(2,3,3)
        temporal_duration_pk = temporal_duration_kernel{sub};
        plot(linspace(0,each_trial_duration,length(temporal_duration_pk)),temporal_duration_pk,'-ob','LineWidth',2)
        hold on;
        ylabel('Temporal weights','fontsize',20);
        xlabel('Time stamp of fixation','fontsize',20);
        xticks([0 each_trial_duration/2 each_trial_duration]);
        xline(0.0,'k','LineWidth',0.75);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(2,3,5)
        llo_mn = [];
        choice_mn = [];
        err_ch = [];
        err_ch_mn = [];
        llo_mean = [];
        choice_mean = [];
        bin_num = 50;
        [sorted_llo,order_llo] = sort(log_bernoulli{sub});
        choice_used = choice{sub}(order_llo);
        bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
        for bn=1:length(bin_edges)-1
            llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
            choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
            err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
        end
        [llo_mn, ord] = sort(llo_mean);
        choice_mn = choice_mean(ord);
        err_ch_mn = err_ch(ord);
        errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
        %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
        ylabel('Probability chose vertical','Fontsize',20);
        xlabel('Log likelihood odds','Fontsize',20);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([0.0 1]);
        xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
        yline(0.5,'-k','linewidth',2);
        hold on;
        xline(0.0,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(2,3,4)
        ratios = linspace(0, 1, 101);
        uniq_vals=linspace(-0.05,1.05,12);
        for tt=1:(length(uniq_vals)-1)
            subj_resp(sub,tt)=mean(choice{sub}(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)));
            ntrial_subj(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
        end
        uniq_vals=linspace(0,1,11);
        errorbar(uniq_vals,subj_resp(sub,:),sqrt((subj_resp(sub,:)).*(1-subj_resp(sub,:))./ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
        uniq_vals=linspace(-0.05,1.05,12);
        for tt=1:(length(uniq_vals)-1)
            pred_resp = lapse(1)*0.5+(1-lapse(1))*sigmoid(log_bernoulli{sub}(((tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)))));
            subj_resp_pred(sub,tt)=mean(pred_resp);
            ntrial_subj_pred(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
            subj_resp_pred_err(sub,tt) = sqrt(subj_resp_pred(sub,tt)*(1-subj_resp_pred(sub,tt))/ntrial_subj_pred(sub,tt));
        end
        uniq_vals=linspace(0,1,11);
        hold on;
        errorbar(uniq_vals,subj_resp_pred(sub,:),subj_resp_pred_err(sub,:),'--or','Linewidth',2);
        ylabel('Probability chose vertical','Fontsize',20);
        xlabel({'Ratio of vertical ellipses'},'Fontsize',20);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([0.0 1]);
        yline(0.5,'-k','linewidth',2);
        hold on;
        xline(0.5,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(2,3,6)
        bar(1, prob_choice_favor{1}(sub),'FaceColor','b','EdgeColor','k','LineWidth',0.75);
        hold on;
        bar(2, prob_choice_favor_random{1}(sub),'FaceColor','r','EdgeColor','k','LineWidth',0.75);
        hold on;
        errorbar([1 2],[prob_choice_favor{1}(sub) prob_choice_favor_random{1}(sub)],[squeeze(prob_choice_favor_err{1}(sub)) squeeze(prob_choice_favor_random_err{1}(sub))],'ok','LineWidth',2,'linestyle','none');
        xlim([0.25 2.75]);
        text(1,0.9,['Trial Num: ' num2str(trial_num(sub))],'Fontsize',12);
        hold on;
        text(1,0.85,['Sacc Num: ' num2str(sum(fixations{sub}))],'Fontsize',12);
        xticks([1 2]);
        xticklabels({'Subject'; 'Baseline'});
        ylabel('Probability of chose in favor','Fontsize',20);
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([0.4 1]);
        yline(0.5,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        sgtitle(['Subject ' num2str(sub)],'fontsize',30);
        
        %         pause;
        close all;
    end
end
%%
disp('Plotting stuff may take a while ...............................');
for sub=1:num_sub
    fig = figure(sub);
    set(fig,'defaultLegendAutoUpdate','off');
    
    subplot(2,4,5)
    %     [floor(sub), thresh(sub), ~] = getThresholdWindow(data,1, 0.5, 0.7);
    plot(ratio_used,'-ob');
    xlabel('Trial number','Fontsize',20);
    ylabel({'Ratio of vertical ellipses'},'Fontsize',20);
    yline(0.5,'k','LineWidth',1.5);
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.5 1.0]);
    xlim([1 length(ratio_used)]);
    %     yline(thresh(sub),'-k','linewidth',2)
    %     hold on;
    %     yline(floor(sub),'-k','linewidth',2)
    %     yline(0.5,'-k','linewidth',2)
    %     ylim([0 1]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Staircase','fontsize',20);
    
    subplot(2,4,1)
    pk = spatial_kernel{sub};
    imagesc(pk);set(gca,'YDir','normal'); h = colorbar; xlabel(h,'Spatial Weights','fontsize',8);h.LineWidth=1.5; h.FontSize=15;
    hold on;
    scatter(screen_mid(1), screen_mid(2), 50, '+', 'r' );
    ylabel('Screen y-axis','fontsize',20);
    xlabel('Screen x-axis','fontsize',20);
    xticks([1 fix(1920/2) 1920]);
    yticks([1 fix(1080/2) 1080]);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Spatial Weights','fontsize',20);
    
    
    subplot(2,4,2)
    temporal_fixation_pk = temporal_fixation_kernel{sub};
    plot(1:length(temporal_fixation_pk),temporal_fixation_pk,'-ob','LineWidth',2)
    hold on;
    ylabel('Fixation weights','fontsize',20);
    xlabel('Fixation Number','fontsize',20);
    xticks([1 fix(max_fixations(sub)/2) max_fixations(sub)]);
    xlim([1 max_fixations(sub)]);
    xline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Fixation weights','Fontsize',20);
    
    subplot(2,4,3)
    temporal_duration_pk = temporal_duration_kernel{sub};
    plot(linspace(0,each_trial_duration,length(temporal_duration_pk)),temporal_duration_pk,'-ob','LineWidth',2)
    hold on;
    ylabel('Temporal weights','fontsize',20);
    xlabel('Timestamps','fontsize',20);
    xticks([0 each_trial_duration/2 each_trial_duration]);
    xline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Timestamp weights','Fontsize',20)
    
    subplot(2,4,4)
    temporal_pk = temporal_kernel_across_trials{sub};
    plot(1:length(temporal_pk),temporal_pk,'-ob','LineWidth',2)
    hold on;
    ylabel('Combined temporal weights','fontsize',20);
    xlabel('Fixation Number','fontsize',20);
    xticks([1 fix(max_fixations(sub)/2) max_fixations(sub)]);
    xlim([1 max_fixations(sub)]);
    xline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Combined temporal','Fontsize',20);
    
    subplot(2,4,6)
    llo_mn = [];
    choice_mn = [];
    err_ch = [];
    err_ch_mn = [];
    llo_mean = [];
    choice_mean = [];
    bin_num = 20;
    [sorted_llo,order_llo] = sort(log_bernoulli{sub});
    choice_used = choice{sub}(order_llo);
    bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
    for bn=1:length(bin_edges)-1
        llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
    end
    [llo_mn, ord] = sort(llo_mean);
    choice_mn = choice_mean(ord);
    err_ch_mn = err_ch(ord);
    errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
    %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
    ylabel('Probability chose vertical','Fontsize',20);
    xlabel('Log likelihood odds','Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('LLO Sanity','fontsize',20);
    
    subplot(2,4,7)
    ratios = linspace(0, 1, 101);
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp(sub,tt)=mean(choice{sub}(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)));
        ntrial_subj(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    errorbar(uniq_vals,subj_resp(sub,:),sqrt((subj_resp(sub,:)).*(1-subj_resp(sub,:))./ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        pred_resp = lapse(1)*0.5+(1-lapse(1))*sigmoid(log_bernoulli{sub}(((tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)))));
        subj_resp_pred(sub,tt) = mean(pred_resp);
        ntrial_subj_pred(sub,tt) = sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
        subj_resp_pred_err(sub,tt) = sqrt(subj_resp_pred(sub,tt)*(1-subj_resp_pred(sub,tt))/ntrial_subj_pred(sub,tt));
    end
    uniq_vals=linspace(0,1,11);
    hold on;
    errorbar(uniq_vals,subj_resp_pred(sub,:),subj_resp_pred_err(sub,:),'--or','Linewidth',2);
    ylabel('Probability chose vertical','Fontsize',20);
    xlabel({'Ratio of vertical ellipses'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('Match to data','fontsize',20);
    
    subplot(2,4,8)
    hold on;
    LH3(1)=bar(1, prob_choice_favor{1}(sub),'FaceColor','b','EdgeColor','k','LineWidth',0.75);
    L3{1} = 'Subject';
    hold on;
    LH3(2)=bar(2, prob_choice_favor_random{1}(sub),'FaceColor','r','EdgeColor','k','LineWidth',0.75);
    L3{2} = 'Baseline';
    hold on;
    errorbar([1 2],[prob_choice_favor{1}(sub) prob_choice_favor_random{1}(sub)],[squeeze(prob_choice_favor_err{1}(sub)) squeeze(prob_choice_favor_random_err{1}(sub))],'ok','LineWidth',2,'linestyle','none');
    xlim([0.25 2.75]);
    text(1,0.9,['Trial Num: ' num2str(trial_num(sub))],'Fontsize',12);
    hold on;
    text(1,0.85,['Sacc Num: ' num2str(sum(fixations{sub}))],'Fontsize',12);
    xticks([1 2]);
    xticklabels({'Subject'; 'Baseline'});
    ylabel('Probability of chose in favor','Fontsize',20);
    hold on;
    legend(LH3,L3, 'Fontsize',15, 'Box','off');
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title('External CB','fontsize',20);
    
    sgtitle(['All analysis for Subject ' num2str(sub)],'fontsize',30);
end

%%
figure();
for sub=1:num_sub
    subplot(2,5,sub)
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp(sub,tt)=mean(choice{sub}(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)));
        ntrial_subj(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    LH1(1) = errorbar(uniq_vals,subj_resp(sub,:),sqrt((subj_resp(sub,:)).*(1-subj_resp(sub,:))./ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
    L1{1} = 'Data';
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        pred_resp = lapse(1)*0.5+(1-lapse(1))*sigmoid(log_bernoulli{sub}(((tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)))));
        subj_resp_pred(sub,tt)=mean(pred_resp);
        subj_resp_pred_err(sub,tt) = sqrt(subj_resp_pred(sub,tt)*(1-subj_resp_pred(sub,tt))/ntrial_subj_pred(sub,tt));
        ntrial_subj_pred(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    hold on;
    LH1(2) = errorbar(uniq_vals,subj_resp_pred(sub,:),subj_resp_pred_err(sub,:),'-or','Linewidth',2);
    L1{2} = 'Fit with PK';
    if sub==1 || sub==6
        ylabel('Probability chose vertical','Fontsize',20);
    end
    if sub==3 || sub==8
        xlabel({'Ratio of vertical ellipses'},'Fontsize',20);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    if sub==1
        legend(LH1,L1, 'Fontsize',12, 'Box','off','Location','northwest');
    end
    title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('PS fit w.r.t learned Weights/PK compared to real data','fontsize',30)
%%
bin_num=20;
figure();
llo_mn = [];
choice_mn = [];
err_ch = [];
err_ch_mn = [];
llo_mean = [];
choice_mean = [];
for sub=1:num_sub
    subplot(2,5,sub)
    %     llo = llo( ~any( isnan( llo ) | isinf( llo ), 2 ),: );
    [sorted_llo,order_llo] = sort(log_bernoulli{sub});
    choice_used = choice{sub}(order_llo);
    bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
    for bn=1:length(bin_edges)-1
        llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
    end
    [llo_mn, ord] = sort(llo_mean);
    choice_mn = choice_mean(ord);
    err_ch_mn = err_ch(ord);
    errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
    %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
    if sub==1 || sub==6
        ylabel('Probability chose vertical','Fontsize',20);
    end
    if sub==3 || sub==8
        xlabel('Log likelihood odds','Fontsize',20);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Check how good weights predict behavior','fontsize',30);
%%
figure();
ratios = linspace(0, 1, 101);
for sub=1:num_sub
    plot(ratios, subject_pm_curve_psig(sub,:), 'color','k','LineWidth', 0.5);
    hold on;
end
yline(0.5,'-k','linewidth',2);
hold on;
xline(0.5,'-k','linewidth',2);
hold on;
ylabel('Probability chose vertical','Fontsize',20);
xlabel({'Ratio of vertical ellipses'},'Fontsize',20);
hold on;
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ylim([0.0 1]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
plot(ratios, mean(subject_pm_curve_psig,1), 'k','LineWidth', 2);
title('Psychometric Curves for all subjects','fontsize',30);
%%
figure();
max_pk = 0;
min_pk = 100000;
max_fx_pk = 0;
min_fx_pk = 100000;
max_dur_pk = 0;
min_dur_pk = 100000;
max_comb_pk = 0;
min_comb_pk = 100000;
for sub=1:num_sub
    max_pk = max(max_pk,max(spatial_kernel{sub}(:)));
    min_pk = min(min_pk,min(spatial_kernel{sub}(:)));
    max_fx_pk = max(max_fx_pk,max(temporal_fixation_kernel{sub}(:) * scaling(sub,1)));
    min_fx_pk = min(min_fx_pk,min(temporal_fixation_kernel{sub}(:) * scaling(sub,1))) ;
    max_dur_pk = max(max_dur_pk,max(temporal_duration_kernel{sub}(:)));% * scaling(sub);
    min_dur_pk = min(min_dur_pk,min(temporal_duration_kernel{sub}(:)));% * scaling(sub);
    max_comb_pk = max(max_comb_pk,max(temporal_kernel_across_trials{sub}(:)));
    min_comb_pk = min(min_comb_pk,min(temporal_kernel_across_trials{sub}(:)));
end
for sub=1:num_sub
    subplot(3,4,sub)
    pk = spatial_kernel{sub};%./scaling(sub,1);
    imagesc(pk);set(gca,'YDir','normal'); h = colorbar; xlabel(h,'Spatial Weights','fontsize',8);h.LineWidth=1.5; h.FontSize=15;
    axis image;
    hold on;
    scatter(screen_mid(1), screen_mid(2), 50, '+', 'r' );
    if sub==1 || sub==5 || sub==9
        ylabel('Screen y-axis','fontsize',20);
    end
    if sub==2 || sub==6 || sub==10
        xlabel('Screen x-axis','fontsize',20);
    end
    xticks([1 fix(1920/2) 1920]);
    yticks([1 fix(1080/2) 1080]);
    hold on;
    caxis([min_pk, max_pk]);
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Spatial weights w.r.t fixation in center of screen','fontsize',30);

figure();
for sub=1:num_sub
    %     subplot(3,4,sub)
    temporal_fixation_pk = temporal_fixation_kernel{sub} * scaling(sub,1);
    plot(1:length(temporal_fixation_pk),temporal_fixation_pk,'-ob','LineWidth',2)
    hold on;
    if sub==1 || sub==5 || sub==9
        ylabel('Temporal Weights','fontsize',20);
    end
    if sub==2 || sub==6 || sub==10
        xlabel('Fixation Number','fontsize',20);
    end
    %     if sub==7
    %         ylim([0 150]);
    %     else
    %         ylim([0 5]);
    %     end
    %     xticks([1 fix(max_fixations(sub)/2) max_fixations(sub)]);
    xlim([1 max([max_fixations(1) max_fixations(3:end)])]);
    xticks(1:max([max_fixations(1) max_fixations(3:end)]));
    ylim([0.0 max_fx_pk]);
    %     ylim([min_fx_pk max_fx_pk]);
    xline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    %     title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Temporal fixation based weights','fontsize',30);

figure();
for sub=1:num_sub
    %     subplot(3,4,sub)
    temporal_duration_pk = temporal_duration_kernel{sub};% * scaling(sub);
    plot(linspace(0,each_trial_duration,length(temporal_duration_pk)),temporal_duration_pk,'-ob','LineWidth',2)
    hold on;
    if sub==1 || sub==5 || sub==9
        ylabel('Temporal weights','fontsize',20);
    end
    if sub==2 || sub==6 || sub==10
        xlabel('Time stamp of fixation','fontsize',20);
    end
    xticks([0 each_trial_duration/2 each_trial_duration]);
    xlim([0 each_trial_duration]);
    ylim([0.0 max_dur_pk]);
    %     ylim([min_dur_pk max_dur_pk]);
    xline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    %     title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Time stamp based weights','fontsize',30)

figure();
for sub=1:num_sub
    subplot(3,4,sub)
    temporal_pk = temporal_kernel_across_trials{sub};
    plot(1:length(temporal_pk),temporal_pk,'-ob','LineWidth',2)
    hold on;
    if sub==1 || sub==5 || sub==9
        ylabel('Combined temporal weights','fontsize',20);
    end
    if sub==2 || sub==6 || sub==10
        xlabel('Fixation Number','fontsize',20);
    end
    xticks([1 fix(max_fixations(sub)/2) max_fixations(sub)]);
    xline(0.0,'k','LineWidth',0.75);
    xlim([1 max_fixations(sub)]);
    ylim([min_comb_pk max_comb_pk]);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end
sgtitle('Combined fixation and time stamp based weights','fontsize',30);

%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, prob_choice_favor{1}(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
        L6{1} = 'Subject';
    else
        bar(k1, prob_choice_favor{1}(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, prob_choice_favor_random{1}(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
        L6{2} = 'Baseline';
    else
        bar(k2, prob_choice_favor_random{1}(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    errorbar([k1 k2],[prob_choice_favor{1}(sub) prob_choice_favor_random{1}(sub)],[prob_choice_favor_err{1}(sub) prob_choice_favor_random_err{1}(sub)],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.25,0.9,[num2str(trial_num(sub)) ' trials'],'Fontsize',15);
    text(k1-0.25,0.8,[num2str(saccades_per_subj{1}(sub)) ' sac.'],'Fontsize',15);
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of chose';' in favor of accumulated evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 2.5;
    k2 = k2 + 2.5;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k+0.5);
legend(LH6,L6, 'Fontsize',20, 'Box','off');
sgtitle('All saccades for all subjects','fontsize',30);

%%
%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    
    bar(k1, prob_choice_favor{1}(sub)-prob_choice_favor_random{1}(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    hold on;
    errorbar([k1],[prob_choice_favor{1}(sub)-prob_choice_favor_random{1}(sub)],[prob_choice_favor_err{1}(sub)+prob_choice_favor_random_err{1}(sub)],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.5,0.175,[num2str(trial_num(sub)) ' trials'],'Fontsize',15);
    text(k1-0.5,0.15,[num2str(saccades_per_subj{1}(sub)) ' sac.'],'Fontsize',15);
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of chose';' in favor of accumulated evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([-0.05 0.2]);
    yline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 1;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k);
% legend(LH6,L6, 'Fontsize',20, 'Box','off');
sgtitle('All saccades for all subjects','fontsize',30);


%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, prob_choice_favor_sanity(sub,1),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
        L6{1} = 'Subject';
    else
        bar(k1, prob_choice_favor_sanity(sub,1),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, prob_choice_favor_sanity_random(sub,1),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
        L6{2} = 'Baseline';
    else
        bar(k2, prob_choice_favor_sanity_random(sub,1),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    errorbar([k1 k2],[prob_choice_favor_sanity(sub,1) prob_choice_favor_sanity_random(sub,1)],[prob_choice_favor_sanity(sub,2) prob_choice_favor_sanity_random(sub,2)],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.25,0.9,[num2str(trial_num_sanity(sub)) ' trials'],'Fontsize',15)
    text(k1-0.25,0.8,[num2str(saccades_per_subj_sanity(sub)) ' sac.'],'Fontsize',15);
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of chose';' in favor of accumulated evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 2.5;
    k2 = k2 + 2.5;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k+0.5);
legend(LH6,L6, 'Fontsize',20, 'Box','off');
sgtitle(['Sanity check figure with ratio > ' num2str(min(sanity_bnd(1),sanity_bnd(2))) ' and < ' num2str(max(sanity_bnd(1),sanity_bnd(2)))],'fontsize',30);

%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, prob_choice_favor_thresh(sub,1),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
        L6{1} = 'Subject';
    else
        bar(k1, prob_choice_favor_thresh(sub,1),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, prob_choice_favor_thresh_random(sub,1),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
        L6{2} = 'Baseline';
    else
        bar(k2, prob_choice_favor_thresh_random(sub,1),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    errorbar([k1 k2],[prob_choice_favor_thresh(sub,1) prob_choice_favor_thresh_random(sub,1)],[prob_choice_favor_thresh(sub,2) prob_choice_favor_thresh_random(sub,2)],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.25,0.9,[num2str(trial_num_thresh(sub)) ' trials'],'Fontsize',15)
    text(k1-0.25,0.75,[num2str(saccades_per_subj_thresh(sub)) ' sac.'],'Fontsize',15);
    text(k1-0.25,0.85,'Thresh ratio: ','Fontsize',15);
    text(k1-0.25,0.82, [num2str(floor(sub),'%1.2f') '     ' num2str(thresh(sub),'%1.2f')],'Fontsize',15);
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of chose';' in favor of accumulated evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 2.5;
    k2 = k2 + 2.5;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k+0.5);
legend(LH6,L6, 'Fontsize',20, 'Box','off');
sgtitle('Saccade bias of all subjects for threshold trials','fontsize',30);

%%
for bn=2:length(saccade_types)
    fig=figure();
    set(fig,'defaultLegendAutoUpdate','off');
    for sub=1:num_sub
        subplot(2,5,sub)
        if sub==1
            %             LH9(1) = errorbar(1:saccade_types(bn), prob_choice_favor{bn}(sub,:),prob_choice_favor_err{bn}(sub,:),'-ob','LineWidth',2);
            LH9(1) = errorbar(saccade_dist_per_subj{bn}(sub,:), prob_choice_favor{bn}(sub,:),prob_choice_favor_err{bn}(sub,:),'-ob','LineWidth',2);
            L9{1} = 'Subject';
        else
            %             errorbar(1:saccade_types(bn), prob_choice_favor{bn}(sub,:),prob_choice_favor_err{bn}(sub,:),'-ob','LineWidth',2);
            errorbar(saccade_dist_per_subj{bn}(sub,:), prob_choice_favor{bn}(sub,:),prob_choice_favor_err{bn}(sub,:),'-ob','LineWidth',2);
        end
        hold on;
        if sub==1
            %             LH9(2) = errorbar(1:saccade_types(bn), prob_choice_favor_random{bn}(sub,:),prob_choice_favor_random_err{bn}(sub,:),'-or','LineWidth',2);
            LH9(2) = errorbar(saccade_dist_per_subj{bn}(sub,:), prob_choice_favor_random{bn}(sub,:),prob_choice_favor_random_err{bn}(sub,:),'-or','LineWidth',2);
            L9{2} = 'Baseline';
            legend(LH9,L9, 'Fontsize',20, 'Box','off');
        else
            errorbar(saccade_dist_per_subj{bn}(sub,:), prob_choice_favor_random{bn}(sub,:),prob_choice_favor_random_err{bn}(sub,:),'-or','LineWidth',2);        
        end
        hold on;
        if sub==3 || sub==8
            xlabel('Saccade length','Fontsize',30);
        end
        if sub==1 || sub==6
            ylabel({'Probability of chose in favor';' of accumulated evidence'},'Fontsize',30);
        end
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([0.4 1]);
        %         xlim([1 saccade_types(bn)]);
        yline(0.5,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        title(['Subject ' num2str(sub)],'Fontsize',20);
    end
    sgtitle(['Bias w.r.t saccade lengths of ' num2str(saccade_types(bn)) ' bins for all trials'],'fontsize',30);
end

%%
for bn=2:length(saccade_types)
    fig=figure();
    set(fig,'defaultLegendAutoUpdate','off');
    for sub=1:num_sub
        subplot(2,5,sub)
        bar(1:saccade_types(bn), prob_choice_favor{bn}(sub,:)-prob_choice_favor_random{bn}(sub,:),'b');  
        hold on;
        errorbar(1:saccade_types(bn), prob_choice_favor{bn}(sub,:)-prob_choice_favor_random{bn}(sub,:),sqrt(prob_choice_favor_err{bn}(sub,:).^2+prob_choice_favor_err{bn}(sub,:).^2),'ok','LineWidth',2);
        hold on;
        if sub==3 || sub==8
            xlabel('Saccade length','Fontsize',30);
        end
        if sub==1 || sub==6
            ylabel({'Probability of chose in favor'},'Fontsize',30);
        end
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([-0.5 0.5]);
        xlim([0.25 saccade_types(bn)+0.75]);
        yline(0.0,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        title(['Subject ' num2str(sub)],'Fontsize',20);
    end
    sgtitle(['Bias w.r.t saccade lengths of ' num2str(saccade_types(bn)) ' bins for all trials difference plot'],'fontsize',30);
end


%%
for ev_bn=1:length(evidence_types)
    fig10=figure();
    set(fig10,'defaultLegendAutoUpdate','off');
    for sub=1:num_sub
        subplot(2,5,sub)
        if sub==1
            LH9(1) = errorbar(binned_acc_evidence{ev_bn}(sub,:), binned_prob_bias{ev_bn}(sub,:),binned_prob_bias_err{ev_bn}(sub,:),'-ob','LineWidth',2);
            L9{1} = 'Subject';
        else
            errorbar(binned_acc_evidence{ev_bn}(sub,:), binned_prob_bias{ev_bn}(sub,:),binned_prob_bias_err{ev_bn}(sub,:),'-ob','LineWidth',2);
        end
        hold on;
        if sub==1
            LH9(2) = errorbar(binned_acc_evidence{ev_bn}(sub,:), binned_baseline{ev_bn}(sub,:),binned_baseline_err{ev_bn}(sub,:),'-or','LineWidth',2);
            L9{2} = 'Baseline';
            legend(LH9,L9, 'Fontsize',20, 'Box','off');
        else
            errorbar(binned_acc_evidence{ev_bn}(sub,:), binned_baseline{ev_bn}(sub,:),binned_baseline_err{ev_bn}(sub,:),'-or','LineWidth',2);
        end
        hold on;
        if sub==3 || sub==8
            xlabel('Evidence Accumulated','Fontsize',30);
        end
        if sub==1 || sub==6
            ylabel({'Probability of chose in favor';' of accumulated evidence'},'Fontsize',30);
        end
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([0.4 1]);
        yline(0.5,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        title(['Subject ' num2str(sub)],'Fontsize',20);
    end
    sgtitle(['Bias w.r.t evidence accumulated of ' num2str(evidence_types(ev_bn)) ' bins for all trials'],'fontsize',30);
end

%%
for ev_bn=1:length(evidence_types)
    fig10=figure();
    set(fig10,'defaultLegendAutoUpdate','off');
    for sub=1:num_sub
        subplot(2,5,sub)
        bar(1:evidence_types(ev_bn), binned_prob_bias{ev_bn}(sub,:)-binned_baseline{ev_bn}(sub,:),'b');
        hold on;
        errorbar(1:evidence_types(ev_bn), binned_prob_bias{ev_bn}(sub,:)-binned_baseline{ev_bn}(sub,:),sqrt(binned_prob_bias_err{ev_bn}(sub,:).^2+binned_prob_bias_err{ev_bn}(sub,:).^2),'ok','LineWidth',2);
        hold on;
        if sub==3 || sub==8
            xlabel('Evidence Accumulated','Fontsize',30);
        end
        if sub==1 || sub==6
            ylabel({'Probability of chose in favor';' of accumulated evidence'},'Fontsize',30);
        end
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ylim([-0.5 0.5]);
        yline(0.0,'-k','linewidth',2);
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        title(['Subject ' num2str(sub)],'Fontsize',20);
    end
    sgtitle(['Bias w.r.t evidence accumulated of ' num2str(evidence_types(ev_bn)) ' bins for all trials difference plot'],'fontsize',30);
end




%%
weight_bins = 2;
mn_weights_thresh = [];
prob_chose_in_favor_thresh_weights = [];
prob_chose_in_favor_thresh_weights_err = [];
prob_random_thresh_weights = [];
prob_random_thresh_weights_err = [];

mn_weights = [];
prob_chose_in_favor_weights = [];
prob_chose_in_favor_weights_err = [];
prob_random_weights = [];
prob_random_weights_err = [];

for sub=1:num_sub
    [sorted_thresh_saccade_weights, order_thresh_sort] = sort(thresh_saccade_weights{sub}(:));
    sorted_thresh_choice_saccades = thresh_choice_saccades{sub}(order_thresh_sort);
    sorted_thresh_random_landing_weights = thresh_random_landing_weights{sub}(order_thresh_sort);
    [sorted_saccade_weights, order_sort] = sort(saccade_weights{sub}(:));
    sorted_choice_saccades = choice_saccades{sub}(order_sort);
    sorted_random_landing_weights = random_landing_weights{sub}(order_sort);
    bin_edges_weights_thresh = linspace(sorted_thresh_saccade_weights(1),sorted_thresh_saccade_weights(end),weight_bins+1);
    bin_edges_weights = linspace(sorted_saccade_weights(1),sorted_saccade_weights(end),weight_bins+1);
    
    for bn=1:weight_bins
        indx_thresh = (sorted_thresh_saccade_weights>=bin_edges_weights_thresh(bn) & sorted_thresh_saccade_weights<bin_edges_weights_thresh(bn+1));
        mn_weights_thresh(sub,bn) = (bin_edges_weights_thresh(bn) + bin_edges_weights_thresh(bn+1))/2;
        prob_chose_in_favor_thresh_weights(sub,bn) = mean(sorted_thresh_choice_saccades(indx_thresh));
        prob_chose_in_favor_thresh_weights_err(sub,bn) = sqrt((prob_chose_in_favor_thresh_weights(bn)*(1-prob_chose_in_favor_thresh_weights(bn)))/length(sorted_thresh_choice_saccades(indx_thresh)));
        prob_random_thresh_weights(sub,bn) = mean(sorted_thresh_random_landing_weights(indx_thresh));
        prob_random_thresh_weights_err(sub,bn) = sqrt((prob_random_thresh_weights(bn)*(1-prob_random_thresh_weights(bn)))/length(sorted_thresh_random_landing_weights(indx_thresh)));
        
        mn_weights(sub,bn) = (bin_edges_weights(bn) + bin_edges_weights(bn+1))/2;
        indx = (sorted_saccade_weights>=bin_edges_weights(bn) & sorted_saccade_weights<bin_edges_weights(bn+1));
        prob_chose_in_favor_weights(sub,bn) = mean(sorted_choice_saccades(indx));
        prob_chose_in_favor_weights_err(sub,bn) = sqrt((prob_chose_in_favor_weights(bn)*(1-prob_chose_in_favor_weights(bn)))/length(sorted_choice_saccades(indx_thresh)));
        prob_random_weights(sub,bn) = mean(sorted_random_landing_weights(indx));
        prob_random_weights_err(sub,bn) = sqrt((prob_random_weights(bn)*(1-prob_random_weights(bn)))/length(sorted_random_landing_weights(indx)));
    end
end

mn_sp_weights_thresh = [];
prob_chose_in_favor_thresh_sp_weights = [];
prob_chose_in_favor_thresh_sp_weights_err = [];
prob_random_thresh_sp_weights = [];
prob_random_thresh_sp_weights_err = [];

mn_sp_weights = [];
prob_chose_in_favor_sp_weights = [];
prob_chose_in_favor_sp_weights_err = [];
prob_random_sp_weights = [];
prob_random_sp_weights_err = [];

for sub=1:num_sub
    [sorted_thresh_saccade_sp_weights, sp_order_thresh_sort] = sort(thresh_saccade_weights{sub}(:));
    sp_sorted_thresh_choice_saccades = thresh_choice_saccades{sub}(sp_order_thresh_sort);
    sp_sorted_thresh_random_landing_weights = thresh_random_landing_weights{sub}(sp_order_thresh_sort);
    [sp_sorted_saccade_weights, sp_order_sort] = sort(saccade_sp_weights{sub}(:));
    sp_sorted_choice_saccades = choice_saccades{sub}(sp_order_sort);
    sp_sorted_random_landing_weights = random_landing_weights{sub}(sp_order_sort);
    bin_edges_sp_weights_thresh = linspace(sorted_thresh_saccade_sp_weights(1),sorted_thresh_saccade_sp_weights(end),weight_bins+1);
    bin_edges_sp_weights = linspace(sp_sorted_saccade_weights(1),sp_sorted_saccade_weights(end),weight_bins+1);
    
    for bn=1:weight_bins
        sp_indx_thresh = (sorted_thresh_saccade_sp_weights>=bin_edges_sp_weights_thresh(bn) & sorted_thresh_saccade_sp_weights<bin_edges_sp_weights_thresh(bn+1));
        mn_sp_weights_thresh(sub,bn) = (bin_edges_sp_weights_thresh(bn) + bin_edges_sp_weights_thresh(bn+1))/2;
        prob_chose_in_favor_thresh_sp_weights(sub,bn) = mean(sp_sorted_thresh_choice_saccades(sp_indx_thresh));
        prob_chose_in_favor_thresh_sp_weights_err(sub,bn) = sqrt((prob_chose_in_favor_thresh_sp_weights(bn)*(1-prob_chose_in_favor_thresh_sp_weights(bn)))/length(sp_sorted_thresh_choice_saccades(sp_indx_thresh)));
        prob_random_thresh_sp_weights(sub,bn) = mean(sp_sorted_thresh_random_landing_weights(sp_indx_thresh));
        prob_random_thresh_sp_weights_err(sub,bn) = sqrt((prob_random_thresh_sp_weights(bn)*(1-prob_random_thresh_sp_weights(bn)))/length(sp_sorted_thresh_random_landing_weights(sp_indx_thresh)));
        
        mn_sp_weights(sub,bn) = (bin_edges_sp_weights(bn) + bin_edges_sp_weights(bn+1))/2;
        sp_indx = (sp_sorted_saccade_weights>=bin_edges_sp_weights(bn) & sp_sorted_saccade_weights<bin_edges_sp_weights(bn+1));
        prob_chose_in_favor_sp_weights(sub,bn) = mean(sp_sorted_choice_saccades(sp_indx));
        prob_chose_in_favor_sp_weights_err(sub,bn) = sqrt((prob_chose_in_favor_sp_weights(bn)*(1-prob_chose_in_favor_sp_weights(bn)))/length(sp_sorted_choice_saccades(sp_indx_thresh)));
        prob_random_sp_weights(sub,bn) = mean(sp_sorted_random_landing_weights(sp_indx));
        prob_random_sp_weights_err(sub,bn) = sqrt((prob_random_sp_weights(bn)*(1-prob_random_sp_weights(bn)))/length(sp_sorted_random_landing_weights(sp_indx)));
    end
end


%%
fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    if sub==1
        LH9(1) = errorbar(mn_weights_thresh(sub,:), prob_chose_in_favor_thresh_weights(sub,:),prob_chose_in_favor_thresh_weights_err(sub,:),'-ob','LineWidth',2);
        L9{1} = 'Subject';
    else
        errorbar(mn_weights_thresh(sub,:), prob_chose_in_favor_thresh_weights(sub,:),prob_chose_in_favor_thresh_weights_err(sub,:),'-ob','LineWidth',2);
    end
    hold on;
    if sub==1
        LH9(2) = errorbar(mn_weights_thresh(sub,:), prob_random_thresh_weights(sub,:),prob_random_thresh_weights_err(sub,:),'-or','LineWidth',2);
        L9{2} = 'Baseline';
        legend(LH9,L9, 'Fontsize',20, 'Box','off');
    else
        errorbar(mn_weights_thresh(sub,:), prob_random_thresh_weights(sub,:),prob_random_thresh_weights_err(sub,:),'-or','LineWidth',2);
    end
    hold on;
    if sub==3 || sub==8
        xlabel('Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t weights for threshold trials'],'fontsize',30);

%%
fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    if sub==1
        LH9(1) = errorbar(mn_sp_weights(sub,:), prob_chose_in_favor_sp_weights(sub,:),prob_chose_in_favor_sp_weights_err(sub,:),'-ob','LineWidth',2);
        L9{1} = 'Subject';
    else
        errorbar(mn_sp_weights(sub,:), prob_chose_in_favor_sp_weights(sub,:),prob_chose_in_favor_sp_weights_err(sub,:),'-ob','LineWidth',2);
    end
    hold on;
    if sub==1
        LH9(2) = errorbar(mn_sp_weights(sub,:), prob_random_sp_weights(sub,:),prob_random_sp_weights_err(sub,:),'-or','LineWidth',2);
        L9{2} = 'Baseline';
        legend(LH9,L9, 'Fontsize',20, 'Box','off');
    else
        errorbar(mn_sp_weights(sub,:), prob_random_sp_weights(sub,:),prob_random_sp_weights_err(sub,:),'-or','LineWidth',2);
    end
    hold on;
    if sub==3 || sub==8
        xlabel('Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t spatial weights for all trials'],'fontsize',30);
%%
fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    if sub==1
        LH9(1) = errorbar(mn_sp_weights_thresh(sub,:), prob_chose_in_favor_thresh_sp_weights(sub,:),prob_chose_in_favor_thresh_sp_weights_err(sub,:),'-ob','LineWidth',2);
        L9{1} = 'Subject';
    else
        errorbar(mn_sp_weights_thresh(sub,:), prob_chose_in_favor_thresh_sp_weights(sub,:),prob_chose_in_favor_thresh_sp_weights_err(sub,:),'-ob','LineWidth',2);
    end
    hold on;
    if sub==1
        LH9(2) = errorbar(mn_sp_weights_thresh(sub,:), prob_random_thresh_sp_weights(sub,:),prob_random_thresh_sp_weights_err(sub,:),'-or','LineWidth',2);
        L9{2} = 'Baseline';
        legend(LH9,L9, 'Fontsize',20, 'Box','off');
    else
        errorbar(mn_sp_weights_thresh(sub,:), prob_random_thresh_sp_weights(sub,:),prob_random_thresh_sp_weights_err(sub,:),'-or','LineWidth',2);
    end
    hold on;
    if sub==3 || sub==8
        xlabel('Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t spatial weights for threshold trials'],'fontsize',30);

%%
fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    if sub==1
        LH9(1) = errorbar(mn_weights(sub,:), prob_chose_in_favor_weights(sub,:),prob_chose_in_favor_weights_err(sub,:),'-ob','LineWidth',2);
        L9{1} = 'Subject';
    else
        errorbar(mn_weights(sub,:), prob_chose_in_favor_weights(sub,:),prob_chose_in_favor_weights_err(sub,:),'-ob','LineWidth',2);
    end
    hold on;
    if sub==1
        LH9(2) = errorbar(mn_weights(sub,:), prob_random_weights(sub,:),prob_random_weights_err(sub,:),'-or','LineWidth',2);
        L9{2} = 'Baseline';
        legend(LH9,L9, 'Fontsize',20, 'Box','off');
    else
        errorbar(mn_weights(sub,:), prob_random_weights(sub,:),prob_random_weights_err(sub,:),'-or','LineWidth',2);
    end
    hold on;
    if sub==3 || sub==8
        xlabel('Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t weights for all trials'],'fontsize',30);
%%
fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    bar(1:weight_bins, prob_chose_in_favor_thresh_weights(sub,:)-prob_random_thresh_weights(sub,:),'b');
    hold on;
    errorbar(1:weight_bins, prob_chose_in_favor_thresh_weights(sub,:)-prob_random_thresh_weights(sub,:),sqrt(prob_chose_in_favor_thresh_weights_err(sub,:).^2+prob_random_thresh_weights_err(sub,:).^2),'ok','LineWidth',2);
    hold on;
    if sub==3 || sub==8
        xlabel('Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    yline(0.0, 'k', 'Linewidth',2);
    ylim([-0.3 0.3]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t weights for threshold trials difference plot'],'fontsize',30);

%%
fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    bar(1:weight_bins, prob_chose_in_favor_sp_weights(sub,:)-prob_random_sp_weights(sub,:),'b');
    hold on;
    errorbar(1:weight_bins, prob_chose_in_favor_sp_weights(sub,:)-prob_random_sp_weights(sub,:),sqrt(prob_chose_in_favor_sp_weights_err(sub,:).^2+prob_random_sp_weights_err(sub,:).^2),'ok','LineWidth',2);
    hold on;
    if sub==3 || sub==8
        xlabel('Spatial Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    yline(0.0,'k','Linewidth',2)
    ylim([-0.25 0.25]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t spatial weights for all trials difference plot'],'fontsize',30);

%%
fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    bar(1:weight_bins, prob_chose_in_favor_weights(sub,:)-prob_random_weights(sub,:),'b');
    hold on;
    errorbar(1:weight_bins, prob_chose_in_favor_weights(sub,:)-prob_random_weights(sub,:),sqrt(prob_random_weights_err(sub,:).^2+prob_chose_in_favor_weights_err(sub,:).^2),'ok','LineWidth',2);
    hold on;
    if sub==3 || sub==8
        xlabel('Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    yline(0.0,'k','Linewidth',2);
    ylim([-0.25 0.25]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t weights for all trials difference plot'],'fontsize',30);

%%
fig=figure();
set(fig,'defaultLegendAutoUpdate','off');
for sub=1:num_sub
    subplot(2,5,sub)
    bar(1:weight_bins, prob_chose_in_favor_thresh_sp_weights(sub,:)-prob_random_thresh_sp_weights(sub,:),'b');
    hold on;
    errorbar(1:weight_bins,prob_chose_in_favor_thresh_sp_weights(sub,:)-prob_random_thresh_sp_weights(sub,:),sqrt(prob_chose_in_favor_thresh_sp_weights_err(sub,:).^2+prob_random_thresh_sp_weights_err(sub,:).^2),'ok','LineWidth',2);
    if sub==3 || sub==8
        xlabel('Weights','Fontsize',30);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose in favor'},'Fontsize',30);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    yline(0.0,'k','Linewidth',2);
    ylim([-0.4 0.4]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'Fontsize',20);
end
sgtitle(['Bias w.r.t spatial weights for threshold trials difference plot'],'fontsize',30);