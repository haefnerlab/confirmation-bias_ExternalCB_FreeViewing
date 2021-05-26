close all; clear all; clc;
warning('off');
fix_cond = 'saccade_code';%'eyelink';
boots = 500;
fix_cluster_dist = 240; % value used to cluster fixations closeby as one point!! prevents bias from artifacts!!
hpr1 = 0.0;
hpr2 = logspace(1, 5, 5);  %logspace(-1, 5, 7);
boots_thresh = 1; %bootstraps to get threshold, takes too much time so set to 1
subjects = {'OvalEyeMovement-subject04';'OvalEyeMovement-subject03';'OvalEyeMovement-subject07';'OvalEyeMovement-subject06';'OvalEyeMovement-subject08';'OvalEyeMovement-subject09'; 'OvalEyeMovement-subject10';'OvalEyeMovement-subject11';'OvalEyeMovement-subject12';'OvalEyeMovement-subject13'};
[num_sub,~] = size(subjects);
bin = [1 2 3]; % bins for evidence accumulated
saccade_dist_types = [1 0.5 0.25 0.17]; % saccade_dist_types x max saccade length used
expt_type = 1; % ratio change
format short g;

for sub=1:num_sub
    disp(['Starting analysis for subject ' num2str(sub) ' ...']);
    [params_boot,sobl,abbl,best_hprs,data,accuracy,num_image,mean_oval_signal,...
        num_trial,mean_bin_index, prob, mean_bin_index_random, prob_random, acc_evidence,choice_in_fav,random_landing,fixations_num,fixations_num_cases,eyelink,saccade_code,saccade_dist_cases,fixation_dist]...
        = analysis_across_allsaccades(subjects{sub}, expt_type, boots, hpr1, hpr2, bin,fix_cond,saccade_dist_types,fix_cluster_dist);
    
    disp(['Doing sanity analysis for subject ' num2str(sub)  ' ...']);
    [~,~,~,~,~,...
        num_trial_sanity,mean_bin_index_sanity, prob_sanity, mean_bin_index_random_sanity, prob_random_sanity, ~,~,~,fixations_num_sanity]...
        = analysis_across_allsaccades_sanity(subjects{sub}, expt_type, boots, params_boot, bin(1),fix_cond,fix_cluster_dist);
    
    alpha(sub,:) = [prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
    bias(sub) = prctile(params_boot(:, end-1), 50);
    fixations_per_trial{sub} = fixations_num;
    temporal_kernel{sub} = prctile(params_boot(:,1:end-2), 50);
    lo_temporal_kernel{sub} = prctile(params_boot(:,1:end-2), 50) - prctile(params_boot(:, 1:end-2), 16);
    hi_temporal_kernel{sub} = prctile(params_boot(:,1:end-2), 84) - prctile(params_boot(:, 1:end-2), 50);
    beta(sub) = prctile(squeeze(abbl(:,2)),50);
    norm_all_linear(sub,:,:) = [sobl(:,1)/mean(temporal_kernel{sub}(1:end-2)) sobl(:,2)/mean(temporal_kernel{sub}(1:end-2)) sobl(:,3) sobl(:,4)];%sobl_norm;
    norm_slope_all(sub,:) = norm_all_linear(sub,:,1);
    norm_slope(sub) = prctile(squeeze(norm_all_linear(sub,:,1)),50);
    hprs_used(sub,:) = best_hprs;
    trial_num(sub) = num_trial;
    trial_num_sanity(sub) = num_trial_sanity;
    fixations{sub} = fixations_num;
    fixations_sanity{sub} = fixations_num_sanity;
    max_fixations{sub} = max(fixations_num);
    eyelink_per_subj{sub} = eyelink;
    saccade_code_per_subj{sub} = saccade_code;
    fixation_dist_per_subj{sub} = fixation_dist;
    
    for bn=1:length(bin)
        for sc=1:length(saccade_dist_types)
            fixations_total(sub,sc) = fixations_num_cases(sc);
            saccade_dist_per_subj(sub,sc) = saccade_dist_cases(sc);
            
            acc_evi_all_mid{bn}{sc}(sub,:) = prctile(mean_bin_index{bn}{sc}(:,:),50);
            pro_choice_all_mid{bn}{sc}(sub,:) = prctile(prob{bn}{sc}(:,:),50);
            pro_choice_all_low{bn}{sc}(sub,:) = prctile(prob{bn}{sc}(:,:), 50) - prctile(prob{bn}{sc}(:,:),16);
            pro_choice_all_high{bn}{sc}(sub,:) = prctile(prob{bn}{sc}(:,:),84) - prctile(prob{bn}{sc}(:,:), 50);
            
            acc_evi_all_mid_random{bn}{sc}(sub,:) = prctile(mean_bin_index_random{bn}{sc}(:,:),50);
            pro_choice_all_mid_random{bn}{sc}(sub,:) = prctile(prob_random{bn}{sc}(:,:),50);
            pro_choice_all_low_random{bn}{sc}(sub,:) = prctile(prob_random{bn}{sc}(:,:), 50) - prctile(prob_random{bn}{sc}(:,:),16);
            pro_choice_all_high_random{bn}{sc}(sub,:) = prctile(prob_random{bn}{sc}(:,:),84) - prctile(prob_random{bn}{sc}(:,:), 50);
            
        end
        
    end
    acc_evi_all_mid_sanity(sub,:) = prctile(mean_bin_index_sanity(:,:),50);
    pro_choice_all_mid_sanity(sub,:) = prctile(prob_sanity(:,:),50);
    pro_choice_all_low_sanity(sub,:) = prctile(prob_sanity(:,:), 50) - prctile(prob_sanity(:,:),16);
    pro_choice_all_high_sanity(sub,:) = prctile(prob_sanity(:,:),84) - prctile(prob_sanity(:,:), 50);
    
    acc_evi_all_mid_random_sanity(sub,:) = prctile(mean_bin_index_random_sanity(:,:),50);
    pro_choice_all_mid_random_sanity(sub,:) = prctile(prob_random_sanity(:,:),50);
    pro_choice_all_low_random_sanity(sub,:) = prctile(prob_random_sanity(:,:), 50) - prctile(prob_random_sanity(:,:),16);
    pro_choice_all_high_random_sanity(sub,:) = prctile(prob_random_sanity(:,:),84) - prctile(prob_random_sanity(:,:), 50);
    
    tr_ratio = data.true_ratio;
    ratio_used = data.ratio;
    choice = data.choice;
    
    disp(['Analysis for subject ' num2str(sub) ' complete!']);
    disp('-----------------------------------------------------------------------------------------------------');

    
    %%
    fig = figure(sub);
    set(fig,'defaultLegendAutoUpdate','off');
    
    subplot(1,5,1)
    [floor(sub), thresh(sub), ~] = getThresholdWindow(data,1, 0.5, 0.7);
    %     thresh(sub,:) = getBootstrapTheshold(data,0.7,boots_thresh);
    plot(ratio_used,'-ob');
    xlabel('Trial number','Fontsize',20);
    ylabel({'Ratio of';' vertical ellipses'},'Fontsize',20);
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
    
    
    subplot(1,5,2)
    LH(1) = errorbar(1:max_fixations{sub},squeeze(temporal_kernel{sub}(1:max_fixations{sub})),squeeze(lo_temporal_kernel{sub}(1:max_fixations{sub})),squeeze(hi_temporal_kernel{sub}(1:max_fixations{sub})),'-ob','LineWidth',2);
    L{1} = 'Closest fix.';
    hold on;
    LH(2) = errorbar(1:max_fixations{sub},squeeze(temporal_kernel{sub}(max_fixations{sub}+1:end)),squeeze(lo_temporal_kernel{sub}(max_fixations{sub}+1:end)),squeeze(hi_temporal_kernel{sub}(max_fixations{sub}+1:end)),'-or','LineWidth',2);
    L{2} = '2nd closest fix.';
    hold on;
    yline(0.0,'-k','LineWidth',2);
    ylabel('Weight on each fixation','fontsize',20);
    xlim([1 max_fixations{sub}]);
    xlabel('Fixations','fontsize',20);
    legend(LH,L, 'Fontsize',10, 'Box','off');
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    
    
    subplot(1,5,3)
    ratios = linspace(0, 1, 101);
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp(sub,tt)=mean(choice(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)));
        ntrial_subj(sub,tt)=sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    LH1(1) = errorbar(uniq_vals,subj_resp(sub,:),(subj_resp(sub,:)).*(1-subj_resp(sub,:))./sqrt(ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
    L1{1} = 'Data';
    subject_pm_curve = (1./(1+exp(-(mean_oval_signal*temporal_kernel{sub}'+bias(sub)))))*( 1-(alpha(sub,1)))+(alpha(sub,1)/2);
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp_pred(sub,tt)=mean(subject_pm_curve(((tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)))));
        ntrial_subj_pred(sub,tt)=sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    hold on;
    LH1(2) = plot(uniq_vals,subj_resp_pred(sub,:),'r','Linewidth',2);
    L1{2} = 'Fit with weights';
    ylabel('Probability chose vertical','Fontsize',20);
    xlabel({'Ratio of';' vertical ellipses'},'Fontsize',20);
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
    legend(LH1,L1, 'Fontsize',10, 'Box','off');
    
    
    subplot(1,5,4)
    LH2(1) = errorbar(uniq_vals,subj_resp(sub,:),(subj_resp(sub,:)).*(1-subj_resp(sub,:))./sqrt(ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
    L2{1} = 'Data';
    hold on;
    ratios = linspace(0, 1, 101);
    avg_pm_curve = zeros(size(ratios));
    [pm_fit, uniq_vals, yvals, stderrs] = GaborPsychometric(data, 1);
    subject_pm_curve_psig(sub,:) = (1-pm_fit.Fit(3)-pm_fit.Fit(4))*arrayfun(@(x) pm_fit.options.sigmoidHandle(x,pm_fit.Fit(1),pm_fit.Fit(2)), ratios)+pm_fit.Fit(4);
    LH2(2) = plot(ratios, subject_pm_curve_psig(sub,:), 'r','LineWidth', 2);
    L2{2} = 'PS-fit curve';
    legend(LH2,L2, 'Fontsize',10, 'Box','off');
    hold on;
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.5,'-k','linewidth',2);
    hold on;
    ylabel('Probability chose vertical','Fontsize',20);
    xlabel({'Ratio of';' vertical ellipses'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    
    
    subplot(1,5,5)
    hold on;
    LH3(1)=bar(1, pro_choice_all_mid{1}{1}(sub),'FaceColor','b','EdgeColor','k','LineWidth',0.75);
    L3{1} = 'Subject';
    hold on;
    LH3(2)=bar(2, pro_choice_all_mid_random{1}{1}(sub),'FaceColor','r','EdgeColor','k','LineWidth',0.75);
    L3{2} = 'Random';
    hold on;
    errorbar([1 2],[pro_choice_all_mid{1}{1}(sub) pro_choice_all_mid_random{1}{1}(sub)],[squeeze(pro_choice_all_low{1}{1}(sub)) squeeze(pro_choice_all_low_random{1}{1}(sub))],[squeeze(pro_choice_all_high{1}{1}(sub)) squeeze(pro_choice_all_high_random{1}{1}(sub))],'ok','LineWidth',2,'linestyle','none');
    xlim([0.25 2.75]);
    text(1,0.9,['Trial Num: ' num2str(trial_num(sub))],'Fontsize',12);
    hold on;
    text(1,0.8,['Sacc Num: ' num2str(fixations_total(sub,1))],'Fontsize',12);
    xlabel('Accumulated Evidence','Fontsize',20);
    ylabel('Probability of chose in favor of accumulated evidence','Fontsize',20);
    hold on;
    legend(LH3,L3, 'Fontsize',10, 'Box','off');
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    suptitle(['All analysis for subject ' num2str(sub)]);
end
%%

figure();
subplot(2,2,1)
for sub=1:num_sub
    ratios = linspace(0, 1, 101);
    plot(ratios, subject_pm_curve_psig(sub,:), 'color',[0 .5 .5],'LineWidth', 0.75);
    hold on;
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
end
plot(ratios, mean(subject_pm_curve_psig,1), 'k','LineWidth', 2);

subplot(2,2,2)
for sub=1:num_sub
    if sub==num_sub
        LH7(1)=plot(1:max_fixations{sub},squeeze(temporal_kernel{sub}(1:max_fixations{sub})),'-ob','Linewidth',0.75);
        L7{1} = 'Closest to fixation';
        hold on;
    else
        plot(1:max_fixations{sub},squeeze(temporal_kernel{sub}(1:max_fixations{sub})),'-ob','Linewidth',0.75);
        hold on;
    end
    if sub==num_sub
        LH7(2)=plot(1:max_fixations{sub},squeeze(temporal_kernel{sub}(max_fixations{sub}+1:end)),'-.r','Linewidth',0.75);
        L7{2} = '2nd Closest to fixation';
        hold on;
    else
        plot(1:max_fixations{sub},squeeze(temporal_kernel{sub}(max_fixations{sub}+1:end)),'-.r','Linewidth',0.75);
        hold on;
    end
    yline(0.0,'-k','LineWidth',2);
    ylabel('Weight on each fixation','fontsize',20);
    %     ylim([-1 1]);
    xlabel('Fixations','fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end
legend(LH7,L7, 'Fontsize',20, 'Box','off','Location','northwest');
xlim([1 max_fixations{sub}]);

subplot(2,1,2)
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, pro_choice_all_mid{1}{1}(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
        L6{1} = 'Subject';
    else
        bar(k1, pro_choice_all_mid{1}{1}(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, pro_choice_all_mid_random{1}{1}(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
        L6{2} = 'Baseline';
    else
        bar(k2, pro_choice_all_mid_random{1}{1}(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    errorbar([k1 k2],[pro_choice_all_mid{1}{1}(sub) pro_choice_all_mid_random{1}{1}(sub)],[squeeze(pro_choice_all_low{1}{1}(sub)) squeeze(pro_choice_all_low_random{1}{1}(sub))],[squeeze(pro_choice_all_high{1}{1}(sub)) squeeze(pro_choice_all_high_random{1}{1}(sub))],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.25,0.85,[num2str(trial_num(sub)) ' trials'],'Fontsize',15);
    text(k1-0.25,0.8,[num2str(fixations_total(sub,1)) ' sac.'],'Fontsize',15);
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
xticks(k+0.5)
legend(LH6,L6, 'Fontsize',20, 'Box','off');


%%
for i=1:length(saccade_dist_types)
    figure();
    k1 = 1;
    k2 = 2;
    k = [];
    for sub=1:num_sub
        if sub==num_sub
            LH6(1) = bar(k1, pro_choice_all_mid{1}{i}(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
            L6{1} = 'Subject';
        else
            bar(k1, pro_choice_all_mid{1}{i}(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
        end
        hold on;
        if sub==num_sub
            LH6(2) = bar(k2, pro_choice_all_mid_random{1}{i}(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
            L6{2} = 'Baseline';
        else
            bar(k2, pro_choice_all_mid_random{1}{i}(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
        end
        hold on;
        errorbar([k1 k2],[pro_choice_all_mid{1}{i}(sub) pro_choice_all_mid_random{1}{i}(sub)],[squeeze(pro_choice_all_low{1}{i}(sub)) squeeze(pro_choice_all_low_random{1}{i}(sub))],[squeeze(pro_choice_all_high{1}{i}(sub)) squeeze(pro_choice_all_high_random{1}{i}(sub))],'ok','LineWidth',2,'linestyle','none');
        text(k1-0.25,0.9,[num2str(trial_num(sub)) ' trials'],'Fontsize',15);
        text(k1-0.25,0.85,[num2str(round((fixations_total(sub,i)/trial_num(sub)))) ' avg sac.'],'Fontsize',15);
        text(k1-0.25,0.8,[num2str(fixations_total(sub,i)) ' sac.'],'Fontsize',15);
        text(k1-0.25,0.75,['<=' num2str(round(saccade_dist_per_subj(sub,i))) ' dist'],'Fontsize',12);
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
    if saccade_dist_types(i)==1
        suptitle('All saccades of each subject across all ratio trials');
    else
        suptitle(['Short saccades of each subject across all ratio trials (<=' num2str(saccade_dist_types(i))  ' x max. fixation length)']);
    end
end
%%
for i=1:length(saccade_dist_types)
    if saccade_dist_types(i)==1
        xlabel_txt{i} = 'All sacc';
    else
        xlabel_txt{i} = ['<=' num2str(saccade_dist_types(i))  ' x max. sacc len'];
    end
    prob_saccade_dist_types(i,:) = squeeze(pro_choice_all_mid{1}{i}(:));
    prob_lo_err_saccade_dist_types(i,:) = squeeze(pro_choice_all_low{1}{i}(:));
    prob_hi_err_saccade_dist_types(i,:) = squeeze(pro_choice_all_high{1}{i}(:));
    
    prob_random_saccade_dist_types(i,:) = squeeze(pro_choice_all_mid_random{1}{i}(:));
    prob_random_lo_err_saccade_dist_types(i,:) = squeeze(pro_choice_all_low_random{1}{i}(:));
    prob_random_hi_err_saccade_dist_types(i,:) = squeeze(pro_choice_all_high_random{1}{i}(:));
end
figure();
for sub=1:num_sub
    subplot(2,5,sub)
    if sub==1
        LH12(1) = errorbar(1:length(saccade_dist_types),prob_saccade_dist_types(:,sub), prob_lo_err_saccade_dist_types(:,sub) , prob_hi_err_saccade_dist_types(:,sub),'-ob','LineWidth',2);
        L12{1} = 'Subject';
        hold on;
        LH12(2) = errorbar(1:length(saccade_dist_types),prob_random_saccade_dist_types(:,sub), prob_random_lo_err_saccade_dist_types(:,sub) , prob_random_hi_err_saccade_dist_types(:,sub),'-or','LineWidth',2);
        L12{2} = 'Baseline';
    else
        errorbar(1:length(saccade_dist_types),prob_saccade_dist_types(:,sub), prob_lo_err_saccade_dist_types(:,sub) , prob_hi_err_saccade_dist_types(:,sub),'-ob','LineWidth',2);
        hold on;
        errorbar(1:length(saccade_dist_types),prob_random_saccade_dist_types(:,sub), prob_random_lo_err_saccade_dist_types(:,sub) , prob_random_hi_err_saccade_dist_types(:,sub),'-or','LineWidth',2);
    end
    if sub==3 || sub==8
        xlabel('Value x Max. saccade length','Fontsize',20);
    end
    if sub==1 || sub==6
        ylabel({'Probability of chose';' in favor'},'Fontsize',20);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 0.75]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    xticks(linspace(1,length(saccade_dist_types),length(saccade_dist_types)));
    xticklabels(saccade_dist_types);
%     xticklabels(xlabel_txt);
    if sub==1
        legend(LH12,L12, 'Fontsize',15, 'Box','off');
    end
    title(['Subject ' num2str(sub)],'Fontsize',15)
end
%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, pro_choice_all_mid_sanity(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
        L6{1} = 'Subject';
    else
        bar(k1, pro_choice_all_mid_sanity(sub),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, pro_choice_all_mid_random_sanity(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
        L6{2} = 'Baseline';
    else
        bar(k2, pro_choice_all_mid_random_sanity(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    errorbar([k1 k2],[pro_choice_all_mid_sanity(sub) pro_choice_all_mid_random_sanity(sub)],[squeeze(pro_choice_all_low_sanity(sub)) squeeze(pro_choice_all_low_random_sanity(sub))],[squeeze(pro_choice_all_high_sanity(sub)) squeeze(pro_choice_all_high_random_sanity(sub))],'ok','LineWidth',2,'linestyle','none');
    text(k1-0.25,0.9,[num2str(trial_num_sanity(sub)) ' trials'],'Fontsize',15)
    text(k1-0.25,0.85,[num2str(round(mean(fixations_sanity{sub}))) ' avg sac.'],'Fontsize',15);
    text(k1-0.25,0.8,[num2str(sum(fixations_sanity{sub})-trial_num_sanity(sub)) ' sac.'],'Fontsize',15);
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
suptitle('Sanity check figure with ratio>0.5 and<0.6');

%%
for bn=2:length(bin)
    for sc=1:length(saccade_dist_types)
        fig = figure();
        set(fig,'defaultLegendAutoUpdate','off');
        for sub=1:num_sub
            subplot(2,5,sub)
            if sub==1
                LH9(1) = errorbar(acc_evi_all_mid{bn}{sc}(sub,:),pro_choice_all_mid{bn}{sc}(sub,:),pro_choice_all_low{bn}{sc}(sub,:), pro_choice_all_high{bn}{sc}(sub,:),'-ob','LineWidth',2);
                L9{1} = 'Subject';
            else
                errorbar(acc_evi_all_mid{bn}{sc}(sub,:),pro_choice_all_mid{bn}{sc}(sub,:) ,pro_choice_all_low{bn}{sc}(sub,:),pro_choice_all_high{bn}{sc}(sub,:),'-ob','LineWidth',2);
            end
            hold on;
            if sub==1
                LH9(2) = errorbar(acc_evi_all_mid_random{bn}{sc}(sub,:),pro_choice_all_mid_random{bn}{sc}(sub,:),pro_choice_all_low_random{bn}{sc}(sub,:),pro_choice_all_high_random{bn}{sc}(sub,:),'-or','LineWidth',2);
                L9{2} = 'Baseline';
                legend(LH9,L9, 'Fontsize',20, 'Box','off');
            else
                errorbar(acc_evi_all_mid_random{bn}{sc}(sub,:),pro_choice_all_mid_random{bn}{sc}(sub,:),pro_choice_all_low_random{bn}{sc}(sub,:),pro_choice_all_high_random{bn}{sc}(sub,:),'-or','LineWidth',2);
            end
            hold on;
            text(acc_evi_all_mid_random{bn}{sc}(sub,1)+0.1,0.9,[num2str(trial_num(sub)) ' trials'],'Fontsize',15);
            if sub==3 || sub==8
                xlabel('Accumulated evidence','Fontsize',30);
            end
            if sub==1 || sub==6
                ylabel({'Probability of chose in favor';' of accumulated evidence'},'Fontsize',30);
            end
            hold on;
            ax = gca;
            ax.LineWidth=2;
            set(ax, 'box','off');
            ylim([0.4 1]);
            xlim([acc_evi_all_mid_random{bn}{sc}(sub,1)-1 acc_evi_all_mid_random{bn}{sc}(sub,end)+1]);
            yline(0.5,'-k','linewidth',2);
            ax.XAxis.FontSize = 20;
            ax.YAxis.FontSize = 20;
            title(['Subject ' num2str(sub)],'Fontsize',20);
        end
        if sc==1
            suptitle([num2str(bn) ' bins for all trials']);
        else
            suptitle([num2str(bn) ' bins across all ratio trials (<=' num2str(saccade_dist_types(sc))  'x max. fixation length)']);
        end
    end
end