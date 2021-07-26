close all; clear all; clc;
warning('off');
fix_cond = 'saccade_code';%'eyelink';
boots = 500;
fix_cluster_dist = 240; % value used to cluster fixations closeby as one point!! prevents bias from artifacts!!
hpr1 = 0.0;
hpr2 = logspace(-5, 5, 100);
boots_thresh = 1; %bootstraps to get threshold, takes too much time so set to 1

subjects = ...
    {'OvalEyeMovement-subject04';'OvalEyeMovement-subject03';...
    'OvalEyeMovement-subject07'; 'OvalEyeMovement-subject06';...
    'OvalEyeMovement-subject08';'OvalEyeMovement-subject09';...
    'OvalEyeMovement-subject10';'OvalEyeMovement-subject11';...
    'OvalEyeMovement-subject12'; 'OvalEyeMovement-subject13'};
% subjects = {'OvalEyeMovement-subject09'};
[num_sub,~] = size(subjects);
bin = [1 2 3]; % bins for evidence accumulated
bnd = [0.5 0.55];
saccade_dist_bounds = [0.0]; % saccade_dist_bounds controls which saccdes to consider setting lower bound on peripheral weights
expt_type = 1; % ratio change
peripheryPKbound = 13; % number of peripheral images to compute PK on based on distance from fixation
format short g;
max_pk = 0;
min_pk = 100000;
disp_fig = 1;

for sub=1:num_sub
    tic;
    disp(['Starting analysis for subject ' num2str(sub) ' ...']);
    [params_boot, sobl, abbl, best_hprs, data, accuracy, num_image, mean_oval_signal,...
        num_trial, mean_bin_index, prob, mean_bin_index_random, prob_random, acc_evidence, choice_in_fav, random_landing, fixations_num, eyelink, saccade_code, fixation_dist,saccades_used,index_cut{sub}]...
        = analysis_across_allsaccades_spotlight(subjects{sub}, expt_type, boots, hpr1, hpr2, bin,fix_cond,saccade_dist_bounds,fix_cluster_dist,peripheryPKbound);
    
    disp(['Doing sanity analysis for subject ' num2str(sub)  ' ...']);
    [~, ~, ~, ~, ~, ...
        num_trial_sanity, mean_bin_index_sanity, prob_sanity, mean_bin_index_random_sanity, prob_random_sanity, ~,~,~,fixations_num_sanity]...
        = analysis_across_allsaccades_sanity(subjects{sub}, expt_type, boots, params_boot, bin(1),fix_cond,fix_cluster_dist,peripheryPKbound,bnd);
    
    alpha(sub,:) = [prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
    bias(sub) = prctile(params_boot(:, end-1), 50);
    fixations_per_trial{sub} = fixations_num;
    temporal_kernel{sub} = prctile(params_boot(:,1:end-2), 50);
    max_pk = max(max_pk,max(temporal_kernel{sub}));
    min_pk = min(min_pk,min(temporal_kernel{sub}));
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
    oval_signal{sub} = mean_oval_signal;
    
    for bn=1:length(bin)
        for sc=1:length(saccade_dist_bounds)
            fixations_total(sub,sc) = saccades_used(sc);
            %             saccade_dist_per_subj(sub,sc) = saccade_dist_cases(sc);
            
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
    
    tr_ratio{sub} = data.true_ratio;
    ratio_used = data.ratio;
    choice{sub} = data.choice;
    
    disp(['Psychometric analysis for subject ' num2str(sub) ' ...']);
    ratios = linspace(0, 1, 101);%used in plots later
    [pm_fit, uniq_vals, yvals, stderrs] = GaborPsychometric(data, 1);
    subject_pm_curve_psig(sub,:) = (1-pm_fit.Fit(3)-pm_fit.Fit(4))*arrayfun(@(x) pm_fit.options.sigmoidHandle(x,pm_fit.Fit(1),pm_fit.Fit(2)), ratios)+pm_fit.Fit(4);
    
    
    disp(['All Analysis for subject ' num2str(sub) ' complete!']);
    disp('-----------------------------------------------------------------------------------------------------');
    toc;
    
    
    if disp_fig==1
        figure();
        subplot(1,3,1)
        pk = reshape(temporal_kernel{sub},peripheryPKbound,max_fixations{sub});
        imagesc(pk);set(gca,'YDir','normal'); h = colorbar; xlabel(h,'Weight','fontsize',8);h.LineWidth=1.5; h.FontSize=15;
        for pp=1:max_fixations{sub}
            for kk=peripheryPKbound:-1:1
                text(pp-0.25,kk,num2str(pk(kk,pp),'%4.2f'),'Fontsize',5);
            end
        end
        hold on;
        plot(linspace(1,max_fixations{sub},max_fixations{sub}),index_cut{sub}+0.25,'or','Linewidth',3);
        ylabel('Distance from fixation','fontsize',20);
        xlabel('Fixations','fontsize',20);
        xticks([1 fix(max_fixations{sub}/2) max_fixations{sub}]);
        yticks([1 7 13]);
        for xl=1:max_fixations{sub}
            xline(xl+0.5,'k','LineWidth',0.75);
        end
        for yl=1:13
            yline(yl+0.5,'k','LineWidth',0.75);
        end
        hold on;
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        
        subplot(1,3,2)
        llo_mn = [];
        choice_mn = [];
        err_ch = [];
        err_ch_mn = [];
        llo_mean = [];
        choice_mean = [];
        epsilon = 1e-35;
        bin_num = 20;
        %     subject_pm{sub} = (1./(1+exp(-(oval_signal{sub}*temporal_kernel{sub}'+bias(sub)))))*( 1-(alpha(sub,1)))+(alpha(sub,1)/2);
        llo = log(sigmoid(oval_signal{sub}*temporal_kernel{sub}'+bias(sub))+epsilon) - log(1 - sigmoid(oval_signal{sub}*temporal_kernel{sub}'+bias(sub))+epsilon);
        llo = llo( ~any( isnan( llo ) | isinf( llo ), 2 ),: );
        [sorted_llo,order_llo] = sort(llo);
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
        
        subplot(1,3,3)
        ratios = linspace(0, 1, 101);
        uniq_vals=linspace(-0.05,1.05,12);
        for tt=1:(length(uniq_vals)-1)
            subj_resp(sub,tt)=mean(choice{sub}(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)));
            ntrial_subj(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
        end
        uniq_vals=linspace(0,1,11);
        errorbar(uniq_vals,subj_resp(sub,:),(subj_resp(sub,:)).*(1-subj_resp(sub,:))./sqrt(ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
        subject_pm_curve = (1./(1+exp(-(oval_signal{sub}*temporal_kernel{sub}'+bias(sub)))))*( 1-(alpha(sub,1)))+(alpha(sub,1)/2);
        uniq_vals=linspace(-0.05,1.05,12);
        for tt=1:(length(uniq_vals)-1)
            subj_resp_pred(sub,tt)=mean(subject_pm_curve(((tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)))));
            ntrial_subj_pred(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
        end
        uniq_vals=linspace(0,1,11);
        hold on;
        plot(uniq_vals,subj_resp_pred(sub,:),'r','Linewidth',2);
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
        
        
%         pause;
%         close all;
    end
end
%%
disp('Plotting stuff may take a while ...............................');
for sub=1:num_sub
    fig = figure(sub);
    set(fig,'defaultLegendAutoUpdate','off');
    
    subplot(1,4,1)
    [floor(sub), thresh(sub), ~] = getThresholdWindow(data,1, 0.5, 0.7);
    %     thresh(sub,:) = getBootstrapTheshold(data,0.7,boots_thresh);
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
    
    
    subplot(1,4,2)
    pk = reshape(temporal_kernel{sub},peripheryPKbound,max_fixations{sub});
    imagesc(pk);set(gca,'YDir','normal'); h = colorbar; xlabel(h,'Weight','fontsize',8);h.LineWidth=1.5; h.FontSize=15;
    for pp=1:max_fixations{sub}
        for kk=peripheryPKbound:-1:1
            text(pp-0.25,kk,num2str(pk(kk,pp),'%4.2f'),'Fontsize',5);
        end
    end
    hold on;
    plot(linspace(1,max_fixations{sub},max_fixations{sub}),index_cut{sub}+0.25,'or','Linewidth',3);
    ylabel('Distance from fixation','fontsize',20);
    xlabel('Fixations','fontsize',20);
    xticks([1 fix(max_fixations{sub}/2) max_fixations{sub}]);
    yticks([1 7 13]);
    for xl=1:max_fixations{sub}
        xline(xl+0.5,'k','LineWidth',0.75);
    end
    for yl=1:13
        yline(yl+0.5,'k','LineWidth',0.75);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    
    
    subplot(1,4,3)
    ratios = linspace(0, 1, 101);
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp(sub,tt)=mean(choice{sub}(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)));
        ntrial_subj(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    LH2(1) = errorbar(uniq_vals,subj_resp(sub,:),(subj_resp(sub,:)).*(1-subj_resp(sub,:))./sqrt(ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
    L2{1} = 'Data';
    hold on;
    LH2(2) = plot(ratios, subject_pm_curve_psig(sub,:), 'r','LineWidth', 2);
    L2{2} = 'PS-fit curve';
    legend(LH2,L2, 'Fontsize',15, 'Box','off','Location','northwest');
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
    
    
    subplot(1,4,4)
    hold on;
    LH3(1)=bar(1, pro_choice_all_mid{1}{1}(sub),'FaceColor','b','EdgeColor','k','LineWidth',0.75);
    L3{1} = 'Subject';
    hold on;
    LH3(2)=bar(2, pro_choice_all_mid_random{1}{1}(sub),'FaceColor','r','EdgeColor','k','LineWidth',0.75);
    L3{2} = 'Baseline';
    hold on;
    errorbar([1 2],[pro_choice_all_mid{1}{1}(sub) pro_choice_all_mid_random{1}{1}(sub)],[squeeze(pro_choice_all_low{1}{1}(sub)) squeeze(pro_choice_all_low_random{1}{1}(sub))],[squeeze(pro_choice_all_high{1}{1}(sub)) squeeze(pro_choice_all_high_random{1}{1}(sub))],'ok','LineWidth',2,'linestyle','none');
    xlim([0.25 2.75]);
    text(1,0.9,['Trial Num: ' num2str(trial_num(sub))],'Fontsize',12);
    hold on;
    text(1,0.85,['Sacc Num: ' num2str(fixations_total(sub,1))],'Fontsize',12);
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
    
end
sgtitle(['All analysis for Subject ' num2str(sub)],'fontsize',30);
%%
figure();
for sub=1:num_sub
    subplot(2,5,sub)
    ratios = linspace(0, 1, 101);
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp(sub,tt)=mean(choice{sub}(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)));
        ntrial_subj(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    LH1(1) = errorbar(uniq_vals,subj_resp(sub,:),(subj_resp(sub,:)).*(1-subj_resp(sub,:))./sqrt(ntrial_subj(sub,:)),'ob','linewidth',2,'linestyle','none');
    L1{1} = 'Data';
    subject_pm_curve = (1./(1+exp(-(oval_signal{sub}*temporal_kernel{sub}'+bias(sub)))))*( 1-(alpha(sub,1)))+(alpha(sub,1)/2);
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp_pred(sub,tt)=mean(subject_pm_curve(((tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1)))));
        ntrial_subj_pred(sub,tt)=sum(tr_ratio{sub}>uniq_vals(tt)&tr_ratio{sub}<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    hold on;
    LH1(2) = plot(uniq_vals,subj_resp_pred(sub,:),'r','Linewidth',2);
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
end
sgtitle('PS fit w.r.t learned Weights/PK compared to real data','fontsize',30)
%%
bin_num=50;
figure();
llo_mn = [];
choice_mn = [];
err_ch = [];
err_ch_mn = [];
llo_mean = [];
choice_mean = [];
epsilon = 1e-100;
for sub=1:num_sub
    subplot(2,5,sub)
    llo = log(sigmoid(oval_signal{sub}*temporal_kernel{sub}'+bias(sub))+epsilon) - log(1 - sigmoid(oval_signal{sub}*temporal_kernel{sub}'+bias(sub))+epsilon);
    %     llo = llo( ~any( isnan( llo ) | isinf( llo ), 2 ),: );
    [sorted_llo,order_llo] = sort(llo);
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
for sub=1:num_sub
    nx = max_fixations{sub};
    ny = 13;
    subplot(2,5,sub)
    pk = reshape(temporal_kernel{sub},peripheryPKbound,max_fixations{sub});
    imagesc(pk);set(gca,'YDir','normal'); h = colorbar; h.FontSize=12;%xlabel(h,'Weight','fontsize',8);h.LineWidth=1.5;
    hold on;
    if sub==1 || sub==6
        ylabel('Distance from fixation','fontsize',20);
    end
    if sub==3 || sub==8
        xlabel('Fixations','fontsize',20);
    end
    xticks([1 fix(max_fixations{sub}/2) max_fixations{sub}]);
    yticks([1 7 13]);
    for xl=1:max_fixations{sub}
        xline(xl+0.5,'k','LineWidth',0.75);
    end
    for yl=1:13
        yline(yl+0.5,'k','LineWidth',0.75);
    end
    hold on;
%     caxis([min_pk, max_pk]);
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
end
sgtitle('Weights learned for subjects','fontsize',30);

%%
max_fixation_across_subj = 0;
for sub=1:num_sub
    if max_fixations{sub}>=max_fixation_across_subj
        max_fixation_across_subj = max_fixations{sub};
    end
end
pk_all_subj = zeros(peripheryPKbound,max_fixation_across_subj);
normalized_pk_all_subj = zeros(peripheryPKbound,max_fixation_across_subj);
count_fix_subj = zeros(1,max_fixation_across_subj);
for sub=1:num_sub
    pk = reshape(temporal_kernel{sub},peripheryPKbound,max_fixations{sub});
    
    % Normalise values of an array to be between -1 and 1
    % original sign of the array values is maintained.
    if abs(min(pk(:))) > max(pk(:))
        max_range_value = abs(min(pk(:)));
        min_range_value = min(pk(:));
    else
        max_range_value = max(pk(:));
        min_range_value = -max(pk(:));
    end
    norm_value = 2 .* pk(:) ./ (max_range_value - min_range_value);
    normalized_pk = reshape(norm_value,peripheryPKbound,max_fixations{sub});
    for fx=1:max_fixations{sub}
        pk_all_subj(:,fx) = pk_all_subj(:,fx) + pk(:,fx);
        normalized_pk_all_subj(:,fx) = normalized_pk_all_subj(:,fx) + normalized_pk(:,fx);
        count_fix_subj(fx) = count_fix_subj(fx) + 1;
    end
end
for fx=1:max_fixation_across_subj
    mean_pk_all_subj(:,fx) = squeeze(pk_all_subj(:,fx))/count_fix_subj(fx);
    mean_normalized_pk_all_subj(:,fx) = squeeze(normalized_pk_all_subj(:,fx))/count_fix_subj(fx);
end



figure()
subplot(1,2,1)
imagesc(mean_normalized_pk_all_subj);set(gca,'YDir','normal'); h = colorbar; h.FontSize=12;xlabel(h,'Weight','fontsize',15);
hold on;
for xl=1:max_fixation_across_subj
    if xl==8
        xline(xl+0.5,'r','LineWidth',5);
    else
        xline(xl+0.5,'k','LineWidth',0.75);
    end
end
for yl=1:peripheryPKbound
    yline(yl+0.5,'k','LineWidth',0.75);
end
hold on;
for pp=1:max_fixation_across_subj
    for kk=peripheryPKbound:-1:1
        text(pp-0.35,kk,num2str(mean_normalized_pk_all_subj(kk,pp),'%4.2f'),'Fontsize',10);
    end
end
hold on;
caxis([-1, 1]);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
xticks([1 fix(max_fixation_across_subj/2) max_fixation_across_subj]);
yticks([1 7 13]);
xlabel('Fixations','fontsize',20);
ylabel('Distance from fixation','fontsize',20);
title('Normalized between -1 and 1 before avg.','fontsize',20);

subplot(1,2,2)
imagesc(mean_pk_all_subj);set(gca,'YDir','normal'); h = colorbar; h.FontSize=12;xlabel(h,'Weight','fontsize',15);
hold on;
for xl=1:max_fixation_across_subj
    if xl==8
        xline(xl+0.5,'r','LineWidth',5);
    else
        xline(xl+0.5,'k','LineWidth',0.75);
    end
end
for yl=1:peripheryPKbound
    yline(yl+0.5,'k','LineWidth',0.75);
end
hold on;
for pp=1:max_fixation_across_subj
    for kk=peripheryPKbound:-1:1
        text(pp-0.35,kk,num2str(mean_pk_all_subj(kk,pp),'%4.2f'),'Fontsize',10);
    end
end
hold on;
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
xticks([1 fix(max_fixation_across_subj/2) max_fixation_across_subj]);
yticks([1 7 13]);
xlabel('Fixations','fontsize',20);
ylabel('Distance from fixation','fontsize',20);
title('Not normalized before avg.','fontsize',20);

sgtitle('Average weights across all subjects!','fontsize',30)

%%
closest_limit = 2;
if mod(closest_limit,2)==0
    addendum = 0;
else
    addendum = 1;
end
figure();
for nn=1:closest_limit
    subplot(2,fix(closest_limit/2) + addendum,nn)
    for sub=1:num_sub
        pk = reshape(temporal_kernel{sub},peripheryPKbound,max_fixations{sub});
        LH13 = plot(linspace(1,max_fixations{sub},max_fixations{sub}),pk(nn,:),'-ob','Linewidth',1);
        if nn>=4
            L13 = [num2str(nn) 'th closest from fixation'];
        elseif nn==1
            L13 = [num2str(nn) 'st closest from fixation'];
        elseif nn==2
            L13 = [num2str(nn) 'nd closest from fixation'];
        elseif nn==3
            L13 = [num2str(nn) 'rd closest from fixation'];
        end
        hold on;
        if nn==1 || nn==(fix(closest_limit/2) + addendum + 1)
            ylabel('Weights','fontsize',20);
        end
        if nn==(fix(closest_limit/4) + addendum) || nn==(fix(closest_limit/2) + addendum + fix(closest_limit/4) + addendum)
            xlabel('Fixations','fontsize',20);
        end
        hold on;
        yline(0.0,'k','LineWidth',0.75);
        hold on;
        legend(LH13,L13, 'Fontsize',12, 'Box','off','Location','northwest');
        ax = gca;
        ax.LineWidth=2;
        set(ax, 'box','off');
        ax.XAxis.FontSize = 15;
        ax.YAxis.FontSize = 15;
    end
end
sgtitle('Weights learned for subjects w.r.t distance from fixation','fontsize',20);
%%
% max_fixation_across_subj = 0;
% for sub=1:num_sub
%     if max_fixations{sub}>=max_fixation_across_subj
%         max_fixation_across_subj = max_fixations{sub};
%     end
% end
% pk_all_subj = zeros(peripheryPKbound,max_fixation_across_subj);
% normalized_pk_all_subj = zeros(peripheryPKbound,max_fixation_across_subj);
% count_fix_subj = zeros(1,max_fixation_across_subj);
% for sub=1:num_sub
%     pk = reshape(temporal_kernel{sub},peripheryPKbound,max_fixations{sub});
% 
%     % Normalise values of an array to be between -1 and 1
%     % original sign of the array values is maintained.
%     if abs(min(pk(:))) > max(pk(:))
%         max_range_value = abs(min(pk(:)));
%         min_range_value = min(pk(:));
%     else
%         max_range_value = max(pk(:));
%         min_range_value = -max(pk(:));
%     end
%     norm_value = 2 .* pk(:) ./ (max_range_value - min_range_value);
%     normalized_pk = reshape(norm_value,peripheryPKbound,max_fixations{sub});
%     for fx=1:max_fixations{sub}
%         pk_all_subj(:,fx) = pk_all_subj(:,fx) + pk(:,fx);
%         normalized_pk_all_subj(:,fx) = normalized_pk_all_subj(:,fx) + normalized_pk(:,fx);
%         count_fix_subj(fx) = count_fix_subj(fx) + 1;
%     end
% end
% for fx=1:max_fixation_across_subj
%     mean_pk_all_subj(:,fx) = squeeze(pk_all_subj(:,fx))/count_fix_subj(fx);
%     mean_normalized_pk_all_subj(:,fx) = squeeze(normalized_pk_all_subj(:,fx))/count_fix_subj(fx);
% end

figure();
closest_limit = 5;
subplot(2,2,1);
rng_clr = linspace(0.1,1,closest_limit);
for nn=1:closest_limit
    if nn==1
        LH14(1) = plot(linspace(1,max_fixation_across_subj,max_fixation_across_subj),mean_pk_all_subj(nn,:),'-o','color',[0 rng_clr(nn) 1],'Linewidth',2);
        L14{1} = [num2str(nn) 'st closest from fixation'];
    elseif nn==closest_limit
        if nn==2
            LH14(2) = plot(linspace(1,max_fixation_across_subj,max_fixation_across_subj),mean_pk_all_subj(nn,:),'-o','color',[0 rng_clr(nn) 1],'Linewidth',2);
            L14{2} = [num2str(nn) 'nd closest from fixation'];
        elseif nn==3
            LH14(2) = plot(linspace(1,max_fixation_across_subj,max_fixation_across_subj),mean_pk_all_subj(nn,:),'-o','color',[0 rng_clr(nn) 1],'Linewidth',2);
            L14{2} = [num2str(nn) 'rd closest from fixation'];
        else
            LH14(2) = plot(linspace(1,max_fixation_across_subj,max_fixation_across_subj),mean_pk_all_subj(nn,:),'-o','color',[0 rng_clr(nn) 1],'Linewidth',2);
            L14{2} = [num2str(nn) 'th closest from fixation'];
        end
    else
        plot(linspace(1,max_fixation_across_subj,max_fixation_across_subj),mean_pk_all_subj(nn,:),'-o','color',[0 rng_clr(nn) 1],'Linewidth',2);
    end
    hold on;
    ylabel('Weights','fontsize',20);
    xlabel('Fixations','fontsize',20);
    hold on;
    yline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
end
xline(8,'--r','LineWidth',5);
xlim([1 max_fixation_across_subj]);
title('Not normalized weights across fixations','fontsize',15);
legend(LH14,L14, 'Fontsize',12, 'Box','off','Location','northeast');

subplot(2,2,2);
fixations_chosen = 6;
rng_clr2 = linspace(0.1,1,fixations_chosen);
for nn=1:fixations_chosen
    if nn==1
        LH15(1) = plot(linspace(1,peripheryPKbound,peripheryPKbound),mean_pk_all_subj(:,nn),'-o','color',[0 rng_clr2(nn) 1],'Linewidth',2);
        L15{1} = [num2str(nn) 'st fixation'];
    elseif nn==fixations_chosen
        if nn==2
            LH15(2) = plot(linspace(1,peripheryPKbound,peripheryPKbound),mean_pk_all_subj(:,nn),'-o','color',[0 rng_clr2(nn) 1],'Linewidth',2);
            L15{2} = [num2str(nn) 'nd fixation'];
        elseif nn==3
            LH15(2) = plot(linspace(1,peripheryPKbound,peripheryPKbound),mean_pk_all_subj(:,nn),'-o','color',[0 rng_clr2(nn) 1],'Linewidth',2);
            L15{2} = [num2str(nn) 'rd fixation'];
        else
            LH15(2) = plot(linspace(1,peripheryPKbound,peripheryPKbound),mean_pk_all_subj(:,nn),'-o','color',[0 rng_clr2(nn) 1],'Linewidth',2);
            L15{2} = [num2str(nn) 'th fixation'];
        end
    else
        plot(linspace(1,peripheryPKbound,peripheryPKbound),mean_pk_all_subj(:,nn),'-o','color',[0 rng_clr2(nn) 1],'Linewidth',2);
    end
    hold on;
    ylabel('Weights','fontsize',20);
    xlabel('Order of distance from fixation','fontsize',20);
    hold on;
    yline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
end
xlim([1 peripheryPKbound]);
title('Not normalized weights across distance from fixation','fontsize',15);
legend(LH15,L15, 'Fontsize',12, 'Box','off','Location','northeast');

subplot(2,2,3);
rng_clr = linspace(0.1,1,closest_limit);
for nn=1:closest_limit
    if nn==1
        LH16(1) = plot(linspace(1,max_fixation_across_subj,max_fixation_across_subj),mean_normalized_pk_all_subj(nn,:),'-o','color',[0 rng_clr(nn) 1],'Linewidth',2);
        L16{1} = [num2str(nn) 'st closest from fixation'];
    elseif nn==closest_limit
        if nn==2
            LH16(2) = plot(linspace(1,max_fixation_across_subj,max_fixation_across_subj),mean_normalized_pk_all_subj(nn,:),'-o','color',[0 rng_clr(nn) 1],'Linewidth',2);
            L16{2} = [num2str(nn) 'nd closest from fixation'];
        elseif nn==3
            LH16(2) = plot(linspace(1,max_fixation_across_subj,max_fixation_across_subj),mean_normalized_pk_all_subj(nn,:),'-o','color',[0 rng_clr(nn) 1],'Linewidth',2);
            L16{2} = [num2str(nn) 'rd closest from fixation'];
        else
            LH16(2) = plot(linspace(1,max_fixation_across_subj,max_fixation_across_subj),mean_normalized_pk_all_subj(nn,:),'-o','color',[0 rng_clr(nn) 1],'Linewidth',2);
            L16{2} = [num2str(nn) 'th closest from fixation'];
        end
    else
        plot(linspace(1,max_fixation_across_subj,max_fixation_across_subj),mean_normalized_pk_all_subj(nn,:),'-o','color',[0 rng_clr(nn) 1],'Linewidth',2);
    end
    hold on;
    ylabel('Weights','fontsize',20);
    xlabel('Fixations','fontsize',20);
    hold on;
    yline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
end
xline(8,'--r','LineWidth',5);
xlim([1 max_fixation_across_subj]);
title('Normalized weights across fixations','fontsize',15);
legend(LH16,L16, 'Fontsize',12, 'Box','off','Location','northeast');

subplot(2,2,4);
fixations_chosen = 6;
rng_clr2 = linspace(0.1,1,fixations_chosen);
for nn=1:fixations_chosen
    if nn==1
        LH17(1) = plot(linspace(1,peripheryPKbound,peripheryPKbound),mean_normalized_pk_all_subj(:,nn),'-o','color',[0 rng_clr2(nn) 1],'Linewidth',2);
        L17{1} = [num2str(nn) 'st fixation'];
    elseif nn==fixations_chosen
        if nn==2
            LH17(2) = plot(linspace(1,peripheryPKbound,peripheryPKbound),mean_normalized_pk_all_subj(:,nn),'-o','color',[0 rng_clr2(nn) 1],'Linewidth',2);
            L17{2} = [num2str(nn) 'nd fixation'];
        elseif nn==3
            LH17(2) = plot(linspace(1,peripheryPKbound,peripheryPKbound),mean_normalized_pk_all_subj(:,nn),'-o','color',[0 rng_clr2(nn) 1],'Linewidth',2);
            L17{2} = [num2str(nn) 'rd fixation'];
        else
            LH17(2) = plot(linspace(1,peripheryPKbound,peripheryPKbound),mean_normalized_pk_all_subj(:,nn),'-o','color',[0 rng_clr2(nn) 1],'Linewidth',2);
            L17{2} = [num2str(nn) 'th fixation'];
        end
    else
        plot(linspace(1,peripheryPKbound,peripheryPKbound),mean_normalized_pk_all_subj(:,nn),'-o','color',[0 rng_clr2(nn) 1],'Linewidth',2);
    end
    hold on;
    ylabel('Weights','fontsize',20);
    xlabel('Order of distance from fixation','fontsize',20);
    hold on;
    yline(0.0,'k','LineWidth',0.75);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
end
xlim([1 peripheryPKbound]);
title('Normalized weights across distnace from fixations','fontsize',15);
legend(LH17,L17, 'Fontsize',12, 'Box','off','Location','northeast');


sgtitle('Mean weights learned for subjects w.r.t distance from fixation','fontsize',20);



%%
for i=1:length(saccade_dist_bounds)
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
        %         text(k1-0.25,0.75,['<=' num2str(round(saccade_dist_per_subj(sub,i))) ' dist'],'Fontsize',12);
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
    sgtitle(['Saccades length <= distance of periphery weights =' num2str(saccade_dist_bounds(i))],'fontsize',30);
end
%%
% for i=1:length(saccade_dist_bounds)
%
%     xlabel_txt{i} = ['sacc. <= dist. of pri wt >=' num2str(saccade_dist_bounds(i))];
%
%     prob_saccade_dist_bounds(i,:) = squeeze(pro_choice_all_mid{1}{i}(:));
%     prob_lo_err_saccade_dist_bounds(i,:) = squeeze(pro_choice_all_low{1}{i}(:));
%     prob_hi_err_saccade_dist_bounds(i,:) = squeeze(pro_choice_all_high{1}{i}(:));
%
%     prob_random_saccade_dist_bounds(i,:) = squeeze(pro_choice_all_mid_random{1}{i}(:));
%     prob_random_lo_err_saccade_dist_bounds(i,:) = squeeze(pro_choice_all_low_random{1}{i}(:));
%     prob_random_hi_err_saccade_dist_bounds(i,:) = squeeze(pro_choice_all_high_random{1}{i}(:));
% end
% figure();
% for sub=1:num_sub
%     subplot(2,5,sub)
%     if sub==1
%         LH12(1) = errorbar(1:length(saccade_dist_bounds),prob_saccade_dist_bounds(:,sub), prob_lo_err_saccade_dist_bounds(:,sub) , prob_hi_err_saccade_dist_bounds(:,sub),'-ob','LineWidth',2);
%         L12{1} = 'Subject';
%         hold on;
%         LH12(2) = errorbar(1:length(saccade_dist_bounds),prob_random_saccade_dist_bounds(:,sub), prob_random_lo_err_saccade_dist_bounds(:,sub) , prob_random_hi_err_saccade_dist_bounds(:,sub),'-or','LineWidth',2);
%         L12{2} = 'Baseline';
%     else
%         errorbar(1:length(saccade_dist_bounds),prob_saccade_dist_bounds(:,sub), prob_lo_err_saccade_dist_bounds(:,sub) , prob_hi_err_saccade_dist_bounds(:,sub),'-ob','LineWidth',2);
%         hold on;
%         errorbar(1:length(saccade_dist_bounds),prob_random_saccade_dist_bounds(:,sub), prob_random_lo_err_saccade_dist_bounds(:,sub) , prob_random_hi_err_saccade_dist_bounds(:,sub),'-or','LineWidth',2);
%     end
%     if sub==3 || sub==8
%         xlabel('Saccade length','Fontsize',20);
%     end
%     if sub==1 || sub==6
%         ylabel({'Probability of chose';' in favor'},'Fontsize',20);
%     end
%     hold on;
%     ax = gca;
%     ax.LineWidth=2;
%     set(ax, 'box','off');
%     ylim([0.4 0.75]);
%     yline(0.5,'-k','linewidth',2);
%     ax.XAxis.FontSize = 20;
%     ax.YAxis.FontSize = 20;
%     xticks(linspace(1,length(saccade_dist_bounds),length(saccade_dist_bounds)));
%     xticklabels(saccade_dist_bounds);
%     %     xticklabels(xlabel_txt);
%     if sub==1
%         legend(LH12,L12, 'Fontsize',15, 'Box','off');
%     end
%     title(['Subject ' num2str(sub)],'Fontsize',12);
% end
% sgtitle('Relation between bias and saccade length',,'fontsize',30)
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
sgtitle(['Sanity check figure with ratio > ' num2str(min(bnd(1),bnd(2))) ' and < ' num2str(max(bnd(1),bnd(2)))],'fontsize',30);

%%
for bn=2:length(bin)
    for sc=1:length(saccade_dist_bounds)
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
            sgtitle([num2str(bn) ' bins for all trials'],'fontsize',30);
        else
            sgtitle([num2str(bn) ' bins across all ratio trials (<=' num2str(saccade_dist_bounds(sc))  'x max. fixation length)'],'fontsize',30);
        end
    end
end

