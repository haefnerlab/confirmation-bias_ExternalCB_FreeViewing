close all; clear all; clc;
% close all;
boots = 500;
hpr1 = 0.0;
hpr2 = logspace(1, 4, 4);
subjects = {'OvalEyeMovement-subject04';'OvalEyeMovement-subject03';'OvalEyeMovement-subject07';'OvalEyeMovement-subject06';'OvalEyeMovement-subject08';'OvalEyeMovement-subject09'; 'OvalEyeMovement-subject10';'OvalEyeMovement-subject11';'OvalEyeMovement-subject12';'OvalEyeMovement-subject13'};
[num_sub,~] = size(subjects);
split_parameter = 2;
bin = 1; 
boots_thresh = 1;
choose_abs_PK = 0;
learn_num_PK = 2;
image_parameter = min(learn_num_PK,5);
format short g
for sub=1:num_sub
    disp(['Starting analysis for Subject ' num2str(sub) ' ...']);
    [params_boot,sobl,abbl,best_hprs,data,accuracy,num_image,mean_oval_signal,...
        acc_evi_all,pro_choice_all,acc_evi_ideal1,pro_choice_ideal1,...
        acc_evi_ideal2,pro_choice_ideal2,acc_evi_ideal_random2,pro_choice_ideal_random2,...
        acc_evi_ideal1avg,pro_choice_ideal1avg,...
        num_trial,tr_complete,tr_ratio,choice,fixations_num]...
        = analysis_compute_all_kernels_new(subjects{sub}, 1, boots, hpr1, hpr2, split_parameter, bin, choose_abs_PK, learn_num_PK);
    alpha(sub,:) = [prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
    bias(sub) = prctile(params_boot(:, end-1), 50);
    fixations_per_trial{sub} = fixations_num;
    trials_done(sub) = tr_complete;
    temporal_kernel(sub,:) = prctile(params_boot(:,1:end-2), 50);
    norm_temporal_kernel(sub,:) = temporal_kernel(sub,1:end-2)/mean(temporal_kernel(sub,1:end-2));
    lo_temporal_kernel(sub,:) = prctile(params_boot(:,1:end-2), 50) - prctile(params_boot(:, 1:end-2), 16);
    hi_temporal_kernel(sub,:) = prctile(params_boot(:,1:end-2), 84) - prctile(params_boot(:, 1:end-2), 50);
    beta(sub) = prctile(squeeze(abbl(:,2)),50);
    norm_all_linear(sub,:,:) = [sobl(:,1)/mean(temporal_kernel(sub,1:end-2)) sobl(:,2)/mean(temporal_kernel(sub,1:end-2)) sobl(:,3) sobl(:,4)];%sobl_norm;
    norm_slope_all(sub,:) = norm_all_linear(sub,:,1);
    norm_slope(sub) = prctile(squeeze(norm_all_linear(sub,:,1)),50);
    hprs_used(sub,:) = best_hprs;
    trial_num(sub) = num_trial;
    
    acc_evi_all_mid(sub,:) = prctile(acc_evi_all(:,:),50);
    pro_choice_all_mid(sub,:) = prctile(pro_choice_all(:,:),50);
    pro_choice_all_low(sub,:) = prctile(pro_choice_all(:,:), 50) - prctile(pro_choice_all(:,:),16);
    pro_choice_all_high(sub,:) = prctile(pro_choice_all(:,:),84) - prctile(pro_choice_all(:,:), 50);
    
    acc_evi_ideal1_mid(sub,:) = prctile(acc_evi_ideal1(:,:),50);
    pro_choice_ideal1_mid(sub,:) = prctile(pro_choice_ideal1(:,:),50);
    pro_choice_ideal1_low(sub,:) = prctile(pro_choice_ideal1(:,:), 50) - prctile(pro_choice_ideal1(:,:),16);
    pro_choice_ideal1_high(sub,:) = prctile(pro_choice_ideal1(:,:),84) - prctile(pro_choice_ideal1(:,:), 50);
    
    acc_evi_ideal2_mid(sub,:) = prctile(acc_evi_ideal2(:,:),50);
    pro_choice_ideal2_mid(sub,:) = prctile(pro_choice_ideal2(:,:),50);
    pro_choice_ideal2_low(sub,:) = prctile(pro_choice_ideal2(:,:), 50) - prctile(pro_choice_ideal2(:,:),16);
    pro_choice_ideal2_high(sub,:) = prctile(pro_choice_ideal2(:,:),84) - prctile(pro_choice_ideal2(:,:), 50);
    
    acc_evi_ideal2_random_mid(sub,:) = prctile(acc_evi_ideal_random2(:,:),50);
    pro_choice_ideal2_random_mid(sub,:) = prctile(pro_choice_ideal_random2(:,:),50);
    pro_choice_ideal2_random_low(sub,:) = prctile(pro_choice_ideal_random2(:,:), 50) - prctile(pro_choice_ideal_random2(:,:),16);
    pro_choice_ideal2_random_high(sub,:) = prctile(pro_choice_ideal_random2(:,:),84) - prctile(pro_choice_ideal_random2(:,:), 50);
    
    acc_evi_ideal1avg_mid(sub,:) = prctile(acc_evi_ideal1avg(:,:),50);
    pro_choice_ideal1avg_mid(sub,:) = prctile(pro_choice_ideal1avg(:,:),50);
    pro_choice_ideal1avg_low(sub,:) = prctile(pro_choice_ideal1avg(:,:), 50) - prctile(pro_choice_ideal1avg(:,:),16);
    pro_choice_ideal1avg_high(sub,:) = prctile(pro_choice_ideal1avg(:,:),84) - prctile(pro_choice_ideal1avg(:,:), 50);
   
    disp(['Completed analysis for Subject ' num2str(sub) '!!']);
    %%
    % for sub=1:num_sub
    figure();
    subplot(1,5,1)
    imagesc(reshape(temporal_kernel(sub,:),13,split_parameter));
    ylabel('Order of Ovals(closest to farthest w.r.t saccade)')
    xlabel('Temporal Frames')
    title('PK')
    colorbar();
    
    
    subplot(1,5,2)
    for j=1:image_parameter
        errorbar(1:split_parameter,squeeze(temporal_kernel(sub,([1:split_parameter]-1)*num_image+j)),squeeze(lo_temporal_kernel(sub,([1:split_parameter]-1)*num_image+j)),squeeze(hi_temporal_kernel(sub,([1:split_parameter]-1)*num_image+j)),'-o','LineWidth',2)
        hold on;
        Legend{j}=strcat('Oval', num2str(j));
    end
    yline(0.0,'--k','LineWidth',1);
    title(['Temporal PK for Each Frame'])
    xlim([1 split_parameter])
    %ylim([-2 2]);
    xLimits = get(gca,'XLim');  %# Get the range of the x axis
    yLimits = get(gca,'YLim');
    %     text(xLimits(1),yLimits(2)-0.07,['Trial Num: ' num2str(trial_num(sub))]);
    xticks([1:1:split_parameter]);
    xlabel('Temporal Frame');
    ylabel('PK');
    legend(Legend)
    
    
    subplot(1,5,3)
    ratios = linspace(0, 1, 101);
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp(sub,tt)=mean(choice(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)));
        ntrial_subj(sub,tt)=sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    errorbar(uniq_vals,subj_resp(sub,:),(subj_resp(sub,:)).*(1-subj_resp(sub,:))./sqrt(ntrial_subj(sub,:)),'o');
    %     [pm_fit, uniq_vals, yvals, stderrs] = GaborPsychometric(data, 1);
    %     alpha_test(1)/2 + (1-alpha_test(i,j,1))./(1+exp(-(temporal_kernel_test(i,j,k))*x-bias_test(i,j)))
    subject_pm_curve = (1./(1+exp(-(mean_oval_signal*temporal_kernel(sub,:)'+bias(sub)))))*( 1-(alpha(sub,1)))+(alpha(sub,1)/2);
    %     ((alpha(sub,1)/2) +(( 1-(alpha(sub,1)))/(1+exp(-(mean_oval_signal*temporal_kernel(sub,:)'+bias(sub))))));
    uniq_vals=linspace(-0.05,1.05,12);
    for tt=1:(length(uniq_vals)-1)
        subj_resp_pred(sub,tt)=mean(subject_pm_curve(((tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)))));
        ntrial_subj_pred(sub,tt)=sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
    end
    uniq_vals=linspace(0,1,11);
    hold on;
    plot(uniq_vals,subj_resp_pred(sub,:),'Linewidth',2);
    title(['Fitted Psychometric Curve'])
    
    
    subplot(1,5,4)
    ratios = linspace(0, 1, 101);
    avg_pm_curve = zeros(size(ratios));
    [pm_fit, uniq_vals, yvals, stderrs] = GaborPsychometric(data, 1);
    thresh(sub,:) = getBootstrapTheshold(data,0.7,boots_thresh);
    subject_pm_curve_psig(sub,:) = (1-pm_fit.Fit(3)-pm_fit.Fit(4))*arrayfun(@(x) pm_fit.options.sigmoidHandle(x,pm_fit.Fit(1),pm_fit.Fit(2)), ratios)+pm_fit.Fit(4);
    plot(ratios, subject_pm_curve_psig(sub,:), 'LineWidth', 3);
    xlabel('True ratio');
    ylabel('Percent chose left');
    axis('tight')
    hold on;
    plot(0.5*ones(1,10),linspace(0,1,10),'--k','LineWidth', 1);
    hold on;
    plot(linspace(0,1,10),0.5*ones(1,10),'--k','LineWidth', 1);
    % Plot data.
    %     errorbar(uniq_vals, yvals, stderrs, 'bs');
    title(['Psychometric Curve'])
    
    
    subplot(1,5,5)
    hold on;
    errorbar(acc_evi_all_mid(sub,:),pro_choice_all_mid(sub,:),squeeze(pro_choice_all_low(sub,:)),squeeze(pro_choice_all_high(sub,:)),'-o','LineWidth',2);
    Legend1{1} = 'All weights per frame';
    hold on;
    errorbar(acc_evi_ideal2_mid(sub,:),pro_choice_ideal2_mid(sub,:),squeeze(pro_choice_ideal2_low(sub,:)),squeeze(pro_choice_ideal2_high(sub,:)),'-o','LineWidth',2);
    Legend1{2} = 'Ideal:Closest two per frame';
    hold on;
    errorbar(acc_evi_ideal1_mid(sub,:),pro_choice_ideal1_mid(sub,:),squeeze(pro_choice_ideal1_low(sub,:)),squeeze(pro_choice_ideal1_high(sub,:)),'-o','LineWidth',2);
    Legend1{3} = 'Ideal:Closest only per frame';
    hold on;
    errorbar(acc_evi_ideal1avg_mid(sub,:),pro_choice_ideal1avg_mid(sub,:),squeeze(pro_choice_ideal1avg_low(sub,:)),squeeze(pro_choice_ideal1avg_high(sub,:)),'-o','LineWidth',2);
    Legend1{4} = 'Ideal:Average of closest across frames';
    hold on;
    title('Choice vs. Evidence')
    ylim([0.0 1]);
    text(xLimits(1),yLimits(2)-0.02,['Trial Num: ' num2str(trial_num(sub))])
    yline(pro_choice_ideal2_random_mid(sub),'--k','LineWidth',1);
    xlabel('Accumulated Evidence');
    ylabel('Probability of chose in favor of accumulated evidence');
    legend(Legend1);
    
    suptitle(['Summary plot for subject ' num2str(sub) ' with ' num2str(trial_num(sub)) ' trials'])

end
%%
figure()
axis image;
subplot(2,3,3)
j=2;
for sub=1:num_sub
    if sub==4
        plot(1:split_parameter,squeeze(temporal_kernel(sub,([1:split_parameter]-1)*num_image+j)),'-dm','LineWidth',0.2);
    else
        plot(1:split_parameter,squeeze(temporal_kernel(sub,([1:split_parameter]-1)*num_image+j)),'-om','LineWidth',0.2);
    end
    hold on;
end
LH(1) = scatter(1:split_parameter,mean(squeeze(temporal_kernel(:,([1:split_parameter]-1)*num_image+j)),1),100,'m','filled');
L{1} = '2nd closest to fixation';
hold on;
plot(1:split_parameter,mean(squeeze(temporal_kernel(:,([1:split_parameter]-1)*num_image+j)),1),'-m','LineWidth',2);
% yline(0.0,'--k','LineWidth',1);
xlim([1-0.1 split_parameter+0.1])
ylim([0.0 3.5]);
yLimits = get(gca,'YLim');
xticks([1:1:split_parameter]);
set(gca,'Fontsize',20)
xlabel('Temporal division of trial','Fontsize',20);
ylabel({'Regression weight to'; 'predict final decision'},'Fontsize',20);
% ylabel('Weights for signal 2^{nd} closest to saccade landing','Fontsize',10');
hold on;
j=1;
for sub=1:num_sub
    if sub==4
        plot(1:split_parameter,squeeze(temporal_kernel(sub,([1:split_parameter]-1)*num_image+j)),'-dk','LineWidth',0.2);
    else
        plot(1:split_parameter,squeeze(temporal_kernel(sub,([1:split_parameter]-1)*num_image+j)),'-ok','LineWidth',0.2);
    end
    hold on;
end
LH(2) = scatter(1:split_parameter,mean(squeeze(temporal_kernel(:,([1:split_parameter]-1)*num_image+j)),1),100,'k','filled');
L{2} = 'Closest to fixation';
hold on;
plot(1:split_parameter,mean(squeeze(temporal_kernel(:,([1:split_parameter]-1)*num_image+j)),1),'k','LineWidth',2);
% yline(0.0,'k','LineWidth',1);
xlim([1-0.1 split_parameter+0.1])
xticks([1 2])
ylim([0.0 6]);
yLimits = get(gca,'YLim');
xticks([1:1:split_parameter]);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
% xlabel('Temporal Frames','Fontsize',10);
% ylabel('Weights for signal closest to saccade landing','Fontsize',10);
legend(LH,L, 'Fontsize',20, 'Box','off');
hold on;


subplot(2,3,1)
ratios = linspace(0, 1, 101);
for sub=1:num_sub
    if sub==4
        plot(ratios, subject_pm_curve_psig(sub,:), '--k','LineWidth', 0.25);
    else
        plot(ratios, subject_pm_curve_psig(sub,:), 'k','LineWidth', 0.25);
    end
    hold on;
end
plot(ratios, mean(subject_pm_curve_psig(:,:),1), 'k','LineWidth', 3);
xlabel('True ratio');
ylabel('Percent chose left');
axis('tight')
hold on;
plot(0.5*ones(1,10),linspace(0,1,10),'k','LineWidth', 1.5);
hold on;
plot(linspace(0,1,10),0.5*ones(1,10),'k','LineWidth', 1.5);
set(gca,'Fontsize',20)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ylabel('Probaility of horizontal choice','Fontsize',20);
xlabel('Ratio of horizontal to vertical','Fontsize',20)


subplot(2,3,2)
bar(1:num_sub,pro_choice_all_mid,'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
hold on;
bar(4,pro_choice_all_mid(4),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
hold on;
yline(0.5,'k','LineWidth',1.5);
hold on;
errorbar(1:num_sub,pro_choice_all_mid(:),squeeze(pro_choice_all_low(:)),squeeze(pro_choice_all_high(:)),'ok','LineWidth',2);
hold on;
xlabel('Subject Number','Fontsize',20);
ylabel({'Probability of saccading to ';'evidence confirming current belief'},'Fontsize',20);
yline(0.5,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ylim([0.4 1]);
xlim([0.5 num_sub+0.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = linspace(1,num_sub,num_sub);
hold on;

%%
figure()
subplot(1,2,1)
bar(1:num_sub,pro_choice_ideal2_mid,'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
hold on;
bar(4,pro_choice_ideal2_mid(4),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
hold on;
yline(0.5,'k','LineWidth',1.5);
hold on;
errorbar(1:num_sub,pro_choice_ideal2_mid(:),squeeze(pro_choice_ideal2_low(:)),squeeze(pro_choice_ideal2_high(:)),'ok','LineWidth',2);
hold on;
xlabel('Subject Number','Fontsize',20);
ylabel({'Probability of saccading to ';'evidence confirming current belief'},'Fontsize',20);
yline(0.5,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ylim([0.0 1]);
xlim([0.5 num_sub+0.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = linspace(1,num_sub,num_sub);
hold on;


subplot(1,2,2)
bar(1:num_sub,pro_choice_ideal2_random_mid,'FaceColor','b','EdgeColor','k','LineWidth',0.75);
hold on;
bar(4,pro_choice_ideal2_random_mid(4),'FaceColor','b','EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
hold on;
yline(0.5,'k','LineWidth',1.5);
hold on;
errorbar(1:num_sub,pro_choice_ideal2_random_mid(:),squeeze(pro_choice_ideal2_random_low(:)),squeeze(pro_choice_ideal2_random_high(:)),'ok','LineWidth',2);
hold on;
xlabel('Subject Number','Fontsize',20);
ylabel({'Baseline probability'},'Fontsize',20);
yline(0.5,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ylim([0.0 1]);
xlim([0.5 num_sub+0.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = linspace(1,num_sub,num_sub);
hold on;


%%
figure();
for i=1:num_sub
    model_series(i,:) = [pro_choice_ideal2_mid(i) pro_choice_ideal2_random_mid(i)];
    model_error_low(i,:) = [pro_choice_ideal2_low(i) pro_choice_ideal2_random_low(i)];
    model_error_high(i,:) = [pro_choice_ideal2_high(i) pro_choice_ideal2_random_high(i)];
end
b = bar(model_series, 'grouped');
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error_low(:,i), model_error_high(:,i),'ok', 'linestyle', 'none');
end
title('Comparison: Blue: Subject, Red: Baseline','fontsize',30);