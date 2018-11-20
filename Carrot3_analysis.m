clear all;
load Carrot3_data_piped;
% load data_averaged;


%% averaging RMSE of multiple iterations of fitting

source = results(1).test_set;
subject = [source.subject]';
d0 = [source.d0]';
v0 = [source.v0]';
dv = [source.dv]';
RMSE_v = [source.RMSE_v]';
colNames = {'subject', 'd0', 'v0', 'dv', 'RMSE_v'};

RMSE_table = table(subject, d0, v0, dv, RMSE_v, 'VariableNames', colNames);

if length(results) > 1
    for i = 2:length(results)
        RMSE_table.RMSE_v = RMSE_table.RMSE_v + [results(i).test_set.RMSE_v]';
    end
end
RMSE_table.RMSE_v = RMSE_table.RMSE_v / length(results);
%%
b = [repmat(1,156,1);repmat(2,156,1);repmat(3,156,1);repmat(4,156,1);repmat(5,156,1);repmat(6,156,1);repmat(7,156,1);repmat(8,156,1)];
title_ = {'model'};
b = table(b,'VariableNames',title_);
%%
save('data_RMSE_table_LOO4_RMSE_v.mat', 'c');
%% Mixed effect ANOVA on RMSE
lme = fitlme(RMSE_table,'RMSE_v ~ d0*v0*dv + (1|subject) + (1|subject:d0) + (1|subject:v0) + (1|subject:dv)',...
'DummyVarCoding','effects');
a = anova(lme);
%% Mixed effect ANOVA on RMSE
lme = fitlme(RMSE_table,'RMSE_v ~ d0*v0*dv + (1|subject)',...
'DummyVarCoding','effects');
a = anova(lme);

%% aggregate data by conditions for final state analysis

sub = [];
finalSpd = [];
finalDv = [];
finalDist = [];

count = 1;
for i = 1:length(following)
    
    subject = following(i).subject;
    trial = following(i).trial;
    d0 = following(i).d0;
    v0 = following(i).v0;
    dv = following(i).dv;
    dump = following(i).dump;
    
    
    if dump == 0
        
        sub(count,1) = subject;
        if d0 == 1
            sub(count,2) = 1;
        elseif d0 == 4
            sub(count,2) = 2;
        elseif d0 == 8
            sub(count,2) = 3;
        end

        if v0 == 0.8
            sub(count,3) = 1;
        elseif v0 == 1.2
            sub(count,3) = 2;
        end    


        if dv == -0.3
            sub(count,4) = 1;
        elseif dv == 0
            sub(count,4) = 2;
        elseif dv == 0.3
            sub(count,4) = 3;
        end
        
        finalSpd(count) = following(i).finalSpd;
        finalDv(count) = following(i).finalDv;
        finalDist(count) = following(i).finalDist;
        
        
        count = count + 1;
    end
    
end

aggregated_finalSpd = accumarray(sub(:,2:4), finalSpd, [], @mean);
aggregated_finalDv = accumarray(sub(:,2:4), finalDv, [], @mean);
aggregated_finalDist = accumarray(sub(:,2:4), finalDist, [], @mean);

%% plot ending speed
figure;
x = 1:3; % d0 = 1,4,8
hold on;
plot(x,aggregated_finalSpd(:,1,1),'b--'); % v0 = 0.8 dv = -0.3
plot(x,aggregated_finalSpd(:,2,1),'r--'); % v0 = 1.2 dv = -0.3
plot(x,aggregated_finalSpd(:,1,3),'b'); % v0 = 0.8 dv = 0.3
plot(x,aggregated_finalSpd(:,2,3),'r'); % v0 = 1.2 dv = 0.3

%% plot ending distance
figure;
x = 1:2; % d0 = 1,4,8
hold on;
plot(x,aggregated_finalSpd(1,:,1),'r--'); % d0 = 1 dv = -0.3
plot(x,aggregated_finalSpd(2,:,1),'b--'); % d0 = 4 dv = -0.3
plot(x,aggregated_finalSpd(3,:,1),'g--'); % d0 = 8 dv = -0.3
plot(x,aggregated_finalSpd(1,:,3),'r'); % d0 = 1 dv = 0.3
plot(x,aggregated_finalSpd(2,:,3),'b'); % d0 = 4 dv = 0.3
plot(x,aggregated_finalSpd(3,:,3),'g'); % d0 = 8 dv = 0.3


%% Analysis

% modelspec = 'y ~ x1 + x2 + x3 + x1*x2 + x1*x3 + x2*x3 + x1*x2*x3';
% 
% rm = fitrm(t,modelspec); % fit a repeated measures model
% 
% ranovatbl = ranova(rm,'WithinModel',WM);

%% cross-correlation
Hz = 60;
data_set = following(1);
for trial = following
    if trial.dump == 0 && trial.t_total >= 8.5 && trial.dv ~= 0 
        data_set(end+1) = trial;
    end
end
data_set = data_set(2:end);



for i = 1:length(data_set)
    t_start = int32((data_set(i).manipOnset-0.5)*Hz)+1;
    t_end = t_start + 6.0*Hz;
    for j = 1:Hz*3
        j_t_start = t_start + j;
        j_t_end = t_end - j;
        data_set(i).xCorrs(j,1) = corr(data_set(i).data(t_start:j_t_end,3), data_set(i).data(j_t_start:t_end,4));
        data_set(i).xCorrs(j,2) = j;
    end
    [M, I] = max(data_set(i).xCorrs(:,1));
    data_set(i).opt_delay = I;
    data_set(i).opt_xCorr = M;
    data_set(i).xCorr_SD = std(data_set(i).xCorrs(:,1));
    data_set(i).xCorr_range = range(data_set(i).xCorrs(:,1));

end
% save('data_xCorrs.mat', 'data_set');

%% Cross RMSE analysis for finding optimal delay
Hz = 60;
data_set = following(1);
for trial = following
    if trial.dump == 0 && trial.t_total >= 8.5 && trial.dv ~= 0
        data_set(end+1) = trial;
    end
end
data_set = data_set(2:end);



for i = 1:length(data_set)
    t_start = int32((data_set(i).manipOnset-0.5)*Hz)+1;
    t_end = t_start + 6.0*Hz;
    for j = 1:Hz*5
        j_t_start = t_start + j;
        j_t_end = t_end - j;
        data_set(i).xRMSEs(j,1) = sqrt(mean((data_set(i).data(t_start:j_t_end,3) - data_set(i).data(j_t_start:t_end,4)).^2));
        data_set(i).xRMSEs(j,2) = j;
    end
    [M, I] = min(data_set(i).xRMSEs(:,1));
    data_set(i).opt_delay = I;
    data_set(i).opt_xRMSE = M;
    data_set(i).xRMSE_SD = std(data_set(i).xRMSEs(:,1));
    data_set(i).xRMSE_range = range(data_set(i).xRMSEs(:,1));

end
% save('data_xRMSEs.mat', 'data_set');

%% Histogram of optimal delay
delay = [];
for i = 1:length(data_set)
    if data_set(i).opt_delay < 120 && data_set(i).opt_delay > 3 && data_set(i).d0 == 8 && data_set(i).dv == 0.3
        delay(end+1) = data_set(i).opt_delay/60;
    end
end

% ['Mean = ', num2str(mean(delay))]
% ['Median = ', num2str(median(delay))]
% ['Mode(0.1) = ', num2str(mode(round(delay, 1)))]
edges = -0.05:0.1:2.05;
histogram(delay,edges);
ylim([0 25]);
% histogram(delay,500);
% xticks(0:0.05:3);

%% xRMSE in each trial
path = 'C:\Users\jbai5\OneDrive\First year project\Data analysis\delay\';
for i = 1:length(data_set)
    fig = plot(data_set(i).xRMSEs(:,2),data_set(i).xRMSEs(:,1));   
    d0 = num2str(data_set(i).d0);
    v0 = num2str(data_set(i).v0);
    dv = num2str(data_set(i).dv);
    title([num2str(i), '  ', d0, ' ', v0, ' ', dv]);
    axis([0 400 0 0.3]);
%     saveas(fig, [path, 'xRMSEs', '_', num2str(i), '(d0=', d0, ' v0=', v0, ' dv=', dv, ').png']);
    pause(0.2);
    
end

%% cross correlations in each trial
path = 'C:\Users\jbai5\OneDrive\First year project\Data analysis\delay\';
for i = 1:length(data_set)
    fig = plot(data_set(i).xCorrs(:,2),data_set(i).xCorrs(:,1));   
    d0 = num2str(data_set(i).d0);
    v0 = num2str(data_set(i).v0);
    dv = num2str(data_set(i).dv);
    title([num2str(i), '  ', d0, ' ', v0, ' ', dv]);
    axis([0 400 0 1]);
%     saveas(fig, [path, 'xCorrs', '_', num2str(i), '(d0=', d0, ' v0=', v0, ' dv=', dv, ').png']);
    pause(0.2);
    
end

%% plot leader and follower on each trial
Hz = 60;
for i = 1:length(data_set)
    t_start = int32((data_set(i).manipOnset-0.5)*Hz)+1;
    t_end = t_start + 6*Hz;
    plot(1:length(data_set(i).data(t_start:t_end,3)),data_set(i).data(t_start:t_end,5));
    hold on;
    plot(1:length(data_set(i).data(t_start:t_end,4)),data_set(i).data(t_start:t_end,6));
    d0 = num2str(data_set(i).d0);
    v0 = num2str(data_set(i).v0);
    dv = num2str(data_set(i).dv);
    title([num2str(i), '  ', d0, ' ', v0, ' ', dv]);
    axis([0 400 -0.6 0.6]);
%     saveas(fig, [path, 'xRMSEs', '_', num2str(i), '(d0=', d0, ' v0=', v0, ' dv=', dv, ').png']);
    pause(0.2);
    hold off;
    
end


%% Preparing lme (BIC matrix) for BMS, init
lme = [];

%% Preparing lme (BIC matrix) for BMS
n = 361*3*2*2; % num of data point in one sim * num of trials = n (total num of data point per subject)
p = length(initParams);
lme(:,end+1) = 0;
for i = 1:13  
    sum = [];
    for j = 1:length(results)
        for k = 1:length(results(j).test_set)
            if results(j).test_set(k).subject == i
                sum = [sum, (results(j).test_set(k).RMSE_v)^2];
            end
        end
    end
    lme(i,end) = n*log(mean(sum)) + p*log(n);  % p is the num of parameters of a model
end



%% Preparing lme for BMS, save lme
save('spm_BMS_lme.mat', 'lme');

%% Bayesian Model Selection (BMS)
% lme columns are: distanceModel speedModel expansionModel sDistanceModel
% ratioModel linearModel limercierModel bruneauModel

% Nsamp = 10^6;
% do_plot = 1;
% sampling = 1;
% ecp = 1;

[alpha,exp_r,xp,pxp,bor] = spm_BMS (lme);
% Bayesian model selection for group studies
% FORMAT [alpha,exp_r,xp,pxp,bor] = spm_BMS (lme, Nsamp, do_plot, sampling, ecp, alpha0)
% 
% INPUT:
% lme      - array of log model evidences (-1/2*BIC)
%              rows: subjects
%              columns: models (1..Nk)
% Nsamp    - number of samples used to compute exceedance probabilities
%            (default: 1e6)
% do_plot  - 1 to plot p(r|y)
% sampling - use sampling to compute exact alpha
% ecp      - 1 to compute exceedance probability
% alpha0   - [1 x Nk] vector of prior model counts
% 
% OUTPUT:
% alpha   - vector of model probabilities
% exp_r   - expectation of the posterior p(r|y)
% xp      - exceedance probabilities
% pxp     - protected exceedance probabilities
% bor     - Bayes Omnibus Risk (probability that model frequencies 
%           are equal)

%% Plot average data and model by condition with confidence interval
x = 1/60:1/60:6;
x = x';
for i = 1:18
    figure;
    l = sim_ave(i).data_ave(:,2);
    f = sim_ave(i).data_ave(:,5);
    m = sim_ave(i).sim_ave(:,2);
    fCI = sim_ave(i).data_CI(:,2);
    mCI = sim_ave(i).sim_CI(:,2);
    fLow = f - fCI;
    fHigh = f + fCI;
    mLow = m - mCI;
    mHigh = m + mCI;

    % confidence intervals
    patch([x;fliplr(x')'], [fLow;fliplr(fHigh')'], 'r', 'FaceColor', [1 0.8 0.8], 'EdgeColor', 'none');
    patch([x;fliplr(x')'], [mLow;fliplr(mHigh')'], 'b', 'FaceColor', [0.8 0.8 1], 'EdgeColor', 'none');

    hold on;
    plot(x, l, 'k'); % plot the leader
    plot(x, f, 'r'); % plot the follower
    plot(x, m, 'b'); % plot the model
    axis([0 6 0 1.6]);
end
%% plot 1 iteration fitting of 1 model to all data
figure;
trial = 200;
plot(results.test_set(trial).tests(:,4));
hold on;
plot(results.test_set(trial).tests(:,3));

%% plot many trials
figure;
hold on;
for i = 1:length(results.test_set)
    if results.test_set(i).d0 == 8
        plot(results.test_set(i).tests(:,4));
%         plot(results.test_set(i).tests(:,3));
    end
end



%% Pearson's r
Pearson_r = tanh(mean([results.test_set.z_v]));
Pearson_r

%% BIC (wikipedia) form the error of 1 iteration fitting of 1 model to all data
k = length(initParams); % number of parameters
n = 143; % number of 
RMSE = mean([results.test_set.RMSE_v].^2);
BIC_wiki = n*log(RMSE) + k*log(n);
BIC_wiki
%% BIC (book) form the error of 1 iteration fitting of 1 model to all data
k = length(initParams); % number of parameters
test_length = 361;
n = test_length * n_trials;
RMSE = mean([results.test_set.RMSE_v]);
X = [];
for i = 1:length(results.test_set)
    X = [X; results.test_set(i).tests(:,3)];
end
VAR = var(X);
VARs = zeros(length(results.test_set),1);
for i = 1:length(results.test_set)
    VARs(i) = var(results.test_set(i).tests(:,3));
end
VAR_mean = mean(VARs);

BIC_book = n*log(RMSE) + k*log(n) - (n - k)*log(1 - k/n) + k*log((VAR_mean/RMSE-1)/k);
BIC_book

%% Ploting d, v at the onset of perturbation
distance = [];
velocity = [];
Hz = 60;
for i = 1:length(following)
    if following(i).dv ~= 0 && following(i).dump == 0&&following(i).d0 == 8
        manipOnset = int32(following(i).manipOnset * Hz);
        curLPos = following(i).data(manipOnset,1);
        curFPos = following(i).data(manipOnset,2);
        curDist = curLPos - curFPos;
        curFSpd = following(i).data(manipOnset,4);
        distance(end+1) = curDist;
        velocity(end+1) = curFSpd;
    end
end
figure;
histogram(distance);
figure;
histogram(velocity);
