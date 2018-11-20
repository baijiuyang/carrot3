% This m file test a model with a fixed parameter(s) and gives Pearson_r,
% Fisher_z, and RMSE of each trial. It also gives the average simulation of
% each condition, including non-perturbation trials.

% The length of simulation for testing is 6 seconds ( from 0.5 s before
% perturbation to 5.5 s after).

% 1/24/2018 Jiuyang Bai baijiuyang@hotmail.com

%% load data_set for testing

load data_piped;
model = 'fit_expansionModel'; % 'fit_distanceModel'
                              % 'fit_speedModel'
                              % 'fit_expansionModel'
                              % 'fit_sDistanceModel'
                              % 'fit_ratioModel'
                              % 'fit_linearModel'
                              % 'fit_lemercierModel'
                              % 'fit_bruneauModel'
p = [3];
delay = 0;
Hz = 60;
time = 6; % time length of simulation
min_length = 8.5; % the smallest acceptable trial length
t_critical = 2; % df = 60-64, confidence level = 95%
%% select trials for simulation
simulation = following(1);
for i = 1:length(following)
    if following(i).dump == 0 && following(i).t_total >= min_length %&& following(i).d0 ~= 8 
        simulation(end+1) = following(i);       
    end
end
simulation = simulation(2:end);
%% set up output structure
sim_ave = struct;
for d0 = [1 4 8]
    for v0 = [0.8 1.2]
        for dv = [-0.3 0 0.3]
            sim_ave(end+1).d0 = d0;
            sim_ave(end).v0 = v0;
            sim_ave(end).dv = dv;
            sim_ave(end).data_ave = zeros(time*Hz,6);
            sim_ave(end).sim_ave = zeros(time*Hz,3);
            sim_ave(end).lPos = [];
            sim_ave(end).lSpd = [];
            sim_ave(end).lAcc = [];
            sim_ave(end).fPos = [];
            sim_ave(end).fSpd = [];
            sim_ave(end).fAcc = [];
            sim_ave(end).mPos = [];
            sim_ave(end).mSpd = [];
            sim_ave(end).mAcc = [];
            sim_ave(end).n = 0;
        end
    end
end
sim_ave = sim_ave(2:end);

%% simulation
for i = 1:length(simulation)
    % select data for simulation
    data = simulation(i).data;
    manipOnset = simulation(i).manipOnset;
    t_start = int32((manipOnset-0.5)*Hz)+1;
    t_end = t_start + time*Hz - 1;
    simulation(i).t_start = t_start;
    simulation(i).t_end = t_end;
    lPos = data(t_start:t_end,1);
    lSpd = data(t_start:t_end,3);
    lAcc = data(t_start:t_end,5);
    fPos = data(t_start:t_end,2);
    fSpd = data(t_start:t_end,4);
    fAcc = data(t_start:t_end,6);
    inputHz = 60;
    outputHz = 60;
    x_start = fPos(1);
    v_start = fSpd(1);
    
    % simulate
    [mPos, mSpd, mAcc] = models(model, p, delay, lPos, lSpd, x_start, v_start, inputHz, outputHz);
    simulation(i).simulation = [mPos mSpd mAcc];
    
    % test
    simulation(i).r_v = corr(fSpd, mSpd, 'type','pearson');
    simulation(i).r_a = corr(fAcc, mAcc, 'type','pearson');
    simulation(i).z_v = atanh(simulation(i).r_v);
    simulation(i).z_a = atanh(simulation(i).r_a);
    simulation(i).RMSE_x = sqrt(mean((fPos - mPos).^2));
    simulation(i).RMSE_v = sqrt(mean((fSpd - mSpd).^2));
    simulation(i).RMSE_a = sqrt(mean((fAcc - mAcc).^2));
    
    % average simulation for each condition
    for j = 1:length(sim_ave)
        if sim_ave(j).d0 == simulation(i).d0 && sim_ave(j).v0 == simulation(i).v0 && sim_ave(j).dv == simulation(i).dv
            sim_ave(j).data_ave = sim_ave(j).data_ave + [lPos lSpd lAcc fPos fSpd fAcc];
            sim_ave(j).sim_ave = sim_ave(j).sim_ave + [mPos mSpd mAcc];
            sim_ave(j).lPos = [sim_ave(j).lPos lPos];
            sim_ave(j).lSpd = [sim_ave(j).lSpd lSpd];
            sim_ave(j).lAcc = [sim_ave(j).lAcc lAcc];
            sim_ave(j).fPos = [sim_ave(j).fPos fPos];
            sim_ave(j).fSpd = [sim_ave(j).fSpd fSpd];
            sim_ave(j).fAcc = [sim_ave(j).fAcc fAcc];
            sim_ave(j).mPos = [sim_ave(j).mPos mPos];
            sim_ave(j).mSpd = [sim_ave(j).mSpd mSpd];
            sim_ave(j).mAcc = [sim_ave(j).mAcc mAcc];
            sim_ave(j).n = sim_ave(j).n + 1;
        end
    end
end



% average simulation and calculate confidence interval
for i = 1:length(sim_ave)
    sim_ave(i).data_CI = t_critical/sqrt(sim_ave(i).n) * [std(sim_ave(i).fPos,0,2),...
                                                          std(sim_ave(i).fSpd,0,2),...
                                                          std(sim_ave(i).fAcc,0,2)];
    sim_ave(i).sim_CI = t_critical/sqrt(sim_ave(i).n) * [std(sim_ave(i).mPos,0,2),...
                                                          std(sim_ave(i).mSpd,0,2),...
                                                          std(sim_ave(i).mAcc,0,2)];
    sim_ave(i).data_ave = sim_ave(i).data_ave ./ sim_ave(i).n;
    sim_ave(i).sim_ave = sim_ave(i).sim_ave ./ sim_ave(i).n;
end
time = string(clock);
save(['sim_test_', model(5:end), '[', num2str(p), ']', char(join(time(1:5),'-')), '.mat'], 'model','p','delay','Hz','time','simulation','sim_ave');

%% print result

['Test results of ' model(5:end) ' p = [' num2str(p) ']:' newline ...
    'The Pearson''s r on v = ' num2str(mean([simulation.r_v])) newline ...
    'The Pearson''s r on a = ' num2str(mean([simulation.r_a])) newline ...
    'The RMSE on x = ' num2str(mean([simulation.RMSE_x])) newline ...
    'The RMSE on v = ' num2str(mean([simulation.RMSE_v])) newline ...
    'The RMSE on a = ' num2str(mean([simulation.RMSE_a]))]


%% Functions
function [mPos, mSpd, mAcc] = models(model, p, delay, lPos, lSpd, p_start, v_start, inputHz, outputHz)

    if strcmp(model,'fit_speedModel')
        c = p(1);
        [mPos, mSpd, mAcc] = fit_speedModel(p_start, v_start, inputHz,outputHz, lSpd, c, delay);

    elseif strcmp(model, 'fit_distanceModel')
        c = p(1);
        [mPos, mSpd, mAcc] = fit_distanceModel(p_start, v_start, inputHz,outputHz, lPos, c, delay);

    elseif strcmp(model, 'fit_expansionModel')
        b = p(1);
        w = 0.6; % the width of the target pole is 0.6 meter
        [mPos, mSpd, mAcc] = fit_expansionModel(p_start, v_start,inputHz,outputHz, lPos, lSpd, b, w, delay);

    elseif strcmp(model,'fit_speedModel+damping')
        c = p(1);
        d = p(2);
        [mPos, mSpd, mAcc] = fit_speedModel_damping(p_start, v_start, inputHz,outputHz, lSpd, c, d, delay);

    elseif strcmp(model,'fit_speedModel+delay')
        c = p(1);
        RT = p(2);
        [mPos, mSpd, mAcc] = fit_speedModel_delay(p_start, v_start, inputHz,outputHz, lSpd, c, RT);    

    elseif strcmp(model,'fit_distanceModel+damping')
        c = p(1);
        d = p(2);
        [mPos, mSpd, mAcc] = fit_distanceModel_damping(p_start, v_start, inputHz,outputHz, lPos, c, d, delay);

    elseif strcmp(model, 'fit_sDistanceModel')
        c = p(1);
        alpha = p(2);
        beta = p(3);
        [mPos, mSpd, mAcc] = fit_sDistanceModel(p_start, v_start, inputHz,outputHz, lPos, c, alpha, beta, delay);

    elseif strcmp(model, 'fit_ratioModel')
        c = p(1);
        M = p(2);
        L = p(3);
        [mPos, mSpd, mAcc] = fit_ratioModel(p_start, v_start,inputHz,outputHz, lPos, lSpd, c, M, L, delay);

    elseif strcmp(model, 'fit_linearModel')
        c1 = p(1);
        c2 = p(2);
        alpha = p(3);
        beta = p(4);
        [mPos, mSpd, mAcc] = fit_linearModel(p_start, v_start,inputHz,outputHz, lPos, lSpd, c1, c2, alpha, beta, delay); 

    elseif strcmp(model, 'fit_bruneauModel')
        ttr = p(1);
        df = p(2);
        [mPos, mSpd, mAcc] = fit_bruneauModel(p_start, v_start,inputHz,outputHz, lPos, ttr, df);

    elseif strcmp(model, 'fit_lemercierModel')
        RT = p(1);
        C = p(2);
        gamma = p(3);
        [mPos, mSpd, mAcc] = fit_lemercierModel(p_start, v_start,inputHz,outputHz, lPos, lSpd, RT, C, gamma);

    end
end

function [fPos, fSpd, fAcc] = fit_speedModel(p_start, v_start,inputHz,outputHz, lSpd, c, delay)
% This function can simulate the trajectory of follower given that of the
% leader. 
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% RT: the delay before the follower react to the speed change of the
% leader.
% c: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                           x..f = c * [x.l - x.f]
% 
% x..f: current acceleration of follower
% c: free parameter
% x.l: current speed of leader
% x.f: current speed of follower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lSpd = following(1).data(:,3);
% c = 1.6;
% p_start = 0;
% v_start = 0;
% inputHz = 60;
% outputHz = 60;
Hz = 6000;
time = length(lSpd)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lSpd_extrap = interp1(x,lSpd,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = zeros(nFrame,1);
for i = 1:nFrame
    fPos(i) = p_start;
end
fSpd = zeros(nFrame,1);
for i = 1:nFrame
    fSpd(i) = v_start;
end
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = c * (lSpd_extrap(i-int32(delay*Hz)) - fSpd(i));
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end

function [fPos, fSpd, fAcc] = fit_distanceModel(p_start, v_start,inputHz,outputHz, lPos, c, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% RT: the delay before the follower react to the distance change between
% the follower and the leader.
% c: free parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                     x..f = c * [delta_x - delta_x0]
% 
% x..f: the current acceleration of follower
% c: free parameter
% delta_x: the current distance between leader and follower
% delta_x0 = d0: initial distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPos_extrap = interp1(x,lPos,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = zeros(nFrame,1);
for i = 1:nFrame
    fPos(i) = p_start;
end
fSpd = zeros(nFrame,1);
for i = 1:nFrame
    fSpd(i) = v_start;
end
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!    
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = c * (lPos_extrap(i-int32(delay*Hz)) - fPos(i) - lPos_extrap(1));
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end

function [fPos, fSpd, fAcc] = fit_expansionModel( p_start, v_start,inputHz,outputHz, lPos, lSpd, b, w, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% b: free parameter
% w, h: width and height of the target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                x..f = -b * alpha.
% 
% x..f(t): the current acceleration of follower
% alpha.: the changing rate of visual angle
% b: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPos_extrap = interp1(x,lPos,xq,'linear','extrap');
lSpd_extrap = interp1(x,lSpd,xq,'linear','extrap');

fPos = zeros(nFrame,1);
for i = 1:nFrame
    fPos(i) = p_start;
end
fSpd = zeros(nFrame,1);
for i = 1:nFrame
    fSpd(i) = v_start;
end
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = b * w * (lSpd_extrap(i-int32(delay*Hz)) - fSpd(i))/((lPos_extrap(i-int32(delay*Hz)) - fPos(i))^2 + w^2/4);    
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end

function [fPos, fSpd, fAcc] = fit_sDistanceModel (p_start, v_start,inputHz,outputHz, lPos, c, alpha, beta, delay)
% This function can simulate the trajectory of follower given that of the
% leader. 
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% RT: the delay before the follower react to the speed change of the
% leader.
% c, alpha, beta: free parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                  x..f = c * [delta_x - alpha - beta*x.f]
% 
% x..f: the current acceleration of follower
% c, alpha, beta: free parameter
% delta_x: the current distance between leader and follower
% x.f: the current speed of the follower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if beta < 0
    beta = 0;
end
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPos_extrap = interp1(x,lPos,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = zeros(nFrame,1);
for i = 1:nFrame
    fPos(i) = p_start;
end
fSpd = zeros(nFrame,1);
for i = 1:nFrame
    fSpd(i) = v_start;
end
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = c * (lPos_extrap(i-int32(delay*Hz)) - fPos(i) - alpha - beta * fSpd(i));
end
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end

function [fPos, fSpd, fAcc] = fit_ratioModel(p_start, v_start,inputHz,outputHz, lPos, lSpd, c, M, L, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% RT: the delay before the follower react to the distance change between
% the follower and the leader.
% c, M, L: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                x..f(t) = c * x.f^M * delta_x. / delta_x^L
% 
% x..f(t): the current acceleration of follower
% x.f: the speed of the follower
% delta_x: the current distance between leader and follower
% delta_x.: the rate of change of the current distance between leader and follower
% c, M, L: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPos_extrap = interp1(x,lPos,xq,'linear','extrap');
lSpd_extrap = interp1(x,lSpd,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = zeros(nFrame,1);
for i = 1:nFrame
    fPos(i) = p_start;
end 
fSpd = zeros(nFrame,1);
for i = 1:nFrame
    fSpd(i) = v_start;
end
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = c * fSpd(i)^M * (lSpd_extrap(i-int32(delay*Hz)) - fSpd(i)) / (lPos_extrap(i-int32(delay*Hz)) - fPos(i))^L; 
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end

function [fPos, fSpd, fAcc] = fit_linearModel(p_start, v_start,inputHz,outputHz, lPos, lSpd, c1, c2, alpha, beta, delay)
% This function can simulate the trajectory of follower given that of the
% leader. 
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% c1, c2, alpha, beta: free parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%            x..f = c1*[delta_x.] + c2*[delta_x - alpha - beta*x.f]
% 
% x..f: current acceleration of follower
% delta_x.: relative speed
% delta_x: relative distance
% x.f: follower speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if beta < 0
    beta = 0;
end
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPos_extrap = interp1(x,lPos,xq,'linear','extrap');
lSpd_extrap = interp1(x,lSpd,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = zeros(nFrame,1);
for i = 1:nFrame
    fPos(i) = p_start;
end 
fSpd = zeros(nFrame,1);
for i = 1:nFrame
    fSpd(i) = v_start;
end
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = c1 * (lSpd_extrap(i-int32(delay*Hz)) - fSpd(i)) + c2 * (lPos_extrap(i-int32(delay*Hz)) - fPos(i) - alpha - beta * fSpd(i));
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end


function [fPos, fSpd, fAcc] = fit_lemercierModel(p_start, v_start,inputHz,outputHz, lPos, lSpd, RT, C, gamma)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% RT: the delay before the follower react to the distance change between
% the follower and the leader.
% C, gamma, RT: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%          x..f(t) = C * delta_x.(t - RT) * (1/delta_x(t))^gamma
% 
% x..f(t): the current acceleration of follower
% x.f: the speed of the follower
% delta_x: the current distance between leader and follower
% delta_x.: the rate of change of the current distance between leader and follower
% C, gamma: free parameter
% RT: reaction time, range [0 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RT < 0
    RT = 0;
end
if RT > 1
    RT = 1;
end

Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPos_extrap = interp1(x,lPos,xq,'linear','extrap');
lSpd_extrap = interp1(x,lSpd,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = zeros(nFrame,1);
for i = 1:nFrame
    fPos(i) = p_start;
end 
fSpd = zeros(nFrame,1);
for i = 1:nFrame
    fSpd(i) = v_start;
end
fAcc = zeros(nFrame,1);

for i = 2+int32(RT*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = C * (lSpd_extrap(i-int32(RT*Hz)) - fSpd(i-int32(RT*Hz))) * (1/(lPos_extrap(i) - fPos(i)))^gamma;
end
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end



function [fPos, fSpd, fAcc] = fit_bruneauModel(p_start, v_start,inputHz,outputHz, lPos, ttr, df)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% delta_t, df, ttr: free parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%          x.f(t) = (x_l(t+delta_t) - x_f(t) - df) / (delta_t + ttr)
% 
% x.f: follower's speed (m/s)
% delta_t: time step (s)
% ttr: time to react (s)
% df: sum of body size and personal distance (m)
% x_f: follower's position (m)
% x_l: leader's position (m)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ttr < 0
    ttr = 0;
end
if ttr > 1
    ttr = 1;
end
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPos_extrap = interp1(x,lPos,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = zeros(nFrame,1);
for i = 1:nFrame
    fPos(i) = p_start;
end
fSpd = zeros(nFrame,1);
for i = 1:nFrame
    fSpd(i) = v_start;
end
fAcc = zeros(nFrame,1);

for i = 2:nFrame-1
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i-1) + fSpd(i-1)/Hz;
    fSpd(i) = (lPos_extrap(i+1) - fPos(i) - df)/(1/Hz + ttr);    
    fAcc(i) = (fSpd(i) - fSpd(i-1))*Hz;
end
fPos(nFrame) = fPos(nFrame-1);
fSpd(nFrame) = fSpd(nFrame-1);
fAcc(nFrame) = fAcc(nFrame-1);
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end
