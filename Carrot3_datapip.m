clear all;
leader_traj_source = 'trajectoryGenerator'; % 'trajectoryGenerator'/'input'/'output'

%% loading
n_subject = 13;
n_trial = 90;
subFolder = 'Carrot3_rawData\';

% load test trials
loading1 = dir(strcat(subFolder,'position*')); % load all txt files in directory of code
following = struct;


for i_trial = 1:length(loading1) % iterate through all trials
    
    if ~contains(loading1(i_trial).name, 'freewalk')
        
        s_ = strfind(loading1(i_trial).name, 'subj');
        following(i_trial).subject = str2double(loading1(i_trial).name(s_+5 : s_+6));

        if strcmp(loading1(i_trial).name(end-5), '_')
            following(i_trial).trial = str2double(loading1(i_trial).name(end-4));
        else
            following(i_trial).trial = str2double(loading1(i_trial).name(end-5:end-4));
        end
        
        following(i_trial).d0 = str2double(loading1(i_trial).name(10));
        following(i_trial).v0 = str2double(loading1(i_trial).name(13:15));

        if strcmp(loading1(i_trial).name(19), '0') 
            following(i_trial).dv = -0.3;

        elseif strcmp(loading1(i_trial).name(19), '.') 
            following(i_trial).dv = 0.3;

        elseif strcmp(loading1(i_trial).name(19), '_')
            following(i_trial).dv = 0;
        end
        
        fileName_ = strcat(subFolder, 'Carrot3_input\Subject_', num2str(following(i_trial).subject,'%02d'),...
            '\conditions.csv');
        conditionFile = csvread(fileName_,1,0); % read data start from row 2 column 1
        if following(i_trial).d0 == conditionFile(following(i_trial).trial,2) &&...
           following(i_trial).v0 == conditionFile(following(i_trial).trial,3) &&...
           following(i_trial).dv == conditionFile(following(i_trial).trial,4)
            following(i_trial).manipOnset = conditionFile(following(i_trial).trial,5);
        end
        
        fileName_ = strcat(subFolder, loading1(i_trial).name); % find name of i_trial file
        raw_traj = load(fileName_);
        following(i_trial).t_total = raw_traj(end,5);
        following(i_trial).raw_traj = raw_traj; 
        
    end
   
end

% load freewalk trials
loading2 = dir(strcat(subFolder,'position_freewalk*')); % load all txt files in directory of code
freewalk = struct;

for i_trial = 1:length(loading2) % iterate through all trials
    
    s_ = strfind(loading2(i_trial).name, 'subj');
    freewalk(i_trial).subject = str2double(loading2(i_trial).name(s_+5 : s_+6));

    freewalk(i_trial).session = str2double(loading2(i_trial).name(20));
        
    if strcmp(loading2(i_trial).name(end-5), '_')
        freewalk(i_trial).trial = str2double(loading2(i_trial).name(end-4));
    else
        freewalk(i_trial).trial = str2double(loading2(i_trial).name(end-5:end-4));
    end
    
    fileName_ = strcat(subFolder, loading2(i_trial).name); % find name of i_trial file
    freewalk(i_trial).raw_traj = load(fileName_); % load in specific file
    
end


% load error notes
loading3 = dir(strcat(subFolder,'Dump*')); % load all txt files in directory of code
dump = struct;

for i_trial = 1:length(loading3)    
    dump(i_trial).subject = str2double(loading3(i_trial).name(14:15));
    dump(i_trial).trial = str2double(loading3(i_trial).name(23:24));
    dump(i_trial).error = loading3(i_trial).name(26:end-4);   
end


%% processing following data
theta = atand(9/11);
Hz = 60;
t = 2; % used the last 2 second to calculate ending speed
sample = Hz*t; % sample of data points used to compute ending speed


for i_trial = 1:length(following)
    
    subject = following(i_trial).subject;
    trial = following(i_trial).trial;
    d0 = following(i_trial).d0;
    v0 = following(i_trial).v0;
    dv = following(i_trial).dv;
    manipOnset = following(i_trial).manipOnset;
    
    % make initial position zero
    pos_0 = following(i_trial).raw_traj(1,3:4);
    following(i_trial).zero_traj(:, 1:2) = following(i_trial).raw_traj(:, 1:2) - pos_0;
    following(i_trial).zero_traj(:, 3:4) = following(i_trial).raw_traj(:, 3:4) - pos_0;
    following(i_trial).zero_traj(:, 5) = following(i_trial).raw_traj(:, 5); % copy time stamp
    

    % rotate the trajectory to 1 dimension (y)
    following(i_trial).rot_traj(:, 1:2) = rotate2D(theta, following(i_trial).zero_traj(:, 1:2), [0, 0]);
    following(i_trial).rot_traj(:, 3:4) = rotate2D(theta, following(i_trial).zero_traj(:, 3:4), [0, 0]);
    % unify two directions of walking
    if following(i_trial).rot_traj(end,4) < 0
       following(i_trial).rot_traj(:, 1:4) = following(i_trial).rot_traj(:, 1:4)*(-1); 
    end
    following(i_trial).rot_traj(:, 5) = following(i_trial).raw_traj(:, 5); % copy time stamp
    
 
    % interpolate the trajectory to 60 Hz data with equal time interval 
    x = following(i_trial).raw_traj(:, 5);
    v = following(i_trial).rot_traj(:, 1:4);
    xp = 1/Hz:1/Hz:following(i_trial).raw_traj(end,5);
    vq = interp1(x,v,xp,'linear','extrap');
    following(i_trial).inter_traj(:, 1:4) = vq;
    following(i_trial).inter_traj(:, 5) = xp;
    

    

    
    if strcmp(leader_traj_source, 'input')
        % load and use leader's trajectory from input files

        fileName_ = strcat(subFolder, 'Inputs\Subject_', num2str(subject,'%02d'),...
            '\trial',num2str(trial,'%03d'),'.csv');
        inputFile = load(fileName_);
        input_traj = inputFile(:,2);

        trial_length = 12;
        Hz = 60;
        inputHz = 90;

        % interpolate the trajectory to 60 Hz data with equal time interval 
        x = 1/inputHz:1/inputHz:trial_length;
        x = x';
        v = input_traj;
        xp = 1/Hz:1/Hz:trial_length;
        xp = xp';
        vq = interp1(x,v,xp,'linear','extrap');
        inter_traj = [vq ; repmat(vq(end),[60 1])];

        % replace leader's trajectory from vizard by that from the input file
        len = size(following(i_trial).inter_traj,1);
        following(i_trial).inter_traj(:,2) = inter_traj(1:len);
   
    end
    
    
    if strcmp(leader_traj_source, 'trajectoryGenerator')
        
        x0 = 0;
        y0 = 0;
        nDuration = 13;
        a = 1;
        heading1 = 0;
        heading2 = 0;
        startupDuration = 0; % how fast the pole start to move

        [x, y, spd, hdn, manipOnset] = Carrot3_trajectoryGenerator(x0,y0,nDuration,...
            d0,v0,dv,a,heading1,heading2,startupDuration,manipOnset,Hz);
        
        len_ = size(following(i_trial).inter_traj,1);
        following(i_trial).inter_traj(:,2) = y(1:len_);
    end
    
    % Butterworth filter
    following(i_trial).filtered_traj(:, 2) = Carrot3_filter_butter(Hz, following(i_trial).inter_traj(:, 4));
    following(i_trial).filtered_traj(:, 1) = following(i_trial).inter_traj(:, 2); % copy leader's traj
    following(i_trial).filtered_traj(:, 3) = following(i_trial).inter_traj(:, 5); % copy time stamp after interpolation
    
    
    % calculate speed and acceleration
    % position
    following(i_trial).data(:, 1) = following(i_trial).filtered_traj(:, 1);% leader
    following(i_trial).data(:, 2) = following(i_trial).filtered_traj(:, 2);% follower
    
    % speed
    d = diff(following(i_trial).data(:, 1));
    following(i_trial).data(:, 3) = [d; d(end)] * Hz;% leader
    d = diff(following(i_trial).data(:, 2));
    following(i_trial).data(:, 4) = [d; d(end)] * Hz;% follower
    
    
    
    % use neighbor points on both sides to calcuate speed
%     for i = 1:length(following(i_trial).data(:, 2))
%         if i == 1
%             p_pre = 0;
%         else
%             p_pre = following(i_trial).filtered_traj(i-1, 2);
%         end
%         
%         if i == length(following(i_trial).data(:, 2))
%             p_post = following(i_trial).filtered_traj(end, 2);
%         else
%             p_post = following(i_trial).filtered_traj(i+1, 2);
%         end
%         
%         following(i_trial).data(i, 4) = (p_post - p_pre) * Hz / 2;
%     end
    

    % acceleration
    d = diff(following(i_trial).data(:, 3));
    following(i_trial).data(:, 5) = [d;d(end)] * Hz;% leader
    d = diff(following(i_trial).data(:, 4));
    following(i_trial).data(:, 6) = [d; d(end)] * Hz;% follower
    
    % time stamp
    following(i_trial).data(:, 7) = following(i_trial).filtered_traj(:, 3);
    
    t_start = int32((manipOnset-0.5)*Hz)+1;
    t_end = size(following(i_trial).data,1) - 1.0*Hz;
    
    % mark thrown trials
    % mark all trials as good first
    following(i_trial).dump = 0;
    
    % mark thrown trials by error notes
    for j_trial = 1:length(dump)
        if following(i_trial).subject == dump(j_trial).subject && following(i_trial).trial == dump(j_trial).trial
            following(i_trial).dump = 1;  
        end
    end
     
    % ending speed
    following(i_trial).finalSpd = mean(following(i_trial).data(t_end-sample+1:t_end, 4));
    following(i_trial).finalDv = mean(following(i_trial).data(t_end-sample+1:t_end, 3)...
                                - following(i_trial).data(t_end-sample+1:t_end,4));
    following(i_trial).finalDist = mean(following(i_trial).data(t_end-sample+1:t_end,1)...
                                - following(i_trial).data(t_end-sample+1:t_end,2));
end





%% Processing freewalk

for i_trial = 1:length(freewalk)
    
    % make initial position zero
    pos_0 = freewalk(i_trial).raw_traj(1,3:4);
    freewalk(i_trial).zero_traj(:, 1:2) = freewalk(i_trial).raw_traj(:, 3:4) - pos_0;
    freewalk(i_trial).zero_traj(:, 3) = freewalk(i_trial).raw_traj(:, 5); % copy time stamp
    
    % rotate the trajectory to 1 dimension (y)
    freewalk(i_trial).rot_traj(:, 1:2) = rotate2D(theta, freewalk(i_trial).zero_traj(:, 1:2), [0, 0]);
    if freewalk(i_trial).rot_traj(end,2) < 0
       freewalk(i_trial).rot_traj(:, 1:2) = freewalk(i_trial).rot_traj(:, 1:2)*(-1); 
    end
    freewalk(i_trial).rot_traj(:, 3) = freewalk(i_trial).raw_traj(:, 5); % copy time stamp
    
    % interpolate the trajectory to 60 Hz data with equal time interval 
    x = freewalk(i_trial).raw_traj(:, 5);
    v = freewalk(i_trial).rot_traj(:, 1:2);
    xp = 1/Hz:1/Hz:freewalk(i_trial).raw_traj(end,5);
    vq = interp1(x,v,xp,'linear','extrap');
    freewalk(i_trial).inter_traj(:, 1:2) = vq;
    freewalk(i_trial).inter_traj(:, 3) = xp;
    
    % Butterworth filter
    freewalk(i_trial).filtered_traj(:, 1) = Carrot3_filter_butter(Hz, freewalk(i_trial).inter_traj(:, 2));
    freewalk(i_trial).filtered_traj(:, 2) = freewalk(i_trial).inter_traj(:, 3); % copy time stamp after interpolation
    
    
    % calculate speed and acceleration
    % position
    freewalk(i_trial).data(:, 1) = freewalk(i_trial).filtered_traj(:, 1);
    
    % speed
    d = diff(freewalk(i_trial).data(:, 1));
    freewalk(i_trial).data(:, 2) = [d; d(end)] * Hz;
    
    % acceleration
    d = diff(freewalk(i_trial).data(:, 2));
    freewalk(i_trial).data(:, 3) = [d; d(end)] * Hz;
    
    % time stamp
    freewalk(i_trial).data(:, 4) = freewalk(i_trial).filtered_traj(:, 2);
    
    % ending speed
    t_end = size(freewalk(i_trial).data,1) - 1.0*Hz;
    freewalk(i_trial).finalSpd = mean(freewalk(i_trial).data(t_end-sample+1:t_end, 2));
    
end



%% save
time = string(clock);
save(['Carrot3_data_piped(4th1Hz)',char(join(time(1:5),'-')),'.mat'],'following','freewalk');



