% selected data range is 0.5s before perturbation. Because the onset of
% perturbation is between 2-3 s, the begining of selected data is between 
% 1.5s to 2.5s. Because the selected data should end no later than 8.5s the
% length should be no longer than 8.5s - 2.5s = 6s.
clear all;
load data_piped;

%% set up the structures


Hz = 60;
time = 6; %seconds


subject_ave = struct;
for i = 1:13
    for j = [1 4 8]
        for k = [0.8 1.2]
            for l = [-0.3 0 0.3]
                subject_ave(end+1).subject = i;
                subject_ave(end).d0 = j;
                subject_ave(end).v0 = k;
                subject_ave(end).dv = l;
                subject_ave(end).data = zeros(time*Hz,6);
                subject_ave(end).n = 0;
            end
        end
    end
end
subject_ave = subject_ave(2:end);

condition_ave = struct;
for j = [1 4 8]
    for k = [0.8 1.2]
        for l = [-0.3 0 0.3]
            condition_ave(end+1).d0 = j;
            condition_ave(end).v0 = k;
            condition_ave(end).dv = l;
            condition_ave(end).data = zeros(time*Hz,6);
            condition_ave(end).n = 0;
        end
    end
end
condition_ave = condition_ave(2:end);


%% loop through all trials to get the sum
for i = 1:length(following)
    subject = following(i).subject;
    trial = following(i).trial;
    d0 = following(i).d0;
    v0 = following(i).v0;
    dv = following(i).dv;
    dump = following(i).dump;
    data = following(i).data;
    t_total = following(i).t_total;
    manipOnset = following(i).manipOnset;

    if dump == 0 && t_total > 8.5
        
        t_start = int32((manipOnset-0.5)*Hz) + 1;
        t_end = int32((manipOnset-0.5)*Hz) + time * Hz;
        

        for j = 1:length(subject_ave)
            if subject == subject_ave(j).subject &&...
                    d0 == subject_ave(j).d0 &&...
                    v0 == subject_ave(j).v0 &&...
                    dv == subject_ave(j).dv
                
                subject_ave(j).data = subject_ave(j).data + data(t_start:t_end,1:6);
                subject_ave(j).n = subject_ave(j).n + 1;                
            end       
        end
        
        for j = 1:length(condition_ave)
            if d0 == condition_ave(j).d0 &&...
               v0 == condition_ave(j).v0 &&...
               dv == condition_ave(j).dv
                
                condition_ave(j).data = condition_ave(j).data + data(t_start:t_end,1:6);
                condition_ave(j).n = condition_ave(j).n + 1;                
            end       
        end
    end 
end


%% divid by n to get the average

for i = 1:length(subject_ave)
    subject_ave(i).data = subject_ave(i).data ./ subject_ave(i).n;
end

for i = 1:length(condition_ave)
    condition_ave(i).data = condition_ave(i).data ./ condition_ave(i).n;
end

save('data_averaged_fullLength.mat','subject_ave','condition_ave');