% load trials
n_subject = 13;
n_trial = 90;
subFolder = 'Carrot3\';

loading = dir(strcat(subFolder,'orientation*')); % load all txt files in directory of code

orientation = struct;

for i = 1:length(loading) % iterate through all trials
    
    if ~contains(loading(i).name, 'freewalk')
        
        fileName_ = strcat(subFolder, loading(i).name); % find name of i_trial file
        orientation(i).orientation = load(fileName_);
    end
   
end

save('data_orientation.mat', 'orientation');