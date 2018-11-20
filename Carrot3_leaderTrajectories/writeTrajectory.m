function [] = writeTrajectory(Nfo, Pos, Spd, frameRate)
% This function takes an array and write it out to text file for Vizard.

nTrial = length(Pos);

fileID = fopen('conditions.csv', 'w');
fprintf(fileID ,'%s\n','Trial,d0,v0,dv,manipStartTime');

for iTrial = 1:nTrial
    % Export positions to text file for Vizard.
    % filename = ['input/trial',num2str(iTrial),'.txt']

    filename = ['trial',num2str(iTrial,'%03d'),'.csv'];

    combo = horzcat(Pos{iTrial,1},Spd{iTrial,1});
    csvwrite(filename, combo);
    

    fileID = fopen('conditions.csv', 'a');
    fprintf(fileID ,'%s\n', strcat(num2str(iTrial), ',', num2str(Nfo{iTrial,1}.('d0')), ',',...
                                num2str(Nfo{iTrial,1}.('v0')), ',',...
                                num2str(Nfo{iTrial,1}.('dv')), ',',...
                                num2str(double(Nfo{iTrial,1}.('manipOnset'))/frameRate, 3), ','));
end
