clear all;
load Carrot3_data_piped;
load Carrot3_data_averaged;

%% plot filtered freewalk 
for i = 1:length(freewalk)
    hold on;
    plot(freewalk(i).data(:,4), freewalk(i).data(:,2));
end

%% plot input traj
figure;
for i = 1:length(inputs)
    if inputs(i).subject == 1
        hold on;
        speed = diff(inputs(i).traj);
        speed = [speed; speed(end)];
        speed = speed * 90;
        plot(1/90:1/90:12, speed);
    %         axis([0 12 -5 5]);
    end

end

%% plot rotated data

figure;
for i = 1:length(following)
    hold on;
    if following(i).dump == 0 && following(i).subject >= 1 
        speed = diff(following(i).rot_traj(:,4));
        speed = [speed; speed(end)];
        speed = speed * 90;
        plot(following(i).rot_traj(:,5), speed);
%         axis([0 12 -5 5]);
    end
end

%% plot interpolated data
figure;
for i = 1:length(following)
    hold on;
    if following(i).dump == 0 && following(i).subject >= 1 
        speed = diff(following(i).inter_traj(:,4));
        speed = [speed; speed(end)];
        speed = speed * 60;
        plot(following(i).inter_traj(:,5), speed);
%         axis([0 12 -5 5]);
    end
end

%% plot filtered data
red = false;
blue = false;
green = false;
figure;
for i = 1:length(following)
    subject = following(i).subject;
    trial = following(i).trial;
    d0 = following(i).d0;
    v0 = following(i).v0;
    dv = following(i).dv;
    dump = following(i).dump;
    data = following(i).data;
    t = following(i).t_total;
    manipOnset = following(i).manipOnset;
    t_start = int32(manipOnset*60) - 30;
    
    hold on;
    if dump == 0 && subject >= 1 && dv ~= 0 %&& d0 == 8
        if d0 == 1
            if ~red
                pr = plot(data(:,7), data(:,4), 'r','DisplayName','1 m');
                red = true;
            else
                plot(data(:,7), data(:,4), 'r');
            end
        end
        if d0 == 4
            if ~blue
                pb = plot(data(:,7), data(:,4), 'b','DisplayName','4 m');
                blue = true;
            else
                plot(data(:,7), data(:,4), 'b');
            end
        end
        if d0 == 8
            if ~green            
                pg = plot(data(:,7), data(:,4), 'g','DisplayName','8 m'); 
                green = true;
            else
                plot(data(:,7), data(:,4), 'g'); 
            end
        end
        
        
%         plot([2 2],[-0.5 3],'--','Color',[0.3 0.3 0.3]),
%         plot([3 3],[-0.5 3],'--','Color',[0.3 0.3 0.3]),

        
        %         axis([0 12 0.4 1.6]);
        
%         if t < 9        
%             plot(data(1:end-60,7), data(1:end-60,4),'r');
%         end
%         
%         if t > 9
%             plot(data(1:end-60,7), data(1:end-60,4), 'b');
%         end
        
        
        
      
%         if d0 == 1 && dv == 0.3
% 
%             plot(data(:,7), data(:,6),'r');
%         end
%         
%         if d0 == 4 && dv == 0.3
%             plot(data(:,7), data(:,6),'b');
%         end
%         
%         if d0 == 8 && dv == 0.3
%             plot(data(:,7), data(:,6),'g');
%         end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%         
%         if d0 == 1 && dv == 0
%             plot(data(:,7), data(:,4),'r--');
%         end
%         
%         if d0 == 4 && dv == 0
%             plot(data(:,7), data(:,4),'b--');
%         end
%         
%         if d0 == 8 && dv == 0
%             plot(data(:,7), data(:,4),'g--');
%         end
%         
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%         
%         if d0 == 1 && dv == -0.3
%             plot(data(:,7), data(:,6),'r:');
%         end
%         
%         if d0 == 4 && dv == -0.3
%             plot(data(:,7), data(:,6),'b:');
%         end
%         
%         if d0 == 8 && dv == -0.3
%             plot(data(:,7), data(:,6),'g:');
%         end
% %         
    end
    
end
hold off;
legend([pr pb pg],'1 m','4 m','8 m');
legend('boxoff')
axis([0 12 0 2.5]);
xlabel('Time(s)');
ylabel('Speed(m/s)');
% title('Acceleration trials');
title('Perturbation trials');
% title('Deceleration trials');
ax = gca;
ax.FontSize = 20;
%% Plot freewalk data
figure;
hold on;
for i=1:length(freewalk)
    plot(freewalk(i).data(:,4), freewalk(i).data(:,2));
end


%% plot for all subject in individual graph

for i = 1:13
    figure;
    for j = 1:length(following)
        hold on;
        if following(j).dump == 0 && following(j).subject == i
            plot(following(j).data(1:end-30,7), following(j).data(1:end-30,1),'.');
        end
    end
end



%% filter testing

Hz = 60;
figure;
for i = 1:length(following)
    hold on;
    if following(i).dump == 0 && following(i).subject >= 1 %&& following(i).trial >= 1 && following(i).raw_traj(end,5)<= 12

        vInput = following(i).inter_traj(:,4);
        cutoff = (0.6); % standard VENLab parameters

        [B,A] = butter(3,cutoff/(Hz/2));

        % length of extrapolation of one side
        pad = 120;
        
%         extension = vInput(end-pad+1:end);
%         extension1 = extension + (vInput(end)-vInput(end-pad));
%         extension2 = extension + (vInput(end)-vInput(end-pad))*2;
%         extendedInput = vertcat(vInput, extension1, extension2);
        % Pad data with timesteps to extrapolate into
        tot_time = 1:length(vInput)+pad*2;

        extrap_pre = 1:pad;
        extrap_post = length(tot_time)-pad+1:length(tot_time);

        exist_data_location = length(extrap_pre)+1:extrap_post(1)-1;
        exist_data_location = exist_data_location'; % ' means transpose

        toExtrap = vertcat(extrap_pre',extrap_post');

%         Extrapolate
        outExtrap = interp1(exist_data_location,vInput,toExtrap,'linear','extrap');
        outExtended = vertcat(outExtrap(1:pad),vInput,outExtrap(pad+1:end));

%         Apply butterworth filter and return input to original size
        outFiltered = filtfilt(B,A,outExtended); 
        vOutput = outFiltered(exist_data_location);

        
        
        speed = diff(vInput);
        speed = [speed; speed(end)];
        speed = speed * Hz;
        
        x = 1:length(vInput);
        plot(x, vInput);
%         % horizontal line
%         plot([0 840], [12.3 12.3],'k');
%         plot([0 840], [0 0],'k');
%         % vertical line
%         plot([840 840], [0 12.3], 'k');
%         plot([120 120], [0 12.3], 'k');
        
    end
end

%% histogram of trial length
t = [];
hold on;
for i = 1:length(following)
    if following(i).v0 == 1.2 && following(i).d0 == 8 && following(i).dv == 0.3
        t(end+1) = length(following(i).data);
    end
end
t = t/60;
histogram(t,10);
axis([7 13 0 500]);

%% trial info check (trial length analysis)
count = 0;
for i = 1:length(following)
    subject = following(i).subject;
    trial = following(i).trial;
    d0 = following(i).d0;
    v0 = following(i).v0;
    dv = following(i).dv;
    dump = following(i).dump;
    data = following(i).data;
    t = following(i).t_total;
    
    if dump == 0 && dv == 0.3
        count = count + 1;
%         [d0 v0 dv]
        
        
    end
end
count

%% plot subject_ave data
x = 1/60:1/60:6;
x = x-0.5;
for s = 1:13
    figure;
    for i = 1:length(subject_ave)       
        if subject_ave(i).subject == s
            hold on;
            fig = plot(x, subject_ave(i).data(:,4)),
            plot([0 0],[0 2],'--k'),
            axis([-0.5 5.5 0 2]);
        end
    end
    saveas(fig, [num2str(s),'.png']);
end
%% plot condition_ave data
x = 1/60:1/60:6;
x = x-0.5;
for i = 1:length(condition_ave)       
    hold on;
    fig = plot(x, condition_ave(i).data(:,4)),
    plot([0 0],[0 2],'--k'),
    axis([-0.5 5.5 0 2]);
end
saveas(fig,'average_across_subjects.png');

%% plot speed data against one model
x = 1/60:1/60:6;
x = x-0.5;
for i = 1:18
    d0 = num2str(model_ave(i).d0);
    v0 = num2str(model_ave(i).v0);
    dv = num2str(model_ave(i).dv);
    condition = ['Initial distance ' d0 '(m) initial leader speed ' v0 '(m/s) leader speed change ' dv '(m/s)'];
    figure;
    hold on;
    fig = plot(x, condition_ave(i).data(:,3), 'k', 'LineWidth',2),
    plot(x, condition_ave(i).data(:,4), 'r', 'LineWidth',2),
    plot(x, model_ave(i).data(:,2), '--r', 'LineWidth',2),
    plot([0 0],[0 2],'--','Color',[0.3 0.3 0.3]),
    text(0, 1.8, '\leftarrow leader speed change'),
    legend('leader','follower',model(5:end)),
    axis([-0.5 5.5 0 2]),
    xlabel('Time(s)'),
    ylabel('Speed(m/s)'),
    title(condition);
    saveas(fig,['data_vs_', model(5:end), 'p[', num2str(p), ']', '(d0=', d0, ' v0=', v0, ' dv=', dv, ').png']);
end

%% plot speed data against two models
x = 1/60:1/60:6;
x = x-0.5;
for i = 1:18
    d0 = num2str(model_ave(i).d0);
    v0 = num2str(model_ave(i).v0);
    dv = num2str(model_ave(i).dv);
    condition = ['Initial distance ' d0 '(m) initial leader speed ' v0 '(m/s) leader speed change ' dv '(m/s)'];
    if model_ave(i).d0 == 1 && model_ave(i).v0 == 1.2 && model_ave(i).dv == 0.3
        figure;
        hold on;
        fig = plot(x, condition_ave(i).data(:,3), 'k', 'LineWidth',2),...
            plot(x, condition_ave(i).data(:,4), 'r', 'LineWidth',2),...
            plot(x, model_ave1(i).data(:,2), '--r', 'LineWidth',2),...
            plot(x, model_ave2(i).data(:,2), '--', 'Color',[0 0.4470 0.7410], 'LineWidth',2),...
            plot([0 0],[0 2],'--','Color',[0.3 0.3 0.3]),...
            txt = text(0, 1.8, '\leftarrow perturbation'),...
            txt.FontSize = 12,...
            led = legend('leader','follower','speed model','visual model'),...
            led.FontSize = 14,...
            axis([-0.5 5.5 0 2]),...
            xlabel('Time(s)'),...
            ylabel('Speed(m/s)'),...
            title('Accelerating Trials');
            ax = gca,...
            ax.FontSize = 20;
%         saveas(fig,['data_vs_Speed&Expansion model', 'p[', num2str(p), ']', '(d0=', d0, ' v0=', v0, ' dv=', dv, ').png']);
    end
end

%% plot distance of two models

x = 1/60:1/60:6;
x = x-0.5;
for i = 1:18
    d0 = num2str(model_ave(i).d0);
    v0 = num2str(model_ave(i).v0);
    dv = num2str(model_ave(i).dv);
    condition = ['Initial distance ' d0 '(m) initial leader speed ' v0 '(m/s) leader speed change ' dv '(m/s)'];
    if model_ave(i).d0 == 1 && model_ave(i).v0 == 0.8 && model_ave(i).dv == -0.3
        figure;
        hold on;
        fig = plot(x, condition_ave(i).data(:,1) - condition_ave(i).data(:,2), 'r', 'LineWidth',2),...
            plot(x, condition_ave(i).data(:,1) - model_ave1(i).data(:,1), '--r', 'LineWidth',2),...
            plot(x, condition_ave(i).data(:,1) - model_ave2(i).data(:,1), '--', 'Color',[0 0.4470 0.7410], 'LineWidth',2),...
            plot([0 0],[0 12],'--','Color',[0.3 0.3 0.3]),...
            txt = text(0, 1.8, '\leftarrow perturbation'),...
            txt.FontSize = 12,...
            led = legend('follower','speed model','expansion model'),...
            led.FontSize = 14,...
            axis([-0.5 5.5 0 12]),...
            xlabel('Time(s)'),...
            ylabel('Distance(m)'),...
            title('Accelerating Trials');
            ax = gca,...
            ax.FontSize = 20;
%         saveas(fig,['data_vs_Speed&Expansion model', 'p[', num2str(p), ']', '(d0=', d0, ' v0=', v0, ' dv=', dv, ').png']);
    end
end


%% plot orientation

% load data_orientation;
figure;
hold on;
for i = 1:length(orientation)
    % plot yaw
    y = orientation(i).orientation(:,1);
    x = 1:length(y);
    x = x';
    plot(x,y);
    title('yaw');
    ax = gca;
    ax.FontSize = 20;
end

figure;
hold on;
for i = 1:length(orientation)
    % plot yaw
    y = orientation(i).orientation(:,2);
    x = 1:length(y);
    x = x';
    plot(x,y);
    title('pitch');
    ax = gca;
    ax.FontSize = 20;  
end
figure;
hold on;
for i = 1:length(orientation)
    % plot yaw
    y = orientation(i).orientation(:,3);
    x = 1:length(y);
    x = x';
    plot(x,y);
    title('roll');
    ax = gca;
    ax.FontSize = 20;
end

%% Check frame rate
FR = [];
for i = 1:length(following)
    FR = [FR;(following(i).raw_traj(2:end,5) - following(i).raw_traj(1:end-1,5)).^(-1)];
end

histogram(FR);
xticks(0:5:100);
axis([20 100 0 200000]),...
ax = gca;
ax.FontSize = 20;


%% histogram of effective trial length
t = [];
Hz = 60;
figure;
hold on;
for j = 1:length(following)
    if following(j).dump ~= 1 %ExperimentalTrials(i).w == 0.2 && ExperimentalTrials(i).dv == -0.3
        t(end+1) = following(j).t_total - following(j).manipOnset;
    end
end
histogram(t,10);
axis([0 13 0 500]);

%% trial info check (trial length analysis)
count = 0;
for j = 1:length(following)
    subject = following(j).subject;
    trial = following(j).trial;
    d0 = following(j).d0;
    v0 = following(j).v0;
    dv = following(j).dv;
    dump = following(j).dump;
    t = following(j).t_total;
    manipOnset = following(j).manipOnset;
    
    if dump ~= 1 && t-manipOnset<6
        count = count + 1;
        [d0 v0 dv]
    end
end
count
