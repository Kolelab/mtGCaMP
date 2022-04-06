%% Before you start

% Download the ephys traces and save them in the same folder as the tiff files
% Save ephys traces as ephys_traces!!

% Add path to where this script is located on your device
addpath('C:\Users\...');
% Add the functions to the path
addpath('C:\Users\...');

% Provide the following parameters:
stimulation_start = 2100; % Window in which to detect calcium responses;
stimulation_stop = 2800; % Window in which to detect calcium responses
stimulation_length = 700;
ephys_min = -1.5; % minimum voltage of ephys traces (according to MES)
ephys_max = 1.5; % maximum voltage of ephys traces (according to MES)
AP_threshold = 0; % above which voltage is it considered an AP?
SD_threshold = 4.5; % From how many times the baseline SD is a peak considere a response?
size_subplots = 1; % Number of rows for subplots

%% Check if titles are right

D = pwd; % directory path
S = dir(fullfile(D,'F*UG*.tif*')); % get list of files in directory

% Ask the user whether the plot titles are correct
title_list = cell(1,length(S));
prompt = cell(1,length(S));
definput = cell(1,length(S));
j = 0;
for i = 1:length(S)
    j = j+1;
    my_title = ('Name your recording');
    title_list(j) = {my_title};
    prompt(j) = {['Title subplot ', num2str(j)]};
    definput(j) = title_list(j);
end
       
dlgtitle = 'Check the titles';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims,definput);

for i = 1:length(title_list)
    title_list(i) = answer(i)
end

%% Correct drift
disp('Correcting drift for every stack...')

% for GCaMP analysis we only use the UG file
D = pwd; % directory path
S = dir(fullfile(D,'F*UG*.tif*')); % get list of files in directory
[~,ndx] = natsortfiles({S.name}); % indices of correct order
S = S(ndx);

datetime=datetime;   % Create new folder for today's analysis
ymd = yyyymmdd(datetime);
mkdir('Analyses', [num2str(ymd)]);
ymd = num2str(ymd);

for i = 1:length(S)
    [my_rows, my_columns] = size(imread(S(i).name,1));
    fname = S(i).name;
    [pathstr,name,ext] = fileparts(fname);
    numimgs = size(imfinfo(fname),1);
    I_matrix = zeros(my_rows, my_columns);
    number_images(i) = numimgs;
    for j = 1:numimgs
        I = imread(fname,j);
        [rows,columns] = size(I);
        row_dif = rows - my_rows;
        if rows > my_rows
            I((rows - row_dif +1) : rows, :) = [];
        elseif rows < my_rows
            I((rows + 1) : (rows - row_dif),:) = zeros(1, columns);
        end
        col_dif = columns - my_columns;
        if columns > my_columns
            I(:, (columns - col_dif + 1) : columns) = [];
        elseif columns < my_columns
            I(:,(columns + 1) : (columns - col_dif)) = zeros(rows, 1);
        end
        I_matrix = cat(3,I_matrix,I);
    end
    I_matrix(:,:,1) = [];

    % Drift correction
    Y = single(I_matrix);
    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2));
    
    % Calculate drift in UR stack
    [M_final,shifts_g,template,options,col_shift] = normcorre(Y, options_rigid);

    cd (['Analyses\',ymd])
    
    for i = 1:size(M_final,3)
    cor_image = M_final(:,:,i);
    cor_image = uint16(cor_image);
    imwrite(cor_image, [name '_drift_cor.tif'], 'WriteMode','append');
    end
    cd .. ;
    cd .. ;
end

cd (['Analyses\',ymd])

%% Draw ROI of GCaMP (UG stack)
% You can draw a maximum of 5 mitos per image!

disp('Draw the first ROI on the mitochondrion, the second in the background. 7 seconds left...')

D = pwd; % directory path
S = dir(fullfile(D,'F*UG*_drift_cor.tif*')); % get list of files in directory
[~,ndx] = natsortfiles({S.name}); % indices of correct order
S = S(ndx);

i = 1;
j = 1;
mito_counter = 1;
while j <= length(S)
    UG_image_corr = S(j).name;
    
    % Draw ROI on first frame
    [X,map] = imread(UG_image_corr, 2);
    imshow(imadjust(X))
    
    % plot image 1 & draw a rectangular ROI on images
    subplot(2,1,1); imshow(imadjust(X));    % the contrast is adjusted to better visualize mitochondria
    title('Select ROI around mitochondria')
    roi = drawrectangle('LineWidth',2,'Color','white');
    subplot(2,1,2); imshow(imadjust(X));
    title('Move ROI to background location')
    roi2 = drawrectangle(gca,'Position',roi.Position);
    % set up listeners for ROI moving events
    addlistener(roi,'MovingROI',@(r1,evt) allevents(r1,evt,roi2));
    addlistener(roi,'ROIMoved',@(r1,evt) allevents(r1,evt,roi2));
    
    pause(5) % enter the amount of seconds you need between drawing the ROI on the mito and moving it to the background in the second image
    
    Stimulation(i).mito_number = mito_counter
    Stimulation(i).mask_mito = createMask(roi);
    Stimulation(i).mask_background = createMask(roi2);
    Stimulation(i).numimgs = number_images(j)
    
    saveas(figure(1), ['ROI ', num2str(i) '.png'])

    i = i+1;
    
    answer1 = questdlg('Do you want to draw another ROI in this stack?', "Mitos", "Yes", "No", "No");
    
    if answer1 == "Yes"
        j = j;
        mito_counter = mito_counter + 1;
    elseif answer1 == "No"
        j = j + 1;
        mito_counter = 1;
    end
           
end
%% Apply ROIs on stacks
disp('Creating masks from ROIs...')
disp('Applying masks on drift corrected stacks...')

D = pwd; % directory path
S = dir(fullfile(D,'F*UG*_drift_cor.tif*')); % get list of files in directory
[~,ndx] = natsortfiles({S.name}); % indices of correct order
S = S(ndx);

i = 1;
j = 0;
while i <= length(Stimulation)
    if Stimulation(i).mito_number == 1;
    j = j + 1;
    UG_image_corr = S(j).name;
    elseif Stimulation(i).mito_number ~= 1;
    UG_image_corr = S(j).name;
    end
    
    % Create mask from ROIs
    mask_mito = (Stimulation(i).mask_mito);
    mask_background = (Stimulation(i).mask_background);
    
    % Apply mito mask on UG stack
    Stimulation(i).UG_mean_values = [];
    
    for k = 1:size(imfinfo(UG_image_corr),1)
        I = imread(UG_image_corr,k);
        average_intensity = mean(I(mask_mito));
        Stimulation(i).UG_mean_values(k) = average_intensity;
    end
    
    % Apply background mask on UG stack
    Stimulation(i).UG_mean_background = [];
    
    for k = 1:size(imfinfo(UG_image_corr),1)
        I = imread(UG_image_corr,k);
        average_intensity = mean(I(mask_background));
        Stimulation(i).UG_mean_background(k) = average_intensity;
    end
    
    i = i + 1;
end

  
%% Automated start and stop times; and interval
% The time of the first frame as well as the interval between frames is
% calculated from the information in the metadata. This information is used
% to make an x-axis in seconds instead of frame numbers

cd ..
cd ..
D = pwd; % directory path
M = dir(fullfile(D,'*metadata.txt')); % get list of files in directory
[~,ndx] = natsortfiles({M.name}); % indices of correct order
M = M(ndx);

k = 0;
for i = 1:length(Stimulation)
    %search for start time
    if Stimulation(i).mito_number == 1
        k = k + 1;
        fid = fopen(M(k).name);
        ftell(fid);
        for j = 1:52
            fseek(fid, 1, 'eof');
            start_time = fgetl(fid);
            start_time = start_time(12:end);
            start_time = str2num(start_time);
        end
        Stimulation(i).start_time = start_time;
    elseif Stimulation(i).mito_number ~= 1
        Stimulation(i).start_time = start_time;
    end
    
    %search for sample interval
    if Stimulation(i).mito_number == 1
        fid = fopen(M(k).name);
        ftell(fid);
        for j = 1:48
            fseek(fid, 1, 'eof');
            sample_interval = fgetl(fid);
            sample_interval = sample_interval(10:end);
            sample_interval = str2num(sample_interval);
        end
        Stimulation(i).sample_interval = sample_interval;
    elseif Stimulation(i).mito_number ~= 1
        Stimulation(i).sample_interval = sample_interval;
    end
    
    % For every frame, determine its time (for x-axis of figures)
    if Stimulation(i).mito_number == 1
        x_inseconds = zeros(1, Stimulation(i).numimgs);
        x_inseconds(1) = start_time;
        for o = 2:length(x_inseconds)
            x_inseconds(o) = start_time + ((o-1) * sample_interval);
        end
        Stimulation(i).x_inseconds = x_inseconds;
    elseif Stimulation(i).mito_number ~= 1
        Stimulation(i).x_inseconds = x_inseconds;
    end
end



%% Adjust titles (add the mito number after the title name)
j = 0;
for i = 1:length(Stimulation)
    if Stimulation(i).mito_number == 1
        j = j+1;
        Stimulation(i).title = join(string([title_list(j), 'mito' Stimulation(i).mito_number]));
    elseif Stimulation(i).mito_number ~= 1
        Stimulation(i).title = join(string([title_list(j), 'mito', Stimulation(i).mito_number]));
    end
end

%% Substract the UG background values from the UG ROI mean values
disp('Substracting background values...')
for i = 1:length(Stimulation)
    Stimulation(i).UG_final = Stimulation(i).UG_mean_values - Stimulation(i).UG_mean_background;
end

%% Correct photobleaching until you are satisfied
disp('Correcting bleach...')

for i = 1:length(Stimulation)
    Stimulation(i).OnsetStimulus = round(stimulation_start/Stimulation(i).sample_interval);                                            
end
    
figure(2),
while 1
    % Ask the user to provide a value for 't' (= rate of decay)
    prompt = {'Enter value for t','Enter maximum value for t', 'Enter minimum value for t'};
    dlgtitle = 'Provide t value (speed of decay)';
    dims = [1 35];
    definput = {'50','t+40','t-40'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    t = str2double(answer{1});
    tmax = eval(answer{2});
    tmin = eval(answer{3});
    
    for i = 1:length(Stimulation)
        UG_data = Stimulation(i).UG_final;
        
        FitTo = nan(size(UG_data))';
        FitTo(1:Stimulation(i).OnsetStimulus) = UG_data(1:Stimulation(i).OnsetStimulus);                
        xfit = (1:1:length(FitTo))';
        idxValid = ~isnan(FitTo);
        
        a = UG_data(end) - UG_data(1);
        c = UG_data(end);
        
        ft = fittype('a * exp(-x / t) + c');
        fo = fitoptions('Method', 'NonlinearLeastSquares',...
        'Start', [a, c, t],...
        'Upper', [inf, inf, tmax],...          %if 'Exiting due to infeasibility' error, try changing upper-c to inf and lower-c to -inf; additionally min(UG_data) and 0 can be changed to inf and -inf
        'Lower', [-inf, -inf, tmin]);              % you can change this 0 to -inf
        
        BleachFit = fit(xfit(idxValid), FitTo(idxValid), ft, fo);
        BleachFitted = feval(BleachFit,xfit);                            % if you bleach fit to part of the trace, feval will create the bleach fit for the entire trace
        
        
        cd (['Analyses\',ymd])
        b = size(dir('*.png'),1);                    % counts the number of ROIs in the directory (png files)
        size_subplots =  ceil(b/3);                   % takes the number of ROIs, and calculates the required number of rows (with 3 subplots per row); the number is rounded up
        cd ..
        cd ..
        
        
        subplot(size_subplots, 3,i)
        hold on
        plot(UG_data, 'r')
        plot(FitTo, 'y')
        plot(BleachFitted, 'k')
        
        BleachFittedNorm = BleachFitted ./ BleachFitted(1,:);
        BleachFittedNormM = permute(BleachFittedNorm,[2 1]);      
        
        DataBC = UG_data ./ BleachFittedNormM;
        Stimulation(i).DataBC = DataBC;
        Stimulation(i).dFF_percent = (DataBC/(mean(DataBC(1:Stimulation(i).OnsetStimulus))) - 1)*200;         % dFF in percentage
        plot(Stimulation(i).dFF_percent, 'b')
        title(Stimulation(i).title, 'FontSize', 5);
        xlabel('Frame');
        ylabel('dF/F (%)');
    end
    
    pause(5)
    
    answer = questdlg('Do you accept this correction?', 'Bleach Correction', 'Yes', 'No, do it again','No, do it again');
    switch answer
        case 'Yes'                                          
        disp('Continuing with this correction!')
        break                                                                   % break --> the loop will stop and you will continue with the script
        case 'No, do it again'
        clf(figure(2))
    end
end

set(figure(2),'PaperOrientation','landscape');
set(figure(2),'PaperUnits','normalized');
set(figure(2),'PaperPosition', [0 0 1 1]);

%% dF/F using baseline averaging
disp('Calculating dF/F...')
% dF/F in percent is available under Stimulation(i).dFF_percent
% But, we just want dF/F (not in percentage)

% Calculate F0 as the average of the frames before stimulation (mean of baseline)
for i = 1:length(Stimulation)
    F0 = mean(Stimulation(i).DataBC(1:(round(stimulation_start/Stimulation(i).sample_interval)-1))); % until 1 frame before the stimuli start
        for j = 1:Stimulation(i).numimgs
        Stimulation(i).dFF(j) = ((Stimulation(i).DataBC(j) - F0) / F0);
        end
end

%% Signal smoothing filter
disp('Smoothing signal...')
% Average factor: 3
for i=1:length(Stimulation)
    Stimulation(i).dFF_smoothed(1:2) = Stimulation(i).dFF(1:2);
    for j = 3:Stimulation(i).numimgs
        Stimulation(i).dFF_smoothed(j) = mean(Stimulation(i).dFF(j-2:j));
    end
end

%% Rest data -> dFF smoothed of Prestim (1 until 1 frame before start stimulation)
disp('Calculating baseline parameters...')

for i = 1:length(Stimulation)
    Prestim(i).dFF_smoothed = Stimulation(i).dFF_smoothed(1:round(stimulation_start / Stimulation(i).sample_interval)-1);
end

baseline_SD = [];
% Prestim SD, average, peak, trendline, AUC
for i = 1:length(Prestim)
    Prestim(i).SD = std(Prestim(i).dFF_smoothed);
    baseline_SD = [baseline_SD, Prestim(i).SD];
    Prestim(i).Average = mean(Prestim(i).dFF_smoothed);
    Prestim(i).Peak = max(Prestim(i).dFF_smoothed);
    % Trendline
    f2 = fit(Stimulation(i).x_inseconds(1:(round(stimulation_start/Stimulation(i).sample_interval)-1))', Prestim(i).dFF_smoothed(1:(round(stimulation_start/Stimulation(i).sample_interval)-1))', 'poly1');
    Prestim(i).Trendline = f2.p1;
    % AUC
    for j = 1:length(Prestim(i).dFF_smoothed)
        if Prestim(i).dFF_smoothed(j) >= 0
            Prestim(i).AUC(j) = Prestim(i).dFF_smoothed(j);
        else
            Prestim(i).AUC(j) = 0;
        end
    end
    Prestim(i).AUC_total = trapz(1:(round(stimulation_start/Stimulation(i).sample_interval)-1), Prestim(i).AUC(1:(round(stimulation_start/Stimulation(i).sample_interval)-1)));
end

%% Prepare plotting
% Calculate minimum and maximum y-value to give all plots the same axes
ymin = inf;
ymax = -inf;
for i = 1:length(Stimulation)
    if min(Stimulation(i).dFF_smoothed) < ymin
        ymin = min(Stimulation(i).dFF_smoothed);
    end
    if max(Stimulation(i).dFF_smoothed) > ymax
        ymax = max(Stimulation(i).dFF_smoothed);
    end
end


%% Frequency APs during stimulation
disp('Calculating AP frequency...')
AP_frequency = [];
figure(12),

load('ephys_traces.mat')
ephys_traces = struct(ephys_traces);

j = 0
for i = 1:length(Stimulation)
    subplot(size_subplots,3,i);
    title(Stimulation(i).title, 'FontSize', 5); %if you have less subplots, you can increase the font size of titles
    ylabel('Ephys');
    % Peaks only during stimulation
    if Stimulation(i).mito_number == 1
        j = j+1
    end
    freq_ephys_sampling = length(ephys_traces(j).y) / Stimulation(i).x_inseconds(end);
    x = ephys_traces(j).y(2000 * freq_ephys_sampling : 3000 * freq_ephys_sampling);
    findpeaks(x, 'MinPeakHeight',AP_threshold);
    ylim([ephys_min ephys_max]);
    [pks, locs] =findpeaks(x, 'MinPeakHeight',AP_threshold);
    APfreq = length(locs) / (stimulation_length / 1000); % frequency of APs in AP/sec
    Stimulation(i).num_APs = length(pks);
    Stimulation(i).APfreq = APfreq;
    AP_frequency = [AP_frequency, APfreq];
end

set(figure(12),'PaperOrientation','landscape');
set(figure(12),'PaperUnits','normalized');
set(figure(12),'PaperPosition', [0 0 1 1]);

% If MinPeakHeight does not work properly, try MinPeakProminence. Change AP
% threshold to for example 1 and play around with this number.

%% Find peak of Calcium response
disp('Finding peaks...')
peak_values = [];
peak_loc = [];
time_to_peak = [];

figure(9),
for i = 1:length(Stimulation)
    subplot(size_subplots,3,i);
    plot(Stimulation(i).x_inseconds, Stimulation(i).dFF_smoothed);
    title(Stimulation(i).title, 'FontSize', 5);
    if i == length(Stimulation)
    xlabel('Time (ms)');
    end
    xlim([0 Stimulation(i).x_inseconds(end)]);
    if i == 1
    ylabel('dF/F');
    end
    ylim([ymin ymax]);
    hold on;
    x = [stimulation_start, stimulation_start, stimulation_stop, stimulation_stop];
    y = [ymin, ymax, ymax, ymin];
    a = fill(x,y,'blue');
    a.FaceAlpha = 0.1;
    set(a, 'EdgeColor', 'none');
    [max_val,idx_val]=max(Stimulation(i).dFF_smoothed(round(stimulation_start/Stimulation(i).sample_interval):(stimulation_stop + 200)/Stimulation(i).sample_interval));
    Stimulation(i).peak = max_val;
    Stimulation(i).peak_loc = idx_val;
    peak_values = [peak_values, max_val];
    peak_loc = [peak_loc, idx_val];
    time_to_peak = [time_to_peak, Stimulation(i).x_inseconds(idx_val + round(stimulation_start / Stimulation(i).sample_interval)-1) - stimulation_start];
    Stimulation(i).time_to_peak = time_to_peak(i);
    if Stimulation(i).peak > SD_threshold * Prestim(i).SD
        plot(Stimulation(i).x_inseconds(idx_val + round(stimulation_start / Stimulation(i).sample_interval)-1), max_val, 'og');
        Stimulation(i).response = 'yes';
    else
        plot(Stimulation(i).x_inseconds(idx_val + round(stimulation_start / Stimulation(i).sample_interval)-1), max_val, 'or');
        Stimulation(i).response = 'no response';
    end
    hold off;
end
sgtitle('Peak Detection, green = response, red = no response')

set(figure(9),'PaperOrientation','landscape');
set(figure(9),'PaperUnits','normalized');
set(figure(9),'PaperPosition', [0 0 1 1]);


%% Plot Ephys trace over dF/F
disp('Plotting Ephys traces over dF/F...')

figure(4),
j = 0;
for i = 1:length(Stimulation)
    if Stimulation(i).mito_number == 1
        j = j + 1;
    end
    subplot(size_subplots,3,i);
    title(Stimulation(i).title, 'FontSize', 5);
    x1 = 1:length(ephys_traces(j).y);
    y1 = ephys_traces(j).y;
    l = line(x1,y1);
    l.Color = [0 0 1 0.05];
    ax1 = gca;
    ax1.XLim = [0 length(ephys_traces(j).y)];
    ax1.XColor = 'k';
    ax1.XAxis.Visible = 'off';
    ax1.YAxis.Visible = 'off';
    ylim([ephys_min ephys_max])
    ax1_pos = ax1.Position;
    ax2 = axes('Position', ax1_pos, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left','Color', 'none');
    ax2.XColor = 'k';
    xlim([0 Stimulation(i).x_inseconds(end)]);
    ylim([ymin ymax]);
    xticks([2000 4000]);
    if i == 1
    ylabel('dF/F');
    end
    if i == length(Stimulation)
    xlabel('Time (ms)');
    end
    x2 = Stimulation(i).x_inseconds;
    y2 = Stimulation(i).dFF_smoothed;
    line(x2, y2, 'Parent', ax2);
    l2 = line(x2,y2);
    if Stimulation(i).peak > SD_threshold * Prestim(i).SD
    l2.Color = [0 0 0];
    else
    l2.Color = [0.5 0.5 0.5];
    end
%     xline(stimulation_start, '--') Use these if you are not sure whether the ephys traces are aligned well
%     xline(stimulation_stop, '--')
end

set(figure(4),'PaperOrientation','landscape');
set(figure(4),'PaperUnits','normalized');
set(figure(4),'PaperPosition', [0 0 1 1]);

%% Area under curve (without correction) -> vanaf stimulatie
disp('Calculating area under curve...')
for i = 1:length(Stimulation)
    for j = 1:Stimulation(i).numimgs
        if Stimulation(i).dFF_smoothed(j) >= 0
            Stimulation(i).AUC(j) = Stimulation(i).dFF_smoothed(j);
        else
            Stimulation(i).AUC(j) = 0;
        end
    end
    Stimulation(i).AUC_total = trapz(round(stimulation_start/Stimulation(i).sample_interval):round((stimulation_stop + 200)/Stimulation(i).sample_interval), Stimulation(i).AUC(round(stimulation_start/Stimulation(i).sample_interval):round((stimulation_stop + 200)/Stimulation(i).sample_interval)));
end

disp('Plotting AUC...')
% visualisation corrected AUC
figure(5),
for k = 1:length(Stimulation)
    subplot(size_subplots,3,k);
    plot(Stimulation(k).dFF_smoothed)
    title(Stimulation(k).title, 'FontSize', 5);
    ylim([ymin ymax]);
    xline(round(stimulation_start/Stimulation(k).sample_interval));
    hold on
    x = [round(stimulation_start/Stimulation(k).sample_interval): round((stimulation_stop + 200)/Stimulation(k).sample_interval)];
    Stimulation(k).AUC_vis = Stimulation(k).AUC;
    y = Stimulation(k).AUC_vis(round(stimulation_start/Stimulation(k).sample_interval): round((stimulation_stop + 200)/Stimulation(k).sample_interval));
    yline(0)
    if k == 1
    ylabel('dF/F');
    end
    if Stimulation(k).peak > SD_threshold * Prestim(k).SD
        a = area(x,y,'FaceColor', 'g', 'EdgeColor', 'none');
        a.FaceAlpha = 0.1;
    else
        a = area(x,y,'FaceColor', 'r', 'EdgeColor', 'none');
        a.FaceAlpha = 0.1;
    end
    hold off;
end
sgtitle('AUC, green = response, red = no response')

set(figure(5),'PaperOrientation','landscape');
set(figure(5),'PaperUnits','normalized');
set(figure(5),'PaperPosition', [0 0 1 1]);



%% Plot sorted APfreq over peak height
% if the cell responded, the trial is indicated by a black square
% if no response, the trial is indicated by a grey circle

peak_values = [];

for i = 1:length(Stimulation)
    peak_values(i) = Stimulation(i).peak;
end

[APfrequency_sorted, idx] = sort(AP_frequency);
peak_height_sorted = peak_values(idx);
prestim_SD_sorted = baseline_SD(idx);

bin1 = 0:89;
bin2 = 90:110;
bin3 = 111:130;
bin4 = 131:150;
bin5 = 151:400;

bin1_data = [];
bin2_data = [];
bin3_data = [];
bin4_data = [];
bin5_data = [];

index_array = cell(1, length(peak_height_sorted));
for i = 1:length(peak_height_sorted)
    if ismember(round(APfrequency_sorted(i)), bin1)
        bin1_data = [bin1_data, peak_height_sorted(i)];
        index_array{i} = '< 90';
    elseif ismember(round(APfrequency_sorted(i)), bin2)
         bin2_data = [bin2_data, peak_height_sorted(i)];
        index_array{i} = '91-110';
    elseif ismember(round(APfrequency_sorted(i)), bin3)
         bin3_data = [bin3_data, peak_height_sorted(i)];
        index_array{i} = '111-130';
    elseif ismember(round(APfrequency_sorted(i)), bin4)
         bin4_data = [bin4_data, peak_height_sorted(i)];
        index_array{i} = '131-150';
    elseif ismember(round(APfrequency_sorted(i)), bin5)
         bin5_data = [bin5_data, peak_height_sorted(i)];
        index_array{i} = '> 150';
    end
end

titles = {'< 90', '91-110', '111-130', '131-150', '> 150'};
x = categorical(titles);
x = reordercats(x, titles);
y = [mean(bin1_data), mean(bin2_data), mean(bin3_data), mean(bin4_data), mean(bin5_data)];

figure(10),
bar(x,y)
hold on
for i = 1:length(peak_height_sorted)
    x_var = categorical(index_array(i));
    y_var = peak_height_sorted(i);
    if y_var > SD_threshold * prestim_SD_sorted(i)
        plot(x_var, y_var, 's', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
    else
        plot(x_var, y_var, 'o', 'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerEdgeColor', [0.6 0.6 0.6]);
    end
end
title('Peak height per AP frequency');
xlabel('AP frequency');
ylabel('Peak height');

set(figure(10),'PaperOrientation','landscape');
set(figure(10),'PaperUnits','normalized');
set(figure(10),'PaperPosition', [0 0 1 1]);

%% Plot sorted APfreq over time to peak

peak_loc = [];
for i = 1:length(Stimulation)
    peak_loc(i) = Stimulation(i).x_inseconds(Stimulation(i).peak_loc + round(stimulation_start / Stimulation(i).sample_interval)-1);
end

time_to_peak_sorted = time_to_peak(idx);

bin1_data = [];
bin2_data = [];
bin3_data = [];
bin4_data = [];
bin5_data = [];

index_array = cell(1, length(time_to_peak_sorted));
for i = 1:length(time_to_peak_sorted)
    if ismember(round(APfrequency_sorted(i)), bin1)
        bin1_data = [bin1_data, time_to_peak_sorted(i)];
        index_array{i} = '< 90';
    elseif ismember(round(APfrequency_sorted(i)), bin2)
         bin2_data = [bin2_data, time_to_peak_sorted(i)];
        index_array{i} = '91-110';
    elseif ismember(round(APfrequency_sorted(i)), bin3)
         bin3_data = [bin3_data, time_to_peak_sorted(i)];
        index_array{i} = '111-130';
    elseif ismember(round(APfrequency_sorted(i)), bin4)
         bin4_data = [bin4_data, time_to_peak_sorted(i)];
        index_array{i} = '131-150';
    elseif ismember(round(APfrequency_sorted(i)), bin5)
         bin5_data = [bin5_data, time_to_peak_sorted(i)];
        index_array{i} = '> 150';
    end
end

titles = {'< 90', '91-110', '111-130', '131-150', '> 150'};
x = categorical(titles);
x = reordercats(x, titles);
y = [mean(bin1_data), mean(bin2_data), mean(bin3_data), mean(bin4_data), mean(bin5_data)];

figure(11),
bar(x,y)
hold on
for i = 1:length(time_to_peak_sorted)
    x_var = categorical(index_array(i));
    y_var2 = time_to_peak_sorted(i);
    y_var = peak_height_sorted(i);
    if y_var > SD_threshold * prestim_SD_sorted(i)
        plot(x_var, y_var2, 's', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
    else
        plot(x_var, y_var2, 'o', 'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerEdgeColor', [0.6 0.6 0.6]);
    end
end
yline(stimulation_stop - stimulation_start, '-')
title('Peak location per AP frequency');
xlabel('AP frequency');
ylabel('Time to Peak');

set(figure(11),'PaperOrientation','landscape');
set(figure(11),'PaperUnits','normalized');
set(figure(11),'PaperPosition', [0 0 1 1]);

%% Plot sorted APfreq over AUC
AUC_values = zeros(size(Stimulation));

for i = 1:length(Stimulation)
    AUC_values(i) = Stimulation(i).AUC_total;
end

AUC_sorted = AUC_values(idx);

bin1_data = [];
bin2_data = [];
bin3_data = [];
bin4_data = [];
bin5_data = [];

index_array = cell(1, length(AUC_sorted));
for i = 1:length(AUC_sorted)
    if ismember(round(APfrequency_sorted(i)), bin1)
        bin1_data = [bin1_data, AUC_sorted(i)];
        index_array{i} = '< 90';
    elseif ismember(round(APfrequency_sorted(i)), bin2)
         bin2_data = [bin2_data, AUC_sorted(i)];
        index_array{i} = '91-110';
    elseif ismember(round(APfrequency_sorted(i)), bin3)
         bin3_data = [bin3_data, AUC_sorted(i)];
        index_array{i} = '111-130';
    elseif ismember(round(APfrequency_sorted(i)), bin4)
         bin4_data = [bin4_data, AUC_sorted(i)];
        index_array{i} = '131-150';
    elseif ismember(round(APfrequency_sorted(i)), bin5)
         bin5_data = [bin5_data, AUC_sorted(i)];
        index_array{i} = '> 150';
    end
end

titles = {'< 90', '91-110', '111-130', '131-150', '> 150'};
x = categorical(titles);
x = reordercats(x, titles);
y = [mean(bin1_data), mean(bin2_data), mean(bin3_data), mean(bin4_data), mean(bin5_data)];

figure(6),
bar(x,y)
hold on
for i = 1:length(AUC_values)
    x_var = categorical(index_array(i));
    y_var2 = AUC_sorted(i);
    y_var = peak_height_sorted(i);
    if y_var > SD_threshold * prestim_SD_sorted(i)
        plot(x_var, y_var2, 's', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
    else
        plot(x_var, y_var2, 'o', 'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerEdgeColor', [0.6 0.6 0.6]);
    end
end

title('AUC per AP frequency');
xlabel('AP frequency');
ylabel('AUC');

set(figure(6),'PaperOrientation','landscape');
set(figure(6),'PaperUnits','normalized');
set(figure(6),'PaperPosition', [0 0 1 1]);

%% Quantify calcium response to stimulus using slope trendline
disp('Drawing trendline...')
trendline_values = [];
figure(8),
for i = 1:length(Stimulation)
    subplot(size_subplots,3,i);
    if Stimulation(i).peak > SD_threshold * Prestim(i).SD
        plot(Stimulation(i).x_inseconds, Stimulation(i).dFF_smoothed, 'g');
    else
        plot(Stimulation(i).x_inseconds, Stimulation(i).dFF_smoothed, 'r');
    end
    
    f = fit(Stimulation(i).x_inseconds(round(stimulation_start/Stimulation(i).sample_interval):round((stimulation_start + (stimulation_stop - stimulation_start))/Stimulation(i).sample_interval))', Stimulation(i).dFF_smoothed(round(stimulation_start/Stimulation(i).sample_interval):round((stimulation_start + (stimulation_stop - stimulation_start))/Stimulation(i).sample_interval))', 'poly1');
    Stimulation(i).fit_response = f.p1;
    hold on;
    plot(f);
    title(Stimulation(i).title, 'FontSize', 5);
    if i == length(Stimulation)
    xlabel('Time (ms)');
    end
    xlim([0 Stimulation(i).x_inseconds(end)]);
    xticks([2000 4000]);
    format short;
    if i == 1
    ylabel('dF/F');
    end
    ylim([ymin ymax]);
    x = [stimulation_start, stimulation_start, stimulation_stop, stimulation_stop];
    y = [ymin, ymax, ymax, ymin];
    a = fill(x,y,'blue');
    a.FaceAlpha = 0.1;
    set(a, 'EdgeColor', 'none');
    legend('hide');
    hold off;
    trendline_values = [trendline_values, f.p1];
end
sgtitle('Trendline, green = response, red = no response')

set(figure(7),'PaperOrientation','landscape');
set(figure(7),'PaperUnits','normalized');
set(figure(7),'PaperPosition', [0 0 1 1]);

%% Plot sorted APfreq over slope trendline
trendline_values = zeros(size(Stimulation));

for i = 1:length(Stimulation)
    trendline_values(i) = Stimulation(i).fit_response;
end

trendline_sorted = trendline_values(idx);

bin1_data = [];
bin2_data = [];
bin3_data = [];
bin4_data = [];
bin5_data = [];

index_array = cell(1, length(trendline_sorted));
for i = 1:length(trendline_sorted)
    if ismember(round(APfrequency_sorted(i)), bin1)
        bin1_data = [bin1_data, trendline_sorted(i)];
        index_array{i} = '< 90';
    elseif ismember(round(APfrequency_sorted(i)), bin2)
         bin2_data = [bin2_data, trendline_sorted(i)];
        index_array{i} = '91-110';
    elseif ismember(round(APfrequency_sorted(i)), bin3)
         bin3_data = [bin3_data, trendline_sorted(i)];
        index_array{i} = '111-130';
    elseif ismember(round(APfrequency_sorted(i)), bin4)
         bin4_data = [bin4_data, trendline_sorted(i)];
        index_array{i} = '131-150';
    elseif ismember(round(APfrequency_sorted(i)), bin5)
         bin5_data = [bin5_data, trendline_sorted(i)];
        index_array{i} = '> 150';
    end
end

titles = {'< 90', '91-110', '111-130', '131-150', '> 150'};
x = categorical(titles);
x = reordercats(x, titles);
y = [mean(bin1_data), mean(bin2_data), mean(bin3_data), mean(bin4_data), mean(bin5_data)];

figure(8),
bar(x,y)
hold on
for i = 1:length(trendline_sorted)
    x_var = categorical(index_array(i));
    y_var2 = trendline_sorted(i);
    y_var = peak_height_sorted(i);
    if y_var > SD_threshold * prestim_SD_sorted(i)
        plot(x_var, y_var2, 's', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
    else
        plot(x_var, y_var2, 'o', 'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerEdgeColor', [0.6 0.6 0.6]);
    end
end
title('Trendline per AP frequency');
xlabel('AP frequency');
ylabel('Slope trendline');

set(figure(8),'PaperOrientation','landscape');
set(figure(8),'PaperUnits','normalized');
set(figure(8),'PaperPosition', [0 0 1 1]);

%% Repetitions (calculate the average of repetitions)

title_list_cell = cell(1,1);
for i = 1:length(Stimulation)
    title_list_cell(i) = cellstr(Stimulation(i).title);
end

Stimulation(1).AvgAUC = "";
Stimulation(1).AvgTrendline = "";
Stimulation(1).AvgPeakHeight = "";
Stimulation(1).AvgTimeToPeak = "";
Stimulation(1).Repetition = "No";

Prestim(1).AvgDFF = "";
Prestim(1).AvgTrendline = "";
Prestim(1).AvgPeakHeight = "";
Prestim(1).AvgAUC = "";
Prestim(1).Repetition = "No";

double_list = cell(1,1);
double_list{1} = char(Stimulation(1).title);
for i = 2:length(Stimulation)
    if Stimulation(i).mito_number == 1
        matches = strfind(double_list(~cellfun('isempty',double_list)), title_list_cell{i});
        if any(horzcat(matches{:}))
            Stimulation(i).Repetition = "Yes";
            Prestim(i).Repetition = "Yes";
            % Calculate average of the repetition
            if i ==2 & Stimulation(2).mito_number == 1 
                row_shift = 1;
            elseif i <= 2
                temp_array = [];
                temp_array = [temp_array, Stimulation(1:i+2).mito_number];
                row_shift = max(temp_array);
            elseif i > 2 & i <= length(Stimulation) - 2
                temp_array = [];
                temp_array = [temp_array, Stimulation(i-2:i+2).mito_number];
                row_shift = max(temp_array);
            elseif i >= length(Stimulation) - 2
                temp_array = [];
                temp_array = [temp_array, Stimulation(i-2:end).mito_number];
                row_shift = max(temp_array);
            end
            Stimulation(i).AvgAUC = (Stimulation(i).AUC_total + Stimulation(i-row_shift).AUC_total) / 2;
            Stimulation(i).AvgTrendline = (Stimulation(i).fit_response + Stimulation(i-row_shift).fit_response) / 2;
            Stimulation(i).AvgPeakHeight = (Stimulation(i).peak + Stimulation(i-row_shift).peak) / 2;
            Stimulation(i).AvgTimeToPeak = (Stimulation(i).time_to_peak + Stimulation(i-row_shift).time_to_peak) / 2;
            Prestim(i).AvgDFF = (Prestim(i).Average + Prestim(i-row_shift).Average) / 2;
            Prestim(i).AvgTrendline = (Prestim(i).Trendline + Prestim(i-row_shift).Trendline) / 2;
            Prestim(i).AvgPeakHeight = (Prestim(i).Peak + Prestim(i-row_shift).Peak) / 2;
            Prestim(i).AvgAUC = (Prestim(i).AUC_total + Prestim(i-row_shift).AUC_total) / 2;
            
        else
            Stimulation(i).Repetition = "No";
            Stimulation(i).AvgAUC = "";
            Stimulation(i).AvgTrendline = "";
            Stimulation(i).AvgPeakHeight = "";
            Stimulation(i).AvgTimeToPeak = "";
            Prestim(i).Repetition = "No";
            Prestim(i).AvgDFF = "";
            Prestim(i).AvgTrendline = "";
            Prestim(i).AvgPeakHeight = "";
            Prestim(i).AvgAUC = "";
        end
    elseif Stimulation(i).mito_number ~= 1
        Stimulation(i).Repetition = Stimulation(i-Stimulation(i).mito_number+1).Repetition;
        if Stimulation(i).Repetition == "Yes"
            Prestim(i).Repetition = "Yes";
            Stimulation(i).AvgAUC = (Stimulation(i).AUC_total + Stimulation(i-row_shift).AUC_total) / 2;
            Stimulation(i).AvgTrendline = (Stimulation(i).fit_response + Stimulation(i-row_shift).fit_response) / 2;
            Stimulation(i).AvgPeakHeight = (Stimulation(i).peak + Stimulation(i-row_shift).peak) / 2;
            Stimulation(i).AvgTimeToPeak = (Stimulation(i).time_to_peak + Stimulation(i-row_shift).time_to_peak) / 2;
            Prestim(i).AvgDFF = (Prestim(i).Average + Prestim(i-row_shift).Average) / 2;
            Prestim(i).AvgTrendline = (Prestim(i).Trendline + Prestim(i-row_shift).Trendline) / 2;
            Prestim(i).AvgPeakHeight = (Prestim(i).Peak + Prestim(i-row_shift).Peak) / 2;
            Prestim(i).AvgAUC = (Prestim(i).AUC_total + Prestim(i-row_shift).AUC_total) / 2;
        else
            Stimulation(i).Repetition = "No";
            Stimulation(i).AvgAUC = "";
            Stimulation(i).AvgTrendline = "";
            Stimulation(i).AvgPeakHeight = "";
            Stimulation(i).AvgTimeToPeak = "";
            Prestim(i).Repetition = "No";
            Prestim(i).AvgDFF = "";
            Prestim(i).AvgTrendline = "";
            Prestim(i).AvgPeakHeight = "";
            Prestim(i).AvgAUC = "";
        end
        
    else
        Stimulation(i).Repetition = "No";
        Stimulation(i).AvgAUC = "";
        Stimulation(i).AvgTrendline = "";
        Stimulation(i).AvgPeakHeight = "";
        Stimulation(i).AvgTimeToPeak = "";
        Prestim(i).Repetition = "No";
        Prestim(i).AvgDFF = "";
        Prestim(i).AvgTrendline = "";
        Prestim(i).AvgPeakHeight = "";
        Prestim(i).AvgAUC = "";
    end
    double_list{i} = char(Stimulation(i).title);
end


%% Repetitions: trendline over stimulation strength

dif_titles = unique(title_list_cell,'stable');
count_titles = cellfun(@(x) sum(ismember(title_list_cell,x)),dif_titles,'un',0);

x3 = cell(1,length(Stimulation));
y3 = [];
for i = 1:length(Stimulation)
    x3{i} = convertStringsToChars(Stimulation(i).title);
    y3 = [y3, Stimulation(i).fit_response];
end

y4 = [];
for i = 1:length(dif_titles)
    y4 = [y4, mean(y3(strcmp(dif_titles{i}, x3)))];
end

x4 = categorical(dif_titles);
x4 = reordercats(x4, dif_titles);

figure(13),
bar(x4, y4, 'red');
hold on
for i = 1:length(Stimulation)
    x_var = categorical(cellstr(Stimulation(i).title));
    y_var = Stimulation(i).fit_response;
    if Stimulation(i).peak > SD_threshold * Prestim(i).SD
        plot(x_var, y_var, 's', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
    else
        plot(x_var, y_var, 'o', 'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerEdgeColor', [0.6 0.6 0.6]);
    end
end
title('Averaged trendline per stimulation strength');
xlabel('Stimulation strength (pA)');
ylabel('Slope trendline');

set(figure(13),'PaperOrientation','landscape');
set(figure(13),'PaperUnits','normalized');
set(figure(13),'PaperPosition', [0 0 1 1]);


%% Repetitions: AUC over stimulation strength

y5 = [];
for i = 1:length(Stimulation)
    y5 = [y5, Stimulation(i).AUC_total];
end

y6 = [];
for i = 1:length(dif_titles)
    y6 = [y6, mean(y5(strcmp(dif_titles{i}, x3)))];
end

figure(14),
bar(x4, y6, 'red');
hold on
for i = 1:length(Stimulation)
    x_var = categorical(cellstr(Stimulation(i).title));
    y_var = Stimulation(i).AUC_total;
    if Stimulation(i).peak > SD_threshold * Prestim(i).SD
        plot(x_var, y_var, 's', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
    else
        plot(x_var, y_var, 'o', 'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerEdgeColor', [0.6 0.6 0.6]);
    end
end
hold off
title('Averaged AUC per stimulation strength');
xlabel('Stimulation strength');
ylabel('AUC');

set(figure(14),'PaperOrientation','landscape');
set(figure(14),'PaperUnits','normalized');
set(figure(14),'PaperPosition', [0 0 1 1]);


%% Repetitions: peak over stimulation strength

y7 = [];
for i = 1:length(Stimulation)
    y7 = [y7, Stimulation(i).peak];
end

y8 = [];
for i = 1:length(dif_titles)
    y8 = [y8, mean(y7(strcmp(dif_titles{i}, x3)))];
end

figure(15),
bar(x4, y8, 'red');
hold on
for i = 1:length(Stimulation)
    x_var = categorical(cellstr(Stimulation(i).title));
    y_var = Stimulation(i).peak;
    if Stimulation(i).peak > SD_threshold * Prestim(i).SD
        plot(x_var, y_var, 's', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
    else
        plot(x_var, y_var, 'o', 'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerEdgeColor', [0.6 0.6 0.6]);
    end
end
hold off
title('Averaged peak height per stimulation strength');
xlabel('Stimulation strength');
ylabel('Peak Height');

set(figure(15),'PaperOrientation','landscape');
set(figure(15),'PaperUnits','normalized');
set(figure(15),'PaperPosition', [0 0 1 1]);


%% Repetitions: time to peak over stimulation strength

y9 = [];
for i = 1:length(Stimulation)
    y9 = [y9, time_to_peak(i)];
end

y10 = [];
for i = 1:length(dif_titles)
    y10 = [y10, mean(y9(strcmp(dif_titles{i}, x3)))];
end

figure(16),
bar(x4, y10, 'red');
hold on
for i = 1:length(Stimulation)
    x_var = categorical(cellstr(Stimulation(i).title));
    y_var = Stimulation(i).time_to_peak;
    if Stimulation(i).peak > SD_threshold * Prestim(i).SD
        plot(x_var, y_var, 's', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
    else
        plot(x_var, y_var, 'o', 'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerEdgeColor', [0.6 0.6 0.6]);
    end
end
hold off
title('Average time to peak per stimulation strength (lower is quicker peak; y=0 is start stimulation)');
xlabel('Stimulation strength');
ylabel('Time to Peak (ms)');
yline(stimulation_stop - stimulation_start)

set(figure(16),'PaperOrientation','landscape');
set(figure(16),'PaperUnits','normalized');
set(figure(16),'PaperPosition', [0 0 1 1]);


%% Divide data per subcellular region

% Per subcellular region, calculate the mean parameter of 'first trials'
% only (so, no repetitions)

% Divide data into dendrite, proximal axon, distal axon
Dendrite_titles = cell(1,1);
Dendrite_AUC = [];
Dendrite_PeakHeight = [];
Dendrite_trendline = [];
Dendrite_PeakLoc = [];

ProxAx_titles = cell(1,1);
ProxAx_AUC = [];
ProxAx_PeakHeight = [];
ProxAx_trendline = [];
ProxAx_PeakLoc = [];

DistAx_titles = cell(1,1);
DistAx_AUC = [];
DistAx_PeakHeight = [];
DistAx_trendline = [];
DistAx_PeakLoc = [];

Soma_titles = cell(1,1);
Soma_AUC = [];
Soma_PeakHeight = [];
Soma_trendline = [];
Soma_PeakLoc = [];

% If a stimulation is not named correctly, it will end up in 'other'(as a
% check).
Other_titles = cell(1,1);
Other_AUC = [];
Other_PeakHeight = [];
Other_trendline = [];
Other_PeakLoc = [];


% stat_titles = cell(1,1); % Enable if doing within cell statistics

for i = 1:length(Stimulation)
    if Stimulation(i).Repetition == "No"
        if ismember('endrite', convertStringsToChars(Stimulation(i).title))
            Dendrite_titles{i} = Stimulation(i).title;
            Dendrite_AUC = [Dendrite_AUC, Stimulation(i).AUC_total];
            Dendrite_PeakHeight = [Dendrite_PeakHeight, Stimulation(i).peak];
            Dendrite_trendline = [Dendrite_trendline, Stimulation(i).fit_response];
            Dendrite_PeakLoc = [Dendrite_PeakLoc, Stimulation(i).x_inseconds(Stimulation(i).peak_loc + round(stimulation_start / Stimulation(i).sample_interval)-1)];
            
        elseif ismember('Soma', convertStringsToChars(Stimulation(i).title))
            Soma_titles{i} = Stimulation(i).title;
            Soma_AUC = [Soma_AUC, Stimulation(i).AUC_total];
            Soma_PeakHeight = [Soma_PeakHeight, Stimulation(i).peak];
            Soma_trendline = [Soma_trendline, Stimulation(i).fit_response];
            Soma_PeakLoc = [Soma_PeakLoc, Stimulation(i).x_inseconds(Stimulation(i).peak_loc + round(stimulation_start / Stimulation(i).sample_interval)-1)];
            
        elseif ismember('unknown', convertStringsToChars(Stimulation(i).title)) | ismember('passant', convertStringsToChars(Stimulation(i).title))
            DistAx_titles{i} = Stimulation(i).title;
            DistAx_AUC = [DistAx_AUC, Stimulation(i).AUC_total];
            DistAx_PeakHeight = [DistAx_PeakHeight, Stimulation(i).peak];
            DistAx_trendline = [DistAx_trendline, Stimulation(i).fit_response];
            DistAx_PeakLoc = [DistAx_PeakLoc, Stimulation(i).x_inseconds(Stimulation(i).peak_loc + round(stimulation_start / Stimulation(i).sample_interval)-1)];
            
        elseif ~isnan(str2double(extractBetween(Stimulation(i).title, 8,8)))
            ProxAx_titles{i} = Stimulation(i).title;
            ProxAx_AUC = [ProxAx_AUC, Stimulation(i).AUC_total];
            ProxAx_PeakHeight = [ProxAx_PeakHeight, Stimulation(i).peak];
            ProxAx_trendline = [ProxAx_trendline, Stimulation(i).fit_response];
            ProxAx_PeakLoc = [ProxAx_PeakLoc, Stimulation(i).x_inseconds(Stimulation(i).peak_loc + round(stimulation_start / Stimulation(i).sample_interval)-1)];
            
        else
            Other_titles{i} = Stimulation(i).title;
            Other_AUC = [Other_AUC, Stimulation(i).AUC_total];
            Other_PeakHeight = [Other_PeakHeight, Stimulation(i).peak];
            Other_trendline = [Other_trendline, Stimulation(i).fit_response];
            Other_PeakLoc = [Other_PeakLoc, Stimulation(i).x_inseconds(Stimulation(i).peak_loc + round(stimulation_start / Stimulation(i).sample_interval)-1)];
            
        end
    end
end

DenAvgAUC = mean(Dendrite_AUC);
DenAvgTrendline = mean(Dendrite_trendline);
DenAvgPeakHeight = mean(Dendrite_PeakHeight);
DenAvgPeakLoc = mean(Dendrite_PeakLoc);

ProxAxAvgAUC = mean(ProxAx_AUC);
ProxAxAvgTrendline = mean(ProxAx_trendline);
ProxAxAvgPeakHeight = mean(ProxAx_PeakHeight);
ProxAxAvgPeakLoc = mean(ProxAx_PeakLoc);

DistAxAvgAUC = mean(DistAx_AUC);
DistAxAvgTrendline = mean(DistAx_trendline);
DistAxAvgPeakHeight = mean(DistAx_PeakHeight);
DistAxAvgPeakLoc = mean(DistAx_PeakLoc);

SomaAvgAUC = mean(Soma_AUC);
SomaAvgTrendline = mean(Soma_trendline);
SomaAvgPeakHeight = mean(Soma_PeakHeight);
SomaAvgPeakLoc = mean(Soma_PeakLoc);

T6 = [SomaAvgAUC, SomaAvgTrendline, SomaAvgPeakHeight, SomaAvgPeakLoc;...
    DenAvgAUC, DenAvgTrendline, DenAvgPeakHeight, DenAvgPeakLoc;...
    ProxAxAvgAUC, ProxAxAvgTrendline, ProxAxAvgPeakHeight, ProxAxAvgPeakLoc;...
    DistAxAvgAUC, DistAxAvgTrendline, DistAxAvgPeakHeight, DistAxAvgPeakLoc];
T6 = array2table(T6, 'Rownames', {'Soma', 'Dendrite', 'Proximal Axon', 'Distal Axon'}, 'VariableNames', {'AUC', 'Slope_Trendline', 'Peak_Height', 'Peak_Location'});

%% Save data
disp('Saving data...')

cd(['Analyses\',ymd]);

print(figure(2),'-dpdf','PDF Bleach Correction.pdf')

% find largest array of dF/F values to preallocate the table
max_length = 0;
for i = 1:length(Stimulation)
    length_dff = length(Stimulation(i).dFF_smoothed);
    if length_dff > max_length
        max_length = length_dff;
    end
end
dFF_values = NaN(length(Stimulation), max_length);

stim_response = cell(1,1);
AvgAUC = [];
AvgTrendline = [];
AvgPeakHeight = [];
AvgTimeToPeak = [];
Repetition = cell(1,1);
PreAvgDFF = [];
PreAvgTrendline = [];
PreAvgPeakHeight = [];
PreAvgAUC = [];
PreRepetition = cell(1,1);
% only save the final dFF values because of the excel limitations
for i = 1:length(Stimulation)
    for j = 1:Stimulation(i).numimgs
        dFF_values(i,j) = Stimulation(i).dFF_smoothed(j);
    end
    stim_response{i} = Stimulation(i).response;
    imaging_freq(i) = 1000 / Stimulation(i).sample_interval;
    AvgAUC = [AvgAUC, Stimulation(i).AvgAUC];
    AvgTrendline = [AvgTrendline, Stimulation(i).AvgTrendline];
    AvgPeakHeight = [AvgPeakHeight, Stimulation(i).AvgPeakHeight];
    AvgTimeToPeak = [AvgTimeToPeak, Stimulation(i).AvgTimeToPeak];
    Repetition{i} = Stimulation(i).Repetition;
    PreAvgDFF = [PreAvgDFF, Prestim(i).AvgDFF];
    PreAvgTrendline = [PreAvgTrendline, Prestim(i).AvgTrendline];
    PreAvgPeakHeight = [PreAvgPeakHeight, Prestim(i).AvgPeakHeight];
    PreAvgAUC = [PreAvgAUC, Prestim(i).AvgAUC];
    PreRepetition{i} = Prestim(i).Repetition;
end

% Make a title list to save all data (names cannot occur more than once, so
% add a counter before the names)
title_list_for_save = cell(1, length(Stimulation));
for i = 1:length(Stimulation)
    temp_title = ['(', num2str(i), ') ', Stimulation(i).title];
    temp_title = join(temp_title);
    temp_title = convertStringsToChars(temp_title);
    title_list_for_save{i} = temp_title;
end

my_T = cell2table(stim_response', 'VariableNames', {'Response'});
T1 = array2table(dFF_values');
T1 = rows2vars(T1);
T2 = table(Repetition', AUC_values', str2double(AvgAUC)', trendline_values', str2double(AvgTrendline)', peak_values', str2double(AvgPeakHeight)', time_to_peak', str2double(AvgTimeToPeak)', AP_frequency', imaging_freq', 'VariableNames', {'Repetition', 'AUC', 'Average_AUC_from_repetitions', 'Trendline_Slope', 'Average_Trendline_slope_from_repetitions', 'Peak_Height', 'Average_peak_height_from_repetitions', 'Time_to_Peak', 'Average_time_to_peak_from_repetitions', 'AP_Frequency', 'Imaging_Frequency'}, 'RowNames', title_list_for_save);
T3 = horzcat(T2, my_T, T1);

% Save prestim data in table
for i = 1:length(Prestim)
    base_Average_dFF(i) = Prestim(i).Average;
    base_peak(i) = Prestim(i).Peak;
    base_slope_trendline(i) = Prestim(i).Trendline;
    base_AUC(i) = Prestim(i).AUC_total;
end
T4 = table(PreRepetition', base_Average_dFF', str2double(PreAvgDFF)', base_peak', str2double(PreAvgPeakHeight)', base_slope_trendline', str2double(PreAvgTrendline)', base_AUC', str2double(PreAvgAUC)', 'VariableNames', {'Repetition', 'Average_dF_F', 'Average_dF_F_from_repetitions', 'Peak', 'Average_Peak_from_repetitions', 'Slope_trendline', 'Average_Slope_Trendline_from_repetitions', 'AUC', 'Average_AUC_from_repetitions'}, 'RowNames', title_list_for_save);

save('data.mat')
writetable(T3, 'stim_data.xlsx', 'WriteRowNames', true)
writetable(T4, 'prestim_data.xlsx', 'WriteRowNames', true)
% writetable(T5, 'statistics.xlsx', 'WriteRowNames', true)
writetable(T6, 'Average Response per Subcellular Region.xlsx', 'WriteRowNames', true)

print(figure(4), '-dpdf', 'PDF Ephys + Calcium response.pdf');
print(figure(5), '-dpdf', 'PDF AUC.pdf');
print(figure(6), '-dpdf', 'PDF AP frequency - AUC .pdf');
print(figure(7), '-dpdf', 'PDF Trendlines.pdf');
print(figure(8), '-dpdf', 'PDF AP frequency - trendline.pdf');
print(figure(9), '-dpdf', 'PDF Peak Detection.pdf');
print(figure(10), '-dpdf', 'PDF AP frequency - Peak Height.pdf');
print(figure(11), '-dpdf', 'PDF AP frequency - Peak Location.pdf');
print(figure(12), '-dpdf', 'PDF AP frequency.pdf');
print(figure(13), '-dpdf', 'PDF Repetitions Trendline.pdf');
print(figure(14), '-dpdf', 'PDF Repititions AUC.pdf');
print(figure(15), '-dpdf', 'PDF Repititions Peak Height.pdf');
print(figure(16), '-dpdf', 'PDF Repetitions Peak Location.pdf');

disp('Done!')
close all;

%% Warning if some trials are referred to as 'Other'

if numel(Other_titles(~cellfun('isempty',Other_titles))) > 0
    warning('Some trials have not been divided into the right subcellular region!')
end

clear ();

%% Use the code below to make heatmaps

% % Heatmap z-stack (select lines and press ctr + t to uncomment)
% D = pwd; % directory path
% S = dir(fullfile(D,'*G_drift_cor*.tif*')); % get list of files in directory
% [~,ndx] = natsortfiles({S.name}); % indices of correct order
% S = S(ndx);
% 
% % Pick a trial number of which you want a heatmap-stack
% trial_numb = 1;
% trial = S(trial_numb).name;
% 
% for i = 1:size(imfinfo(trial),1);
%     I = imread(trial,i);
%     clf(figure(1));
%     figure(1),
%     h = pcolor(I);
%     h.FaceColor = 'interp';
%     colormap jet;
%     caxis([250 2500]);
%     arrayfun(@(h) set(h,'EdgeColor','none'), findobj(gcf,'type','surface'));
%     set(gca,'XColor', 'none','YColor','none')
%     c = colorbar
%     c.Ticks = [350, 800];
%     c.TickLabels = {'\bf Low Ca^{2+}', '\bf High Ca^{2+}'};
%     saveas(figure(1),sprintf('heatmap%d.tif',i)); % will create heatmap1, heatmap2,...
% end
% 
% V = dir(fullfile(D,'heatmap*.tif*')); % get list of files in directory
% [~,ndx] = natsortfiles({V.name}); % indices of correct order
% V = V(ndx);
% 
% for i = 1:length(V)
%     I = imread(V(i).name);
%     imwrite(I, 'AxonProx_heatmap.tif', 'WriteMode', 'append');
% end
% 
% delete heatmap*.tif


