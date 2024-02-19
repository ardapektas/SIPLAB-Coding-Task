clearvars
close all

%% Import Data

% Import csv file
data_file= "ptt_dataset\gt01.csv";
opts = detectImportOptions(data_file);
data= readmatrix(data_file,opts);
% Create variables in column header names and assign the data in each column to these variables
for i = 1:length(opts.VariableNames)
    assignin('base', opts.VariableNames{i}, data(:,i));
end

%% Pan-Tompkins algorithm

% Sampling frequency
fs= 2000;

% Bandpass Butterworth filter
[b, a]=butter(1,[5, 15]*2/fs,"bandpass");
filtered_ecg= filtfilt(b,a, chestSternumECG);
filtered_ecg= filtered_ecg/max(abs(filtered_ecg));

% Diffrentiator
differentiated_ecg= filtfilt([1 2 0 -2 -1]*fs/8,1, filtered_ecg);
differentiated_ecg= differentiated_ecg/max(abs(differentiated_ecg));

% Squaring operation
squared_ecg = differentiated_ecg.^2;

% Moving filter integrator
moved_integrated_ecg = movsum(squared_ecg, 15);

% Thresholding
threshold_ecg=0.005*max(moved_integrated_ecg);
[peak_values, peak_indices]= findpeaks(moved_integrated_ecg, 'MinPeakDistance', round(0.2*fs),'MinPeakHeight', threshold_ecg);

% R points
r_points= zeros(size(chestSternumECG));
for i=1:length(peak_indices)
    r_points(peak_indices(i))= abs(chestSternumECG(peak_indices(i)));
end

% Smoothing BCG and PPG
filtered_bcg = medfilt1(biopacBCG,10);
filtered_ppg = medfilt1(biopacPPG,10);


% Finding session intervals
i0 = find(session==0);
i1 = find(session==1);
i2 = find(session==2);
i3 = find(session==3);
i4 = find(session==4);
i5 = find(session==5);
i6 = find(session==6);
i7 = find(session==7);
i8 = find(session==8);
i9 = find(session==9);
i10 = find(session==10);

% Finding R points for each session 
peak_indices0= [];
peak_indices1= [];
peak_indices2= [];
peak_indices3= [];
peak_indices4= [];
peak_indices5= [];
peak_indices6= [];
peak_indices7= [];
peak_indices8= [];
peak_indices9= [];
peak_indices10= [];

for i = 1:length(peak_values)
    if peak_indices(i)<= i0(end) 
        peak_indices0 = [peak_indices0; peak_indices(i)];
    elseif peak_indices(i)<= i1(end) 
        peak_indices1 = [peak_indices1; peak_indices(i)];
    elseif peak_indices(i)<= i2(end) 
        peak_indices2 = [peak_indices2; peak_indices(i)];
    elseif peak_indices(i)<= i3(end) 
        peak_indices3 = [peak_indices3; peak_indices(i)];
    elseif peak_indices(i)<= i4(end) 
        peak_indices4 = [peak_indices4; peak_indices(i)];
    elseif peak_indices(i)<= i5(end) 
        peak_indices5 = [peak_indices5; peak_indices(i)];
    elseif peak_indices(i)<= i6(end) 
        peak_indices6 = [peak_indices6; peak_indices(i)];
    elseif peak_indices(i)<= i7(end) 
        peak_indices7 = [peak_indices7; peak_indices(i)];
    elseif peak_indices(i)<= i8(end) 
        peak_indices8 = [peak_indices8; peak_indices(i)];
    elseif peak_indices(i)<= i9(end) 
        peak_indices9 = [peak_indices9; peak_indices(i)];
    elseif peak_indices(i)<= i10(end) 
        peak_indices10 = [peak_indices10; peak_indices(i)];
    end
end

% Create correlation vecors
correlation_diastole = zeros(11,1);
correlation_systole  = zeros(11,1);

% Set window size
window_size= 200;

% Use 2-FIR filter to predict values instead of NaN in MSE sense. Use
% sigmoid function to stabilize
nan_diastolic_indices= find(isnan(finapresDiastolic)==1);
non_nan_diastolic_indices= find(isnan(finapresDiastolic)==0);
rd= autocorr(finapresDiastolic(non_nan_diastolic_indices));
opt_d=inv([2*rd(1)-2*rd(2), 2*rd(2)-rd(3)-rd(1);2*rd(2)-rd(3)-rd(1),2*rd(1)-2*rd(2)])*[rd(2)-rd(3),rd(3)-rd(4)]';
opt_d=opt_d/norm(opt_d);
for i=1:length(nan_diastolic_indices)
    finapresDiastolic(nan_diastolic_indices(i))= 1/(1+exp(-1*([ finapresDiastolic(nan_diastolic_indices(i)-1), finapresDiastolic(nan_diastolic_indices(i)-2)] *opt_d)));
end

nan_systolic_indices= find(isnan(finapresSystolic)==1);
non_nan_systolic_indices= find(isnan(finapresSystolic)==0);
rs= autocorr(finapresSystolic(non_nan_systolic_indices));
opt_s=inv([2*rs(1)-2*rs(2), 2*rs(2)-rs(3)-rs(1);2*rs(2)-rs(3)-rs(1),2*rs(1)-2*rs(2)])*[rs(2)-rs(3),rs(3)-rs(4)]';
opt_s=opt_s/norm(opt_s);
for i=1:length(nan_systolic_indices)
    finapresSystolic(nan_systolic_indices(i))= 1/(1+exp(-1*([ finapresSystolic(nan_systolic_indices(i)-1), finapresSystolic(nan_systolic_indices(i)-2)] *opt_s)));;
end

if(isempty(peak_indices0)==0)
%% Session 0

    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices0));
    diastole_peaks = zeros(window_size,length(peak_indices0));
    systole_peaks = zeros(window_size,length(peak_indices0));
    bcg_peaks = zeros(window_size,length(peak_indices0));
    ppg_peaks =zeros(window_size,length(peak_indices0));
    ppg_peak_indices = zeros(size(peak_indices0));
    bcg_peak_indices = zeros(size(peak_indices0));
    
    for i = 1:length(peak_indices0)
        start_point = max(1, peak_indices0(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices0(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices0), 1);
    
    for i = 1:length(peak_indices0)
        peak_index = peak_indices0(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(1) = corr(ptt_values, diastole_peaks);
    correlation_systole(1) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 0: ', num2str(correlation_diastole(1))]);
    disp(['Correlation coefficient between PTT and Systole for Session 0: ', num2str(correlation_systole(1))]);
end
if(isempty(peak_indices1)==0)
    %% Session 1
    
    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices1));
    diastole_peaks = zeros(window_size,length(peak_indices1));
    systole_peaks = zeros(window_size,length(peak_indices1));
    bcg_peaks = zeros(window_size,length(peak_indices1));
    ppg_peaks =zeros(window_size,length(peak_indices1));
    ppg_peak_indices = zeros(size(peak_indices1));
    bcg_peak_indices = zeros(size(peak_indices1));
    
    for i = 1:length(peak_indices1)
        start_point = max(1, peak_indices1(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices1(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices1), 1);
    
    for i = 1:length(peak_indices1)
        peak_index = peak_indices1(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(2) = corr(ptt_values, diastole_peaks);
    correlation_systole(2) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 1: ', num2str(correlation_diastole(2))]);
    disp(['Correlation coefficient between PTT and Systole for Session 1: ', num2str(correlation_systole(2))]);
end
if(isempty(peak_indices2)==0)
    %% Session 2
    
    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices2));
    diastole_peaks = zeros(window_size,length(peak_indices2));
    systole_peaks = zeros(window_size,length(peak_indices2));
    bcg_peaks = zeros(window_size,length(peak_indices2));
    ppg_peaks =zeros(window_size,length(peak_indices2));
    ppg_peak_indices = zeros(size(peak_indices2));
    bcg_peak_indices = zeros(size(peak_indices2));
    
    for i = 1:length(peak_indices2)
        start_point = max(1, peak_indices2(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices2(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices2), 1);
    
    for i = 1:length(peak_indices2)
        peak_index = peak_indices2(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(3) = corr(ptt_values, diastole_peaks);
    correlation_systole(3) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 2: ', num2str(correlation_diastole(3))]);
    disp(['Correlation coefficient between PTT and Systole for Session 2: ', num2str(correlation_systole(3))]);
end
if(isempty(peak_indices3)==0)
    %% Session 3
    
    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices3));
    diastole_peaks = zeros(window_size,length(peak_indices3));
    systole_peaks = zeros(window_size,length(peak_indices3));
    bcg_peaks = zeros(window_size,length(peak_indices3));
    ppg_peaks =zeros(window_size,length(peak_indices3));
    ppg_peak_indices = zeros(size(peak_indices3));
    bcg_peak_indices = zeros(size(peak_indices3));
    
    for i = 1:length(peak_indices3)
        start_point = max(1, peak_indices3(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices3(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices3), 1);
    
    for i = 1:length(peak_indices3)
        peak_index = peak_indices3(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(4) = corr(ptt_values, diastole_peaks);
    correlation_systole(4) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 3: ', num2str(correlation_diastole(4))]);
    disp(['Correlation coefficient between PTT and Systole for Session 3: ', num2str(correlation_systole(4))]);
end
if(isempty(peak_indices4)==0)
    %% Session 4
    
    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices4));
    diastole_peaks = zeros(window_size,length(peak_indices4));
    systole_peaks = zeros(window_size,length(peak_indices4));
    bcg_peaks = zeros(window_size,length(peak_indices4));
    ppg_peaks =zeros(window_size,length(peak_indices4));
    ppg_peak_indices = zeros(size(peak_indices4));
    bcg_peak_indices = zeros(size(peak_indices4));
    
    for i = 1:length(peak_indices4)
        start_point = max(1, peak_indices4(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices4(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices4), 1);
    
    for i = 1:length(peak_indices4)
        peak_index = peak_indices4(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(5) = corr(ptt_values, diastole_peaks);
    correlation_systole(5) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 4: ', num2str(correlation_diastole(5))]);
    disp(['Correlation coefficient between PTT and Systole for Session 4: ', num2str(correlation_systole(5))]);
end
if(isempty(peak_indices5)==0)
    %% Session 5
    
    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices5));
    diastole_peaks = zeros(window_size,length(peak_indices5));
    systole_peaks = zeros(window_size,length(peak_indices5));
    bcg_peaks = zeros(window_size,length(peak_indices5));
    ppg_peaks =zeros(window_size,length(peak_indices5));
    ppg_peak_indices = zeros(size(peak_indices5));
    bcg_peak_indices = zeros(size(peak_indices5));
    
    for i = 1:length(peak_indices5)
        start_point = max(1, peak_indices5(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices5(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices5), 1);
    
    for i = 1:length(peak_indices5)
        peak_index = peak_indices5(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(6) = corr(ptt_values, diastole_peaks);
    correlation_systole(6) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 5: ', num2str(correlation_diastole(6))]);
    disp(['Correlation coefficient between PTT and Systole for Session 5: ', num2str(correlation_systole(6))]);
end
if(isempty(peak_indices6)==0)
    %% Session 6
    
    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices6));
    diastole_peaks = zeros(window_size,length(peak_indices6));
    systole_peaks = zeros(window_size,length(peak_indices6));
    bcg_peaks = zeros(window_size,length(peak_indices6));
    ppg_peaks =zeros(window_size,length(peak_indices6));
    ppg_peak_indices = zeros(size(peak_indices6));
    bcg_peak_indices = zeros(size(peak_indices6));
    
    for i = 1:length(peak_indices6)
        start_point = max(1, peak_indices6(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices6(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices6), 1);
    
    for i = 1:length(peak_indices6)
        peak_index = peak_indices6(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(7) = corr(ptt_values, diastole_peaks);
    correlation_systole(7) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 6: ', num2str(correlation_diastole(7))]);
    disp(['Correlation coefficient between PTT and Systole for Session 6: ', num2str(correlation_systole(7))]);
end
if(isempty(peak_indices7)==0)
    %% Session 7
    
    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices7));
    diastole_peaks = zeros(window_size,length(peak_indices7));
    systole_peaks = zeros(window_size,length(peak_indices7));
    bcg_peaks = zeros(window_size,length(peak_indices7));
    ppg_peaks =zeros(window_size,length(peak_indices7));
    ppg_peak_indices = zeros(size(peak_indices7));
    bcg_peak_indices = zeros(size(peak_indices7));
    
    for i = 1:length(peak_indices7)
        start_point = max(1, peak_indices7(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices7(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices7), 1);
    
    for i = 1:length(peak_indices7)
        peak_index = peak_indices7(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(8) = corr(ptt_values, diastole_peaks);
    correlation_systole(8) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 7: ', num2str(correlation_diastole(8))]);
    disp(['Correlation coefficient between PTT and Systole for Session 7: ', num2str(correlation_systole(8))]);
end
if(isempty(peak_indices8)==0)
    %% Session 8
    
    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices8));
    diastole_peaks = zeros(window_size,length(peak_indices8));
    systole_peaks = zeros(window_size,length(peak_indices8));
    bcg_peaks = zeros(window_size,length(peak_indices8));
    ppg_peaks =zeros(window_size,length(peak_indices8));
    ppg_peak_indices = zeros(size(peak_indices8));
    bcg_peak_indices = zeros(size(peak_indices8));
    
    for i = 1:length(peak_indices8)
        start_point = max(1, peak_indices8(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices8(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices8), 1);
    
    for i = 1:length(peak_indices8)
        peak_index = peak_indices8(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(9) = corr(ptt_values, diastole_peaks);
    correlation_systole(9) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 8: ', num2str(correlation_diastole(9))]);
    disp(['Correlation coefficient between PTT and Systole for Session 8: ', num2str(correlation_systole(9))]);
end
if(isempty(peak_indices9)==0)
    %% Session 9
    
    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices9));
    diastole_peaks = zeros(window_size,length(peak_indices9));
    systole_peaks = zeros(window_size,length(peak_indices9));
    bcg_peaks = zeros(window_size,length(peak_indices9));
    ppg_peaks =zeros(window_size,length(peak_indices9));
    ppg_peak_indices = zeros(size(peak_indices9));
    bcg_peak_indices = zeros(size(peak_indices9));
    
    for i = 1:length(peak_indices9)
        start_point = max(1, peak_indices9(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices9(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices9), 1);
    
    for i = 1:length(peak_indices9)
        peak_index = peak_indices9(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(10) = corr(ptt_values, diastole_peaks);
    correlation_systole(10) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 9: ', num2str(correlation_diastole(10))]);
    disp(['Correlation coefficient between PTT and Systole for Session 9: ', num2str(correlation_systole(10))]);
end
if(isempty(peak_indices10)==0)
    %% Session 10
    
    % Extract beats around each R point within selected window size and find peak indices
    bp_peaks =  zeros(window_size,length(peak_indices10));
    diastole_peaks = zeros(window_size,length(peak_indices10));
    systole_peaks = zeros(window_size,length(peak_indices10));
    bcg_peaks = zeros(window_size,length(peak_indices10));
    ppg_peaks =zeros(window_size,length(peak_indices10));
    ppg_peak_indices = zeros(size(peak_indices10));
    bcg_peak_indices = zeros(size(peak_indices10));
    
    for i = 1:length(peak_indices10)
        start_point = max(1, peak_indices10(i) - window_size/2);
        end_point = min(length(finapresBP),peak_indices10(i) + window_size/2 - 1);
        bp_peaks(1:end_point-start_point+1,i) = finapresBP(start_point:end_point);
        diastole_peaks(1:end_point-start_point+1,i) = finapresDiastolic(start_point:end_point);
        systole_peaks(1:end_point-start_point+1,i) =  finapresSystolic(start_point:end_point);
        bcg_peaks(1:end_point-start_point+1,i) =filtered_bcg(start_point:end_point);
        [~,I]= max(filtered_bcg(start_point:end_point));
        bcg_peak_indices(i)= start_point+I;
        ppg_peaks(1:end_point-start_point+1,i) = filtered_ppg(start_point:end_point);
        [~,I]= max(filtered_ppg(start_point:end_point));
        ppg_peak_indices(i)= start_point+I;
    end
    
    bp_peaks = max(bp_peaks)';
    diastole_peaks = max(diastole_peaks)';
    systole_peaks = max(systole_peaks)';
    bcg_peaks = max(bcg_peaks)';
    ppg_peaks = max(ppg_peaks)';
    
    % Find PTT values
    ptt_values = zeros(length(peak_indices10), 1);
    
    for i = 1:length(peak_indices10)
        peak_index = peak_indices10(i);
        [~, ppg_peak_index] = min(abs(ppg_peak_indices - peak_index));
        [~, bcg_peak_index] = min(abs(bcg_peak_indices - peak_index));
        ppg_peak_time = ppg_peak_index / fs;
        bcg_peak_time = bcg_peak_index / fs;
        r_peak_time = peak_index / fs;
        ptt_values(i) = (abs(ppg_peak_time - r_peak_time)+abs(bcg_peak_time - r_peak_time))/2;
    end
    
    % Calculate correlation coefficients
    correlation_diastole(11) = corr(ptt_values, diastole_peaks);
    correlation_systole(11) = corr(ptt_values, systole_peaks);
    
    disp(['Correlation coefficient between PTT and Diastole for Session 10: ', num2str(correlation_diastole(11))]);
    disp(['Correlation coefficient between PTT and Systole for Session 10: ', num2str(correlation_systole(11))]);
end

%% R Points

% Plot ECG and R points
figure
plot(chestSternumECG)
hold on
plot(r_points,'ro')
title('ECG and R points')

%% Overall Correlation Coefficients

figure
stem([0:10],correlation_diastole)
hold on
stem([0:10],correlation_systole)
legend
legend('Diastole','Systole')
xlabel('Session')
ylabel('Correlation Coefficient')
title('Correlation Coefficients for each Session')