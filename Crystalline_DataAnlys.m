%------------------------------------------------------------------------------------------------
% Application: Crystalline Image and Data Analysis Framework
% Subroutine: Crystalline_DataAnlys
% 
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2020b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/FrederikDoerr


function[MDB_Data,R_DataAnlys] = Crystalline_DataAnlys(Opt,Opt_DS,MDB_File,MDB_Img)

fprintf('%s - %s (%s) INITIATED (Elapsed time: %.0f sec)\n',Opt.ProjectShorthand,Opt_DS.ExpShorthand,mfilename(),toc(Opt.tic))

%% Graphs - File Processing
fig = figure('units','pixel','position',[100 100 Opt_DS.Graph.width Opt_DS.Graph.height]);
ax = axes('Parent',fig,...
    'Position',[Opt_DS.Graph.ax_1,Opt_DS.Graph.ax_2,Opt_DS.Graph.ax_1_width,Opt_DS.Graph.ax_2_height]);
hold(ax,'on');
box(ax,'on');
ax.XAxis.Exponent = 0;
xtickformat('%.0f')

yyaxis left
ax.YColor = 'r';
plot(MDB_File.file_time_list,MDB_File.Temperature,'r-','DisplayName','Temperature')

xlabel('Time [s]')
ylabel('Temperature')
ax = gca;
ax.XLim = [min(MDB_File.file_time_list) max(MDB_File.file_time_list)];

yyaxis right

plot(MDB_File.file_time_list,MDB_File.Transmissivity,'g-','DisplayName','Transmissivity')
ax.YColor = 'g'; 
ylabel('Transmissivity')

print(fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_File_Temp_Trans',Opt_DS.ExpShorthand)),'-djpeg',Opt_DS.Graph.print_rS)


%% Find heating/cooling sequences from temperature profile
time = MDB_File.file_time_list(1:end);
Temperature = MDB_File.Temperature(1:end);
Temp_diff = movsum(gradient(Temperature,1),Opt_DS.DataPrc.k_movsum_Temp_diff);
% compensate for changing sampling rate
Temp_diff = Temp_diff./gradient(time);
Temp_diff_sm = sgolayfilt(Temp_diff, 1, Opt_DS.DataPrc.av_StepSize);

% Temp_diff for plotting [degC/min]
Temp_diff_sm_p = Temp_diff_sm./Opt_DS.DataPrc.k_movsum_Temp_diff*60;

fig = figure('units','pixel','position',[100 100 Opt_DS.Graph.width Opt_DS.Graph.height]);
ax = axes('Parent',fig,...
    'Position',[Opt_DS.Graph.ax_1,Opt_DS.Graph.ax_2,Opt_DS.Graph.ax_1_width,Opt_DS.Graph.ax_2_height]);
hold(ax,'on');
box(ax,'on');
ax.XAxis.Exponent = 0;
xtickformat('%.0f')
ax.YColor = Opt_DS.Graph.Color_Seq_grey{1};
hold on

yyaxis left
plot(time,Temperature,'LineStyle','-','Color',Opt_DS.Graph.Color_Seq_grey{9},'LineWidth',1.5,'DisplayName','Temperature')
yyaxis right
plot(time,Temp_diff_sm_p,'LineStyle','-','Color',Opt_DS.Graph.Color_Seq_grey{5},'LineWidth',1.5,'DisplayName','Heating Rate')

ax = gca;
ax.XLim = [min(time) max(time)];
yyaxis left
ax.YColor = 'k';
yyaxis right
ax.YColor = 'k';

% Constant period 
idx_Temp_const = find(abs(Temp_diff_sm) < 1e-4*Opt_DS.DataPrc.k_movsum_Temp_diff);

time_const = time(idx_Temp_const);
Temp_const = Temperature(idx_Temp_const);
Temp_const_diff = gradient(Temp_const,1);

Temp_const_diff_TF = isoutlier(Temp_const_diff,'movmean',length(Temp_const_diff));
idx_TF = find(Temp_const_diff_TF);
% Filter through temperatures that meet Temp_const_thres
idx_TF_DEL = [];
for i_TF = 1:length(idx_TF)
    if abs(Temp_const_diff(idx_TF(i_TF))) < Opt_DS.DataPrc.Temp_const_thres
         idx_TF_DEL = [idx_TF_DEL,i_TF]; %#ok<AGROW>
    end
end    
idx_TF(idx_TF_DEL) = [];
idx_TF = [1;idx_TF;length(Temp_const_diff_TF)];

% linear fit Temperature and error estimation
LinReg_Offset = 0.1;
p_list = nan(length(idx_TF)-1,2);
T_mean_list = nan(length(idx_TF)-1,1);
T_std_list = nan(length(idx_TF)-1,1);
idx_TF_list = nan(length(idx_TF)-1,2);
for i_TF = 1:length(idx_TF)-1
    idx1 = idx_TF(i_TF);
    idx2 = idx_TF(i_TF+1);
    d_idx = idx2-idx1;
    idx1 = floor(idx1 +d_idx*LinReg_Offset);
    idx2 = ceil(idx2 -d_idx*LinReg_Offset);

    if idx2-idx1 > 1
        [p_list(i_TF,:),~,~] = polyfit(time_const(idx1:idx2),Temp_const(idx1:idx2),1);
        T_mean_list(i_TF) = mean(Temp_const(idx1:idx2));
        T_std_list(i_TF) = std(Temp_const(idx1:idx2));
        idx_TF_list(i_TF,:) = [idx_Temp_const(idx1),idx_Temp_const(idx2)];

        yyaxis left
        plot(time(idx_TF_list(i_TF,1):idx_TF_list(i_TF,2)), ...
            Temperature(idx_TF_list(i_TF,1):idx_TF_list(i_TF,2)),'r--','DisplayName','Temp Constant')
    end
end

% Reove nan
idx_row_SEL = find(~isnan(T_mean_list));
idx_TF_list = idx_TF_list(idx_row_SEL,:);
T_mean_list = T_mean_list(idx_row_SEL,:);
T_std_list = T_std_list(idx_row_SEL,:);


% Find cooling/heating
Temp_Det_cooling = nan(length(T_mean_list)-1,1);
Temp_Det_heating = nan(length(T_mean_list)-1,1);
idx_Det_TempRamp_init = nan(length(T_mean_list)-1,1);
idx_Det_TempRamp_end = nan(length(T_mean_list)-1,1);
for i_TF = 1:length(T_mean_list)-1
    idx2 = idx_TF_list(i_TF,2);
    if i_TF < length(T_mean_list)
        idx_nT = idx_TF_list(i_TF+1,1);
    else
        idx_nT = length(time);
    end

    idx = idx2;
    while idx < idx_nT
        Temperature(idx);
        Res = abs(mean(Temperature(idx-5:idx+5)) - T_mean_list(i_TF));
        if (Res > 5*T_std_list(i_TF) || Res > Opt_DS.DataPrc.Prc_Tthres) && Temperature(idx) > T_mean_list(i_TF)+4
            Temp_Det_cooling(i_TF) = false;
            Temp_Det_heating(i_TF) = true;

            idx_Det_TempRamp_init(i_TF) = idx;
            while idx < idx_nT && median(Temperature(idx-2:idx+2)) < T_mean_list(i_TF+1)-0.2
                idx = idx +1;
            end
            idx_Det_TempRamp_end(i_TF) = idx;
            
            % Refine transition points
            numIdx_Temp = idx_Det_TempRamp_end(i_TF) - idx_Det_TempRamp_init(i_TF);
            numIdx_Temp_offset = round(numIdx_Temp*0.25);
            xi_time = time(idx_Det_TempRamp_init(i_TF)+numIdx_Temp_offset:idx_Det_TempRamp_end(i_TF)-numIdx_Temp_offset) - time(idx_Det_TempRamp_init(i_TF)+numIdx_Temp_offset);
            yi_Temperature = Temperature(idx_Det_TempRamp_init(i_TF)+numIdx_Temp_offset:idx_Det_TempRamp_end(i_TF)-numIdx_Temp_offset);
            [p_TempRamp,~] = polyfit(xi_time,yi_Temperature,1);
            
            t_TempRamp_init_Corr =  (T_mean_list(i_TF) - p_TempRamp(2))/p_TempRamp(1);
            t_TempRamp_end_Corr = (T_mean_list(i_TF+1) - p_TempRamp(2))/p_TempRamp(1);
            time_Det_TempRamp_init_fit = time(idx_Det_TempRamp_init(i_TF)+numIdx_Temp_offset) + t_TempRamp_init_Corr;
            time_Det_TempRamp_init_end = time(idx_Det_TempRamp_init(i_TF)+numIdx_Temp_offset) + t_TempRamp_end_Corr;
            idx_Det_TempRamp_init(i_TF) = find(abs(time-time_Det_TempRamp_init_fit) == min(abs(time-time_Det_TempRamp_init_fit)),1,'first');
            idx_Det_TempRamp_end(i_TF) = find(abs(time-time_Det_TempRamp_init_end) == min(abs(time-time_Det_TempRamp_init_end)),1,'first');
            
            
            line([time(idx_Det_TempRamp_init(i_TF)),time(idx_Det_TempRamp_init(i_TF))],[min(Temperature),max(Temperature)],'Color','red')
            line([time(idx_Det_TempRamp_end(i_TF)),time(idx_Det_TempRamp_end(i_TF))],[min(Temperature),max(Temperature)],'Color','red')
            time_Det_grad_init = time(idx_Det_TempRamp_init(i_TF));
            time_Det_grad_end = time(idx_Det_TempRamp_end(i_TF)); 
            Heat_ramp = (Temperature(idx_Det_TempRamp_end(i_TF)) - Temperature(idx_Det_TempRamp_init(i_TF)))/(time_Det_grad_end-time_Det_grad_init);
            text(time(idx)+time(idx)*0.03,max(Temperature)+(max(Temperature)-min(Temperature))*-0.1, ...
                sprintf('Heating Seq:\nT = %.2f > %.2f\nHeatRate %.4f K/s\nHeatRate %.2f K/min',Temperature(idx_Det_TempRamp_init(i_TF)),Temperature(idx_Det_TempRamp_end(i_TF)), ...
                Heat_ramp,Heat_ramp*60),'Color','r','FontSize',4)
            break
        elseif (Res > 5*T_std_list(i_TF) || Res > Opt_DS.DataPrc.Prc_Tthres) && Temperature(idx) < T_mean_list(i_TF)-4
            Temp_Det_cooling(i_TF) = true;
            Temp_Det_heating(i_TF) = false;

            idx_Det_TempRamp_init(i_TF) = idx;
            while idx < idx_nT && median(Temperature(idx-2:idx+2)) > T_mean_list(i_TF+1)+0.2
                idx = idx +1;
            end
            idx_Det_TempRamp_end(i_TF) = idx;

            % Refine transition points
            numIdx_Temp = idx_Det_TempRamp_end(i_TF) - idx_Det_TempRamp_init(i_TF);
            numIdx_Temp_offset = round(numIdx_Temp*0.35);
            xi_time = time(idx_Det_TempRamp_init(i_TF)+numIdx_Temp_offset:idx_Det_TempRamp_end(i_TF)-numIdx_Temp_offset) - time(idx_Det_TempRamp_init(i_TF)+numIdx_Temp_offset);
            yi_Temperature = Temperature(idx_Det_TempRamp_init(i_TF)+numIdx_Temp_offset:idx_Det_TempRamp_end(i_TF)-numIdx_Temp_offset);
            [p_TempRamp,~] = polyfit(xi_time,yi_Temperature,1);
            
            t_TempRamp_init_Corr =  (T_mean_list(i_TF) - p_TempRamp(2))/p_TempRamp(1);
            t_TempRamp_end_Corr = (T_mean_list(i_TF+1) - p_TempRamp(2))/p_TempRamp(1);
            time_Det_TempRamp_init_fit = time(idx_Det_TempRamp_init(i_TF)+numIdx_Temp_offset) + t_TempRamp_init_Corr;
            time_Det_TempRamp_init_end = time(idx_Det_TempRamp_init(i_TF)+numIdx_Temp_offset) + t_TempRamp_end_Corr;
            idx_Det_TempRamp_init(i_TF) = find(abs(time-time_Det_TempRamp_init_fit) == min(abs(time-time_Det_TempRamp_init_fit)),1,'first');
            idx_Det_TempRamp_end(i_TF) = find(abs(time-time_Det_TempRamp_init_end) == min(abs(time-time_Det_TempRamp_init_end)),1,'first');
            
            line([time(idx_Det_TempRamp_init(i_TF)),time(idx_Det_TempRamp_init(i_TF))],[min(Temperature),max(Temperature)],'Color','blue')
            line([time(idx_Det_TempRamp_end(i_TF)),time(idx_Det_TempRamp_end(i_TF))],[min(Temperature),max(Temperature)],'Color','blue')
            time_Det_grad_init = time(idx_Det_TempRamp_init(i_TF));
            time_Det_grad_end = time(idx_Det_TempRamp_end(i_TF));
            Heat_ramp = (Temperature(idx_Det_TempRamp_end(i_TF)) - Temperature(idx_Det_TempRamp_init(i_TF)))/(time_Det_grad_end-time_Det_grad_init);

            text(time(idx)+time(idx)*0.03,min(Temperature)+(max(Temperature)-min(Temperature))*+0.1, ...
                sprintf('Cooling Seq:\nT = %.2f > %.2f\nHeatRate %.4f K/s\nHeatRate %.2f K/min',Temperature(idx_Det_TempRamp_init(i_TF)),Temperature(idx_Det_TempRamp_end(i_TF)), ...
                Heat_ramp,Heat_ramp*60),'Color','b','FontSize',4)
            break
        end
        idx = idx +1;
    end 
end


xlabel('Time [s]')
yyaxis left
ylabel('Temperature [degC]')
yyaxis right
ylabel('Heating Rate [degC/min]')

print(fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_ImgAnlys_Heat_Cool_Det',Opt_DS.ExpShorthand)),'-djpeg',Opt_DS.Graph.print_rS)

%% Image tags: 0 = constant temperature, 1 = heating, -1 = cooling
num_TempStep = sum(~isnan(Temp_Det_cooling));

Temp_mode_Img = zeros(length(MDB_Img.Img_PrcIdx),1);
Temp_mode_File = zeros(length(time),1);
for k =1:num_TempStep
    
    if Temp_Det_heating(k)
        Temp_mode_File(idx_Det_TempRamp_init(k):idx_Det_TempRamp_end(k)) = 1;
    else
        Temp_mode_File(idx_Det_TempRamp_init(k):idx_Det_TempRamp_end(k)) = -1;
    end
    
    idx_Img_grad_t1 = find(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_init(k)))) == min(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_init(k))))),1,'last');
    idx_Img_grad_t2 = find(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_end(k)))) == min(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_end(k))))),1,'first');

    if Temp_Det_heating(k)
        Temp_mode_Img(idx_Img_grad_t1:idx_Img_grad_t2) = 1;
    else
        Temp_mode_Img(idx_Img_grad_t1:idx_Img_grad_t2) = -1;
    end
    
end

%% Detect in-sample PLS calibration for crsytal suspension density estimation

% Conditions: 
%   1) series of heating steps (no cooling), min 3 consecutive  heating stages
%   2) first heating after exp start

PLS_logic = zeros(length(Temp_Det_heating),1);
T_mean_grad = gradient(T_mean_list);
iter = 1;
while iter < length(Temp_Det_heating) && ...
    T_mean_grad(iter) > 0
    
    PLS_logic(iter) = true;
    
    iter = iter+1;
end

% Assuming at least PLS_logic_nStages_min consequtive heating stages to idicate PLSR Data
PLS_logic_Img = zeros(length(MDB_Img.Img_PrcIdx),1);


if sum(PLS_logic) > Opt_DS.DataPrc.PLS_logic_nStages_min
    % Isolate PLS dataset
    i_c = find(PLS_logic,1,'last')+1;
    idx_Img_grad_t1 = find(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_init(1)))) == min(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_init(1))))),1,'last');
    idx_Img_grad_t2 = find(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_init(i_c)))) == min(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_init(i_c))))),1,'first')-1;
    idx_File_grad_t1 = find(abs(time-(time(idx_Det_TempRamp_init(1)))) == min(abs(time-(time(idx_Det_TempRamp_init(1))))),1,'last');
    idx_File_grad_t2 = find(abs(time-(time(idx_Det_TempRamp_init(i_c)))) == min(abs(time-(time(idx_Det_TempRamp_init(i_c))))),1,'first')-1;
    
    % Clip data
    PLS_logic_Img(idx_Img_grad_t1:idx_Img_grad_t2) = 1;
    MDB_Img_PLSR = MDB_Img(idx_Img_grad_t1:idx_Img_grad_t2,:);
    File_PLSR = MDB_File(idx_File_grad_t1:idx_File_grad_t2,:);
    
    save(fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_PLSR_Data_MDB_Img.mat',Opt_DS.ExpShorthand)), 'MDB_Img_PLSR');
    save(fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_PLSR_Data_File.mat',Opt_DS.ExpShorthand)), 'File_PLSR');
end


%% Analysis Heating / Cooling Steps
num_TempStep = sum(~isnan(Temp_Det_cooling));
Opt.num_TempStep = sum(~isnan(Temp_Det_cooling));

% Preallocation Data lists
HeatDet_list = nan(num_TempStep,1);
CoolDet_list = nan(num_TempStep,1);
HeatRate_list = nan(num_TempStep,1);
T_Const_1_list = nan(num_TempStep,1);
T_Const_2_list = nan(num_TempStep,1);
time_File_grad_init = nan(num_TempStep,1);
time_File_grad_end = nan(num_TempStep,1);

numImg_list = nan(num_TempStep,1);

Nuc_Time_abs_File_list = nan(num_TempStep,1);
Nuc_Time_Ind_File_list = nan(num_TempStep,1);
Nuc_Time_Ind_File_Err_list = nan(num_TempStep,1);
Nuc_Img_name_File_list = cell(num_TempStep,1);
Nuc_Time_abs_Img_list = nan(num_TempStep,1);
Nuc_Time_Ind_Img_list = nan(num_TempStep,1);
Nuc_Time_Ind_Img_Err_list = nan(num_TempStep,1);
Nuc_Img_name_Img_list = cell(num_TempStep,1);
Nuc_Temp_Img_list = nan(num_TempStep,1);
Nuc_Temp_File_list = nan(num_TempStep,1);
Nuc_Img_IsoThermal_list = nan(num_TempStep,1);
Nuc_File_IsoThermal_list = nan(num_TempStep,1);

Growth_PSD_3_25_list = nan(num_TempStep,1);
Growth_PSD_3_25_Err_list = nan(num_TempStep,1);
Growth_PSD_3_50_list = nan(num_TempStep,1);
Growth_PSD_3_50_Err_list = nan(num_TempStep,1);
Growth_PSD_3_75_list = nan(num_TempStep,1);
Growth_PSD_3_75_Err_list = nan(num_TempStep,1);
Growth_PSD_3_90_list = nan(num_TempStep,1);
Growth_PSD_3_90_Err_list = nan(num_TempStep,1);

for i_c = 1:num_TempStep
    
    % Clear and re-create result structure
    clearvars R_ic
    R_ic = struct();
    
    R_ic.num_idx_ext = 0;
    % Extend sequence 1) to last image or 2) 30 seconds before following
    if i_c == num_TempStep
        R_ic.Ramp_t_ext = MDB_Img.Img_time_s(end)-time(idx_Det_TempRamp_end(i_c));
    else
        R_ic.Ramp_t_ext = time(idx_Det_TempRamp_init(i_c+1))-time(idx_Det_TempRamp_end(i_c))-30;

    end
    
    % Start sequence 1) at first image, 2) end point of previous  TempRamp
    % or 3) 5 mins before start TempRamp
    ro = -5*60;
    if i_c == 1
        if idx_Det_TempRamp_init(i_c)+ro < 0
            R_ic.Ramp_t_offset = MDB_Img.Img_time_s(1);
        else
            R_ic.Ramp_t_offset = ro;
        end
    elseif idx_Det_TempRamp_init(i_c)+ro < idx_Det_TempRamp_end(i_c-1)
        R_ic.Ramp_t_offset = time(idx_Det_TempRamp_end(i_c))+1;
    else
        R_ic.Ramp_t_offset = ro;
    end

    
    
    HeatDet_list(i_c) = Temp_Det_heating(i_c);
    CoolDet_list(i_c) = Temp_Det_cooling(i_c);

    T_Const_1_list(i_c) = T_mean_list(i_c);
    T_Const_2_list(i_c) = T_mean_list(i_c+1);

    R_ic.idx_Img_grad_t1 = find(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_init(i_c)))) == min(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_init(i_c))))),1,'last');
    R_ic.idx_Img_grad_t2 = find(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_end(i_c)))) == min(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_end(i_c))))),1,'first');
    R_ic.idx_File_grad_t1 = find(abs(time-(time(idx_Det_TempRamp_init(i_c)))) == min(abs(time-(time(idx_Det_TempRamp_init(i_c))))),1,'last');
    R_ic.idx_File_grad_t2 = find(abs(time-(time(idx_Det_TempRamp_end(i_c)))) == min(abs(time-(time(idx_Det_TempRamp_end(i_c))))),1,'first');

    time_File_grad_init(i_c) = time(R_ic.idx_File_grad_t1);
    time_File_grad_end(i_c) = time(R_ic.idx_File_grad_t2);

    [R_ic.p_HeatRate,R_ic.S_HeatRate] = polyfit(time(R_ic.idx_File_grad_t1:R_ic.idx_File_grad_t2),Temperature(R_ic.idx_File_grad_t1:R_ic.idx_File_grad_t2),1);
    HeatRate_list(i_c) = R_ic.p_HeatRate(1); % K/s
    

    %% Identify cooling phase
    R_ic.idx_Img_t1 = find(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_init(i_c))+R_ic.Ramp_t_offset)) == min(abs(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_init(i_c))+R_ic.Ramp_t_offset))),1,'first');
    R_ic.idx_Img_t2_max = find(MDB_Img.Img_time_s-(time(idx_Det_TempRamp_end(i_c))+R_ic.Ramp_t_ext)<0,1,'last');
    R_ic.idx_Img_t2 = find(abs(MDB_Img.Img_time_s(1:R_ic.idx_Img_t2_max)-(time(idx_Det_TempRamp_end(i_c))+R_ic.Ramp_t_ext)) == min(abs(MDB_Img.Img_time_s(1:R_ic.idx_Img_t2_max)-(time(idx_Det_TempRamp_end(i_c))+R_ic.Ramp_t_ext))),1,'first');
    if R_ic.idx_Img_t2+R_ic.num_idx_ext > size(MDB_Img.Img_time_s,1)
        R_ic.idx_Img_t2 = size(MDB_Img.Img_time_s,1)-Opt_DS.DataPrc.l_movAvg-1;
    else
        R_ic.idx_Img_t2 = R_ic.idx_Img_t2+R_ic.num_idx_ext-Opt_DS.DataPrc.l_movAvg-1;
    end
    R_ic.idx_Img = R_ic.idx_Img_t1:1:R_ic.idx_Img_t2;

    % Match with file
    R_ic.idx_File_PSD_t1 = find(abs(MDB_File.file_PSD_time_list-(time(idx_Det_TempRamp_init(i_c))+R_ic.Ramp_t_offset)) == min(abs(MDB_File.file_PSD_time_list-(time(idx_Det_TempRamp_init(i_c))+R_ic.Ramp_t_offset))),1,'last');
    R_ic.idx_File_PSD_t2 = find(abs(MDB_File.file_PSD_time_list-(time(idx_Det_TempRamp_end(i_c))+R_ic.Ramp_t_ext)) == min(abs(MDB_File.file_PSD_time_list-(time(idx_Det_TempRamp_end(i_c))+R_ic.Ramp_t_ext))),1,'first');
    if R_ic.idx_File_PSD_t2+R_ic.num_idx_ext > size(MDB_File.file_PSD_time_list,1)
        R_ic.idx_File_PSD_t2 = size(MDB_File.file_PSD_time_list,1);
    else
        R_ic.idx_File_PSD_t2 = R_ic.idx_File_PSD_t2+R_ic.num_idx_ext;
    end
    R_ic.idx_File_PSD = R_ic.idx_File_PSD_t1:1:R_ic.idx_File_PSD_t2;

    idx_File_t1 = find(abs(MDB_File.file_time_list-(time(idx_Det_TempRamp_init(i_c))+R_ic.Ramp_t_offset)) == min(abs(MDB_File.file_time_list-(time(idx_Det_TempRamp_init(i_c))+R_ic.Ramp_t_offset))),1,'last');
    idx_File_t2 = find(abs(MDB_File.file_time_list-(time(idx_Det_TempRamp_end(i_c))+R_ic.Ramp_t_ext)) == min(abs(MDB_File.file_time_list-(time(idx_Det_TempRamp_end(i_c))+R_ic.Ramp_t_ext))),1,'first');
    if idx_File_t2+R_ic.num_idx_ext > size(MDB_File.file_time_list,1)
        idx_File_t2 = size(MDB_File.file_time_list,1);
    else
        idx_File_t2 = idx_File_t2+R_ic.num_idx_ext;
    end
    idx_File = idx_File_t1:1:idx_File_t2;

    %% Images collected for DataAnlys
    R_ic.numImg = length(R_ic.idx_Img);
    numImg_list(i_c) = R_ic.numImg;
    
    %%% Images within sensible range:
    % timestamp of last image before temperature ramp starts? i.e. ramp not
    % imaged
    if time(idx_Det_TempRamp_init(i_c)) > MDB_Img.Img_time_s(R_ic.idx_Img_t2)
        warning('%s: timestamp of last image before temperature ramp starts! (skipped)',mfilename())   
        continue
    % dataset with less than 51 images?
    elseif R_ic.numImg < 51
        warning('%s: less than 51 images for data anlysis step! (skipped)',mfilename())   
        continue
    end

    
    %% Clip Datatset
    Img_time_c = MDB_Img.Img_time_s(R_ic.idx_Img);
    Img_PSD_10_c = MDB_Img.PSD_10(R_ic.idx_Img);
    Img_PSD_50_c = MDB_Img.PSD_50(R_ic.idx_Img);
    Img_PSD_90_c = MDB_Img.PSD_90(R_ic.idx_Img);
    
    Img_FM_GLLV_c = MDB_Img.FM_GLLV(R_ic.idx_Img);
    Img_FM_GLLV_c_min = min(Img_FM_GLLV_c);
    Img_FM_GLLV_c_max = max(Img_FM_GLLV_c);
    Img_FM_GLLV_c_norm = (Img_FM_GLLV_c-Img_FM_GLLV_c_min)./(Img_FM_GLLV_c_max-Img_FM_GLLV_c_min);

    Img_FM_HELM_c = MDB_Img.FM_HELM_005(R_ic.idx_Img);
    Img_FM_HELM_c_min = min(Img_FM_HELM_c);
    Img_FM_HELM_c_max = max(Img_FM_HELM_c); 
    Img_FM_HELM_c_norm = (Img_FM_HELM_c-Img_FM_HELM_c_min)./(Img_FM_HELM_c_max-Img_FM_HELM_c_min);
    
    Img_FM_GLLV_bk_median = median((MDB_Img.FM_GLLV(Opt_DS.ImgPrc.bk_idx)-Img_FM_GLLV_c_min)./(Img_FM_GLLV_c_max-Img_FM_GLLV_c_min));
    Img_FM_GLLV_bk_std = std((MDB_Img.FM_GLLV(Opt_DS.ImgPrc.bk_idx)-Img_FM_GLLV_c_min)./(Img_FM_GLLV_c_max-Img_FM_GLLV_c_min));
    Img_FM_HELM_bk_median = median((MDB_Img.FM_HELM_005(Opt_DS.ImgPrc.bk_idx)-Img_FM_HELM_c_min)./(Img_FM_HELM_c_max-Img_FM_HELM_c_min));
    Img_FM_HELM_bk_std  = std((MDB_Img.FM_HELM_005(Opt_DS.ImgPrc.bk_idx)-Img_FM_HELM_c_min)./(Img_FM_HELM_c_max-Img_FM_HELM_c_min)); 
    
    if Img_FM_HELM_bk_std < 1e-5
        Img_FM_HELM_bk_std = 1e-5;
    end
    
    Img_Int_median_c = MDB_Img.Int_median(R_ic.idx_Img);

    file_PSD_time_c = MDB_File.file_PSD_time_list(R_ic.idx_File_PSD);
    file_PSD_D10_c = MDB_File.PSD_D10(R_ic.idx_File_PSD);
    file_PSD_D50_c = MDB_File.PSD_D50(R_ic.idx_File_PSD);
    file_PSD_D90_c = MDB_File.PSD_D90(R_ic.idx_File_PSD);

    file_time_c = MDB_File.file_time_list(idx_File);
    file_Transmissivity_c = MDB_File.Transmissivity(idx_File)./100;
    file_Temperature_c = MDB_File.Temperature(idx_File);

    %% Graphs
    if any(~isnan(Img_PSD_10_c))
        fig = figure('units','pixel','position',[100 100 Opt_DS.Graph.width Opt_DS.Graph.height]);
        ax = axes('Parent',fig,...
            'Position',[Opt_DS.Graph.ax_1,Opt_DS.Graph.ax_2,Opt_DS.Graph.ax_1_width,Opt_DS.Graph.ax_2_height]);
        hold(ax,'on');
        box(ax,'on');
        ax.XAxis.Exponent = 0;
        xtickformat('%.0f')
        ax.YColor = Opt_DS.Graph.Color_Seq_grey{1};
        hold on
        
        hold on
        plot(Img_time_c,Img_PSD_10_c, ...
            'LineStyle','-','Color','g','DisplayName','PSD_10 (Img)')
        plot(Img_time_c,Img_PSD_50_c, ...
            'LineStyle','-','Color','b','DisplayName','PSD_50 (Img)')
        plot(Img_time_c,Img_PSD_90_c, ...
             'LineStyle','-','Color','r','DisplayName','PSD_90 (Img)')

        plot(file_PSD_time_c,file_PSD_D10_c, ...
            'LineStyle','--','Color','k','DisplayName','PSD_10 (File)')
        plot(file_PSD_time_c,file_PSD_D50_c, ...
            'LineStyle','--','Color','k','DisplayName','PSD_50 (File)')
        plot(file_PSD_time_c,file_PSD_D90_c, ...
            'LineStyle','--','Color','k','DisplayName','PSD_90 (File)')


        axis([min([Img_time_c;file_PSD_time_c]),max([Img_time_c;file_PSD_time_c]),-Inf,Inf])

        lgd = legend();
        lgd.Interpreter = 'none';
        lgd.Location = 'NorthEast';
        xlabel('Time [s]')
        ylabel('Size [µm]')

        print(fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_ImgAnlys_PSD_Trend_i_c_%.0f',Opt_DS.ExpShorthand,i_c)),'-djpeg',Opt_DS.Graph.print_rS)
    end

    %% Nucleation and Growth Analysis - Local Average

    if Opt_DS.PSD_Anlys_On
        numObj_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        numObj = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        area_Mean_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        d_eqSph_Mean_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        PSD_3_d10_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        PSD_3_d25_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        PSD_3_d50_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        PSD_3_d75_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        PSD_3_d90_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        PSD_3_Span_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        PSD_D_43_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);
        PSD_D_43_Std_movAvg = zeros(R_ic.idx_Img_t2-R_ic.idx_Img_t1+1,1);

        numiter = R_ic.idx_Img_t2 - R_ic.idx_Img_t1;
        k_Fig_Print_On = R_ic.idx_Img_t1+round(0:numiter/10:numiter);
        for k = R_ic.idx_Img_t1:R_ic.idx_Img_t2
            k_idx = k - R_ic.idx_Img_t1+1;
            idx = R_ic.idx_Img(k_idx)-Opt_DS.Img_ObjPara.l_movAvg:R_ic.idx_Img(k_idx)+Opt_DS.Img_ObjPara.l_movAvg;

            area = MDB_SiPa.area(ismember(MDB_SiPa.Img_PrcIdx,idx));
            area_movAvg = area(~isnan(area));
            area_Mean_movAvg(k_idx) = mean(area_movAvg);

            d_eqSph_movAvg = (area_movAvg./pi()).^(0.5).*2;
            d_eqSph_Mean_movAvg(k_idx) = mean(d_eqSph_movAvg);

            numObj(k_idx) = size(d_eqSph_movAvg,1);
            numImg = length(idx);
            numObj_movAvg(k_idx) = numObj(k_idx)/numImg;
            if numObj(k_idx) > Opt.Img_ObjPara.PSD_calc_numThres
                Opt_PSD_Calc = Opt.Img_ObjPara;

                if any(k == k_Fig_Print_On) && Opt_PSD_Calc.Fig_Print_On
                    Opt_PSD_Calc.Fig_Print_On = true;
                    Opt_PSD_Calc.Graph_name = sprintf('%s_ImgAnlys_PSD_Calc_i_c_%.0f_movAvg_k_%.0f',Opt_DS.ExpShorthand,i_c,k);
                else
                    Opt_PSD_Calc.Fig_Print_On = false;
                end
                
                addpath(Opt.path_FileExchange)
                [PSD_3_d25_movAvg(k_idx),PSD_3_d50_movAvg(k_idx),PSD_3_d75_movAvg(k_idx),PSD_3_Span_movAvg(k_idx),R] = Crystalline_PSD_Calc(d_eqSph_movAvg,Opt_PSD_Calc); 
                PSD_3_d10_movAvg(k_idx) = R.PSD_3_d10;
                PSD_3_d90_movAvg(k_idx) = R.PSD_3_d90;
                PSD_D_43_movAvg(k_idx) = R.D_43;
                PSD_D_43_Std_movAvg(k_idx) = R.D_43_StdDev;
            end
        end
    end
    close all

    %% Nucleation and Growth estimation

    % % % % % % % % % % % % % % % % % % % % % 
    % Nucleation Detection:        
    R_ic.Nuc_Time_File_offset = idx_Det_TempRamp_end(i_c);
    R_ic.Nuc_Time_Img_offset = idx_Det_TempRamp_end(i_c);

    Img_FM_GLLV_c_norm_sm = hampel(Img_FM_GLLV_c_norm,Opt_DS.NucDet.k_hampel,Opt_DS.NucDet.k_sigma);
    Img_FM_GLLV_c_norm_sm = sgolayfilt(Img_FM_GLLV_c_norm_sm, 1, Opt_DS.NucDet.k_sgolayfilt);
   
    Img_FM_HELM_c_norm_sm = hampel(Img_FM_HELM_c_norm,Opt_DS.NucDet.k_hampel,Opt_DS.NucDet.k_sigma);
    Img_FM_HELM_c_norm_sm = sgolayfilt(Img_FM_HELM_c_norm_sm, 1, Opt_DS.NucDet.k_sgolayfilt);
    
    file_Transmissivity_c_sm = hampel(file_Transmissivity_c,Opt_DS.NucDet.k_hampel_Transmissivity,Opt_DS.NucDet.k_sigma_Transmissivity);
    file_Transmissivity_c_sm = sgolayfilt(file_Transmissivity_c_sm, 1, Opt_DS.NucDet.k_sgolayfilt_Transmissivity);  
    
    if CoolDet_list(i_c)        
        % Check if:
        % (1) clear at start (Img Feature values(1:15)) and 
        % (2) cloudy at end (Img Feature values(end-15:end))
        if (any(Img_FM_HELM_c_norm_sm(end-15:end) > Img_FM_HELM_bk_median + 4*Img_FM_HELM_bk_std) || ...
                any(Img_FM_HELM_c_norm_sm(end-15:end) < Img_FM_HELM_bk_median - 4*Img_FM_HELM_bk_std)) && ...
                (any(Img_FM_HELM_c_norm_sm(1:15) < Img_FM_HELM_bk_median + 4*Img_FM_HELM_bk_std) && ...
                any(Img_FM_HELM_c_norm_sm(1:15) > Img_FM_HELM_bk_median - 4*Img_FM_HELM_bk_std))
            
            R_ic.idx_Img_FM_GLLV_Nuc = find((MDB_Img.Int_median(R_ic.idx_Img) > Opt_DS.ImgPrc.Img_Int_median_thres)& ...
                ((Img_FM_GLLV_c_norm_sm > Img_FM_GLLV_bk_median + 4*Img_FM_GLLV_bk_std)|(Img_FM_GLLV_c_norm_sm < Img_FM_GLLV_bk_median - 4*Img_FM_GLLV_bk_std)),1,'first');
            
            R_ic.idx_Img_FM_HELM_Nuc = find((MDB_Img.Int_median(R_ic.idx_Img) > Opt_DS.ImgPrc.Img_Int_median_thres)& ...
                ((Img_FM_HELM_c_norm_sm > Img_FM_HELM_bk_median + 4*Img_FM_HELM_bk_std)|(Img_FM_HELM_c_norm_sm < Img_FM_HELM_bk_median - 4*Img_FM_HELM_bk_std)),1,'first');
            
            if ~isempty(R_ic.idx_Img_FM_HELM_Nuc)
                R_ic.idx_Img_Nuc = R_ic.idx_Img_FM_HELM_Nuc;
            else
                R_ic.idx_Img_Nuc = [];
            end
        else
            R_ic.idx_Img_FM_GLLV_Nuc = [];
            R_ic.idx_Img_FM_HELM_Nuc = [];
            R_ic.idx_Img_Nuc = [];
        end
        
        % clear at start (> 0.999)
        if (median(file_Transmissivity_c_sm(1:15)) > 0.999)
            R_ic.idx_file_Transmissivity_Nuc = find((file_Transmissivity_c_sm < 0.999) ,1,'first');
        else
            R_ic.idx_file_Transmissivity_Nuc = [];
        end
    elseif HeatDet_list(i_c)
        % Equal to cloud point detection but with reversed direction
        Img_FM_HELM_c_norm_sm_rev = flip(Img_FM_HELM_c_norm_sm);
        Img_FM_GLLV_c_norm_sm_rev = flip(Img_FM_GLLV_c_norm_sm);
        Img_Int_median_c_rev = flip(Img_Int_median_c);
        
        % Check if:
        % (1) cloudy at start (Img Feature values(1:15)) and 
        % (2) clear at end (Img Feature values rev(1:15))
        if (any(Img_FM_HELM_c_norm_sm(1:15) > Img_FM_HELM_bk_median + 4*Img_FM_HELM_bk_std) || ...
                any(Img_FM_HELM_c_norm_sm(1:15) < Img_FM_HELM_bk_median - 4*Img_FM_HELM_bk_std)) && ...
                (any(Img_FM_HELM_c_norm_sm_rev(1:15) < Img_FM_HELM_bk_median + 4*Img_FM_HELM_bk_std) && ...
                any(Img_FM_HELM_c_norm_sm_rev(1:15) > Img_FM_HELM_bk_median - 4*Img_FM_HELM_bk_std))

            R_ic.idx_Img_FM_GLLV_Nuc = length(Img_FM_GLLV_c_norm_sm_rev) - find((Img_Int_median_c_rev > Opt_DS.ImgPrc.Img_Int_median_thres)& ...
                ((Img_FM_GLLV_c_norm_sm_rev > Img_FM_GLLV_bk_median + 4*Img_FM_GLLV_bk_std)|(Img_FM_GLLV_c_norm_sm_rev < Img_FM_GLLV_bk_median - 4*Img_FM_GLLV_bk_std)),1,'first');
            
            R_ic.idx_Img_FM_HELM_Nuc = length(Img_FM_HELM_c_norm_sm_rev) - find((Img_Int_median_c_rev > Opt_DS.ImgPrc.Img_Int_median_thres)& ...
                ((Img_FM_HELM_c_norm_sm_rev > Img_FM_HELM_bk_median + 4*Img_FM_HELM_bk_std)|(Img_FM_HELM_c_norm_sm_rev < Img_FM_HELM_bk_median - 4*Img_FM_HELM_bk_std)),1,'first');
            
            if ~isempty(R_ic.idx_Img_FM_HELM_Nuc)
                R_ic.idx_Img_Nuc = R_ic.idx_Img_FM_HELM_Nuc;
            else
                R_ic.idx_Img_Nuc = [];
            end
        else
            R_ic.idx_Img_FM_GLLV_Nuc = [];
            R_ic.idx_Img_FM_HELM_Nuc = [];
            R_ic.idx_Img_Nuc = [];
        end

        % clear at start (< 0.999)
        if (median(file_Transmissivity_c_sm(1:15)) < 0.999)
            R_ic.idx_file_Transmissivity_Nuc = find(file_Transmissivity_c_sm > 0.999 ,1,'first');
        else
            R_ic.idx_file_Transmissivity_Nuc = [];
        end
    end
    
    
    %% Save Nucleation Descriptors
    if ~isempty(R_ic.idx_file_Transmissivity_Nuc)
        Nuc_Time_abs_File_list(i_c) = file_time_c(R_ic.idx_file_Transmissivity_Nuc);
        Nuc_Time_Ind_File_list(i_c) = file_time_c(R_ic.idx_file_Transmissivity_Nuc) - R_ic.Nuc_Time_File_offset;
        Nuc_Time_Ind_File_Err_list(i_c) = nan;
        Nuc_Temp_File_list(i_c) = file_Temperature_c(R_ic.idx_file_Transmissivity_Nuc); % degC
        [~,R_ic.idx_Img_Transmissivity_Nuc] = min(abs(file_time_c(R_ic.idx_file_Transmissivity_Nuc)-Img_time_c));
        Nuc_Img_name_File_list{i_c} = MDB_Img.Img_Name{R_ic.idx_Img_t1+R_ic.idx_Img_Transmissivity_Nuc};
        if Nuc_Time_abs_File_list(i_c) > time_File_grad_end(i_c) 
            Nuc_File_IsoThermal_list(i_c) = true;
        else
            Nuc_File_IsoThermal_list(i_c) = false;
        end
    else
        Nuc_Time_abs_File_list(i_c) = nan;
        Nuc_Time_Ind_File_list(i_c) = nan;
        Nuc_Time_Ind_File_Err_list(i_c) = nan;
        Nuc_Temp_File_list(i_c) = nan;
        Nuc_File_IsoThermal_list(i_c) = nan;
    end

    if ~isempty(R_ic.idx_Img_Nuc)
        Nuc_Time_abs_Img_list(i_c) = Img_time_c(R_ic.idx_Img_Nuc);
        Nuc_Time_Ind_Img_list(i_c) = Img_time_c(R_ic.idx_Img_Nuc) - R_ic.Nuc_Time_Img_offset;
        Nuc_Time_Ind_Img_Err_list(i_c) = nan;
        [~,R_ic.idx_file_Nuc] = min(abs(file_time_c-Img_time_c(R_ic.idx_Img_Nuc)));
        Nuc_Temp_Img_list(i_c) = file_Temperature_c(R_ic.idx_file_Nuc); % degC
        Nuc_Img_name_Img_list{i_c} = MDB_Img.Img_Name{R_ic.idx_Img_t1+R_ic.idx_Img_Nuc};
        if Nuc_Time_abs_Img_list(i_c) > time_File_grad_end(i_c)
            Nuc_Img_IsoThermal_list(i_c) = true;
        else
            Nuc_Img_IsoThermal_list(i_c) = false;
        end 
    else
        Nuc_Time_abs_Img_list(i_c) = nan;
        Nuc_Time_Ind_Img_list(i_c) = nan;
        Nuc_Time_Ind_Img_Err_list(i_c) = nan;
        Nuc_Temp_Img_list(i_c) = nan;
        Nuc_Img_IsoThermal_list(i_c) = nan;
    end

    
    %% Graphs
    % % % % % % % % % % % % % % % % % % % % % 
    % Graph Nucleation Detection

    fig = figure('units','pixel','position',[100 100 Opt_DS.Graph.width Opt_DS.Graph.height]);
    ax = axes('Parent',fig,...
        'Position',[Opt_DS.Graph.ax_1,Opt_DS.Graph.ax_2,Opt_DS.Graph.ax_1_width,Opt_DS.Graph.ax_2_height]);
    hold(ax,'on');
    box(ax,'on');
    ax.XAxis.Exponent = 0;
    xtickformat('%.0f')

    plot(Img_time_c,Img_FM_HELM_c_norm,'m-','DisplayName','ImgFeat')
    plot(Img_time_c,Img_FM_HELM_c_norm_sm,'k--','DisplayName','ImgFeat (sm)')

    plot(file_time_c,file_Transmissivity_c,'g-','DisplayName','Transmissivity')
    plot(file_time_c,file_Transmissivity_c_sm,'k--','DisplayName','Transmissivity (sm)')

    % Mark detected nucleation point   
    l1 = line([Nuc_Time_abs_Img_list(i_c) Nuc_Time_abs_Img_list(i_c)],get(ax,'YLim'));
    set(l1,'LineStyle','--','Color','b','LineWidth',1.5,'DisplayName','tNuc ImgFeat')
    l2 = line([Nuc_Time_abs_File_list(i_c) Nuc_Time_abs_File_list(i_c)],get(ax,'YLim'));
    set(l2,'LineStyle',':','Color','b','LineWidth',1.5,'DisplayName','tNuc Transmissivity')

    xlabel('Time [s]')
    ylabel('Intensity')
    ax = gca;
    ax.XLim = [min([file_time_c;Img_time_c]) max([file_time_c;Img_time_c])];
    yyaxis right
    ax.YColor = 'k'; 
    ylabel('Temperature [degC]')
    plot(file_time_c,file_Temperature_c,'r-','DisplayName','Temperature')
    legend('Interpreter','none');

    annotation(fig,'textbox',...
        [ax.Position(1) 0.92 0.8 0.06],...
        'Interpreter', 'Latex', ...
        'String',{sprintf('Nucleation Detetcion: Nuc ImgFeat %.0f s (abs %.0f s, T %.2f degC),\n Nuc Transmissivity %.0f s (abs %.0f s, T %.2f degC)', ...
        Nuc_Time_Ind_Img_list(i_c),Nuc_Time_abs_Img_list(i_c),Nuc_Temp_Img_list(i_c), ...
        Nuc_Time_Ind_File_list(i_c),Nuc_Time_abs_File_list(i_c),Nuc_Temp_File_list(i_c))}, ...
        'FitBoxToText','on');

    print(fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_ImgAnlys_Nuc_Det_%.0f',Opt_DS.ExpShorthand,i_c)),'-djpeg',Opt_DS.Graph.print_rS)

    % % % % % % % % % % % % % % % % % % % % % 
    % Graph Crystal Growth (PSD)
    if Opt_DS.PSD_Anlys_On && any(PSD_3_d25_movAvg > 0)
        fig = figure('units','pixel','position',[100 100 Opt_DS.Graph.width Opt_DS.Graph.height]);
        ax = axes('Parent',fig,...
            'Position',[Opt_DS.Graph.ax_1,Opt_DS.Graph.ax_2,Opt_DS.Graph.ax_1_width,Opt_DS.Graph.ax_2_height]);
        hold(ax,'on');
        box(ax,'on');
        ax.XAxis.Exponent = 0;
        xtickformat('%.0f')
        hold on

        plot(Img_time_c,PSD_3_d25_movAvg, ...
            'LineStyle','-','Color',Opt_DS.Graph.Color_Seq_grey{2},'DisplayName','PSD_3_25 movAvg (Img)')
        plot(Img_time_c,PSD_3_d50_movAvg, ...
            'LineStyle','-','Color',Opt_DS.Graph.Color_Seq_grey{4},'DisplayName','PSD_3_50 movAvg (Img)')
        plot(Img_time_c,PSD_3_d75_movAvg, ...
            'LineStyle','-','Color',Opt_DS.Graph.Color_Seq_grey{6},'DisplayName','PSD_3_75 movAvg (Img)')
        plot(Img_time_c,PSD_3_d90_movAvg, ...
            'LineStyle','-','Color',Opt_DS.Graph.Color_Seq_grey{8},'DisplayName','PSD_3_90 movAvg (Img)')

        if isfield(R_ic,'yi_PSD_3_d25') && ~isempty(R_ic.yi_PSD_3_d25)
            plot(xi_PSD_3,R_ic.yi_PSD_3_d25, ...
                'LineStyle','-','Color','y','DisplayName','PSD_3_25 Fit (Img)')
            plot(xi_PSD_3,R_ic.yi_PSD_3_d50, ...
                'LineStyle','-','Color','g','DisplayName','PSD_3_50 Fit (Img)')
            plot(xi_PSD_3,R_ic.yi_PSD_3_d75, ...
                 'LineStyle','-','Color','b','DisplayName','PSD_3_75 Fit (Img)')
            plot(xi_PSD_3,R_ic.yi_PSD_3_d90, ...
                 'LineStyle','-','Color','r','DisplayName','PSD_3_90 Fit (Img)')  
        end

        lgd = legend('Autoupdate','off');
        lgd.Interpreter = 'none';
        lgd.Location = 'NorthEast';
        xlabel('Time [s]')
        ylabel('Size [µm]')

        if isfield(R_ic,'yi_PSD_3_d25') && ~isempty(R_ic.yi_PSD_3_d25)
            plot(xi_PSD_3,R_ic.yi_PSD_3_d25+2*R_ic.delta_PSD_3_d25,'y--',xi_PSD_3,R_ic.yi_PSD_3_d25-2*R_ic.delta_PSD_3_d25,'y--')
            plot(xi_PSD_3,R_ic.yi_PSD_3_d50+2*R_ic.delta_PSD_3_d50,'g--',xi_PSD_3,R_ic.yi_PSD_3_d50-2*R_ic.delta_PSD_3_d50,'g--')
            plot(xi_PSD_3,R_ic.yi_PSD_3_d75+2*R_ic.delta_PSD_3_d75,'b--',xi_PSD_3,R_ic.yi_PSD_3_d75-2*R_ic.delta_PSD_3_d75,'b--')
            plot(xi_PSD_3,R_ic.yi_PSD_3_d90+2*R_ic.delta_PSD_3_d90,'r--',xi_PSD_3,R_ic.yi_PSD_3_d90-2*R_ic.delta_PSD_3_d90,'r--')

            annotation(fig,'textbox',...
                 [Opt_DS.Graph.ax_1+0.01 Opt_DS.Graph.ax_2_height-0.03 0.10 0.15],...
                'String',{sprintf(['p_PSD_3_d25 = %.2e +/- %.2e um/s\n' ...
                'p_PSD_3_d50 = %.2e +/- %.2e um/s\n' ...
                'p_PSD_3_d75 = %.2e +/- %.2e um/s\n' ...
                'p_PSD_3_d90 = %.2e +/- %.2e um/s'], ...
                R_ic.p_PSD_3_d25(1),mean(R_ic.delta_PSD_3_d25),R_ic.p_PSD_3_d50(1),mean(R_ic.delta_PSD_3_d50), ...
                R_ic.p_PSD_3_d75(1),mean(R_ic.delta_PSD_3_d75),R_ic.p_PSD_3_d90(1),mean(R_ic.delta_PSD_3_d90))},...
                'FitBoxToText','on','Fontsize',9,'Interpreter','none','BackgroundColor','w');
        else
            annotation(fig,'textbox',...
                 [Opt_DS.Graph.ax_1+0.01 Opt_DS.Graph.ax_2_height-0.03 0.10 0.15],...
                'String',{sprintf(['p_PSD_3_d25 = 0 +/- 0 um/s\n' ...
                'p_PSD_3_d50 = 0 +/- 0 um/s\n' ...
                'p_PSD_3_d75 = 0 +/- 0 um/s\n' ...
                'p_PSD_3_d90 = 0 +/- 0 um/s'])},...
                'FitBoxToText','on','Fontsize',9,'Interpreter','none','BackgroundColor','w');

        end

        print(fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_ImgAnlys_Growth_PSD_i_c_%.0f_movAvg',Opt_DS.ExpShorthand,i_c)),'-djpeg',Opt_DS.Graph.print_rS)

        if isfield(R_ic,'yi_PSD_3_d25') && ~isempty(R_ic.yi_PSD_3_d25)
            yi_PSD_max = max([max(R_ic.yi_PSD_3_d25),max(R_ic.yi_PSD_3_d50),max(R_ic.yi_PSD_3_d75),max(R_ic.yi_PSD_3_d90)])*1.1;
            yi_PSD_min = min([min(R_ic.yi_PSD_3_d25),min(R_ic.yi_PSD_3_d50),min(R_ic.yi_PSD_3_d75),min(R_ic.yi_PSD_3_d90)])*0.9;
            if yi_PSD_max > 0
                axis([xi_PSD_3(1),xi_PSD_3(end),yi_PSD_min,yi_PSD_max]) %#ok<FNCOLND>
            end

            print(fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_ImgAnlys_Growth_PSD_i_c_%.0f_movAvg_rS',Opt_DS.ExpShorthand,i_c)),'-djpeg',Opt_DS.Graph.print_rS)                
        end
    end
end

%% Assemble MDB_Data Table Structure
MDB_Data = table();

MDB_Data.HeatDet = HeatDet_list;
MDB_Data.CoolDet = CoolDet_list;
MDB_Data.HeatRate = HeatRate_list;
MDB_Data.T_Const_1 = T_Const_1_list;
MDB_Data.T_Const_2 = T_Const_2_list;
MDB_Data.time_File_grad_init = time_File_grad_init;
MDB_Data.time_File_grad_end = time_File_grad_end;

MDB_Data.Nuc_Time_abs_File = Nuc_Time_abs_File_list;
MDB_Data.Nuc_Time_Ind_File = Nuc_Time_Ind_File_list;
MDB_Data.Nuc_Time_Ind_File_Err = Nuc_Time_Ind_File_Err_list;

MDB_Data.Nuc_Time_abs_Img = Nuc_Time_abs_Img_list;
MDB_Data.Nuc_Time_Ind_Img = Nuc_Time_Ind_Img_list;
MDB_Data.Nuc_Time_Ind_Img_Err = Nuc_Time_Ind_Img_Err_list;

MDB_Data.Nuc_Img_name_File = Nuc_Img_name_File_list;
MDB_Data.Nuc_Img_name_Img = Nuc_Img_name_Img_list;
MDB_Data.numImg = numImg_list;

MDB_Data.Nuc_Temp_Img = Nuc_Temp_Img_list;
MDB_Data.Nuc_Temp_File = Nuc_Temp_File_list;
MDB_Data.Nuc_Img_IsoThermal = Nuc_Img_IsoThermal_list;
MDB_Data.Nuc_File_IsoThermal = Nuc_File_IsoThermal_list;

MDB_Data.Growth_PSD_3_25 = Growth_PSD_3_25_list;
MDB_Data.Growth_PSD_3_25_Err = Growth_PSD_3_25_Err_list;
MDB_Data.Growth_PSD_3_50 = Growth_PSD_3_50_list;
MDB_Data.Growth_PSD_3_50_Err = Growth_PSD_3_50_Err_list;
MDB_Data.Growth_PSD_3_75 = Growth_PSD_3_75_list;
MDB_Data.Growth_PSD_3_75_Err = Growth_PSD_3_75_Err_list;
MDB_Data.Growth_PSD_3_90 = Growth_PSD_3_90_list;
MDB_Data.Growth_PSD_3_90_Err = Growth_PSD_3_90_Err_list;
MDB_Data.PLS_logic = PLS_logic;

if ~isempty(MDB_Data)
    % Save MDB_Data
    writetable(MDB_Data,fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_MDB_Data_%s.csv',Opt_DS.ExpShorthand,datestr(now, 'yyyymmdd'))))
end

%% Additonal DataAnlys result structure
R_DataAnlys.Temp_mode_Img = Temp_mode_Img;
R_DataAnlys.PLS_logic_Img = PLS_logic_Img; 

R_DataAnlys.Temp_mode_File = Temp_mode_File;

%% End of Script
fprintf('%s - %s (%s) COMPLETE (Elapsed time: %.0f sec)\n',Opt.ProjectShorthand,Opt_DS.ExpShorthand,mfilename(),toc(Opt.tic))
