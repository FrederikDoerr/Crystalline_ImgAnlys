%------------------------------------------------------------------------------------------------
% Application: Crystalline Image and Data Analysis Framework
% Main Routine
% 
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2020b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/FrederikDoerr


fclose('all');
clear all %#ok<CLALL>
close all
clc


set(0,'DefaultFigureVisible','on');
Opt.version = 0.5; %


%% Setup
Opt.ProjectShorthand = 'FrDo_Crystalline_Anlys';

% Main folder location
path = matlab.desktop.editor.getActiveFilename;
[Opt.path_ProjectFolder,~,~] = fileparts(path);

Opt.path_MatFolder_lc = fullfile('C:\Users',getenv('USERNAME'),'Documents\MATLAB');
Opt.path_ParamFolder = fullfile(Opt.path_ProjectFolder,'_Parameter_Data');
Opt.path_FileExchange = fullfile(Opt.path_ProjectFolder ,'_FileExchange');

Opt.name_Masterfile = 'CryExp_Masterfile.xlsx';
Opt.MainPrc_T_Crystalline = readtable(fullfile(Opt.path_ProjectFolder,Opt.name_Masterfile),'Sheet','Crystalline');
Opt.ImgPrc_T_Crystalline = readtable(fullfile(Opt.path_ProjectFolder,Opt.name_Masterfile),'Sheet','Prc_ImgAnlys');

Opt.name_ExportFolder = sprintf('%s_Export_%s',Opt.ProjectShorthand,datestr(now,'yyyy-mm-dd'));
Opt.path_ExportFolder = fullfile(Opt.path_ProjectFolder,Opt.name_ExportFolder);
if ~exist(Opt.path_ExportFolder,'dir')
    mkdir(Opt.path_ExportFolder)
end
cd(Opt.path_ExportFolder)
addpath(Opt.path_ProjectFolder)

% Get local computer ID
addpath(fullfile(Opt.path_FileExchange,'getComputerName'))
Opt.localPC_name = getComputerName();

formatOut = 'yyyy-mm-dd';
datestring = datestr(datetime('now'),formatOut);
Opt.idx_MasterDB = length(dir(fullfile(Opt.path_ExportFolder,sprintf('%s_%s_%s_*.csv',Opt.ProjectShorthand,datestring,Opt.localPC_name))))+1;
Opt.filesName_MasterDB = sprintf('%s_%s_%s_%03.0f.csv',Opt.ProjectShorthand,datestring,Opt.localPC_name,Opt.idx_MasterDB);
Opt.filesPath_MasterDB = fullfile(Opt.path_ExportFolder,Opt.filesName_MasterDB);
Opt.numDS = find(cellfun(@isempty,Opt.MainPrc_T_Crystalline.ExpShorthand),1,'first')-1;
if isempty(Opt.numDS)
    Opt.numDS = size(Opt.MainPrc_T_Crystalline.ExpShorthand,1);
end

[folder, name, ext] = fileparts(which('Crystalline_Main'));
Opt.pwd = folder;

fprintf('%s - Setup complete ---\n',Opt.ProjectShorthand)

%% Loop Dataset
for k_DS = 1:Opt.numDS
    
    Opt.k_DS = k_DS;
    
    MDB_Exp_VariableNames = { ...   
        'ExpShorthand','ELN','OperatorShorthand', ...
        'Src_Data_Path','Src_Img','Src_PrcFile','Dest_Data_Path', ...
        'Info', 'Img_RePrc','Platform','SerialNo','Position', ...
        'CompoundShorthand','SolventShorthand','ChemSys_ID', ...
        'Cmp_Mass','Solv_Mass','Cmp_Conc_w_w', ...
        'time_Prc','time_Prc_ImgAvg', ...
            };

    MDB_Exp = array2table(nan(1,length(MDB_Exp_VariableNames)),'VariableNames',MDB_Exp_VariableNames);  


    clearvars -except Opt
    
    
    % User-defined Metadata  
    Opt_DS.ExpShorthand = Opt.MainPrc_T_Crystalline.ExpShorthand{Opt.k_DS};   
    Opt_DS.Dest_Data_Path = Opt.MainPrc_T_Crystalline.Dest_Data_Path{Opt.k_DS};
    if ~exist(Opt_DS.Dest_Data_Path,'dir')
        Opt_DS.Dest_Data_Path = fullfile(Opt.pwd,Opt_DS.Dest_Data_Path);
    end
    Opt_DS.Src_Data_Path = Opt.MainPrc_T_Crystalline.Src_Data_Path{Opt.k_DS};
    if ~exist(Opt_DS.Src_Data_Path,'dir')
        Opt_DS.Src_Data_Path = fullfile(Opt.pwd,Opt_DS.Src_Data_Path);
    end
    Opt_DS.Src_PrcFile = sprintf('%s.csv',Opt.MainPrc_T_Crystalline.Src_PrcFile{Opt.k_DS});
    Opt_DS.Src_PrcFile_path = fullfile(Opt_DS.Src_Data_Path,Opt_DS.Src_PrcFile); 
    Opt_DS.ImgFolder_name = Opt.MainPrc_T_Crystalline.Src_Img{Opt.k_DS};
    Opt_DS.ImgFolder_path = fullfile(Opt_DS.Src_Data_Path,Opt_DS.ImgFolder_name);
    
    Opt_DS.path_ExportFolder_Exp = fullfile(Opt.path_ExportFolder,Opt_DS.ExpShorthand);
    if ~exist(Opt_DS.path_ExportFolder_Exp,'dir')
       mkdir(Opt_DS.path_ExportFolder_Exp) 
    end
    
    Opt_DS.Img_RePrc = Opt.MainPrc_T_Crystalline.Img_RePrc(Opt.k_DS);
    Opt_DS.Prc_ImgAnlys_REF = Opt.MainPrc_T_Crystalline.Prc_ImgAnlys_REF{Opt.k_DS};

    Opt_DS.OperatorShorthand = Opt.MainPrc_T_Crystalline.OperatorShorthand{Opt.k_DS};
    Opt_DS.ELN = Opt.MainPrc_T_Crystalline.ELN{Opt.k_DS};
    Opt_DS.Platform = Opt.MainPrc_T_Crystalline.Platform{Opt.k_DS};
    Opt_DS.SerialNo = Opt.MainPrc_T_Crystalline.SerialNo{Opt.k_DS};
    Opt_DS.Position = Opt.MainPrc_T_Crystalline.Position{Opt.k_DS};
    
    Opt_DS.Info = Opt.MainPrc_T_Crystalline.Info{Opt.k_DS};
    
    % Load Processing Parameters
    Opt_DS = Crystalline_ParamLoader(Opt_DS,Opt);
    
    % log-file (initialise)
    Opt.tic = tic;
    Opt.date_formatOut = 'yyyy-mm-dd';
    Opt.date = datestr(datetime('now'),Opt.date_formatOut);
    diary(sprintf('LOG_%s_%s_%s.txt',Opt.ProjectShorthand,Opt_DS.ExpShorthand,Opt.date))
    diary on
    
    % log-file
    fprintf('Date: %s \n',Opt.date)
    fprintf('ExpShorthand: %s\n',Opt_DS.ExpShorthand)
    fprintf('OperatorShorthand: %s\n',Opt_DS.OperatorShorthand)
    fprintf('ELN: %s\n',Opt_DS.ELN)
    fprintf('Platform: %s\n',Opt_DS.Platform)
    fprintf('SerialNo: %s\n',Opt_DS.SerialNo)
    fprintf('Position: %s\n',Opt_DS.Position)
    
    fprintf('Src_Data_Path: %s\n',Opt.MainPrc_T_Crystalline.Src_Data_Path{Opt.k_DS})
    fprintf('Src_PrcFile: %s\n',Opt_DS.Src_PrcFile)
    fprintf('Src_Img: %s\n',Opt_DS.ImgFolder_name)
    fprintf('Dest_Data_Path: %s\n',Opt_DS.Dest_Data_Path)
    fprintf('Prc_ImgAnlys_REF: %s\n',Opt_DS.Prc_ImgAnlys_REF)
    fprintf('Info: %s \n',Opt_DS.Info)
    
    %% Update MDB_Exp
    MDB_Exp.ExpShorthand = Opt_DS.ExpShorthand;
    MDB_Exp.ELN = Opt_DS.ELN;
    MDB_Exp.Platform = Opt_DS.Platform;
    MDB_Exp.SerialNo = Opt_DS.SerialNo;
    MDB_Exp.Position = Opt_DS.Position;
    
    MDB_Exp.OperatorShorthand = Opt_DS.OperatorShorthand;
    MDB_Exp.Src_Data_Path = Opt_DS.Src_Data_Path;
    MDB_Exp.Src_PrcFile = Opt_DS.Src_PrcFile;
    MDB_Exp.Src_Img = Opt_DS.ImgFolder_name;
    MDB_Exp.Dest_Data_Path = Opt_DS.Dest_Data_Path;
    MDB_Exp.Info = Opt.MainPrc_T_Crystalline.Info{Opt.k_DS};
    
    MDB_Exp.Img_RePrc = Opt_DS.Img_RePrc;

    
    %% Load CSV File
    fread_opts = detectImportOptions(Opt_DS.Src_PrcFile_path, ...
        'NumHeaderLines', 1);
    fread_opts.VariableNamesLine = 1;
    fread_opts = setvaropts(fread_opts,'Var1','InputFormat','dd/MM/yyyy HH:mm:ss');
    if length(fread_opts.VariableNames) == 12
        fread_opts = setvaropts(fread_opts,'Var12','InputFormat','dd/MM/yyyy HH:mm:ss');
    end

    warning('off','MATLAB:table:ModifiedAndSavedVarnames')
    MDB_File =  readtable(Opt_DS.Src_PrcFile_path,fread_opts);
    
    % Problems reading timestamp
    [warnMsg, warnId] = lastwarn;
    if any(isnat(table2array(MDB_File(:,1)))) || strcmp(warnId,'MATLAB:readtable:AllNaTVariable')
        
        fread_opts = detectImportOptions(Src_PrcFile_path, ...
        'NumHeaderLines', 1);
        fread_opts.VariableNamesLine = 1;
        
        MDB_File =  readtable(Src_PrcFile_path,fread_opts);
        MDB_File.Properties.VariableNames{1} = 'Timestamp';
        MDB_File.Timestamp= datetime(MDB_File.Timestamp,'Format','dd/MM/yyyy HH:mm:ss');
        if length(fread_opts.VariableNames) > 6
            MDB_File.Properties.VariableNames{12} = 'PSD_Timestamp';
            MDB_File.PSD_Timestamp= datetime(MDB_File.PSD_Timestamp,'Format','dd/MM/yyyy HH:mm:ss');
        end
    end
    
    % Rename Column Headers
    MDB_File.Properties.VariableNames{1} = 'Timestamp';
    MDB_File.Properties.VariableNames{2} = 'Temperature';
    MDB_File.Properties.VariableNames{4} = 'Transmissivity';
    MDB_File.Properties.VariableNames{5} = 'StirringSpeed_top';
    MDB_File.Properties.VariableNames{6} = 'StirringSpeed_bottom';
    
    if length(fread_opts.VariableNames) > 6
        MDB_File.Properties.VariableNames{7} = 'Temp_Ctrl_1';
        MDB_File.Properties.VariableNames{8} = 'Temp_Ctrl_2';
        MDB_File.Properties.VariableNames{9} = 'Temp_Ctrl_3';
        MDB_File.Properties.VariableNames{10} = 'Temp_Ctrl_4';
        MDB_File.Properties.VariableNames{11} = 'Spacer';
    else
        MDB_File.Temp_Ctrl_1 = nan(size(MDB_File.PSD_Timestamp,1),1);
        MDB_File.Temp_Ctrl_2 = nan(size(MDB_File.PSD_Timestamp,1),1);
        MDB_File.Temp_Ctrl_3 = nan(size(MDB_File.PSD_Timestamp,1),1);
        MDB_File.Spacer = nan(size(MDB_File.PSD_Timestamp,1),1);
    end
    
    if length(fread_opts.VariableNames) > 11
        
        MDB_File.Properties.VariableNames{12} = 'PSD_Timestamp';
        MDB_File.Properties.VariableNames{13} = 'PSD_D10';
        MDB_File.Properties.VariableNames{14} = 'PSD_D50';
        MDB_File.Properties.VariableNames{15} = 'PSD_D90';
        MDB_File.Properties.VariableNames{16} = 'PSD_D10_Shape';
        MDB_File.Properties.VariableNames{17} = 'PSD_D50_Shape';
        MDB_File.Properties.VariableNames{18} = 'PSD_D90_Shape';
        
        Opt_DS.File.CSV_PSD_Check = true;
    else
        MDB_File.PSD_Timestamp = MDB_File.Timestamp;
        MDB_File.PSD_D10 = zeros(size(MDB_File.PSD_Timestamp,1),1);
        MDB_File.PSD_D50 = zeros(size(MDB_File.PSD_Timestamp,1),1);
        MDB_File.PSD_D90 = zeros(size(MDB_File.PSD_Timestamp,1),1);
        MDB_File.PSD_D10_Shape = zeros(size(MDB_File.PSD_Timestamp,1),1);
        MDB_File.PSD_D50_Shape = zeros(size(MDB_File.PSD_Timestamp,1),1);
        MDB_File.PSD_D90_Shape = zeros(size(MDB_File.PSD_Timestamp,1),1);

        Opt_DS.File.CSV_PSD_Check = false;
    end

    %% File Timestamp processing
    
    Opt_DS.File.file_time_init = sum( ...
        day(MDB_File.Timestamp(1))*24*3600 + ...
        hour(MDB_File.Timestamp(1))*3600 + ...
        minute(MDB_File.Timestamp(1))*60 + ...
        second(MDB_File.Timestamp(1)) ...
        );
    
    MDB_File.file_time_list = ( ...
        day(MDB_File.Timestamp(:))*24*3600 + ...
        hour(MDB_File.Timestamp(:))*3600 + ...
        minute(MDB_File.Timestamp(:))*60 + ...
        second(MDB_File.Timestamp(:)) ...
        ) - Opt_DS.File.file_time_init;
    
    Opt_DS.File.file_PSD_time_init = sum( ...
        day(MDB_File.PSD_Timestamp(1))*24*3600 + ...
        hour(MDB_File.PSD_Timestamp(1))*3600 + ...
        minute(MDB_File.PSD_Timestamp(1))*60 + ...
        second(MDB_File.PSD_Timestamp(1)) ...
        );

    MDB_File.file_PSD_time_list = ( ...
        day(MDB_File.PSD_Timestamp(:))*24*3600 + ...
        hour(MDB_File.PSD_Timestamp(:))*3600 + ...
        minute(MDB_File.PSD_Timestamp(:))*60 + ...
        second(MDB_File.PSD_Timestamp(:)) ...
        ) - Opt_DS.File.file_PSD_time_init;
    
    
    %% Reload processed image data
    if ~Opt_DS.Img_RePrc
    
        % Check whether datafile exists (located with original image data)
        if exist(fullfile(Opt_DS.Src_Data_Path,sprintf('%s_MDB.mat',Opt_DS.ExpShorthand)),'file')
            R = load(fullfile(Opt_DS.Src_Data_Path,sprintf('%s_MDB.mat',Opt_DS.ExpShorthand)));
            
            MDB_Img = R.MDB_Img;
            MDB_File = R.MDB_File;
            MDB_Exp = R.MDB_Exp;
            MDB_SiPa = R.MDB_SiPa;
            Opt_DS.ImgPrc.bk_idx = R.Opt_DS.ImgPrc.bk_idx;
            fprintf('--- %s - %s - Processed Image Data Reloaded ---\n\n',Opt.ProjectShorthand,Opt_DS.ExpShorthand)
        else
            warning('%s - %s - Processed Image Data Not Found!\n',Opt.ProjectShorthand,Opt_DS.ExpShorthand)
            Opt_DS.Img_RePrc = true;
        end
    end
    
    if Opt_DS.Img_RePrc
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        addpath(Opt.path_ProjectFolder)
        set(0,'DefaultFigureVisible','off');
        [MDB_Img, MDB_SiPa, Opt_DS] = Crystalline_ImgAnlys(Opt,Opt_DS,MDB_File);
        set(0,'DefaultFigureVisible','on');
        
        %% Update MDB_Exp
        MDB_Exp.time_Prc_Img = toc(Opt.tic);
        MDB_Exp.time_Prc_ImgAvg = toc(Opt.tic)/height(MDB_Img);
    end
    

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %% Data Post-Analysis
    [MDB_Data,R_DataAnlys] = Crystalline_DataAnlys(Opt,Opt_DS,MDB_File,MDB_Img);
    close all
    
    %% Update MDB_Exp
    MDB_Exp.time_Prc_Data = toc(Opt.tic);
    
    %% Update MDB_Img
    MDB_Img.Temp_mode_Img = R_DataAnlys.Temp_mode_Img;
    MDB_Img.PLS_logic_Img = R_DataAnlys.PLS_logic_Img; 
    
    %% Update MDB_File
    MDB_File.Temp_mode_File = R_DataAnlys.Temp_mode_File;
    
    %% Save MDB Files for re-loading    
    save(fullfile(Opt.path_ExportFolder,sprintf('%s_MDB.mat',Opt_DS.ExpShorthand)),'MDB_Data','MDB_Img','MDB_SiPa','MDB_File','MDB_Exp','Opt','Opt_DS');
    
    %% Export MDB_Exp to csv file    
    addpath(fullfile(Opt.path_FileExchange,'struct2csv'))
    struct2csv_FrDo(MDB_Exp,fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_MDB_Exp_%s.csv',Opt_DS.ExpShorthand,datestr(now, 'yyyymmdd'))))
    
    %% Export MDB_Img to csv file
    writetable(MDB_Img,fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_MDB_Img_%s.csv',Opt_DS.ExpShorthand,datestr(now, 'yyyymmdd'))))

    %% Export MDB_File to csv file
    writetable(MDB_File,fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_MDB_File_%s.csv',Opt_DS.ExpShorthand,datestr(now, 'yyyymmdd'))))
    
    %% Export MDB_Data to csv file
    writetable(MDB_Data,fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_MDB_Data_%s.csv',Opt_DS.ExpShorthand,datestr(now, 'yyyymmdd'))))
    
    %% End Sequence
    fprintf('Processed Crystalline Datasets: %.0f / %.0f (%s)\n',Opt.k_DS,Opt.numDS,Opt_DS.ExpShorthand)
    pause(5)
    close all
    diary off
end
