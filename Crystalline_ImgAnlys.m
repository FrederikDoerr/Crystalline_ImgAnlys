%------------------------------------------------------------------------------------------------
% Application: Crystalline Image and Data Analysis Framework
% Subroutine: Crystalline_ImgAnlys
% 
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2020b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/FrederikDoerr


function[MDB_Img, MDB_SiPa, Opt_DS] = Crystalline_ImgAnlys(Opt,Opt_DS,MDB_File)

fprintf('%s - %s (%s) INITIATED (Elapsed time: %.0f sec)\n',Opt.ProjectShorthand,Opt_DS.ExpShorthand,mfilename(),toc(Opt.tic))

%% Crystalline Image processing

if Opt_DS.ImgPrc.Fig_Print_On
    if ~exist(Opt_DS.ImgPrc.path_ExportFolder_DS_Img,'dir') && Opt_DS.ImgPrc.Fig_Print_On
        mkdir(Opt_DS.ImgPrc.path_ExportFolder_DS_Img)
    end
end

if Opt_DS.ImgPrc.BW_Print_On
    if ~exist(Opt_DS.ImgPrc.path_ExportFolder_DS_bw,'dir') && Opt_DS.ImgPrc.BW_Print_On
        mkdir(Opt_DS.ImgPrc.path_ExportFolder_DS_bw)
    end
end
    
if Opt_DS.ImgPrc.BW_Prc_Ext_Print_On
    if ~exist(Opt_DS.ImgPrc.path_ExportFolder_DS_bw_Ext,'dir') && Opt_DS.ImgPrc.BW_Prc_Ext_Print_On
        mkdir(Opt_DS.ImgPrc.path_ExportFolder_DS_bw_Ext)
    end
end

Img_filelist = dir(fullfile(Opt_DS.ImgFolder_path,sprintf('*.%s',Opt_DS.ImgPrc.Img_Format)));
Opt_DS.ImgPrc.numImg = size(Img_filelist,1);



%% Analyse image timestamps and assess time coverage
Img_time_abs_list = nan(Opt_DS.ImgPrc.numImg,1);
Img_time_list = nan(Opt_DS.ImgPrc.numImg,1);
Img_time_logic = zeros(size(MDB_File.file_time_list),'logical');

for k = 1:Opt_DS.ImgPrc.numImg
    img_name_F = Img_filelist(k).name;
    
    % Test wether time stamp changes number of digits
    idx_bk_img_name = find(isspace(img_name_F),1,'last');
    idx_bk_img_name_m = max(strfind(img_name_F(1:end-4),'m'));
    idx_bk_img_name_s = max(strfind(img_name_F(1:end-4),'s'));

    Img_time_abs_list(k) = str2double(img_name_F(idx_bk_img_name+1:idx_bk_img_name_m-1))*60+ ...
    str2double(img_name_F(idx_bk_img_name_m+1:idx_bk_img_name_s-1));

    Img_time_list(k) = Img_time_abs_list(k) - Opt_DS.tstamp_offset;
end

% Update image sorting related to timestamp
[Img_time_list,idx_Img_time_sort] = sort(Img_time_list);
Img_filelist = Img_filelist(idx_Img_time_sort);

for k = 1:Opt_DS.ImgPrc.numImg
    % Match against recorded temperature profile
    Img_time_logic(abs(MDB_File.file_time_list - Img_time_list(k)) == min(abs(MDB_File.file_time_list - Img_time_list(k)))) = 1;
end

% Graph
fig = figure('units','pixel','position',[100 100 Opt_DS.Graph.width Opt_DS.Graph.height]);
ax = axes('Parent',fig,...
    'Position',[Opt_DS.Graph.ax_1,Opt_DS.Graph.ax_2,Opt_DS.Graph.ax_1_width,Opt_DS.Graph.ax_2_height]);
hold(ax,'on');
box(ax,'on');

yyaxis left
ax.YColor = 'r';
plot(MDB_File.file_time_list,MDB_File.Temperature,'r-','DisplayName','Temperature')

xlabel('Time [s]')
ylabel('Temperature')
ax = gca;
ax.XLim = [min(MDB_File.file_time_list) max(MDB_File.file_time_list)];

yyaxis right
scatter(MDB_File.file_time_list,Img_time_logic,'ks','DisplayName','Image Coverage')
ax.YColor = 'k';
ax.YLim = [0 1.1];
ylabel('ImgMatched')

print(fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_File_Temp_ImgageCov',Opt_DS.ExpShorthand)),'-djpeg',Opt_DS.Graph.print_rS)

%% Background Image
Img_bk_filelist = dir(fullfile(Opt_DS.ImgFolder_path,'bk',sprintf('*.%s',Opt_DS.ImgPrc.Img_Format)));
if ~isempty(Img_bk_filelist)
    Opt_DS.ImgPrc.bk_idx = nan(length(Img_bk_filelist),1);
    Opt_DS.ImgPrc.bk_img_name = cell(length(Img_bk_filelist),1);
    
    img_name = Img_filelist(1).name;
    I = imread(fullfile(Img_filelist(1).folder,img_name));
    
    % match against timeseries and average background image
    I_bk = zeros(size(I,1),size(I,2),length(Img_bk_filelist),'double');
    for i_bk = 1:length(Img_bk_filelist)
        Opt_DS.ImgPrc.bk_idx(i_bk) = find(strcmp({Img_filelist(:).name},Img_bk_filelist(i_bk).name));
        Opt_DS.ImgPrc.bk_img_name{i_bk} = Img_filelist(Opt_DS.ImgPrc.bk_idx(i_bk)).name;
        I_bk(:,:,i_bk) = imread(fullfile(Img_filelist(Opt_DS.ImgPrc.bk_idx(i_bk)).folder,Opt_DS.ImgPrc.bk_img_name{i_bk}));
    end
    
    Opt_DS.ImgPrc.bk_img = uint8(mean(I_bk,3));
else
    warning('%s: no background image found (subfolder bk)!',mfilename())    
    Opt_DS.ImgPrc.bk_idx = nan(length(Img_bk_filelist),1);
    Opt_DS.ImgPrc.bk_img_name = cell(length(Img_bk_filelist),1);
    
    img_name = Img_filelist(1).name;
    I = imread(fullfile(Img_filelist(1).folder,img_name));
    
    % Find max temperature range (assuming full dissolution) and match
    % against images
    n_Img = 15;
    Opt_DS.ImgPrc.bk_idx = zeros(n_Img,1);
    idx_maxTemp = find(MDB_File.Temperature > max(MDB_File.Temperature)-2);
        
    [~,idx_Img_matched] = maxk(MDB_File.file_time_list(idx_maxTemp).*Img_time_logic(idx_maxTemp),n_Img);
    Opt_DS.ImgPrc.bk_idx = nan(length(idx_Img_matched),1);
    for k =1:length(idx_Img_matched)
        Opt_DS.ImgPrc.bk_idx(k) = find(abs(Img_time_list - MDB_File.file_time_list(idx_Img_matched(k))) == min(abs(Img_time_list - MDB_File.file_time_list(idx_Img_matched(k)))),1,'first');
    end
    
    if ~exist(fullfile(Opt_DS.ImgFolder_path,'bk'),'dir')
        mkdir(fullfile(Opt_DS.ImgFolder_path,'bk'))
    end
    % average background image
    I_bk = zeros(size(I,1),size(I,2),length(Opt_DS.ImgPrc.bk_idx),'double');
    for i_bk = 1:length(Opt_DS.ImgPrc.bk_idx)
        Opt_DS.ImgPrc.bk_img_name{i_bk} = Img_filelist(Opt_DS.ImgPrc.bk_idx(i_bk)).name;
        I_bk(:,:,i_bk) = imread(fullfile(Img_filelist(Opt_DS.ImgPrc.bk_idx(i_bk)).folder,Opt_DS.ImgPrc.bk_img_name{i_bk}));
        imwrite(uint8(I_bk(:,:,i_bk)),fullfile(Opt_DS.ImgFolder_path,'bk',Opt_DS.ImgPrc.bk_img_name{i_bk}),'bmp')
    end
    Opt_DS.ImgPrc.bk_img = uint8(mean(I_bk,3));
end

Opt_DS.ImgPrc.bk_img = Opt_DS.ImgPrc.bk_img(Opt_DS.ImgPrc.ymin:Opt_DS.ImgPrc.ymax,Opt_DS.ImgPrc.xmin:Opt_DS.ImgPrc.xmax);
Opt_DS.ImgPrc.bk_bw = ones(Opt_DS.ImgPrc.ymax-Opt_DS.ImgPrc.ymin+1,Opt_DS.ImgPrc.xmax-Opt_DS.ImgPrc.xmin+1,'logical');

Opt_DS.ImgPrc.idx_bk_img_name = find(isspace(Opt_DS.ImgPrc.bk_img_name{1}),1,'last');
Opt_DS.ImgPrc.idx_bk_img_name_m = max(strfind(Opt_DS.ImgPrc.bk_img_name{1}(1:end-4),'m'));
Opt_DS.ImgPrc.idx_bk_img_name_s = max(strfind(Opt_DS.ImgPrc.bk_img_name{1}(1:end-4),'s'));


%% Define Data Structures

Img_Name_list = cell(Opt_DS.ImgPrc.numImg,1);
Img_PrcIdx_list = nan(Opt_DS.ImgPrc.numImg,1);

Img_StirringSpeed_top_list = nan(Opt_DS.ImgPrc.numImg,1);
Img_StirringSpeed_bottom_list = nan(Opt_DS.ImgPrc.numImg,1);
Img_Transmissivity_list = nan(Opt_DS.ImgPrc.numImg,1);
Img_Temperature_list = nan(Opt_DS.ImgPrc.numImg,1);


Opt_ImgPrc = Opt_DS.ImgPrc;

% % % % % % % % % % % % % % % % % % % % % 
% Image features

% Crop image border
[y_ROI_sq_max,x_ROI_sq_max] = size(Opt_DS.ImgPrc.bk_img);
ROI_b = round(min([0.1*y_ROI_sq_max,0.1*x_ROI_sq_max]));
Opt_DS.Img_ImgFeat.x_ROI_sq = ROI_b;
Opt_DS.Img_ImgFeat.y_ROI_sq = ROI_b;
Opt_DS.Img_ImgFeat.dx_ROI_sq = x_ROI_sq_max-2*Opt_DS.Img_ImgFeat.x_ROI_sq;
Opt_DS.Img_ImgFeat.dy_ROI_sq = y_ROI_sq_max-2*Opt_DS.Img_ImgFeat.y_ROI_sq;
Opt_DS.Img_ImgFeat.xcoords = [ ...
    Opt_DS.Img_ImgFeat.x_ROI_sq, ...
    Opt_DS.Img_ImgFeat.x_ROI_sq+Opt_DS.Img_ImgFeat.dx_ROI_sq, ...
    Opt_DS.Img_ImgFeat.x_ROI_sq+Opt_DS.Img_ImgFeat.dx_ROI_sq, ...
    Opt_DS.Img_ImgFeat.x_ROI_sq, ...
    Opt_DS.Img_ImgFeat.x_ROI_sq];
Opt_DS.Img_ImgFeat.ycoords = [ ...
    Opt_DS.Img_ImgFeat.y_ROI_sq, ...
    Opt_DS.Img_ImgFeat.y_ROI_sq, ...
    Opt_DS.Img_ImgFeat.y_ROI_sq+Opt_DS.Img_ImgFeat.dy_ROI_sq, ...
    Opt_DS.Img_ImgFeat.y_ROI_sq+Opt_DS.Img_ImgFeat.dy_ROI_sq, ...
    Opt_DS.Img_ImgFeat.y_ROI_sq];

% Compensate for static image background variance
Opt_DS.Img_ImgFeat.I_bk = Opt_DS.ImgPrc.bk_img( ...
    Opt_DS.Img_ImgFeat.y_ROI_sq:Opt_DS.Img_ImgFeat.y_ROI_sq+Opt_DS.Img_ImgFeat.dy_ROI_sq, ...
    Opt_DS.Img_ImgFeat.x_ROI_sq:Opt_DS.Img_ImgFeat.x_ROI_sq+Opt_DS.Img_ImgFeat.dx_ROI_sq);
Opt_DS.Img_ImgFeat.I_bk_pxVar = imsubtract(Opt_DS.Img_ImgFeat.I_bk, ...
    ones(size(Opt_DS.Img_ImgFeat.I_bk),'uint8').*mean(Opt_DS.Img_ImgFeat.I_bk(:)));

Opt_DS.Img_ImgFeat.numImg = Opt_DS.ImgPrc.numImg;
Opt_DS.Img_ImgFeat.T_init = true;
[R_ImgFeat,R_ImgFeat_Prc_time] = Crystalline_ImgAnlys_Feat(I,Opt_DS.Img_ImgFeat);
Opt_DS.Img_ImgFeat.T_init = false;

% % % % % % % % % % % % % % % % % % % % % 
% Object Parameterisation (binary images)
Opt_DS.Img_ObjPara.numImg = Opt_DS.ImgPrc.numImg;
Opt_DS.Img_ObjPara.T_init = true;
[R_ObjPara] = Crystalline_ImgAnlys_ObjPara(Opt_DS.ImgPrc.bk_bw,Opt_DS.ImgPrc.bk_bw,Opt_DS.ImgPrc.bk_img,'bk',Opt_DS.Img_ObjPara);
Opt_DS.Img_ObjPara.T_init = false;


% % % % % % % % % % % % % % % % % % % % % 
% Single particle data structure (length unclear > no indexing, preallocate empty cell structure and assemble table post-loop)
MDB_SiPa_header = { ...
    'Img_PrcIdx','area','solidity_ch','convex_area','convexity_area_frac','eccentricity','maxFeretDiameter', ...
    'PgFit_area','PgFit_perimeter', 'RectFit_area','RectFit_perimeter','RectFit_length','RectFit_width', ...
    'AspectRatio_RectFit','AspectRatio_BB',...
    'numCorner','bw_Frac_CryFit', 'CryFit_perimeter','CryFit_area',...
    };

MDB_SiPa_data = cell(Opt_DS.ImgPrc.numImg,1);


%% Parallel Loop Image Processing

% Progress bar for parfor
if isfield(Opt_DS.ImgPrc,'progressPrompt_On') && Opt_DS.ImgPrc.progressPrompt_On
    addpath(fullfile(Opt.path_FileExchange,'parfor_progress'))
    Opt_DS.ImgPrc.progressPrompt_string = sprintf('%s (DS %.0f/%.0f) - Image Processing',Opt_DS.ExpShorthand,Opt.k_DS,Opt.numDS);
    Opt_parfor_progress.path_fileCom = Opt.path_MatFolder_lc;
    parfor_progress_FrDo(Opt_DS.ExpShorthand,Opt_DS.ImgPrc.progressPrompt_string,Opt_DS.ImgPrc.numImg,Opt_parfor_progress);
end

for k = 1:Opt_DS.ImgPrc.numImg
    
    img_name_F = Img_filelist(k).name;
    img_name = img_name_F(1:end-4);
    
    Img_Name_list{k} = string(img_name);
    Img_PrcIdx_list(k) = k;
    
    I = imread(fullfile(Img_filelist(k).folder,img_name_F));
    % Crop scale-bar section
    I = I(Opt_DS.ImgPrc.ymin:Opt_DS.ImgPrc.ymax,Opt_DS.ImgPrc.xmin:Opt_DS.ImgPrc.xmax);  

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %% Image processing
    R = Crystalline_ImgPrc(I,Opt_ImgPrc);
    
    bw = R.bw;
    bw_ch = R.bw_ch;
    bw_thres = R.bw_thres;
    I_d = R.I_d;
    I_grad = R.I_grad;

    if Opt_DS.ImgPrc.BW_Print_On
        imwrite(bw,fullfile(Opt_DS.ImgPrc.path_ExportFolder_DS_bw,sprintf('%s_bw_k_%04.0f.bmp',...
                img_name,k)))
        if Opt_DS.ImgPrc.BW_Prc_Ext_Print_On
            max_I_d = 255;
            I_d_Img = uint8(I_d./max_I_d.*255);
            I_grad_Img = uint8(I_grad./max(I_grad).*255);

            imwrite(I,fullfile(Opt_DS.ImgPrc.path_ExportFolder_DS_bw_Ext,sprintf('%s_I_k_%04.0f.bmp',...
                img_name,k)))
            imwrite(I_d_Img,fullfile(Opt_DS.ImgPrc.path_ExportFolder_DS_bw_Ext,sprintf('%s_I_d_k_%04.0f.bmp',...
                img_name,k)))
            imwrite(I_grad_Img,fullfile(Opt_DS.ImgPrc.path_ExportFolder_DS_bw_Ext,sprintf('%s_I_grad_k_%04.0f.bmp',...
                img_name,k)))
            imwrite(bw_thres,fullfile(Opt_DS.ImgPrc.path_ExportFolder_DS_bw_Ext,sprintf('%s_bw_thres_k_%04.0f.bmp',...
                img_name,k)))
        end
    end
   
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% Extract Direct Image Features
    [R_ImgAnlys_Feat_data,R_stats] = Crystalline_ImgAnlys_Feat(I,Opt_DS.Img_ImgFeat);
    
    % Transfer to data table
    R_ImgFeat(k,:) = R_ImgAnlys_Feat_data;
    R_ImgFeat_Prc_time(k,:) = R_stats.Prc_time;
    
    FM_LAPV = R_ImgAnlys_Feat_data.FM_LAPV;
    FM_GLLV = R_ImgAnlys_Feat_data.FM_GLLV;
    Int_mean = R_ImgAnlys_Feat_data.Int_mean;
    Int_std = R_ImgAnlys_Feat_data.Int_std;
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% Object Parameterisation (binary images, 1bit)  
    [R,CrystAnlys_R] = Crystalline_ImgAnlys_ObjPara(bw,bw_ch,I,img_name,Opt_DS.Img_ObjPara);
    
    % Transfer to data table
    R_ObjPara(k,:) = R;

    numObj = R.numObj;
    d_eqSph_Mean = R.d_eqSph_Mean;
    d_eqSph_Std = R.d_eqSph_Std;
    d_RectFit_max_Mean = R.d_RectFit_max_Mean;
    d_RectFit_max_Std = R.d_RectFit_max_Std;
    shapeFact_Mean = R.shapeFact_Mean;
    shapeFact_Std = R.shapeFact_Std;
    convexFact_Mean = R.convexFact_Mean;
    convexFact_Std = R.convexFact_Std;
    areaCov_Mean = R.areaCov_Mean*100;
    PSD_10 = R.PSD_10;
    PSD_50 = R.PSD_50;
    PSD_90 = R.PSD_90;
    PSD_span = R.PSD_span;      
        
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %% Retrieve MDB_File/Process Data
    idx_file = find(abs(MDB_File.file_time_list-Img_time_list(k)) == min(abs(MDB_File.file_time_list-Img_time_list(k))),1,'first');
    Img_StirringSpeed_top_list(k) = MDB_File.StirringSpeed_top(idx_file);
    Img_StirringSpeed_bottom_list(k) = MDB_File.StirringSpeed_bottom(idx_file);
    Img_Transmissivity_list(k) = MDB_File.Transmissivity(idx_file);
    Img_Temperature_list(k) = MDB_File.Temperature(idx_file);
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% Show Processed Image
        if Opt_DS.ImgPrc.Fig_Print_On && isfield(CrystAnlys_R,'I')
            % Need to change parfor to normal for loop to run correctly
            minutes = floor(Img_time_list(k)/60);
            seconds = Img_time_list(k)-minutes*60;
            
            x_Dim = size(CrystAnlys_R.I,2);
            y_Dim = size(CrystAnlys_R.I,1);

            Fig_width = x_Dim*2;
            Fig_height = y_Dim;

            axes_Pos = [0 0 0.5 1];
            fig = figure('units','pixel','position',[100 100 Fig_width Fig_height]);

            axes1 = axes('Parent',fig,...
                'Position',axes_Pos);
            hold(axes1,'on');
            box(axes1,'on');
            axis([0,x_Dim,0,y_Dim])

            warning('off', 'images:initSize:adjustingMag')
            imshow(CrystAnlys_R.I);

            mask_1_RGB = cat(3,ones(size(bw)).*Opt_DS.mask_1_RGB_s(1), ...
                ones(size(bw)).*Opt_DS.mask_1_RGB_s(2), ones(size(bw)).*Opt_DS.mask_1_RGB_s(3));
            hold on
            h1 = imshow(mask_1_RGB);
            set(h1, 'AlphaData', bw.*0.2)

            [B,~] = bwboundaries(bw,'noholes');
            for i_B = 1:length(B)
                line('XData',B{i_B}(:,2),'YData',B{i_B}(:,1),'Color','green', ...
                    'LineWidth', 1)
            end
            rectangle('Position',[R_stats.x_ROI_sq R_stats.y_ROI_sq R_stats.dx_ROI_sq R_stats.dy_ROI_sq],'EdgeColor','r')

            annotation(fig,'textbox',...
                [0.52 0.90 0.765071428571429 0.0571428571428568],...
                'String',{sprintf([ ...
                'ID: %s\n' ...
                'Image: %s\n' ...
                'Iteration: %.0f/%.0f\n' ...
                'Rel. Time: %.0f:%02.0f\n' ...
                'Temperature: %.2f degC\n' ...
                'StirringSpeed_top: %.2f rpm\n' ...
                'StirringSpeed_bottom: %.2f rpm\n' ...
                'Transmissivity: %.2f\n' ...
                'FM_LAPV: %.2f\n' ...
                'FM_GLLV: %.2f\n' ...
                'Intensity: %.2f +/- %.2f\n\n' ...
                'numObj: %.0f\n' ...
                'Dia eqSph: %.2f +/- %.2f um\n' ...
                'Length RectFit: %.2f +/- %.2f um\n' ...
                'shapeFact: %.2f +/- %.2f\n' ...
                'convexFact: %.2f +/- %.2f\n' ...
                'AreaCov: %.2f %%\n' ...
                'PSD (Q0 10/50/90): %.2f um / %.2f um / %.2f um (span %.2f)' ...
                ], ...
                Opt_DS.ExpShorthand,img_name,k,Opt_DS.ImgPrc.numImg,...
                minutes,seconds,Img_Temperature_list(k),Img_StirringSpeed_top_list(k), ...
                Img_StirringSpeed_bottom_list(k),Img_Transmissivity_list(k), ...
                FM_LAPV,FM_GLLV,Int_mean,Int_std, ...
                numObj,d_eqSph_Mean,d_eqSph_Std,d_RectFit_max_Mean,d_RectFit_max_Std,...
                shapeFact_Mean,shapeFact_Std,convexFact_Mean,convexFact_Std, ...
                areaCov_Mean,PSD_10,PSD_50,PSD_90,PSD_span)},...
                'FitBoxToText','on','Fontsize',9,'Interpreter','none','BackgroundColor','w');
            
            print(fullfile(Opt_DS.ImgPrc.path_ExportFolder_DS_Img,sprintf('%s_k_%04.0f',img_name,k)),'-djpeg',Opt_DS.Graph.print_rS)
        

        end

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% Update DSDB_Data
    [~,numCol] = size(MDB_SiPa_header);
    numRow = CrystAnlys_R.numObj;
    if numRow < 1
        numRow = 1;
    end
    DSDB_tmp = nan(numRow,numCol);

    if CrystAnlys_R.numObj == 0
        DSDB_idx_Row = 1;
        DSDB_tmp(DSDB_idx_Row, strcmp('Img_PrcIdx',MDB_SiPa_header(1,:))) = k;

        DSDB_tmp(DSDB_idx_Row, strcmp('area',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('solidity_ch',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('convex_area',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('convexity_area_frac',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('eccentricity',MDB_SiPa_header(1,:))) = nan;
        
        DSDB_tmp(DSDB_idx_Row, strcmp('PgFit_area',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('PgFit_perimeter',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('RectFit_area',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('RectFit_perimeter',MDB_SiPa_header(1,:))) = nan;

        DSDB_tmp(DSDB_idx_Row, strcmp('RectFit_length',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('RectFit_width',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('AspectRatio_RectFit',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('AspectRatio_BB',MDB_SiPa_header(1,:))) = nan;

        DSDB_tmp(DSDB_idx_Row, strcmp('numCorner',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('bw_Frac_CryFit',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('CryFit_perimeter',MDB_SiPa_header(1,:))) = nan;
        DSDB_tmp(DSDB_idx_Row, strcmp('CryFit_area',MDB_SiPa_header(1,:))) = nan;
    end


    for i_Obj = 1:CrystAnlys_R.numObj
        DSDB_idx_Row = i_Obj;
        DSDB_tmp(DSDB_idx_Row, strcmp('Img_PrcIdx',MDB_SiPa_header(1,:))) = k;

        DSDB_tmp(DSDB_idx_Row, strcmp('area',MDB_SiPa_header(1,:))) = CrystAnlys_R.area(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('solidity_ch',MDB_SiPa_header(1,:))) = CrystAnlys_R.solidity_ch(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('convex_area',MDB_SiPa_header(1,:))) = CrystAnlys_R.convex_area(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('convexity_area_frac',MDB_SiPa_header(1,:))) = CrystAnlys_R.convexity_area_frac(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('eccentricity',MDB_SiPa_header(1,:))) = CrystAnlys_R.eccentricity(i_Obj);

        DSDB_tmp(DSDB_idx_Row, strcmp('PgFit_area',MDB_SiPa_header(1,:))) = CrystAnlys_R.PgFit_area(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('PgFit_perimeter',MDB_SiPa_header(1,:))) = CrystAnlys_R.PgFit_perimeter(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('RectFit_area',MDB_SiPa_header(1,:))) = CrystAnlys_R.RectFit_area(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('RectFit_perimeter',MDB_SiPa_header(1,:))) = CrystAnlys_R.RectFit_perimeter(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('RectFit_length',MDB_SiPa_header(1,:))) = CrystAnlys_R.RectFit_length(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('RectFit_width',MDB_SiPa_header(1,:))) = CrystAnlys_R.RectFit_width(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('AspectRatio_RectFit',MDB_SiPa_header(1,:))) = CrystAnlys_R.AspectRatio_RectFit(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('AspectRatio_BB',MDB_SiPa_header(1,:))) = CrystAnlys_R.AspectRatio_BB(i_Obj);

        DSDB_tmp(DSDB_idx_Row, strcmp('numCorner',MDB_SiPa_header(1,:))) = CrystAnlys_R.numCorner(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('bw_Frac_CryFit',MDB_SiPa_header(1,:))) = CrystAnlys_R.bw_Frac_CryFit(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('CryFit_perimeter',MDB_SiPa_header(1,:))) = CrystAnlys_R.CryFit_perimeter(i_Obj);
        DSDB_tmp(DSDB_idx_Row, strcmp('CryFit_area',MDB_SiPa_header(1,:))) = CrystAnlys_R.CryFit_area(i_Obj);
    end
    
    MDB_SiPa_data{k} = DSDB_tmp;

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% End Sequence
    close all
    if isfield(Opt_DS.ImgPrc,'progressPrompt_On') && Opt_DS.ImgPrc.progressPrompt_On
        parfor_progress_FrDo(Opt_DS.ExpShorthand,Opt_DS.ImgPrc.progressPrompt_string,-1,Opt_parfor_progress);
    else
        fprintf('%s (DS %.0f/%.0f) - Image Processing: %.0f/%.0f (Elapsed time: %.0f sec)\n',Opt_DS.ExpShorthand,Opt.k_DS,Opt.numDS,k,Opt_DS.ImgPrc.numImg,toc(Opt.tic))
    end
end

% Close waitbar
if isfield(Opt_DS.ImgPrc,'progressPrompt_On') && Opt_DS.ImgPrc.progressPrompt_On
    parfor_progress_FrDo(Opt_DS.ExpShorthand,Opt_DS.ImgPrc.progressPrompt_string,0,Opt_parfor_progress);
end

%% Save performance data
Opt_DS.ImgPrc.ImgFeat_Prc_time = R_ImgFeat_Prc_time;

%% Re-structure DSDB_data
[nrows_DSDB_data,~] = cellfun(@size,MDB_SiPa_data,'uni',false);
numObj = sum([nrows_DSDB_data{:}]);

MDB_SiPa = array2table(nan(numObj,size(MDB_SiPa_header,2)), ...
    'VariableNames', MDB_SiPa_header(1,:));

numObj_cum = 1;
for k =1:Opt_DS.ImgPrc.numImg
    numObj = nrows_DSDB_data{k};
    MDB_SiPa(numObj_cum:numObj_cum+numObj-1,:) = array2table(MDB_SiPa_data{k});

    numObj_cum = numObj_cum+numObj;
end

% Export Crystalline Particle Data
writetable(MDB_SiPa,fullfile(Opt_DS.path_ExportFolder_Exp,sprintf('%s_MDB_SiPa_%s.csv',Opt_DS.ExpShorthand,datestr(now, 'yyyymmdd'))))

%% Assemble Table Structure
Img_time_list_min = Img_time_list./60;

MDB_Img_VariableNames = { ...
    'Img_Name', ...
    'Img_PrcIdx', ...
    'Img_time_s', ...
    'Img_time_min', ...
    'Img_StirringSpeed_top', ...
    'Img_StirringSpeed_bottom', ...
    'Img_Transmissivity', ...
    'Img_Temperature', ...
    };

MDB_Img = array2table(nan(Opt_DS.ImgPrc.numImg,length(MDB_Img_VariableNames)),'VariableNames',MDB_Img_VariableNames);

MDB_Img.Img_Name = [Img_Name_list{:}]';
MDB_Img.Img_PrcIdx = Img_PrcIdx_list;
MDB_Img.Img_time_s = Img_time_list;
MDB_Img.Img_time_min = Img_time_list_min;
MDB_Img.Img_StirringSpeed_top = Img_StirringSpeed_top_list;
MDB_Img.Img_StirringSpeed_bottom = Img_StirringSpeed_bottom_list;
MDB_Img.Img_Transmissivity = Img_Transmissivity_list;
MDB_Img.Img_Temperature = Img_Temperature_list;

% Add R_ImgFeat, R_ObjPara
MDB_Img = [MDB_Img,R_ImgFeat,R_ObjPara];

%% End of Script
fprintf('%s - %s (%s) COMPLETE (Elapsed time: %.0f sec)\n',Opt.ProjectShorthand,Opt_DS.ExpShorthand,mfilename(),toc(Opt.tic))

