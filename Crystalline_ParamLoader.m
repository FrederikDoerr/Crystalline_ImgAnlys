%------------------------------------------------------------------------------------------------
% Application: Crystalline Image and Data Analysis Framework
% Subroutine: Crystalline_ParamLoader
% 
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2020b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/FrederikDoerr


function[Opt_DS] = Crystalline_ParamLoader(Opt_DS,Opt)

%% Graph Opt_DSions

% Figure print resolution
Opt_DS.Graph.print_rS = '-r300';

% Dimensions
Opt_DS.Graph.width = 600;
Opt_DS.Graph.height = 500;
Opt_DS.Graph.height_ext = (Opt_DS.Graph.height*1.1);
Opt_DS.Graph.width_ext = 800;

Opt_DS.Graph.ax_1 = 0.10;
Opt_DS.Graph.ax_1_width = 0.80;
Opt_DS.Graph.ax_2 = 0.125;
Opt_DS.Graph.ax_2_height = 0.75;

Opt_DS.Graph.Color_Seq_blue = cell(1,9);
Opt_DS.Graph.Color_Seq_blue{1} = [236,231,242]./255; %
Opt_DS.Graph.Color_Seq_blue{2} = [208,209,230]./255; %
Opt_DS.Graph.Color_Seq_blue{3} = [166,189,219]./255; %
Opt_DS.Graph.Color_Seq_blue{4} = [116,169,207]./255; %
Opt_DS.Graph.Color_Seq_blue{5} = [54,144,192]./255; %
Opt_DS.Graph.Color_Seq_blue{6} = [5,112,176]./255; %
Opt_DS.Graph.Color_Seq_blue{7} = [4,90,141]./255; %
Opt_DS.Graph.Color_Seq_blue{8} = [2,56,88]./255; %
Opt_DS.Graph.Color_Seq_blue{9} = [0,24,35]./255; %

Opt_DS.Graph.Color_Seq_red = cell(1,9);
Opt_DS.Graph.Color_Seq_red{1} = [255,247,236]./255; %
Opt_DS.Graph.Color_Seq_red{2} = [254,232,200]./255; %
Opt_DS.Graph.Color_Seq_red{3} = [253,212,158]./255; %
Opt_DS.Graph.Color_Seq_red{4} = [253,187,132]./255; %
Opt_DS.Graph.Color_Seq_red{5} = [252,141,89]./255; %
Opt_DS.Graph.Color_Seq_red{6} = [239,101,72]./255; %
Opt_DS.Graph.Color_Seq_red{7} = [215,48,31]./255; %
Opt_DS.Graph.Color_Seq_red{8} = [179,0,0]./255; %
Opt_DS.Graph.Color_Seq_red{9} = [127,0,0]./255; %

Opt_DS.Graph.Color_Seq_green = cell(1,9);
Opt_DS.Graph.Color_Seq_green{1} = [247,252,253]./255; %
Opt_DS.Graph.Color_Seq_green{2} = [229,245,249]./255; %
Opt_DS.Graph.Color_Seq_green{3} = [204,236,230]./255; %
Opt_DS.Graph.Color_Seq_green{4} = [153,216,201]./255; %
Opt_DS.Graph.Color_Seq_green{5} = [102,194,164]./255; %
Opt_DS.Graph.Color_Seq_green{6} = [65,174,118]./255; %
Opt_DS.Graph.Color_Seq_green{7} = [35,139,69]./255; %
Opt_DS.Graph.Color_Seq_green{8} = [0,109,44]./255; %
Opt_DS.Graph.Color_Seq_green{9} = [0,68,27]./255; %

Opt_DS.Graph.Color_Seq_grey = cell(1,9);
Opt_DS.Graph.Color_Seq_grey{1} = [255,255,255]./255; %
Opt_DS.Graph.Color_Seq_grey{2} = [240,240,240]./255; %
Opt_DS.Graph.Color_Seq_grey{3} = [217,217,217]./255; %
Opt_DS.Graph.Color_Seq_grey{4} = [189,189,189]./255; %
Opt_DS.Graph.Color_Seq_grey{5} = [150,150,150]./255; %
Opt_DS.Graph.Color_Seq_grey{6} = [115,115,115]./255; %
Opt_DS.Graph.Color_Seq_grey{7} = [82,82,82]./255; %
Opt_DS.Graph.Color_Seq_grey{8} = [37,37,37]./255; %
Opt_DS.Graph.Color_Seq_grey{9} = [0,0,0]./255; %

%% Processing Parameter

% Progress updates parfor
Opt_DS.ImgPrc.progressPrompt_On = true;

Opt_DS.mask_1_RGB_s = [1,0,0];
Opt_DS.ImgPrc.Img_Format = Opt.MainPrc_T_Crystalline.Img_Format{Opt.k_DS};

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% % % Image Feature Extraction Settings
Opt_DS.Img_ImgFeat.path_FileExchange = Opt.path_FileExchange;
Opt_DS.Img_ImgFeat.subSet_On = true; % only extract defined subset of simple image features to accelerate routine analysis


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% % % Image Processing Parameters
k_ImgPrc = find(strcmp(Opt.ImgPrc_T_Crystalline.Prc_ImgAnlys_ID,Opt_DS.Prc_ImgAnlys_REF),1,'first');

% Image Pixelsize
Opt_DS.imagepixelsize = Opt.ImgPrc_T_Crystalline.ImagePixelSize_um_px(k_ImgPrc); % um/px

% Crop scale bar / timestamp
Opt_DS.ImgPrc.xmin = Opt.ImgPrc_T_Crystalline.Crop_xmin_px(k_ImgPrc);
Opt_DS.ImgPrc.xmax = Opt.ImgPrc_T_Crystalline.Crop_xmax_px(k_ImgPrc);
Opt_DS.ImgPrc.ymin = Opt.ImgPrc_T_Crystalline.Crop_ymin_px(k_ImgPrc);
Opt_DS.ImgPrc.ymax = Opt.ImgPrc_T_Crystalline.Crop_ymax_px(k_ImgPrc);
% Opt_DS.ImgPrc.Img_Cr_dynm = false;

Opt_DS.ImgPrc.I_obj_m_thres = Opt.ImgPrc_T_Crystalline.I_obj_m_thres(k_ImgPrc);
Opt_DS.ImgPrc.I_ob_thres = Opt.ImgPrc_T_Crystalline.I_ob_thres(k_ImgPrc);
Opt_DS.ImgPrc.bw_threSize = Opt.ImgPrc_T_Crystalline.bw_threSize_px(k_ImgPrc);
Opt_DS.ImgPrc.bw_strelSize = Opt.ImgPrc_T_Crystalline.bw_strelSize_px(k_ImgPrc);
Opt_DS.ImgPrc.Method_ConvexHullRefine = Opt.ImgPrc_T_Crystalline.Method_ConvexHullRefine(k_ImgPrc);
Opt_DS.ImgPrc.Info = Opt.ImgPrc_T_Crystalline.Info{k_ImgPrc};

Opt_DS.ImgPrc.Img_Int_median_thres = Opt.ImgPrc_T_Crystalline.Img_Int_median_thres(k_ImgPrc); % intensity theshold to identify image polarity

Opt_DS.ImgPrc.BW_Print_On = false;
Opt_DS.ImgPrc.Fig_Print_On = false;
Opt_DS.ImgPrc.BW_Prc_Ext_Print_On = false;
Opt_DS.ImgPrc.BW_Obj_Print_On = false;

Opt_DS.ImgPrc.path_ExportFolder_DS_Img = fullfile(Opt_DS.Dest_Data_Path,sprintf('%s_Img_MatPrc',Opt_DS.ExpShorthand));
Opt_DS.ImgPrc.path_ExportFolder_DS_bw = fullfile(Opt_DS.Dest_Data_Path,sprintf('%s_Img_MatPrc_bw',Opt_DS.ExpShorthand));
Opt_DS.ImgPrc.path_ExportFolder_DS_bw_Ext = fullfile(Opt_DS.Dest_Data_Path,sprintf('%s_Img_MatPrc_Prc_Ext',Opt_DS.ExpShorthand));

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% % % Data Processing
Opt_DS.DataPrc.Prc_Tthres = 1; % absolute temperature threshold for heating/cooling detection
Opt_DS.DataPrc.l_movAvg = 2;  % Moving average
Opt_DS.tstamp_offset = 0;

Opt_DS.DataPrc.av_StepSize = 21;
Opt_DS.DataPrc.Temp_const_thres = 2; %2K thrshold
Opt_DS.DataPrc.k_movsum_Temp_diff = 30;

Opt_DS.DataPrc.PLS_logic_nStages_min = 3; % Assuming at least three consequtive heating stages to idicate PLSR Data

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% % % Nucleation Detection Filter
Opt_DS.NucDet.k_hampel = 29;
Opt_DS.NucDet.k_sigma = 0.1;
Opt_DS.NucDet.k_sgolayfilt = 11;

Opt_DS.NucDet.k_hampel_Transmissivity = 7;
Opt_DS.NucDet.k_sigma_Transmissivity = 3;
Opt_DS.NucDet.k_sgolayfilt_Transmissivity = 9;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% % % PSD Processing Settings
Opt_DS.PSD_Anlys_On = false;
Opt_DS.Img_ObjPara.l_movAvg = 2;  % Moving average PSD calculation 
Opt_DS.Img_ObjPara.path_ProjectFolder = Opt.path_ProjectFolder;
Opt_DS.Img_ObjPara.path_FileExchange = Opt.path_FileExchange;
Opt_DS.Img_ObjPara.path_ExportFolder = Opt_DS.ImgPrc.path_ExportFolder_DS_Img;
Opt_DS.Img_ObjPara.Fig_Print_On = Opt_DS.ImgPrc.Fig_Print_On;
Opt_DS.Img_ObjPara.BW_Obj_Print_On = Opt_DS.ImgPrc.BW_Obj_Print_On;
Opt_DS.Img_ObjPara.imagepixelsize = Opt_DS.imagepixelsize;
Opt_DS.Img_ObjPara.PSD_calc_scale = 'logarithmic';
Opt_DS.Img_ObjPara.PSD_calc_numThres = 0;
Opt_DS.Img_ObjPara.PSD_calc_numBin = 100;
Opt_DS.Img_ObjPara.PSD_calc_minD = 0; % um
Opt_DS.Img_ObjPara.PSD_calc_maxD = 300; % um
Opt_DS.Img_ObjPara.print_rS = Opt_DS.Graph.print_rS;




