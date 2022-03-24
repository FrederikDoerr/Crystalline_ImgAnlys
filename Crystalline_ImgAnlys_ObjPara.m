function[R,CrystAnlys_R] = Crystalline_ImgAnlys_ObjPara(bw,bw_ch,I,img_name,Opt_ObjPara)
%Image object paramertisation.
%   [R,CrystAnlys_R] = Crystalline_ImgAnlys_ObjPara(bw,bw_ch,img_name,Opt_ObjPara) 
%   function to extract object features related to size and shape. bw is a
%   binary mask of all objects. bw_ch is a binary convex hull of all objects.
%   I is the original 8bit image. Returns R with extracted mean size/shape 
%   descriptors and CrystAnlys_R containing size/shape descriptors of all 
%   objects.
%   Initialise using Opt.T_init = true and Opt.numImg to create empty data
%   table.
%
%   Dependencies on FileExchange packages/functions:
%       > CrystAnlys
%       > Crystalline_PSD_Calc
%
%   OPTIONS:
%      'PSD_calc_numThres'  Minimum number of object for PSD calculations.
%      'T_init'             Initialise empty table creation (requires
%                           numImg)
%      'numImg'             Total number of images in dataset.
% 
%   additional OPTIONS related to CrystAnlys and Crystalline_PSD_Calc!
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2020b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/FrederikDoerr


R = table();


    stats_bw = regionprops(bw,'Area','Perimeter');
    bw_area = sum([stats_bw(:).Area]);  
    bw_perimeter = sum([stats_bw(:).Perimeter]);

    stats_bw_ch = regionprops(bw_ch,'Area','Perimeter');
    bw_ch_Area = sum([stats_bw_ch(:).Area]);  
    bw_ch_Perimeter = sum([stats_bw_ch(:).Perimeter]);
    
    bw_siPa_Anlys_Area = nan;
    bw_siPa_Anlys_Perimeter = nan;


    %% Single Particle Analysis        
    
    R.area = sum(bw(:));
    [~,numObj] = bwlabel(bw);
    
    R.numObj = numObj;
    R.x_Dim = size(bw,2);
    R.y_Dim = size(bw,1);
    
    R.bw_area = bw_area;
    R.bw_perimeter = bw_perimeter;
    R.bw_ch_Area = bw_ch_Area;
    R.bw_ch_Perimeter = bw_ch_Perimeter;
    R.bw_siPa_Anlys_Area = bw_siPa_Anlys_Area;
    R.bw_siPa_Anlys_Perimeter = bw_siPa_Anlys_Perimeter;
    
    CrystAnlys_Opts = struct();
    CrystAnlys_Opts.path_MatFolder = Opt_ObjPara.path_FileExchange;
    CrystAnlys_Opts.path_FileExchange = Opt_ObjPara.path_FileExchange;
    CrystAnlys_Opts.path_ExportFolder = Opt_ObjPara.path_ExportFolder;
    CrystAnlys_Opts.ID = img_name;
    CrystAnlys_Opts.I = I;
    CrystAnlys_Opts.Imagepixelsize = Opt_ObjPara.imagepixelsize;
    CrystAnlys_Opts.CryFit_On = false;
    CrystAnlys_Opts.PgFit_On = false;
    CrystAnlys_Opts.Print_Control_On = Opt_ObjPara.Fig_Print_On;
    CrystAnlys_Opts.Print_imwrite_On = false;
    CrystAnlys_Opts.Obj_Print_Control_On = Opt_ObjPara.BW_Obj_Print_On;
    
    addpath(Opt_ObjPara.path_FileExchange)
    [CrystAnlys_R] = CrystAnlys(bw,CrystAnlys_Opts);

    % Post calculation
    R.d_eqSph_Mean = mean(CrystAnlys_R.d_eqSph);
    R.d_eqSph_Std = std(CrystAnlys_R.d_eqSph);
    R.d_eqSqu_Mean = mean(CrystAnlys_R.d_eqSqu);
    R.d_eqSqu_Std = std(CrystAnlys_R.d_eqSqu);

    RectFit_length = CrystAnlys_R.RectFit_length;
    RectFit_width = CrystAnlys_R.RectFit_width;
    R.d_RectFit_max_Mean = mean(RectFit_length);
    R.d_RectFit_max_Std = std(RectFit_length);
    R.d_RectFit_min_Mean = mean(RectFit_width);
    R.d_RectFit_min_Std = std(RectFit_width);

    R.area_Mean = mean(CrystAnlys_R.area).*(Opt_ObjPara.imagepixelsize)^2;
    R.area_Std = std(CrystAnlys_R.area).*(Opt_ObjPara.imagepixelsize)^2;
    R.area_Min = min(CrystAnlys_R.area).*(Opt_ObjPara.imagepixelsize)^2;
    R.area_Max = max(CrystAnlys_R.area).*(Opt_ObjPara.imagepixelsize)^2;


    R.areaCov_Mean = mean(CrystAnlys_R.area)/(R.x_Dim*R.y_Dim);

    R.shapeFact_Mean = mean(CrystAnlys_R.AspectRatio_RectFit);
    R.shapeFact_Std = std(CrystAnlys_R.AspectRatio_RectFit);

    R.convexFact_Mean = mean(CrystAnlys_R.convexity_area_frac);
    R.convexFact_Std = std(CrystAnlys_R.convexity_area_frac);
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% Calculate PSD
    if numObj > Opt_ObjPara.PSD_calc_numThres
        addpath(Opt_ObjPara.path_FileExchange)
        Opt_PSD_Calc = Opt_ObjPara;
        Opt_PSD_Calc.Graph_name = sprintf('%s_PSD_Calc',img_name);
        [~,~,~,~,R_PSD_Calc] = Crystalline_PSD_Calc(CrystAnlys_R.d_eqSph,Opt_PSD_Calc);
        R.PSD_10 = R_PSD_Calc.PSD_0_d10;
        R.PSD_50 = R_PSD_Calc.PSD_0_d50;
        R.PSD_90 = R_PSD_Calc.PSD_0_d90;
        R.PSD_span = R_PSD_Calc.PSD_0_Span;
    else
        R.PSD_10 = nan;
        R.PSD_50 = nan;
        R.PSD_90 = nan;
        R.PSD_span = nan;
    end

    
% Create empty table     
if Opt_ObjPara.T_init
    R = array2table(nan(Opt_ObjPara.numImg,length(R.Properties.VariableNames)),'VariableNames',R.Properties.VariableNames);
end



