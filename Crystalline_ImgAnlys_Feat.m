function[R_data,R_stats] = Crystalline_ImgAnlys_Feat(I,Opt)
%Image feature extraction.
%   [R_data,R_stats] = Crystalline_ImgAnlys_Feat(I,Opt) function to extract
%   image features from a greyscale image I.  
%   Returns extracted image feature data (R_data) and function statistics
%   (R_stats) e.g. sub-funstion processing times.
%   Initialise using Opt.T_init = true and Opt.numImg to create empty data
%   table.
%
%   Dependencies on FileExchange packages/functions:
%       > fmeasure\fmeasure
%
%   OPTIONS:
%      'PLS_sub_On'         Extract sub-set of image features.
%      'T_init'             Initialise empty table creation (requires
%                           numImg)
%      'numImg'             Total number of images in dataset.
% 
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2020b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/FrederikDoerr

R_data = table();
Prc_time = table();

%% Removing image border regions (optical aberration)
I = I(Opt.y_ROI_sq:Opt.y_ROI_sq+Opt.dy_ROI_sq,Opt.x_ROI_sq:Opt.x_ROI_sq+Opt.dx_ROI_sq);
ROI = ones(Opt.dy_ROI_sq+1,Opt.dx_ROI_sq+1);

if isfield(Opt,'I_bk_pxVar')
    I = imsubtract(I,Opt.I_bk_pxVar );
end


addpath(fullfile(Opt.path_FileExchange,'fmeasure\fmeasure'))
if Opt.subSet_On
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % % % % %  Assess relative degree of focus of an image
    Opt.tic_f = tic; WSize = 15;
    R_data.FM_LAPV = fmeasure_FD(I, 'LAPV',ROI,WSize); % Variance of laplacian (Pech2000)
    Prc_time.FM_LAPV = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 15;
    R_data.FM_GLLV = fmeasure_FD(I, 'GLLV',ROI,WSize); % Graylevel local variance (Pech2000)
    Prc_time.FM_GLLV = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 15;
    R_data.FM_BREN = fmeasure_FD(I, 'BREN',ROI,WSize);
    Prc_time.FM_BREN = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 5;
    R_data.FM_HELM_005 = fmeasure_FD(I, 'HELM',ROI,WSize);
    Prc_time.FM_HELM_005 = toc(Opt.tic_f);
    

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % % % % %  Basic Statistics for Image Intensity
    Opt.tic_f = tic;
    R_data.Int_mean = mean2(I);
    R_data.Int_std = std2(I);
    % Variance
    R_data.Int_var = var(double(I(:)));

    R_data.Int_median = median(double(I(:)));
    
    % Index of dispersion
    % Distribution 	variance-to-mean ratio (VMR) 	
    % constant random variable 	VMR = 0 	not dispersed
    % binomial distribution 	0 < VMR < 1 	under-dispersed
    % Poisson distribution 	VMR = 1 	
    % negative binomial distribution 	VMR > 1 	over-dispersed 

    R_data.Int_IdxD = R_data.Int_var/R_data.Int_mean;
    Prc_time.IntStats_Global = toc(Opt.tic_f);
    

else
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % % % % %  Assess relative degree of focus of an image
    Opt.tic_f = tic; WSize = 15;
    R_data.FM_LAPV = fmeasure_FD(I, 'LAPV',ROI,WSize); % Variance of laplacian (Pech2000)
    Prc_time.FM_LAPV = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 23;
    R_data.FM_LAPV = fmeasure_FD(I, 'LAPV_WSize',ROI,WSize); % Variance of laplacian (Pech2000)
    Prc_time.FM_LAPV = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 5;
    R_data.FM_GLLV_005 = fmeasure_FD(I, 'GLLV',ROI,WSize); % Graylevel local variance (Pech2000)
    Prc_time.FM_GLLV_005 = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 15;
    R_data.FM_GLLV = fmeasure_FD(I, 'GLLV',ROI,WSize); % Graylevel local variance (Pech2000)
    Prc_time.FM_GLLV = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_BREN = fmeasure_FD(I, 'BREN',ROI,WSize);
    Prc_time.FM_BREN = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 15;
    R_data.FM_GDER = fmeasure_FD(I, 'GDER',ROI,WSize);
    Prc_time.FM_GDER = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_GLVA = fmeasure_FD(I, 'GLVA',ROI,WSize);
    Prc_time.FM_GLVA = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_GRAE = fmeasure_FD(I, 'GRAE',ROI,WSize);
    Prc_time.FM_GRAE = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_GRAT = fmeasure_FD(I, 'GRAT',ROI,WSize);
    Prc_time.FM_GRAT = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_GRAS = fmeasure_FD(I, 'GRAS',ROI,WSize);
    Prc_time.FM_GRAS = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 5;
    R_data.FM_HELM_005 = fmeasure_FD(I, 'HELM',ROI,WSize);
    Prc_time.FM_HELM_005 = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 15;
    R_data.FM_HELM = fmeasure_FD(I, 'HELM',ROI,WSize);
    Prc_time.FM_HELM = toc(Opt.tic_f);

    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_HISE = fmeasure_FD(I, 'HISE',ROI,WSize);
    Prc_time.FM_HISE = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_HISR = fmeasure_FD(I, 'HISR',ROI,WSize);
    Prc_time.FM_HISR = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 15;
    R_data.FM_LAPM = fmeasure_FD(I, 'LAPM',ROI,WSize);
    Prc_time.FM_LAPM = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 23;
    R_data.FM_LAPM_023 = fmeasure_FD(I, 'LAPM_WSize',ROI,WSize);
    Prc_time.FM_LAPM_023 = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_SFRQ = fmeasure_FD(I, 'SFRQ',ROI,WSize);
    Prc_time.FM_SFRQ = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 15;
    R_data.FM_TENG = fmeasure_FD(I, 'TENG',ROI,WSize);
    Prc_time.FM_TENG = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = 15;
    R_data.FM_TENV = fmeasure_FD(I, 'TENV',ROI,WSize);
    Prc_time.FM_TENV = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_VOLA = fmeasure_FD(I, 'VOLA',ROI,WSize);
    Prc_time.FM_VOLA = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_WAVS = fmeasure_FD(I, 'WAVS',ROI,WSize);
    Prc_time.FM_WAVS = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_WAVV = fmeasure_FD(I, 'WAVV',ROI,WSize);
    Prc_time.FM_WAVV = toc(Opt.tic_f);
    
    Opt.tic_f = tic; WSize = nan;
    R_data.FM_WAVR = fmeasure_FD(I, 'WAVR',ROI,WSize);
    Prc_time.FM_WAVR = toc(Opt.tic_f);
    
    Opt.tic_f = tic;
    Ix = I;
    Iy = I;
    I_d_2 = diff(I, 15, 2);
    I_d_1 = diff(I, 15, 1);
    Ix(:,1:end-1) = imresize(I_d_2,size(Ix(:,1:end-1)),'nearest');
    Iy(1:end-1,:) = imresize(I_d_1,size(Iy(1:end-1,:)),'nearest');
    R_data.FM_SFRQ_15 = mean2(sqrt(double(Iy.^2+Ix.^2)));
    Prc_time.FM_SFRQ_15 = toc(Opt.tic_f);
    
    Opt.tic_f = tic;
    Ix = I;
    Iy = I;
    I_d_2 = diff(I, 50, 2);
    I_d_1 = diff(I, 50, 1);
    Ix(:,1:end-1) = imresize(I_d_2,size(Ix(:,1:end-1)),'nearest');
    Iy(1:end-1,:) = imresize(I_d_1,size(Iy(1:end-1,:)),'nearest');
    R_data.FM_SFRQ_50 = mean2(sqrt(double(Iy.^2+Ix.^2)));
    Prc_time.FM_SFRQ_50 = toc(Opt.tic_f);

    Opt.tic_f = tic;
    SOB = fspecial('sobel');
    ISOB = imfilter(I, SOB, 'replicate', 'conv');
    R_data.FM_SOBL = std2(ISOB)^2;
    Prc_time.FM_SOBL = toc(Opt.tic_f);
    
    Opt.tic_f = tic;
    PRW = fspecial('prewitt');
    IPRW = imfilter(I, PRW, 'replicate', 'conv');
    R_data.FM_PRWI = std2(IPRW)^2;
    Prc_time.FM_PRWI = toc(Opt.tic_f);
    
    Opt.tic_f = tic;
    PRW = fspecial('prewitt');
    Gx = imfilter(double(I), PRW, 'replicate', 'conv');
    Gy = imfilter(double(I), PRW', 'replicate', 'conv');
    FM = Gx.^2 + Gy.^2;
    R_data.FM_PRWM = std2(FM)^2;
    R_data.FM_PRWV = mean2(FM);
    Prc_time.FM_PRWM = toc(Opt.tic_f);

    
    % % % % %  Local range of image
    Opt.tic_f = tic;
    WSize = 7;
    nhood = ones(WSize,WSize);
    J = rangefilt(I,nhood);
    R_data.RngFilt_007_ent = entropy(J);
    R_data.RngFilt_007_V = std2(J)^2;
    R_data.RngFilt_007_M = mean2(J);
    Prc_time.RngFilt_007 = toc(Opt.tic_f);
    
    Opt.tic_f = tic;
    WSize = 39;
    nhood = ones(WSize,WSize);
    J = rangefilt(I,nhood);
    R_data.RngFilt_039_ent = entropy(J);
    R_data.RngFilt_039_V = std2(J)^2;
    R_data.RngFilt_039_M = mean2(J);
    Prc_time.RngFilt_039 = toc(Opt.tic_f);
    
    Opt.tic_f = tic;
    WSize = 71;
    nhood = ones(WSize,WSize);
    J = rangefilt(I,nhood);
    R_data.RngFilt_071_ent = entropy(J);
    R_data.RngFilt_071_V = std2(J)^2;
    R_data.RngFilt_071_M = mean2(J);
    Prc_time.RngFilt_071 = toc(Opt.tic_f);
    
    Opt.tic_f = tic;
    k_sigma = 1;
    WSize = 7;
    J = imgaussfilt(I,k_sigma,'FilterSize',[WSize,WSize]);
    WSize = 71;
    nhood = ones(WSize,WSize);
    J = rangefilt(J,nhood);
    R_data.RngFilt_Gauss_071_ent = entropy(J);
    R_data.RngFilt_Gauss_071_V = std2(J)^2;
    R_data.RngFilt_Gauss_071_M = mean2(J);
    Prc_time.RngFilt_Gauss_071 = toc(Opt.tic_f);
    
    Opt.tic_f = tic;
    k_sigma = 1;
    WSize = 15;
    J = imgaussfilt(I,k_sigma,'FilterSize',[WSize,WSize]);
    WSize = 157;
    nhood = ones(WSize,WSize);
    J = rangefilt(J,nhood);
    R_data.RngFilt_Gauss_157_ent = entropy(J);
    R_data.RngFilt_Gauss_157_V = std2(J)^2;
    R_data.RngFilt_Gauss_157_M = mean2(J);
    Prc_time.RngFilt_Gauss_157 = toc(Opt.tic_f);
    
    
    % % % % %  wiener2
    Opt.tic_f = tic;
    WSize = 39;
    J = wiener2(I,[WSize,WSize]);
    Id = I-J;
    R_data.Filt_wiener2_ent = entropy(J);
    R_data.Filt_wiener2_V = std2(Id)^2;
    R_data.Filt_wiener2_M = mean2(Id);
    Prc_time.Filt_wiener2 = toc(Opt.tic_f);
    
    % % % % %  imgaussfilt
    Opt.tic_f = tic;
    k_sigma = 1;
    WSize = 7;
    J = imgaussfilt(I,k_sigma,'FilterSize',[WSize,WSize]);
    Id = I-J;
    R_data.Filt_gauss_ent = entropy(J);
    R_data.Filt_gaussV = std2(Id)^2;
    R_data.Filt_gaussM = mean2(Id);
    Prc_time.Filt_gauss = toc(Opt.tic_f);
    
    
    % % % % %  Local entropy of grayscale image
    Opt.tic_f = tic;
    WSize = 7;
    nhood = ones(WSize,WSize);
    J = entropyfilt(I,nhood);
    R_data.ImgStats_EntFilt_007_ent = entropy(J);
    R_data.ImgStats_EntFilt_007_V = std2(J)^2;
    R_data.ImgStats_EntFilt_007_M = mean2(J);
    Prc_time.ImgStats_EntFilt_007 = toc(Opt.tic_f);
    
    Opt.tic_f = tic;
    WSize = 71;
    nhood = ones(WSize,WSize);
    J = entropyfilt(I,nhood);
    R_data.ImgStats_EntFilt_071_ent = entropy(J);
    R_data.ImgStats_EntFilt_071_V = std2(J)^2;
    R_data.ImgStats_EntFilt_071_M = mean2(J);
    Prc_time.ImgStats_EntFilt_071 = toc(Opt.tic_f);
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % % % % %  Basic Statistics for Image Intensity
    Opt.tic_f = tic;
    R_data.Int_mean = mean2(I);
    R_data.Int_std = std2(I);
    % Variance
    R_data.Int_var = var(double(I(:)));

    R_data.Int_median = median(double(I(:)));
    
    % Index of dispersion
    % Distribution 	variance-to-mean ratio (VMR) 	
    % constant random variable 	VMR = 0 	not dispersed
    % binomial distribution 	0 < VMR < 1 	under-dispersed
    % Poisson distribution 	VMR = 1 	
    % negative binomial distribution 	VMR > 1 	over-dispersed 

    R_data.Int_IdxD = R_data.Int_var/R_data.Int_mean;
    Prc_time.IntStats_Global = toc(Opt.tic_f);
    
    

    % % % % %  Properties of gray-level co-occurrence matrix
    Opt.tic_f = tic;
    glcm = graycomatrix(I,'Offset',[5 0;0 5]);
    stats_glcm = graycoprops(glcm,'all');
    
    R_data.GLCM_CONTM = mean(stats_glcm.Contrast);
    R_data.GLCM_CONTV = std(stats_glcm.Contrast)^2;
    
    R_data.GLCM_CORRM = mean(stats_glcm.Correlation);
    R_data.GLCM_CORRV = std(stats_glcm.Correlation)^2;
    
    R_data.GLCM_ENGM = mean(stats_glcm.Energy);
    R_data.GLCM_ENGV = std(stats_glcm.Energy)^2;
        
    R_data.GLCM_HOMM = mean(stats_glcm.Homogeneity);
    R_data.GLCM_HOMV = std(stats_glcm.Homogeneity)^2;
    
    Prc_time.GLCM = toc(Opt.tic_f);

    
    % % % % %  ROI Histogram Parameterisation
    Opt.tic_f = tic;
    edges = 0:1:256;
    bins = 0:1:255;
    stats_hist = histogram(I(:),edges);  
    Q0_Int = cumtrapz(bins,stats_hist.Values);
    Q0_Int = Q0_Int./Q0_Int(end);

    y = bins;
    x = Q0_Int;
    [x, idx] = unique(x); 
    R_data.I_hist_10 = interp1(x, y(idx), 0.10);
    R_data.I_hist_25 = interp1(x, y(idx), 0.25);
    R_data.I_hist_50 = interp1(x, y(idx), 0.50);
    R_data.I_hist_75 = interp1(x, y(idx), 0.75);
    R_data.I_hist_90 = interp1(x, y(idx), 0.90);
    R_data.I_hist_Span = (R_data.I_hist_90-R_data.I_hist_10)/R_data.I_hist_50;

    % Discretization
    R_data.I_hist_Sum_Q4 = sum(stats_hist.Values(bins <= max(bins)*1.00 & bins > max(bins)*0.75));
    R_data.I_hist_Sum_Q3 = sum(stats_hist.Values(bins <= max(bins)*0.75 & bins > max(bins)*0.50));
    R_data.I_hist_Sum_Q2 = sum(stats_hist.Values(bins <= max(bins)*0.50 & bins > max(bins)*0.25));
    R_data.I_hist_Sum_Q1 = sum(stats_hist.Values(bins <= max(bins)*0.25 & bins > max(bins)*0.00));
    R_data.I_hist_Sum_bk = sum(stats_hist.Values(bins == 0));
    Prc_time.IntStats_hist = toc(Opt.tic_f);

    
    % % % % %  Otsu
    Opt.tic_f = tic;
    % T — Global threshold
    % EM — Effectiveness metric
    [R_data.ImgStats_Otsu_T,R_data.ImgStats_Otsu_EM] = graythresh(I(:));
    Prc_time.ImgStats_Otsu = toc(Opt.tic_f);
    
end

R_stats.Prc_time = Prc_time;

    
% Create empty table     
if Opt.T_init
    R_data = array2table(nan(Opt.numImg,length(R_data.Properties.VariableNames)),'VariableNames',R_data.Properties.VariableNames);
    
    R_stats = R_stats.Prc_time;
    R_stats = array2table(nan(Opt.numImg,length(R_stats.Properties.VariableNames)),'VariableNames',R_stats.Properties.VariableNames);
end

