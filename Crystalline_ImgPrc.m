function R = Crystalline_ImgPrc(I,Opt_ImgPrc)
%Image processing.
%   [R] = Crystalline_ImgPrc(I,Opt_ImgPrc) 
%   function to process images. I is the original 8bit image. 
%   Returns R with bw is a binary mask of all objects. bw_ch is a binary
%   convex hull of each object.
%   
%   OPTIONS:
%      'I_obj_m_thres'      Threshold marker regions from gradient map
%      'I_ob_thres'         Threshold mask regions from gradient map
%      'bw_threSize'        Threshold size for despeckle small binary noise
%      'bw_strelSize'       Size of structuring element for  despeckling
% 
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2020b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/FrederikDoerr


% Background subtraction
I_d = imabsdiff(I,Opt_ImgPrc.bk_img);

% Multiply with Image ROI (bk_bw, logic mask)
I_d = immultiply(I_d,Opt_ImgPrc.bk_bw);

% Local low-contrast smoothing (edge aware)
edgeThreshold = 0.05;
amount = -0.75;
I_d = localcontrast(I_d,edgeThreshold,amount);
[I_d,noise_out] = wiener2(I_d,[11 11]);

% Gradient method for edge detection
I_grad = imgradient(I_d);

% Thresholding (marker and mask regions)
bw_obj_m = I_grad > Opt_ImgPrc.I_obj_m_thres;
bw_obj = I_grad > Opt_ImgPrc.I_ob_thres;
bw_thres = imbinarize(I_d);
bw = imreconstruct(bw_obj_m,bw_obj);

% Remove object connected to image border
bw = imclearborder(bw);

% Fill closed object contours
bw = imfill(bw,'holes');

% Convex hull
bw_ch = bwconvhull(bw,'object');

%% Loop over all Objects to refine borders

[L,numObj] = bwlabel(bw);

for k_Obj = 1:numObj
    bw_obj = L == k_Obj;
    bw_obj_ch = bwconvhull(bw_obj,'object');
    bw_obj_ch(bw_obj) = 0;
    
    % Single object median intensity and standard variation
    bw_obj_Int = median(I_d(bw_obj));
    bw_obj_Int_Err = std2(I_d(bw_obj));
    
    % Assess intensity of concave surface areas
    [L_ch,numObj_ch] = bwlabel(bw_obj_ch);
    for k_ch = 1:numObj_ch
        bw_obj_ch_k = L_ch == k_ch;
        bw_obj_ch_Int = median(I_d(bw_obj_ch_k)); 
        
        % Add concave surface areas to object if median intensity smaller
        % than single object median intensity minus standard variation
        if bw_obj_Int-bw_obj_Int_Err < bw_obj_ch_Int
            bw_obj(bw_obj_ch_k) = 1;
        end
    end
    bw(bw_obj) = 1;
end   

% Remove outer edge (partial volume effect)
se = strel('disk',1);
bw = imerode(bw,se);

% Despeckle small binary noise
bw = bwareaopen(bw,Opt_ImgPrc.bw_threSize);

% Despeckle fine binary strcutures (noise) and separate light touching
% structures
if Opt_ImgPrc.bw_strelSize > 0
    se = strel('disk',Opt_ImgPrc.bw_strelSize);
    bw = imopen(bw,se);
end

%% Update return structure
R.bw = bw;
R.bw_ch = bw_ch;
R.bw_thres = bw_thres;
R.I_d = I_d;
R.I_grad = I_grad;


