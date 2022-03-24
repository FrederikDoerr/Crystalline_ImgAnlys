function[CrystAnlys_R] = CrystAnlys(bw,Opt)
%Crystal size and shape characterisation.
%   [CrystAnlys_R] = CrystAnlys(bw,Opt) function to extract object features
%   related to size and shape. bw is a binary mask of all objects.
%   Returns CrystAnlys_R with extracted mean size/shape descriptors
%   containing size/shape descriptors of all objects.
%
%   Dependencies on FileExchange packages/functions:
%       > ShapeFitting_MinBoundSuite
%
%   OPTIONS:
%      'Imagepixelsize'     Imagepixelsize to convert measured results from
%                          	px to um.
%      'PgFit_On'           Parallelogram shape fitting
%      'Print_Control_On'   Print results to figure;
%      'Print_imwrite_On'   Save figure with results
% 
%   additional OPTIONS related to CrystAnlys and Crystalline_PSD_Calc!
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2019b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/frederik-d


if ~isfield(Opt,'PgFit_On') || Opt.PgFit_On
    Opt.PgFit_On = true;
else
    Opt.PgFit_On = false;
end


if ~isfield(Opt,'Imagepixelsize')
    Opt.Imagepixelsize = 1;
    warning('%s: Imagepixelsize not defined (Imagepixelsize = 1)!\n',mfilename())
end

if Opt.Print_Control_On && isfield(Opt,'I')
    figure;
    imshow(Opt.I)
    hold on
elseif Opt.Print_Control_On
    figure;
    imshow(bw)
    hold on
end
if ~isfield(Opt,'Obj_Print_Control_On')
    Opt.Obj_Print_Control_On = false;
end

if ~isfield(Opt,'Print_imwrite_On')
    Opt.Print_imwrite_On = false;
end

if ~isfield(Opt,'Print_textwrite_On')
    Opt.Print_textwrite_On = false;
end

%% Basic regionprops statistics

[L,CrystAnlys_R.numObj] = bwlabeln(bw);
CrystAnlys_R.stats = regionprops(bw,'all');

CrystAnlys_R.area = [CrystAnlys_R.stats.Area].'.*Opt.Imagepixelsize^2;
CrystAnlys_R.d_eqSph = 2*sqrt(CrystAnlys_R.area/pi());
CrystAnlys_R.d_eqSqu = 2*sqrt(CrystAnlys_R.area);

CrystAnlys_R.solidity_ch = [CrystAnlys_R.stats.Solidity].';
CrystAnlys_R.convex_area = [CrystAnlys_R.stats.ConvexArea].'.*Opt.Imagepixelsize^2;
CrystAnlys_R.convexity_area_frac = [CrystAnlys_R.stats.Area].'./[CrystAnlys_R.stats.ConvexArea].';
CrystAnlys_R.eccentricity = [CrystAnlys_R.stats.Eccentricity].';

%% Loop - Preallocation

CrystAnlys_R.bw_cry_fit = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.CryFit_R = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.PgFit_area = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.PgFit_perimeter = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.RectFit_length = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.RectFit_width = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.RectFit_area = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.RectFit_perimeter = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.AspectRatio_RectFit = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.AspectRatio_BB = nan(CrystAnlys_R.numObj,1);

CrystAnlys_R.numCorner = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.CryFit_perimeter = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.CryFit_area = nan(CrystAnlys_R.numObj,1);
CrystAnlys_R.bw_Frac_CryFit = nan(CrystAnlys_R.numObj,1);

% [Myerson2002] p. 102, ka is the area shape factor that relates the area of a particle
% to the square of the characteristic dimension.

%% Loop - Start

if CrystAnlys_R.numObj > 0
    for k = 1:CrystAnlys_R.numObj

        bw_cry = L ==k;
        bw_cry = bwconvhull(bw_cry);
        [B,~] = bwboundaries(bw_cry,'noholes');
        
        if length(unique(B{1}(:,1))) < 2 || length(unique(B{1}(:,2))) < 2
            continue
        end
        

        % PgFit
        if Opt.PgFit_On
            addpath(fullfile(Opt.path_FileExchange,'ShapeFitting_MinBoundSuite'));
            [x_PgFit,y_PgFit,PgFit_area,PgFit_perimeter] = minboundparallelogram(B{1}(:,1),B{1}(:,2));
                x_PgFit = [x_PgFit;x_PgFit(2)];
                y_PgFit = [y_PgFit;y_PgFit(2)];
                CrystAnlys_R.PgFit_area(k)= PgFit_area.*Opt.Imagepixelsize^2;
                CrystAnlys_R.PgFit_perimeter(k) = PgFit_perimeter.*Opt.Imagepixelsize;
        else
            CrystAnlys_R.PgFit_area(k) = nan;
            CrystAnlys_R.PgFit_perimeter(k) = nan;
        end
        
        % RectFit
        addpath(fullfile(Opt.path_FileExchange,'ShapeFitting_MinBoundSuite'));
        [x_RectFit,y_RectFit,RectFit_area,RectFit_perimeter] = minboundrect(B{1}(:,1),B{1}(:,2));
        w = pdist([x_RectFit(1),y_RectFit(1);x_RectFit(2),y_RectFit(2)],'euclidean');
        l = pdist([x_RectFit(3),y_RectFit(3);x_RectFit(2),y_RectFit(2)],'euclidean');
        CrystAnlys_R.AspectRatio_RectFit(k) = max([w,l])/min([w,l]);
        
        CrystAnlys_R.RectFit_area(k) = RectFit_area.*Opt.Imagepixelsize^2;
        CrystAnlys_R.RectFit_perimeter(k) = RectFit_perimeter.*Opt.Imagepixelsize;
        CrystAnlys_R.RectFit_length(k) = l.*Opt.Imagepixelsize;
        CrystAnlys_R.RectFit_width(k) = w.*Opt.Imagepixelsize;
        
        %% Visualisation
        if Opt.Print_Control_On
            for i = 1:length(x_RectFit)-1
                line('XData',y_RectFit,'YData',x_RectFit,'Color','g', ...
                    'LineWidth', 1)
            end

            if Opt.Print_textwrite_On
                text(CrystAnlys_R.stats(k).Centroid(1)-20, CrystAnlys_R.stats(k).Centroid(2)+12, sprintf('Convexity: %.2f',CrystAnlys_R.stats(k).Solidity), 'Color', 'r' )
                text(CrystAnlys_R.stats(k).Centroid(1)-20, CrystAnlys_R.stats(k).Centroid(2)+24, sprintf('AR_RectF: %.2f',CrystAnlys_R.AspectRatio_RectFit(k)), 'Color', 'r' )
                text(CrystAnlys_R.stats(k).Centroid(1)-20, CrystAnlys_R.stats(k).Centroid(2)+36, sprintf('nC: %.0f',CrystAnlys_R.numCorner(k)), 'Color', 'r' )
                text(CrystAnlys_R.stats(k).Centroid(1)-20, CrystAnlys_R.stats(k).Centroid(2)+48, sprintf('nCF: %.2f',CrystAnlys_R.bw_Frac_CryFit(k)), 'Color', 'r' )
            end
            text(CrystAnlys_R.stats(k).Centroid(1)+10, CrystAnlys_R.stats(k).Centroid(2), sprintf('ID: %.0f',k), 'Color','r','FontSize',6)

            addpath(fullfile(Opt.path_FileExchange,'export_fig-5b3965b'))
            [imageData, ~] = export_fig(gcf);
            if Opt.Print_imwrite_On
                imwrite(imageData,fullfile(Opt.path_ExportFolder,sprintf('%s_%s.jpg',...
                    Opt.ID,mfilename())))
            else
                CrystAnlys_R.I = imageData;
            end
        else
            
        end
        %% General CrystAnlys_R.stats
        BB = CrystAnlys_R.stats(k).BoundingBox;
        CrystAnlys_R.AspectRatio_BB(k) = max([BB(3),BB(4)])/min([BB(3),BB(4)]);

    end
else
    if isfield(Opt,'I')
        CrystAnlys_R.I = Opt.I;
    end
            CrystAnlys_R.numObj = 0;
            CrystAnlys_R.bw = bw;
            
            CrystAnlys_R.area = nan;
            CrystAnlys_R.solidity_ch = nan;
            CrystAnlys_R.convex_area = nan;
            CrystAnlys_R.convexity_area_frac = nan;
            CrystAnlys_R.eccentricity = nan;
            CrystAnlys_R.maxFeretDiameter = nan;

            CrystAnlys_R.PgFit_area = nan;
            CrystAnlys_R.PgFit_perimeter = nan;
            CrystAnlys_R.RectFit_length = nan;
            CrystAnlys_R.RectFit_width = nan;
            CrystAnlys_R.RectFit_area = nan;
            CrystAnlys_R.RectFit_perimeter = nan;
            CrystAnlys_R.AspectRatio_RectFit = nan;
            CrystAnlys_R.AspectRatio_BB = nan;
            CrystAnlys_R.numCorner = nan;
            CrystAnlys_R.bw_Frac_CryFit = nan;
            CrystAnlys_R.CryFit_perimeter = nan;
            CrystAnlys_R.CryFit_area = nan;
            CrystAnlys_R.CryFit_length = nan;
            CrystAnlys_R.CryFit_angles = nan;
end

