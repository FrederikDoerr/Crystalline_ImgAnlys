function[PSD_3_d25,PSD_3_d50,PSD_3_d75,PSD_3_Span,R] = Crystalline_PSD_Calc(d_eqSph,Opt_PSD_Calc)        
%Particle size distribution calculations.
%   [PSD_3_d25,PSD_3_d50,PSD_3_d75,PSD_3_Span,R] = Crystalline_PSD_Calc(d_eqSph,Opt_PSD_Calc) 
%   function to transform measured object sizes to particle size
%   distribution. Returns selected 
%
%   OPTIONS:
%      'PSD_calc_scale'   	PSD scale: 'logarithmic' or 'linear'.
%      'PSD_calc_minD'      PSD scale: minimum size.
%      'PSD_calc_maxD'      PSD scale: maximum size.
%      'PSD_calc_numBin'    PSD scale: number of bins.
%      'Fig_Print_On'       Export PSD figure (default: false).
%      'path_ExportFolder'  Folder path for figure export.
%      'Warning_Print_On'   Print out/suppress function warnings (default:
%                           true)
% 
%
%   Author
%   --------
%   Frederik Doerr, Aug 2020 (MATLAB R2020b)
%   frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
%   Github: https://github.com/FrederikDoerr


R = struct();

if ~isfield(Opt_PSD_Calc,'Fig_Print_On')
    Opt_PSD_Calc.Fig_Print_On = false;
end

if ~isfield(Opt_PSD_Calc,'Warning_Print_On')
    Opt_PSD_Calc.Warning_Print_On = true;
end


if (~isempty(d_eqSph) || ~any(isnan(d_eqSph))) && length(d_eqSph) > 2
    switch  lower(Opt_PSD_Calc.PSD_calc_scale)
        case 'logarithmic'
            if Opt_PSD_Calc.PSD_calc_minD == 0
                Opt_PSD_Calc.PSD_calc_minD = 1;
            end
            PSD_calc_minD = log10(Opt_PSD_Calc.PSD_calc_minD);
            PSD_calc_maxD = log10(Opt_PSD_Calc.PSD_calc_maxD);
            edges = logspace(PSD_calc_minD,PSD_calc_maxD,Opt_PSD_Calc.PSD_calc_numBin+1);
        case 'linear'
            edges = linspace(Opt_PSD_Calc.PSD_calc_minD,Opt_PSD_Calc.PSD_calc_maxD,Opt_PSD_Calc.PSD_calc_numBin+1);
    end
    
    numObj_PSD = length(d_eqSph);
    
    if any(d_eqSph <= edges(1))
        R.Warning_SizeMin = true;
        R.Warning_SizeMin_num = sum(d_eqSph <= edges(1));
        
        d_eqSph(d_eqSph < edges(1)) = edges(1)+(edges(2)-edges(1))*0.05;
        if Opt_PSD_Calc.Warning_Print_On
            warning('PSD with sizes smaller than first edge!')
        end

    else
        R.Warning_SizeMin = false;
    end
    
    if any(d_eqSph > edges(end))
        R.Warning_SizeMax = true;
        R.Warning_SizeMax_num = sum(d_eqSph > edges(end));
        
        d_eqSph(d_eqSph > edges(end)) = edges(end)-(edges(end)-edges(end-1))*0.05;
        if Opt_PSD_Calc.Warning_Print_On
            warning('PSD with sizes larger than last edge!')
        end
    else
    	R.Warning_SizeMax = false;
    end
    
    cumVol = sum((d_eqSph/2).^3*4/3*pi());
    
    centerbin = nan(Opt_PSD_Calc.PSD_calc_numBin,1);
    q0_List_bin = nan(Opt_PSD_Calc.PSD_calc_numBin,1);
    Q0_List_bin = nan(Opt_PSD_Calc.PSD_calc_numBin,1);
    q3_List_bin = nan(Opt_PSD_Calc.PSD_calc_numBin,1);
    Q3_List_bin = nan(Opt_PSD_Calc.PSD_calc_numBin,1);
    for i = 1:Opt_PSD_Calc.PSD_calc_numBin
        centerbin(i) = (edges(i) + edges(i+1))/2;
        logicArray = logical((d_eqSph > edges(i)).*(d_eqSph <= edges(i+1)));
        q0_List_bin(i) = sum(logicArray)/numObj_PSD;
        Q0_List_bin(i) = sum(q0_List_bin(1:i));
        
        q3_List_bin(i) = sum((d_eqSph(logicArray)./2).^3*4/3*pi())/cumVol;
        Q3_List_bin(i) = sum(q3_List_bin(1:i));
    end

    % Calculate D1,10, D1,50, D1,90, SPAN1
    x = Q0_List_bin;
    y = centerbin; 
    [x, idx] = unique(x); 
    R.PSD_0_d10 = interp1(x, y(idx), 0.10);
    R.PSD_0_d25 = interp1(x, y(idx), 0.25);
    R.PSD_0_d50 = interp1(x, y(idx), 0.50);
    R.PSD_0_d75 = interp1(x, y(idx), 0.75);
    R.PSD_0_d90 = interp1(x, y(idx), 0.90);
    R.PSD_0_Span = (R.PSD_0_d90-R.PSD_0_d10)/R.PSD_0_d50;
    
    x = Q3_List_bin;
    y = centerbin; 
    [x, idx] = unique(x); 
    R.PSD_3_d10 = interp1(x, y(idx), 0.10);
    R.PSD_3_d25 = interp1(x, y(idx), 0.25);
    R.PSD_3_d50 = interp1(x, y(idx), 0.50);
    R.PSD_3_d75 = interp1(x, y(idx), 0.75);
    R.PSD_3_d90 = interp1(x, y(idx), 0.90);
    R.PSD_3_Span = (R.PSD_3_d90-R.PSD_3_d10)/R.PSD_3_d50;
    
    % Gaussian fit
    counts = q3_List_bin;
    d = centerbin;
    % ft = fittype( 'gauss1' );
    %  probability density function
    ft = fittype( '1/(c1 * sqrt(2*pi()))*exp(-0.5*((x-b1)/c1)^2)' );
    % General model Gauss1:
    %      f(x) =  a1*exp(-((x-b1)/c1)^2)

    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 0];
    opts.Upper = [max(d) Inf];
    opts.StartPoint = [R.PSD_3_d50 R.PSD_3_Span];

    rng default;  % For reproducibility
    [fitresult, gof] = fit(d, counts, ft, opts );
    
    d_fit = 0:0.01:max(d);
    q3_GFit = 1/(fitresult.c1 * sqrt(2*pi()))*exp(-0.5*((d_fit-fitresult.b1)/fitresult.c1).^2);
    Q3_GFit = cumtrapz(d_fit,q3_GFit);
    
    % Graph
    if Opt_PSD_Calc.Fig_Print_On
        fig = figure;
        left_color = [0 0 0];
        right_color = [0 0 0];
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        b = bar(log10(d),counts,'DisplayName','Exp');
        b.FaceColor = [0,0,1];
        b.EdgeColor = [0,0,1];
        hold on
        yyaxis left
        plot(log10(centerbin),q3_List_bin,'k-','LineWidth',1.5,'DisplayName','q3 (Img)')
        yyaxis right
        plot(log10(centerbin),Q3_List_bin,'k--','LineWidth',1.5,'DisplayName','Q3 (Img)')
        yyaxis left
        plot(log10(d_fit),q3_GFit,'r-','LineWidth',1.5,'DisplayName','q3_GFit (Img)')
        yyaxis right
        plot(log10(d_fit),Q3_GFit,'r--','LineWidth',1.5,'DisplayName','Q3_GFit (Img)')
        lgd = legend();
        lgd.Interpreter = 'none';
        lgd.Location = 'NorthEast';
        xlabel('Size [µm]')
        yyaxis left
        ylabel('q3')
        yyaxis right
        ylabel('Q3')
       
%         set(gca,'XScale','log','XLim',[0.1 max(d_fit)])
         xlim([log10(1) max(log10(d_fit))]) 
        set(gca,'Xtick',[-1:3]); %// adjust manually; values in log scale
        set(gca,'Xticklabel',10.^get(gca,'Xtick')); %// use labels with linear values
        
        print(fullfile(Opt_PSD_Calc.path_ExportFolder,sprintf('%s',Opt_PSD_Calc.Graph_name)),'-djpeg',Opt_PSD_Calc.print_rS)
    end

    x = Q3_GFit;
    y = d_fit; 
    [x, idx] = unique(x); 
    R.PSD_3_d10_GFit = interp1(x, y(idx), 0.10);
    R.PSD_3_d25_GFit = interp1(x, y(idx), 0.25);
    R.PSD_3_d50_GFit = interp1(x, y(idx), 0.50);
    R.PSD_3_d75_GFit = interp1(x, y(idx), 0.75);
    R.PSD_3_d90_GFit = interp1(x, y(idx), 0.90);
    R.Span_3_GFit = (R.PSD_3_d90_GFit-R.PSD_3_d10_GFit)/R.PSD_3_d50_GFit;
    
    
    R.centerbin = centerbin;
    R.q0 = q0_List_bin;
    R.Q0 = Q0_List_bin;
    R.q3 = q3_List_bin;
    R.Q3 = Q3_List_bin;
    R.q3_GFit = q3_GFit;
    R.Q3_GFit = Q3_GFit;

    R.fitresult = fitresult;
    R.gof = gof;
    
    
    % Volume Mean Diameter (D_43)
    R.D_43 = mean(d_eqSph);
    R.D_43_StdDev = std(d_eqSph);

else    
    
    R.PSD_0_d10 = nan;
    R.PSD_0_d25 = nan;
    R.PSD_0_d50 = nan;
    R.PSD_0_d75 = nan;
    R.PSD_0_d90 = nan;
    R.PSD_0_Span = nan;

    R.PSD_3_d10 = nan;
    R.PSD_3_d25 = nan;
    R.PSD_3_d50 = nan;
    R.PSD_3_d75 = nan;
    R.PSD_3_d90 = nan;
    R.PSD_3_Span = nan;
    
    R.PSD_3_d10_GFit = nan;
    R.PSD_3_d25_GFit = nan;
    R.PSD_3_d50_GFit = nan;
    R.PSD_3_d75_GFit = nan;
    R.PSD_3_d90_GFit = nan;
    R.Span_3_GFit = nan;
    
    
    R.centerbin = nan;
    R.q0 = nan;
    R.Q0 = nan;
    R.q3 = nan;
    R.Q3 = nan;
    R.q3_GFit = nan;
    R.Q3_GFit = nan;
    R.fitresult = nan;
    R.gof = nan;
    
    R.D_43 = nan;
    R.D_43_StdDev = nan;
end

PSD_3_d25 = R.PSD_3_d25; 
PSD_3_d50 = R.PSD_3_d50;
PSD_3_d75 = R.PSD_3_d75;
PSD_3_Span = R.PSD_3_Span;

