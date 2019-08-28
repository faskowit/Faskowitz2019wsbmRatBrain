
clc
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = strcat(PROJECT_DIR,'/data/processed/',OUTPUT_STR, '_', GRID_RUN,'_evalGenReps.mat') ;
load(loadName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig stuff

FIGURE_NAME = 'figGener' ;
outputdir = strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/');
mkdir(outputdir) 

writeit = 1 ;

fontsize = 16 ;

% general graphics, this will apply to any figure you open
% (groot is the default figure object).
set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesFontSize', 14, ...
'DefaultAxesFontName', 'Arial', ...
'DefaultLineLineWidth', 1, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 16, ...
'DefaultTextFontName', 'Arial', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KS

figure

% calculate num bins
X = [ mean(eval_wsbm_K,2) ; mean(eval_mod_K,2) ] ;
binSize = 3.5*std(X(:))*numel(X)^(-1/3) ;

% if third needed, we can add
cmap = [0    0.4470    0.7410 ;
    0.8500    0.3250    0.0980 ;
    0.8500    0.3250    0.3920 ];

histogram(mean(eval_wsbm_K(:,1:5),2),...
    'normalization','probability',...
    'FaceColor',cmap(1,:),'EdgeAlpha',0.01,...
    'BinWidth',binSize) 
hold 
histogram(mean(eval_mod_K(:,1:5),2),...
    'normalization','probability',...
    'FaceColor',cmap(2,:),'EdgeAlpha',0.01,...
    'BinWidth',binSize)

lg = legend('WSBM','Modular','Location','northwest') ;
lg.FontSize = fontsize ;

axis square
xl = xlabel('Mean KS energy') ;
yl = ylabel('Normalized frequency') ;
xl.FontSize = fontsize ;
yl.FontSize = fontsize ;

%set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.55, 0.5]);
pbaspect([1 1 1])

% tight
%tightfig

if writeit
    fileName = strcat('gen_wsbmNmod.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot each with their null

figure

combovec = [ mean(eval_wsbm_K,2) ; mean(eval_wsbmRand_K,2) ; 
    mean(eval_mod_K,2) ; mean(eval_modRand_K,2) ];

x_lim = [ min(combovec)*0.9 max(combovec)*1.1 ] ;
y_lim = [ 0 0.17 ] ;

% tight subplot
sp = tight_subplot(2,1,.01,[.1 .025],[.15 .15]);

axes(sp(1))

hh = histogram(mean(eval_wsbm_K,2),'normalization','probability','FaceColor',cmap(1,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
hh.BinWidth = hh.BinWidth * 2 ;
hold
hh = histogram(mean(eval_wsbmRand_K,2),'normalization','probability','FaceColor',cmap(1,:),'FaceAlpha',0.25,'EdgeAlpha',0.05)
hh.BinWidth = hh.BinWidth * 2 ;

set(gca,'XTickLabel',[])

xlim(x_lim)
ylim(y_lim)

lg = legend('intact','null') ;
lg.FontSize = fontsize - 2 ;

yl = ylabel('Normalized frequency') ;
yl.FontSize = fontsize ;

%set(gca,'FontSize',fontsize)

axes(sp(2))

hh = histogram(mean(eval_mod_K,2),'normalization','probability','FaceColor',cmap(2,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
hh.BinWidth = hh.BinWidth * 2 ;
hold
hh = histogram(mean(eval_modRand_K,2),'normalization','probability','FaceColor',cmap(2,:),'FaceAlpha',0.25,'EdgeAlpha',0.05)
hh.BinWidth = hh.BinWidth * 2 ;

xlim(x_lim)
ylim(y_lim)

lg = legend('intact','null') ;
lg.FontSize = fontsize - 2 ;

yl = ylabel('Normalized frequency') ;
yl.FontSize = fontsize ;
xl = xlabel('Mean KS energy') ;
xl.FontSize = fontsize ;

%set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.3, 0.5]);

% tight
tightfig

if writeit
    fileName = strcat('gen_wsbmNmod_wNull.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%% plot each statistic

%     evalFuncs{1} = @(A) sum(triu(A,1),2) + sum(tril(A,1),1)';
%     evalFuncs{2} = @(A) sum(triu(A ~= 0,1),2) + sum(tril(A ~= 0,1),1)';
%     evalFuncs{3} = @(A) clustering_coef_wd(A);
%     evalFuncs{4} = @(A) betweenness_wei(1 ./ A);
%     evalFuncs{5} = @(A) eigenvector_centrality_und(A) ;

figure

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.95, 0.50]);

% tight subplot
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(2,5,[ .035 .04 ],[.12 .033],[.15 .1]);

measure_names = {'Strength' 'BDegree' 'Clustering' 'Betwness' 'BetwnessBin'} ;
measure_short_names = { 's' 'bd' 'c' 'b' 'bb' } ;

x_lim = cell(5,1) ;

for idx = 1:5 
   
    combo = [ eval_wsbmRand_K(:,idx) ; eval_modRand_K(:,idx) ; ...
        eval_wsbm_K(:,idx) ; eval_mod_K(:,idx) ] ;
    x_lim{idx} = [ 0 (max(combo) .* 1.02) ] ;
    
end

for idx = 1:5
       
    axes(sp(idx))
    
    hh = histogram(eval_wsbm_K(:,idx),'BinMethod','sturges',...
        'normalization','probability','EdgeAlpha',0.05,'FaceAlpha',0.6,...
        'FaceColor',cmap(1,:),'BinWidth',binSize)
    hh.BinWidth = hh.BinWidth * 2.6 ;
    hold
    hh = histogram(eval_mod_K(:,idx),'BinMethod','sturges',...
        'normalization','probability','EdgeAlpha',0.05,'FaceAlpha',0.6,...
        'FaceColor',cmap(2,:),'BinWidth',binSize)
    hh.BinWidth = hh.BinWidth * 2.6 ;
    axis square
    tl = title(measure_names{idx},'FontWeight','normal') ;
    tl.FontSize = fontsize ;
    
%     if idx == 1
%        ylabel('Normalized frequency') 
%     end
        
    [~,p,ci,stat] = ttest2(eval_mod_K(:,idx),eval_wsbm_K(:,idx),'Vartype','unequal') ;
    
    % xlabel(strcat('EMD({\it ',measure_short_names{idx},'})'))
   set(gca,'XTickLabel',[])
    
    %xlim([0 60])
    xlim(x_lim{idx}) 
    %set(gca,'FontSize',16)

end

for idx = 1:5
   
    axes(sp(idx+5))
    
    hh = histogram(eval_wsbmRand_K(:,idx),'normalization','probability',...
        'FaceColor',cmap(1,:),'FaceAlpha',0.25,'EdgeAlpha',0.01)
    hh.BinWidth = hh.BinWidth * 2.5 ;
    hold
    hh = histogram(eval_modRand_K(:,idx),'normalization','probability',...
        'FaceColor',cmap(2,:),'FaceAlpha',0.25,'EdgeAlpha',0.01)
    hh.BinWidth = hh.BinWidth * 2.5 ;
    axis square
    
    % and plot a line for the empirical
    ylimits = ylim ;
    
    tmp = median(mean(eval_wsbm_K(:,idx)),2) ;
    plot([ tmp tmp ],...
        [ylimits(1) ylimits(2)*0.95],...
        'Color',[cmap(1,:) 0.9],'LineWidth',1.5)
    
    tmp = median(mean(eval_mod_K(:,idx)),2) ;
    plot([ tmp tmp ],...
        [ylimits(1) ylimits(2)*0.95],...
        'Color',[cmap(2,:) 0.9],'LineWidth',1.5)
    
    ylim(ylimits)
    
    %x_lim{idx} = xlim() ;
    xlim(x_lim{idx}) ;
    
    xl = xlabel(strcat('KS({\it ',measure_short_names{idx},'})'))
    xl.FontSize = fontsize ;
    
end

% set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.99, 0.50]);
%tightfig

if writeit
    fileName = strcat('gen_wsbmNmod_wNull_splitout.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    
    
    fileName = strcat('gen_wsbmNmod_wNull_splitout2.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    export_fig(gcf,ff,'-r550')
    
    close(gcf)
end
