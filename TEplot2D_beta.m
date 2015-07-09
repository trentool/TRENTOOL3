function TEplot2D_beta(cfg,data)

% TEPLOT2D: Plots a 2D representation of the data returned by
%
%   - single subject analysis (InteractionDelayReconstruction_calculate, 
%     TEsurrogatestats, TEsurrogatestats_ensemble)
%   - group analysis (TEgroup_stats)
%   - binomial testing of single subject data (TEsurrogate_binomstats)
%   - graph correction (TEgraphanalysis)
%
%
% The cfg structure must contain:
%
%   .plottype     = type of data to be plotted (see above):
%                   
%                   'singlesubj' - single subject analysis, you may provide
%                   cfg.statstype and cfg.linktype
%                   'groupstats' - group statistics, you may provide
%                   cfg.statstype and cfg.linktype
%                   'binomstats' - binomial testing of single subject data,
%                   you may provide cfg.statstype and cfg.linktype
%                   'graphanalysis' - results of graph correction
%
%   .statstype    = links to be plotted
%                   'all':         plot all links in the data set;
%                   'corrected':   significant links after correction for 
%                                  multiple comparison;
%                   'uncorrected': significant links without correction for 
%                                  multiple comparison.
%
%   .linktype     = statistics to be plotted 
%                   'sign':     plot significant links only (according to
%                               cfg.statstype);
%                   'pval':     encode significance 1-pval; 
%                   'rawdist':  encode rawdistance; 
%                   'noccur':   number of occurences (for output from 
%                               TEsurrogate_binomstats.m only)
%                   'graphres': results from graph correction(for output  
%                               from TEgraphanalysis.m only)
%
%   .threshold    = provide a numerical value as a threshold for linktypes
%                   'pval' ([0 1]), 'rawdist' (same range as raw values),
%                   'noccur' ([0 no. subjects])
%
%
% Layout options:
%
%   .arrowheadpos = Position of arrowhead: 0 = end of the line; 1 =
%                   centre of the line (default = 1)
%   .arrowcolor   = Color of arrows (default = [1 0 0])
%   .alinewidth   = Linewidth of arrows (default = 2)
%   
%   .head         = Plot a head cartoon: 'on', 'off'
%   .hcolor       = Color of head cartoon (default = [0 0 0])
%   .hlinewidth   = Linewidth of the drawn head, nose and ears
%                    (default = 2)
%
%   .electrodes   = 'on','off','labels','numbers' (default = 'on')
%   .emarker      = Marker symbol for electrodes (default = 'o')
%   .ecolor       = Marker color for electrodes (default = [0 0 0])
%   .emarkersize  = Marker size for electrodes (default = 2)
%   .efontsize    = Font size of electrode labels/numbers (default = 8 pt)
%                   when cfg.electrodes = 'numbers' or 'labels'
%   .efontcolor   = Font color of electrode labels/numbers when
%                   cfg.electrodes = 'numbers' or 'labels' 
%                   (default = [0 0 0])
%   .layout       = specification of the layout which defines how the 
%                   channels are arranged: 
%
%   .colorbar     = plot a colorbar for arrow values; 1 = yes, 0 = no 
%                   (default = 1)
%   .plotinfo     = add plotting information to figure; 1 = yes, 0 = no 
%                   (default = 1)
%
% The layout defines how the channels are arranged. You can specify the
% layout in two ways:
%  - you can provide a pre-computed layout structure (see prepare_layout)
%  - you can give the name of an ascii layout file with extension *.lay
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 2.0 by Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt 2012
%

% CHANGELOG
% 18/12/2014 PW: new statstype 5 for plotting of the graph analysis output
% 18/12/2014 PW: new option to switch offf/on plotting of head contour
% 01/19/2015 PW: removed cfg option 'alpha', wasn't used anymore


%% set default values

% if cfg.statstype == 3
%     if ~isfield(cfg,'alpha');    cfg.alpha = 2;         end;
% end

if ~isfield(cfg,'plottype'); 
    error('TRENTOOL ERROR: Specify at least the field or key "plottype".'); 
end

if ~isfield(cfg,'statstype'); cfg.statstype = 'corrected';   end;
if ~isfield(cfg,'linktype');  cfg.linktype = 'sign';         end;
if ~isfield(cfg,'threshold'); cfg.threshold = -inf;          end;


if ~isfield(cfg,'arrowheadpos'); cfg.arrowheadpos = 1;       end;
if ~isfield(cfg,'arrowcolor');   cfg.arrowcolor  = [0 0 1];  end;
if ~isfield(cfg,'alinewidth');   cfg.alinewidth = 2;         end;
if ~isfield(cfg,'colormap');     cfg.colormap = 'Autumn';    end;

if ~isfield(cfg,'electrodes');   cfg.electrodes = 'on'; end; % on,off,label,numbers or highlights
if ~isfield(cfg,'emarker');      cfg.emarker = 'o';     end;
if ~isfield(cfg,'ecolor');       cfg.ecolor = [0 0 0];  end;
if ~isfield(cfg,'emarkersize');  cfg.emarkersize = 2;   end;
if ~isfield(cfg,'efontsize');    cfg.efontsize = get(0,'DefaultAxesFontSize');end;
if ~isfield(cfg,'efontcolor');   cfg.efontcolor = [0 0 0];  end;

if ~isfield(cfg, 'head');          cfg.head = 'on';      end;
if ~isfield(cfg, 'hcolor');        cfg.hcolor = [0 0 0]; end;
if ~isfield(cfg, 'hlinewidth');    cfg.hlinewidth = 2;   end;

if ~isfield(cfg, 'colorbar');    cfg.colorbar = 1;       end;
if ~isfield(cfg, 'plotinfo');    cfg.plotinfo = 1;       end;

if ~isfield(cfg,'layout')
    error('TRENTOOL ERROR: Specify at least the field or key "layout".');
end;

%% additional parameters

cfg.efsize = get(0,'DefaultAxesFontSize');
cfg.efcolor = [0 0 0];

%% check layout file

% if cfg.statstype == 1; sl='corrected';
% elseif cfg.statstype == 2; sl='uncorrected';
% elseif cfg.statstype == 3; sl='1-pvalue';
% elseif cfg.statstype == 4; sl='raw effect';
% elseif cfg.statstype == 5; sl='graph correction';
% end
% 
% widthmultiplier = 5; % for arrowlines in case of 1-pvalue

% check whether the entry in cfg.layout is a structure (if not assume ascii
% file)
if isstruct(cfg.layout) && all(isfield(cfg.layout, {'pos';'label'}))
    lay = cfg.layout;
    xcoord = lay.pos(:,1);
    ycoord = lay.pos(:,2);
    layout_labels = lay.label;
    index = 1:size(xcoord,1);
else
    % read .lay file
    filename = strcat(cfg.layout,'.lay');
    [index,xcoord,ycoord,zcoord,v5,layout_labels] = textread(filename,'%d %f %f %f %f %s',-1);
end;
% Scale the data to a circle with x-axis and y-axis: -0.45 to 0.45
xcoord = 0.8*((xcoord-min(xcoord))/(max(xcoord)-min(xcoord))-0.5);
ycoord = 0.8*((ycoord-min(ycoord))/(max(ycoord)-min(ycoord))-0.5);


%% check channel combinations

% find channel index from layout file for all used channel combinations
channelcombi = nan(size(data.sgncmb));
for kk=1:size(channelcombi,1)
    channelcombi(kk,1) = find(strcmp(layout_labels,data.sgncmb{kk,1}));
    channelcombi(kk,2) = find(strcmp(layout_labels,data.sgncmb{kk,2}));
end
if any(isnan(channelcombi))
    error('TRENTOOL ERROR: mismatch between channel labels in layout file and data.sgncmb!');
end
     
% find coordinates and indices (ordered like in the layout file)
usedchannels   = unique(channelcombi);
channelxcoords = xcoord(usedchannels);
channelycoords = ycoord(usedchannels);
channelindex   = index(usedchannels);
channelnames   = layout_labels(usedchannels);
channelnames   = strrep(channelnames,'_',' ');

%% find sigificant link (either corrected or uncorrected)

if strcmp(cfg.statstype, 'all')
    linkind = true(size(data.TEpermvalues(:,3),1),1);
elseif strcmp(cfg.statstype, 'corrected')
    linkind = logical(data.TEpermvalues(:,3));
elseif strcmp(cfg.statstype, 'uncorrected')
    linkind = logical(data.TEpermvalues(:,2));
else
    error('TRENTOOL ERROR: please specify either ''all'', ''corrected'', or ''uncorrected'' as cfg.statstype.');
end

% add links flagged during graph analysis
if strcmp(cfg.linktype, 'graphres')
    linkind = linkind | data.TEpermvalues(:,5)>1;
end

%% find statistic to be plotted, apply threshold (sign, p-value, raw TE value, etc.)

colormap(cfg.colormap);
cmap = colormap;

switch cfg.linktype
    case 'sign'
        acolor = repmat(cfg.arrowcolor, sum(linkind), 1);
        minval = 0;       % scaling for color bar axis
        maxval = 1;
        stepsize = 0.5;
    case 'pval'
        linkind = data.TEpermvalues(linkind,1) > cfg.threshold;
        pval    = data.TEpermvalues(linkind,1);
        acolor  = getColorMap(pval,cmap);
        minval = min(pval);       % scaling for color bar axis
        maxval = max(pval);
        stepsize = ceil(length(unique(pval))/5);
    case 'rawdist'
        linkind = data.TEpermvalues(linkind,4) > cfg.threshold;
        rawdist = data.TEpermvalues(linkind,4); 
        acolor = getColorMap(rawdist,cmap);
        minval = min(rawdist);       % scaling for color bar axis
        maxval = max(rawdist);
        stepsize = ceil(length(unique(rawdist))/5);
    case 'noccur'       % works for output of TEsurrogate_binomstats only
        if ~isfield(data, 'occurrences')
            error('TRENTOOL ERROR: To use ''noccur'' as linktype, data has to come from TEsurrogate_binomstats.m.')
        end;
        linkind = data.occurrences .* linkind > cfg.threshold & linkind;    % threshold data if requested
        noccur  = data.occurrences(linkind); 
        acolor  = getColorMap(noccur,cmap);
        minval = min(noccur);       % scaling for color bar axis
        maxval = max(noccur);
        stepsize = ceil(length(unique(noccur))/5);  % scaling for color bar axis
    case 'graphres'     % works for graph corrected data only
        graphres = data.TEpermvalues(linkind,5); 
        acolor = repmat(cfg.arrowcolor,sum(linkind),1);
        %acolor(graphres == 1,:) = cfg.arrowcolor;
        acolor(graphres == 2,:) = repmat([1 0 0], sum(graphres == 2), 1); % 2 = cascade effect
        acolor(graphres == 3,:) = repmat([0 1 0], sum(graphres == 3), 1); % 3 = cascade effect triangle
        acolor(graphres == 4,:) = repmat([1 1 0], sum(graphres == 4), 1); % 4 = common drive link triangle
        minval = 2;
        maxval = 4;
    otherwise 
        error('TRENTOOL ERROR: Please specify ''.linktype''.');
end   

if minval == maxval;
    minval = maxval - stepsize;
end

nlinks = sum(linkind);
X = zeros(nlinks,2);
Y = zeros(nlinks,2);

X(:,1)=xcoord(channelcombi(linkind,1));
X(:,2)=xcoord(channelcombi(linkind,2));
Y(:,1)=ycoord(channelcombi(linkind,1));
Y(:,2)=ycoord(channelcombi(linkind,2));

if nlinks == 0
    warning('TRENTOOL ERROR: No significant links were found.')
end

%% plotting

% Define the outline of the head, ears and nose
rmax=.5;
l     = 0:2*pi/100:2*pi;
tip   = rmax*1.15; base = rmax-.004;
EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
hold on

% set background to white
% set(gcf, 'Color', 'w');

% plot head, ears, and nose
if strcmp(cfg.head, 'on');
    plot(cos(l).*rmax, sin(l).*rmax, 'color', cfg.hcolor , 'Linestyle', '-', 'LineWidth', cfg.hlinewidth);
    
    plot([0.18*rmax;0;-0.18*rmax], [base;tip;base], 'Color', cfg.hcolor , 'LineWidth', cfg.hlinewidth);
    plot( EarX, EarY, 'color', cfg.hcolor , 'LineWidth', cfg.hlinewidth)
    plot(-EarX, EarY, 'color', cfg.hcolor , 'LineWidth', cfg.hlinewidth)
end

% plot electrodes
if strcmp(cfg.electrodes,'on') || strcmp(cfg.electrodes,'labels') || strcmp(cfg.electrodes,'numbers')
    for ii = 1:size(xcoord,1)
        plot(xcoord(ii),ycoord(ii),cfg.emarker,'Color',cfg.ecolor,'Linewidth',cfg.emarkersize)
        
    end
end

% plot arrows
for ll = 1:nlinks
    TEarrow(...
        X(ll,1), X(ll,2), Y(ll,1), Y(ll,2), ...
        cfg.alinewidth, acolor(ll,:), cfg.arrowheadpos ...
        );
end

% Plot labels of used channels
if strcmp(cfg.electrodes,'labels')  
    for ii=1:size(usedchannels,1)
        h = text(channelxcoords(ii)+0.03,channelycoords(ii)+0.03,channelnames(ii));
        set(h,'FontSize',cfg.efsize,'Color',cfg.efcolor);
    end    
elseif strcmp(cfg.electrodes,'numbers') 
    for ii=1:size(usedchannels,1)
        h = text(channelxcoords(ii)+0.01,channelycoords(ii)+0.01,num2str(channelindex(ii)));
        set(h,'FontSize',cfg.efsize,'Color',cfg.efcolor);
    end   
end


axis off;
set(gca, 'YLim', [-0.5 .6],...
    'XLim', [-0.6 .6],...
    'Visible', 'off',...
    'Box','off');

if cfg.colorbar
    
    if ~strcmp(cfg.linktype,'graphres')
        cbh = colorbar;        
        caxis([minval maxval]);
        set(cbh, ...
            'YTick',minval:stepsize:maxval, ...
            'YTickLabel',minval:stepsize:maxval)
    else
        h=zeros(1,4);
        h(1) = plot(100,100,'Color',cfg.arrowcolor);
        h(2) = plot(100,100,'Color',[1 0 0]);
        h(3) = plot(100,100,'Color',[0 1 0]);
        h(4) = plot(100,100,'Color',[1 1 0]);
        legend(h, ...
            'Significant TE','Cascade Effect','Simple Cascade Effect','Simple Common Drive', ...
            'Location','NorthEast');
    end
end  

if cfg.plotinfo    
    text(-0.8,-0.45, ['linktype: ' num2str(cfg.linktype)])
    if cfg.threshold == -inf
        text(-0.8,-0.5, ['threshold: none'])
    else
        text(-0.8,-0.5, ['threshold: ' num2str(cfg.threshold)])
    end
    text(-0.8,-0.55, ['multiple comp.: ' num2str(cfg.statstype)]);
end

title(gca, ...
    strcat(['Transfer Entropy - ' cfg.plottype ]),'FontSize',20);
hold off
end


function acolor = getColorMap(data,cmap)
    ncolors = size(cmap,1);
    [binning,binind] = histc(data, linspace( min(data), max(data), ncolors ));
    acolor = cmap(binind,:);
end
