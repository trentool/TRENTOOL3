function TEplot2D(cfg,data)

%
%
%
% You can use the function as following:
%         TEplot2D(cfg,Data) 
%
% cfg.statstype    = 1: corrected; 2:uncorrected; 3: 1-pval; 4:
%                     rawdistance; 5:results graph correction
% cfg.arrowpos      = Position of arrowhead: 1 = end of the line; 2 =
%                     centre of the line (default = 2)
% cfg.arrowcolorpos = Color of arrows (default = [1 0 0])
% cfg.arrowcolorneg = Color of arrows in case of negative mean distances
%                     (default = [0 0 1])
% cfg.alinewidth    = Linewidth of arrows in case of significance
%                     (default = 2)
% cfg.electrodes    = 'on','off','labels','numbers','highlights'
%                     (default = 'on')
% cfg.hcolor        = Color of head cartoon (default = [0,0,0])
% cfg.hlinewidth    = number, Linewidth of the drawn head, nose and ears
%                     (default = 2)
% cfg.emarker       = Marker symbol (default = 'o')
% cfg.ecolor        = Marker color (default = [0 0 0] (black))
% cfg.emarkersize   = Marker size (default = 2)
% cfg.efontsize     = Font size of electrode labels/numbers (default = 8 pt)
%                     when cfg.electrodes = 'numbers' or 'labels'
% cfg.efontcolor    = Font color of electrode labels/numbers when
%                     cfg.electrodes = 'numbers' or 'labels'
%                     (default = [0 0 0])
% cfg.hlmarker        = Highlight marker symbol (default = 'o')
% cfg.hlcolor         = Highlight marker color (default = [1 0 0] (red))
% cfg.hlmarkersize    = Highlight marker size (default = 4)
% cfg.layout          = specification of the layout which defines how the 
%                       channels are arranged: 
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



if cfg.statstype == 3
    if ~isfield(cfg,'alpha');    cfg.alpha = 2;         end;
end


if ~isfield(cfg,'arrowpos');     cfg.arrowpos = 2;      end;
if ~isfield(cfg,'arrowcolorpos');   cfg.arrowcolorpos  = [1 0 0];     end;
if ~isfield(cfg,'arrowcolorneg');   cfg.arrowcolorneg  = [0 0 1];     end;
if ~isfield(cfg,'alinewidth');   cfg.alinewidth = 2;     end;

if ~isfield(cfg,'electrodes');   cfg.electrodes = 'on'; end; % on,off,label,numbers or highlights
if ~isfield(cfg,'showlabels'); % for compatibility with OLDSTYLE
    cfg.showlabels = '';
%else
%    cfg.electrodes = '';
end;

if ~isfield(cfg,'emarker');      cfg.emarker = 'o';     end;
if ~isfield(cfg,'ecolor');       cfg.ecolor = [0 0 0];  end;
if ~isfield(cfg,'emarkersize');  cfg.emarkersize = 2;   end;
if ~isfield(cfg,'efontsize');    cfg.efontsize = get(0,'DefaultAxesFontSize');end;
if ~isfield(cfg,'efontcolor');   cfg.efontcolor = [0 0 0];  end;

if ~isfield(cfg,'plothead');     cfg.plothead = 1;      end;
if ~isfield(cfg,'hlmarker');     cfg.hlmarker = 'o';    end;
if ~isfield(cfg,'hlcolor');      cfg.hlcolor = [1 0 0]; end;
if ~isfield(cfg,'hlmarkersize'); cfg.hlmarkersize = 4;  end;

if ~isfield(cfg, 'hcolor');        cfg.hcolor = [0 0 0];     end;
if ~isfield(cfg, 'hlinewidth');    cfg.hlinewidth = 2;       end;


if ~isfield(cfg,'layout')
    error('TRENTOOL ERROR: Specify at least the field or key "layout".');
end;

if isfield(cfg,'electrod')
    cfg.electrodes = lower(cfg.electrod);
    cfg            = rmfield(cfg,'electrod');
end;
if isfield(cfg,'electcolor')
    cfg.ecolor = cfg.electcolor;
    cfg        = rmfield(cfg,'electcolor');
end;
if isfield(cfg,'emsize')
    cfg.emarkersize = cfg.emsize;
    cfg             = rmfield(cfg,'emsize');
end;
if isfield(cfg,'headcolor')
    cfg.hcolor = cfg.headcolor;
    cfg        = rmfield(cfg,'headcolor');
end;
if isfield(cfg,'efontsize')
    cfg.efsize = cfg.efontsize;
    cfg        = rmfield(cfg,'efontsize');
end;
if isfield(cfg,'efontcolor')
    cfg.efcolor = cfg.efontcolor;
    cfg        = rmfield(cfg,'efontcolor');
end;

if cfg.statstype == 1; sl='corrected';
elseif cfg.statstype == 2; sl='uncorrected';
elseif cfg.statstype == 3; sl='1-pvalue';
elseif cfg.statstype == 4; sl='raw effect';
elseif cfg.statstype == 5; sl='graph correction';
end

widthmultiplier = 5; % for arrowlines in case of 1-pvalue

% check whether the entry in cfg.layout is a structure (if not assume ascii
% file)
if isstruct(cfg.layout) && all(isfield(cfg.layout, {'pos';'label'}))
    lay = cfg.layout;
    xcoord = lay.pos(:,1);
    ycoord = lay.pos(:,2);
    allchannels = lay.label;
    index = 1:size(xcoord,1);
else
    % read .lay file
    filename = strcat(cfg.layout,'.lay');
    [index,xcoord,ycoord,zcoord,v5,allchannels] = textread(filename,'%d %f %f %f %f %s',-1);
end;
% Scale the data to a circle with x-axis and y-axis: -0.45 to 0.45
xcoord = 0.8*((xcoord-min(xcoord))/(max(xcoord)-min(xcoord))-0.5);
ycoord = 0.8*((ycoord-min(ycoord))/(max(ycoord)-min(ycoord))-0.5);


channelcombi=zeros(size(data.sgncmb));
for cc = 1:2
    channelcmp=data.sgncmb(:,cc);
    usedcounter = 0;
    for jj=1:size(allchannels,1)
        for kk=1:size(channelcmp,1)
            if strcmp(allchannels{jj},channelcmp{kk})
                channelcombi(kk,cc)=jj;
                usedcounter = usedcounter + 1;
            end
        end
    end
    if usedcounter ~= size(channelcombi,1)
        error('TRENTOOL ERROR: mismatch between cfg.sgncmb and data.label');
    end
end

usedchannels=unique(sort(reshape(channelcombi,size(channelcombi,1)*size(channelcombi,2),1)));
channelxcoords=xcoord(usedchannels);
channelycoords=ycoord(usedchannels);
channelindex=index(usedchannels);
channelnames=allchannels(usedchannels);

flag1=0;flag2=0;

counter=1;
for ii=1:size(data.TEpermvalues,1),
    switch lower(cfg.statstype)
        case 1,
            % corrected
            if data.TEpermvalues(ii,3)==1
                if flag1==0; flag1=1; end;
                X(1,counter)=xcoord(channelcombi(ii,1));
                X(2,counter)=xcoord(channelcombi(ii,2));
                Y(1,counter)=ycoord(channelcombi(ii,1));
                Y(2,counter)=ycoord(channelcombi(ii,2));
                if data.TEpermvalues(ii,4)<=0
                    arrowvalue(counter)=0;
                else
                    arrowvalue(counter)=1;
                end
                counter=counter+1;
            end
        case 2,
            % uncorrected
            if data.TEpermvalues(ii,2)==1
                if flag2==0; flag2=1; end;
                X(1,counter)=xcoord(channelcombi(ii,1));
                X(2,counter)=xcoord(channelcombi(ii,2));
                Y(1,counter)=ycoord(channelcombi(ii,1));
                Y(2,counter)=ycoord(channelcombi(ii,2));
                if data.TEpermvalues(ii,4)<=0
                    arrowvalue(counter)=0;
                else
                    arrowvalue(counter)=1;
                end
                counter=counter+1;
                
            end
        case 3,
            % 1-pvalue
            %if data.TEpermvalues(ii,1)<=cfg.alpha
            X(1,counter)=xcoord(channelcombi(ii,1));
            X(2,counter)=xcoord(channelcombi(ii,2));
            Y(1,counter)=ycoord(channelcombi(ii,1));
            Y(2,counter)=ycoord(channelcombi(ii,2));
            alinewidth(ii) = widthmultiplier*(1.0-data.TEpermvalues(ii,1)+0.0000000001);
            if data.TEpermvalues(ii,4)<=0
                arrowvalue(counter)=0;
            else
                arrowvalue(counter)=1;
            end
            counter=counter+1;
            %end
        case 4,
            % mean dist
            X(1,counter)=xcoord(channelcombi(ii,1));
            X(2,counter)=xcoord(channelcombi(ii,2));
            Y(1,counter)=ycoord(channelcombi(ii,1));
            Y(2,counter)=ycoord(channelcombi(ii,2));
            alinewidth(ii) = widthmultiplier*(data.TEpermvalues(ii,4));
            if data.TEpermvalues(ii,4)<=0
                arrowvalue(counter)=0;
            else
                arrowvalue(counter)=1;
            end
            counter=counter+1;
        case 5,
            % results after graph correction
            if data.TEpermvalues(ii,2)==1 || data.TEpermvalues(ii,5)>1
                X(1,counter)=xcoord(channelcombi(ii,1));
                X(2,counter)=xcoord(channelcombi(ii,2));
                Y(1,counter)=ycoord(channelcombi(ii,1));
                Y(2,counter)=ycoord(channelcombi(ii,2));
                alinewidth(ii) = widthmultiplier*(data.TEpermvalues(ii,4));
                if data.TEpermvalues(ii,5)<=1
                    arrowvalue(counter)=0;
                else
                    arrowvalue(counter)=1;
                end
                counter=counter+1;
            end
    end
end

if flag1 == 0 && cfg.statstype == 1
    error('TRENTOOL error: No corrected significances are found. Try lower cfg.statstype value. (see help)')
elseif flag2 == 0 && cfg.statstype == 2
    error('TRENTOOL error: No uncorrected significances are found. Try lower cfg.statstype value. (see help)')
end


% Define the outline of the head, ears and nose:
rmax=.5;
l     = 0:2*pi/100:2*pi;
tip   = rmax*1.15; base = rmax-.004;
EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
hold on

% Plot head, ears, and nose:
if cfg.plothead
    plot(cos(l).*rmax, sin(l).*rmax, 'color', cfg.hcolor , 'Linestyle', '-', 'LineWidth', cfg.hlinewidth);
    
    plot([0.18*rmax;0;-0.18*rmax], [base;tip;base], 'Color', cfg.hcolor , 'LineWidth', cfg.hlinewidth);
    plot( EarX, EarY, 'color', cfg.hcolor , 'LineWidth', cfg.hlinewidth)
    plot(-EarX, EarY, 'color', cfg.hcolor , 'LineWidth', cfg.hlinewidth)
    %hold off
end

%hold on
% Plot all electrodes
if strcmp(cfg.electrodes,'on')||strcmp(cfg.showlabels,'markers')
    %hold on
    for ii = 1:size(xcoord,1)
        plot(xcoord(ii),ycoord(ii),cfg.emarker,'Color',cfg.ecolor,'Linewidth',cfg.emarkersize)
        
    end
end


%plot used channels
if strcmp(cfg.electrodes,'on')||strcmp(cfg.electrodes,'highlights')||strcmp(cfg.electrodes,'labels')||strcmp(cfg.electrodes,'numbers')||strcmp(cfg.showlabels,'markers')
    %hold on
    for jj = 1:size(channelxcoords,1)
        plot(channelxcoords(jj),channelycoords(jj),cfg.hlmarker,'Color',cfg.hlcolor,'Linewidth',cfg.hlmarkersize)
    end
end

% plot arrows
if cfg.statstype <= 2 || cfg.statstype == 5
    for kk = 1:size(X,2)
        if arrowvalue(kk)==0
            TEarrow(X(1,kk),X(2,kk),Y(1,kk),Y(2,kk),cfg.alinewidth,cfg.arrowcolorneg,cfg.arrowpos);
        else
            TEarrow(X(1,kk),X(2,kk),Y(1,kk),Y(2,kk),cfg.alinewidth,cfg.arrowcolorpos,cfg.arrowpos);
        end
    end
else
    for kk = 1:size(X,2)
        if arrowvalue(kk)==0
            TEarrow(X(1,kk),X(2,kk),Y(1,kk),Y(2,kk),cfg.alinewidth,cfg.arrowcolorneg,cfg.arrowpos);
            TEarrow(0.58,0.58,-0.2,0.2,widthmultiplier,cfg.arrowcolorneg,cfg.arrowpos);
        else
            TEarrow(X(1,kk),X(2,kk),Y(1,kk),Y(2,kk),cfg.alinewidth,cfg.arrowcolorpos,cfg.arrowpos);
            TEarrow(0.58,0.58,-0.2,0.2,widthmultiplier,cfg.arrowcolorpos,cfg.arrowpos);
        end
    end
    
end

% Plot labels of used channels
if strcmp(cfg.electrodes,'labels') || strcmp(cfg.showlabels,'yes')
    for ii=1:size(usedchannels,1)
        %h=text(channelxcoords(ii)+0.01,channelycoords(ii)+0.01,channelnames(ii));
        h=text(channelxcoords(ii)+0.03,channelycoords(ii)+0.03,channelnames(ii));
        set(h,'FontSize',cfg.efsize,'Color',cfg.efcolor);
    end
    
elseif strcmp(cfg.electrodes,'numbers') || strcmp(cfg.showlabels,'numbers')
    for ii=1:size(usedchannels,1)
        h=text(channelxcoords(ii)+0.01,channelycoords(ii)+0.01,num2str(channelindex(ii)));
        set(h,'FontSize',cfg.efsize,'Color',cfg.efcolor);
    end
    
end


hold off
axis off;
set(gca, 'YLim', [-0.5 .6],...
    'XLim', [-0.6 .6],...
    'Visible', 'off',...
    'Box','off');

title(gca,strcat(['Connectivity arrows: ',sl]),'FontSize',20);


end
