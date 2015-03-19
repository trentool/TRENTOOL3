% function PlotTE
function []=InteractionDelayReconstruction_plotting(cfg)

% InteractionDelayReconstruction_plotting
% provides a detailed plot of TE and (1-p) versus perdictiontime u
% created from the intermediate results of InteractionDelayReconstruction_calculate
% which are stored in files with names that follow patterns like
% 'ORIGINALDATASETNAME_FILEIDOUT_u_*_TIMEINFO_TEpermtest_output.mat'
% (when you check the directory where your output is stored look for the
% files that end in 'TEpermtest_output.mat' it shoud be pretty clear how to set up the pattern)
%
% INPUTS:
%
%   cfg                   configuration structure  with fields
%
%   cfg.pattern           the Filename pattern with MATLAB compatible 
%                         wildcards (such as ?,*)
%
%   cfg.directory         The directory where the files conatining the
%                         results of TEpermtest are stored
%   cfg.graphics    =     'pp_ready' or 'view'
%                         The configuration for paper ready/ testview,
%                         default 'view'
%   cfg.scaletype   =     'log' or 'lin' - the predictiontime u scale
%                          (default 'log')
%
% 2012 Michael Wibral based on code by Viola Priesemann


% CHANGELOG:
%
% 2012.02.27 NP: plot for multiple channel combination;
% 2012.07.05 NP: adding scaletype for ploting
% 2014.12.01 PW: added correction for negative TE values (dividing by the
% max abs value doesn't work for negative values)
% 2014.12.02 PW: compatibility with new delay reconstruction
% 2015/01/30 PW: changed the way figures are handled (now called without the indexing),
%		 this caused problems with multiple calls to the function

%set defaluts
if  ~isfield(cfg,'directory')
    cfg.directory=pwd;
end;

if ~isfield(cfg,'scaletype')
    cfg.scaletype = 'log';
end;

if ~strcmp(cfg.directory(end),'/')
    cfg.directory=strcat(cfg.directory,'/');
end;

if ~isfield(cfg, 'graphics' )
    cfg.graphics = 'view';
end;

ddir = dir(strcat(cfg.directory,cfg.pattern));

if isempty(ddir)
    error('TRENTOOL: Didn''t find any entries for the given directory and pattern.');
end

TEmat = [];
for k=1:length(ddir)
   load(strcat(cfg.directory,ddir(k).name)); 
   u(k)=TEpermtest.cfg.predicttime_u(1); 
   %TE{k}=squeeze(TEpermtest.TEpermvalues); %squeeze not necessary in TRENTOOL2
   TEmat=cat(2,TEmat,mean(TEpermtest.TEmat,2));
end

[u,idx] = sort(u,'ascend');
sgs = TEpermtest.sgncmb;
pval  = zeros(size(sgs,1),length(idx));
TEval = zeros(size(sgs,1),length(idx));
TEsig = zeros(size(sgs,1),length(idx));
TEsig2= zeros(size(sgs,1),length(idx));
TEvolC= zeros(size(sgs,1),length(idx));

for uu = 1:length(idx)
    %pval(:,uu)  = TE{uu}(:,1); 
    %TEval(:,uu) = TE{uu}(:,4); 
    TEval(:,uu) = TEmat(:,uu); 
    %TEsig(:,uu) = TE{uu}(:,2);
    %TEsig2(:,uu)= TE{uu}(:,3);
    %TEvolC(:,uu)= TE{uu}(:,5);
end

% ch_nr = length(sgs);
% j = 0;
% for k=1:length(sgs); 
%     ii=1;
%     for ku = idx
%         pval(k,ii)  = TE{ku}(1*ch_nr-(ch_nr-1)+j); % was (k,1)
%         TEval(k,ii) = TE{ku}(4*ch_nr-(ch_nr-1)+j); % ...
%         TEsig(k,ii) = TE{ku}(2*ch_nr-(ch_nr-1)+j);
%         TEsig2(k,ii)= TE{ku}(3*ch_nr-(ch_nr-1)+j);        
%         TEvolC(k,ii)= TE{ku}(5*ch_nr-(ch_nr-1)+j); 
%         ii=ii+1;
%     end
%     j = j+1;
% end

% shift all values if there are TE values smaller than 0 (otherwise
% dividing by max abs doesn't work)
if any(TEval < 0)
    TEval = TEval + abs(min(min(TEval)));
end

if strcmp(cfg.graphics,'view')
    
    cmap = jet(7);
    cc=3; dy = 1.3; df=0;
    % for k=1:6, figure(k+df), hold off, end

    for k=1:length(sgs)
        %figure(k+df), %  hold off
        figure
        if strcmp(cfg.scaletype,'log') 
            semilogx(u,1-pval(k,:),'--','Color',cmap(cc,:)), hold all
        elseif strcmp(cfg.scaletype,'lin')
            plot(u,1-pval(k,:),'--','Color',cmap(cc,:)), hold all
        else error('TRENTOOL error: chose a correct plot scaletype');
        end
        
        plot(u,TEval(k,:)/max(abs(TEval(k,:))),'-o','Color',cmap(cc,:), 'LineWidth',2)
        plot(u,ones(length(u),1)*dy,'.','Color',[0.4,0.4,0.4])
        plot(u,1./TEsig(k,:)*dy,'*k')
        plot(u,1./TEsig2(k,:)*dy,'*g')
        plot(u,1./TEvolC(k,:)*dy,'*r')
        title (strcat('signal comb ', sgs(k,1), ' -> ', sgs(k,2), ' maxTE @ u= ', num2str(u(TEval(k,:) == max(TEval(k,:))))));
        h = legend('1-p','TE/max(TE)','calculated','sign','signB','Volume Conduct');
        set(h,'box','off','Color', 'none','Location','SouthEast')
        xlabel('pred time u in ms')
    end
elseif strcmp(cfg.graphics,'pp_ready')

    cc=3; dy = 1.3; df=0;
    for k=1:length(sgs)
        %figure(k+df), %  hold off
        figure
        if strcmp(cfg.scaletype,'log') 
            semilogx(u,TEval(k,:)/max(abs(TEval(k,:))),'o-','Color','k','LineWidth',2,'MarkerSize',6), hold all
        elseif strcmp(cfg.scaletype,'lin')
            plot(u,TEval(k,:)/max(abs(TEval(k,:))),'o-','Color','k','LineWidth',2,'MarkerSize',6), hold all
        else error('TRENTOOL error: chose a correct plot scaletype');
        end
       
        plot(u,ones(length(u),1)*dy,'.','Color','k','MarkerSize',4)
        plot(u,1./TEsig(k,:)*dy,'*k','MarkerSize',5.5)
        plot(u,1./TEsig2(k,:)*dy,'dk','MarkerSize',6.5)
        plot(u,1./TEvolC(k,:)*dy,'xk','MarkerSize',5.5)
        
%    title (['signal comb ', sgs(k,:), 'maxTE ', num2str(max(TEval(k,:)))])
        %h = legend('1-p','TE/max(TE)','calculated','sign','signB','Volume Conduct');
        %set(h,'box','off','Color', 'none','Location','SouthEast')
        ylabel('TE/max(TE)')
        xlabel('pred time u in ms')
        %grid(gca,'minor');
        %set(gca,'XGrid','on');
    end
else fprintf('no method %s is implemented\n',cfg.graphics);
end
