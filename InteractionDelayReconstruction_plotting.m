function []=InteractionDelayReconstruction_plotting(cfg, TEpermtest)

% FUNCTION INTERACTIONDELAYRECONSTRUCTION_PLOTTING
%
% plots raw TE values versus assumed delay u from the output of 
% InteractionDelayReconstruction_calculate.
%
% INPUTS:
%   TEpermtest            output structure with TE estimates from 
%                         InteractionDelayReconstruction_calculate
%
%   cfg                   configuration structure  with fields
%
%   cfg.standardize =     'yes' or 'no' - shift and z-standardize TE values
%                         for individual signal combinations
%                          (default = 'no')
%
%   cfg.scaletype   =     'log' or 'lin' - the predictiontime u scale
%                          (default = 'lin')
%
%   cfg.ch_per_fig  =     no. channels that are plotted within one Figure,
%                         if the number of channels in the data set is
%                         higher, each channel is plotted into a new Figure
%                         (default = 8)
%
% Version 1.1 by Patricia Wollstadt Michael Wibral, Viola Priesemann
%
% Frankfurt 2015

% CHANGELOG:
%
% 2012.02.27 NP: plot for multiple channel combination;
% 2012.07.05 NP: adding scaletype for ploting
% 2014.12.01 PW: added correction for negative TE values (dividing by the
% max abs value doesn't work for negative values)
% 2014.12.02 PW: compatibility with new delay reconstruction
% 2015/01/30 PW: changed the way figures are handled (now called without the indexing),
%		 this caused problems with multiple calls to the function

%% set defaluts

if ~isfield(cfg,'scaletype');    cfg.scaletype = 'lin'; end;
if ~isfield(cfg,'standardize');  cfg.standardize = 'no'; end;
if ~isfield(cfg,'ch_per_fig');   cfg.ch_per_fig = 8; end;

%%
TEmat = TEpermtest.TEbyU;
uvec  = ...
    TEpermtest.TEprepare.cfg.predicttimemin_u: ...
    TEpermtest.TEprepare.cfg.predicttimestepsize: ...
    TEpermtest.TEprepare.cfg.predicttimemax_u;
optu  = TEpermtest.TEpermvalues(:,end);
n_channels = size(TEmat, 1);
n_u        = size(TEmat, 2);

% shift all values if there are TE values smaller than 0 (otherwise
% dividing by max abs doesn't work)

if strcmp(cfg.standardize, 'yes')
    if any(TEmat < 0)
        TEmat = TEmat + abs(min(min(TEmat)));
    end

    TEmat = TEmat ./ repmat(max(TEmat,[],2), 1, n_u);
    
    maxTE = ones(n_channels,1);
    ylab  = 'TE [z-score]';
else
    maxTE = max(TEmat, [], 2);
    ylab  = 'TE [arb. units]';
end

%% get channel labels for legend

labels = cell(1, n_channels);

for c=1:n_channels
    labels{c} = [...
        strrep(TEpermtest.sgncmb{c,1}, '_', ' ') ... 
        ' - ' ...
        strrep(TEpermtest.sgncmb{c,2}, '_', ' ')];
end

%%
if n_channels < cfg.ch_per_fig
    figure
    hold all
    plot(uvec,TEmat)
    scatter(optu, maxTE, 'r')
    h = gca;
    legend(labels)
    xlabel('u [ms]'); ylabel(ylab);
else
    h = zeros(1, n_channels);
    for c=1:n_channels
        figure
        plot(uvec,TEmat(c,:))
        h(c) = gca;
        legend(labels{c})
        xlabel('u [ms]'); ylabel(ylab);
    end
end

set(h, ...
    'Xlim', [uvec(1)-1 uvec(end)+1], ...
    'Box', 'on' ...
    )

end
