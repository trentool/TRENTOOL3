function handles = TEarrow(x1,x2,y1,y2,linewidth,color,arrowpos)

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation;
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY;
%
% Version 1.0 by Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt 2009
%


% =============================================
% calculate the arrow head coordinates
% =============================================
den         = x2 - x1 + eps;
teta        = atan((y2-y1)/den)+pi*(x2<x1)-pi/2;
cs          = cos(teta);
ss          = sin(teta);
R           = [cs -ss;ss cs];
%linelength  = sqrt((y2-y1)^2+(x2-x1)^2);
headlength  = .05;%min(linelength*alpha,maxlength);
headwidth   = .025;%min(linelength*beta,maxlength);
if arrowpos == 1
    x0          = x2*cs + y2*ss;
    y0          = -x2*ss + y2*cs;
elseif arrowpos == 2
    x0          = 0.5*(x1+x2)*cs + 0.5*(y1+y2)*ss;
    y0          = -0.5*(x1+x2)*ss + 0.5*(y1+y2)*cs;
end
coords      = R*[x0 x0+headwidth/2 x0-headwidth/2; y0 y0-headlength y0-headlength];

% =============================================
% plot arrow  (= line + patch of a triangle)
% =============================================
%xf,yf,0.5*(xt+xf),0.5*(yt+yf)
h1          = plot([x1,x2],[y1,y2],'k','Linewidth',linewidth,'Color',color);
h2          = patch(coords(1,:),coords(2,:),color);

% =============================================
% return handles
% =============================================
handles = [h1 h2];
end
