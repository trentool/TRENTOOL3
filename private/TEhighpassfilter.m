function [datafiltered] = TEhighpassfilter(cfg,dat)

% TEhighpassfilter performs an high pass filter on the data.

% The subfunctions butter, bilinear, postpad, sftrans and filtfilt are used
% from the Fieldtrip toolbox.


[B, A] = butter(6, max((1000/cfg.TR)/2)/(1/cfg.hpfreq), 'high');
datafiltered = filtfilt(B, A, dat')';


    function [a, b, c, d] = butter (n, W, varargin)
        
        if (nargin>4 || nargin<2) || (nargout>4 || nargout<2)
            usage ('[b, a] or [z, p, g] or [a,b,c,d] = butter (n, W [, "ftype"][,"s"])');
        end
        
        % interpret the input parameters
        if (~(length(n)==1 && n == round(n) && n > 0))
            error ('butter: filter order n must be a positive integer');
        end
        
        stop = 0;
        digital = 1;
        for i=1:length(varargin)
            switch varargin{i}
                case 's', digital = 0;
                case 'z', digital = 1;
                case { 'high', 'stop' }, stop = 1;
                case { 'low',  'pass' }, stop = 0;
                otherwise,  error ('butter: expected [high|stop] or [s|z]');
            end
        end
        
        
        [r, c]=size(W);
        if (~(length(W)<=2 && (r==1 || c==1)))
            error ('butter: frequency must be given as w0 or [w0, w1]');
        elseif (~(length(W)==1 || length(W) == 2))
            error ('butter: only one filter band allowed');
        elseif (length(W)==2 && ~(W(1) < W(2)))
            error ('butter: first band edge must be smaller than second');
        end
        
        if ( digital && ~all(W >= 0 & W <= 1))
            error ('butter: critical frequencies must be in (0 1)');
        elseif ( ~digital && ~all(W >= 0 ))
            error ('butter: critical frequencies must be in (0 inf)');
        end
        
        % Prewarp to the band edges to s plane
        if digital
            T = 2;       % sampling frequency of 2 Hz
            W = 2/T*tan(pi*W/T);
        end
        
        % Generate splane poles for the prototype butterworth filter
        % source: Kuc
        C = 1; % default cutoff frequency
        pole = C*exp(1i*pi*(2*[1:n] + n - 1)/(2*n));
        if mod(n,2) == 1, pole((n+1)/2) = -1; end  % pure real value at exp(i*pi)
        zero = [];
        gain = C^n;
        
        % splane frequency transform
        [zero, pole, gain] = sftrans(zero, pole, gain, W, stop);
        
        % Use bilinear transform to convert poles to the z plane
        if digital
            [zero, pole, gain] = bilinear(zero, pole, gain, T);
        end
        
        % convert to the correct output form
        if nargout==2,
            a = real(gain*poly(zero));
            b = real(poly(pole));
        elseif nargout==3,
            a = zero;
            b = pole;
            c = gain;
        else
            % output ss results
            [a, b, c, d] = zp2ss (zero, pole, gain);
        end
        
        
        
        function [Zz, Zp, Zg] = bilinear(Sz, Sp, Sg, T)
            
            
            p = length(Sp);
            z = length(Sz);
            if z > p || p==0
                error('bilinear: must have at least as many poles as zeros in s-plane');
            end
            
            % ----------------  -------------------------  ------------------------
            % Bilinear          zero: (2+xT)/(2-xT)        pole: (2+xT)/(2-xT)
            %      2 z-1        pole: -1                   zero: -1
            % S -> - ---        gain: (2-xT)/T             gain: (2-xT)/T
            %      T z+1
            % ----------------  -------------------------  ------------------------
            Zg = real(Sg * prod((2-Sz*T)/T) / prod((2-Sp*T)/T));
            Zp = (2+Sp*T)./(2-Sp*T);
            if isempty(Sz)
                Zz = -ones(size(Zp));
            else
                Zz = [(2+Sz*T)./(2-Sz*T)];
                Zz = postpad(Zz, p, -1);
            end
            
            
            
            function y = postpad (x, l, c, dim)
                
                nd = ndims(x);
                sz = size(x);
                if nargin < 4
                    % Find the first non-singleton dimension
                    dim  = 1;
                    while dim < nd+1 && sz(dim)==1
                        dim = dim + 1;
                    end
                    if dim > nd
                        dim = 1;
                    elseif ~(isscalar(dim) && dim == round(dim)) && dim > 0 && dim< nd+1
                        error('postpad: dim must be an integer and valid dimension');
                    end
                end
                
                if ~isscalar(l) || l<0
                    error ('second argument must be a positive scalar');
                end
                
                if dim > nd
                    sz(nd+1:dim) = 1;
                end
                
                d = sz(dim);
                
                if d >= l
                    idx = cell(1,nd);
                    for i = 1:nd
                        idx{i} = 1:sz(i);
                    end
                    idx{dim} = 1:l;
                    y = x(idx{:});
                else
                    sz(dim) = l-d;
                    y = cat(dim, x, c * ones(sz));
                end
                
            end
            
            
        end
        
        
    end


%%%%%% sftrans
    function [Sz, Sp, Sg] = sftrans(Sz, Sp, Sg, W, stop)
        
        if (nargin ~= 5)
            usage('[Sz, Sp, Sg] = sftrans(Sz, Sp, Sg, W, stop)');
        end;
        
        C = 1;
        p = length(Sp);
        z = length(Sz);
        if z > p || p == 0
            error('sftrans: must have at least as many poles as zeros in s-plane');
        end
        
        if length(W)==2
            Fl = W(1);
            Fh = W(2);
            if stop
                % ----------------  -------------------------  ------------------------
                % Band Stop         zero: b ? sqrt(b^2-FhFl)   pole: b ? sqrt(b^2-FhFl)
                %        S(Fh-Fl)   pole: ?sqrt(-FhFl)         zero: ?sqrt(-FhFl)
                % S -> C --------   gain: -x                   gain: -1/x
                %        S^2+FhFl   b=C/x (Fh-Fl)/2            b=C/x (Fh-Fl)/2
                % ----------------  -------------------------  ------------------------
                if (isempty(Sz))
                    Sg = Sg * real (1./ prod(-Sp));
                elseif (isempty(Sp))
                    Sg = Sg * real(prod(-Sz));
                else
                    Sg = Sg * real(prod(-Sz)/prod(-Sp));
                end
                b = (C*(Fh-Fl)/2)./Sp;
                Sp = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
                extend = [sqrt(-Fh*Fl), -sqrt(-Fh*Fl)];
                if isempty(Sz)
                    Sz = [extend(1+rem([1:2*p],2))];
                else
                    b = (C*(Fh-Fl)/2)./Sz;
                    Sz = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
                    if (p > z)
                        Sz = [Sz, extend(1+rem([1:2*(p-z)],2))];
                    end
                end
            else
                
                % ----------------  -------------------------  ------------------------
                % Band Pass         zero: b ? sqrt(b^2-FhFl)   pole: b ? sqrt(b^2-FhFl)
                %        S^2+FhFl   pole: 0                    zero: 0
                % S -> C --------   gain: C/(Fh-Fl)            gain: (Fh-Fl)/C
                %        S(Fh-Fl)   b=x/C (Fh-Fl)/2            b=x/C (Fh-Fl)/2
                % ----------------  -------------------------  ------------------------
                Sg = Sg * (C/(Fh-Fl))^(z-p);
                b = Sp*((Fh-Fl)/(2*C));
                Sp = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
                if isempty(Sz)
                    Sz = zeros(1,p);
                else
                    b = Sz*((Fh-Fl)/(2*C));
                    Sz = [b+sqrt(b.^2-Fh*Fl), b-sqrt(b.^2-Fh*Fl)];
                    if (p>z)
                        Sz = [Sz, zeros(1, (p-z))];
                    end
                end
            end
        else
            Fc = W;
            if stop
                % ----------------  -------------------------  ------------------------
                % High Pass         zero: Fc C/x               pole: Fc C/x
                % S -> C Fc/S       pole: 0                    zero: 0
                %                   gain: -x                   gain: -1/x
                % ----------------  -------------------------  ------------------------
                if (isempty(Sz))
                    Sg = Sg * real (1./ prod(-Sp));
                elseif (isempty(Sp))
                    Sg = Sg * real(prod(-Sz));
                else
                    Sg = Sg * real(prod(-Sz)/prod(-Sp));
                end
                Sp = C * Fc ./ Sp;
                if isempty(Sz)
                    Sz = zeros(1,p);
                else
                    Sz = [C * Fc ./ Sz];
                    if (p > z)
                        Sz = [Sz, zeros(1,p-z)];
                    end
                end
            else
                % ----------------  -------------------------  ------------------------
                % Low Pass          zero: Fc x/C               pole: Fc x/C
                % S -> C S/Fc       gain: C/Fc                 gain: Fc/C
                % ----------------  -------------------------  ------------------------
                Sg = Sg * (C/Fc)^(z-p);
                Sp = Fc * Sp / C;
                Sz = Fc * Sz / C;
            end
        end
    end


%%%% filtfilt
    function y = filtfilt(b, a, x)
        if (nargin ~= 3)
            usage('y=filtfilt(b,a,x)');
        end
        
        rotate = (size(x, 1)==1);
        
        if rotate	% a row vector
            x = x(:);			% make it a column vector
        end
        
        lx = size(x,1);
        a = a(:).';
        b = b(:).';
        lb = length(b);
        la = length(a);
        n = max(lb, la);
        lrefl = 3 * (n - 1);
        if la < n, a(n) = 0; end
        if lb < n, b(n) = 0; end
        
        % Compute a the initial state taking inspiration from
        % Likhterov & Kopeika, 2003. "Hardware-efficient technique for
        %     minimizing startup transients in Direct Form II digital filters"
        kdc = sum(b) / sum(a);
        if (abs(kdc) < inf) % neither NaN nor +/- Inf
            si = fliplr(cumsum(fliplr(b - kdc * a)));
        else
            si = zeros(size(a)); % fall back to zero initialization
        end
        si(1) = [];
        
        for c = 1:size(x, 2)	% filter all columns, one by one
            v = [2*x(1,c)-x((lrefl+1):-1:2,c); x(:,c);
                2*x(end,c)-x((end-1):-1:end-lrefl,c)]; % a column vector
            
            % Do forward and reverse filtering
            v = filter(b,a,v,si*v(1));		       % forward filter
            v = flipud(filter(b,a,flipud(v),si*v(end))); % reverse filter
            y(:,c) = v((lrefl+1):(lx+lrefl));
        end
        
        if (rotate)			% x was a row vector
            y = rot90(y);		% rotate it back
        end
    end

end



