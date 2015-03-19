function [c]=TEautocorr(data,maxlag)

% TEXCORR
%
%
%
%
%
%
%
%
%
%
% Version 1.0 by Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt 2011
%


% size of data
[M,N] = size(data);

% find nearest power
m = 2*M-1;
[F,E] = log2(abs(m));

% Check if n is an exact power of 2.
if ~isempty(F) && F == 0.5
    E = E-1;
end

% check infinity
checkinfinity = ~isfinite(F);
E(checkinfinity) = F(checkinfinity);


% Compute autocorrelation via FFT
fft_data = fft(data,2^E);
c = ifft(abs(fft_data).^2);


% Force real data
c = real(c);

% define lags for moving
lags = -maxlag:maxlag;
if maxlag >= M
	c = [zeros(maxlag-M+1,N^2);c(end-M+2:end,:);c(1:M,:);zeros(maxlag-M+1,N^2)];
else
	c = [c(end-maxlag+1:end,:);c(1:maxlag+1,:)];
end

% normalize data 
c = c./c(maxlag+1);


            
            
            
            
            