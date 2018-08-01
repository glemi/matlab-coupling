% function [z, zupper, zlower, f] = deripple(f, z, n)
% 
% Removes the ripples of HBAR spectrum by finding the envelope. Uses a
% modified version of the matlab function 'envelope'.
% 
% usage: 
%           [zcenter, zupper, zlower] = deripple(f, z);
%           
%           % resample to have n equally spaced samples (fnew has length n)
%           [zcenter, zupper, zlower, fnew] = deripple(f, z, n);
%           
% see also: 
%           envelope
%
function [z, zupper, zlower, f] = deripple(f, z, n)
    if ~isreal(z)
        z = abs(z);
    end
    
    shift = 0;
    if min(z) < 0
        minz = min(z);
        maxz = max(z);
        shift = abs(min(z)) + (maxz-minz);
        z = z + shift;
    end
    
    df = mean(diff(f(1:10)));
    [zupper, zlower] = myenvelope(z, round(4e6/df), 'peak');
    inan = isnan(zupper); zupper(inan) = z(inan);
    inan = isnan(zlower); zlower(inan) = z(inan);

    z = sqrt(zupper.*zlower);
    
    if shift ~= 0
        z = z - shift;
        zupper = zupper - shift;
        zlower = zlower - shift;
    end
    
    if nargin >=3 && ~isempty(n)
        N = length(f);
        s = floor(N/n);
        i = 1:s:N;
        m = length(i);
        if m > n
            i = i((m-n+1):end);
        end
            
        f = f(i);
        z = z(i);
        zupper = zupper(i);
        zlower = zlower(i);
    end
end