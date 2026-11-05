function p = cns_expandPath(p)
% cns_expandPath  Expand a leading ~ to the home folder (MATLAB file functions do not).
if ischar(p) && strncmp(p, '~', 1)
    home = getenv('HOME');
    if isempty(home), home = getenv('USERPROFILE'); end
    p = [home p(2:end)];
end
end
