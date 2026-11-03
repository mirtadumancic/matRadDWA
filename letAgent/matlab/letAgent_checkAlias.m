function letAgent_checkAlias(alias)
% letAgent_checkAlias  The alias is the only case label that leaves MATLAB.
% It must be a short study code (e.g. CNS_03), never a name, MRN or date.
if ~ischar(alias) || isempty(regexp(alias, '^[A-Za-z][A-Za-z0-9_-]{0,31}$', 'once'))
    error('letAgent:alias', ['Case alias must start with a letter and use only letters, ' ...
        'digits, _ or - (max 32 characters), e.g. CNS_03.']);
end
if ~isempty(regexp(alias, '\d{6,}', 'once')) || ~isempty(regexp(alias, '(19|20)\d{2}[-_]?\d{2}[-_]?\d{2}', 'once'))
    error('letAgent:alias', ['Case alias looks like it contains an ID number or a date. ' ...
        'Use a study code such as CNS_03 instead.']);
end
end
