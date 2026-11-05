function c = cns_toCell(x)
% cns_toCell  jsondecode returns struct arrays, cells or char; always give a 1xN cell.
if isempty(x)
    c = {};
elseif iscell(x)
    c = x(:)';
elseif isstruct(x)
    c = num2cell(x(:)');
else
    c = {x};
end
end
