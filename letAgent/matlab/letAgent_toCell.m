function c = letAgent_toCell(x)
% letAgent_toCell  jsondecode returns struct arrays or cells; always give a cell.
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
