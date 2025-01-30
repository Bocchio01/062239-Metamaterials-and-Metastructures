function [measured, x, nodes, infos] = read_data_set(filename)
% This function assumes that x vector of acquisition is equal for every node
% observed. If not, one should use the provided `readuff()` function and
% load each node data separately in a cell of structs.

if (nargin == 0)
    filename = uigetfile({'*.unv'}, 'Select experimental dataset', [pwd '\experimental\data\nonreciprocal']);
    if (filename == 0)
        error('No file chosen. Exit.');
    end
end

[UffDataSets, Info] = readuff(['experimental/data/nonreciprocal/' filename]);

assert(all(Info.errcode == 0), 'Error while loading data');

nodes = UffDataSets(Info.dsTypes == 2411);
nodes = rmfield(nodes{1}, {'nodeN', 'defCS', 'dispCS', 'color', 'dsType', 'binary'});

data = UffDataSets(Info.dsTypes == 58);

x = data{1, 1}.x';

infos = struct( ...
    'filename', filename, ...
    'name', data{1, 1}.d1, ...
    'description', data{1, 1}.d2, ...
    'x_unit', data{1, 1}.abscissaUnitsLabel, ...
    'y_unit_num', data{1, 1}.ordinateNumUnitsLabel, ...
    'x_unit_den', data{1, 1}.ordinateDenumUnitsLabel);

measured = zeros(length(x), length(nodes.x));
for node_idx = 1:length(nodes.x)
    measured(:, node_idx) = data{1, node_idx}.measData;
end

end

