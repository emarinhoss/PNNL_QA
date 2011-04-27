function [output, x, y] = load_data_new(run,run_2,t_step)

temp = hdf5info([run '_' int2str(t_step)]);

name = [temp.GroupHierarchy.Groups(2).Name '/' run_2];

output_temp = hdf5read([run '_' int2str(t_step)],name);

dims = size(output_temp);

name_3 = [temp.GroupHierarchy.Groups(2).Name '/grid/'];

lower = hdf5read([run '_' int2str(t_step)],name_3, 'vsLowerBounds');
upper = hdf5read([run '_' int2str(t_step)],name_3, 'vsUpperBounds');
cells = hdf5read([run '_' int2str(t_step)],name_3, 'vsNumCells');

% lower = hdf5read([run '_' int2str(t_step)],[name_3 'lowerBounds']);
% upper = hdf5read([run '_' int2str(t_step)],[name_3 'upperBounds']);
% cells = hdf5read([run '_' int2str(t_step)],[name_3 'numPhysCells']);

% temp_2 = temp.GroupHierarchy.Groups(2).Datasets(2);
% 
% output_temp = hdf5read(temp_2);
% 
% dims = temp_2.Dims;
% 
% lower = temp.GroupHierarchy.Groups(2).Groups(1).Attributes(3).Value;
% upper = temp.GroupHierarchy.Groups(2).Groups(1).Attributes(4).Value;
% cells = temp.GroupHierarchy.Groups(2).Groups(1).Attributes(5).Value;

dx = (upper(1)-lower(1))/double(cells(1));
dy = (upper(2)-lower(2))/double(cells(2));

x = linspace(lower(1)+dx/2,upper(1)-dx/2,cells(1));
y = linspace(lower(2)+dy/2,upper(2)-dy/2,cells(2));

output = zeros(dims(2),dims(3),dims(1));

for k = 1:dims(2)
    for l = 1:dims(3)
        for m = 1:dims(1)
            output(k,l,m) = output_temp(m,k,l);
        end
    end
end
% 
% % output_data = zeros(row,column,depth);
