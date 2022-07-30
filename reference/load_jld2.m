function [data_2] = load_jld2(data, step)
%load_mat Load jld2 version of data .jld2 at given step
%     path_2 = 'C:\Users\ich\outTest\';
    path_2 = '/Users/z7717/test/';
    name_2 = 'M_1078_';
    root_file_2 = [path_2 name_2 num2str(step) '.jld2'];
    data_2 = h5read(root_file_2, ['/' data]);
end