function [correct, eps_total, f] = check(grid1,grid2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
assert(all(size(grid1)==size(grid2)));
[ysize, xsize] = size(grid1);
maximum = max(max(max(grid1)), max(max(grid2)));
minimum = min(min(min(grid1)), min(min(grid2)));
max_min = [maximum minimum];
eps_at_max_min = [eps(maximum) eps(minimum)];
[eps_total, idx_eps] = max(eps_at_max_min);
assert(all(size(grid1) == size(grid2)));
atol = 1e-4;
correct = sum(sum(bsxfun(@(x,y) isapprox(x,y,atol),grid1,grid2))) == numel(grid1);
L1_errors = abs(grid1-grid2);
avg_L1_error = sum(sum(L1_errors))/numel(grid1);
med_L1_error = median(L1_errors, 'all');
rel_errors = abs((grid1-grid2)./grid1);
rel_errors(isinf(rel_errors)) = NaN;
rel_errors(rel_errors>100.0) = NaN;
[max_L1_error, I_L1] = max(L1_errors(:));
[I_row_L1, I_col_L1] = ind2sub(size(L1_errors), I_L1);
[max_rel_error, I_rel] = max(rel_errors(:));
[I_row_rel, I_col_rel] = ind2sub(size(rel_errors), I_rel);
disp(strcat("eps(", num2str(max_min(idx_eps)),") = ", num2str(eps_total)));
disp(strcat("median L1 error: ", num2str(med_L1_error)));
disp(strcat("mean L1 error: ", num2str(avg_L1_error)));
disp(strcat("max L1 error: ", num2str(max_L1_error), " at: ", num2str([I_row_L1, I_col_L1])));
disp([grid1(I_row_L1, I_col_L1) grid2(I_row_L1, I_col_L1)]);
disp(strcat("max rel error [%]: ", num2str(max_rel_error*100), " at: ", num2str([I_row_rel, I_col_rel])));
disp([grid1(I_row_rel, I_col_rel) grid2(I_row_rel, I_col_rel)]);
% ymin = max(1, I_row_L1-32);
% ymax = min(ysize, I_row_L1+32);
% xmin = max(1, I_col_L1-32);
% xmax = min(xsize, I_col_L1+32);
ymin = 1;
ymax = ysize;
xmin = 1;
xmax = xsize;
hm = [];
f = figure;
if (1<ysize && 1<xsize)
    f.Position =[10 300 1800 800];
    p1=subplot(1,2,1);
else
    f.Position =[10 300 895 800];
    p1=subplot(1,1,1);
end
hs=histogram(L1_errors);
axis('square');
xlabel('L1 error');
ylabel('node count');
title('L1 error histogram');
if (1<ysize && 1<xsize)
    p2=subplot(1,2,2);
    hm=heatmap(abs(grid1(ymin:ymax,xmin:xmax)-grid2(ymin:ymax,xmin:xmax)),...
        'XData', xmin:xmax, 'YData', ymin:ymax);
    i_x=false(size(hm.XDisplayLabels)); i_x(1)=true; i_x(end)=true;
    i_y=false(size(hm.YDisplayLabels)); i_y(1)=true; i_x(end)=true;
    hm.XDisplayLabels(~i_x) = {''};
    hm.YDisplayLabels(~i_x) = {''};
    xlabel('x [nodes]');
    ylabel('y [nodes]');
    title('L1 error heatmap');
end
end