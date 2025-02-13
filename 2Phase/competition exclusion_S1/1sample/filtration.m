clc
clear
% 读取文件
fileID = fopen('data.txt', 'r');
data = fscanf(fileID, '%f', [3, Inf])'; % 假设每行有3列
fclose(fileID);

% 初始化新数据
newData = [];

% 遍历每一行数据
for i = 1:size(data, 1)
    % 如果行中没有大于1000的数，则保留该行
    if ~any(data(i, :) > 1000)
        newData = [newData; data(i, :)];
    end
end

% 将结果写入新文件
fileID = fopen('filtered_data.txt', 'w');
fprintf(fileID, '%f %f %f\n', newData');
fclose(fileID);
