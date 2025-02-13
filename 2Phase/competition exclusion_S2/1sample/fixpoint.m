% 读取文件
fid = fopen('filtered_data.txt', 'r');
if fid == -1
    error('无法打开文件。');
end

% 逐行读取数据并存储
data = [];
line = fgetl(fid);
while ischar(line)
    line_data = sscanf(line, '%f')'; % 从字符行转换为数值
    if numel(line_data) == 3
        data = [data; line_data]; % 将数据添加到矩阵中
    end
    line = fgetl(fid); % 读取下一行
end
fclose(fid);

% 找出重复的行
[unique_data, ~, idx] = unique(data, 'rows');
counts = accumarray(idx, 1);
duplicates = find(counts > 1);

% 输出重复的行数据
if ~isempty(duplicates)
    disp('重复的行数据：');
    for i = 1:numel(duplicates)
        duplicate_idx = find(idx == duplicates(i));
        disp(data(duplicate_idx, :));
    end
else
    disp('没有重复的行数据。');
end
