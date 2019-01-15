%% Find id of data1 elements in data2
%   Author:         Xiaochen Qiu from Beihang Univ.
%   Announcement:   Not debug free. Feel free to do any modifications and
%                   use it at will
%-Inputs:
% @data1: datas that might share data at same time as @data2
% @data2: datas that might share data at same time as @data1
% @th: threshold, if difference between two quantity is less than @th then we say they are the same
%-Outputs:
% @Ids1: array that includes all found element id in @data1
% @Ids2: array that includes all id in @data2 of founded element in @data1
function [Ids1, Ids2] = findIds (data1, data2, th)
Ids1 = -1;
Ids2 = -1;
for i = 1:length(data1)
    up = data1(i) + th;
    down = data1(i) - th;
    id = find( down<data2 & data2<up );
    if ~isempty(id)
        Ids1 = [Ids1; i];
        Ids2 = [Ids2; id];
    end
end
Ids1 = Ids1(2:length(Ids1));
Ids2 = Ids2(2:length(Ids2));
end