function [num] = getSymNum(s)
%getSymNum You pass it a symbolic variable, it gives you the number
%   Detailed explanation goes here
    symChar = char(s);
    numChar = symChar(2:end);
    num = str2num(numChar);
end

