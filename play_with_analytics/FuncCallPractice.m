funcName = 'addy';
nums = 1:2;
for i = 1:length(nums)
    funstr = strcat(funcName, num2str(nums(i)),'(',num2str(1),')');
    eval(funstr)
end