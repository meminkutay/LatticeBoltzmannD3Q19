function j=NumToStr(n)

if n<10, j=['00' num2str(n)];, end
if n>=10, j=['0' num2str(n)];, end
if n>=100, j=num2str(n);, end

