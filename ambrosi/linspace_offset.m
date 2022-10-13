function x = linspace_offset(x1,x2,N)


N1 = N + 1;
x1 = linspace(x1,x2,N1);
x = x1(1:end-1) + diff(x1);