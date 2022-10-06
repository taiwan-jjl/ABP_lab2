%%
%
semilogx(array_size_O0, bw_O0, "-x", array_size_cublas, bw_cublas, "-x");
xlabel("# of N");
ylabel("GBytes/sec");
title("Task-2c, cuda matrix-vector product");
L1 = 65536;
L2 = 3145728/4;
%solve x^2 + x == memory
%xline(5.5536);
%xline(3.99976);

yline(264);

legend("O0", "cublas", "GDDR6-192bit=264GB/s");
%xlim([L2/3*0.8 L3/3*1.3]);
%xlim([L1/3*0.8 L2/3*1.3]);
%xlim([0/3*0.5 L1/3*1.3]);
%xlim([0/3*0.5 L1/3*0.1]);