%%
%
semilogx(array_size_O0, gflop_O0, "-x", array_size_cublas, gflop_cublas, "-x", array_size_cpu, gflop_cpu, "-x");
%loglog(array_size_O0, gflop_O0, "-x", array_size_cublas, gflop_cublas, "-x", array_size_cpu, gflop_cpu, "-x");
xlabel("# of N");
ylabel("GFlop/sec");
title("Task-4, cuda matrix-matrix product");
L1 = 65536;
L2 = 3145728/4;
%solve x^2 + x == memory
xline(147.802);
xline(512);

yline(12.78*1000);

legend("O0", "cublas", "CPU", "reg. # = 65536", "L2=3MB", "FP32 peak perf.=12.78 TFLOPS");
%xlim([L2/3*0.8 L3/3*1.3]);
%xlim([L1/3*0.8 L2/3*1.3]);
%xlim([0/3*0.5 L1/3*1.3]);
%xlim([0/3*0.5 L1/3*0.1]);