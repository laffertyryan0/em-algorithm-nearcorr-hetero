addpath("../src");

nc_false = load('../artifacts/rho_est_fisherinv_nearcorrfalse.mat');
nc_true = load('../artifacts/rho_est_fisherinv_nearcorrtrue.mat');

falseeig = eig(vecLInverse(nc_false.rho_est_fisherinv{1}));
trueeig = eig(vecLInverse(nc_true.rho_est_fisherinv{1}));

figure,
hold on
histogram(falseeig,'Normalization','pdf')
histogram(trueeig,'Normalization','pdf')
hold off