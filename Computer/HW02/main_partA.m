clc;clear;close all
addpath("../../gspbox"); gsp_start
%% Q01
W = [0 1 1 0 0 0 0 1
     1 0 1 1 1 0 0 1
     1 1 0 1 0 0 0 0 
     0 1 1 0 1 1 0 1
     0 1 0 1 0 1 1 1
     0 0 0 1 1 0 1 0
     0 0 0 0 1 1 0 0
     1 1 0 1 1 0 0 0];
coords = [0 0; 0.2 0.5; 0.1 1; 
     .6 1; 0.7 0.5; 0.9 0.8; 0.8 .35;
     0.5 0];
G = gsp_graph(W, coords);
 
figure("Color", "w")
gsp_plot_graph(G); title("G")
%% Q02
if not(isfield(G, 'U'))
G = gsp_compute_fourier_basis(G);
end
x = 2 * G.U(:,1) + G.U(:,2);
%% Q03
% % snr = 10*log10(sp / sn)
snr = 10;
% sigma2 = 10^-(snr/10) * mean(abs(x).^2);
% n = sqrt(sigma2) * randn(size(x));
% xn = x + n;
% rng(4)
xn = awgn(x, snr, 'measured');
disp(10*log10(sum(xn.^2) / sum((xn-x).^2)))

figure("Color","w")
subplot(121); title("x");
gsp_plot_signal(G,x); 
caxis([-0.1 1.8])
subplot(122); title("x+n");
gsp_plot_signal(G,xn); 
caxis([-0.1 1.8])

figure("Color","w"); hold on
stem(x) 
stem(xn, "--");
legend("x", "x+n");
xlim([0 size(x,1)]+0.5)
% close all
%% 

sigma = svd(G.W);
W_norm = G.W / max(sigma);
[U_norm,e_norm] = eig(W_norm, 'vector');

xhat = G.U'*x;
xnhat = G.U'*xn;

xhat_norm = U_norm'*x;
xnhat_norm = U_norm'*xn;

figure("Color","w")
subplot(2,2,1)
stem(G.e, xhat); ylim([-1 3]); xlim([-1 7])
title("$\hat{x}$ using $L$",'interpreter','latex')
subplot(2,2,3)
stem(G.e, xnhat); ylim([-1 3]); xlim([-1 7])
title("$\hat{x_n}$ using $L$",'interpreter','latex')
subplot(2,2,2); 
stem(e_norm, xhat_norm); ylim([-3 1.5]); xlim([-.75 1.25])
title("$\hat{x}$ using $W_{norm}$",'interpreter','latex')
subplot(2,2,4)
stem(e_norm, xnhat_norm); ylim([-3 1.5]); xlim([-.75 1.25])
title("$\hat{x_n}$ using $W_{norm}$",'interpreter','latex')

%%
h = [1 1 0 0 0 0 0 0].';
h_norm = [0 0 0 0 0 0 1 1].';
figure("Color","w")
subplot(1,2,1)
stem(G.e, h); ylim([-0.5 1.5]); xlim([-1 7])
title("$h_L$ filter",'interpreter','latex')
subplot(1,2,2)
stem(e_norm, h_norm); ylim([-0.5 1.5]); xlim([-.75 1.25])
title("$h_{W_{norm}}$ filter",'interpreter','latex')
%% Q06
xd = G.U*(xnhat.*h);
xd_norm = U_norm*(xnhat_norm.*h_norm);
%%
figure("Color","w"); hold on

figure("Color","w")
subplot(221); title("clean x");
gsp_plot_signal(G,x); 
caxis([-0.2 1.8])

subplot(222); title("noisy x");%, "denoised x using L", "denoisted x using w_{norm}"));
gsp_plot_signal(G,xn); 
caxis([-0.2 1.8])

subplot(223); title("denoised x using L");%, ", "denoisted x using w_{norm}"));
gsp_plot_signal(G,xd); 
caxis([-0.2 1.8])

subplot(224); title("denoisted x using w_{norm}");%, "denoised x using L", ));
gsp_plot_signal(G,xd_norm); 
caxis([-0.2 1.8])
%%
figure("Color","w"); hold on
a = 0.1;
stem((1:8)+0*a, x);
stem((1:8)+1*a, xn);
stem((1:8)+2*a, xd);
stem((1:8)+3*a, xd_norm);
xlim([0.5 8.5])
legend("clean x", "noisy x", "denoised x using L", "denoisted x using w_{norm}")
%% Q07
disp(10*log10(sum(xn.^2) / sum((xn-x).^2)))
disp(10*log10(sum(xd.^2) / sum((xd-x).^2)))
disp(10*log10(sum(xd_norm.^2) / sum((xd_norm-x).^2)))
%% Q08
func_L =@(x) G.U*((G.U'*x).*h);
y = func_L(xn);
Ly = func_L(G.L*xn);
disp(mean(abs(G.L*y - Ly)))
disp(mean(abs(func_L(x + xn)-func_L(x)-func_L(xn))))


func_Wnorm =@(x) U_norm*((U_norm'*x).*h_norm);
y = func_Wnorm(xn);
Ly = func_Wnorm(W_norm*xn);
disp(mean(abs(W_norm*y - Ly)))
disp(mean(abs(func_Wnorm(x + xn)-func_Wnorm(x)-func_Wnorm(xn))))

%% Q09
V = G.e.^(0:2);
h_fir = pinv(V)*h;

V_norm = e_norm.^(0:2);
h_fir_norm = pinv(V_norm)*h_norm;

%% Q10
h_new = V*h_fir;
h_new_norm = V_norm*h_fir_norm;

figure("Color","W");
subplot(211); hold on
stem(G.e, h);
stem(G.e, h_new);
legend("ideal filter", "FIR filter")
title("using L")

subplot(212); hold on
stem(e_norm, h_norm);
stem(e_norm, h_new_norm);
legend("ideal filter", "FIR filter")
title("using W_{norm}")

%% Q11
xd_new = G.U*(xnhat.*h_new);
xd_new_norm = U_norm*(xnhat_norm.*h_new_norm);

figure("Color","w"); hold on
a = 0.1;
stem((1:8)+0*a, x);
stem((1:8)+1*a, xn);
stem((1:8)+2*a, xd);
stem((1:8)+3*a, xd_norm);
stem((1:8)+4*a, xd_new);
stem((1:8)+5*a, xd_new_norm);
xlim([0.5 8.5+5*a])
legend("clean x", "noisy x", "denoised x using L", "denoisted x using w_{norm}", ...
    "new denoised x using L", "new denoisted x using w_{norm}")

disp(10*log10(sum(xn.^2) / sum((xn-x).^2)))
disp(10*log10(sum(xd_new.^2) / sum((xd_new-x).^2)))
disp(10*log10(sum(xd_new_norm.^2) / sum((xd_new_norm-x).^2)))

%%
clc
M = 2;
V = G.e.^(0:(M-1));
h_fir = pinv(V)*h;
V_norm = e_norm.^(0:(M-1));
h_fir_norm = pinv(V_norm)*h_norm;

h_new = V*h_fir;
h_new_norm = V_norm*h_fir_norm;

xd_new = G.U*(xnhat.*h_new);
xd_new_norm = U_norm*(xnhat_norm.*h_new_norm);


figure("Color","w"); hold on
a = 0.1;
stem((1:8)+0*a, x);
stem((1:8)+1*a, xn);
stem((1:8)+2*a, xd);
stem((1:8)+3*a, xd_norm);
stem((1:8)+4*a, xd_new);
stem((1:8)+5*a, xd_new_norm);
xlim([0.5 8.5+5*a])
legend("clean x", "noisy x", "denoised x using L", "denoisted x using w_{norm}", ...
    "new denoised x using L", "new denoisted x using w_{norm}")

disp(10*log10(sum(xn.^2) / sum((xn-x).^2)))
disp(10*log10(sum(xd_new.^2) / sum((xd_new-x).^2)))
disp(10*log10(sum(xd_new_norm.^2) / sum((xd_new_norm-x).^2)))
%% Q13
snr_L = zeros(7,1);
snr_Wnorm = zeros(7,1);
for M=2:8
V = G.e.^(0:(M-1));
h_fir = pinv(V)*h;
V_norm = e_norm.^(0:(M-1));
h_fir_norm = pinv(V_norm)*h_norm;

h_new = V*h_fir;
h_new_norm = V_norm*h_fir_norm;

xd_new = G.U*(xnhat.*h_new);
xd_new_norm = U_norm*(xnhat_norm.*h_new_norm);

snr_L(M-1) = (10*log10(sum(xd_new.^2) / sum((xd_new-x).^2)));
snr_Wnorm(M-1) = (10*log10(sum(xd_new_norm.^2) / sum((xd_new_norm-x).^2)));
end
figure("Color","w"); hold on
plot(2:8, snr_L, "-o");
plot(2:8, snr_Wnorm,"-o");
xlim([1 9])
xlabel("M"); ylabel("SNR(db)")
legend("using L", "using W_{norm}")