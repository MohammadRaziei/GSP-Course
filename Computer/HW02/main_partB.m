clc;clear;close all
addpath("../../gspbox"); gsp_start
%% Q1
n = 8;
Wstar = [0 ones(1,n-1); ones(n-1,1) zeros(n-1,n-1)];
coords = [0 0 ;cos(2*pi/(n-1)*(0:n-2)).' sin(2*pi/(n-1)*(0:n-2)).'];
Gstar = gsp_graph(Wstar, coords);
if not(isfield(Gstar, 'U'))
Gstar = gsp_compute_fourier_basis(Gstar);
end
Xstar = [100;zeros(n-1,1)];
figure("Color", "w")
% gsp_plot_graph(Gstar); 
gsp_plot_signal(Gstar,Xstar); 
title("Gـ{star}")
%%
n = 8;
Wloop = zeros(n,n);
for i = 1:n
    Wloop(i,:) = circshift([1 0 1 zeros(1,n-3)],i-2);
end
coords = [cos(2*pi/(n)*(0:n-1)).' sin(2*pi/(n)*(0:n-1)).'];
Gloop = gsp_graph(Wloop, coords);
if not(isfield(Gloop, 'U'))
Gloop = gsp_compute_fourier_basis(Gloop);
end
Xloop = [100;zeros(n-1,1)];
figure("Color", "w")
% gsp_plot_graph(Gstar); 
gsp_plot_signal(Gloop,Xloop); 
title("Gـ{loop}")
%%
n = 8;
Wcomp = ones(n,n)-eye(n,n);
coords = [cos(2*pi/(n)*(0:n-1)).' sin(2*pi/(n)*(0:n-1)).'];
Gcomp = gsp_graph(Wcomp, coords);
if not(isfield(Gcomp, 'U'))
Gcomp = gsp_compute_fourier_basis(Gcomp);
end
Xcomp = [100;zeros(n-1,1)];
figure("Color", "w")
% gsp_plot_graph(Gstar); 
gsp_plot_signal(Gcomp,Xcomp); 
title("Gـ{complete}")
%%
Gstar = calc_wnorm(Gstar);
Gloop = calc_wnorm(Gloop);
Gcomp = calc_wnorm(Gcomp);
%%
Xstar1_1 = shift_vertex(Xstar, 1, Gstar.U);
Xstar1_2 = shift_vertex(Xstar, 2, Gstar.U);
Xstar1_3 = shift_vertex(Xstar, 3, Gstar.U);
Xstar1_4 = shift_vertex(Xstar, 4, Gstar.U);
Xstar1_5 = shift_vertex(Xstar, 5, Gstar.U);
Xstar1_6 = shift_vertex(Xstar, 6, Gstar.U);
Xstar1_7 = shift_vertex(Xstar, 7, Gstar.U);
Xstar1_8 = shift_vertex(Xstar, 8, Gstar.U);

Xloop1_1 = shift_vertex(Xloop, 1, Gloop.U);
Xloop1_2 = shift_vertex(Xloop, 2, Gloop.U);
Xloop1_3 = shift_vertex(Xloop, 3, Gloop.U);
Xloop1_4 = shift_vertex(Xloop, 4, Gloop.U);
Xloop1_5 = shift_vertex(Xloop, 5, Gloop.U);
Xloop1_6 = shift_vertex(Xloop, 6, Gloop.U);
Xloop1_7 = shift_vertex(Xloop, 7, Gloop.U);
Xloop1_8 = shift_vertex(Xloop, 8, Gloop.U);

Xcomp1_1 = shift_vertex(Xcomp, 1, Gcomp.U);
Xcomp1_2 = shift_vertex(Xcomp, 2, Gcomp.U);
Xcomp1_3 = shift_vertex(Xcomp, 3, Gcomp.U);
Xcomp1_4 = shift_vertex(Xcomp, 4, Gcomp.U);
Xcomp1_5 = shift_vertex(Xcomp, 5, Gcomp.U);
Xcomp1_6 = shift_vertex(Xcomp, 6, Gcomp.U);
Xcomp1_7 = shift_vertex(Xcomp, 7, Gcomp.U);
Xcomp1_8 = shift_vertex(Xcomp, 8, Gcomp.U);

Xstar2_1 = shift_vertex(Xstar, 1, Gstar.U_Wnorm);
Xstar2_2 = shift_vertex(Xstar, 2, Gstar.U_Wnorm);
Xstar2_3 = shift_vertex(Xstar, 3, Gstar.U_Wnorm);
Xstar2_4 = shift_vertex(Xstar, 4, Gstar.U_Wnorm);
Xstar2_5 = shift_vertex(Xstar, 5, Gstar.U_Wnorm);
Xstar2_6 = shift_vertex(Xstar, 6, Gstar.U_Wnorm);
Xstar2_7 = shift_vertex(Xstar, 7, Gstar.U_Wnorm);
Xstar2_8 = shift_vertex(Xstar, 8, Gstar.U_Wnorm);

Xloop2_1 = shift_vertex(Xloop, 1, Gloop.U_Wnorm);
Xloop2_2 = shift_vertex(Xloop, 2, Gloop.U_Wnorm);
Xloop2_3 = shift_vertex(Xloop, 3, Gloop.U_Wnorm);
Xloop2_4 = shift_vertex(Xloop, 4, Gloop.U_Wnorm);
Xloop2_5 = shift_vertex(Xloop, 5, Gloop.U_Wnorm);
Xloop2_6 = shift_vertex(Xloop, 6, Gloop.U_Wnorm);
Xloop2_7 = shift_vertex(Xloop, 7, Gloop.U_Wnorm);
Xloop2_8 = shift_vertex(Xloop, 8, Gloop.U_Wnorm);

Xcomp2_1 = shift_vertex(Xcomp, 1, Gcomp.U_Wnorm);
Xcomp2_2 = shift_vertex(Xcomp, 2, Gcomp.U_Wnorm);
Xcomp2_3 = shift_vertex(Xcomp, 3, Gcomp.U_Wnorm);
Xcomp2_4 = shift_vertex(Xcomp, 4, Gcomp.U_Wnorm);
Xcomp2_5 = shift_vertex(Xcomp, 5, Gcomp.U_Wnorm);
Xcomp2_6 = shift_vertex(Xcomp, 6, Gcomp.U_Wnorm);
Xcomp2_7 = shift_vertex(Xcomp, 7, Gcomp.U_Wnorm);
Xcomp2_8 = shift_vertex(Xcomp, 8, Gcomp.U_Wnorm);
%%


    
figure("Color", "w")
sgtitle('shifting "star graph" using "L"') 
subplot(3,3,1);
gsp_plot_signal(Gstar,Xstar); 
title("raw data")
for i=1:8
   subplot(3,3,i+1)
   gsp_plot_signal(Gstar,eval("Xstar1_"+i)); 
   title("shifted data to vertex "+i)
end

figure("Color", "w")
sgtitle('shifting "loop graph" using "L"') 
subplot(3,3,1);
gsp_plot_signal(Gloop,Xloop); 
title("raw data")
for i=1:8
   subplot(3,3,i+1)
   gsp_plot_signal(Gloop,eval("Xloop1_"+i)); 
   title("shifted data to vertex "+i)
end

figure("Color", "w")
sgtitle('shifting "complete graph" using "L"') 
subplot(3,3,1);
gsp_plot_signal(Gcomp,Xcomp); 
title("raw data")
for i=1:8
   subplot(3,3,i+1)
   gsp_plot_signal(Gcomp,eval("Xcomp1_"+i)); 
   title("shifted data to vertex "+i)
end



figure("Color", "w")
sgtitle('shifting "star graph" using "W_{norm}"') 
subplot(3,3,1);
gsp_plot_signal(Gstar,Xstar); 
title("raw data")
for i=1:8
   subplot(3,3,i+1)
   gsp_plot_signal(Gstar,eval("Xstar2_"+i)); 
   title("shifted data to vertex "+i)
end

figure("Color", "w")
sgtitle('shifting "loop graph" using "W_{norm}"') 
subplot(3,3,1);
gsp_plot_signal(Gloop,Xloop); 
title("raw data")
for i=1:8
   subplot(3,3,i+1)
   gsp_plot_signal(Gloop,eval("Xloop2_"+i)); 
   title("shifted data to vertex "+i)
end

figure("Color", "w")
sgtitle('shifting "complete graph" using "W_{norm}"') 
subplot(3,3,1);
gsp_plot_signal(Gcomp,Xcomp); 
title("raw data")
for i=1:8
   subplot(3,3,i+1)
   gsp_plot_signal(Gcomp,eval("Xcomp2_"+i)); 
   title("shifted data to vertex "+i)
end



%%
function xm = shift_vertex(x, m, U)
xm = U*diag(U'*x)*U(m,:)';
end

%%
function G = calc_wnorm(G)
sigma = svd(G.W);
W_norm = G.W / max(sigma);
[U_norm,e_norm] = eig(W_norm, 'vector');
G.U_Wnorm = U_norm;
G.e_Wnorm = e_norm;
G.Wnorm = W_norm;
end