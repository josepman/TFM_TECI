close all

%% Evolución temporal
close all
figure(1)
% Matrices
subplot(241),pcolor(G.Hpos)
%subplot(5,4,9),pcolor(G.GWpos)
title('Hpos')
xlabel('Tiempo (w)')
ylabel('Hpos.')

subplot(242),pcolor(G.particip)  % == subplot(5,4,10),pcolor(G.Ppos)
title('Participación Local')
xlabel('Tiempo (w)')
ylabel('Particip.')

subplot(243),pcolor(G.node_str)  %==subplot(554),pcolor(G.avg_str) == subplot(5,4,18),pcolor(G.C)   % ==  subplot(5,5,23),pcolor(G.Eloc)
title('Eficiencia Local')
xlabel('Tiempo (w)')
ylabel('E_(loc)')

subplot(244),pcolor(G.GWpos)
title('GWpos')
xlabel('Tiempo (w)')
ylabel('BC')

% Series temporales
subplot(245),plot(G.maximized_mod),hold on %== subplot(551),plot(G.optimized_B)
plot(smooth(smooth(G.maximized_mod)),'r--')
title('Modularidad')
xlabel('Tiempo (w)')
ylabel('Q_(opt)')

subplot(246),plot(G.char_path_len), hold on
plot(smooth(smooth(G.char_path_len)),'r--')
title('Characteristic Path Length')
xlabel('Tiempo (w)')
ylabel('L')

subplot(247),plot(G.global_efficiency), hold on
plot(smooth(smooth(G.global_efficiency)),'r--')
title('Eficiencia Global')
xlabel('Tiempo (w)')
ylabel('E_(glob)')

subplot(248),plot(G.transitv), hold on  % ==  subplot(5,5,22),plot(G.Eglob) = subplot(5,5,11),plot(G.vpos)/450 == G.Eglob
plot(smooth(smooth(G.transitv)),'r--')
title('Transitividad Global')
xlabel('Tiempo (w)')
ylabel('T_glob')

% Nodo-nodo
figure(2)

subplot(241),plot(G.Z)
title('Z rand coeff')
xlabel('Nodos')
ylabel('Z')
axis([1 26 -2.5 2.2])

subplot(245),plot(G.large_eig)  % ==  subplot(5,5,10),plot(G.Sp)
title('Eigenvalue')
xlabel('Nodos')
ylabel('Eign')
axis([1 26 0.12 0.25])
subplot(242),pcolor(G.match_ind{1, 1})
title('Match ind w=1')
xlabel('Nodos')
ylabel('Nodos')

subplot(246),pcolor(G.match_ind{1, 25})
title('Match ind w=25')
xlabel('Nodos')
ylabel('Nodos')

subplot(243),pcolor(G.edgeBC_matrix{1, 1})
title('BC w=1')
xlabel('Nodos')
ylabel('Nodos')

subplot(247),pcolor(G.edgeBC_matrix{1, 25})
title('BC w=25')
xlabel('Nodos')
ylabel('Nodos')

subplot(244),pcolor(G.D{1, 1})
title('D w=1')
xlabel('Nodos')
ylabel('Nodos')

subplot(248),pcolor(G.D{1,25})
title('D w=25')
xlabel('Nodos')
ylabel('Nodos')






