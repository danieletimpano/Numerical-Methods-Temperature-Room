close all
clear all 
clc

%% Dati del problema
%geometria
H_I= 2.4; %m
H_A= 2.7;
H_T= 4;
L= 8;
H_tr= ((L/2)^2 + H_T^2)^(0.5);

%superfici
S_s= L^2 ; %superficie tetto
S_tr= 4*0.5*H_tr*L ; %superficie esterna piramide
S_Aext= 4*L*H_A ; %superficie esterna livello abitabile
S_Iext= 4*L*H_I ; %superficie esterna sottotetto

%volumi
V_T= (L^2*H_T)/3;
V_A= L^2*H_A;
V_I= L^2*H_I; 

%aria
rho= 1.225 %kg/m^3
c= 1000 %J/Kg°C

%trasmittanze
h_t= 0.3; %W/m^2 K
h_m= 0.35;
h_s= 0.4;

%temperatura esterna
T0= 5;
C1= 5;
C11= 1;
C2= ((2*pi)/(24*60*60));

%sorgente
q= 0.8*10^3 %kW
Te1= @(t) T0 - C1*cos(C2.*t); 
Te2= @(t) T0 - C11*cos(C2.*t); 

%% Definizione del sistema di equazioni differenziali
f= @(t,x) [(h_t*S_tr*(Te1(t) - x(1)) + h_t*S_s*(x(2) - x(1)))/(rho*c*V_T); ... %x1 punto
           (h_t*S_s*(x(1) - x(2)) + h_m*S_Aext*(Te1(t) - x(2)) + h_s*S_s*(x(3) - x(2)))/(rho*c*V_A); ... %x2 punto
           (h_s*S_s*(x(2) - x(3)) + h_m*S_Iext*(Te2(t) - x(3)) + h_s*S_s*(Te2(t) - x(3)))/(rho*c*V_I)];%x2 punto
t0=0;
t1=86400;
x0= [18; 18;18];
h=10; %passo di 10 secondi

       
[t_hav, u_hav] = eulero_avanti(f, t0, t1, x0, h)
[t_hind, u_hind, iter_pf] = eulero_indietro_vett(f, t0, t1, x0, h)

figure(1)
plot(t_hav, u_hav(1,:), 'linewidth', 3)
hold on 
plot(t_hav, u_hav(2,:), 'linewidth', 3)
plot(t_hav, u_hav(3,:), 'linewidth', 3)
grid on
title('Andamento della temperatura in 24h [eulero avanti]')
legend('tetto', 'zona abitabile', 'seminterrato')

figure(2)
plot(t_hind, u_hind(1,:), 'linewidth', 3)
hold on 
plot(t_hind, u_hind(2,:), 'linewidth', 3)
plot(t_hind, u_hind(3,:), 'linewidth', 3)
grid on
title('Andamento della temperatura in 24h [eulero indietro]')
legend('tetto', 'zona abitabile', 'seminterrato')

%% Introduzione di una sorgente di calore
g= @(t,x) [(h_t*S_tr*(Te1(t) - x(1)) + h_t*S_s*(x(2) - x(1)))/(rho*c*V_T); ... %x1 punto
           (h_t*S_s*(x(1) - x(2)) + h_m*S_Aext*(Te1(t) - x(2)) + h_s*S_s*(x(3) - x(2)) +q)/(rho*c*V_A); ... %x2 punto
           (h_s*S_s*(x(2) - x(3)) + h_m*S_Iext*(Te2(t) - x(3)) + h_s*S_s*(Te2(t) - x(3)))/(rho*c*V_I)];%x2 punto
t0=0;
t1=86400;
x0= [18; 18;18];
h=2; %passo di 10 secondi

       
[t_hav_sorg, u_hav_sorg] = eulero_avanti(g, t0, t1, x0, h)
[t_hind_sorg, u_hind_sorg, iter_pf_sorg] = eulero_indietro_vett(g, t0, t1, x0, h)

figure(3)
plot(t_hav_sorg, u_hav_sorg(1,:), 'linewidth', 3)
hold on 
plot(t_hav_sorg, u_hav_sorg(2,:), 'linewidth', 3)
plot(t_hav_sorg, u_hav_sorg(3,:), 'linewidth', 3)
grid on
title('Andamento della temperatura in 24h con termine di generazione [eulero avanti]')
legend('tetto', 'zona abitabile', 'seminterrato')

figure(4)
plot(t_hind_sorg, u_hind_sorg(1,:), 'linewidth', 3)
hold on 
plot(t_hind_sorg, u_hind_sorg(2,:), 'linewidth', 3)
plot(t_hind_sorg, u_hind_sorg(3,:), 'linewidth', 3)
grid on
title('Andamento della temperatura in 24h con termine di generazione [eulero indietro]')
legend('tetto', 'zona abitabile', 'seminterrato')

%% Tentativi di isolamento termico
for h_m= [0.35 : -0.02 : 0.05]
    
    y= @(t,x) [(h_t*S_tr*(Te1(t) - x(1)) + h_t*S_s*(x(2) - x(1)))/(rho*c*V_T); ... %x1 punto
           (h_t*S_s*(x(1) - x(2)) + h_m*S_Aext*(Te1(t) - x(2)) + h_s*S_s*(x(3) - x(2)) +q)/(rho*c*V_A); ... %x2 punto
           (h_s*S_s*(x(2) - x(3)) + h_m*S_Iext*(Te2(t) - x(3)) + h_s*S_s*(Te2(t) - x(3)))/(rho*c*V_I)];%x2 punto
t0=0;
t1=86400;
x0= [18; 18;18];
h=2; %passo di 10 secondi

       
[t_hav_is, u_hav_is] = eulero_avanti(y, t0, t1, x0, h)

if u_hav_is(2,:) >= 17
    break
end

end

figure(5)
plot(t_hav_is, u_hav_is(1,:), 'linewidth', 3)
hold on 
plot(t_hav_is, u_hav_is(2,:), 'linewidth', 3)
plot(t_hav_is, u_hav_is(3,:), 'linewidth', 3)
grid on
title('Andamento della temperatura in 24h con isolamento [eulero avanti]')
legend('tetto', 'zona abitabile', 'seminterrato')


