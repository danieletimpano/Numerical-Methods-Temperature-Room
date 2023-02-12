function [succ, it] = ptofis_new(x0, phi, nmax, toll)
%
% [succ, it] = ptofis(x0, phi, nmax, toll)
% Metodo di punto fisso x = phi(x) 
%
% --------Parametri di ingresso:
% x0      Punto di partenza
% phi     Funzione di punto fisso (definita inline o anonimous)
% nmax    Numero massimo di iterazioni
% toll    Tolleranza sul test d'arresto
%
% --------Parametri di uscita:
% succ   Vett. contenente tutte le iterate calcolate
%         (l'ultima componente e' la soluzione)
% it      Iterazioni effettuate
err   = 1 + toll;
it    = 0;
succ  = x0;
xv    = x0;
while (it < nmax && err > toll)
   xn    = phi(xv);
   % adattato al caso vettoriale
   err   = norm(xn - xv, inf);
   succ = [succ, xn];
   it    = it + 1;
   xv    = xn;
end
%fprintf(' \n Numero di Iterazioni    : %d \n',it);
%fprintf(' Punto fisso calcolato   : %12.13f \n',succ(end));
