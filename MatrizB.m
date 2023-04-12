%Programa para hallar matriz de pérdidas de un sistema de potencia
%Se necesitan como inputs, Ybarra, valores de flujo de potencia, matriz de
%combinaciones de la red
%Una limitacion de este código es que requiere Ibase(Ib_g) igual para todo  
%generador y cargas

% Variables de entrada necesarias
    % Ybarra: Matriz de admitancia ya reducida de la red
    % C12: Matriz de conexión entre barras y generadores y cargas
    % I_g: Magnitudes de corriente de generadores (Unidades de corriente 
    % deben ser iguales)
    % I_s: Magnitudes de corrientes de carga (Unidades de corriente
    % iguales)
    % Fp_g: Factores de potencia de generadores
    % Fp_s: Factores de potencia de cargas
    % Alfa_g: Ángulo de voltaje en barras asociadas a generadores (Grados)
    % Alfa_s: Ángulo de voltaje en barras asociadas a cargas (Grados)
    % V: Voltajes en barras asociadas a generadores (%)

Ib_g=83673.95206;             %Corriente base generadores. Constante para 
                              % todas las combinaciones en estudio

% Ajustando unidades y signos
i=sqrt(-1);
I_s=-I_s;
Fp_g=Fp_g/100;
Fp_s=Fp_s/100;
V=V/100;
Alfa_g=deg2rad(Alfa_g);
Alfa_s=deg2rad(Alfa_s);

% Inicialización
Ng=length(I_g);                         %Numero de generadores
Ns=length(I_s);                         %Numero de cargas

% Proceso
for j=1:Ng
  Beta_g(j)=-acos(Fp_g(j))+Alfa_g(j);   %Ángulo de corriente generadores
  CI_g(j)=I_g(j)*exp(Beta_g(j)*i);      %Forma compleja corrientes gens
end

for j=1:Ns
  Beta_s(j)=-acos(Fp_s(j))+Alfa_s(j);   %Ángulo de corriente cargas
  CI_s(j)=I_s(j)*exp(Beta_s(j)*i);      %Forma compleja corrientes cargas
end

R_b=real(inv(Ybarra));
C23=eye(Ng);
for j=1:Ns
  K(j)=CI_s(j)/sum(CI_s);
  C23(Ng+j,Ng+1)=K(j);
end
C34=eye(Ng);
for j=1:(Ng+1)
  C34(Ng+1,j)=-1;
end
C34(Ng+1,Ng+1)=1;
It=(sum(CI_g)+sum(CI_s))/Ib_g;          %Corriente de gens a pu
for j=1:Ng
  m(j)=1*exp(-Beta_g(j)*i)*V(j)*cos(Alfa_g(j)-Beta_g(j));
  C45(j,j)=m(j)^-1;
end
C45(Ng+1,Ng+1)=It;
C_t=C12*C23*C34*C45;
B=real(transpose(conj(C_t))*R_b*C_t);
