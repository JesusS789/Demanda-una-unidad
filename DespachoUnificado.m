% Rutina para hallar despacho economico real
syms Lambda;
Ng=length(a);
dX=zeros(Ng+1,1)+10;
iter=0;
Pg=sym('Pg',[1,Ng]);
X=Pg;
X(Ng+1)=Lambda;
for k=1:Ng
Serie=0;
  for j=1:Ng
    Serie=Serie+Pg(j)*B(j,k);
  end
F(k)=Lambda*(1-2*(Serie+B(k,Ng+1)))-a(k)*Pg(k)-b(k);
end
F(Ng+1)=sum(Pg)-Pc;
F=vpa(transpose(F),10);
J=jacobian(F,X);
J=vpa(J,10) ;

Var_sim=X;
while (abs(max(dX)) >=Error)
iter=iter+1;
Jn=vpa(subs(J,Var_sim,Xs'),10);
Fn=vpa(subs(F,Var_sim,Xs'),10);
dX=vpa(Jn^-1*-Fn,10);
Xs=vpa(Xs+dX,7);
end
%Hace primero el loop para asegurar no caer en un fuera de limites
%momentaneo
for k=1:Ng
  Pgd(k)=double(Xs(k));
end
disp([iter , Xs'])
clearvars F J Jn Fn;
