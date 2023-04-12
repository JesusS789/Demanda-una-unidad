% Programa para hallar el despacho ecónomico considerando perdidas
% Se necesitan como entradas los coeficientes asociados a la funcion de
% costo de los generadores involucrados

    %Variables de entrada necesarias
        % B: Matriz de pérdidas del sistema de potencia
        % Sbase: Potencia base del sistema
        % a: Coeficientes del factor cuadratico del costo en generadores
        % b: Coeficientes del factor lineal del costo en generadores
        % Pc: Potencia real consumida en el intervalo evaluado
        % Error: Precisión deseada en la solución
        % Xs: Punto inicial del metodo iterativo
        % Pgmin: Potencias minimas de los generadores
        % Pgmax: Potencias maximas de los generadores

% Inicialización
a=a*Sbase^2;
b=b*Sbase;
Pc=Pc/Sbase;
Error=1e-3;
Pgmin=Pgmin/Sbase;
Pgmax=Pgmax/Sbase;

DespachoUnificado;      % Primera solución a encontrar
Resul=zeros(1,Ng);
Ng0=Ng;
n=0;

% Proceso para verificar los limites de los generadores
while (any(Pgd>Pgmax) || any(Pgd<Pgmin))
  Xs=double(Xs);
  for i=1:Ng
      if (Xs(i)>Pgmax(i))
        Xs(i)=[];
        Resul(i+n)=Pgmax(i);
        Pc=Pc-Pgmax(i);
        Pgd(i)=[]; Pgmax(i)=[]; Pgmin(i)=[];
        a(i)=[]; b(i)=[];
        B(:,i)=[];
        B(i,:)=[];
        fprintf("El generador %d se ajusto a su máx\n", i+n); 
        DespachoUnificado;
        n=n+1;
        break;
 elseif (Xs(i)<Pgmin(i))
        Xs(i)=[];
        Resul(i+n)=Pgmin(i);
        Pc=Pc-Pgmin(i);
        Pgd(i)=[]; Pgmax(i)=[]; Pgmin(i)=[];
        a(i)=[]; b(i)=[];
        B(:,i)=[];
        B(i,:)=[];
        fprintf("El generador %d se ajusto a su min\n", i+n); 
        DespachoUnificado;
        n=n+1;
        break;
    end
  end
end
Xs=double(Xs);
j=1;
for i=1:Ng0
   if (Resul(i)==0)
     Resul(i)=Xs(j); j=j+1;
   end
end
Potencias=Resul*Sbase
