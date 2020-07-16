clear
X=[0.164 0.079 0.005;
   0.149 0.254 0.255;
   0.126 0.215 0.252;
   0.509 0.302 0.267;
   0.321 0.426 0.367;
   0.577 0.565 0.291;
   0.667 0.655 0.291;
   0.704 0.741 0.470;
   0.779 0.687 0.646;
   0.886 1.070 0.717];
prom=mean(X);
desv=std(X)*sqrt(9/10); %Se calcula la desviación estandar poblacional (matlab calcula la muestral)
[n,m]=size(X);
X_norm=X; %Datos a normalizar
%Normalizar datos
for i=1:n
    for j=1:m
        X_norm(i,j)=(X_norm(i,j)-prom(j))/desv(j);
    end
end
%Matriz de varianza covarianza con datos normalizados
S=(1/(n-1))*(X_norm')*(X_norm);
%Vectores y valores propios de la matriz var-covar
[VecP,ValP]=eig(S);
[Va_p,in]= sort(diag(ValP),'descend'); %Se ordenan los valores propios de mayor a menor. in es el indice de orden
P=VecP(:,in); %Se ordenan los vectores propios siguiendo el orden de los valores propios de mayor a menor
P=P(:,1:2); %Dos componentes principales (2 primeras columnas de VecP)
Y=X_norm*P; %Matriz de proyeccion de los datos a las dos componentes principales

%Score Plot
figure(1)
hold on
scatter(Y(:,1),Y(:,2))
title 'Score Plot de 2 primeras componentes principales'
xlabel '1° componente principal'
ylabel '2° componente principal'
legend({'Puntos generados por 1° y 2° C.P.'},'location','northwest')
hold off

%%%%Elipse de Control%%%%
%Calcular forma cuadratica sobre las dos C.P
ValP=ValP(in,in); %Ordenar valores propios de mayor a menor
sigmacuad_1=ValP(1,1); %Primera varianza=Primer Valor Propio
sigmacuad_2=ValP(2,2); %Segunda varianza=Segundo Valor Propio
%Se obtienen las distancias cuadradaticas para cada una de las mediciones
%usando la formula de la forma cuadratica
for i=1:1:n
   distancia_cuad(i,1)=X_norm(i,:)*P*(ValP(1:2,1:2))^(-1)*P'*X_norm(i,:)'; %Se obtiene la distancia cuadrada  de cada punto
end

%Graficar curvas de nivel 
t=0:0.001:2*pi; %Intervalo de tiempo para seno y coseno
figure(2)
hold on
for i=1:1:n
    %Parametrizando la elipse para cada par ordenado:
    x_eli_1=(sqrt(distancia_cuad(i)*sigmacuad_1))*cos(t); %Componente en eje x
    y_eli_1=(sqrt(distancia_cuad(i)*sigmacuad_2))*sin(t); %Componente en eje y
    plot(x_eli_1,y_eli_1)
    scatter(Y(:,1),Y(:,2))
end
title 'Curvas de nivel de la distribución de los datos'
xlabel '1° componente principal'
ylabel '2° componente principal'
hold off

%%%%Test de Hotelling%%%%
a=2;
f_1=4.737; %alpha=0.05
f_2=9.547; %alpha=0.01
f_3=3.257; %alpha=0.1

eta_1=(((n-1)^2)*(a/(n-a-1))*f_1)/(n*(1+(a*f_1/(n-a-1))));
eta_2=(((n-1)^2)*(a/(n-a-1))*f_2)/(n*(1+(a*f_2/(n-a-1))));
eta_3=(((n-1)^2)*(a/(n-a-1))*f_3)/(n*(1+(a*f_3/(n-a-1))));

figure(3)
hold on
    %%%%elipse para eta_1 (95%)
x_eta1=(sqrt(eta_1*sigmacuad_1))*cos(t);
y_eta1=(sqrt(eta_1*sigmacuad_2))*sin(t);

    %%%%elipse para eta_2 (99%)
x_eta2=(sqrt(eta_2*sigmacuad_1))*cos(t);
y_eta2=(sqrt(eta_2*sigmacuad_2))*sin(t);

    %%%%elipse para eta_3 (90%)
x_eta3=(sqrt(eta_3*sigmacuad_1))*cos(t);
y_eta3=(sqrt(eta_3*sigmacuad_2))*sin(t);
plot(x_eta3,y_eta3)
plot(x_eta1,y_eta1)
plot(x_eta2,y_eta2)
scatter(Y(:,1),Y(:,2),'+')
legend('Elipse 90%','Elipse 95%','Elipse 99%','Datos')
title 'Test de Hotelling para umbral 90%, 95% y 99%'
xlabel '1° componente principal'
ylabel '2° componente principal'
hold off

