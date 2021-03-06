
%Valores parametros
    ro=12/100 ; %radio olla (12 cm) OK
    ho=27/100; %altura olla (27 cm) OK
    Ao=2*pi*ro*ho+2*pi*(ro^2); %Area superficial olla OK
    kma=225; %Coeficiente de transferencia aluminio OK
    pa= 997; %densidad agua kg/m3 OK
    Va= pi*(ro^2)*ho; %volumen de agua contenido en la olla (volumen de la olla) OK OK
    ca= 4.175*1000; % calor especifico agua J/(kg*K)OK
    rl= 3.5/100; % radio tarro de leche condensada (3.5 cm) OK
    hl=0.1 ;%altura lata de leche condensada (10 cm)OK
    gl= 0.001; %grosor aluminio de la lata (1mm) OK
    Al=2*pi*rl*hl+2*pi*(rl^2); %area manto y tapas lata (m2)OK
    kla=225; %Coeficiente de transferencia aluminio OK
    pl= 2700; %densidad aluminio (kg/m3)OK
    Vl=2*pi*(rl^2)*gl+hl*2*pi*rl*gl; %volumen metal lata de aluminio OK
    ml=pl*Vl; %masa lata de aluminio OK
    cl=880; % calor especifico aluminio J/(kg*K)OK
    kcl=225; %Coeficiente de transferencia aluminio OK
    pc=1300; %densidad leche condensada (kg/m3)OK
    Vc=pi*(rl^2)*hl; %volumen lata (volumen de leche condensada)OK
    cc=3.984*1000; % calor especifico leche cond. J/(kg*K) OK
    T0=20; %temperatura ambienteOK
    
    
    %Calculo coeficientes de cada expresion
    K0=1/(pa*Va*ca);
    K1=Ao*kma/(pa*Va*ca);
    K2=Al*kla/(pa*Va*ca);
    K3=Al*kla/(ml*cl);
    K4=Al*kcl/(ml*cl);
    K5=Al*kcl/(pc*Vc*cc);
    
    %Matrices Simulacion
    A = [-K1-K2,K2,0; K3, -K3-K4, K4;0,K5,-K5]; %Matriz con constantes
    B = [K0; 0;0];
    C = [1, 1,1];
    E=[K1*T0;0;0];
    
    x0 = [20,20,20];
    u=6000; %Calor de entrada
    %Grafico
    temp=simout;
    t=tout;
    
    plot(t,temp)
    xlabel 'Tiempo [seg]';
    ylabel 'Temperatura [°C]';
    grid;
    title 'Temperatura de las distintas partes del sistema en funcion del tiempo'
    legend ('Temperatura agua T_a','Temperatura lata T_l','Temperatura manjar T_m')
    legend('Location','southeast')
    
 
