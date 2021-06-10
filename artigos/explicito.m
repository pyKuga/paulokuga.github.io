%Versao dia 09/06 -21h

clc; clear; close all;

type="inf";       %tipo de parede(Finita="fin", infinita="inf")
inicial="cte";

h = 100; %coeficiente de convecção do ar
k = 401; %coeficiente de condução do cobre [J/s*m*K]
rho = 8960; %massa específica do cobre [kg/m^3] 
cp= 376.56; %calor específico do cobre[J/Kg*K]
alpha = k/(rho*cp); %indice de difusividade termica
Tinf= 300; %temperatura do fluido
tempo = 85000; %quantidade de temp0
nt = 101; %numero de pontos de tempo
dt = tempo/(nt-1) ; %tamanho do passo do tempo
vt=0:dt:tempo; %vetor com os pontos no tempo

L = 100; %comprimento da parede [m] (PARA CASO PAREDE INFINITA, USE UM NUMERO GRANDE ex: 3000)
switch type
  case "fin"
    L = L/2
endswitch
 

ns = 101; %numero de pontos de espaço
dx = L/(ns-1) ; %tamanho do passo do espaço
vs=0:dx:L; %vetor com os pontos no espaço

switch inicial
     case "cte"
         T0 = 500; %temperatura inicial
         Ti = T0*ones(1,ns); %temperatura constante inicial no dominio
     case "sin"
         m = 50;
         T1 = 500;
         T2 = 500 ;
         Ti = sin(m*vs) + T1 + (T2-T1)*(vs/(2*L));
endswitch

switch type
    case "fin"
        Fo = (alpha*dt)/(dx^2); %numero de Fourier
        Bi = h*dx/k; %numero de Biot

        Tp = zeros(1,ns);
        Ttotal = zeros(nt,ns);
        Ttotal(1,:)=Ti;
        tic
        for i=2:nt      %loop temporal
            for j=2:ns-1 
              Tp(j)=Fo*Ti(j+1) + (1-2*Fo)*Ti(j) + Fo*Ti(j-1); %tempratura nos nós internos a cada instante
            end
            Tp(1)=Fo*Ti(2) + (1-Fo)*Ti(1); %tempratura na superfície 1 a cada instante
            %Tp(1)=;
            Tp(ns)=2*Fo*(Ti(ns-1)+Bi*Tinf) + (1-2*Fo - 2*Bi*Fo)*Ti(ns); %tempratura na superfície 2 a cada instante
            %Tp(ns)=;
            Ttotal(i,:)=Tp;
            Ti=Tp;
        end
        tempo_de_execucao = toc
        for t=1:nt
            plot(vs,Ttotal(t,:),'r');
            hold on;
            pause(.01)
            clf
        endfor 
        
        
    case "inf"
        Fo = (alpha*dt)/(dx^2); %numero de Fourier
        Bi = h*dx/k; %numero de Biot

        Tp = zeros(1,ns);
        Ttotal = zeros(nt,ns);
        Ttotal(1,:)=Ti;
        tic
        for i=2:nt      %loop temporal
            for j=2:ns-1 
              Tp(j)=Fo*Ti(j+1) + (1-2*Fo)*Ti(j) + Fo*Ti(j-1); %tempratura nos nós internos a cada instante
            endfor
            Tp(1)=2*Fo*(Ti(2)+Bi*Tinf) + (1-2*Fo - 2*Bi*Fo)*Ti(1); %tempratura na superfície 1 a cada instante
            %Tp(1)=;
            %Tp(ns)=2*Fo*(Ti(ns-1)+Bi*Tinf) + (1-2*Fo - 2*Bi*Fo)*Ti(ns); %tempratura na superfície 2 a cada instante
            Tp(ns)=Ttotal(1,end);
            Ttotal(i,:)=Tp;
            Ti=Tp;
        endfor
        tempo_de_execucao = toc

        for t=1:nt
            plot(vs,Ttotal(t,:),'r');
            hold on;
            pause(.01)
            clf
        endfor 
        %soluçao analítica para o caso infinito
        Ta_total=zeros(nt,ns);
        Theta=zeros(1,ns);
        for j=1:nt
            for i=1:ns
                Theta(i)=erfc(vs(i)/(2*((alpha*vt(j))^0.5)))-(exp(h*vs(i)/k + (h^2)*alpha*vt(j)/(k^2)))*(erfc(vs(i)/(2*((alpha*vt(j))^0.5)) + h*((alpha*vt(j))^0.5)/k));
            end
            for s=1:ns
            Ta(s)=Theta(s)*(Tinf-Ttotal(1,s)) + Ttotal(1,s);
            end
            Ta_total(j,:)=Ta;
        endfor
        Ta_total(1,1) = Ti(1); 
        erro = norm(Ttotal-Ta_total)/norm(Ta_total)
      
endswitch    
      
switch type
  case "fin"
         function ksi =ksisol(Bi,n)
          if Bi <= 4
            f = @(x) cosh(pi*x/2) -1 - Bi;
            ini = fsolve(f,1);
          else 
            ini = 1.3;
          endif
          x = [ini (1:1:n)*pi]';
          f = @(x) x*tan(x) - Bi;
          ksi = []
          for i=1:rows(x)
            ksi(i,1) = fsolve(f,x(i));
          endfor
        endfunction

        x_p=vs/L;
        Bi2=h*L/(k);

        Fo2 = alpha*vt/L^2;
        n = 10; #tem uma hora q começa a ficar sem sentido colocar mais do que ísso
        ksi = ksisol(Bi2,n);
        theta = zeros(nt,ns);
        for i=1:n
          z = ksi(i); 
          Cn=(4*sin(z))/(2*z + sin(2*z));
          theta = theta + Cn*exp((-z**2)*Fo2')*cos(z*x_p);
        endfor
        Ta_total =Tinf + (T0-Tinf)*theta;
        erro = norm(Ttotal-Ta_total)/norm(Ta_total) #note que o erro é aceitável (e ainda nem é a solução exata, é apenas a aproximada)
               
        
endswitch

plot(vs,Ttotal(end,:),'r',vs,Ta_total(end,:),'b')
legend('Numérico','Analítico')
