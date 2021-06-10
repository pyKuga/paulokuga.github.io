clc; close all; clear all
L = 10; #[m] tamanho total da parede (PARA CASO PAREDE INFINITA USE UM NUMERO GRANDE ex: 3000)
ns = 1001; #numero de pontos de espaço

tempo = 100; #[s] quantidade de tempo
nt = 3e3+1; #numero de pontos de tempo


h = 100; %coeficiente de convecção do ar [J/(s*K*m^2)]
k = 401; %coeficiente de condução do cobre [J/(s*m*K)]
rho = 8960; %massa especí­fica do cobre [kg/m^3] 
cp= 376.56; %calor específico do cobre[J/(Kg*K)]

alpha = k/(rho*cp); %indice de difusividade termica
dt = tempo/(nt-1) ; #tamanho do passo do tempo
dx = L/(ns-1) ; #tamanho do passo do espaço
Fo = alpha*dt/dx**2 #numero de fourier
Bil = h*dx/k; #biot local (serve p/ reduzir o numero de termos dentro da convecção)

T = 500*ones(ns,1);
Tinf = 300;
Tminf = [Bil*Fo*Tinf; zeros(ns-2,1); Bil*Fo*Tinf];

K = zeros(ns,ns);

#convec

K(1,1:2) = [(1-Fo*(1+Bil)) Fo];
K(end,end-1:end) = [Fo (1-Fo*(1+Bil))];

savem = [T';];
k1 = Fo; 
k2 = (1-2*Fo); 
k3 = Fo;

for i=2:ns-1
	K(i, i-1) = k1;
	K(i, i) = k2;
	K(i, i+1) = k3;
end

tic;
savem(1,:) = (K*T+Tminf);
for i=1:nt-1
  savem(i+1,:) = (K*savem(i,:)'+Tminf);
endfor
tempo_matricial= toc

tic;
Tp = zeros(1,ns);
Ttotal = zeros(nt,ns);
Ttotal(1,:)=T;

for i=2:nt      %loop temporal
  for j=2:ns-1 
    Tp(j)=Fo*T(j+1) + (1-2*Fo)*T(j) + Fo*T(j-1); %tempratura nos nós internos a cada instante
  end
  Tp(1)=2*Fo*(T(2)+Bil*Tinf) + (1-2*Fo - 2*Bil*Fo)*T(1); %tempratura na superfície 1 a cada instante
  %Tp(1)=;
  %Tp(ns)=2*Fo*(Ti(ns-1)+Bi*Tinf) + (1-2*Fo - 2*Bi*Fo)*Ti(ns); %tempratura na superfície 2 a cada instante
  Tp(ns)=Ttotal(1,end);
  Ttotal(i,:)=Tp;
  T=Tp;
end
tempo_array_edit = toc