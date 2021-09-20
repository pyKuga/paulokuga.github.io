clc; clear all; close all
format short
          
titulo = "  z   y     b     h       p "
infor = [   0   0     10    150    1;
            0   80    200   10     1;
            0   -80   200   10    1]
AVISO = "CUIDADO COM AS UNIDADES QUE VOCÊ COLOCOU! SERÃO CONVERTIDAS EM METROS"
unidade = "mm"
infor(:,1:end-1)
unidade = "m"
infor(:,1:end-1) = infor(:,1:end-1)*1e-3
pause
clc

A = infor(:,5).*infor(:,3).*infor(:,4);

zn = (infor(:,1)'*A)/sum(A); #centroide em z
yn = (infor(:,2)'*A)/sum(A); #centroide em y

dz = (infor(:,1)-zn); #z em relacao ao centroide
dy = (infor(:,2)-yn); #y em relacao ao centroide


Izz = (infor(:,5).*infor(:,3).*infor(:,4).^3)/12; #momento de inercia em z
Iyy = (infor(:,5).*infor(:,4).*infor(:,3).^3)/12; #momento de inercia em y

Adz2 = A.*dz.^2;
Ady2 = A.*dy.^2;
format short rat

Izzt = Izz + Ady2; #dy é a distancia pro eixo z
Iyyt = Iyy + Adz2; #dz é a distancia pro eixo y

AVISO = "CUIDADO COM AS UNIDADES QUE VOCÊ COLOCOU"

titulo = "          z         y          b          h    existe?          A         dz         dy        Izz        Iyy      Adz^2      Ady^2       Izzt       Iyyt";
infor = [infor A dz dy Izz Iyy Ady2 Adz2 Izzt Iyyt];
data = num2str(rats(infor));
tabela = [titulo; " "; data ]


format short

IyyV = sum(Iyyt);
IzzV = sum(Izzt) ;

resultados = [" " ;
 "centroide z = " num2str(zn) " " unidade ;
 "centroide y = " num2str(yn) " " unidade; 
 "momento de inercia Izz da viga = " num2str(IzzV) " " unidade "^4";
 "momento de inercia Iyy da viga = " num2str(IyyV) " " unidade "^4" ;]

##pause
##clc
##
##Mz = 5*sin(pi/6)*1500
##My = -5*cos(pi/6)*1500
##
##alpha = rad2deg(atan((IzzV/IyyV)*(My/Mz)))
##
##z = 0
##y = 0
##function res = sigma_xx(Mz, My, IzzV, IyyV, z, y)
##  res = -Mz*y/IzzV + My*z/IyyV;
##endfunction
##tensao_em_xx = sigma_xx(Mz, My, IzzV, IyyV, z, y)
