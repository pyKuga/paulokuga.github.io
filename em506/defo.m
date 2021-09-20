clear all; close all; clc
#TRANSFORMAÇÃO DE DEFORMAÇÕES

ex = 0; #deformacao x
ey = 0; #deformação y
exy = 5.03e-3; #deformação cisalhante

E = [ex exy/2;
     exy/2 ey] #tensor de deformações. geralmente em mm


angulo = deg2rad(45);
P = [cos(angulo) -sin(angulo); 
     sin(angulo) cos(angulo)]; #matriz de rotação

E_trans = P'*E*P #transformação de deformação

ex_linha = E_trans(1,1)
ey_linha = E_trans(2,2)
exy_linha = E_trans(1,2)+E_trans(2,1)

def_prin = eig(E)
def_cis_max = (max(def_prin)-min(def_prin)) #lembr-se que ele é dividido por 2
def_med = (ex+ey)/2
