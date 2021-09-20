clear all; close all; clc
#TRANSFORMA��O DE DEFORMA��ES

ex = -2.435e-3; #deformacao x
ey = -0.972e-3; #deforma��o y
ez = 2.353e-3;   #deforma��o z
exy = 0; #deforma��o cisalhante xy
exz = 0; #cis xz
eyz = 0; #cis yz
E = [ex exy/2 exz/2;
     exy/2 ey eyz/2;
     exz/2 eyz/2 ez] #tensor de deforma��es. geralmente em mm

     
def_prin = eig(E)
