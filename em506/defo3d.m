clear all; close all; clc
#TRANSFORMAÇÃO DE DEFORMAÇÕES

ex = -2.435e-3; #deformacao x
ey = -0.972e-3; #deformação y
ez = 2.353e-3;   #deformação z
exy = 0; #deformação cisalhante xy
exz = 0; #cis xz
eyz = 0; #cis yz
E = [ex exy/2 exz/2;
     exy/2 ey eyz/2;
     exz/2 eyz/2 ez] #tensor de deformações. geralmente em mm

     
def_prin = eig(E)
