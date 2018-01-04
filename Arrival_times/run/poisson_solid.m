clear all
close all
clc

syms('rho');
syms('mu');

vp = input('Vp: ');
vs = input('Vs: ');

eqn1 = vp == sqrt(3*(vs^2*rho)/rho)


