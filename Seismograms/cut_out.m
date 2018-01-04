clear all
close all
clc

obs_ew=load('obs_velo_ew');
obs_ns=load('obs_velo_ns');
obs_zz=load('obs_velo_zz');

len = length(obs_ew);

nsta = 30;

samps = len/nsta;

for i=1:nsta
 tsampi(i) =  1 + (i-1)*samps;
 tsampf(i) = samps + (i-1)*samps;
end


dt=0.1;
t=0:dt:dt*(samps-1);
nsta=30

for k=1:nsta
%   e=sprintf('dat/obs_S%03d_C1.a',k);
%   n=sprintf('dat/obs_S%03d_C2.a',k);
%   v=sprintf('dat/obs_S%03d_C3.a',k);
%   eobs2=load(e);
%   nobs2=load(n);
%   vobs2=load(v);
    eobs2(:,2)=obs_ew(tsampi(k):tsampf(k));
    nobs2(:,2)=obs_ns(tsampi(k):tsampf(k));
    vobs2(:,2)=obs_zz(tsampi(k):tsampf(k));
    eobs2(:,1)=t';
    nobs2(:,1)=t';
    vobs2(:,1)=t';

   eobs=eobs2(:,2);
   nobs=nobs2(:,2);
   vobs=vobs2(:,2);

  e=sprintf('obs_S%03d_C1',k);
  n=sprintf('obs_S%03d_C2',k);
  v=sprintf('obs_S%03d_C3',k);
  save('-ascii',e,'eobs');
  save('-ascii',n,'nobs');
  save('-ascii',v,'vobs');

  e=sprintf('obs_S%03d_C1.a',k);
  n=sprintf('obs_S%03d_C2.a',k);
  v=sprintf('obs_S%03d_C3.a',k);
  save('-ascii',e,'eobs2');
  save('-ascii',n,'nobs2');
  save('-ascii',v,'vobs2');

end


