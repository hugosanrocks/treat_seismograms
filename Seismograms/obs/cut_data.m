clear all
close all
clc


nsta=30;
stan(1,:)='KMMH01';
stan(2,:)='KMMH02';
stan(3,:)='KMMH03';
stan(4,:)='KMMH06';
stan(5,:)='KMMH09';
%stan(6,:)='KMMH10';
stan(6,:)='KMMH11';
%stan(8,:)='KMMH12';
%stan(9,:)='KMMH13';
stan(7,:)='KMMH14';
%stan(11,:)='KMMH15';
stan(8,:)='KMMH16';
stan(9,:)='KMM003';
stan(10,:)='KMM005';
stan(11,:)='KMM012';
stan(12,:)='KMM013';
stan(13,:)='KMM014';
stan(14,:)='KMM018';
stan(15,:)='MYZ001';
stan(16,:)='TKD001';
stan(17,:)='TMC001';
stan(18,:)='KMMH12';
stan(19,:)='KMMH13';
stan(20,:)='KMMH15';
stan(21,:)='OITH08';
stan(22,:)='FKOH10';
stan(23,:)='SAGH04';
stan(24,:)='FKOH03';
stan(25,:)='SAGH02';
stan(26,:)='MYZH04';
stan(27,:)='MYZH08';
stan(28,:)='NGSH06';
stan(29,:)='KGSH04';
stan(30,:)='KGSH08';


comp(1,1:2)='EW';
comp(2,1:2)='NS';
comp(3,1:2)='UD';

fig=1;
for ista = 1:nsta
for icom = 1:3

  %Name of files
%  file=input('File name: ')
  file = sprintf('%s.%s.out.velo.ascii',stan(ista,:),comp(icom,:))

  %Load data and header info
  data=load(file);

  %Header information
  tini=data(1,1);
  tfin=data(end,1);
  duration = tfin - tini;
  dt = data(2,1) - tini;

  %torig=1*3600+25*60+5.47;

option = 1;
if (option == 1)
if (ista == 1)
  torig = 5112.690;
elseif (ista == 2)
  torig = 5114.392;
elseif (ista == 3)
  torig = 5110.850;
elseif (ista == 4)
  torig = 5111.432;
elseif (ista == 5)
  torig = 5111.566;
elseif (ista == 6)
  torig = 5115.094;
elseif (ista == 7)
  torig = 5108.786;
elseif (ista == 8)
  torig = 5108.049;
elseif (ista == 9)
  torig = 5110.875;
elseif (ista == 10)
  torig = 5109.250;
elseif (ista == 14)
  torig = 5114.579;
elseif (ista == 15)
  torig = 5114.547;

elseif (ista == 11)
  torig = 5111.440;
elseif (ista == 12)
  torig = 5114.325;
elseif (ista == 13)
  torig = 5112.848;


elseif (ista == 18)
  torig = 5116.195;
elseif (ista == 19)
  torig = 5116.173;
elseif (ista == 20)
  torig = 5118.428;
elseif (ista == 21)
  torig = 5118.012;
elseif (ista == 22)
  torig = 5115.892;
elseif (ista == 23)
  torig = 5118.504;
elseif (ista == 24)
  torig = 5121.077;
elseif (ista == 25)
  torig = 5122.448;
elseif (ista == 26)
  torig = 5115.915;
elseif (ista == 27)
  torig = 5121.462;
elseif (ista == 28)
  torig = 5119.965;
elseif (ista == 29)
  torig = 5123.909;
elseif (ista == 30)
  torig = 5127.862;
  torig = torig + 2;
end
end

  cont = 1;
  t=tini;
  while (t <= torig)
   t = t+dt;
   cont = cont + 1;
  end 
  fin = cont + 3999;

if ( (ista <= 15) || (ista > 17) )
  data_cut = data(cont:fin,:);
  data_cut(:,1) = data_cut(:,1) - torig;

  %Resample every 0.1 sec
  data_samp1(:,:) = data_cut(1:end,:);
  data_samp(:,:) = data_samp1(1:end,:);

else

  if (ista == 16) 

   t = 1*3600+25*60;
   cont = 1;
   torig = 5114.55;%5115.550;
   while (t <= torig)
    t = t+0.03;
    cont = cont + 1;
   end
   clear data_samp1 data_samp;
   data_samp1(:,:) = data(cont:end,:);
   data_2 = data_samp1(1:end,:);
   data_2(:,1) = data_2(:,1) - torig;
   data_samp(:,:) = data_2(1:end,:);
  end
  if (ista == 17)
   t = 1*3600+25*60;
   cont = 1;
   torig = 5109.901;
   while (t <= torig)
    t = t+0.03;
    cont = cont + 1;
   end
   clear data_samp1 data_samp;
   data_samp1(:,:) = data(cont:end,:);
   data_2 = data_samp1(1:end,:);
   data_2(:,1) = data_2(:,1) - torig;
   data_samp(:,:) = data_2(1:end,:);
  end


end

%  if (ista > 17)
%   clear data_samp data_cut data_samp1
%   torig = 5105.47;
%   cont = 1;
%   t=tini;
%   while (t <= torig)
%    t = t+dt;
%    cont = cont + 1;
%   end
%   fin = cont + 3999;
%   data_cut = data(cont:fin,:);
%   data_cut(:,1) = data_cut(:,1) - torig;
%   data_samp1(:,:) = data_cut(1:end,:);
%   data_samp(:,:) = data_samp1(1:end,:);
%  end

%  if ( icom == 1 )
%    data_samp(:,2) = data_samp(:,2).*-1;
%   if (ista == 17 )
%     data_samp(:,2) = data_samp(:,2).*-1;
%   end
%   if (ista == 21 )
%     data_samp(:,2) = data_samp(:,2).*-1;
%   end

%  end


% CHANGE SIGN ACCORDING TO AXITRA COORDINATE SYSTEM
  if ( icom == 1 )
    data_samp(:,2) = data_samp(:,2).*-1;
  end


  row1 = data_samp(:,2);
  %Write the file to process the data
  fileout=sprintf('obs_S%03d_C%01d.a',ista,icom)
  save('-ascii',fileout,'row1');%'data_samp');


  %figure(fig)
  %plot(data_cut(:,1),data_cut(:,2)),xlim([data_cut(1,1),data_cut(end,1)])
  %fig=fig+1;
  clear row1 data_samp1 data_samp
end
end
