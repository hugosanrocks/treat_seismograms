close all
clear all

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT DATA
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ntot = 651;%751;%1024;             % Number of samples per output seismogram
dtot = 0.1;%0.1;%0.064;            % Output sampling rate for output seismograms (s) 

bpf=3;                   % Apply filter (0=no; 1=lowpass ; 2=highpass ; 3=bandpass)
highcut=0.25;            % High cutoff frequency (Hz)
%highcut=1.0;            % High cutoff frequency (Hz)
%highcut=0.2;            % High cutoff frequency (Hz)
%lowcut=0.0125;             % Low cutoff frequency (Hz)
%lowcut=0.0166;             % Low cutoff frequency (Hz)
lowcut=0.025;             % Low cutoff frequency (Hz)

ifield=1;                % Output field (1=velo ; 2=disp)
wdata=1;                 % Write output files (1=yes ; 0=no)

% Open output files
fid1 = fopen('screen.out', 'wt');

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Read arrival times
tempo=importdata('arrival_times.dat');
statn=tempo.textdata;
arrtime=tempo.data;
clear tempo

% Number of etations
nstat=length(statn);

if wdata == 1
    delete('obs_velo_ns');
    delete('obs_velo_ew');
    delete('obs_velo_zz');
end

% Read and process data
for istat=1:nstat
%for istat=1:1

    % +++++++++++++++++++++
    % North-south component
    % +++++++++++++++++++++
    
    fprintf(1,'%s\n','')
    fprintf(1,'%s\n','             North-South Components...')
    
%     fns=[char(statn(istat)) '.an'];
    fns=[char(statn(istat)) '.vn'];
    tmp=importdata(fns);
    accnsori=tmp.data;
    tmpc=char(tmp.textdata);    
    
    % Original time increment (s)
    dt(istat)=str2num(tmpc(14:22));  
    
    % Number of samples, window length and time serie
    ntori=length(accnsori);
    tlori=(ntori-1)*dt(istat);
    timeori=(0:dt(istat):tlori);
    timeint=(0:dtot:tlori);
    tl=(ntot-1)*dtot;
    time=(0:dtot:tl);
    ntint=length(timeint);
    
    % Detrend and remove base line
%     accnsori=detrend(accnsori,'constant');
%     nsp(1)=nearest((arrtime(istat,2)-arrtime(istat,1))/dtot);
%     nsp(2)=nearest((arrtime(istat,2)+6-arrtime(istat,1))/dtot);
%     accnsori=detrend(accnsori,'linear',nsp);
%     accnsori=detrend(accnsori,'linear');
%     % Taper window
%     accnsori=accnsori.*tukeywin(ntori,0.2);
    
    % Interpolate seismogram (dtot)
    accnsint=interp1(timeori,accnsori,timeint);
    
    % Taper window
    accnsint=accnsint'.*tukeywin(ntint,0.1);
    
    % Add zeros to match the theoretical P-wave arrival time
    accns(1:ntot)=0.0;
    nzero=nearest(arrtime(istat,1)/dtot);    
    if nzero+ntint < ntot
        accns(nzero+1:nzero+ntint)=accnsint(1:ntint);
    else
        accns(nzero+1:ntot)=accnsint(1:ntot-nzero);
    end

    % Taper window
%     accns=accns'.*tukeywin(ntot,0.1);
    % Detrend and remove base line
%     accns=detrend(accns,'constant');
%     accns=detrend(accns,'linear');
    % Taper window
%     accns=accns.*tukeywin(ntori,0.05);
     
    % Compute velocity (integrate acceleration)
    velnsold = 0;
    for it=1:ntot               
        % Compute displacements (integrate velocities)
        if it > 1
            velns(it) = accns(it) * dtot + velnsold;
        else
            velns(it) = 0.0;
        end
        velnsold = velns(it);
    end
    
    % Write output seismograms
    if wdata == 1
        fprintf(1,'%s\n','             Write sliprate ascii files...')
        if ifield == 1
            fileoutns=fopen('obs_velo_ns','a');
            fprintf(fileoutns,'%14.6e\n',accns);
        else
            fileoutns=fopen('obs_disp_ns','a');
            fprintf(fileoutns,'%14.6e\n',velns);
        end
    end
    
    % Filter observed seismograms
    if bpf == 1
        fprintf(1,'              > Lowpass filter seismograms (f < %4.2f Hz) \n',highcut)
        fprintf(fid1,'              > Lowpass filter seismograms (f < %4.2f Hz) \n',highcut);
    end
    if bpf == 2
        fprintf(1,'              > Highpass filter seismograms (f > %4.2f Hz) \n',lowcut')
        fprintf(fid1,'              > Highpass filter seismograms (f > %4.2f Hz) \n',lowcut');
    end
    if bpf == 3
        fprintf(1,'              > Bandpass filter seismograms (%4.2f < f < %4.2f Hz) \n',lowcut,highcut')
        fprintf(fid1,'              > Bandpass filter seismograms (%4.2f < f < %4.2f Hz) \n',lowcut,highcut');
    end
    
    if ifield == 1
        field=accns;
    else
        field=velns;
    end
    
    if bpf ~= 0
        
        % Compute lowpass filter
        ratbobs=highcut/(1/(2*dtot));
        [Bfd,Afd]=butter(4,ratbobs,'low');
        Lowobs=impz(Bfd,Afd);

        % Compute highpass filter
        ratbobs=lowcut/(1/(2*dtot));
        [Bfd,Afd]=butter(4,ratbobs,'high');
        Highobs=impz(Bfd,Afd);
        %
        if bpf == 1  % Lowpass filter
            % Filter NS component
            tempo=conv(field,Lowobs);
            field(1:ntot)=tempo(1:ntot);
        end
        if bpf == 2  % Highpass filter
            % Filter NS component
            tempo=conv(field,Highobs);
            field(1:ntot)=tempo(1:ntot);
        end
        if bpf == 3  % Bandpass filter
            % Filter NS component
            tempo1=conv(field,Lowobs);
            tempo=conv(tempo1,Highobs);
            field(1:ntot)=tempo(1:ntot);
        end

        if ifield == 1
            accns=field;
        else
            velns=field;
        end
    end
    clear tempo
    
    % P and S waves arrival time ticks
    natp=nearest(arrtime(istat,1)/dtot)+1;
    nats=nearest(arrtime(istat,2)/dtot)+1;
    
    if ifield == 1
        vmax=max(accns);
        vmin=min(accns);
        absmax=max(abs(accns));
        refp=accns(natp);
        refs=accns(nats);
    else
        vmax=max(velns);
        vmin=min(velns);
        absmax=max(abs(velns));
        refp=velns(natp);
        refs=velns(nats);
    end
    
    ltick=0.15*absmax;
    tickpx=[(natp-1)*dtot (natp-1)*dtot];
    tickpy=[refp+ltick*0.5 refp-ltick*0.5];
    ticksx=[(nats-1)*dtot (nats-1)*dtot];
    ticksy=[refs+ltick*0.5 refs-ltick*0.5];

    h=figure(istat);
    %     plot(timeint,accnsint,'r'); hold on
    %     plot(timeori,accnsori,'g'); hold on
    if ifield == 1
        plot(time,accns,'b'); hold on
        xlabel('Time (s)');
        ylabel('Velocity (cm/s)');
    else
        plot(time,velns,'b'); hold on
        xlabel('Time (s)');
        ylabel('Displacement (cm)');
    end
    plot(tickpx,tickpy,'r'); hold on
    text(tickpx(1),tickpy(1)+0.1*ltick,'P','FontSize',12); hold on
    plot(ticksx,ticksy,'r'); hold on
    text(ticksx(1),ticksy(1)+0.1*ltick,'S','FontSize',12); hold on
    text(tl*0.8,0.8*vmax,char(statn(istat)),'FontSize',12); hold on
    grid on
    axis([0 tl vmin-0.1*absmax vmax+0.1*absmax]);
    title('North-South Component');
    
    if ifield == 1
        fname=['./figure/' char(statn(istat)) '.vns.ps'];
    else
        fname=['./figure/' char(statn(istat)) '.dns.ps'];
    end
%    print(h,'-depsc2',fname)
    
%     figure(100)
%     shift=30;
%     if ifield == 1
%         plot(time,accns+shift*(istat-1),'b'); hold on
%         plot(tickpx,tickpy+shift*(istat-1),'r'); hold on
%         text(tickpx(1),tickpy(1)+0.1*ltick+shift*(istat-1),'P','FontSize',12); hold on
%         plot(ticksx,ticksy+shift*(istat-1),'r'); hold on
%         text(ticksx(1),ticksy(1)+0.1*ltick+shift*(istat-1),'S','FontSize',12); hold on
%     else
%         plot(time,velns,'b'); hold on
%     end
%     text(tl*0.8,5+shift*(istat-1),char(statn(istat)),'FontSize',12); hold on
%     grid on
%     xlim([0 tl]);
%     xlabel('Time (s)');
%     ylabel('Velocity (cm/s)');
%     title('North-South Components');

    clear accnsori tmp tmpc accns

    % +++++++++++++++++++++
    % East-west component
    % +++++++++++++++++++++
    
    fprintf(1,'%s\n','')
    fprintf(1,'%s\n','             East-west Components...')
    
%     fns=[char(statn(istat)) '.ae'];
    few=[char(statn(istat)) '.ve'];
    tmp=importdata(few);
    accewori=tmp.data;
    tmpc=char(tmp.textdata);    
    
    % Original time increment (s)
    dt(istat)=str2num(tmpc(14:22));  
    
    % Number of samples, window length and time serie
    ntori=length(accewori);
    tlori=(ntori-1)*dt(istat);
    timeori=(0:dt(istat):tlori);
    timeint=(0:dtot:tlori);
    tl=(ntot-1)*dtot;
    time=(0:dtot:tl);
    ntint=length(timeint);
    
    % Detrend and remove base line
%     accewori=detrend(accewori,'constant');
%     ewp(1)=nearest((arrtime(istat,2)-arrtime(istat,1))/dtot);
%     ewp(2)=nearest((arrtime(istat,2)+6-arrtime(istat,1))/dtot);
%     accewori=detrend(accewori,'linear',ewp);
    % Taper window
%     accewori=accewori.*tukeywin(ntori,0.05);
%     accewori=detrend(accewori,'linear');
%     % Taper window
%     accewori=accewori.*tukeywin(ntori,0.2);
    
    % Interpolate seismogram (dtot)
    accewint=interp1(timeori,accewori,timeint);
    
    % Taper window
    accewint=accewint'.*tukeywin(ntint,0.1);
    
    % Add zeros to match the theoretical P-wave arrival time
    accew(1:ntot)=0.0;
    nzero=nearest(arrtime(istat,1)/dtot);    
    if nzero+ntint < ntot
        accew(nzero+1:nzero+ntint)=accewint(1:ntint);
    else
        accew(nzero+1:ntot)=accewint(1:ntot-nzero);
    end

    % Taper window
%     accew=accew'.*tukeywin(ntot,0.1);
    
    % Compute velocity (integrate acceleration)
    velewold = 0;
    for it=1:ntot               
        % Compute displacements (integrate velocities)
        if it > 1
            velew(it) = accew(it) * dtot + velewold;
        else
            velew(it) = 0.0;
        end
        velewold = velew(it);
    end
    
    % Write output seismograms
    if wdata == 1
        fprintf(1,'%s\n','             Write sliprate ascii files...')
        if ifield == 1
            fileoutns=fopen('obs_velo_ew','a');
            fprintf(fileoutns,'%14.6e\n',accew);
        else
            fileoutns=fopen('obs_disp_ew','a');
            fprintf(fileoutns,'%14.6e\n',velew);
        end
    end
    
    % Filter observed seismograms
    if bpf == 1
        fprintf(1,'              > Lowpass filter seismograms (f < %4.2f Hz) \n',highcut)
        fprintf(fid1,'              > Lowpass filter seismograms (f < %4.2f Hz) \n',highcut);
    end
    if bpf == 2
        fprintf(1,'              > Highpass filter seismograms (f > %4.2f Hz) \n',lowcut')
        fprintf(fid1,'              > Highpass filter seismograms (f > %4.2f Hz) \n',lowcut');
    end
    if bpf == 3
        fprintf(1,'              > Bandpass filter seismograms (%4.2f < f < %4.2f Hz) \n',lowcut,highcut')
        fprintf(fid1,'              > Bandpass filter seismograms (%4.2f < f < %4.2f Hz) \n',lowcut,highcut');
    end
    
    if ifield == 1
        field=accew;
    else
        field=velew;
    end
    
    if bpf ~= 0
        
        % Compute lowpass filter
        ratbobs=highcut/(1/(2*dtot));
        [Bfd,Afd]=butter(4,ratbobs,'low');
        Lowobs=impz(Bfd,Afd);

        % Compute highpass filter
        ratbobs=lowcut/(1/(2*dtot));
        [Bfd,Afd]=butter(4,ratbobs,'high');
        Highobs=impz(Bfd,Afd);
        %
        if bpf == 1  % Lowpass filter
            % Filter NS component
            tempo=conv(field,Lowobs);
            field(1:ntot)=tempo(1:ntot);
        end
        if bpf == 2  % Highpass filter
            % Filter NS component
            tempo=conv(field,Highobs);
            field(1:ntot)=tempo(1:ntot);
        end
        if bpf == 3  % Bandpass filter
            % Filter NS component
            tempo1=conv(field,Lowobs);
            tempo=conv(tempo1,Highobs);
            field(1:ntot)=tempo(1:ntot);
        end

        if ifield == 1
            accew=field;
        else
            velew=field;
        end
    end
    clear tempo
    
    % P and S waves arrival time ticks
    natp=nearest(arrtime(istat,1)/dtot)+1;
    nats=nearest(arrtime(istat,2)/dtot)+1;
    
    if ifield == 1
        vmax=max(accew);
        vmin=min(accew);
        absmax=max(abs(accew));
        refp=accew(natp);
        refs=accew(nats);
    else
        vmax=max(velew);
        vmin=min(velew);
        absmax=max(abs(velew));
        refp=velew(natp);
        refs=velew(nats);
    end
    
    ltick=0.15*absmax;
    tickpx=[(natp-1)*dtot (natp-1)*dtot];
    tickpy=[refp+ltick*0.5 refp-ltick*0.5];
    ticksx=[(nats-1)*dtot (nats-1)*dtot];
    ticksy=[refs+ltick*0.5 refs-ltick*0.5];

    h=figure(istat+10);
    %     plot(timeint,accewint,'r'); hold on
    %     plot(timeori,accewori,'g'); hold on
    if ifield == 1
        plot(time,accew,'b'); hold on
        xlabel('Time (s)');
        ylabel('Velocity (cm/s)');
    else
        plot(time,velew,'b'); hold on
        xlabel('Time (s)');
        ylabel('Displacement (cm)');
    end
    plot(tickpx,tickpy,'r'); hold on
    text(tickpx(1),tickpy(1)+0.1*ltick,'P','FontSize',12); hold on
    plot(ticksx,ticksy,'r'); hold on
    text(ticksx(1),ticksy(1)+0.1*ltick,'S','FontSize',12); hold on
    text(tl*0.8,0.8*vmax,char(statn(istat)),'FontSize',12); hold on
    grid on
    axis([0 tl vmin-0.1*absmax vmax+0.1*absmax]);
    title('East-West Component');
    
    if ifield == 1
        fname=['./figure/' char(statn(istat)) '.vew.ps'];
    else
        fname=['./figure/' char(statn(istat)) '.dew.ps'];
    end
%    print(h,'-depsc2',fname)

    clear accewori tmp tmpc accew
    
    % +++++++++++++++++++++
    % Vertical component
    % +++++++++++++++++++++
    
    fprintf(1,'%s\n','')
    fprintf(1,'%s\n','             Vertical Components...')
    
%     fns=[char(statn(istat)) '.ae'];
    fzz=[char(statn(istat)) '.vz'];
    tmp=importdata(fzz);
    acczzori=tmp.data;
    tmpc=char(tmp.textdata);    
    
    % Original time increment (s)
    dt(istat)=str2num(tmpc(14:22));  
    
    % Number of samples, window length and time serie
    ntori=length(acczzori);
    tlori=(ntori-1)*dt(istat);
    timeori=(0:dt(istat):tlori);
    timeint=(0:dtot:tlori);
    tl=(ntot-1)*dtot;
    time=(0:dtot:tl);
    ntint=length(timeint);
    
    % Detrend and remove base line
%     acczzori=detrend(acczzori,'constant');
%     zzp(1)=nearest((arrtime(istat,2)-arrtime(istat,1))/dtot);
%     zzp(2)=nearest((arrtime(istat,2)+6-arrtime(istat,1))/dtot);
%     acczzori=detrend(acczzori,'linear',zzp);
    % Taper window
%     acczzori=acczzori.*tukeywin(ntori,0.05);
%     acczzori=detrend(acczzori,'linear');
%     % Taper window
%     acczzori=acczzori.*tukeywin(ntori,0.);
    
    % Interpolate seismogram (dtot)
    acczzint=interp1(timeori,acczzori,timeint);
    
    % Taper window
    acczzint=acczzint'.*tukeywin(ntint,0.1);
    
    % Add zeros to match the theoretical P-wave arrival time
    acczz(1:ntot)=0.0;
    nzero=nearest(arrtime(istat,1)/dtot);    
    if nzero+ntint < ntot
        acczz(nzero+1:nzero+ntint)=acczzint(1:ntint);
    else
        acczz(nzero+1:ntot)=acczzint(1:ntot-nzero);
    end

    % Taper window
%     acczz=acczz'.*tukeywin(ntot,0.1);
    
    % Compute velocity (integrate acceleration)
    velzzold = 0;
    for it=1:ntot               
        % Compute displacements (integrate velocities)
        if it > 1
            velzz(it) = acczz(it) * dtot + velzzold;
        else
            velzz(it) = 0.0;
        end
        velzzold = velzz(it);
    end
    
    % Write output seismograms
    if wdata == 1
        fprintf(1,'%s\n','             Write sliprate ascii files...')
        if ifield == 1
            fileoutns=fopen('obs_velo_zz','a');
            fprintf(fileoutns,'%14.6e\n',acczz);
        else
            fileoutns=fopen('obs_disp_zz','a');
            fprintf(fileoutns,'%14.6e\n',velzz);
        end
    end
    
    % Filter observed seismograms
    if bpf == 1
        fprintf(1,'              > Lowpass filter seismograms (f < %4.2f Hz) \n',highcut)
        fprintf(fid1,'              > Lowpass filter seismograms (f < %4.2f Hz) \n',highcut);
    end
    if bpf == 2
        fprintf(1,'              > Highpass filter seismograms (f > %4.2f Hz) \n',lowcut')
        fprintf(fid1,'              > Highpass filter seismograms (f > %4.2f Hz) \n',lowcut');
    end
    if bpf == 3
        fprintf(1,'              > Bandpass filter seismograms (%4.2f < f < %4.2f Hz) \n',lowcut,highcut')
        fprintf(fid1,'              > Bandpass filter seismograms (%4.2f < f < %4.2f Hz) \n',lowcut,highcut');
    end
    
    if ifield == 1
        field=acczz;
    else
        field=velzz;
    end
    
    if bpf ~= 0
        
        % Compute lowpass filter
        ratbobs=highcut/(1/(2*dtot));
        [Bfd,Afd]=butter(4,ratbobs,'low');
        Lowobs=impz(Bfd,Afd);

        % Compute highpass filter
        ratbobs=lowcut/(1/(2*dtot));
        [Bfd,Afd]=butter(4,ratbobs,'high');
        Highobs=impz(Bfd,Afd);
        %
        if bpf == 1  % Lowpass filter
            % Filter NS component
            tempo=conv(field,Lowobs);
            field(1:ntot)=tempo(1:ntot);
        end
        if bpf == 2  % Highpass filter
            % Filter NS component
            tempo=conv(field,Highobs);
            field(1:ntot)=tempo(1:ntot);
        end
        if bpf == 3  % Bandpass filter
            % Filter NS component
            tempo1=conv(field,Lowobs);
            tempo=conv(tempo1,Highobs);
            field(1:ntot)=tempo(1:ntot);
        end

        if ifield == 1
            acczz=field;
        else
            velzz=field;
        end
    end
    clear tempo
    
    % P and S waves arrival time ticks
    natp=nearest(arrtime(istat,1)/dtot)+1;
    nats=nearest(arrtime(istat,2)/dtot)+1;
    
        if ifield == 1
        vmax=max(acczz);
        vmin=min(acczz);
        absmax=max(abs(acczz));
        refp=acczz(natp);
        refs=acczz(nats);
    else
        vmax=max(velzz);
        vmin=min(velzz);
        absmax=max(abs(velzz));
        refp=velzz(natp);
        refs=velzz(nats);
    end
    
    ltick=0.15*absmax;
    tickpx=[(natp-1)*dtot (natp-1)*dtot];
    tickpy=[refp+ltick*0.5 refp-ltick*0.5];
    ticksx=[(nats-1)*dtot (nats-1)*dtot];
    ticksy=[refs+ltick*0.5 refs-ltick*0.5];

    h=figure(istat+20);
    %     plot(timeint,acczzint,'r'); hold on
    %     plot(timeori,acczzori,'g'); hold on
    if ifield == 1
        plot(time,acczz,'b'); hold on
        xlabel('Time (s)');
        ylabel('Velocity (cm/s)');
    else
        plot(time,velzz,'b'); hold on
        xlabel('Time (s)');
        ylabel('Displacement (cm)');
    end
    plot(tickpx,tickpy,'r'); hold on
    text(tickpx(1),tickpy(1)+0.1*ltick,'P','FontSize',12); hold on
    plot(ticksx,ticksy,'r'); hold on
    text(ticksx(1),ticksy(1)+0.1*ltick,'S','FontSize',12); hold on
    text(tl*0.8,0.8*vmax,char(statn(istat)),'FontSize',12); hold on
    grid on
    axis([0 tl vmin-0.1*absmax vmax+0.1*absmax]);
    title('Vertical Component');
    
    if ifield == 1
        fname=['./figure/' char(statn(istat)) '.vzz.ps'];
    else
        fname=['./figure/' char(statn(istat)) '.dzz.ps'];
    end
 %   print(h,'-depsc2',fname)

    clear acczzori tmp tmpc acczz
end

if wdata == 1; fclose(fileoutns); end
%if wdata == 1; fclose(fileoutew); end
%if wdata == 1; fclose(fileoutzz); end
