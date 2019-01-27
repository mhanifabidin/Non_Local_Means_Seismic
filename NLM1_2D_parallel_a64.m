clear
addpath('C:\Temp\Seismic_filters\NLMeans\NLM_code_2D','C:\Temp\Seismic_filters\NLMeans\NLM_code_2D\SegyMAT','C:\Temp\Seismic_filters\NLMeans\NLM_code_2D\SegyMAT\GUI');
parpool('local',12);
seisfile = '1445810_8601G_Reg_onshore_Regional';
[data] = ReadSegyFast(strcat(seisfile,'.sgy'));
u_scale = max(max(abs(data))); 
data_in = data./u_scale;
clearvars data;
clearvars i;
t0 = clock;
parfor i = 1:12
    searchw = 9;
    simw = 19;
    rangew = searchw + simw;
    [nx ny]=size(data_in);
    ny_part = floor(ny/12);
    pointer1=((i-1)*ny_part)+1-searchw;
    pointer2=(i*ny_part)+searchw;
    h_filt = 0.0015;
    if i == 1
        pointer1=1;
    end
    if i == 12
        pointer2=ny;
    end
    ny = ((pointer2-pointer1)+1);
    temp = data_in(1:nx,pointer1:pointer2);
    [w] = mexnlm2D(temp,searchw,simw,h_filt,nx,ny);
    if i == 12
        w_sel2(i,:,:) = w(1:nx,1+searchw:ny);
    elseif i == 1
        w_sel(i,:,:) = w(1:nx,1:ny-searchw);
    else
        w_sel(i,:,:) = w(1:nx,1+searchw:ny-searchw);
    end
end
t1 = clock;
et1 = etime(t1,t0);
clearvars data_in temp w;
for i = 1:11
    if i == 1
        w_combi=squeeze(w_sel(i,:,:));
    else
        w_combi=cat(2,w_combi,squeeze(w_sel(i,:,:)));
    end
end
w_combi=cat(2,w_combi,squeeze(w_sel2(12,:,:)));
clearvars w_sel w_sel2;
w_out = w_combi.*u_scale;
clearvars w_combi;
% save(strcat(seisfile,'_NLmeans_9_19_0_001.mat'),'w_out','-v7.3');
[data,Trhead,SegyHead] = ReadSegy(strcat(seisfile,'.sgy'));
WriteSegyStructure(strcat(seisfile,'_NLmeans_9_19_0_001.sgy'),SegyHead,Trhead,w_out,'revision',0,'dsf',1);
clearvars w_out
quit force
