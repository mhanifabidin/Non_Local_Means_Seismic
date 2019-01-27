%[data,Trhead,SegyHead] = ReadSegy('99_GNL91-112_PSTM_r6000_noAA_AGC.sgy');
%[data2,Trhead2,SegyHead2] = ReadSegy('610135_803007_h_repro.sgy');
%data = data(250:650,2200:2800);

data = read_segy_file('99_GNL91-112_PSTM_r6000_noAA.sgy', ...
    {'format','ibm'},{'times',2000,2400},{'traces','cdp >= 1100 & cdp <=1800'}, ...
    {'ignoreshift',true})
data_NLM = read_segy_file('99_GNL91-112_PSTM_r6000_noAA_NLmeans_50_3_0_001.sgy', ...
    {'format','ibm'},{'times',2000,2400},{'traces','cdp >= 1100 & cdp <=1800'}, ...
    {'ignoreshift',true})
diff_NLM = data-data_NLM.traces

subplot(2,2,1)
s_cplot(data,{'annotation','cdp'},{'title','PSTM'},{'fontsize',11},...
    {'figure','old'},{'time_lines',[]},{'annotation','time'},{'interpol','v5cubic'})
subplot(2,2,2)
s_cplot(data_NLM,{'annotation','cdp'},{'title','sim win 3'},{'fontsize',11},...
    {'figure','old'},{'time_lines',[]},{'annotation','time'},{'interpol','no'})
subplot(2,2,3)
s_cplot(diff_NLM,{'annotation','cdp'},{'title','sim win 3'},{'fontsize',11},...
    {'figure','old'},{'time_lines',[]},{'annotation','time'},{'interpol','no'})
% subplot(2,1,3)
% s_cplot(diff_NLM,{'annotation','cdp'},{'title','difference'},{'fontsize',11},...
%    {'figure','old'},{'time_lines',[]},{'annotation','time'})


% [data_agc500,tgain]=s_gain(data,{'type','trace'},{'average','mean'},...
%     {'wlength',500});
% [data_agc250,tgain]=s_gain(data,{'type','trace'},{'average','mean'},...
%     {'wlength',250});
% [data_NLM_agc500,tgain]=s_gain(data_NLM,{'type','trace'},{'average','mean'},...
%     {'wlength',500});
% [data_NLM_agc250,tgain]=s_gain(data_NLM,{'type','trace'},{'average','median'},...
%     {'wlength',250});
% [diff_NLM_agc250,tgain]=diff_NLM(data_NLM,{'type','trace'},{'average','median'},...
%     {'wlength',250});


% subplot(1,3,1)
% s_cplot(data_agc250,{'annotation','cdp'},{'title','PSTM'},{'fontsize',11},...
%     {'figure','old'},{'time_lines',[]},{'annotation','time'})
% subplot(1,3,2)
% s_cplot(data_NLM_agc250,{'annotation','cdp'},{'title','NLM sw 50 nw 3 h 0.001'},{'fontsize',20},...
%     {'figure','old'},{'time_lines',[]},{'annotation','time'})
% subplot(1,3,3)
%s_cplot(diff_NLM_agc250,{'annotation','cdp'},{'title','NLM sw 50 nw 3 h 0.001'},{'fontsize',20},...
%    {'figure','old'},{'time_lines',[]},{'annotation','time'})


% subplot(1,3,3)
% s_cplot(gdata_NLM,{'annotation','cdp'},{'title','NLM sw 50 nw 3 h 0.001'},{'fontsize',20},...
%     {'figure','old'},{'time_lines',[]},{'annotation','time'})
%       


   
% [nx ny]=size(data)
% [data_agc,tgain] = s_gain(data,'type','trace','average','median','wlength',500);
% figure(1);
% imagesc(data); colormap(redblue); caxis([-2.25 2.25])
% figure(2);
% imagesc(data_agc); colormap(redblue); caxis([-2.25 2.25])
