
%% load xlsx
[data,text]=xlsread('..\raw\Cluster12_0.707.xlsx','Sheet2','A2:B4097');
%% convert
tic
fid=fopen('..\cache\cluster.csv','w');
fclose(fid);

mzresult=data(:,1);
fid=fopen('..\cache\cluster.csv','a');
likelyrhomzresult=diag(likelyrho(diag(mzresult)));
%likelyrhomzresult=mzresult;
fprintf(fid,sprintf('%0.20f,',likelyrhomzresult));
fprintf(fid,sprintf('\n'));
fclose(fid);
toc
mzresult=data(:,2);
fid=fopen('..\cache\cluster.csv','a');
likelyrhomzresult=diag(likelyrho(diag(mzresult)));
%likelyrhomzresult=mzresult;
fprintf(fid,sprintf('%0.20f,',likelyrhomzresult));
fclose(fid);
toc