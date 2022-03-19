


load('..\raw\circuits_filenames.mat')

%% print fidelity
if 1
tic
clc
disp(size(big_circuit3_filenames))
count=0;
tomopack={};
fidelitypack=[];
for ii = big_circuit3_filenames
    count=count+1;
    if (count >48 )&& (count < 1200-48+1)
        continue
    end
    circuitname=split(ii{1},'\');
    circuitname=['..\raw\circuits\' char(circuitname(end))];
    load(circuitname);
    disp([num2str(count) ': ' num2str(fidelity) ','])
    tomopack{end+1}={'description,count,tomoResult,rho0,rho,purity,rhoideal',count,tomoResult,rho0,rho,purity,finalstate*finalstate'};
    fidelitypack(end+1)=fidelity;
end


fid=fopen('..\cache\big_circuit3_fidelity.txt','w');
fprintf(fid,sprintf('%0.20f,',fidelitypack));
fclose(fid);
toc
end
%% write likelyrho data
tic
fid=fopen('..\cache\big_circuit3_check_temp_likelyrho_1.txt','w');
fclose(fid);
fid=fopen('..\cache\big_circuit3_check_temp_likelyrho_2.txt','w');
fclose(fid);


count=0;
for ii = big_circuit3_filenames
    count=count+1;
    circuitname=split(ii{1},'\');
    circuitname=['..\raw\circuits\' char(circuitname(end))];
    load(circuitname);
    % -> mzresult
    %disp(fidelity)
    %disp([mzresult;abs(finalstate').^2])
    
    % mzresultraw=mzresult;
    % mzresultraw(:)=0;
    % size(singleShotEvents,2);%40000  ,1 4 or 3
    % for eventii = singleShotEvents
    %     %[0;0;0;1]
    %     bb=eventii(1)+eventii(2)*2+eventii(3)*4;
    %     if size(eventii,1)==4
    %         bb=bb+eventii(4)*8;
    %     end
    %     mzresultraw(bb+1)=mzresultraw(bb+1)+1;
    % end
    % mzresultraw=mzresultraw/size(singleShotEvents,2);

    if big_circuit3_value.index_measure<=36
        fid=fopen('..\cache\big_circuit3_check_temp_likelyrho_1.txt','a');
    else
        fid=fopen('..\cache\big_circuit3_check_temp_likelyrho_2.txt','a');
    end
    fprintf(fid,sprintf('===%d\n',big_circuit3_value.index_measure));
    likelyrhomzresult=diag(likelyrho(diag(mzresult)));
    likelyrhomzresult=likelyrhomzresult';
    %size(mzresult)
    %size(likelyrhomzresult)
    if big_circuit3_value.index_measure<=36
        fprintf(fid,sprintf('%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n',[likelyrhomzresult;abs(finalstate').^2]'));
    else
        fprintf(fid,sprintf('%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n',[likelyrhomzresult;abs(finalstate').^2]'));
    end
    
    fclose(fid);
    
    

    
end
toc
