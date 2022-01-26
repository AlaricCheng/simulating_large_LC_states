

if 1
load('..\raw\big_circuit3_filenames.mat')

%% print fidelity
clc
disp(size(big_circuit3_filenames))
count=0;
tomopack={};
for ii = big_circuit3_filenames
    count=count+1;
    if (count >48 )&& (count < 1200-48+1)
        continue
    end
    load(ii{1});
    disp([num2str(count) ': ' num2str(fidelity) ','])
    tomopack{end+1}={'description,count,tomoResult,rho0,rho,purity,rhoideal',count,tomoResult,rho0,rho,purity,finalstate*finalstate'};
end
end
%% define observable
I=[1 0;0 1];
Z=[1 0;0 -1];
X=[0 1;1 0];
Y=[0 -1i;1i 0];

O=struct();
O.M1=[1 0;0 0];
O.M2=[0 0;0 1];
% O.M3=X;
% O.M4=X;
% O.M5=Y;
% O.M6=Y;
O.M3=Z;
O.M4=Z;
O.M5=Z;
O.M6=Z;

c=[1 1 1/2 -1/2 1/2 -1/2];

E1_o=struct();

for ii = 1:6
    % E1_o.(['M' num2str(ii)])=kron(O.(['M' num2str(ii)]),kron(I,Z));
    E1_o.(['M' num2str(ii)])=kron(kron(Z,I),O.(['M' num2str(ii)]));
end

E2_o=kron(Z,kron(Z,Z));

%%
% clc

allTurn=struct();
allTurn.count=0;
allTurn.value=[];
allTurn.trace=[];
allTurn.E1=zeros(1,6);
allTurn.E2=zeros(1,6);
allTurn.mzresult=struct();
for index_measure = [1,2,3,4,5,6,10,11,12]
    allTurn.mzresult.(['Z' num2str(index_measure)])=zeros(1,8);
end

currentTurn=struct();
currentTurn.E1=zeros(1,6);
currentTurn.E2=zeros(1,6);

fid=fopen('..\cache\big_circuit3_check_all25turn_likelyrho.txt','w');
fclose(fid);


figureid=ceil(2147483646*rand());
count=0;
for ii = big_circuit3_filenames
    count=count+1;
    %if (count >= 1200-48+1)
    %    continue %放弃最后一轮数据, 因为程序bug只跑了两个
    %end
    load(ii{1}); % -> mzresult
    disp(fidelity)
    disp([mzresult;abs(finalstate').^2])
    
    mzresultraw=mzresult;
    mzresultraw(:)=0;
    size(singleShotEvents,2);%40000  ,1 4 or 3
    for eventii = singleShotEvents
        %[0;0;0;1]
        bb=eventii(1)+eventii(2)*2+eventii(3)*4;
        if size(eventii,1)==4
            bb=bb+eventii(4)*8;
        end
        mzresultraw(bb+1)=mzresultraw(bb+1)+1;
    end
    mzresultraw=mzresultraw/size(singleShotEvents,2);

    fid=fopen('big_circuit3_check_all25turn_likelyrho.txt','a');
    fprintf(fid,sprintf('\n===%d',big_circuit3_value.index_measure));
    likelyrhomzresult=diag(lib.likelyrho(diag(mzresult)));
    likelyrhomzresult=likelyrhomzresult';
    size(mzresult)
    size(likelyrhomzresult)
    if big_circuit3_value.index_measure<=36
        fprintf(fid,sprintf('\n%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n',[likelyrhomzresult;abs(finalstate').^2]'));
    else
        fprintf(fid,sprintf('\n%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n',[likelyrhomzresult;abs(finalstate').^2]'));
    end
    
    fclose(fid);
    
    

    % % 25 turns * 9 circuits * 40000

    % index_measure=big_circuit3_value.index_measure;
    % % index_measure in [1,2,3,4,5,6,10,11,12]
    % % ~ rho 1-6, O X Y Z
    % index_measure_to_index=1:12;
    % index_measure_to_index(10:12)=[3,5,1];

    % index = index_measure_to_index(index_measure);

    % allTurn.mzresult.(['Z' num2str(index_measure)])=allTurn.mzresult.(['Z' num2str(index_measure)])+mzresult;

    % if index_measure <=8
    %     % circuit2
    %     currentTurn.E2(index)=trace(E2_o*diag(mzresult));
    % else
    %     % circuit1
    %     currentTurn.E1(index)=trace(E1_o.(['M' num2str(index)])*diag(mzresult));
    %     currentTurn.E1(index+1)=trace(E1_o.(['M' num2str(index+1)])*diag(mzresult));
    % end

    % if mod(count,9)~=0
    %     continue
    % end

    % allTurn.count=allTurn.count+1;
    % allTurn.trace(end+1)=sum(c.*currentTurn.E1.*currentTurn.E2);

    % for index_measure = [1,2,3,4,5,6,10,11,12]
    %     index = index_measure_to_index(index_measure);
    %     mzresult=allTurn.mzresult.(['Z' num2str(index_measure)])/allTurn.count;
    %     if index_measure <=8
    %         % circuit2
    %         allTurn.E2(index)=trace(E2_o*diag(mzresult));
    %     else
    %         % circuit1
    %         allTurn.E1(index)=trace(E1_o.(['M' num2str(index)])*diag(mzresult));
    %         allTurn.E1(index+1)=trace(E1_o.(['M' num2str(index+1)])*diag(mzresult));
    %     end
    % end

    % allTurn.value(end+1)=sum(c.*allTurn.E1.*allTurn.E2);

    % figure(figureid);
    % ax1=subplot(2,1,1);
    % plot(allTurn.value)
    % ax2=subplot(2,1,2);
    % plot(allTurn.trace)
    
end

