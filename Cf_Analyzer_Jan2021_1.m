%This program will iterate over many pairings of Cf for the Atchafalaya and
%Mississippi River side and return the partitioning in the R_pre model. 


load('AtchMiss_v12.mat')
load('netstruct50M_v13.mat');
load('netstruct16M_v12.mat');
options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunEvals',1500)
load('ValidatorData_v1.mat');

Qset=[18300]';%[10000:5000:25000]';
Cf=0.01;%This is just to fill the voice in v6Cf, it is a dummy variable.

AtchCfset=linspace(0.002,0.004,9)';
MRCfset=linspace(0.002,0.004,9)';
tic;
for h=1:size(Qset,1);
    for i=1:size(AtchCfset,1);
        for j=1:size(MRCfset,1);
            Q=Qset(h);
            AtchCf=AtchCfset(i);
            MRCf=MRCfset(j);


            netstruct16M(5).Cf=MRCf;
            for k=1:4;
                netstruct16M(k).Cf=AtchCf;
            end
            netstruct50M(13).Cf=MRCf;
            for k=[1:12 14];
                netstruct50M(k).Cf=AtchCf;
            end
            [poptim16M,fval16M,exitflag] = fminunc(@(p)arbbwnpo_v10(netstruct16M,AtchMiss,Q,p,delzmatrix16M,1),[0.3 0.5],options);
            [delz2, printout16M]=arbbwnpo_v10(netstruct16M,AtchMiss,Q,poptim16M,delzmatrix16M,2)
            Atchplot16M=[printout16M(3).simout(2:end,:);printout16M(2).simout(2:end,:);printout16M(1).simout(2:end,:)];
            
%             [poptim50M,fval50M,exitflag] = fminunc(@(p)arbbwnpo_v10(netstruct50M,AtchMiss,Q,p,delzmatrix50M,1),[0.3 0.6 0.20 0.1 0.3 0.6],options);
%             [delz2, printout50M]=arbbwnpo_v10(netstruct50M,AtchMiss,Q,poptim50M,delzmatrix50M,2)
%             Atchplot50M=[printout50M(7).simout(1:end,:);printout50M(12).simout(2:end,:);printout50M(11).simout(2:end,:);printout50M(5).simout(2:end,:);printout50M(9).simout(2:end,:);printout50M(2).simout(2:end,:);printout50M(1).simout(2:end,:)];


            Cfdelz16(i,j,h)=delz2;
            AtchMisspartition16(i,j,h)=poptim16M(1);
            
            Cfdelz50(i,j,h)=delz2;
            AtchMisspartition50(i,j,h)=poptim50M(1);
            
            Qrr16=Qset(h);%RedRiverLanding
            Zrr16=interp1(ValidatorData(1).QZ16(:,1),ValidatorData(1).QZ16(:,2),Qrr16,'linear','extrap');
            %plot(sum(Atchplot16M(1:9,3))/1000,Zsp16,'ko')
            dz(i,j,1,h)=printout16M(3).simout(1,6)-Zrr16;

            Qsp16=printout16M(3).reachQ;%Simmesport
            Zsp16=interp1(ValidatorData(2).QZ16(:,1),ValidatorData(2).QZ16(:,2),Qsp16,'linear','extrap');
            %plot(sum(Atchplot16M(1:9,3))/1000,Zsp16,'ko')
            dz(i,j,2,h)=Atchplot16M(9,6)-Zsp16;

            Qmv16=printout16M(3).reachQ;%Melville
            Zmv16=interp1(ValidatorData(3).QZ16(:,1),ValidatorData(3).QZ16(:,2),Qmv16,'linear','extrap');
            %plot(sum(Atchplot16M(1:20,3))/1000,Zmv16,'ko')
            dz(i,j,3,h)=Atchplot16M(20,6)-Zmv16;

            Qks16=printout16M(3).reachQ;%Krotz Springs
            Zks16=interp1(ValidatorData(4).QZ16(:,1),ValidatorData(4).QZ16(:,2),Qks16,'linear','extrap');
            %plot(sum(Atchplot16M(1:26,3))/1000,Zks16,'ko')
            dz(i,j,4,h)=Atchplot16M(26,6)-Zks16;

            Qaa16=printout16M(3).reachQ;%Atchafalaya %AtchMiss(97)
            Zaa16=interp1(ValidatorData(5).QZ16(:,1),ValidatorData(5).QZ16(:,2),Qaa16,'linear','extrap');
            %plot(sum(Atchplot16M(1:37,3))/1000,Zaa16,'ko')
            dz(i,j,5,h)=Atchplot16M(37,6)-Zaa16;

            Qkp16=printout16M(2).reachQ;%Keel Boat Pass %AtchMiss(5903)
            Zkp16=interp1(ValidatorData(6).QZ16(:,1),ValidatorData(6).QZ16(:,2),Qkp16,'linear','extrap');
            %plot(sum(Atchplot16M(1:54,3))/1000,Zkp16,'ko')
            dz(i,j,6,h)=Atchplot16M(54,6)-Zkp16;

            Qmc16=printout16M(1).reachQ;%Morgan City %AtchMiss(55)
            Zmc16=interp1(ValidatorData(7).QZ16(:,1),ValidatorData(7).QZ16(:,2),Qmc16,'linear','extrap');
            %plot(sum(Atchplot16M(1:63,3))/1000,Zmc16,'ko')
            dz(i,j,7,h)=Atchplot16M(63,6)-Zmc16;
            
            
%             %Now for R50.
% 
%             Qrr50=Qset(h);%RedRiverLanding
%             Zrr50=interp1(ValidatorData(1).QZ50(:,1),ValidatorData(1).QZ50(:,2),Qrr50,'linear','extrap');
%             %plot(sum(Atchplot50M(1:9,3))/1000,Zsp50,'ko')
%             dz50(i,j,1,h)=printout50M(13).simout(1,6)-Zrr50;
% 
%             Qsp50=printout50M(3).reachQ;%Simmesport
%             Zsp50=interp1(ValidatorData(2).QZ50(:,1),ValidatorData(2).QZ50(:,2),Qsp50,'linear','extrap');
%             %plot(sum(Atchplot50M(1:9,3))/1000,Zsp50,'ko')
%             dz50(i,j,2,h)=Atchplot50M(9,6)-Zsp50;
% 
%             Qmv50=printout50M(3).reachQ;%Melville
%             Zmv50=interp1(ValidatorData(3).QZ50(:,1),ValidatorData(3).QZ50(:,2),Qmv50,'linear','extrap');
%             %plot(sum(Atchplot50M(1:20,3))/1000,Zmv50,'ko')
%             dz50(i,j,3,h)=Atchplot50M(20,6)-Zmv50;
% 
%             Qks50=printout50M(3).reachQ;%Krotz Springs
%             Zks50=interp1(ValidatorData(4).QZ50(:,1),ValidatorData(4).QZ50(:,2),Qks50,'linear','extrap');
%             %plot(sum(Atchplot50M(1:26,3))/1000,Zks50,'ko')
%             dz50(i,j,4,h)=Atchplot50M(26,6)-Zks50;
% 
%             Qaa50=printout50M(3).reachQ;%Atchafalaya %AtchMiss(97)
%             Zaa50=interp1(ValidatorData(5).QZ50(:,1),ValidatorData(5).QZ50(:,2),Qaa50,'linear','extrap');
%             %plot(sum(Atchplot50M(1:37,3))/1000,Zaa50,'ko')
%             dz50(i,j,5,h)=Atchplot50M(37,6)-Zaa50;
% 
%             Qkp50=printout50M(2).reachQ;%Keel Boat Pass %AtchMiss(5903)
%             Zkp50=interp1(ValidatorData(6).QZ50(:,1),ValidatorData(6).QZ50(:,2),Qkp50,'linear','extrap');
%             %plot(sum(Atchplot50M(1:54,3))/1000,Zkp50,'ko')
%             dz50(i,j,6,h)=Atchplot50M(54,6)-Zkp50;
% 
%             Qmc50=printout50M(1).reachQ;%Morgan City %AtchMiss(55)
%             Zmc50=interp1(ValidatorData(7).QZ50(:,1),ValidatorData(7).QZ50(:,2),Qmc50,'linear','extrap');
%             %plot(sum(Atchplot50M(1:63,3))/1000,Zmc50,'ko')
%             dz50(i,j,7,h)=Atchplot50M(63,6)-Zmc50;

            disp(toc)
            disp(['h=' num2str(h)])
            disp(['i=' num2str(i)])
            disp(['j=' num2str(j)])
        end
    end
end
save('Cf_Analyzer_v9Q18300.mat') 

figure
subplot(1,3,1)
set(gca,'Units','inches','Position',[0.5 0.5 2.25 2],'FontSize',12)
image(MRCfset,AtchCfset,abs(AtchMisspartition16(:,:,1)),'CDataMapping','scaled')
hold on
plot(0.0025,0.00325,'ko')
daspect([1 1 1])
colorbar
caxis([0.15 0.35])
xlabel('C_{f Miss} (-)')
ylabel('C_{f Atch} (-)')
title('f_A')
set(gca,'Ydir','normal')

subplot(1,3,2)
set(gca,'Units','inches','Position',[3 0.5 2.25 2],'FontSize',12)
image(MRCfset,AtchCfset,sqrt(mean(dz(:,:,:,1).^2,3)),'CDataMapping','scaled')
hold on
plot(0.0025,0.00325,'ko')
daspect([1 1 1])
colorbar
caxis([0 2])
xlabel('C_{f Miss} (-)')
ylabel('C_{f Atch} (-)')
title('RMSE gauges (m)')
set(gca,'Ydir','normal')


% for i=1:10;
%     for j=1:10;
%         dzavg(i,j)=
%         