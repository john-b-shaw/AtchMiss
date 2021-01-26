%This mfile is written to display the key figures in a manuscript
%tentatively titled "Natural and human controls on the
%Atchafalaya-Mississippi avulsion. Written by John Shaw with help from
%Kashauna Mason.

load('AtchMiss_v12.mat')%This is a library of all bathymetric transects used here.
%Here are the network structure files. Each structure contains a separate
%transect (trans), discharge function Qfun, a downstream condition (dsc)
%with -999 showing that the transect ends in another transect, and 0
%indicating that the system ends at 0 elevation (sea level).
load('netstruct16M_v12.mat'); 
load('netstruct50M_v13.mat');
load('netstruct16MD_v8.mat');
load('netstruct16MA_v5.mat');
load('netstruct16MM_v5.mat');
%load('netstruct50GL_v2.mat');
load('netstruct16GL_v2.mat');
load('ValidatorData_v1.mat'); %This is a database of 

Q=20000;%Here is the upstream discharge in m^3/s.

options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunEvals',1500)
%that last step set the options for the optimization


%Here is where the friction factors (C_f) are set for the Mississippi and
%Atchafalaya individually. Routines for setting these friction factors are
%described in a companion mfile. AtchCf and MRCf were varied through a
%range, and (a) the error associated with the gauges and (b) the partitioning fraction at Old RIver was assessed for each
%combination.
AtchCf=0.00325;
MRCf=0.0025;


netstruct16M(5).Cf=MRCf;
for i=1:4;
    netstruct16M(i).Cf=AtchCf;
end
%netstruct16M(2).Cf=0.00000;

netstruct50M(13).Cf=MRCf;
for i=[1:12 14];
    netstruct50M(i).Cf=AtchCf;
end
%netstruct50M(9).Cf=0.00000;

netstruct16MD(13).Cf=MRCf;
for i=[1:12 14];
    netstruct16MD(i).Cf=AtchCf;
end

netstruct16MM(5).Cf=MRCf;
for i=1:4;
    netstruct16MM(i).Cf=AtchCf;
end

netstruct16MA(5).Cf=MRCf;
for i=1:4;
    netstruct16MA(i).Cf=AtchCf;
end

% netstruct50GL(13).Cf=MRCf;
% for i=[1:12 14];
%     netstruct50GL(i).Cf=AtchCf;
% end

netstruct16GL(5).Cf=MRCf;
for i=1:4;
    netstruct16GL(i).Cf=AtchCf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now, simulate flow through each model. For each model, the first line
%optimizes flow, the second line prints out the optimized solution, and the
%third line builds a flow through the Atchafalaya Network for comparison
%and visualization.

[poptim16M,fval16M,exitflag] = fminunc(@(p)arbbwnpo_v10(netstruct16M,AtchMiss,Q,p,delzmatrix16M,1),[0.3 0.5],options);
[delz2, printout16M]=arbbwnpo_v10(netstruct16M,AtchMiss,Q,poptim16M,delzmatrix16M,2)
Atchplot16M=[printout16M(3).simout(1:end,:);printout16M(2).simout(2:end,:);printout16M(1).simout(2:end,:)];


[poptim50M,fval50M,exitflag] = fminunc(@(p)arbbwnpo_v10(netstruct50M,AtchMiss,Q,p,delzmatrix50M,1),[0.3 0.6 0.25 0.12 0.3 0.25],options);
[delz2, printout50M]=arbbwnpo_v10(netstruct50M,AtchMiss,Q,poptim50M,delzmatrix50M,2)
Atchplot50M=[printout50M(7).simout(1:end,:);printout50M(12).simout(2:end,:);printout50M(11).simout(2:end,:);printout50M(5).simout(2:end,:);printout50M(9).simout(2:end,:);printout50M(2).simout(2:end,:);printout50M(1).simout(2:end,:)];
%Atchplot50M=[printout50M(7).simout(2:end,:);printout50M(12).simout(2:end,:);printout50M(7).simout(2:end,:);printout50M(12).simout(2:end,:);printout50M(10).simout(2:end,:);printout50M(2).simout(2:end,:);printout50M(1).simout(2:end,:)];




[poptim16MD,fval16MD,exitflag] = fminunc(@(p)arbbwnpo_v10(netstruct16MD,AtchMiss,Q,p,delzmatrix16MD,1),[0.3 0.6 0.25 0.12 0.3 0.35],options);
[delz2, printout16MD]=arbbwnpo_v10(netstruct16MD,AtchMiss,Q,poptim16MD,delzmatrix16MD,2)
Atchplot16MD=[printout16MD(7).simout(1:end,:);printout16MD(12).simout(2:end,:);printout16MD(11).simout(2:end,:);printout16MD(5).simout(2:end,:);printout16MD(9).simout(2:end,:);printout16MD(2).simout(2:end,:);printout16MD(1).simout(2:end,:)];

[poptim16MA,fval16MA,exitflag] = fminunc(@(p)arbbwnpo_v10(netstruct16MA,AtchMiss,Q,p,delzmatrix16MA,1),[0.3 0.5],options);
[delz2, printout16MA]=arbbwnpo_v10(netstruct16MA,AtchMiss,Q,poptim16MA,delzmatrix16MA,2)
Atchplot16MA=[printout16MA(3).simout(1:end,:);printout16MA(2).simout(2:end,:);printout16MA(1).simout(2:end,:)];

[poptim16MM,fval16MM,exitflag] = fminunc(@(p)arbbwnpo_v10(netstruct16MM,AtchMiss,Q,p,delzmatrix16MM,1),[0.3 0.5],options);
[delz2_16MM, printout16MM]=arbbwnpo_v10(netstruct16MM,AtchMiss,Q,poptim16MM,delzmatrix16MM,2)
Atchplot16MM=[printout16MM(3).simout(1:end,:);printout16MM(2).simout(2:end,:);printout16MM(1).simout(2:end,:)];

[poptim16GL,fval16GL,exitflag] = fminunc(@(p)arbbwnpo_v10(netstruct16GL,AtchMiss,Q,p,delzmatrix16GL,1),[0.3 0.5],options);
[delz2, printout16GL]=arbbwnpo_v10(netstruct16GL,AtchMiss,Q,poptim16GL,delzmatrix16GL,2)
Atchplot16GL=[printout16GL(3).simout(1:end,:);printout16GL(2).simout(2:end,:);printout16GL(1).simout(2:end,:)];

%%%%%%%%%%%%%%%%%%%%%
%Now build validation figure.

%We need to validate 16M and 50M.

%RedRiverLanding NOTE THAT THIS ACTUALLY TARBERT LANDING DATA, WHICH IS
%UPSTREAM OF THE BIFURCATION
Qrr16=Q;
Zrr16=interp1(ValidatorData(1).QZ16(:,1),ValidatorData(1).QZ16(:,2),Qrr16,'linear','extrap');
%plot(0,Zrr16,'ko')
Qrr50=Q;
Zrr50=interp1(ValidatorData(1).QZ50(:,1),ValidatorData(1).QZ50(:,2),Qrr50,'linear','extrap');
%plot(0,Zrr50,'ro')
V16(1,:)=[0 Zrr16];
V50(1,:)=[0 Zrr50];


Qsp16=printout16M(3).reachQ;%Simmesport
Zsp16=interp1(ValidatorData(2).QZ16(:,1),ValidatorData(2).QZ16(:,2),Qsp16,'linear','extrap');
%plot(sum(Atchplot16M(1:9,3))/1000,Zsp16,'ko')

Qsp50=printout50M(7).reachQ;%Simmesport
Zsp50=interp1(ValidatorData(2).QZ50(:,1),ValidatorData(2).QZ50(:,2),Qsp50,'linear','extrap');
%plot(sum(Atchplot50M(1:9,3))/1000,Zsp50,'ro')

V16(2,:)=[sum(Atchplot16M(1:9,3))/1000 Zsp16];
V50(2,:)=[sum(Atchplot50M(1:9,3))/1000 Zsp50];


Qmv16=printout16M(3).reachQ;%Melville
Zmv16=interp1(ValidatorData(3).QZ16(:,1),ValidatorData(3).QZ16(:,2),Qmv16,'linear','extrap');
%plot(sum(Atchplot16M(1:20,3))/1000,Zmv16,'ko')

Qmv50=printout50M(7).reachQ;%Melville
Zmv50=interp1(ValidatorData(3).QZ50(:,1),ValidatorData(3).QZ50(:,2),Qmv50,'linear','extrap');
%plot(sum(Atchplot50M(1:20,3))/1000,Zmv50,'ro')
V16(3,:)=[sum(Atchplot16M(1:20,3))/1000 Zmv16];
V50(3,:)=[sum(Atchplot50M(1:20,3))/1000 Zmv50];

Qks16=printout16M(3).reachQ;%Krotz Springs
Zks16=interp1(ValidatorData(4).QZ16(:,1),ValidatorData(4).QZ16(:,2),Qks16,'linear','extrap');
%plot(sum(Atchplot16M(1:26,3))/1000,Zks16,'ko')

Qks50=printout50M(7).reachQ;%Krotz Springs
Zks50=interp1(ValidatorData(4).QZ50(:,1),ValidatorData(4).QZ50(:,2),Qks50,'linear','extrap');
%plot(sum(Atchplot50M(1:26,3))/1000,Zks50,'ro')

V16(4,:)=[sum(Atchplot16M(1:26,3))/1000 Zks16];
V50(4,:)=[sum(Atchplot50M(1:26,3))/1000 Zks50];

Qaa16=printout16M(3).reachQ;%Krotz Springs %AtchMiss(97)
Zaa16=interp1(ValidatorData(5).QZ16(:,1),ValidatorData(5).QZ16(:,2),Qaa16,'linear','extrap');
%plot(sum(Atchplot16M(1:37,3))/1000,Zaa16,'ko')

Qaa50=printout50M(7).reachQ;%Krotz Springs %AtchMiss(98)
Zaa50=interp1(ValidatorData(5).QZ50(:,1),ValidatorData(5).QZ50(:,2),Qaa50,'linear','extrap');
%plot(sum(Atchplot50M(1:37,3))/1000,Zaa50,'ro')

V16(5,:)=[sum(Atchplot16M(1:37,3))/1000 Zaa16];
V50(5,:)=[sum(Atchplot50M(1:37,3))/1000 Zaa50];

Qkp16=printout16M(2).reachQ;%Atchafalaya %AtchMiss(5903)
Zkp16=interp1(ValidatorData(6).QZ16(:,1),ValidatorData(6).QZ16(:,2),Qkp16,'linear','extrap');
%plot(sum(Atchplot16M(1:54,3))/1000,Zkp16,'ko')

Qkp50=printout50M(12).reachQ;%Atchafalaya %AtchMiss(5903)
Zkp50=interp1(ValidatorData(6).QZ50(:,1),ValidatorData(6).QZ50(:,2),Qkp50,'linear','extrap');
%plot(sum(Atchplot50M(1:54,3))/1000,Zkp50,'ro')

V16(6,:)=[sum(Atchplot16M(1:54,3))/1000 Zkp16];
V50(6,:)=[sum(Atchplot50M(1:54,3))/1000 Zkp50];

Qmc16=printout16M(1).reachQ;%Morgan City %AtchMiss(55)
Zmc16=interp1(ValidatorData(7).QZ16(:,1),ValidatorData(7).QZ16(:,2),Qmc16,'linear','extrap');
%plot(sum(Atchplot16M(1:63,3))/1000,Zmc16,'ko')

Qmc50=printout16M(1).reachQ;%Morgan City %AtchMiss(55)
Zmc50=interp1(ValidatorData(7).QZ50(:,1),ValidatorData(7).QZ50(:,2),Qmc50,'linear','extrap');
%plot(sum(Atchplot50M(1:63,3))/1000,Zmc50,'ro')

V16(7,:)=[sum(Atchplot16M(1:63,3))/1000 Zmc16];
V50(7,:)=[sum(Atchplot50M(1:63,3))/1000 Zmc50];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Look at the 50 km averages of width and elevation for these rivers.
mrb50=printout50M(13).simout;
mrb16=printout16M(5).simout;
rc50=cumsum(mrb50(:,3))/1000;
rc16=cumsum(mrb16(:,3))/1000;
bins=[0:50:500];
clear mrouts
for i=1:length(bins)-1;
    ind50=find(rc50>bins(i) & rc50<bins(i+1));
    ind16=find(rc16>bins(i) & rc16<bins(i+1));
    eta50=mean(mrb50(ind50,4));
    w50=mean(mrb50(ind50,5));
    csa50=mean(mrb50(ind50,11));
    slope50=mean(mrb50(ind50,10));
    eta16=mean(mrb16(ind16,4));
    w16=mean(mrb16(ind16,5));
    csa16=mean(mrb16(ind16,11));
    slope16=mean(mrb16(ind16,10));
    x=(bins(i)+bins(i+1))/2;
    mrouts(i,:)=[x eta16 eta50 eta50/eta16 w16 w50 w50/w16 csa16 csa50 slope16, slope50];
end


pgon1 = polyshape([0 98 98 0],[-50 -50 20 20]);
pgon2 = polyshape([149 171 171 149],[-50 -50 20 20]);
pgon11 = polyshape([0 98 98 0],[-50 -50 2E4 2E4]);
pgon21 = polyshape([149 171 171 149],[-50 -50 2E4 2E4]);
pgon12 = polyshape([0 98 98 0],[.1 .1 10 10]);
pgon22 = polyshape([149 171 171 149],[.1 .1 10 10]);

figure;
subplot(2,2,1)
set(gca,'Units','Inches','Position',[0.5 3 2.75 2],'FontSize',9)
hold on
plot(pgon1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
plot(pgon2,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
p1=plot(cumsum(Atchplot50M(:,3))/1000,Atchplot50M(:,6),'Color',[230 115 114]/256,'LineWidth',2)
p2=plot(cumsum(Atchplot16M(:,3))/1000,Atchplot16M(:,6),'k','LineWidth',2)
p3=plot(cumsum(Atchplot50M(:,3))/1000,Atchplot50M(:,4),'Color',[230 115 114]/256,'LineWidth',2)
p4=plot(cumsum(Atchplot16M(:,3))/1000,Atchplot16M(:,4),'k','LineWidth',2)
p5=plot(V16(:,1),V16(:,2),'ko')
p6=plot(V50(:,1),V50(:,2),'o','Color',[230 115 114]/256)
axis([0 240 -40 15])
%plot([0 98],[14 14],'k')%This is the leveed Atchafalaya Range
text(49,13.5,'Upper A.R.','HorizontalAlignment','center','FontSize',9)
%plot([98 149],[13.5 13.5],'k')% Unleveed, multichannel Atchafalaya
text(123.5,-38,'Lower A.R.','HorizontalAlignment','center','FontSize',9)
%plot([149 171],[14 14],'k')%The Delta
text(160,13.5,'Delta','HorizontalAlignment','center','FontSize',9)
%plot([171 235],[13.5 13.5],'k')%The Lakes
%text(203,-38,'Lakes and Below','HorizontalAlignment','center','FontSize',9)
ylabel('Elevation (m MSL)')
legend([p2 p1 p5 p6],{'R16','R50','R16 Data','R50 Data'},'Location','southwest')

subplot(2,2,3)
set(gca,'Units','Inches','Position',[0.5 0.5 2.75 2],'FontSize',9)
hold on
plot(pgon11,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
plot(pgon21,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
plot(cumsum(Atchplot50M(:,3))/1000,Atchplot50M(:,11),'Color',[230 115 114]/256,'LineWidth',2)
plot(cumsum(Atchplot16M(:,3))/1000,Atchplot16M(:,11),'k','LineWidth',2)
axis([0 240 0 2E4])
ylabel('A (m^2)')
xlabel('Distance Downstream of RRL (km)')
% subplot(2,1,2)
% set(gca,'Units','Inches','Position',[1 1 5.5 2],'FontSize',12)
% hold on
% plot(cumsum(Atchplot50M(:,3))/1000,Atchplot50M(:,5),'Color',[230 115 114]/256,'LineWidth',2)
% plot(cumsum(Atchplot16M(:,3))/1000,Atchplot16M(:,5),'k','LineWidth',2)
% axis([0 240 0 2E3])
% ylabel('Surface Width (m)')
% xlabel('Distance Downstream of RRL (km)')


%figure
subplot(2,2,2)
set(gca,'Units','Inches','Position',[4 3 2.75 2],'FontSize',9)
hold on
plot(cumsum(printout50M(13).simout(:,3))/1000,printout50M(13).simout(:,4),'-','Color',[230 115 114]/256)
plot(cumsum(printout50M(13).simout(:,3))/1000,printout50M(13).simout(:,6),'Color',[230 115 114]/256,'LineWidth',2)
plot(cumsum(printout16M(5).simout(:,3))/1000,printout16M(5).simout(:,4),'k')
plot(cumsum(printout16M(5).simout(:,3))/1000,printout16M(5).simout(:,6),'k','LineWidth',2)
plot(mrouts(:,1),mrouts(:,2),'d-','Color',[0.5 0.5 0.5],'LineWidth',2)
plot(mrouts(:,1),mrouts(:,3),'d-','Color',[234 186 186]/256,'LineWidth',2)
axis([0 475 -40 15])
ylabel('Elevation (m MSL)')
% subplot(2,1,2)
% set(gca,'Units','Inches','Position',[1 1 5.5 2],'FontSize',12)
% hold on
% plot(cumsum(printout50M(13).simout(:,3))/1000,printout50M(13).simout(:,5),'.','Color',[230 115 114]/256)
% plot(cumsum(printout16M(5).simout(:,3))/1000,printout16M(5).simout(:,5),'k.')
% plot(mrouts(:,1),mrouts(:,5),'o-','Color','k','LineWidth',2)
% plot(mrouts(:,1),mrouts(:,6),'o-','Color',[230 115 114]/256,'LineWidth',2)
% axis([0 475 200 1000])
% ylabel('Surface Width (m)')
% xlabel('Distance Downstream of RRL (km)')
subplot(2,2,4)
set(gca,'Units','Inches','Position',[4 0.5 2.75 2],'FontSize',9)
hold on
plot(cumsum(printout50M(13).simout(:,3))/1000,printout50M(13).simout(:,11),'-','Color',[230 115 114]/256)
plot(cumsum(printout16M(5).simout(:,3))/1000,printout16M(5).simout(:,11),'k-')
plot(mrouts(:,1),mrouts(:,8),'d-','Color',[0.5 0.5 0.5],'LineWidth',2)
plot(mrouts(:,1),mrouts(:,9),'d-','Color',[234 186 186]/256,'LineWidth',2)
axis([0 475 0 20000])
ylabel('A (m^2)')
xlabel('Distance Downstream of RRL (km)')

figure
subplot(2,1,1)
set(gca,'Units','Inches','Position',[1 4 5.5 2],'FontSize',9)
hold on
plot(pgon1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
plot(pgon2,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
pR50=plot(cumsum(Atchplot50M(:,3))/1000,Atchplot50M(:,6),'Color',[230 115 114]/256,'LineWidth',2)
pR16=plot(cumsum(Atchplot16M(:,3))/1000,Atchplot16M(:,6),'k','LineWidth',2)
pR16D=plot(cumsum(Atchplot16MD(:,3))/1000,Atchplot16MD(:,6),'Color',[255 205 101]/256,'LineWidth',2)%Yellowish
pR16A=plot(cumsum(Atchplot16MA(:,3))/1000,Atchplot16MA(:,6),'Color',[138 198 74]/256,'LineWidth',2)%Greenish
pR16M=plot(cumsum(Atchplot16MM(:,3))/1000,Atchplot16MM(:,6),'Color',[116 208 237]/256,'LineWidth',2)%Light Blue
pR16GL=plot(cumsum(Atchplot16GL(:,3))/1000,Atchplot16GL(:,6),'Color',[112 89 165]/256,'LineWidth',2)
p5=plot(V16(:,1),V16(:,2),'ko')
p6=plot(V50(:,1),V50(:,2),'o','Color',[230 115 114]/256)
axis([0 240 0 15])
%plot([0 98],[14 14],'k')%This is the leveed Atchafalaya Range
text(49,13.5,'Upper A.R.','HorizontalAlignment','center')
%plot([98 149],[13.5 13.5],'k')% Unleveed, multichannel Atchafalaya
text(123.5,13.5,'Lower A.R.','HorizontalAlignment','center')
%plot([149 171],[14 14],'k')%The Delta
text(160,13.5,'Delta','HorizontalAlignment','center')
%plot([171 235],[13.5 13.5],'k')%The Lakes
%text(203,14.5,'Lakes and Below','HorizontalAlignment','center')
legend([pR16 pR50 pR16D pR16A pR16M pR16GL],{'R\_pre','R\_post','R\_pre\_D','R\_pre\_A','R\_pre\_M','R\_pre\_GL'},'Location','east')

ylabel('Elevation (m MSL)')
xlabel('Distance Downstream of RRL (km)')


%Look at shear stress variations
subplot(2,1,2)
set(gca,'Units','Inches','Position',[1 1 5.5 2],'FontSize',9)
tau16=1000*Atchplot16M(2:end,8).^2.*Atchplot16M(2:end,9);
tau50=1000*Atchplot50M(2:end,8).^2.*Atchplot50M(2:end,9);
tau16D=1000*Atchplot16MD(2:end,8).^2.*Atchplot16MD(2:end,9);
tau16A=1000*Atchplot16MA(2:end,8).^2.*Atchplot16MA(2:end,9);
tau16M=1000*Atchplot16MM(2:end,8).^2.*Atchplot16MM(2:end,9);
tau16G=1000*Atchplot16GL(2:end,8).^2.*Atchplot16GL(2:end,9);

hold on
plot(pgon12,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
plot(pgon22,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
% tR50=plot(cumsum(Atchplot50M(2:end,3))/1000,1000*Atchplot50M(2:end,8).^2.*Atchplot50M(2:end,9),'Color',[230 115 114]/256,'LineWidth',2)
% tR16=plot(cumsum(Atchplot16M(2:end,3))/1000,1000*Atchplot16M(2:end,8).^2.*Atchplot16M(2:end,9),'Color','k','LineWidth',2)
% tR16D=plot(cumsum(Atchplot16MD(2:end,3))/1000,1000*Atchplot16MD(2:end,8).^2.*Atchplot16MD(2:end,9),'Color',[255 205 101]/256,'LineWidth',2)
% tR16A=plot(cumsum(Atchplot16MA(2:end,3))/1000,1000*Atchplot16MA(2:end,8).^2.*Atchplot16MA(2:end,9),'Color',[138 198 74]/256,'LineWidth',2)

tR50=plot(cumsum(Atchplot50M(2:end,3))/1000,tau50./tau16,'Color',[230 115 114]/256,'LineWidth',2)
tR16=plot(cumsum(Atchplot16M(2:end,3))/1000,ones(size(Atchplot16M(2:end,3))),'Color','k','LineWidth',2)
tR16D=plot(cumsum(Atchplot16MD(2:end,3))/1000,tau16D./tau16,'Color',[255 205 101]/256,'LineWidth',2)
tR16A=plot(cumsum(Atchplot16MA(2:end,3))/1000,tau16A./tau16,'Color',[138 198 74]/256,'LineWidth',2)
tR16M=plot(cumsum(Atchplot16MM(2:end,3))/1000,tau16M./tau16,'Color',[116 208 237]/256,'LineWidth',2)
tR16G=plot(cumsum(Atchplot16GL(2:end,3))/1000,tau16G./tau16,'Color',[112 89 165]/256,'LineWidth',2)
set(gca,'Yscale','log')
axis([0 240 0.1 10])


%axis([0 240 0 15])
ylabel('tau_b relative to R16 (-)')
xlabel('Distance Downstream of RRL (km)')

