%This is a function to simulate backwater flow through a network with an
%arbitrary network structure. 

function varargout=arbbwnpo_v8(netstruct,AtchMiss,Qin,p,delzmatrix,nargs);

%The network structure will have this order
%The first network networkstruct(1).transects will be the call numbers of
%the transects in the Atch library that go from the ocean to the upstream
%most point. The downstream condition networkstruct(i).dsc is the water
%surface elevation at the downstream node. If networkstruct(i).dsc=0, then
%it is at the downstream end. If networkstruct(i).dsc=-999, it's downstream
%end is within the network and we will search for the elevation from
%previously mapped transects.

%netstruct(i).Qfun is a function handle describing the discharge through
%the network link. The coefficients in .Qfun will be 1 if the reach uses
%the a given p fraction, -1 if the reach uses 1-p fraction and 0 if the
%reaches has nothing to do with the p fraction.

%p is a vector with one less entry than the number of network links. 

%netstruct(i).trans is a nx2 matrix, the first column is the call number of
%the transect in the Atch library, and the second column is the distance
%from the transect to the transect upstream.
p=abs(p);
printout=struct('reachno',{},'simout',{},'reachQ',{});
for i = 1:length(netstruct);
    printout(i).reachno=i;
    for j=1:length(netstruct(i).trans);
        printout(i).simout(j,1)=j;
        printout(i).simout(j,2:3)=netstruct(i).trans(j,:); %This is the transect library and the distance.
        %printout(i).simout(j,4)=AtchMiss(netstruct(i).trans(j,1)).avgeta;
        %printout(i).simout(j,5)=AtchMiss(netstruct(i).trans(j,1)).effwidth;
    end
    if netstruct(i).dsc~=-999;
        printout(i).simout(end,6)=netstruct(i).dsc;
    else %Search the previously modelled transects for the starting point of this transect
        %Then use it as a downstream boundary. Note that program will break
        %if the first transect i=1 doesn't have a downstream boundary set.
        searchind=printout(i).simout(end,2);
        for k=1:i-1;
            match=find(printout(k).simout(:,2)==searchind);
            if length(match)>0;
            printout(i).simout(end,6)=printout(k).simout(match,6);
            end
        end               
        
    end
    
    Qcoeff=netstruct(i).Qfun;
    [r c]=size(Qcoeff);
    reachQ=Qin*ones(r,1);
    for j=1:r;
        for k=1:c;
            if Qcoeff(j,k)==1;
                reachQ(j)=reachQ(j)*p(k);
            elseif Qcoeff(j,k)==0;
                reachQ(j)=reachQ(j);
            elseif Qcoeff(j,k)==-1;
                reachQ(j)=reachQ(j)*(1-p(k));
            else
                break
            end
            
        end
        
    end
    reachQ=sum(reachQ');
    printout(i).reachQ=reachQ;
    %Now calculate Bw
    [r1 c1]=size(printout(i).simout);
    for k=1:r1-1;
        zeta=printout(i).simout(end-k+1,6);
        trans=AtchMiss(printout(i).simout(end-k+1,2)).tdata_meters;
        %[CSA,eta,H,effwidth]=tranalyze_v1(trans(:,1:2),zeta);
        transtatsEH=AtchMiss(printout(i).simout(end-k+1,2)).transtatsEH;
        %printout(i).simout(end-k+1,4)=min(trans(:,2));
        printout(i).simout(end-k+1,4)=nanmin(transtatsEH(:,3));
        eta=printout(i).simout(end-k+1,4);
        printout(i).simout(end-k+1,5)=interp1(transtatsEH(:,1),transtatsEH(:,5),zeta,'linear','extrap');%surface_width
        H=zeta-eta;
        CSA=interp1(transtatsEH(:,1),transtatsEH(:,2),zeta,'linear','extrap');
        hy_radius=interp1(transtatsEH(:,1),transtatsEH(:,7),zeta,'linear','extrap');
        U=printout(i).reachQ/CSA;
        Fr2=U^2/9.8/(hy_radius);
        
%             Fr2=0.1;
%         end
        %Cf=[1/.407*log(11*H/(3*D90))].^-2;
        Cf=netstruct(i).Cf;
        printout(i).simout(end-k+1,7)=H;
        printout(i).simout(end-k+1,8)=U;
        printout(i).simout(end-k+1,9)=Cf;
        printout(i).simout(end-k+1,11)=CSA;
               
    %     if Fr2>0.5;
    %         Fr2=0.5;
    %     end
        
        %nSlopes=100;
        %Slopes=logspace(-10,-1,nSlopes);
        %newzetas=printout(i).simout(end-k+1,3).*Slopes+zeta;
        %newtrans=AtchMiss(printout(i).simout(end-k,2)).tdata_meters;
        uptrans=AtchMiss(printout(i).simout(end-k,2)).tdata_meters;
        uptranstatsEH=AtchMiss(printout(i).simout(end-k,2)).transtatsEH;
        %etaup=min(uptrans(:,2));;
        etaup=nanmin(uptranstatsEH(:,3));
        etachange=10;
%         while etachange>0;
%             detadxtst=(eta-etatst)/printout(i).simout(end-k+1,3);
%             dhdxtst=(-1*detadxtst-Cf*Fr2)/(1-Fr2);
%             wss=(detadxtst+dhdxtst);
%             if wss>0;
%                 wss=-1E-8;
%             end
%             zetanew=zeta+wss*printout(i).simout(end-k+1,3)*-1;
%             newtrans=AtchMiss(printout(i).simout(end-k,2)).tdata_meters;
%             [CSAtst,etatst1,Htst,effwidthtst]=tranalyze_v1(newtrans,zetanew);
%             etachange=abs(etatst1-etatst);
%             etatst=etatst1;
%         end
        %for l=1:2;
            detadxtst=(eta-etaup)/printout(i).simout(end-k+1,3);
            dhdxtst=(-1*detadxtst-Cf*Fr2)/(1-Fr2);
            wss=(detadxtst+dhdxtst);
             if wss>10E-4;
                 wss=10E-4;
                 %wss=-1*Cf*Fr2;
             elseif wss<-10E-4;
                wss=-10E-4;
            end
            zetanew=zeta+wss*printout(i).simout(end-k+1,3)*-1;
            %newtrans=AtchMiss(printout(i).simout(end-k,2)).tdata_meters;
%             if zetanew>30;
%                 22;
%             end
            %etatst=interp1(uptranstats(:,1),uptranstats(:,3),zetanew,'linear','extrap');
            
            %[CSAtst,etatst,Htst,effwidthtst]=tranalyze_v1(newtrans,zetanew);
        %end
        %UtstMB=printout(i).reachQ./CSAtst;%Velocity just from mass balance
        %Cftst=[1/.407*log(11*Htst/3*D90)].^-2;
        %detadxtst=(eta-etatst)/printout(i).simout(end-k+1,3);
        %dhdxtst=(-1*detadxtst-Cf*Fr2)/(1-Fr2);
        %realdzetadx=interp1(Slopes+dhdxtst+detadxtst,Slopes,0,'linear','extrap');
        %zetanew=zeta+realdzetadx*printout(i).simout(end-k+1,3);
        printout(i).simout(end-k+1,10)=wss;
        printout(i).simout(end-k,4)=etaup;
        printout(i).simout(end-k,6)=zetanew;
%         
%         detadx=(printout(i).simout(end-k,4)-printout(i).simout(end-k+1,4))/printout(i).simout(end-k+1,3);
%         dHdx=(-1*detadx-Cf*Fr2)/(1-Fr2);
%         dzeta=(dHdx+detadx)*printout(i).simout(end-k,3)*-1/2;%Only integrate halfway through the step!
%         if dzeta<0;
%             dzeta=10^-5;
%         end
%         newzeta=printout(i).simout(end-k+1,6)+dzeta;
% %         if newzeta<printout(i).simout(end-k,4);
% %             newzeta=printout(i).simout(end-k,4)+0.3;
% %         end
%         H2=newzeta-mean(printout(i).simout(end-k:end-k+1,4));
%         
%     %      if H2<0.3;
%     %          H2=0.3;
%     %      end
%          U2=printout(i).reachQ/mean(printout(i).simout(end-k:end-k+1,5))/H2;
%           Fr22=U2^2/9.8/H2;
%     %      if Fr22>0.5;
%     %          Fr22=0.5;
%     %      end
%          dHdx2=(-1*detadx-Cf*Fr22)/(1-Fr22);
%          dzeta2=(dHdx2+detadx)*printout(i).simout(end-k,3)*-1/2;
%          if dzeta2<0;
%              dzeta2=10^-5;
%          %elseif dzeta2>
%          end
%          printout(i).simout(end-k,6)=printout(i).simout(end-k+1,6)+dzeta+dzeta2;
    end
end


%OK, we simulated flow in all of the channels, now to calculate the
%elevation differences at the bifurcations. For channels 1 to n-1, we will
%look at all channels with larger index numbers for a matching transect
%call number at the top of the section. To do this right, the transect that
%goes furtest upstream (Mississippi in this case) will be the last
%transect.

% for i=1:length(printout)-1;
%     bifind=printout(i).simout(1,2);
%     bifzeta=printout(i).simout(1,6);
%     for k=i+1:length(printout);
%         match2=find(printout(k).simout(:,2)==bifind);
%         if length(match)>0;
%         bifzeta2=printout(k).simout(match2,6);
%         break
%         end
%     end
%     delz2(i)=(bifzeta-bifzeta2)^2;
% end      
[r c]=size(delzmatrix);
for i=1:r;
    zeta1=printout(delzmatrix(i,1)).simout(delzmatrix(i,2),6);
    zeta2=printout(delzmatrix(i,3)).simout(delzmatrix(i,4),6);
    delz(i)=(zeta2-zeta1)^2;
end
    delz2=sum(delz);
    if max(abs(p))>1;
        delz2=99999;
    end
    if nargs==1;
        varargout{1}=delz2;
        %varargout{2}=0;
    elseif nargs==2;
        varargout{1}=delz2;
        varargout{2}=printout;
    end
