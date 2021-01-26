%This is a function to simulate backwater flow through a network with an
%arbitrary network structure. 

function varargout=arbbwnpo_v10(netstruct,AtchMiss,Qin,p,delzmatrix,nargs);

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

%Nov. 20, 2020, This script will be updated to include the spatial area
%variation. 
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
    
    
    %Build a set of reasonable slopes to solve for.
    wssguess=[-1 linspace(-1E-3,1E-3,200)];
        for k=1:r1-1;
            %recall water surface at this transect, and find possible water
            %surfaces at the upstream transect.
            zeta=printout(i).simout(end-k+1,6);
            zetaupguess=zeta+wssguess*printout(i).simout(end-k+1,3)*-1;
            %find this and the upstream transect... maybe not necessary.
            trans=AtchMiss(printout(i).simout(end-k+1,2)).tdata_meters;
            uptrans=AtchMiss(printout(i).simout(end-k,2)).tdata_meters;
            %Find the transstats for each transect.
            transtatsEH=AtchMiss(printout(i).simout(end-k+1,2)).transtatsEH;
            uptranstatsEH=AtchMiss(printout(i).simout(end-k,2)).transtatsEH;
            
            printout(i).simout(end-k+1,4)=nanmin(transtatsEH(:,3));
            eta=printout(i).simout(end-k+1,4);
            printout(i).simout(end-k+1,5)=interp1(transtatsEH(:,1),transtatsEH(:,5),zeta,'linear','extrap');%surface_width
            H=zeta-eta;%This is not used in modeling, but is recorded... for posterity.
            CSA=interp1(transtatsEH(:,1),transtatsEH(:,2),zeta,'linear','extrap');
            hy_radius=interp1(transtatsEH(:,1),transtatsEH(:,7),zeta,'linear','extrap');
            printout(i).simout(end-k+1,7)=H;
            

            %etaup=nanmin(uptranstatsEH(:,3));

            
            %S_bed=(eta-etaup)/printout(i).simout(end-k+1,3);
            CSA_up=interp1(uptranstatsEH(:,1),uptranstatsEH(:,2),zetaupguess,'linear','extrap');
            CSA_up(CSA_up<0)=1;
            CSA_use=(CSA_up+CSA)/2;
            
            hy_radius_up=interp1(uptranstatsEH(:,1),uptranstatsEH(:,7),zetaupguess,'linear','extrap');
            hy_radius_up(hy_radius_up<0)=1;
            hy_radius_use=(hy_radius+hy_radius_up)/2;
            
            U_use=printout(i).reachQ./CSA_use;
            Fr2=U_use.^2/9.8./(hy_radius_use);
            Cf=netstruct(i).Cf;
            Sf_use=Cf*Fr2;
            
            
            %U_up=printout(i).reachQ./CSA_up;
            %Fr2_up=U_up.^2/9.8./(hy_radius_up);
            %Sf_up=Cf*Fr2_up;
            
            %Sf_use=(Sf_up+Sf)/2;
            
            Q2overA=(printout(i).reachQ.^2)/CSA;
            upQ2overA=(printout(i).reachQ.^2)./CSA_up;
            Wchange=1./(9.8*CSA_use).*(Q2overA-upQ2overA)/printout(i).simout(end-k+1,3);
            
            if isnan(sum(Wchange))==1;
                break
            end
            wsschoices=-1*Wchange-Sf_use-wssguess;

            [jnk,ind]=find(wsschoices<0,1);
            if size(ind,2)==0;
                ind=20;
            end
            ind2=min(ind+5,length(wsschoices));
            wss=interp1(wsschoices(1:ind2),wssguess(1:ind2),0,'linear','extrap');
            zetanew=interp1(wssguess,zetaupguess,wss,'linear','extrap');
            U=interp1(wssguess,U_use,wss,'linear','extrap');
            CSA=interp1(wssguess,CSA_use,wss,'linear','extrap');

        
            printout(i).simout(end-k+1,10)=wss;
            %printout(i).simout(end-k,4)=etaup;
            printout(i).simout(end-k,6)=zetanew;
            printout(i).simout(end-k+1,8)=U;
            printout(i).simout(end-k+1,9)=Cf;
            printout(i).simout(end-k+1,11)=CSA;
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
