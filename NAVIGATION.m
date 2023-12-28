[n1,n2]=uigetfile('*.txt','Select Export matches');
T=strcat(n2,n1);
T=readtable(T);
format long g
%%=====================================================
Data=T.Var2;
Data(14:14:end)=[];
Data=reshape(Data,[],31);
e=Data(3,:);
t0=Data(4,:);
i=Data(5,:);
Rate_of_Right_Ascen=Data(6,:);
a=Data(7,:).^2;
Right_Ascen_at_Week=Data(8,:);
Argument_of_Perigee=Data(9,:);
Mean_Anom=Data(10,:);
GM=3.986004418e14;
we=7.2921151467e-5;
%----------------------DATUM------------------------
% WGS-84(1984)
A=6378137;	  B=6356752.3142;  E2=(A^2-B^2)/(A^2);
% h=input('input h = ');
phi=input('input phi  d= ');
Landa=input('input Landa d= ');
t=([phi;Landa]).*pi/180;
 phi=t(1,1);LAndA=t(2,1);PP=phi;LL=LAndA;
%---------------Ex------------------------------
h=0;
% phi=32*pi/180;LAndA=51*pi/180;

%------------------------------------------------
N=A^2/sqrt(A^2*cos(phi)^2+B^2*sin(phi)^2);
X_CT=(N+h)*cos(phi)*cos(LAndA);
Y_CT=(N+h)*cos(phi)*sin(LAndA);
Z_CT=(N*(B^2/A^2)+h)*sin(phi);
%---------------------------------------------------
% Greenwich True Sidereal Time at Longitude 0.0°
% 13h 43m 30.759s = 13.7252107292 h
GAST=205.8781609381*pi/180;
xp=0;yp=0;
S=Rotation(2,-xp)*Rotation(1,-yp)*Rotation(3,GAST);
%---------------------------------------------------
ch=menu('','SKYPLOT','PRNJ GROUND TRACK')
switch ch
    case 1
        Day_of_week=input('input Day_of_week =')-1;
        T1=input('input start time =')*3600+Day_of_week*24*3600;
        T2=input('input End time =')*3600+Day_of_week*24*3600;
        t=T1:20:T2;
        for j=1:31
            tk=t-t0(j);
            Mk=Mean_Anom(j)+(sqrt(GM/a(j)^3) ).*tk;
            omegak=Right_Ascen_at_Week(j)+(Rate_of_Right_Ascen(j)-we).*tk-we.*t0(j);
            w=Argument_of_Perigee(j);
            E0=Mk;
            for iii=1:5
                Ek=Mk+e(j)*sin(E0);
                E0=Ek;
            end
            Ek=wrapTo2Pi(Ek);
            r_orb=[a(j)*(cos(Ek)-e(j));a(j)*sin(Ek)*sqrt(1-e(j)^2);zeros(1,max(size(Ek)))];
            r=Rotation(1,-i(j))*Rotation(3,-w)*r_orb;
            r_ct=[];
            c=[];
            for k=1:max(size(Ek));
                r_ct(:,k)=S*Rotation(3,-omegak(k))*r(:,k);
           
            end
            dr=r_ct-[X_CT;Y_CT;Z_CT]*ones(1,k);
            for k=1:max(size(Ek));
                c(k)=norm(dr(:,k));
            end
            r_LA=[1 0 0;0 -1 0;0 0 1]*Rotation(2,phi-pi/2)*Rotation(3,LAndA-pi)*dr;
            Az=atan2(r_LA(2,:),r_LA(1,:));
            v=asin(r_LA(3,:)./c);
            Az=Az(find(v>0))*180/pi;
            v=v(find(v>0))*180/pi;max(v);
            v=-v./90 +1;
            hold on
            plot(v.*sind(Az),v.*cosd(Az),'.')
        end
        plotEllipseRotated([1],[1],[0,0],[0])
    case 2
        
        
        j=input('input PRN number =');
        for jj=1:max(size(j))
            w=Argument_of_Perigee(j(jj));
            Ek=0:0.01:10*pi;
            Mk=Ek-e(j(jj))*sin(Ek);
            tk=(Mk-Mean_Anom(j(jj)))/(sqrt(GM/a(j(jj))^3) );
            t=tk+t0(j(jj));
            omegak=Right_Ascen_at_Week(j(jj))+(Rate_of_Right_Ascen(j(jj))-we).*tk-we.*t0(j(jj));
            r_orb=[a(j(jj))*(cos(Ek)-e(j(jj)));a(j(jj))*sin(Ek)*sqrt(1-e(j(jj))^2);zeros(1,max(size(Ek)))];
            r=Rotation(1,-i(j(jj)))*Rotation(3,-w)*r_orb;
            r_ct=[];
            for k=1:max(size(Ek));
                r_ct(:,k)=S*Rotation(3,-omegak(k))*r(:,k);
            end
            for k=1:size(r_ct,2)
                X=r_ct(1,k);
                Y=r_ct(2,k);
                Z=r_ct(3,k);
                %kartezian be monhani
                p=sqrt(X^2+Y^2);
                phi0=atan((Z/p)/(1-E2));
                N0=A/(1-E2*sin(phi0)^2)^.5;
                h0=p/cos(phi0)-N0;
                Landa=wrapTo2Pi(atan2(Y,X));
                repitition=0;delta_phi=1;delta_h=1;
                while abs(delta_phi)>10^-5 && abs(delta_h)>10^-5
                    N=A/(1-E2*sin(phi0)^2)^.5;
                    h=(p/cos(phi0))-N;
                    phi=atan((Z/p)*((1-((E2*N)/(N+h)))^-1));
                    delta_phi=phi-phi0;delta_h=h-h0;phi0=phi;h0=h;
                    repitition=repitition+1;
                end
                PHI(k)=phi;LANDA(k)=Landa;
            end
            %wgs84
            Landa=[LANDA LL];
            Phi=[PHI PP];
            q=log(tan(pi/4 +Phi./2).*((1-sqrt(E2).*sin(Phi))./(1+sqrt(E2).*sin(Phi))).^(sqrt(E2)/2));
            x=A*(Landa);
            y=A*q;
            hold on
            plot(x,y,'.');
            title('in Mercator')
        end   
end