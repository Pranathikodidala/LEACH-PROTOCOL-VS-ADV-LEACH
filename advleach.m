%%%%%%Advanced LEACH%%%%%%
clear;
xm=300;
ym=300;
sink.x=0.5*xm;
sink.y=0.5*ym;
n = 200;
p=0.1;
Eo=0.5;
ch=n/10;
ETX=50*0.000000001;
ERX=50*0.000000001;
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
Efs1=Efs/10;
Emp1=Emp/10;
EDA=5*0.000000001;
a=Eo/2;
rmax=500;
h=100;
s=2;
do=sqrt(Efs/Emp);
do1=sqrt(Efs1/Emp1);
for h=1:1
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
Et=0;
end
for i=1:1:n
S(i).xd=rand(1,1)*xm;
S(i).yd=rand(1,1)*ym;
S(i).G=0;
S(i).type='N';
S(i).E=Eo;
Et=Et+S(i).E;
figure(h*1)
plot(S(i).xd,S(i).yd,'bo');
text(S(i).xd+1,S(i).yd-0.5,num2str(i));
hold on;
end
plot(S(i).xd,S(i).yd,'b^','MarkerSize', 5, 'MarkerFaceColor', 'b');
text(S(n+1).xd+1,S(n+1).yd-0.5,num2str(n+1));
title("Nodes formation")
hold off ;
countCHs=0;
cluster=1;
flag_first_DEAD=0;
flag_teenth_DEAD=0;
flag_all_DEAD=0;
DEAD=0;
first_DEAD=0;
teenth_DEAD=0;
all_DEAD=0;
ALLIVE=n;
for r=0:1:rmax
if(mod(r, round(1/p) )==0)
for i=1:1:n
S(i).G=0;
end
end
DEAD=0;
for i=1:1:n
if (S(i).E<=0)
DEAD=DEAD+1;
if (DEAD==1)
if(flag_first_DEAD==0)
first_DEAD=r;
flag_first_DEAD=1;
end
end
if(DEAD==0.1*n)
if(flag_teenth_DEAD==0)
teenth_DEAD=r;
flag_teenth_DEAD=1;
end
end
if(DEAD==n)
if(flag_all_DEAD==0)
all_DEAD=r;
flag_all_DEAD=1;
end
end
end
if S(i).E>0
S(i).type='N';
end
end
STATS.DEAD(r+1)=DEAD;
STATS.ALLIVE(r+1)=ALLIVE-DEAD;
countCHs=0;
cluster=1;
if S(i).type=='C' && S(i).E>a
for j=1:1:ch
countCHs=countCHs+1;
S(i).type='C';
S(i).G=round(1/p)-1;
C(cluster).xd=S(i).xd;
C(cluster).yd=S(i).yd;
distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
C(cluster).distance=distance;
C(cluster).id=i;
X(cluster)=S(i).xd;
Y(cluster)=S(i).yd;
cluster=cluster+1;
distance;
if (distance>do)
S(i).E=S(i).E- ( (ETX+EDA)*(4000) +Emp*4000*(distance*distance*distance*distance ));
end
if (distance<=do)
S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Efs*4000*(distance * distance ));
end
end
else
Etotal=0;
for i=1:n
if S(i).E>0
Etotal=Etotal+S(i).E;
end
end
Eavg=Etotal/n;
STATS.TotalENERGY(r+1)=Etotal;
STATS.AvgENERGY(r+1)=Eavg;
for i=1:1:n
if(S(i).E>0)
temp_rand=rand;
if ( (S(i).G)<=0)
if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
if(S(i).E>Eavg)
countCHs=countCHs+1;
S(i).type='C';
S(i).G=round(1/p)-1;
C(cluster).xd=S(i).xd;
C(cluster).yd=S(i).yd;
distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
C(cluster).distance=distance;
C(cluster).id=i;
X(cluster)=S(i).xd;
Y(cluster)=S(i).yd;
cluster=cluster+1;
distance;
if (distance>do)
S(i).E=S(i).E- ( (ETX+EDA)*(4000) +Emp*4000*(distance*distance*distance*distance ));
end
if (distance<=do)
S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Efs*4000*(distance * distance ));
end
end
end
end
end
end
end
STATS.COUNTCHS(r+1)=countCHs;
for i=1:1:n
if ( S(i).type=='N' && S(i).E>0 )
if(cluster-1>=1)
min_distance=Inf;
min_distance_cluster=0;
for c=1:1:cluster-1
temp=min(min_distance,(sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 )+sqrt((S(n+1).xd-C(c).xd)^2 + (S(n+1).yd-C(c).yd)^2 )));
if ( temp<min_distance )
min_distance=temp;
min_distance_cluster=c;
end
end
min_distance;
if (min_distance>do1)
S(i).E=S(i).E- ( ETX*(4000) + Emp1*4000*( min_distance *min_distance *min_distance * min_distance));
end
if (min_distance<=do1)
S(i).E=S(i).E- ( ETX*(4000) + Efs1*4000*( min_distance * min_distance));
end
S(C(min_distance_cluster).id).E =S(C(min_distance_cluster).id).E- ( (ERX +EDA)*4000 );
S(i).min_distance=min_distance;
S(i).min_distance_cluster=min_distance_cluster;
else
min_distance=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
if (min_distance>do)
S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_distance *min_distance *min_distance * min_distance));
end
if (min_distance<=do)
S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_distance * min_distance));
end
end
end
end
En=0;
for i=1:n
    if S(i).E<=0
continue;
end
En=En+S(i).E;
end
ENERGY(r+1)=En;
STATS.ENERGY(h,r+1)=En;
end
for r=0:rmax
STATS.DEAD(h+1,r+1)=sum(STATS.DEAD(:,r+1))/h;
STATS.ALLIVE(h+1,r+1)=sum(STATS.ALLIVE(:,r+1))/h;
STATS.ENERGY(h+1,r+1)=sum(STATS.ENERGY(:,r+1))/h;
STATS.COUNTCHS(h+1,r+1)=sum(STATS.COUNTCHS(:,r+1))/h;
end
%%%%% ACTUAL LEACH %%%%%
xm=300;
ym=300;
sink.x=xm;
sink.y=ym;
n=200;
p=0.1;
Eo=0.5;
ETX=50*0.000000001;
ERX=50*0.000000001;
Efs=10e-12;
Emp=0.0013e-12;
EDA=5*0.000000001;
rmax=500;
do=sqrt(Efs/Emp);
Et=0;
for h=1:1
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
Et=0;
for i=1:1:n
S(i).xd=rand(1,1)*xm;
XR(i)=S(i).xd;
S(i).yd=rand(1,1)*ym;
YR(i)=S(i).yd;
distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
S(i).distance=distance;
S(i).G=0;
S(i).type='N';
S(i).E=Eo;
Et=Et+S(i).E;
figure(h*2)
plot(S(i).xd,S(i).yd,'b^','MarkerSize', 5, 'MarkerFaceColor', 'b');
text(S(i).xd+1,S(i).yd-0.5,num2str(i));
hold on;
end
plot(S(n+1).xd,S(n+1).yd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(S(n+1).xd+1,S(n+1).yd-0.5,num2str(n+1));
title("Nodes formation");
hold off ;
countCHs=0;
cluster=1;
flag_first_DEAD=0;
flag_half_DEAD=0;
flag_all_DEAD=0;
first_DEAD=0;
half_DEAD=0;
all_DEAD=0;
ALLIVE=n;
for r=0:1:rmax
if(mod(r, round(1/p) )==0)
for i=1:1:n
S(i).G=0;
S(i).cl=0;
end
end
DEAD=0;
for i=1:1:n
if (S(i).E<=0)
DEAD=DEAD+1;
if (DEAD==1)
if(flag_first_DEAD==0)
first_DEAD=r;
flag_first_DEAD=1;
end
end
if(DEAD==0.5*n)
if(flag_half_DEAD==0)
half_DEAD=r;
flag_half_DEAD=1;
end
end
if(DEAD==n)
if(flag_all_DEAD==0)
all_DEAD=r;
flag_all_DEAD=1;
end
end
end
if S(i).E>0
S(i).type='N';
end
end
STATISTICS(2).DEAD(h,r+1)=DEAD;
STATISTICS(2).ALLIVE(h,r+1)=ALLIVE-DEAD;
countCHs=0;
cluster=1;
for i=1:1:n
if(S(i).E>0)
temp_rand=rand;
if ( (S(i).G)<=0)
if ( temp_rand <= ( p/ ( 1 - p * mod(r,round(1/p)) )))
countCHs=countCHs+1;
S(i).type='C';
S(i).G=round(1/p)-1;
C(cluster).xd=S(i).xd;
C(cluster).yd=S(i).yd;
distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
C(cluster).distance=distance;
C(cluster).id=i;
X(cluster)=S(i).xd;
Y(cluster)=S(i).yd;
cluster=cluster+1;
distance;
if (distance>do)
S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*(distance*distance*distance*distance ));
end
if (distance<=do)
S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Efs*4000*( distance * distance ));
end
end
end
end
end
STATISTICS(2).COUNTCHS(h,r+1)=countCHs;
for i=1:1:n
if ( S(i).type=='N' && S(i).E>0 )
if(cluster-1>=1)
min_distance=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
min_distance_cluster=0;
for c=1:1:cluster-1
temp=min(min_distance,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
if ( temp<min_distance )
min_distance=temp;
min_distance_cluster=c;
end
end
if(min_distance_cluster~=0)
min_distance;
if (min_distance>do)
S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_distance * min_distance*min_distance * min_distance));
end
if (min_distance<=do)
S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_distance * min_distance));
end
S(C(min_distance_cluster).id).E = S(C(min_distance_cluster).id).E- ( (ERX +EDA)*4000 );
else
min_distance;
if (min_distance>do)
S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_distance * min_distance *min_distance * min_distance));
end
if (min_distance<=do)
S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_distance * min_distance));
end
end
S(i).min_distance=min_distance;
S(i).min_distance_cluster=min_distance_cluster;
else
min_distance=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
if (min_distance>do)
S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_distance * min_distance *min_distance * min_distance));
end
if (min_distance<=do)
S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_distance * min_distance));
end
end
end
end
En=0;
for i=1:n
if S(i).E<=0
continue;
end
En=En+S(i).E;
end
ENERGY(r+1)=En;
STATISTICS(2).ENERGY(h,r+1)=En;
end
first_DEAD_LEACH(h)=first_DEAD;
half_DEAD_LEACH(h)=half_DEAD;
all_DEAD_LEACH(h)=all_DEAD;
end
for r=0:rmax
STATISTICS(2).DEAD(h+1,r+1)=sum(STATISTICS(2).DEAD(:,r+1))/h;
STATISTICS(2).ALLIVE(h+1,r+1)=sum(STATISTICS(2).ALLIVE(:,r+1))/h;
STATISTICS(2).ENERGY(h+1,r+1)=sum(STATISTICS(2).ENERGY(:,r+1))/h;
STATISTICS(2).COUNTCHS(h+1,r+1)=sum(STATISTICS(2).COUNTCHS(:,r+1));
end
first_DEAD=sum(first_DEAD_LEACH)/h;
half_DEAD=sum(half_DEAD_LEACH)/h;
all_DEAD=sum(all_DEAD_LEACH)/h;
r=0:rmax;
figure(3)
plot(r,STATS.DEAD(h+1,r+1))
xlabel('Rounds');
ylabel('DEAD Nodes');
title('Rounds VS DEAD Nodes');
hold on
plot(r,STATISTICS(2).DEAD(h+1,r+1));
legend('Advanced LEACH','LEACH')
hold off
figure(4)
plot(r,STATS.ALLIVE(h+1,r+1))
xlabel('Rounds');
ylabel('Live Nodes')
title('Rounds VS Live Nodes');
hold on
plot(r,STATISTICS(2).ALLIVE(h+1,r+1))
legend('Advanced LEACH','LEACH')
hold off
figure(5)
plot(r,STATS.ENERGY(h+1,r+1))
xlabel('Rounds')
ylabel('Average ENERGY')
title('Average Residual ENERGY')
hold on
plot(r,STATISTICS(2).ENERGY(h+1,r+1))
legend('Advanced LEACH','LEACH')
hold off
