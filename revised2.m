clc
clear
%%moteghayerat

MaxMin=input('-1 For Minimization/ 1 For Maximazation: ');
C=input('Coefficients of your problem:  ');
C=-MaxMin*C;
A1=input('Constraint coefficients (<=) : ');
A2=input('Constraint coefficients (>=) : ');
A3=input('Constraint coefficients (=) : ');
RHS1=input('Right hand side (<=) : ' );
RHS2=input('Right hand side (>=) : ');
RHS3=input('Right hand side (=) : ');
b=[RHS1;RHS3;RHS2];
b1=b;
SC=size(C);
SA1=size(A1);
SA2=size(A2);
SA3=size(A3);
A=[A1;A3;A2];
SA=size(A);
I=eye(SA(1));
B=I;
%syms M;
M(1,1:SA2(1)+SA3(1))=100;
CB=[zeros(1,SA1(1)) M];
P=[A B];
CB2=zeros(1,SA2(1));
B1=[];
%%jadval avalie
if ~isempty(A2)
for i=SA3(1)+SA1(1)+1:SA3(1)+SA2(1)+SA1(1)
    B1(:,i-SA3(1)-SA1(1))=-B(:,i);
end
CAT=[C CB CB2 0;A B B1 b];
else 
    CAT=[C CB CB2 0;A B b];
end
SCAT=size(CAT);
T=table();
for i=1:SCAT(2)
    T{:,i}=CAT(:,i);
end
V1=cell(1,SCAT(2));
for i=1:SA(2)
    V1{1,i}=strcat('X',num2str(i));
end
for i=SA(2)+1:SA(2)+SA1(1)
    V1{1,i}=strcat('S',num2str(i-SA(2)));
end
for i=SA(2)+SA1(1)+1:SA(2)+SA1(1)+SA2(1)+SA3(1)
    V1{1,i}=strcat('R',num2str(i-SA(2)));
end
for i=SA(2)+SA1(1)+SA2(1)+SA3(1)+1:SA(2)+SA1(1)+SA3(1)+(SA2(1)*2)
    V1{1,i}=strcat('S',num2str(i-SA(2)-SA2(1)));
end
V1{1,end}=strcat('RHS');
T.Properties.VariableNames=V1;
V2=cell(SCAT(1),1);
if MaxMin==-1
    V2{1,1}=strcat('W');
else
V2{1,1}=strcat('Z');
end
for i=2:SA1(1)+1
    V2{i,1}=strcat('S',num2str(i-1));
end
for i=SA1(1)+2:SA(1)+1
    V2{i,1}=strcat('R',num2str(i-1));
end
T.Properties.RowNames=V2;
disp(T)
if ~isempty(A2)|| ~isempty(A3)
     for i=SA1(1)+2:SCAT(1)
 CAT(1,:)=-100*CAT(i,:)+CAT(1,:);
     end
end
CB=CAT(1,SA(2)+1:SCAT(2)-SA2(1)-1);
C=CAT(1,1:SA(2));
CT=CAT(1,1:end-1);
CTsabet=-CT;
XB=b;
Z1=CAT(1,end);
P2=P;
P=[A B B1];
%P=CAT(2:end,1:end-1);
 T{(1:end),(1:end)}=CAT(:,:);
 disp(T);
%%Mohasebat Revised
%%
%%moteghayer vorodi
while 1
Mn=min(CT);
if Mn==0
    disp('Optimal solution')
    break
end
F=find(CT==Mn);
P1=max(P2(:,F),0);
if sum(P1)==0
    disp('Unbounded')
    break
end
D=XB./P1;
%Moteghayer khorogi
Mn2=min(D);
F2=find(D==Mn2);
SF2=size(F2);
%Tabahgeni
if SF2(1)>1
    disp('Degeneracy')
end
%amaliat
B(:,F2(1,1))=P(:,F);
XB=B\b;
CB(1,F2(1,1))=CTsabet(1,F);
Z=CB*XB+Z1;
C1=(CB/B)*A+C;
CBn=CB/B;
A1=B\A;
%vared kardan slack variable dar jadval
if ~isempty(A2)
    for i=SA3(1)+SA1(1)+1:SA3(1)+SA2(1)+SA1(1)
    CB2(1,i-SA3(1)-SA1(1))=CBn(1,i);
    end
    CB1=-CB2+M(1,1);
    CT=[C1 CBn CB1];
    for i=SA3(1)+SA1(1)+1:SA3(1)+SA2(1)+SA1(1)
    B2=inv(B);
    B1(:,i-SA3(1)-SA1(1))=-B2(:,i);
    end
    P2=[A1 inv(B) B1];
else
    CT=[C1 CBn];
    P2=[A1 inv(B)];
    MOS=find(CT==0);
    SMOS=size(MOS);
    if SMOS(2)>size(b,1)
        disp('Multiple optimal solution')
    end
end
CAT=[CT Z;P2 XB];
T.Properties.VariableNames=V1;
V2{F2(1,1)+1,1}=V1{1,F};
T.Properties.RowNames=V2;
T{(1:end),(1:size(CAT,2))}=CAT(:,:);
disp(T)
end
if ~isempty(A2)|| ~isempty(A3)
for i=1+SA(2)+SA1(1):SA(2)+SA1(1)+SA3(1)+SA2(1)
    CC(1,i-SA(2)-SA1(1))=CT(1,i);
    MNC=min(CC);
    if MNC==0
        disp('infeasible');
    end
end
end
if MaxMin==-1
    T1=table();
    T1{1,1}=-Z;
    V3=cell(1,1);
    V3{1,1}=strcat('Z = ');
    T1.Properties.RowNames=V3;
    V4=cell(1,1);
    V4{1,1}=strcat('Meghdar_Z');
    T1.Properties.VariableNames=V4;
    disp(T1)
end

%%


%%Sensivity Analysis b

G=input('show sensetivity analysis?yes(1)/no(2)');
if G==1
    S_b=size(b);
if S_b(1)==1
    
% sensivity analysis b1

for i = 0 : 1 : 200000
    ba = b;
    ba(1) = (i/100) - 1000;
    sen1 = B\ba;
    if sen1(1)>=0
                sen2(i+1,1) = (i/100) - 1000;
    else
        sen2(i+1,1) = NaN;
    end 
end

lb1 = min(sen2);
ub1 = max(sen2);

if lb1 == -1000
    ddd = ['sensivity analysis b1: ','-inf' , '<b1<=' , num2str(ub1)];
    disp(ddd)
elseif ub1 == 1000
     ddd = ['sensivity analysis b1: ', num2str(lb1) , '<=b1<' , 'inf'];
     disp(ddd)
elseif lb1 == -1000 && ub1 == 1000
     ddd = ['sensivity analysis b1: ', '-inf' , '<b1<' , 'inf'];
     disp(ddd)
else
   ddd= ['sensivity analysis b1: ', num2str(lb1) , '<=b1<=' , num2str(ub1)];
    disp(ddd)
end
end





if S_b(1)==2
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b1

for i = 0 : 1 : 200000
    ba = b;
    ba(1) = (i/100) - 1000;
    sen1 = B\ba;
    if sen1(1)>=0
        if sen1(2)>=0
                sen2(i+1,1) = (i/100) - 1000;
            else
                sen2(i+1,1) = NaN;
        end
    else
        sen2(i+1,1) = NaN;
    end 
end

lb1 = min(sen2);
ub1 = max(sen2);

if lb1 == -1000
    ddd = ['sensivity analysis b1: ','-inf' , '<b1<=' , num2str(ub1)];
    disp(ddd)
elseif ub1 == 1000
     ddd = ['sensivity analysis b1: ', num2str(lb1) , '<=b1<' , 'inf'];
     disp(ddd)
elseif lb1 == -1000 && ub1 == 1000
     ddd = ['sensivity analysis b1: ', '-inf' , '<b1<' , 'inf'];
     disp(ddd)
else
   ddd= ['sensivity analysis b1: ', num2str(lb1) , '<=b1<=' , num2str(ub1)];
    disp(ddd)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b2
for i = 0 : 1 : 200000
    ba = b;
    ba(2) = (i/100) - 1000;
    sen3 = B\ba;
    if sen3(1)>=0
        if sen3(2)>=0
                sen4(i+1,1) = (i/100) - 1000;
            else
                sen4(i+1,1) = NaN;
        end
    else
        sen4(i+1,1) = NaN;
    end  
end

lb2 = min(sen4);
ub2 = max(sen4);

if lb2 == -1000
    ddd2 = ['sensivity analysis b2: ','-inf' , '<b2<=' , num2str(ub2)];
    disp(ddd2)
elseif ub2 == 1000
     ddd2 = ['sensivity analysis b2: ', num2str(lb2) , '<=b2<' , 'inf'];
     disp(ddd2)
elseif lb2 == -1000 && ub2 == 1000
     ddd2 = ['sensivity analysis b2: ', '-inf' , '<b2<' , 'inf'];
     disp(ddd2)
else
   ddd2= ['sensivity analysis b2: ', num2str(lb2) , '<=b2<=' , num2str(ub2)];
    disp(ddd2)
end
end






if S_b(1)==3
    
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b1

for i = 0 : 1 : 200000
    ba = b;
    ba(1) = (i/100) - 1000;
    sen1 = B\ba;
    if sen1(1)>=0
        if sen1(2)>=0
            if sen1(3)>=0
                sen2(i+1,1) = (i/100) - 1000;
            else
                sen2(i+1,1) = NaN;
            end
        else
            sen2(i+1,1) = NaN;
        end
    else
        sen2(i+1,1) = NaN;
    end 
end

lb1 = min(sen2);
ub1 = max(sen2);

if lb1 == -1000
    ddd = ['sensivity analysis b1: ','-inf' , '<b1<=' , num2str(ub1)];
    disp(ddd)
elseif ub1 == 1000
     ddd = ['sensivity analysis b1: ', num2str(lb1) , '<=b1<' , 'inf'];
     disp(ddd)
elseif lb1 == -1000 && ub1 == 1000
     ddd = ['sensivity analysis b1: ', '-inf' , '<b1<' , 'inf'];
     disp(ddd)
else
   ddd= ['sensivity analysis b1: ', num2str(lb1) , '<=b1<=' , num2str(ub1)];
    disp(ddd)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b2
for i = 0 : 1 : 200000
    ba = b;
    ba(2) = (i/100) - 1000;
    sen3 = B\ba;
    if sen3(1)>=0
        if sen3(2)>=0
            if sen3(3)>=0
                sen4(i+1,1) = (i/100) - 1000;
            else
                sen4(i+1,1) = NaN;
            end
        else
            sen4(i+1,1) = NaN;
        end
    else
        sen4(i+1,1) = NaN;
    end  
end

lb2 = min(sen4);
ub2 = max(sen4);

if lb2 == -1000
    ddd2 = ['sensivity analysis b2: ','-inf' , '<b2<=' , num2str(ub2)];
    disp(ddd2)
elseif ub2 == 1000
     ddd2 = ['sensivity analysis b2: ', num2str(lb2) , '<=b2<' , 'inf'];
     disp(ddd2)
elseif lb2 == -1000 && ub2 == 1000
     ddd2 = ['sensivity analysis b2: ', '-inf' , '<b2<' , 'inf'];
     disp(ddd2)
else
   ddd2= ['sensivity analysis b2: ', num2str(lb2) , '<=b2<=' , num2str(ub2)];
    disp(ddd2)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b3
for i = 0 : 1 : 200000
    ba = b;
    ba(3) = (i/100) - 1000;
    sen5 = B\ba;
    if sen5(1)>=0
        if sen5(2)>=0
            if sen5(3)>=0
                sen6(i+1,1) = (i/100) - 1000;
            else
                sen6(i+1,1) = NaN;
            end
        else
            sen6(i+1,1) = NaN;
        end
    else
        sen6(i+1,1) = NaN;
    end
end

lb3 = min(sen6);
ub3 = max(sen6);

if lb3 == -1000
    ddd3 = ['sensivity analysis b3: ','-inf' , '<b3<=' , num2str(ub3)];
    disp(ddd3)
elseif ub3 == 1000
     ddd3 = ['sensivity analysis b3: ', num2str(lb3) , '<=b3<' , 'inf'];
     disp(ddd3)
elseif lb3 == -1000 && ub3 == 1000
     ddd3 = ['sensivity analysis b3: ', '-inf' , '<b3<' , 'inf'];
     disp(ddd3)
else
   ddd3= ['sensivity analysis b3: ', num2str(lb3) , '<=b3<=' , num2str(ub3)];
    disp(ddd3)
end
end








if S_b(1)==4
    
    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b1

for i = 0 : 1 : 200000
    ba = b;
    ba(1) = (i/100) - 1000;
    sen1 = B\ba;
    if sen1(1)>=0
        if sen1(2)>=0
            if sen1(3)>=0
                if sen1(4)>=0
                sen2(i+1,1) = (i/100) - 1000;
            else
                sen2(i+1,1) = NaN;
            end
        else
            sen2(i+1,1) = NaN;
        end
    else
        sen2(i+1,1) = NaN;
        end 
    else
        sen2(i+1,1) = NaN;
    end 
end

lb1 = min(sen2);
ub1 = max(sen2);

if lb1 == -1000
    ddd = ['sensivity analysis b1: ','-inf' , '<b1<=' , num2str(ub1)];
    disp(ddd)
elseif ub1 == 1000
     ddd = ['sensivity analysis b1: ', num2str(lb1) , '<=b1<' , 'inf'];
     disp(ddd)
elseif lb1 == -1000 && ub1 == 1000
     ddd = ['sensivity analysis b1: ', '-inf' , '<b1<' , 'inf'];
     disp(ddd)
else
   ddd= ['sensivity analysis b1: ', num2str(lb1) , '<=b1<=' , num2str(ub1)];
    disp(ddd)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b2
for i = 0 : 1 : 200000
    ba = b;
    ba(2) = (i/100) - 1000;
    sen3 = B\ba;
    if sen3(1)>=0
        if sen3(2)>=0
            if sen3(3)>=0
                if sen3(4)>=0
                sen4(i+1,1) = (i/100) - 1000;
            else
                sen4(i+1,1) = NaN;
            end
        else
            sen4(i+1,1) = NaN;
        end
    else
        sen4(i+1,1) = NaN;
        end  
    else
        sen4(i+1,1) = NaN;
    end  
end

lb2 = min(sen4);
ub2 = max(sen4);

if lb2 == -1000
    ddd2 = ['sensivity analysis b2: ','-inf' , '<b2<=' , num2str(ub2)];
    disp(ddd2)
elseif ub2 == 1000
     ddd2 = ['sensivity analysis b2: ', num2str(lb2) , '<=b2<' , 'inf'];
     disp(ddd2)
elseif lb2 == -1000 && ub2 == 1000
     ddd2 = ['sensivity analysis b2: ', '-inf' , '<b2<' , 'inf'];
     disp(ddd2)
else
   ddd2= ['sensivity analysis b2: ', num2str(lb2) , '<=b2<=' , num2str(ub2)];
    disp(ddd2)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b3
for i = 0 : 1 : 200000
    ba = b;
    ba(3) = (i/100) - 1000;
    sen5 = B\ba;
    if sen5(1)>=0
        if sen5(2)>=0
            if sen5(3)>=0
                if sen5(4)>=0
                sen6(i+1,1) = (i/100) - 1000;
            else
                sen6(i+1,1) = NaN;
            end
        else
            sen6(i+1,1) = NaN;
        end
    else
        sen6(i+1,1) = NaN;
        end
    else
        sen6(i+1,1) = NaN;
    end
end

lb3 = min(sen6);
ub3 = max(sen6);

if lb3 == -1000
    ddd3 = ['sensivity analysis b3: ','-inf' , '<b3<=' , num2str(ub3)];
    disp(ddd3)
elseif ub3 == 1000
     ddd3 = ['sensivity analysis b3: ', num2str(lb3) , '<=b3<' , 'inf'];
     disp(ddd3)
elseif lb3 == -1000 && ub3 == 1000
     ddd3 = ['sensivity analysis b3: ', '-inf' , '<b3<' , 'inf'];
     disp(ddd3)
else
   ddd3= ['sensivity analysis b3: ', num2str(lb3) , '<=b3<=' , num2str(ub3)];
    disp(ddd3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b4
for i = 0 : 1 : 200000
    ba = b;
    ba(4) = (i/100) - 1000;
    sen7 = B\ba;
    if sen7(1)>=0
        if sen7(2)>=0
            if sen7(3)>=0
                if sen7(4)>=0
                sen8(i+1,1) = (i/100) - 1000;
            else
                sen8(i+1,1) = NaN;
            end
        else
            sen8(i+1,1) = NaN;
        end
    else
        sen8(i+1,1) = NaN;
        end
    else
        sen8(i+1,1) = NaN;
    end
end

lb4 = min(sen8);
ub4 = max(sen8);

if lb4 == -1000
    ddd4 = ['sensivity analysis b4: ','-inf' , '<b4<=' , num2str(ub4)];
    disp(ddd4)
elseif ub4 == 1000
     ddd4 = ['sensivity analysis b4: ', num2str(lb4) , '<=b4<' , 'inf'];
     disp(ddd4)
elseif lb4 == -1000 && ub4 == 1000
     ddd4 = ['sensivity analysis b4: ', '-inf' , '<b4<' , 'inf'];
     disp(ddd4)
else
   ddd4= ['sensivity analysis b4: ', num2str(lb4) , '<=b4<=' , num2str(ub4)];
    disp(ddd4)
end
end








if S_b(1)==5
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b1

for i = 0 : 1 : 200000
    ba = b;
    ba(1) = (i/100) - 1000;
    sen1 = B\ba;
    if sen1(1)>=0
        if sen1(2)>=0
            if sen1(3)>=0
                if sen1(4)>=0
                    if sen1(5)>=0
                
                sen2(i+1,1) = (i/100) - 1000;
            else
                sen2(i+1,1) = NaN;
            end
        else
            sen2(i+1,1) = NaN;
        end
    else
        sen2(i+1,1) = NaN;
    end 
            else
            sen2(i+1,1) = NaN;
        end
    else
        sen2(i+1,1) = NaN;
    end 
end

lb1 = min(sen2);
ub1 = max(sen2);

if lb1 == -1000
    ddd = ['sensivity analysis b1: ','-inf' , '<b1<=' , num2str(ub1)];
    disp(ddd)
elseif ub1 == 1000
     ddd = ['sensivity analysis b1: ', num2str(lb1) , '<=b1<' , 'inf'];
     disp(ddd)
elseif lb1 == -1000 && ub1 == 1000
     ddd = ['sensivity analysis b1: ', '-inf' , '<b1<' , 'inf'];
     disp(ddd)
else
   ddd= ['sensivity analysis b1: ', num2str(lb1) , '<=b1<=' , num2str(ub1)];
    disp(ddd)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b2
for i = 0 : 1 : 200000
    ba = b;
    ba(2) = (i/100) - 1000;
    sen3 = B\ba;
    if sen3(1)>=0
        if sen3(2)>=0
            if sen3(3)>=0
                if sen3(4)>=0
                    if sen3(5)>=0
                sen4(i+1,1) = (i/100) - 1000;
            else
                sen4(i+1,1) = NaN;
            end
        else
            sen4(i+1,1) = NaN;
        end
    else
        sen4(i+1,1) = NaN;
            end  
            else
            sen4(i+1,1) = NaN;
        end
    else
        sen4(i+1,1) = NaN;
    end  
end

lb2 = min(sen4);
ub2 = max(sen4);

if lb2 == -1000
    ddd2 = ['sensivity analysis b2: ','-inf' , '<b2<=' , num2str(ub2)];
    disp(ddd2)
elseif ub2 == 1000
     ddd2 = ['sensivity analysis b2: ', num2str(lb2) , '<=b2<' , 'inf'];
     disp(ddd2)
elseif lb2 == -1000 && ub2 == 1000
     ddd2 = ['sensivity analysis b2: ', '-inf' , '<b2<' , 'inf'];
     disp(ddd2)
else
   ddd2= ['sensivity analysis b2: ', num2str(lb2) , '<=b2<=' , num2str(ub2)];
    disp(ddd2)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b3
for i = 0 : 1 : 200000
    ba = b;
    ba(3) = (i/100) - 1000;
    sen5 = B\ba;
    if sen5(1)>=0
        if sen5(2)>=0
            if sen5(3)>=0
                if sen5(4)>=0
                    if sen5(5)>=0
                sen6(i+1,1) = (i/100) - 1000;
            else
                sen6(i+1,1) = NaN;
            end
        else
            sen6(i+1,1) = NaN;
        end
    else
        sen6(i+1,1) = NaN;
    end
            else
            sen6(i+1,1) = NaN;
        end
    else
        sen6(i+1,1) = NaN;
    end
end

lb3 = min(sen6);
ub3 = max(sen6);

if lb3 == -1000
    ddd3 = ['sensivity analysis b3: ','-inf' , '<b3<=' , num2str(ub3)];
    disp(ddd3)
elseif ub3 == 1000
     ddd3 = ['sensivity analysis b3: ', num2str(lb3) , '<=b3<' , 'inf'];
     disp(ddd3)
elseif lb3 == -1000 && ub3 == 1000
     ddd3 = ['sensivity analysis b3: ', '-inf' , '<b3<' , 'inf'];
     disp(ddd3)
else
   ddd3= ['sensivity analysis b3: ', num2str(lb3) , '<=b3<=' , num2str(ub3)];
    disp(ddd3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b4
for i = 0 : 1 : 200000
    ba = b;
    ba(4) = (i/100) - 1000;
    sen7 = B\ba;
    if sen7(1)>=0
        if sen7(2)>=0
            if sen7(3)>=0
                if sen7(4)>=0
                    if sen7(5)>=0
                sen8(i+1,1) = (i/100) - 1000;
            else
                sen8(i+1,1) = NaN;
            end
            
        else
            sen8(i+1,1) = NaN;
        end
    else
        sen8(i+1,1) = NaN;
    end
            else
            sen8(i+1,1) = NaN;
        end
    else
        sen8(i+1,1) = NaN;
    end
end

lb4 = min(sen8);
ub4 = max(sen8);

if lb4 == -1000
    ddd4 = ['sensivity analysis b4: ','-inf' , '<b4<=' , num2str(ub4)];
    disp(ddd4)
elseif ub4 == 1000
     ddd4 = ['sensivity analysis b4: ', num2str(lb4) , '<=b4<' , 'inf'];
     disp(ddd4)
elseif lb4 == -1000 && ub4 == 1000
     ddd4 = ['sensivity analysis b4: ', '-inf' , '<b4<' , 'inf'];
     disp(ddd4)
else
   ddd4= ['sensivity analysis b4: ', num2str(lb4) , '<=b4<=' , num2str(ub4)];
    disp(ddd4)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensivity analysis b5
for i = 0 : 1 : 200000
    ba = b;
    ba(5) = (i/100) - 1000;
    sen9 = B\ba;
    if sen9(1)>=0
        if sen9(2)>=0
            if sen9(3)>=0
                if sen9(4)>=0
                    if sen9(5)>=0
                sen10(i+1,1) = (i/100) - 1000;
            else
                sen10(i+1,1) = NaN;
            end
        else
            sen10(i+1,1) = NaN;
        end
    else
        sen1(i+1,1) = NaN;
            end
            else
            sen10(i+1,1) = NaN;
        end
    else
        sen1(i+1,1) = NaN;
    end
end

lb5 = min(sen10);
ub5 = max(sen10);

if lb5 == -1000
    ddd5 = ['sensivity analysis b5: ','-inf' , '<b5<=' , num2str(ub5)];
    disp(ddd5)
elseif ub5 == 1000
     ddd5 = ['sensivity analysis b5: ', num2str(lb5) , '<=b5<' , 'inf'];
     disp(ddd5)
elseif lb5 == -1000 && ub5 == 1000
     ddd5 = ['sensivity analysis b5: ', '-inf' , '<b5<' , 'inf'];
     disp(ddd5)
else
   ddd5= ['sensivity analysis b5: ', num2str(lb5) , '<=b5<=' , num2str(ub5)];
    disp(ddd5)
end
end
end
%%
%jadval nahaii doal
%%moarefi moteghayer
Dual=input('Show Dual Final Table? yes(1)/no(2)');
if Dual==1
if MaxMin==1
    T2=table();
    for i=1:SA1(1)+SA2(1)+2*SA3(1)+2*SA(2)+1
        V3=cell(1,i);
    end
    for i=1:SA1(1)
        V3{1,i}=strcat('Y',num2str(i));
    end
    for i=1+SA1(1):SA1(1)+SA3(1)
        V3{1,i}=strcat('Yp',num2str(i));
        V3{1,i+1}=strcat('Yn',num2str(i));
    end
    for i=1+SA1(1)+SA3(1):SA1(1)+SA3(1)+SA2(1)
        V3{1,i+SA3(1)}=strcat('Y^',num2str(i));
    end
    for i=1+SA1(1)+SA2(1)+2*SA3(1):SA1(1)+SA2(1)+2*SA3(1)+SA(2)
        V3{1,i}=strcat('S',num2str(i-SA1(1)-SA2(1)-2*SA3(1)));
        V3{1,i+SA(2)}=strcat('R',num2str(i-SA1(1)-SA2(1)-2*SA3(1)));
    end
    V3{1,end}=strcat('RHS');
    F3=find(CT(1,1:SA(2)+SA1(1)+SA2(1)+SA3(1)));
    F4=find(CT==0);
    for i=1:size(F4,2)
        F5(1,i)=find(P2(:,F4(1,i))==1);
    end
    for i=1:size(F3,2)
        V4=cell(i,1);
    end
    V4{1,1}=strcat('W');
    for i=1:size(F3,2)
        for j=1:SA(2)
            if F3(1,i)==j
                V4{i+1,1}=strcat('S',num2str(j));
                U1=CT(1,j);
            end
        end
        for j=1+SA(2):SA(2)+SA1(1)
            if F3(1,i)==j
                V4{i+1,1}=strcat('Y',num2str(j-SA(2)));
            end
        end
        for j=1+SA(2)+SA1(1):SA(2)+SA1(1)+SA3(1)
           if F3(1,i)==j
               U(1,j-SA(2)-SA1(1))=CT(1,j)-100;
               if U(1,j-SA(2)-SA1(1)) >= 0
                   V4{i+1,1}=strcat('Yp',num2str(j-SA(2)));
               else
                   V4{i+1,1}=strcat('Yn',num2str(j-SA(2)));
               end
           end
        end
           for j=1+SA(2)+SA1(1)+SA3(1):SA(2)+SA1(1)+SA3(1)+SA2(1)
               if F3(1,i)==j
                   V4{i+1,1}=strcat('Y^',num2str(j-SA(2)-SA3(1)));
               end
           end
    end
    for i=1:size(F3,2)
        if F3(1,i)>SA(2)+SA1(1) && F3(1,i)<=SA(2)+SA1(1)+SA3(1)+SA2(1)
            b22(1,i)=abs(CT(1,F3(1,i))-100);
        else
            b22(1,i)=CT(1,F3(1,i));
        end
    end
    %%martesi jadval
    b2=transpose(b22);
    for i=1:size(F4,2)
        if F4(1,i)<=SA(2)
            CT2(1,F4(1,i)+SA1(1)+SA2(1)+2*SA3(1))=XB(F5(1,i),1);
            CT2(1,F4(1,i)+SA1(1)+SA2(1)+2*SA3(1)+SA(2))=100-XB(F5(1,i),1);
        elseif F4(1,i)>SA(2) && F4(1,i)<=SA(2)+SA1(1)
            CT2(1,F4(1,i)-SA(2))=XB(F5(1,i),1);
        elseif F4(1,i)>SA(2)+SA1(1) && F4(1,i)<=SA(2)+SA1(1)+SA(3)
            disp('infeasible')
        elseif F4(1,i)>SA(2)+SA1(1)+SA(3)+SA(2) && F4(1,i)<=SA(2)+SA1(1)+SA(3)+2*SA(2)
            CT2(1,F4(1,i)-SA(2)-SA2(1))=XB(F5(1,i),1);
        end
    end
    for i=SA(1)+2*SA3(1)+1:SA(1)+2*SA3(1)+SA(2)
        if CT2(1,i)==0
            CT2(1,i+SA(2))=100;
        end
    end
    F3=find(CT(1,1:end-SA2(1)));
    for i=1:size(F3,2)
        for j=1:SA(2)
            if F3(1,i)==j
                P3(i,j+SA1(1)+2*SA3(1)+SA2(1))=1;
                P3(i,j+SA1(1)+2*SA3(1)+SA(2)+SA2(1))=-1;
            else
                P3(i,j+SA1(1)+2*SA3(1))=0;
            end
        end
        for j=1+SA(2):SA(2)+SA1(1)
            if F3(1,i)==j
                P3(i,j-SA(2))=1;
            else
                P3(i,j-SA(2))=0;
            end
        end
        for j=1+SA(2)+SA1(1):SA(2)+SA1(1)+SA3(1)
           if F3(1,i)==j
               P3(i,j-SA(2))=1; 
               P3(i,j-SA(2)+1)=-1;
           else
               P3(i,j-SA(2))=0;
               P3(i,j-SA(2)+1)=0;
           end
        end
        for j=1+SA(2)+SA1(1)+SA3(1):SA(2)+SA1(1)+SA3(1)+SA2(1)
            if F3(1,i)==j
                P3(i,j-SA(2)+SA3(1))=1; 
            else
                P3(i,j-SA(2)+SA3(1))=0;
            end
        end
    end
    %for i=1:size(F4,2)
        P2(:,F4)=[];
       if  ~isempty(A2)
        for i=1:size(b2,1)
            h(1,i)=P2(1,i);
            P2=h;
        end
        end 
    %end
    for i=1:size(XB,1)
        if XB(i,1)==0
        else
        F6=find(CT2==XB(i,1));
        P3(:,F6)=-transpose(P2(i,:));
        P3(:,F6+SA(2))=transpose(P2(i,:));
        end
    end
    
for i=1:size(F3,2)
    for j=1+SA(2)+SA1(1):SA(2)+SA1(1)+SA3(1)
           if F3(1,i)==j
               if U(1,j-SA(2)-SA1(1)) < 0
                  P3(i,:)=-P3(i,:);
               end
           end
    end
end
CAT2=[CT2 -Z;P3 b2];
for i=1:size(CAT2,2)
    T2{:,i}=CAT2(:,i);
end
T2.Properties.VariableNames=V3;
T2.Properties.RowNames=V4;
disp(T2)
disp('Dual Final Table')
%%braye min
%%

else
    T2=table();
    for i=1:SA1(1)+SA2(1)+2*SA3(1)+SA(2)+1
        V3=cell(1,i);
    end
    for i=1:SA1(1)
        V3{1,i}=strcat('Y^',num2str(i));
    end
    for i=1+SA1(1):SA1(1)+SA3(1)
        V3{1,i}=strcat('Yp',num2str(i));
        V3{1,i+1}=strcat('Yn',num2str(i));
    end
    for i=1+SA1(1)+SA3(1):SA1(1)+SA3(1)+SA2(1)
        V3{1,i+SA3(1)}=strcat('Y',num2str(i));
    end
    for i=1+SA1(1)+SA2(1)+2*SA3(1):SA1(1)+SA2(1)+2*SA3(1)+SA(2)
        V3{1,i}=strcat('S',num2str(i-SA1(1)-SA2(1)-2*SA3(1)));
    end
    V3{1,end}=strcat('RHS');
    F3=find(CT(1,1:SA(2)+SA1(1)+SA2(1)+SA3(1)));
    F4=find(CT==0);
    for i=1:size(F4,2)
        F5(1,i)=find(P2(:,F4(1,i))==1);
    end
    for i=1:size(F3,2)
        V4=cell(i+1,1);
    end
    V4{1,1}=strcat('Z');
    for i=1:size(F3,2)
        for j=1:SA(2)
            if F3(1,i)==j
                V4{i+1,1}=strcat('S',num2str(j));
                U1=CT(1,j);
            end
        end
        for j=1+SA(2):SA(2)+SA1(1)
            if F3(1,i)==j
                V4{i+1,1}=strcat('Y^',num2str(j-SA(2)));
            end
        end
        for j=1+SA(2)+SA1(1):SA(2)+SA1(1)+SA3(1)
           if F3(1,i)==j
               U(1,j-SA(2)-SA1(1))=CT(1,j)-100;
               if U(1,j-SA(2)-SA1(1)) >= 0
                   V4{i+1,1}=strcat('Yp',num2str(j-SA(2)));
               else
                   V4{i+1,1}=strcat('Yn',num2str(j-SA(2)));
               end
           end
        end
           for j=1+SA(2)+SA1(1)+SA3(1):SA(2)+SA1(1)+SA3(1)+SA2(1)
               if F3(1,i)==j
                   V4{i+1,1}=strcat('Y',num2str(j-SA(2)-SA3(1)));
               end
           end
    end
    for i=1:size(F3,2)
        if F3(1,i)>SA(2)+SA1(1) && F3(1,i)<=SA(2)+SA1(1)+SA3(1)+SA2(1)
            b22(1,i)=abs(CT(1,F3(1,i))-100);
        else
            b22(1,i)=CT(1,F3(1,i));
        end
    end
    %%martesi jadval
    b2=transpose(b22);
    for i=1:size(F4,2)
        if F4(1,i)<=SA(2)
            CT2(1,F4(1,i)+SA1(1)+SA2(1)+SA3(1))=XB(F5(1,i),1);
            
        elseif F4(1,i)>SA(2) && F4(1,i)<=SA(2)+SA1(1)
            CT2(1,F4(1,i)-SA(2))=XB(F5(1,i),1);
        elseif F4(1,i)>SA(2)+SA1(1) && F4(1,i)<=SA(2)+SA1(1)+SA(3)
            disp('infeasible')
        elseif F4(1,i)>SA(2)+SA1(1)+SA(3)+SA(2) && F4(1,i)<=SA(2)+SA1(1)+SA(3)+2*SA(2)
            CT2(1,F4(1,i)-SA(2)-SA2(1))=XB(F5(1,i),1);
        end
    end

    F3=find(CT(1,1:end-SA2(1)));
    for i=1:size(F3,2)
        for j=1:SA(2)
            if F3(1,i)==j
                P3(i,j+SA1(1)+SA3(1)+SA2(1))=1;
                
            else
                P3(i,j+SA1(1)+SA3(1))=0;
            end
        end
        for j=1+SA(2):SA(2)+SA1(1)
            if F3(1,i)==j
                P3(i,j-SA(2))=1;
            else
                P3(i,j-SA(2))=0;
            end
        end
        for j=1+SA(2)+SA1(1):SA(2)+SA1(1)+SA3(1)
           if F3(1,i)==j
               P3(i,j-SA(2))=1; 
               P3(i,j-SA(2)+1)=-1;
           else
               P3(i,j-SA(2))=0;
               P3(i,j-SA(2)+1)=0;
           end
        end
        for j=1+SA(2)+SA1(1)+SA3(1):SA(2)+SA1(1)+SA3(1)+SA2(1)
            if F3(1,i)==j
                P3(i,j-SA(2)+SA3(1))=1; 
            else
                P3(i,j-SA(2)+SA3(1))=0;
            end
        end
    end
    %for i=1:size(F4,2)
        P2(:,F4)=[];
        
        if  ~isempty(A2)
        for i=1:size(b2,1)
            h(1,i)=P2(1,i);
            P2=h;
        end
        end 
    %end
    for i=1:size(XB,1)
        if XB(i,1)==0
        else
        F6=find(CT2==XB(i,1));
        P3(:,F6)=-transpose(P2(i,:));
        
        end
    end
   
CAT2=[CT2 Z;P3 b2];
for i=1:size(CAT2,2)
    T2{:,i}=CAT2(:,i);
end
T2.Properties.VariableNames=V3;
T2.Properties.RowNames=V4;
disp(T2)
disp('Dual Final Table')    

end
end
