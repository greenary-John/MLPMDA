%based on dataset:AMVML
eta=0.01;
dis_qz1=0.4;
dis_qz2=0.1;
dis_z=1;
mi_qz1=1;
mi_qz2=0.7;
mi_z=10;
top_n=11;
%第一层视图
dis_qz1_v1=0.1;
dis_qz2_v1=0.0;
dis_z_v1=1;
mi_qz1_v1=1;
mi_qz2_v1=0;
mi_z_v1=30;
top_n_v1=12;
%第二层视图
dis_qz1_v2=1;
dis_qz2_v2=1;
dis_z_v2=16;
mi_qz1_v2=0.1;
mi_qz2_v2=0.6;
mi_z_v2=16;
top_n_v2=25;

MD=load("dataset/dataset_AMVML/RDmat2.txt");

RS_source_go=load("dataset/dataset_AMVML/mi_go_sim.csv");
RS_source_fun=load("dataset/dataset_AMVML/mi_fun_sim.csv");

RS_source=(RS_source_go*0.089+RS_source_fun*0.911); 

y = MD;
y_train = y;
[m,n]=size(y_train);
[rna_guass,dis_guass]=gaussiansimilarity(y_train,m,n);
DS=(DS_source*dis_qz1+dis_guass*dis_qz2)*dis_z;
RS=(RS_source*mi_qz1+rna_guass*mi_qz2)*mi_z;
MD_update=preprocess_WKNKN(MD,RS,DS,top_n,1);
yy=y;
yy(yy==0)=-1;
H1=[RS,MD_update;MD_update',DS]; 
C1 = (1/eta*eye(size(H1'*H1)) + ...
     H1' * H1)\H1'* H1; 
P_temp_1 = H1 * C1;
P1 = P_temp_1(size(RS,1)+1:end,1:size(MD',2))';  

[m,n]=size(y_train);
[rna_guass2,dis_guass2]=gaussiansimilarity(P1,m,n);
DS2=(DS_source*dis_qz1+dis_guass2*dis_qz2)*dis_z;
RS2=(RS_source*mi_qz1+rna_guass2*mi_qz2)*mi_z;
P1_add=max(P1,MD);
H2=[RS2,P1_add;P1_add',DS2];
C2 = (1/eta*eye(size(H2'*H2)) + ...
        H2' * H2)\H2'* H2; 
P_temp_2 = H2 * C2;
P2 = P_temp_2(size(RS,1)+1:end,1:size(MD',2))';  
        
[m,n]=size(P2);
[rna_guass_v1,dis_guass_v1]=gaussiansimilarity(P2,m,n);
DS_v1=(DS_source*dis_qz1_v1+dis_guass_v1*dis_qz2_v1)*dis_z_v1;
RS_v1=(RS_source_fun*mi_qz1_v1+rna_guass_v1*mi_qz2_v1)*mi_z_v1;
train_v1=preprocess_WKNKN(P2,RS_v1,DS_v1,top_n_v1,1);
       
H_v1=[RS_v1,train_v1;train_v1',DS_v1]; 

C_v1 =(1/eta*eye(size(H_v1'*H_v1)) + ...
                    H_v1' * H_v1)\H_v1'* H_v1; 
P_temp_v1 = H_v1 * C_v1;
P_v1 = P_temp_v1(size(RS,1)+1:end,1:size(MD',2))';  

[m,n]=size(P2);
[rna_guass_v2,dis_guass_v2]=gaussiansimilarity(P2,m,n);
DS_v2=(DS_source*dis_qz1_v2+dis_guass_v2*dis_qz2_v2)*dis_z_v2;
RS_v2=(RS_source_go*mi_qz1_v2+rna_guass_v2*mi_qz2_v2)*mi_z_v2;
train_v2=preprocess_WKNKN(P2,RS_v2,DS_v2,top_n_v2,1);

H_v2=[RS_v2,train_v2;train_v2',DS_v2]; 

C_v2 =(1/eta*eye(size(H_v2'*H_v2)) + ...
   H_v2' * H_v2)\H_v2'* H_v2; 
P_temp_v2 = H_v2 * C_v2;
P_v2 =P_temp_v2(size(RS,1)+1:end,1:size(MD',2))';  

P_v_all=max(P_v1,P_v2);
MD_out=max(P_v_all,P2);

[Breast,  index_Breast] = sort(-(MD(:,5)-1).*MD_out(:,5),'descend');
C_breast=join(table(index_Breast),miRNANames,'keys',1);
string_brest = strrep(table2array(C_breast(1:50,2)),'"','')


[HCC,  index_HCC] = sort(-(MD(:,7)-1).*MD_out(:,7),'descend');
C_HCC=join(table(index_HCC),miRNANames,'keys',1);
string_HCC = strrep(table2array(C_HCC(1:50,2)),'"','');

[Colon,  index_Colon] = sort(-(MD(:,1)-1).*MD_out(:,1),'descend');
C_Colon=join(table(index_Colon),miRNANames,'keys',1);
string_Colon = strrep(table2array(C_Colon(1:50,2)),'"','');