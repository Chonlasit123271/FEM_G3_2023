close all;clear all;clc
tic
% % ------------------------------------START_INPUT-SECTION--------------------------------------
nu=0.2;                         %Poisson ratio
E0=22.36*10^9;                         %Young's Modulus
SigmaMax=20*10^6;                       %Compressive Strength of Concrete             
%CHANGE THE MESH SIZE HERE
inputfile    =  ('80nodes_1.msh');
[Coord,Quads,ELEMCon] = autoread(inputfile);
Plot_Element01(ELEMCon,Coord)
exportgraphics(gcf, ['Node_and_Connectivity_2.jpg']); % Save the figure as a JPEG image
%%%%__________Define support__________%%%%
%%____matrix form [Node_no ux uy uz]------------sp --> support matrix
sp = [1 1 1 1;
      2 1 1 1;
      3 1 1 1;
      4 1 1 1;
      16 1 1 1;
      17 1 1 1;
      19 1 1 1;
      20 1 1 1];   %% sp --> support matrix
nsp = length(sp(:,1)) ; %% number of support nodes

%%%%__________Load conditions__________%%%%
%%____matrix form [Node_no Fx Fy Fz]------------load --> load matrix
load=[10 696000 0 0;
      18 696000 0 0];
nload = length(load(:,1)) ; %% number of support nodes.
% NINC ----> No. of an increment
NINC=15;
% % ------------------------------------FINISH_INPUT-SECTION--------------------------------------
% % ----------------------------------------------------------------------------------------------

% % ----------------------------------------------------------------------------------------------
% % ---------------------------------------Processing_Part----------------------------------------
% % ----------------------------------------------------------------------------------------------
disp("-------------------------------------------------------------------------------------------")
disp("-------------------------------------Processing--------------------------------------------")
disp("-------------------------------------------------------------------------------------------")
savejpg(Coord,Quads,'Before')
NE=length(ELEMCon(:,1));                %Number of elements
nNode=length(Coord(:,1));             %Number of nodes
nDof=nNode*3;                            %Number of degrees of freedom%% 
%%%%__________Node and Element container__________%%%%
for p=1:nNode
    NODE(p).X=Coord(p,1); NODE(p).Y=Coord(p,2); NODE(p).Z=Coord(p,3);
end
for eNo=1:NE
    ELEMENT(eNo).con=ELEMCon(eNo,:);
end

%%%%__________Shape functions__________%%%%
syms xi eta zeta
N(1)=(1/8)*(1-xi)*(1-eta)*(1-zeta);
N(2)=(1/8)*(1+xi)*(1-eta)*(1-zeta);
N(3)=(1/8)*(1+xi)*(1+eta)*(1-zeta);
N(4)=(1/8)*(1-xi)*(1+eta)*(1-zeta);
N(5)=(1/8)*(1-xi)*(1-eta)*(1+zeta);
N(6)=(1/8)*(1+xi)*(1-eta)*(1+zeta);
N(7)=(1/8)*(1+xi)*(1+eta)*(1+zeta);
N(8)=(1/8)*(1-xi)*(1+eta)*(1+zeta);
nN=length(N);

%%%%__________Prepare array container of Young's modulus and displacement__________%%%%
for eNo=1:NE
    ELEMENT(eNo).Ex=E0;
    ELEMENT(eNo).Ey=E0;
    ELEMENT(eNo).Ez=E0;
    ELEMENT(eNo).Exy=E0;
    ELEMENT(eNo).Eyz=E0;
    ELEMENT(eNo).Exz=E0;
end
uDisp=zeros(nDof,1); 

%%%%__________Prepare array container for an incremental__________%%%
uDel_inc=zeros(NINC,nload);
uLoad_inc=zeros(NINC,nload);
% ------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------Each element--------------------------------------------------
% ------------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------------
for nin = 1:NINC
    disp("START NINC:    " +nin)
    for eNo=1:NE
        %Element Coordinates
        for i=1:nN
            xp(i)=Coord(ELEMCon(eNo,i),1);
            yp(i)=Coord(ELEMCon(eNo,i),2);
            zp(i)=Coord(ELEMCon(eNo,i),3);
            ELEMENT(eNo).X(i)=NODE(ELEMENT(eNo).con(i)).X;
            ELEMENT(eNo).Y(i)=NODE(ELEMENT(eNo).con(i)).Y;
            ELEMENT(eNo).Z(i)=NODE(ELEMENT(eNo).con(i)).Z;
        end
    
        %Coordinate Mapping
        %Mapping
        x=0;y=0;z=0;
        for i=1:nN
            x=x+xp(i)*N(i);
            y=y+yp(i)*N(i);
            z=z+zp(i)*N(i);
        end
    
        %3D Jacobian matrix 
        J=[ diff(x,xi) diff(y,xi) diff(z,xi);
            diff(x,eta) diff(y,eta) diff(z,eta);
            diff(x,zeta) diff(y,zeta) diff(z,zeta)];
        
        %Differentiating the N matrix to get doMat, dN
        for i=1:nN
            dN(:,i)=inv(J)*[diff(N(i),xi);diff(N(i),eta); diff(N(i),zeta)];
        end
        
        %Arranging the terms to get the B matrix
        for i=1:nN
            j=3*i-2;k=3*i;
            B(1:6,j:k)=[dN(1,i) , 0         , 0 ;
                        0       , dN(2,i)   , 0 ;
                        0       , 0         , dN(3,i);
                        dN(2,i) , dN(1,i)   , 0 ;
                        0       , dN(3,i)   , dN(2,i) ;
                        dN(3,i) , 0         , dN(1,i)];
        end                 
       
        ELEMENT(eNo).B=B;
        D=(1/(1-nu^2))*[ELEMENT(eNo).Ex    ELEMENT(eNo).Ex*nu  ELEMENT(eNo).Ex*nu  0           0           0;
                                   ELEMENT(eNo).Ey*nu   ELEMENT(eNo).Ey*1   ELEMENT(eNo).Ey*nu  0           0           0;
                                   ELEMENT(eNo).Ez*nu   ELEMENT(eNo).Ez*nu  ELEMENT(eNo).Ez*1   0           0           0;
                                   0    0   0   ELEMENT(eNo).Exy*(1-nu)/2    0           0;
                                   0    0   0   0           ELEMENT(eNo).Eyz*(1-nu)/2    0;
                                   0    0   0   0           0           ELEMENT(eNo).Exz*(1-nu)/2];
        %Integrand is
        f= B'*D*B*det(J);
        [frows, fcols] = size(f);
        
        % Gauss points
        [xg,wg]=lgwt(2,-1,1);
        
        % Stiffness matrix by guassian integration points
        Ke=zeros(frows, fcols);
        Ke=Ke+wg(1)*wg(1)*wg(1)*subs(f,[xi,eta,zeta],[xg(1),xg(1),xg(1)])+wg(2)*wg(1)*wg(1)*subs(f,[xi,eta,zeta],[xg(2),xg(1),xg(1)])+wg(1)*wg(2)*wg(1)*subs(f,[xi,eta,zeta],[xg(1),xg(2),xg(1)])+wg(2)*wg(2)*wg(1)*subs(f,[xi,eta,zeta],[xg(2),xg(2),xg(1)])+wg(1)*wg(1)*wg(2)*subs(f,[xi,eta,zeta],[xg(1),xg(1),xg(2)])+wg(2)*wg(1)*wg(2)*subs(f,[xi,eta,zeta],[xg(2),xg(1),xg(2)])+wg(1)*wg(2)*wg(2)*subs(f,[xi,eta,zeta],[xg(1),xg(2),xg(2)])+wg(2)*wg(2)*wg(2)*subs(f,[xi,eta,zeta],[xg(2),xg(2),xg(2)]);
        ELEMENT(eNo).stiffness=Ke;
    end
    disp(['      Completed Jacobian Matrix: J      ', num2str(toc), ' sec'])
    disp(['      Completed B Matrix: B      ', num2str(toc), ' sec'])
    disp(['      Completed Element Stiffness Matrix: Ke      ', num2str(toc), ' sec'])
    %-------------------------------------------------------------------------%
    %                  Global Stiffness Matrix Calculation                    %
    %-------------------------------------------------------------------------%
    % This part calculates the global stiffness matrix. Basically; for each
    % element it takes 12x12 part from the element stiffness matrix and puts to
    % the correct spot on the global stiffness matrix. This process loops until
    % all elements all parts placed in to the global stiffness matrix.
    KG = zeros(nNode*3,nNode*3);
    for eNo=1:NE
        kElem=  ELEMENT(eNo).stiffness;
        for j=1:nN
            for i=1:nN
                n = ELEMCon(eNo,i);
                m = ELEMCon(eNo,j);
                KG(3*n-2,3*m-2) = KG(3*n-2,3*m-2)+kElem(3*i-2,3*j-2);
                KG(3*n-2,3*m-1) = KG(3*n-2,3*m-1)+kElem(3*i-2,3*j-1);
                KG(3*n-2,3*m) = KG(3*n-2,3*m)+kElem(3*i-2,3*j);
                KG(3*n-1,3*m-2) = KG(3*n-1,3*m-2)+kElem(3*i-1,3*j-2);
                KG(3*n-1,3*m-1) = KG(3*n-1,3*m-1)+kElem(3*i-1,3*j-1);
                KG(3*n-1,3*m) = KG(3*n-1,3*m)+kElem(3*i-1,3*j);
                KG(3*n,3*m-2) = KG(3*n,3*m-2)+kElem(3*i,3*j-2);
                KG(3*n,3*m-1) = KG(3*n,3*m-1)+kElem(3*i,3*j-1);
                KG(3*n,3*m) = KG(3*n,3*m)+kElem(3*i,3*j);
            end
        end
    end
    disp(['      Completed Global Stiffness Matrix: KG      ', num2str(toc), ' sec'])
    KG= sparse(KG);
    %-------------------------------------------------------------------------%
    %           Apply Support Conditions to the Global Stiffnes Matrix        %
    %-------------------------------------------------------------------------%
    % This part makes zeros all columns and rows where the supports are except
    % the diagonal element of the matrix. Diagonal element set to 1. I choose
    % this method because with this way sort of displacement evaulated are not
    % changes. And later it is easier to use evaluated values. Only negative
    % side of this approach is we have to be careful not to put force where the
    % support is fixed.
    for i=1:nsp
        n = sp(i,1);
        if (sp(i,2) == 1)  %%For ux
            KG(3*n-2,:) = 0; 
            KG(:,3*n-2) = 0;
            KG(3*n-2,3*n-2) = 1;
        end
        if (sp(i,3) == 1)  %%For uy
            KG(3*n-1,:) = 0;
            KG(:,3*n-1) = 0;
            KG(3*n-1,3*n-1) = 1;
        end
        if (sp(i,4) == 1)  %%For uy
            KG(3*n,:) = 0;
            KG(:,3*n) = 0;
            KG(3*n,3*n) = 1;
        end
    end
    disp(['      Applied Supports Condition to Global Stiffness Matrix: KG      ', num2str(toc), ' sec'])
    %-------------------------------------------------------------------------%
    %                       Load Vector Computation                           %
    %-------------------------------------------------------------------------%
    % In this part load vector created. If there is a load vector get the value
    % from load matrix. If not not load value set to zero.
    fvec = zeros(nDof,1);
    for i=1:nload
        n = load(i,1);
        fvec(3*n-2) = load(i,2)/NINC;       %%%%____Fx
        fvec(3*n-1) = load(i,3)/NINC;       %%%%____Fy
        fvec(3*n) = load(i,4)/NINC;         %%%%____Fz
    end
    disp(['      Created Force Matrix: fvec      ', num2str(toc), ' sec'])
    %-------------------------------------------------------------------------%
    %                       Displacement Computation                          %
    %-------------------------------------------------------------------------%
    uDisp =uDisp+ KG\fvec;
    for i=1:(nload)          % Display all displacements on the free end
        uDel_inc(nin,i)=uDisp(load(i,1)*3-2);
        if(nin==1)
            uLoad_inc(nin,i)=fvec(load(i,1)*3-2);
        else
            uLoad_inc(nin,i)=uLoad_inc(nin-1,i)+fvec(load(i,1)*3-2);
        end
    end
    disp(['      Completed Displacement Computation: uDisp      ', num2str(toc), ' sec'])
    %Calculation of stresses
    
    %Nodal displacements
    for i=1:nNode
        NODE(i).u=[uDisp(3*i-2); uDisp(3*i-1); uDisp(3*i)];
    end
    %Element displacements
    for i=1:NE
        temp=[];
        for j=1:nN
            temp=[temp;
                  NODE(ELEMCon(i,j)).u];
        end
                ELEMENT(i).u=temp;
    end
    %Element Stresses
    SCP = [-1 -1 -1;
           1 -1 -1;
           1 1 -1;
           -1 1 -1;
           -1 -1 1;
           1 -1 1;
           1 1 1;
           -1 1 1];                % Stress calculation points
    for i=1:                                NE
       uelem = ELEMENT(i).u;
       Belem=ELEMENT(i).B;
       Stress = D*Belem*uelem;               %Stress=D*B*u
       strain=Belem*uelem;
       for j=1:size(SCP,1)
           ELEMENT(i).Stress(:,j) =double(subs(Stress,[xi,eta,zeta],SCP(j,:))); 
           ELEMENT(i).strain(:,j) =double(subs(strain,[xi,eta,zeta],SCP(j,:))); 
       end
      
       Strain_x=mean(ELEMENT(i).strain(1,:));
       Stress_x=mean(ELEMENT(i).Stress(1,:));
       [Stress_x, ELEMENT(i).Ex]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_x));
       ELEMENT(i).Actual_stress_x = Stress_x*sign(Strain_x);

       Strain_y=mean(ELEMENT(i).strain(2,:));
       Stress_y=mean(ELEMENT(i).Stress(2,:));
       [Stress_y, ELEMENT(i).Ey]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_y));
        ELEMENT(i).Actual_stress_y = Stress_y*sign(Strain_y);

       Strain_z=mean(ELEMENT(i).strain(3,:));
       Stress_z=mean(ELEMENT(i).Stress(3,:));
       [Stress_z, ELEMENT(i).Ez]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_z));
        ELEMENT(i).Actual_stress_z = Stress_z*sign(Strain_z);

       Strain_xy=mean(ELEMENT(i).strain(4,:));
       Stress_xy=mean(ELEMENT(i).Stress(4,:));
       [Stress_xy, ELEMENT(i).Exy]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_xy));
        ELEMENT(i).Actual_stress_xy = Stress_xy*sign(Strain_xy);

       Strain_yz=mean(ELEMENT(i).strain(5,:));
       Stress_yz=mean(ELEMENT(i).Stress(5,:));
       [Stress_yz, ELEMENT(i).Eyz]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_yz));
        ELEMENT(i).Actual_stress_yz = Stress_yz*sign(Strain_yz);

       Strain_xz=mean(ELEMENT(i).strain(6,:));
       Stress_xz=mean(ELEMENT(i).Stress(6,:));
       [Stress_xz, ELEMENT(i).Exz]=GetHDStressStiffness(E0,SigmaMax,abs(Strain_xz));
       ELEMENT(i).Actual_stress_xz = Stress_xz*sign(Strain_xz);

       Actual_stress_matrix = [ELEMENT(i).Actual_stress_x; ELEMENT(i).Actual_stress_y; ELEMENT(i).Actual_stress_z; ELEMENT(i).Actual_stress_xy; ELEMENT(i).Actual_stress_yz; ELEMENT(i).Actual_stress_xz];
       strain_matrix = [Strain_x; Strain_y; Strain_z; Strain_xy; Strain_yz; Strain_xz];
       if nin == NINC
            writematrix(Actual_stress_matrix, ['ELEMENT_' num2str(i) '_Actual_stress_x.txt']);
            writematrix(strain_matrix, ['ELEMENT_' num2str(i) '_strain_matrix.txt']);
       end
    end
    disp(['      Completed Stress & Strain Computation      ', num2str(toc), ' sec'])
    disp("Finished Loadstep:      " +nin                            )
end

disp("------------------POST-Processing------------------")
%Force Displacement Graph
figure; % Create a new figure for each plot
plot(abs(uDel_inc),uLoad_inc);
xlabel('Displacement')
ylabel('Force')
exportgraphics(gcf,'Force-Displacement_Graph.jpg')

for i = 1:size(Coord, 1)
    Deformed_Node(i,:) = [uDisp(3*i-2,1) uDisp(3*i-1,1) uDisp(3*i,1)];
end
Deformed_structure = Coord-Deformed_Node;
savejpg(Deformed_structure,Quads,'After Deformed')

save('myWorkspace.mat');
disp("-------------------------------------------------------------------------------------------")
disp("----------------------------------------Finished-------------------------------------------")
disp("-------------------------------------------------------------------------------------------")
toc
