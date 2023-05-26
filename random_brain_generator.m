%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%                MRE RANDOM BRAIN GENERATOR               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                         %%%
%%%                       _---~~(~~-_.                      %%%
%%%                     _{        )   )                     %%%
%%%                   ,   ) -~~- ( ,-' )_                   %%%
%%%                  (  `-,_..`., )-- '_,)                  %%%
%%%                 ( ` _)  (  -~( -_ `,  }                 %%%
%%%                 (_-  _  ~_-~~~~`,  ,' )                 %%%
%%%                   `~ -^(    __;-,((()))                 %%%
%%%                         ~~~~ {_ -_(())                  %%%
%%%                                `\  }                    %%%
%%%                                  { }                    %%%
%%%                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dhrubo Jyoti, PhD & Matt McGarry, PhD. Edited 10/14/22  %%%
%%%  Thayer School of Engineering, Dartmouth College, USA   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: Study MRE NLI reconstruction performance in an ensemble of
% simulated brains with randomized property variations of WM tracts.
%
% This script generates run files and submit files for both AP and LR
% actuations for N=20 random brains. Set N to desired number below.
% Three distinct brain models are generated:
% (1) isotropic brain (abbreviated: ISO)
% (2) transversely-isotropic brain with isotropic damping (TIID), and finally
% (3) transversely-isotropic brain with anisotropic damping (TIAD)
% Simply run this script in Matlab 2022a or later
% to output 20 x 3 x 2 = 120 run files and 120 submit files.
% More details in published paper [DOI: ...]
tic;
addpath('inputs');
load('invivo_TIID.mat'); %Use isotropic damping in vivo data
N=20; %Number of brains to create
prefix='102822-brain'; %Use any desired filename. Good to use date.
phi=zeros(18,1); %Index 1, and 13 thru 18, are gray matter, thus zero.
zeta=zeros(18,1); %Brain segmented into 18 segments for this simulation.
phi_all=squeeze(mean(Realphireg,4));
zeta_all=squeeze(mean(Realzetareg,4));
nonans=~isnan(phi_all) & ~isnan(zeta_all); %Just skip over any NaNs in the data.

%% Firstly, GATHER property ranges from available in vivo data
% necessary for building up our brain simulation

% ESTIMATE SHEAR AND TENSILE ANISOTROPY FOR EACH WM TRACT
% Unspecified white matter
jj=sum(tractsLABELS(:,:,:,[ 1:6 9:28  33 34 37:40 43 44 47 48]),4)>0 & nonans ;
phi(2)=mean(phi_all(jj)); zeta(2)=mean(zeta_all(jj));

% Anterior thalamic radiation
jj=sum(tractsLABELS(:,:,:,[29 30]),4)>0 & nonans;
phi(3)=mean(phi_all(jj)); zeta(3)=mean(zeta_all(jj));

% Corticospinal tract
jj=sum(tractsLABELS(:,:,:,[7 8]),4)>0 & nonans;
phi(4)=mean(phi_all(jj)); zeta(4)=mean(zeta_all(jj));

% Cingulate gyrus
jj=sum(tractsLABELS(:,:,:,[35 36]),4)>0 & nonans;
phi(5)=mean(phi_all(jj)); zeta(5)=mean(zeta_all(jj));

% Cingulate hippocampus
jj=sum(tractsLABELS(:,:,:,[37 38]),4)>0 & nonans;
phi(6)=mean(phi_all(jj)); zeta(6)=mean(zeta_all(jj));

% Forceps
jj=sum(tractsLABELS(:,:,:,[5]),4)>0 & nonans;
phi(7)=mean(phi_all(jj)); zeta(7)=mean(zeta_all(jj));

% Inferior fronto-occipital fasciculus
jj=sum(tractsLABELS(:,:,:,[31 32]),4)>0 & nonans;
phi(8)=mean(phi_all(jj)); zeta(8)=mean(zeta_all(jj));

% inferior longitudinal fasciculus
jj=sum(tractsLABELS(:,:,:,[31 32]),4)>0 & nonans;
phi(9)=mean(phi_all(jj)); zeta(9)=mean(zeta_all(jj));

% superior longitudinal fasciculus (regular)
jj=sum(tractsLABELS(:,:,:,[41 42]),4)>0 & nonans;
phi(10)=mean(phi_all(jj)); zeta(10)=mean(zeta_all(jj));

% Uncinate fasciculus
jj=sum(tractsLABELS(:,:,:,[45 46]),4)>0 & nonans;
phi(11)=mean(phi_all(jj)); zeta(11)=mean(zeta_all(jj));

% superior longitudinal fasciculus (temporal)
jj=sum(tractsLABELS(:,:,:,[41 42]),4)>0 & nonans;
phi(12)=mean(phi_all(jj)); zeta(12)=mean(zeta_all(jj)); clear jj;

%% Use moduli and damping ratio population range data from Lucy Hiscox et. al.
% See Excel file [...]
Realmu=zeros(18,2);
Realmu(1,:)=[2.01 2.45]*10^3; %Cortical gray matter
Realmu(2,:)=[2.38 2.80]*10^3; %Unspecified White matter
Realmu(3,:)=[2.72 3.58]*10^3; %Anterior thalamic radiation
Realmu(4,:)=[2.53 3.40]*10^3; %Corticospinal tract
Realmu(5,:)=[2.98 3.86]*10^3; %Cingulate gyrus
Realmu(6,:)=[2.50 3.40]*10^3; %Cingulate hippocampus
Realmu(7,:)=[2.43 3.11]*10^3; %Forceps
Realmu(8,:)=[2.75 3.38]*10^3; %Inferior fronto-occipital fasciulus
Realmu(9,:)=[2.54 3.16]*10^3; %Inferior longitudinal fasciulus
Realmu(10,:)=[2.36 3.19]*10^3; %Superior longitudinal fasciulus (regular)
Realmu(11,:)=[2.56 3.38]*10^3; %Uncinate fasciulus
Realmu(12,:)=[2.36 3.19]*10^3; %Superior longitudinal fasciulus (temporal)
Realmu(13,:)=[2.57 3.63]*10^3; %Thalamus (gray matter)
Realmu(14,:)=[2.28 3.38]*10^3; %Caudate
Realmu(15,:)=[3.15 4.08]*10^3; %Putamen
Realmu(16,:)=[3.06 4.20]*10^3; %Pallidum
Realmu(17,:)=[1.98 3.21]*10^3; %Hippocampus
Realmu(18,:)=[1.92 3.55]*10^3; %Amygdala

Imagmu=zeros(18,2);
Imagmu(1,:)=[.682 1.07]*10^3; %Cortical gray matter
Imagmu(2,:)=[.477 .634]*10^3; %Unspecified White matter
Imagmu(3,:)=[.555 .762]*10^3; %Anterior thalamic radiation
Imagmu(4,:)=[.509 .697]*10^3; %Corticospinal tract
Imagmu(5,:)=[.319 .592 ]*10^3; %Cingulate gyrus
Imagmu(6,:)=[.198 .347]*10^3; %Cingulate hippocampus
Imagmu(7,:)=[.554 .766]*10^3; %Forceps
Imagmu(8,:)=[.545 .775]*10^3; %Inferior fronto-occipital fasciulus
Imagmu(9,:)=[.484 .802 ]*10^3; %Inferior longitudinal fasciulus
Imagmu(10,:)=[.530 .769 ]*10^3; %Superior longitudinal fasciulus
Imagmu(11,:)=[.497 .830 ]*10^3; %Uncinate fasciulus
Imagmu(12,:)=[.530 .769]*10^3; %Superior longitudinal fasciulus
Imagmu(13,:)=[.977 1.44]*10^3; %Thalamus
Imagmu(14,:)=[.998 1.49]*10^3; %Caudate
Imagmu(15,:)=[1.06 1.61]*10^3; %Putamen
Imagmu(16,:)=[1.08 1.61 ]*10^3; %Pallidum
Imagmu(17,:)=[.635 1.16 ]*10^3; %Hippocampus
Imagmu(18,:)=[.423 1.30 ]*10^3; %Amygdala

%% Estimate phi and zeta population range by taking 80% and 120% of single indiviudal in vivo data
phi(:,2)=phi+0.2*abs(phi); %DJ: Generalized to accomodate negative values
zeta(:,2)=zeta+0.2*abs(zeta);
phi(:,1)=phi(:,1)-0.2*abs(phi(:,1));
zeta(:,1)=zeta(:,1)-0.2*abs(zeta(:,1));
%phi(:,2)=phi*1.2;
%zeta(:,2)=zeta*1.2;
%phi(:,1)=phi(:,1)*0.8;
%zeta(:,1)=zeta(:,1)*0.8;
Realmus=Realmu.*(1+phi);
Realmut=Realmu.*(3+4*zeta);

%% Not much is known about aniso damping, so just use full range values of iso damp
Imagmus=Realmus.*Imagmu./Realmu; %Convolve ranges of Realmus and DR from Lucy Excel file
Imagmut=Realmut.*Imagmu./Realmu; 

%% Gathered all necessary data. Now create brains and save fwd run files!
for j=1:N
    %% TIAD model
    nameAP=['TIADfwdrunfile_v9p33_1p0xDTI_AP_',prefix,'-',num2str(j),'.dat'];
    nameLR=['TIADfwdrunfile_v9p33_1p0xDTI_LR_',prefix,'-',num2str(j),'.dat'];
    brain=zeros(18,11); brain(:,8:11)=1; %gray matter
    brain(1:18,1)=1:18; %this column is just numbering
    for row=1:18;
        %uniformly distributed. Increase ranges by 50% to compensate for
        %incomplete contrast recovery
        randRealmu=Realmu(row,1)-((Realmu(row,2)-Realmu(row,1))/4)+rand()*(Realmu(row,2)-Realmu(row,1))*1.5;
        randImagmu=Imagmu(row,1)-((Imagmu(row,2)-Imagmu(row,1))/4)+rand()*(Imagmu(row,2)-Imagmu(row,1))*1.5;
        randRealmus=Realmus(row,1)-((Realmus(row,2)-Realmus(row,1))/4)+rand()*(Realmus(row,2)-Realmus(row,1))*1.5;
        randImagmus=Imagmus(row,1)-((Imagmus(row,2)-Imagmus(row,1))/4)+rand()*(Imagmus(row,2)-Imagmus(row,1))*1.5;
        randRealmut=Realmut(row,1)-((Realmut(row,2)-Realmut(row,1))/4)+rand()*(Realmut(row,2)-Realmut(row,1))*1.5;
        randImagmut=Imagmut(row,1)-((Imagmut(row,2)-Imagmut(row,1))/4)+rand()*(Imagmut(row,2)-Imagmut(row,1))*1.5;

        brain(row,2)=randRealmu;
        brain(row,3)=randImagmu;
        brain(row,4)=randRealmus;
        brain(row,5)=randImagmus;
        brain(row,6)=randRealmut;
        brain(row,7)=randImagmut;
    end
    writematrix(brain,'temp.dat','Delimiter','\t');
    f1=fileread('TIADfwdrunfile_v9p33_1p0xDTI_AP_template.dat'); % TEMPLATE file! Need it.
    f2=fileread('temp.dat');
    [fid,msg]=fopen(['output/',nameAP],'wt');
    fprintf(fid,'%s\n%s',f1,f2); 
    fclose(fid);
    S = readlines(['output/',nameAP]);
    S{13}=regexprep(S{13},'50Hz',['50Hz_',prefix,'-',num2str(j)]);
    [fid, msg] = fopen(['output/',nameAP], 'w');
    if fid < 1;
        error('could not write output file because "%s"', msg);
    end
    fwrite(fid, strjoin(S, '\n')); fclose(fid);

    %Do same thing for LR
    f3=fileread('TIADfwdrunfile_v9p33_1p0xDTI_LR_template.dat'); % TEMPLATE file for LR.
    [fid,msg]=fopen(['output/',nameLR],'wt');
    fprintf(fid,'%s\n%s',f3,f2); 
    fclose(fid);
    S = readlines(['output/',nameLR]);
    S{13}=regexprep(S{13},'50Hz',['50Hz_',prefix,'-',num2str(j)]);
    [fid, msg] = fopen(['output/',nameLR], 'w');
    if fid < 1;
        error('could not write output file because "%s"', msg);
    end
    fwrite(fid, strjoin(S, '\n'));
    fclose(fid);
    delete('temp.dat');

    % Now create the DSLURM files. First AP:
    fid = fopen('inputs/DSLURM_TIADfwd_v9p33_brain_AP','r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);
    % Change cell A
    A{3} = ['#SBATCH --job-name=ADAP',num2str(j)];
    A{27} = ['rfile=',nameAP];

    % Write cell A into txt
    fid = fopen(['output/DSLURM_',nameAP], 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);

    % Same thing for LR
    fid = fopen('inputs/DSLURM_TIADfwd_v9p33_brain_LR','r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);
    % Change cell A
    A{3} = ['#SBATCH --job-name=ADLR',num2str(j)];
    A{27} = ['rfile=',nameLR];

    % Write cell A into txt
    fid = fopen(['output/DSLURM_',nameLR], 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);

    %%  TIID model
    nameAP=['TIIDfwdrunfile_v9p33_1p0xDTI_AP_',prefix,'-',num2str(j),'.dat'];
    nameLR=['TIIDfwdrunfile_v9p33_1p0xDTI_LR_',prefix,'-',num2str(j),'.dat'];
    %brain=zeros(18,11);
    %brain(:,[2 3 8:11])=1;
    brain(:,4:7)=0;
    for row=1:18;
        %uniformly distributed. Increase ranges by 50% to compensate for
        %incomplete contrast recovery
        %randRealmu=Realmu(row,1)-((Realmu(row,2)-Realmu(row,1))/4)+rand()*(Realmu(row,2)-Realmu(row,1))*1.5;
        %randImagmu=Imagmu(row,1)-((Imagmu(row,2)-Imagmu(row,1))/4)+rand()*(Imagmu(row,2)-Imagmu(row,1))*1.5;
        %Use same values as previous model!

        randphi=phi(row,1)-((phi(row,2)-phi(row,1))/4)+rand()*(phi(row,2)-phi(row,1))*1.5;
        randzeta=zeta(row,1)-((zeta(row,2)-zeta(row,1))/4)+rand()*(zeta(row,2)-zeta(row,1))*1.5;
        %brain(row,2)=randRealmu;
       % brain(row,3)=randImagmu;
        brain(row,4)=randphi;
        brain(row,6)=randzeta;
    end
    writematrix(brain,'temp.dat','Delimiter','\t');
    f1=fileread('TIIDfwdrunfile_v9p33_1p0xDTI_AP_template.dat'); % TEMPLATE file! Need it.
    f2=fileread('temp.dat');
    [fid,msg]=fopen(['output/',nameAP],'wt');
    fprintf(fid,'%s\n%s',f1,f2);
    fclose(fid);
    S = readlines(['output/',nameAP]);
    S{13}=regexprep(S{13},'50Hz',['50Hz_',prefix,'-',num2str(j)]);
    [fid, msg] = fopen(['output/',nameAP], 'w');
    if fid < 1;
        error('could not write output file because "%s"', msg);
    end
    fwrite(fid, strjoin(S, '\n'));
    fclose(fid);

    %Do same thing for LR
    f3=fileread('TIIDfwdrunfile_v9p33_1p0xDTI_LR_template.dat'); % TEMPLATE file for LR.
    [fid,msg]=fopen(['output/',nameLR],'wt');
    fprintf(fid,'%s\n%s',f3,f2);
    fclose(fid);
    S = readlines(['output/',nameLR]);
    S{13}=regexprep(S{13},'50Hz',['50Hz_',prefix,'-',num2str(j)]);
    [fid, msg] = fopen(['output/',nameLR], 'w');
    if fid < 1;
        error('could not write output file because "%s"', msg);
    end
    fwrite(fid, strjoin(S, '\n'));
    fclose(fid);
    delete('temp.dat');

    % Now create the DSLURM files. First AP:
    fid = fopen('inputs/DSLURM_TIIDfwd_v9p33_brain_AP','r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);
    % Change cell A
    A{3} = ['#SBATCH --job-name=IDAP',num2str(j)];
    A{27} = ['rfile=',nameAP];

    % Write cell A into txt
    fid = fopen(['output/DSLURM_',nameAP], 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);

    % Same thing for LR
    fid = fopen('inputs/DSLURM_TIIDfwd_v9p33_brain_LR','r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);
    % Change cell A
    A{3} = ['#SBATCH --job-name=IDLR',num2str(j)];
    A{27} = ['rfile=',nameLR];

    % Write cell A into txt
    fid = fopen(['output/DSLURM_',nameLR], 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);

    %%  ISO model
    nameAP=['ISOfwdrunfile_v9p33_1p0xDTI_AP_',prefix,'-',num2str(j),'.dat'];
    nameLR=['ISOfwdrunfile_v9p33_1p0xDTI_LR_',prefix,'-',num2str(j),'.dat'];
    %brain=zeros(18,11);
    %brain(:,[2 3 8:11])=1;
    brain(:,4:7)=0;
    %for row=1:18;
        %uniformly distributed. Increase ranges by 50% to compensate for
        %incomplete contrast recovery
        %randRealmu=Realmu(row,1)-((Realmu(row,2)-Realmu(row,1))/4)+rand()*(Realmu(row,2)-Realmu(row,1))*1.5;
        %randImagmu=Imagmu(row,1)-((Imagmu(row,2)-Imagmu(row,1))/4)+rand()*(Imagmu(row,2)-Imagmu(row,1))*1.5;
        %Use same values as previous model!
        %brain(row,2)=randRealmu;
        %brain(row,3)=randImagmu;
    %end
    writematrix(brain,'temp.dat','Delimiter','\t');
    f1=fileread('ISOfwdrunfile_v9p33_1p0xDTI_AP_template.dat'); % TEMPLATE file! Need it.
    f2=fileread('temp.dat');
    [fid,msg]=fopen(['output/',nameAP],'wt');
    fprintf(fid,'%s\n%s',f1,f2);
    fclose(fid);
    S = readlines(['output/',nameAP]);
    S{13}=regexprep(S{13},'50Hz',['50Hz_',prefix,'-',num2str(j)]);
    [fid, msg] = fopen(['output/',nameAP], 'w');
    if fid < 1;
        error('could not write output file because "%s"', msg);
    end
    fwrite(fid, strjoin(S, '\n'));
    fclose(fid);

    %Do same thing for LR
    f3=fileread('ISOfwdrunfile_v9p33_1p0xDTI_LR_template.dat'); % TEMPLATE file for LR.
    [fid,msg]=fopen(['output/',nameLR],'wt');  fprintf(fid,'%s\n%s',f3,f2);
    fclose(fid);
    S = readlines(['output/',nameLR]);
    S{13}=regexprep(S{13},'50Hz',['50Hz_',prefix,'-',num2str(j)]);
    [fid, msg] = fopen(['output/',nameLR], 'w');
    if fid < 1;
        error('could not write output file because "%s"', msg);
    end
    fwrite(fid, strjoin(S, '\n'));
    fclose(fid);
    delete('temp.dat');

    % Now create the DSLURM files. First AP:
    fid = fopen('inputs/DSLURM_ISOfwd_v9p33_brain_AP','r');
    i = 1;  tline = fgetl(fid);  A{i} = tline;
    while ischar(tline)
        i = i+1;  tline = fgetl(fid);   A{i} = tline;
    end
    fclose(fid);
    % Change cell A
    A{3} = ['#SBATCH --job-name=ISOAP',num2str(j)];
    A{27} = ['rfile=',nameAP];

    % Write cell A into txt
    fid = fopen(['output/DSLURM_',nameAP], 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);

    % Same thing for LR
    fid = fopen('inputs/DSLURM_ISOfwd_v9p33_brain_LR','r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);
    % Change cell A
    A{3} = ['#SBATCH --job-name=ISOLR',num2str(j)];
    A{27} = ['rfile=',nameLR];

    % Write cell A into txt
    fid = fopen(['output/DSLURM_',nameLR], 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);
    disp(['Generated brain ',num2str(j)]);
    fclose('all');

    pause(1);
end
toc;

